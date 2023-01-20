#!/usr/bin/env python3
'''
  extract time series of points@xyz or transects@xy from SCHISM outputs
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/sciclone/data10/wangzg/CBP/RUN04a'
svars=('elev','salt','hvel','NO3') #variables to be extracted
bpfile='./station.bp'  #station file
sname='./icm'

#optional
#itype=1         #0: time series of points @xyz;  1: time series of trasects @xy
#ifs=0           #0: refer to free surface; 1: fixed depth
#stacks=[1,3]    #outputs stacks to be extracted
#nspool=12       #sub-sampling frequency within each stack (1 means all)
#rvars=['elev','salt','hvel','NO3'] #rname the varibles 
#prj=['epsg:26918','epsg:4326']  #projections used to transform coord. in station.bp

#resource requst 
walltime='00:10:00' 
#qnode='x5672'; nnode=2; ppn=8       #hurricane, ppn=8
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=2; ppn=12    #vortex, ppn=12
qnode='femto'; nnode=1; ppn=32     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='frontera'; nnode=1; ppn=56  #frontera, ppn=56 (flex,normal)
#qnode='levante'; nnode=1; ppn=36   #levante, ppn=128
#qnode='stampede2'; nnode=1; ppn=48 #stampede2, ppn=48 (skx-normal,skx-dev,normal,etc)

#additional information:  frontera,levante,stampede2
qname='flex'                        #partition name
account='TG-OCE140024'              #stampede2: NOAA_CSDL_NWI,TG-OCE140024; levante: gg0028 

brun=os.path.basename(run); jname='Rd_'+brun #job name
ibatch=1; scrout='screen_{}.out'.format(brun); bdir=os.path.abspath(os.path.curdir)
#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt=fmt,qname=qname)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
bdir=os.getenv('bdir'); os.chdir(bdir) #enter working dir
odir=os.path.dirname(os.path.abspath(sname))
comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()
if myrank==0 and (not fexist(odir)): os.mkdir(odir)

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
sdir=run+'/outputs'                                              #output directory
if 'itype' not in locals(): itype=0                              #time series or transect
if 'ifs' not in locals(): ifs=0                                  #refer to free surface
if 'nspool' not in locals(): nspool=1                            #subsample
if 'rvars' not in locals(): rvars=svars                          #rename variables
if 'prj' not in locals(): prj=None                               #projections
modules, outfmt, dstacks, dvars, dvars_2d = get_schism_output_info(sdir,1)     #schism outputs information
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks #check stacks

#read model grid
gd=loadz(run+'/grid.npz').hgrid if fexist(run+'/grid.npz') else read_schism_hgrid(run+'/hgrid.gr3')
gd.compute_bnd(); sys.stdout.flush()

#extract results
irec=0; oname=odir+'/.schout'
for svar in svars: 
   ovars=get_schism_var_info(svar,modules,fmt=outfmt)
   if ovars[0][1] not in dvars: continue 
   for istack in stacks:
       fname='{}_{}_{}'.format(oname,svar,istack); irec=irec+1; t00=time.time()
       if irec%nproc==myrank: 
          read_schism_output(run,svar,bpfile,istack,ifs,nspool,fname=fname,grid=gd,fmt=itype,prj=prj)
          dt=time.time()-t00; print('finishing reading {}_{}.nc on myrank={}: {:.2f}s'.format(svar,istack,myrank,dt)); sys.stdout.flush()

#combine results
comm.Barrier()
if myrank==0:
   S=zdata(); S.time=[]; fnames=[]
   for i,[k,m] in enumerate(zip(svars,rvars)):
       data=[]; mtime=[]
       for istack in stacks:
           fname='{}_{}_{}.npz'.format(oname,k,istack)
           if not fexist(fname): continue
           C=loadz(fname); datai=C.__dict__[k]; fnames.append(fname)
           data.extend(datai.transpose([1,0,*arange(2,datai.ndim)])); mtime.extend(C.time)
       if len(data)>0: S.__dict__[m]=array(data).transpose([1,0,*arange(2,array(data).ndim)])
       if len(mtime)>len(S.time): S.time=array(mtime)
   S.bp=read_schism_bpfile(bpfile)
   for pn in ['param','icm','sediment','cosine','wwminput']:
       if fexist('{}/{}.nml'.format(run,pn)): S.__dict__[pn]=read_schism_param('{}/{}.nml'.format(run,pn),3)
   savez(sname,S)
   for i in fnames: os.remove(i)

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
