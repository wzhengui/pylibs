#!/usr/bin/env python3
'''
  extract time series of points@xyz or transects@xy from SCHISM outputs
'''
from pylib import *
from mpi4py import MPI

#-----------------------------------------------------------------------------
#Input
#hpc: kuro, femto, bora, potomac, james, frontera, levante, stampede2
#ppn:  64,   32,    20,    12,     20,     56,      128,      48
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
#mdt=1           #time window (day) for averaging output
#rvars=['elev','salt','hvel','NO3'] #rname the varibles 
#prj=['epsg:26918','epsg:4326']  #projections used to transform coord. in station.bp

#hpc resource requst
walltime='00:10:00'; nnode=1;  ppn=4

#optional: (frontera,levante,stampede2,etc.)
ibatch     =1              #0: serial mode;  1: parallel mode
qnode      =None           #specify node name, or default qnode based on HOST will be used
qname      =None           #partition name
account    =None           #account name
reservation=None           #reservation information

#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
brun=os.path.basename(run); jname='Rd_'+brun; scrout='screen_{}.out'.format(brun); bdir=os.path.abspath(os.path.curdir)
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]; os.environ['qnode']=get_qnode(qnode)
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt,'param',qname,account,reservation)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
bdir=os.getenv('bdir'); os.chdir(bdir) #enter working dir
odir=os.path.dirname(os.path.abspath(sname))
if ibatch==0: nproc=1; myrank=0
if ibatch==1: comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()
if myrank==0 and (not fexist(odir)): os.mkdir(odir)

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
add_var(['itype','ifs','nspool','rvars','prj','mdt'],[0,0,1,svars,None,None],locals()) #add default values
modules, outfmt, dstacks, dvars, dvars_2d = schout_info(run+'/outputs',1)   #schism outputs information
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks #check stacks
gd,vd=grd(run,fmt=2); gd.compute_bnd()                                      #read model grid

#extract results
irec=0; oname=odir+'/.schout_'+os.path.basename(os.path.abspath(sname))
for svar in svars: 
   ovars=schvar_info(svar,modules,fmt=outfmt)
   if ovars[0][1] not in dvars: continue 
   for istack in stacks:
       fname='{}_{}_{}'.format(oname,svar,istack); irec=irec+1; t00=time.time()
       if irec%nproc==myrank: 
          try:
             read_schism_output(run,svar,bpfile,istack,ifs,nspool,fname=fname,hgrid=gd,vgrid=vd,fmt=itype,prj=prj,mdt=mdt)
             dt=time.time()-t00; print('finishing reading {}_{}.nc on myrank={}: {:.2f}s'.format(svar,istack,myrank,dt)); sys.stdout.flush()
          except:
             pass

#combine results
if ibatch==1: comm.Barrier()
if myrank==0:
   S=zdata(); S.bp=read(bpfile); S.time=[]; fnss=[]
   for i,[k,m] in enumerate(zip(svars,rvars)):
       fns=['{}_{}_{}.npz'.format(oname,k,n) for n in stacks]; fnss.extend(fns)
       data=[read(fn,k).astype('float32') for fn in fns if fexist(fn)]; mtime=[read(fn,'time') for fn in fns if fexist(fn)]
       if len(data)>0: S.attr(m,concatenate(data,axis=1)); mtime=concatenate(mtime)
       if len(mtime)>len(S.time): S.time=array(mtime)
   [S.attr(pn,read('{}/{}.nml'.format(run,pn),3)) for pn in ['param','icm','sediment','cosine','wwminput'] if fexist('{}/{}.nml'.format(run,pn))]
   S.save(sname); [os.remove(fn) for fn in fnss if fexist(fn)] 

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
