#!/usr/bin/env python3
'''
  extract SCHISM slab outputs
'''
from pylib import *

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/home/g/g260135/work/wangzg/DSP/RUN08a'  #run dir containing outputs
svars=('elev','hvel','GEN_1')                 #variables to be extracted
levels=[1,3,]        #schism level indices (1-nvrt: surface-bottom; (>nvrt): kbp level)
sname='RUN08a/slab'  #name for saving the resutls

#optional
#stacks=[1,5]   #outputs stacks to be extracted
#nspool=12      #sub-sampling frequency within each stack (1 means all)
#mdt=1          #time window (day) for averaging output
#rvars=['elev','hvel','G1'] #rname the varibles 
#reg=None       #region for subsetting reslts (*.reg, or *.bp, or gd_subgrid)

#resource requst 
walltime='00:10:00'; nnode=1;  ppn=4
#hpc: femto, hurricane, bora, vortex, potomac, james, frontera, levante, stampede2
#ppn:   32,       8,     8,    12,       12,     20,     56,      128,      48

#optional: (frontera,levante,stampede2)
qname   ='compute'         #partition name
account ='TG-OCE140024'    #stampede2: NOAA_CSDL_NWI,TG-OCE140024; levante: gg0028 
qnode   =None              #specify node name, or default qnode based on HOST will be used

brun=os.path.basename(run); jname='Rd_'+brun #job name
ibatch=1; scrout='screen_{}.out'.format(brun); bdir=os.path.abspath(os.path.curdir)
#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]; os.environ['qnode']=get_qnode(qnode)
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
sdir=run+'/outputs'                            #output directory
if 'nspool' not in locals(): nspool=1          #subsample
if 'rvars' not in locals(): rvars=svars        #rename variables
if 'mdt' not in locals(): mdt=None             #averaging
if 'reg' not in locals(): reg=None             #region
modules, outfmt, dstacks, dvars, dvars_2d = get_schism_output_info(sdir,1)  #schism outputs information
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks #check stacks

#extract results
irec=0; oname=odir+'/.schout'
for svar in svars: 
   ovars=get_schism_var_info(svar,modules,fmt=outfmt)
   if ovars[0][1] not in dvars: continue 
   for istack in stacks:
       fname='{}_{}_{}_slab'.format(oname,svar,istack); irec=irec+1; t00=time.time()
       if irec%nproc==myrank: 
          read_schism_slab(run,svar,levels,istack,nspool,mdt,fname=fname,reg=reg)
          dt=time.time()-t00; print('finishing reading {}_{}.nc on myrank={}: {:.2f}s'.format(svar,istack,myrank,dt)); sys.stdout.flush()

#combine results
comm.Barrier()
if myrank==0:
   S=zdata(); S.time=[]; fnames=[]
   for i,[k,m] in enumerate(zip(svars,rvars)):
       data=[]; mtime=[]
       for istack in stacks:
           fname='{}_{}_{}_slab.npz'.format(oname,k,istack)
           if not fexist(fname): continue
           C=loadz(fname); data.extend(C.__dict__[k]); mtime.extend(C.time); fnames.append(fname)
       if len(data)>0: S.__dict__[m]=array(data)
       if len(mtime)>len(S.time): S.time=array(mtime)
   savez(sname,S)
   for i in fnames: os.remove(i)

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
