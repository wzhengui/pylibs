#!/usr/bin/env python3
'''
Template for running MPI job on sciclone/ches 
Note: For MPI jobs demanding large memory, use small ppn
'''
from pylib import *
from mpi4py import MPI

#-----------------------------------------------------------------------------
#Input
#hpc: kuro, femto, bora, potomac, james, frontera, levante, stampede2
#ppn:  64,   32,    20,    12,     20,     56,      128,      48
#-----------------------------------------------------------------------------
walltime='00:10:00'; nnode=1;  ppn=4

#optional: (frontera,levante,stampede2)
#ibatch     =1              #0: serial mode;  1: parallel mode
#qnode      =None           #specify node name, or default qnode based on HOST will be used
#qname      =None           #partition name
#account    =None           #account name
#reservation=None           #reservation information
#jname      ='mpi4py'       #job name
#scrout     ='screen.out'   #fname for outout and error

#-----------------------------------------------------------------------------
#on front node: 1). submit jobs  2) running parallel jobs
#-----------------------------------------------------------------------------
bdir=os.path.abspath(os.path.curdir)
add_var(['ibatch','qnode','qname','account','reservation','jname','scrout'],
        [1,None,None,None,None,'mpi4py','screen.out'],locals()) #add default values
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
if ibatch==0: nproc=1; myrank=0
if ibatch==1: comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
print('myrank={}, nproc={}, host={}'.format(myrank,nproc,os.getenv('HOST'))); sys.stdout.flush()

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
