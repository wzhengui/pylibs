#!/usr/bin/env python3
'''
Template for running MPI job on sciclone/ches 
Note: For MPI jobs demanding large memory, use small ppn
'''
from pylib import *

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
gpu        =1

#-----------------------------------------------------------------------------
#on front node: 1). submit jobs
#-----------------------------------------------------------------------------
bdir=os.path.abspath(os.path.curdir)
add_var(['ibatch','qnode','qname','account','reservation','jname','scrout','gpu'],
        [1,None,None,None,None,'mpi4py','screen.out',None],locals()) #add default values
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally

if os.getenv('param')==None:
   fmt=0; bcode=sys.argv[0]; os.environ['qnode']=get_qnode(qnode)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt,'param',qname,account,reservation,gpu)
   print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#actually work between this part
#-----------------------------------------------------------------------------
import torch
device=torch.accelerator.current_accelerator().type if torch.accelerator.is_available() else "cpu"
nproc=os.environ['SLURM_CPUS_ON_NODE']; host=os.environ['HOST']
print('nproc={}, device={}: {}'.format(nproc,device,host)); sys.stdout.flush()
