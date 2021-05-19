#!/usr/bin/env python3
'''
Template for running MPI job on sciclone/ches @VIMS
Note: For MPI jobs demanding large memory, use small ppn
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
jname='mpi4py' #job name
walltime='00:10:00' 

#resource requst 
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=10; ppn=12   #vortex, ppn=12
qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='femto'; nnode=1; ppn=2      #femto,ppn=32
#qnode='skylake'; nnode=2; ppn=36    #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2   #viz3,haswell, ppn=24,or 28

ibatch=1
#-----------------------------------------------------------------------------
#pre-processing
#-----------------------------------------------------------------------------
nproc=nnode*ppn; bdir=os.path.abspath(os.path.curdir); scrout='screen.out'

if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run local
#-----------------------------------------------------------------------------
#on front node; submit jobs in this section
#-----------------------------------------------------------------------------
if os.getenv('param')==None and os.getenv('job_on_node')==None:
    args=sys.argv; param=[bdir,args[0]]
    
    #submit job on node
    if qnode=='femto': 
        scode='sbatch --export=param="{} {}" -J {} -N {} -n {} -t {} {}'.format(*param,jname,nnode,nproc,walltime,args[0])
    else:
        scode='qsub {} -v param="{} {}", -N {} -j oe -l nodes={}:{}:ppn={} -l walltime={}'.format(args[0],*param,jname,nnode,qnode,ppn,walltime)
    print(scode); os.system(scode); os._exit(0)

#-----------------------------------------------------------------------------
#still on front node, but in batch mode; running jobs in this section
#-----------------------------------------------------------------------------
if os.getenv('param')!=None and os.getenv('job_on_node')==None:
    param=os.getenv('param').split();
    param=[int(i) if i.isdigit() else i for i in param]
    bdir=param[0]; bcode=param[1]
    os.chdir(bdir)

    if qnode=='femto':
       rcode="srun --export=ALL,job_on_node=1,bdir={} {} >& {}".format(bdir,bcode,scrout)
    elif qnode=='bora':
       rcode="mpiexec -x job_on_node=1 -x bdir='{}' -n {} {} >& {}".format(bdir,nproc,bcode,scrout)
    elif qnode=='x5672' or qnode=='vortex' or qnode=='potomac' or qnode=='james':
       rcode="mvp2run -v -e job_on_node=1 -e bdir='{}' {} >& {}".format(bdir,bcode,scrout)
    elif qnode=='skylake' or qnode=='haswell':
       rcode="mpiexec --env job_on_node 1 --env bdir='{}' -np {} {} >& {}".format(bdir,nproc,bcode,scrout)
    print(rcode); os.system(rcode); sys.stdout.flush(); os._exit(0)

#-----------------------------------------------------------------------------
#on computation node
#-----------------------------------------------------------------------------
#enter working dir
bdir=os.getenv('bdir'); os.chdir(bdir)

#get nproc and myrank
comm=MPI.COMM_WORLD
nproc=comm.Get_size()
myrank=comm.Get_rank()
if myrank==0: t0=time.time()

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
print('myrank={}, nproc={}, host={}'.format(myrank,nproc,os.getenv('HOST'))); sys.stdout.flush()

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
if qnode=='x5672' or qnode=='james':
   os._exit(0)
else:
   sys.exit(0)
