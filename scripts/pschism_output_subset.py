#!/usr/bin/env python3
'''
running MPI job on sciclone/ches to get schism output subset
Note: For MPI jobs demanding large memory, use small ppn
'''
from pylib import *
from shutil import copyfile
from mpi4py import MPI
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='../STOFS3D' #SCHISM run dir
sdir='./' #directory for saving subset data
bxy='CBP.reg' #ACE regfile or c_[x,y]

#optional outputs
#stacks=[1,5]  #outputs stack
fnames=['out2d',] #output files for subset
svars=['time','dryFlagNode','elevation']

#resource requst 
walltime='01:10:00'; nnode=4;  ppn=16
#optional: (frontera,levante,stampede2)
#ibatch     =1              #0: serial mode;  1: parallel mode
#qnode      =None           #specify node name, or default qnode based on HOST will be used
#qname      =None           #partition name
#account    =None           #account name
#reservation=None           #reservation information
#jname      ='mpi4py'       #job name
#scrout     ='screen.out'   #fname for outout and error

#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
bdir=os.path.abspath(os.path.curdir)
add_var(['ibatch','qnode','qname','account','reservation','jname','scrout'],
        [1,None,None,None,None,'subset','screen.out'],locals()) #add default values
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
if myrank==0:
   t0=time.time()
   if not fexist(sdir): os.mkdir(sdir)
   if not fexist(sdir+'/outputs'): os.mkdir(sdir+'/outputs')
if ibatch==1: comm.Barrier()

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
#determine default variables,stacks
odir=run+'/outputs'; dmodules,outfmt,dstacks,dvars,dvars_2D=get_schism_output_info(odir)
if outfmt==1: sys.exit('OLDIO is not support yet')
dvars=r_[array(['out2d']),setdiff1d(dvars,dvars_2D)]
if 'fnames' not in locals(): fnames=dvars
if 'svars' not in locals(): svars=None
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks

#get output subset
gd=grd(run)
fn=0; gdn=None
for m,fname in enumerate(fnames):
    for n,stack in enumerate(stacks):
        bname='{}_{}.nc'.format(fname,stack); oname='{}/{}'.format(odir,bname); sname='{}/outputs/{}'.format(sdir,bname)
        if not fexist(oname): continue
        if fn%nproc==myrank: 
           t00=time.time(); gdn=get_schism_output_subset(oname,sname,bxy,svars,gd)
           print('writing {} on rank {}, {:.1f} s'.format(bname,myrank,time.time()-t00)); sys.stdout.flush()
           if gdn is not None: gd=gdn
        fn=fn+1

#write new grid
if myrank==0:
   vd=zdata(); vd.nvrt=gd.nvrt; vd.kbp=gd.kbp; pn0='{}/param.nml'.format(run); pn='{}/param.nml'.format(sdir)
   c=zdata(); c.hgrid=gdn; c.vgrid=vd;  c.save('{}/grid.npz'.format(sdir))
   if fexist(pn0): copyfile(pn0,pn)

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
