#!/usr/bin/env python3
'''
running MPI job on sciclone/ches to get schism output subset
Note: For MPI jobs demanding large memory, use small ppn
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='./RUN19n' #SCHISM run dir
sdir='/sciclone/pscr/wangzg/tmp' #directory for saving subset data
bxy=array([[110, 130, 130, 110],[-13,-13,7,7]]).T  #ACE regfile or c_[x,y]

#optional outputs
stacks=[1,5]  #outputs stack
#svars=['out2d', 'zCoordinates', 'salinity'] #output files for subset
#grd='grid.npz' #hgrid to speed up

#resource requst 
walltime='00:10:00'; nnode=1;  ppn=4
#hpc: femto, hurricane, bora, vortex, potomac, james, frontera, levante, stampede2
#ppn:   32,       8,     8,    12,       12,     20,     56,      128,      48

#optional: (frontera,levante,stampede2)
qname   ='compute'         #partition name
account ='TG-OCE140024'    #stampede2: NOAA_CSDL_NWI,TG-OCE140024; levante: gg0028 
qnode   =None              #specify node name, or default qnode based on HOST will be used

ibatch=1; scrout='screen.out'; jname='subset'; bdir=os.path.abspath(os.path.curdir)
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
comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: t0=time.time()
if not fexist(sdir): os.mkdir(sdir)

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
#determine default variables,stacks
odir=run+'/outputs'; dmodules,outfmt,dstacks,dvars,dvars_2D=get_schism_output_info(odir)
if outfmt==1: sys.exit('OLDIO is not support yet')
dvars=r_[array(['out2d']),setdiff1d(dvars,dvars_2D)]
if 'svars' not in locals(): svars=dvars
if 'grd' not in locals(): grd=run+'/hgrid.gr3'
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks

#get output subset
gd=loadz(grd).hgrid if grd.endswith('.npz') else read_schism_hgrid(grd)
fn=0; gdn=None
for m,svar in enumerate(svars):
    for n,stack in enumerate(stacks):
        bname='{}_{}.nc'.format(svar,stack); fname='{}/{}'.format(odir,bname); sname='{}/{}'.format(sdir,bname)
        if not fexist(fname): continue
        if fn%nproc==myrank: 
           t00=time.time(); gdn=get_schism_output_subset(fname,sname,bxy,gd)
           print('writing {} on rank {}, {:.1f} s'.format(bname,myrank,time.time()-t00)); sys.stdout.flush()
           if gdn is not None: gd=gdn
        fn=fn+1

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
