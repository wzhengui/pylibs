#!/usr/bin/env python3
'''
  extract time series for SCHISM outputs
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/sciclone/data10/wangzg/BGC_Tests/Test_ICM_SAV_ChesBay_v1'
svars=('elev','salt','temp','PB1','PB2','PB3','DO','DOC','NH4','NO3','PO4') #variables to be extracted
bpfile='./database/station.bp'  #station file
sname='./outputs/icm'

#optional
#stacks=[1,10] #outputs stacks to be extracted
#ifs=0        #=0: refer to free surface; =1: fixed depth
#nspool=6     #sub-sampling frequency within each stack (1 means all)
#modules=['Hydro','ICM'] #SCHISM modules that output variables belong to 
#rvars=['elev','salt','temp','PB1','PB2','PB3','DO','DOC','NH4','NO3','PO4'] #rname the varibles 

#resource requst 
walltime='00:10:00' 
qnode='x5672'; nnode=2; ppn=8       #hurricane, ppn=8
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=2; ppn=12    #vortex, ppn=12
#qnode='femto'; nnode=1; ppn=32     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='frontera'; nnode=1; ppn=56  #frontera, ppn=56 (flex,normal)
#qnode='levante'; nnode=1; ppn=36   #levante, ppn=128
#qnode='stampede2'; nnode=1; ppn=48 #stampede2, ppn=48 (skx-normal,skx-dev,normal,etc)

#additional information:  frontera,levante,stampede2
qname='flex'                        #partition name
account='TG-OCE140024'              #stampede2: NOAA_CSDL_NWI,TG-OCE140024; levante: gg0028 

jname='Rd_{}'.format(os.path.basename(run)) #job name
ibatch=1; scrout='screen.out'; bdir=os.path.abspath(os.path.curdir)
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
comm=MPI.COMM_WORLD; nproc=comm.Get_size(); myrank=comm.Get_rank()
if myrank==0: 
   t0=time.time(); odir=os.path.dirname(os.path.abspath(sname))
   if not fexist(odir): os.system('mkdir -p {}'.format(odir))

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
sdir=run+'/outputs'                                              #output directory
n2d=len([i for i in os.listdir(sdir) if i.startswith('out2d_')]) #scribe IO or OLDIO 
if 'ifs' not in locals(): ifs=0                                  #refer to free surface
if 'nspool' not in locals(): nspool=1                            #subsample
if 'modules' not in locals(): modules=None                       #module (Hydro,ICM, etc.)
if 'rvars' not in locals(): rvars=svars                          #rename variables

#check available stacks
if 'stacks' in locals(): 
    stacks=arange(stacks[0],stacks[1]+1)
else:
    if n2d==0:
       stacks=sort([int(i.split('.')[0].split('_')[-1]) for i in os.listdir(sdir) if (i.startswith('schout_') and len(i.split('_'))==2)]) 
    else:
       stacks=sort([int(i.split('.')[0].split('_')[-1]) for i in os.listdir(sdir) if i.startswith('out2d_')])

#get 2D variables list
if n2d==0:
   C=ReadNC('{}/schout_{}.nc'.format(sdir,stacks[0]),1); svars_2d=[i for i in C.variables if C.variables[i].ndim==2]; C.close()
else:
   C=ReadNC('{}/out2d_{}.nc'.format(sdir,stacks[0]),1); svars_2d=[*C.variables]; C.close()

#read model grid
if os.path.exists(run+'/grid.npz'):
    gd=loadz('{}/grid.npz'.format(run)).hgrid
else:
    gd=read_schism_hgrid(run+'/hgrid.gr3')
gd.compute_bnd(); sys.stdout.flush()

for itype in [1,2]: 
    irec=0; oname=os.path.dirname(os.path.abspath(sname))+'/.schout'
    for svar in svars: 
       ovars=get_schism_output_info(svar,modules) 
       if itype==1 and (ovars[0][1] not in svars_2d): continue #read 2D outputs 
       if itype==2 and (ovars[0][1] in svars_2d): continue #read 3D outputs 
       for istack in stacks:
           fname='{}_{}_{}'.format(oname,svar,istack); irec=irec+1; t00=time.time()
           if irec%nproc==myrank: 
              read_schism_output_xyz(run,svar,bpfile,istack,ifs,nspool,fname=fname,grid=gd)
              dt=time.time()-t00; print('finishing reading {}_{}.nc on myrank={}: {:.2f}s'.format(svar,istack,myrank,dt)); sys.stdout.flush()

#combine results
comm.Barrier()
if myrank==0:
   S=zdata(); S.time=[]
   for m,[svar,rvar] in enumerate(zip(svars,rvars)): 
       exec('S.{}=[]'.format(rvar)) 
       for istack in stacks:
           fname='{}_{}_{}.npz'.format(oname,svar,istack); C=loadz(fname); os.remove(fname) 
           exec('S.{}.extend(C.{}.transpose([1,0,*arange(2,C.{}.ndim)]))'.format(rvar,svar,svar))
           if m==0: S.time.extend(C.time)
       exec('S.{}=array(S.{}); S.{}=S.{}.transpose([1,0,*arange(2,S.{}.ndim)])'.format(rvar,rvar,rvar,rvar,rvar))
       S.time=array(S.time); S.bp=read_schism_bpfile(bpfile)
   savez(sname,S)

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
