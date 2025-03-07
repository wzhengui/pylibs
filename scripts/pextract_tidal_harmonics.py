#!/usr/bin/env python3
'''
perform harmonic analysis for tide in whole domain
'''
from pylib import *
from mpi4py import MPI

#-----------------------------------------------------------------------------
#Input
#hpc: kuro, femto, bora, potomac, james, frontera, levante, stampede2
#ppn:  64,   32,    20,    12,     20,     56,      128,      48
#-----------------------------------------------------------------------------
run='/sciclone/data10/wangzg/CBP/RUN10i'  #run dir containing outputs
tidal_names=['O1','K1','Q1','P1','M2','S2','K2','N2']
sname='elev'  #name for saving the resutls
isave_raw=0    #save original time series on each rank (each rank only for subset of nodes)
#tidal_names=['O1','K1','Q1','P1','M2','S2','K2','N2','M3','M4','M6','M7','M8','M10']

#optional 
#stacks=[1,5]   #outputs stacks to be extracted

#resource requst 
walltime='12:00:00'; nnode=5;  ppn=64

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
sdir=run+'/outputs'                            #output directory
modules, outfmt, dstacks, dvars, dvars_2d = get_schism_output_info(sdir,1)    #schism outputs information
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks   #check stacks

#--------------------------------------------------------
#extract elev on each node and do HA
#--------------------------------------------------------
#distribute jobs on each node
gd=grd(run); npt=int(gd.np/nproc); mtime=[]; elev=[]
sindps=[arange(i*npt,gd.np) if i==(nproc-1) else arange(i*npt,(i+1)*npt) for i in arange(nproc)]; sindp=sindps[myrank]

#read elev on myrank=0, and distribute 
for istack in stacks: 
    if myrank==0:
       t00=time.time()
       C=read_schism_slab(run,'elev',[1,],istack); es=[C.elev[:,i] for i in sindps]; mtime.extend(C.time)
    else:
       es=None

    #collect data on each rank
    elevi=comm.scatter(es,root=0) if ibatch==1 else es[0];  elev.extend(elevi)
    if myrank==0: dt=time.time()-t00; print('reading stack: {}, dt={:0.3f}'.format(istack,dt)); sys.stdout.flush()
mtime=comm.bcast(array(mtime),root=0) if ibatch==1 else array(mtime)
elev=array(elev).T; dt=mean(diff(mtime))

#HA for each pt
C=zdata(); C.time=array(mtime); C.elev=elev;  C.sindp=sindp
C.stime=arange(int(C.time.min()),int(C.time.max())).astype('float'); C.amplitude,C.phase,C.selev=[],[],[]
for i,y0 in enumerate(C.elev):
    tn='.tidal_const_{}'.format(myrank); fn='.tidal_series_{}'.format(myrank); sn='.tidal_consit_{}'.format(myrank)
    H=harmonic_analysis(y0,dt,tidal_names=tidal_names,tname=tn,fname=fn,sname=sn)

    #construct tidal signals, and get subtidal signal
    fy=zeros(len(C.time))
    for k,tname in enumerate(H.tidal_name): 
        if tname=='Z0': continue
        fy=fy+H.amplitude[k]*cos(H.freq[k]*(C.time-C.time[0])*86400-H.phase[k])
    sy=interp(C.time,y0-fy,C.stime)
    if ibatch==0 and i%100==0: print('HA on node {}'.format(i))

    #save HA
    C.amplitude.append(H.amplitude); C.phase.append(H.phase); C.selev.append(sy)
C.to_array('amplitude','phase'); C.to_array('selev',dtype='float32'); C.tidal_name=H.tidal_name; C.freq=H.freq
if isave_raw==1: C.save('{}_{}'.format(sname,myrank)) #save raw data

#combine results
C.delattr('time','elev')  #remove raw data
CS=comm.gather(C,root=0) if ibatch==1 else [C]
if myrank==0:
   svars=['amplitude','phase','selev','sindp']
   S=zdata(); sdict=S.__dict__; [S.attr(i,[]) for i in svars]
   for C in CS: [sdict[i].extend(C.attr(i)) for i in svars] #collect data
   sind=argsort(S.sindp); [S.attr(i,array(S.attr(i))[sind].T) for i in svars]
   [S.attr(i,C.attr(i)) for i in ['freq','stime','tidal_name']] #no combine
   S.delattr('sindp'); S.save(sname)

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
if ibatch==1: comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora','levante'] else os._exit(0)
