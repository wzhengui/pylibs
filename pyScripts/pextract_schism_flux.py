#!/usr/bin/env python3
'''
Compute fluxes based on SCHISM node information
'''
from pylib import *

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/home/g/g260135/work/wangzg/DSP/RUN08a'
svars=['GEN_1','GEN_2','GEN_3',]
txy=[[[630597.229,630752.001], [4257510.76,4257544.77]],      #1st transect: [xi,yi]
     'transect_2.bp',                                         #2nd transect; ts2.bp
     [[519390, 522149, 525018],[4186127, 4180971, 4177533]] ] #3rd transect
sname='RUN08a/flux'

#optional
#stacks=[1,9]    #output stacks
#nspool=12       #sub-sampling frequency within each stack (1 means all)
#dx=10          #interval of sub-section, used to divide transect
#prj='cpp'      #projection that convert lon&lat to local project when ics=2
#rvars=['g1','g2','g3',] #rname the varibles

#resource requst
walltime='00:10:00'; nnode=1;  ppn=4
#hpc: femto, hurricane, bora, vortex, potomac, james, frontera, levante, stampede2
#ppn:   32,       8,     8,    12,       12,     20,     56,      128,      48

#optional: (frontera,levante,stampede2)
qname   ='compute'         #partition name
account ='TG-OCE140024'    #stampede2: NOAA_CSDL_NWI,TG-OCE140024; levante: gg0028
qnode   =None              #specify node name, or default qnode based on HOST will be used

brun=os.path.basename(run); jname='Rd_'+brun #job name
ibatch=1; scrout='screen.out'; bdir=os.path.abspath(os.path.curdir)
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

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
if 'nspool' not in locals(): nspool=1       #subsample
if 'rvars' not in locals(): rvars=svars     #rename variables
modules, outfmt, dstacks, dvars, dvars_2d=get_schism_output_info(run+'/outputs',1) #schism outputs info
stacks=arange(stacks[0],stacks[1]+1) if ('stacks' in locals()) else dstacks #check stacks

#check format of transects
if isinstance(txy,str) or array(txy[0]).ndim==1: txy=[txy]
rdp=read_schism_bpfile; txy=[[rdp(i).x,rdp(i).y] if isinstance(i,str) else i for i in txy]

#read grid
fgz=run+'/grid.npz'; fgd=run+'/hgrid.gr3'; fvd=run+'/vgrid.in'
gd=loadz(fgz,'hgrid') if fexist(fgz) else read_schism_hgrid(fgd)
vd=loadz(fgz,'vgrid') if fexist(fgz) else read_schism_vgrid(fvd)

#compute transect information
nps,dsa=[],[]; sx,sy,sinds,angles=[],[],[],[]; ns=len(txy); pxy=ones(ns).astype('O'); ipt=0
for m,[x0,y0] in enumerate(txy):
    #compute transect pts
    x0=array(x0); y0=array(y0)
    if 'dx' in locals():  #divide transect evenly
       ds=abs(diff(x0+1j*y0)); s=cumsum([0,*ds]); npt=int(s[-1]/dx)+1; ms=linspace(0,s[-1],npt)
       xi=interpolate.interp1d(s,x0)(ms); yi=interpolate.interp1d(s,y0)(ms)
    else:
       xi,yi=x0,y0;
    npt=len(xi); ds=abs(diff(xi+1j*yi))
    if sum(gd.inside_grid(c_[xi,yi])==0)!=0: sys.exit('pts outside of domain: {}'.format(m))
    if 'prj' in locals(): pxi,pyi=proj_pts(xi,yi,'epsg:4326',prj); ds=abs(diff(pxi+1j*pyi))

    #transect property
    angle=array([arctan2(yi[i+1]-yi[i],xi[i+1]-xi[i]) for i in arange(npt-1)])   #angle for each subsection
    #pie,pip,pacor=gd.compute_acor(c_[xi,yi]); sigma=(vd.sigma[pip]*pacor[...,None]).sum(axis=1)  #sigma coord.

    nps.append(npt); pxy[m]=c_[xi,yi].T; dsa.append(ds); sx.extend(xi); sy.extend(yi)
    sinds.append(arange(ipt,ipt+npt)); angles.append(angle); ipt=ipt+npt

#compute flux
S=zdata(); S.time=[]; S.flux=[[] for i in txy]; S.tflux=[[[] for i in txy] for i in svars]
for istack in stacks:
    if istack%nproc!=myrank: continue
    t00=time.time(); C=read_schism_output(run,['zcor','hvel',*svars],c_[sx,sy],istack,nspool=nspool,hgrid=gd,vgrid=vd,fmt=1) #read profile
    for m,npt in enumerate(nps): #for each transect
        sind=sinds[m]; angle=angles[m][:,None,None]; ds=dsa[m][:,None,None]
        dz=diff(C.zcor[sind],axis=2); dz=(dz[:-1]+dz[1:])/2; dz[dz<0]=0
        U=C.hvel[sind]; U=(U[:-1,:,:-1]+U[:-1,:,1:]+U[1:,:,:-1]+U[1:,:,1:])/4; u,v=U[...,0],U[...,1]

        #volume flux, and tracers flux
        flx=((sin(angle)*u-cos(angle)*v)*dz*ds).sum(axis=0).sum(axis=1); S.flux[m].extend(flx)
        for n,svar in enumerate(svars):
            T=C.__dict__[svar][sind]; tr=(T[:-1,:,:-1]+T[:-1,:,1:]+T[1:,:,:-1]+T[1:,:,1:])/4
            flx=((sin(angle)*u-cos(angle)*v)*dz*ds*tr).sum(axis=0).sum(axis=1); S.tflux[n][m].extend(flx)
    S.time.extend(C.time); C=None
    print('reading stack {} on rank {}: {:0.2f}'.format(istack,myrank,time.time()-t00)); sys.stdout.flush()
S.time,S.flux,S.tflux=array(S.time),array(S.flux).T,array(S.tflux).T

#gather flux for all ranks
data=comm.gather(S,root=0)
C=zdata(); C.nps=array(nps); C.xy=pxy; C.xy0=txy; C.time,C.flux=[],[]; tflux=[]
if myrank==0:
   for i in data: C.time.extend(i.time); C.flux.extend(i.flux); tflux.extend(i.tflux)
   it=argsort(C.time); C.time=array(C.time)[it]; C.flux=array(C.flux)[it].T.astype('float32')
   for i,rvar in enumerate(rvars): C.__dict__['flux_'+rvar]=array(tflux)[it][...,i].T.astype('float32')

   #save
   sdir=os.path.dirname(os.path.abspath(sname))
   if not fexist(sdir): os.system('mkdir -p '+sdir)
   savez(sname,C)

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora'] else os._exit(0)
