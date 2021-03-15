#!/usr/bin/env python3
'''
Compute fluxes based on SCHISM node information
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='RUN01b'
jname='Rd_{}'.format(run) #job name
stacks=[1,74]
sname='flux'
svars=['GEN_1','GEN_2','GEN_3','GEN_4','GEN_5']

tnames=['SCR_1','SCR_2'] #trasect names
txys=array([[[630597.229, 4257510.76],[630752.001, 4257544.77]],
            [[630168.278, 4234802.51],[630272.344, 4234842.0]]  ]);
dx=10  #interval (m) of the transects

walltime='2:00:00'

#resource requst 
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=10; ppn=12   #vortex, ppn=12
#qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
qnode='james'; nnode=16; ppn=1     #james, ppn=20
#qnode='femto'; nnode=1; ppn=2      #femto,ppn=32, not working yet
#qnode='skylake'; nnode=2; ppn=36    #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2   #viz3,haswell, ppn=24,or 28

#-----------------------------------------------------------------------------
#pre-processing
#-----------------------------------------------------------------------------
nproc=nnode*ppn
bdir=os.path.abspath(os.path.curdir)

#os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run local
#-----------------------------------------------------------------------------
#on front node; submit jobs in this section
#-----------------------------------------------------------------------------
if os.getenv('param')==None and os.getenv('job_on_node')==None:
    args=sys.argv
    param=[bdir,args[0]]
    
    #submit job on node
    if qnode=='femto': 
        scode='sbatch --export=param="{} {}" -J {} -N {} -n {} -t {} {}'.format(*param,jname,nnode,nproc,walltime,args[0])
    else:
        scode='qsub {} -v param="{} {}", -N {} -j oe -l nodes={}:{}:ppn={} -l walltime={}'.format(args[0],*param,jname,nnode,qnode,ppn,walltime)
    print(scode); os.system(scode)
    os._exit(0)

#-----------------------------------------------------------------------------
#still on front node, but in batch mode; running jobs in this section
#-----------------------------------------------------------------------------
if os.getenv('param')!=None and os.getenv('job_on_node')==None:
    param=os.getenv('param').split();
    param=[int(i) if i.isdigit() else i for i in param]
    bdir=param[0]; bcode=param[1]
    os.chdir(bdir)

    if qnode=='bora':
       rcode="mpiexec -x job_on_node=1 -x bdir='{}' -n {} {} >& screen.out".format(bdir,nproc,bcode)
    elif qnode=='femto':
       pypath='/sciclone/home10/wangzg/bin/pylibs/Scripts/:/sciclone/home10/wangzg/bin/pylibs/Utility/'
       rcode="srun --export=job_on_node=1,bdir='{}',PYTHONPATH='{}' {} >& screen.out".format(bdir,pypath,bcode)
    elif qnode=='x5672' or qnode=='vortex' or qnode=='potomac' or qnode=='james':
       rcode="mvp2run -v -e job_on_node=1 -e bdir='{}' {} >& screen.out".format(bdir,bcode)
    elif qnode=='skylake' or qnode=='haswell':
       rcode="mpiexec --env job_on_node 1 --env bdir='{}' -np {} {} >& screen.out".format(bdir,nproc,bcode)
    print(rcode); os.system(rcode); sys.stdout.flush()
    os._exit(0)

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
#link results
if myrank==0:
    if not os.path.exists(run): os.mkdir(run)
    rnames_1=['schout_?.nc','schout_??.nc','schout_???.nc']; [os.system('cd {}; ln -sf ../../{}/outputs/{} ./'.format(run,run,i)) for i in rnames_1]
    rnames_2=['hgrid.npz','vgrid.npz']; [os.system('cd {}; ln -sf ../../{}/{} ./'.format(run,run,i)) for i in rnames_2]
comm.Barrier()

#-----------------------------------------------------------------------------
#compute transect information
#-----------------------------------------------------------------------------
if myrank==0:
   gd=loadz('{}/hgrid.npz'.format(run)).hgrid; gd.compute_ctr()
   vd=loadz('{}/vgrid.npz'.format(run)).vgrid

   #for each transect
   T=npz_data(); T.nt=len(tnames); T.tnames=tnames; T.transects=dict(); T.npt=0; T.sind=[];
   T.x=[]; T.y=[]; T.z=[]; T.dist=[]; T.dx=[]; T.angle=[]; T.ie=[]; T.ip=[]; T.acor=[]; T.sigma=[]
   for m,tname in enumerate(tnames):
       
       #check pts are inside grid
       sxy=txys[m]; sind=gd.inside_grid(sxy)
       if sum(sind==-1)!=0: sys.exit('pts outside of domain: {},{} {}'.format(m,sxy,sind))
       
       #compute xy of transect
       sx,sy=sxy.T
       tdist=abs(diff(sx)+1j*diff(sy))[0]; npt=int(around(tdist/dx))+1; dx=tdist/npt
       sxi=linspace(*sx,npt); syi=linspace(*sy,npt); dist=linspace(0,tdist,npt)
       
       #compute angle for each subsection
       angle=array([arctan2((syi[i+1]-syi[i]),(sxi[i+1]-sxi[i])) for i in arange(npt-1)])
       
       #compute area coordinate
       pie,pip,pacor=gd.compute_acor(c_[sxi,syi]);
       szi=(gd.dp[pip]*pacor).sum(axis=1); sigma=(vd.sigma[pip]*pacor[...,None]).sum(axis=1)
       
       #save transect information
       S=npz_data(); S.npt=npt; S.x=sxi; S.y=syi; S.z=szi; S.dist=dist; S.dx=dx
       S.angle=angle; S.ie=pie; S.ip=pip; S.acor=pacor; S.sigma=sigma

       T.transects[tname]=S; T.sind.append([*arange(T.npt,T.npt+npt)]); T.npt=T.npt+npt 
       T.x.extend(sxi); T.y.extend(syi); T.z.extend(szi); T.dist.extend(dist); T.dx.append(dx)
       T.angle.append(angle); T.ie.extend(pie); T.ip.extend(pip); T.acor.extend(pacor); T.sigma.extend(sigma)  
   
   #save information
   save_npz('{}/{}'.format(run,sname),T)
   
   #plot for check
   #gd.plot_bnd(); [plot(transects[i].x,transects[i].y,'r-') for i in tnames]; 
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
comm.Barrier()
    
#read hgrid and vgrid
os.chdir(run); gd=loadz('hgrid.npz').hgrid; vd=loadz('vgrid.npz').vgrid

#distribute jobs
istacks=[i for i in arange(stacks[0],stacks[1]+1) if i%nproc==myrank]

#read in trasect information
S=loadz('{}.npz'.format(sname))

#compute flux
S.time=[]; S.flux=[]; [exec('S.flux_{}=[]'.format(i)) for i in svars]
for n,istack in enumerate(istacks):
    if not os.path.exists('schout_{}.nc'.format(istack)): continue
    C=ReadNC('schout_{}.nc'.format(istack),1)
    mti=array(C.variables['time'][:])/86400; nt=len(mti)

    #compute zcor and vel, and volume flux
    zis=[]; uis=[]; vis=[]; flux=[]; t00=time.time()
    for i in arange(nt):
        ei=(array(C.variables['elev'][i][S.ip])*S.acor).sum(axis=1)
        uvi=(array(C.variables['hvel'][i][S.ip])*S.acor[...,None,None]).sum(axis=1)
        fpn=uvi>1e10; uvi[fpn]=0.0
        zi=compute_zcor(S.sigma,S.z,ei); ui=uvi[:,:,0]; vi=uvi[:,:,1]
        zis.append(zi); uis.append(ui); vis.append(vi)

        #compute flux
        fluxi=[]
        for k,tname in enumerate(tnames):
            sind=S.sind[k]
            dzi=diff(zi[sind],axis=1); fpz=dzi<0; dzi[fpz]=0.0
            uii=ui[sind]; vii=vi[sind]; angle=S.angle[k]; dx=S.dx[k]
            fluxii=(-sin(angle)[:,None]*dzi[:-1,:]*(uii[:-1,:-1]+uii[:-1,1:]) \
                   -cos(angle)[:,None]*dzi[:-1,:]*(vii[:-1,:-1]+vii[:-1,1:]) \
                   -sin(angle)[:,None]*dzi[1:,:]*(uii[1:,:-1]+uii[1:,1:])    \
                   -cos(angle)[:,None]*dzi[1:,:]*(vii[1:,:-1]+vii[1:,1:])).sum()*dx/4

            fluxi.append(fluxii)
        flux.append(fluxi)
    S.time.extend(mti); S.flux.extend(flux)
    dt00=time.time()-t00; print('finish reading stack={},var=hvel; time={:.2f}s'.format(istack,dt00)); sys.stdout.flush()
   
    #compute fluxes for  each variables
    for m, svar in enumerate(svars):
        tflux=[]; t00=time.time()
        for i in arange(nt):
            ui=uis[i]; vi=vis[i]; zi=zis[i]
            #get tracer values 
            tri=(array(C.variables[svar][i][S.ip])*S.acor[...,None]).sum(axis=1)
            fpn=tri>1e10; tri[fpn]=0.0
  
            #compute flux
            tfluxi=[]
            for k,tname in enumerate(tnames):
                sind=S.sind[k]
                tui=tri[sind]*ui[sind]; tvi=tri[sind]*vi[sind]; angle=S.angle[k]; dx=S.dx[k]
                dzi=diff(zi[sind],axis=1); fpz=dzi<0; dzi[fpz]=0.0 
                tfluxii=(-sin(angle)[:,None]*dzi[:-1,:]*(tui[:-1,:-1]+tui[:-1,1:]) \
                       -cos(angle)[:,None]*dzi[:-1,:]*(tvi[:-1,:-1]+tvi[:-1,1:]) \
                       -sin(angle)[:,None]*dzi[1:,:]*(tui[1:,:-1]+tui[1:,1:])    \
                       -cos(angle)[:,None]*dzi[1:,:]*(tvi[1:,:-1]+tvi[1:,1:])).sum()*dx/4

                tfluxi.append(tfluxii)
            tflux.append(tfluxi)
        
        #save results
        exec('S.flux_{}.extend(tflux)'.format(svar))
        dt00=time.time()-t00; print('finish reading stack={},var={}; time={:0.2f}s'.format(istack,svar,dt00)); sys.stdout.flush()

#combine
comm.Barrier()
sdata=comm.gather(S,root=0)

if myrank==0:
   S=loadz('{}.npz'.format(sname))
   S.time=[]; S.flux=[]; [exec('S.flux_{}=[]'.format(i)) for i in svars]
   for i in arange(nproc):
       Si=sdata[i]
       S.time.extend(Si.time); S.flux.extend(Si.flux)
       for m,svar in enumerate(svars):
           exec('S.flux_{}.extend(Si.flux_{})'.format(svar,svar))

   #save data        
   S.time=array(S.time); sind=argsort(S.time); S.time=S.time[sind]; S.flux=array(S.flux)[sind].T
   [exec('S.flux_{}=array(S.flux_{})[sind].T'.format(i,i)) for i in svars]
   save_npz('{}'.format(sname),S)

   #clean
   #[os.system('rm {}'.format(i)) for i in [*rnames_1,*rnames_2]]
   [os.system('rm {}'.format(i)) for i in [*rnames_1]]

print('myrank={}, nproc={}, host={}'.format(myrank,nproc,os.getenv('HOST'))); sys.stdout.flush()

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
MPI.Finalize()
if qnode=='x5672' or qnode=='james':
   os._exit(0)
else:
   sys.exit(0)
