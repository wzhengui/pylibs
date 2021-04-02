#!/usr/bin/env python3
'''
Extract (xyz) values from station.bp
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/ches/data10/wangzg/DSP/RUN01a'     #run dir containing outputs
stacks=[1,1]           #stacks of schout_*.nc 
sname='./test'         #name for results
svars=['salt','temp','hvel']        #variable to be extracted
bpfile='./station.bp'  #file name of station.bp
ifs=1                  #ifs=1: depth relative to surface; ifs=0: fixed depth (z coordiante) 
icmb=0                 #icmb=0: work on uncombined; icmb=1: work on combined schout_*.nc

#optional
grid='grid.npz'        #saved grid info; use hgrid.gr3 and vgrid.in if not exist in run
igather=1              #igather=1: save data on each rank,then combine; igather=0: use MPI  

#resource requst 
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=10; ppn=12   #vortex, ppn=12
#qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
qnode='potomac'; nnode=1; ppn=2    #ches, ppn=12
#qnode='james'; nnode=1; ppn=2     #james, ppn=20
#qnode='femto'; nnode=1; ppn=2      #femto,ppn=32, not working yet
#qnode='skylake'; nnode=2; ppn=36    #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2   #viz3,haswell, ppn=24,or 28

walltime='1:00:00'
#-----------------------------------------------------------------------------
#pre-processing
#-----------------------------------------------------------------------------
nproc=nnode*ppn
bdir=os.path.abspath(os.path.curdir)
jname='Rd_{}'.format(os.path.basename(run)) #job name

os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run local
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

#-----------------------------------------------------------------------------
#compute grid and bpfile information
#-----------------------------------------------------------------------------
#read grid information
gname='{}/{}'.format(run,grid); hgrid='{}/hgrid.gr3'.format(run); vgrid='{}/vgrid.in'.format(run)
if os.path.exists(gname):
   gd=loadz(gname).hgrid; vd=loadz(gname).vgrid
else:
   gd=read_schism_hgrid(hgrid); vd=read_schism_vgrid(vgrid)

#read bpfile
bp=read_schism_bpfile(bpfile)

#check pts are inside grid
sind=gd.inside_grid(c_[bp.x,bp.y]); sindn=nonzero(sind==-1)[0]
if sum(sindn)!=0: sys.exit('pts outside of domain: {}'.format(sindn+1))

#compute area coordinate
bp.ie,bp.ip,bp.acor=gd.compute_acor(c_[bp.x,bp.y]); #bp.ne,bp.np=gd.ne,gd.np
bp.dp=gd.dp[bp.ip]; bp.dp0=(bp.dp*bp.acor).sum(axis=1)
if vd.ivcor==1: bp.sigma=vd.sigma[bp.ip]; bp.kbp=vd.kbp[bp.ip]

#read subdomain info
if icmb==0:
   subs=gd.read_prop('{}/outputs/global_to_local.prop'.format(run)).astype('int')[bp.ie]
   isub=unique(subs); sbps=[]; sindes=[]
   for i, isubi in enumerate(isub):
       sinde=nonzero(subs==isubi)[0] #elem index of stations 
     
       #build the iegl and ipgl 
       T=read_schism_local_to_global('{}/outputs/local_to_global_{:04}'.format(run,isubi))
       iegl=dict(zip(T.ielg,arange(T.ne))); ipgl=dict(zip(T.iplg,arange(T.np)))
     
       #compute subdomain ie,ip and acor,dp,z,sigma,kbp
       sbp=npz_data(); #sbp.ne,sbp.np=T.ne,T.np
       sbp.ie=array([iegl[k] for k in bp.ie[sinde]])
       sbp.ip=array([[ipgl[k] for k in n ] for n in bp.ip[sinde]])
       sbp.acor=bp.acor[sinde]; sbp.dp=bp.dp[sinde]; sbp.z=bp.z[sinde]; sbp.nsta=len(sinde) 
       if vd.ivcor==1: sbp.sigma=bp.sigma[sinde]; sbp.kbp=bp.kbp[sinde]
       sbps.append(sbp); sindes.extend(sinde)
   sinds=argsort(array(sindes)) #indices to sort station order
else: 
   isub=[None]; sbps=[bp]; sinds=arange(bp.nsta)
#-----------------------------------------------------------------------------
#extract data on each processor
#-----------------------------------------------------------------------------
#distribute jobs
istacks=[i for i in arange(stacks[0],stacks[1]+1) if i%nproc==myrank]

#initilize data capsule
S=npz_data(); S.time=[]; S.bp=bp; [exec('S.{}=[]'.format(i)) for i in svars]

#extract (x,y,z) value for each stack and each subdomain
for n,istack in enumerate(istacks):
    Si=npz_data(); [exec('Si.{}=[]'.format(m)) for m in svars]
    t00=time.time()
    for m,isubi in enumerate(isub):
        #open schout_*.nc
        if icmb==0: fname='{}/outputs/schout_{:04}_{}.nc'.format(run,isubi,istack)
        if icmb==1: fname='{}/outputs/schout_{}.nc'.format(run,istack)
        if (not os.path.exists(fname)) and icmb==0: sys.exit('not exist: {}'.format(fname))
        C=ReadNC(fname,1); bp=sbps[m]
        
        #read time
        mti=array(C.variables['time'][:])/86400; nt=len(mti); 
        if m==0: S.time.extend(mti)

        #extract elevation -> compute zcor -> vertical interploate
        eis=[]; k1s=[]; k2s=[]; rats=[]  
        for i in arange(nt):
            eii=array(C.variables['elev'][i][bp.ip]) if ('elev' in C.variables) else 0*bp.dp

            #compute zcor
            zii=[]
            for k in arange(3):
                if vd.ivcor==1: ziii=vd.compute_zcor(bp.dp[:,k],eii[:,k],sigma=bp.sigma[:,k,:],kbp=bp.kbp[:,k])
                if vd.ivcor==2: ziii=vd.compute_zcor(bp.dp[:,k],eii[:,k])
                zii.append(ziii)
            zi=(array(zii)*bp.acor.T[...,None]).sum(axis=0).T
            ei=(eii*bp.acor).sum(axis=1); eis.append(ei)
 
            #station depth
            mzi=bp.z.copy()
            if ifs==1: mzi=-mzi+ei

            #interpolation in the vertical
            k1=ones(bp.nsta)*nan; k2=ones(bp.nsta)*nan; rat=ones(bp.nsta)*nan
            fp=mzi<=zi[0];  k1[fp]=0; k2[fp]=0; rat[fp]=0   #bottom
            fp=mzi>=zi[-1]; k1[fp]=(vd.nvrt-1); k2[fp]=(vd.nvrt-1); rat[fp]=1  #surface
            for k in arange(vd.nvrt-1):
                fp=(mzi>=zi[k])*(mzi<zi[k+1])
                k1[fp]=k; k2[fp]=k+1
                rat[fp]=(mzi[fp]-zi[k][fp])/(zi[k+1][fp]-zi[k][fp])
            if sum(isnan(r_[k1,k2,rat]))!=0: sys.exit('check vertical interpolation')
            k1s.append(k1); k2s.append(k2); rats.append(rat)
        eis=array(eis); k1s=array(k1s).astype('int'); k2s=array(k2s).astype('int'); rats=array(rats)

        #compute (x,y,z) for each variables
        Sii=npz_data(); [exec('Sii.{}=[]'.format(mm)) for mm in svars]
        for mm, svar in enumerate(svars):
            ndim=C.variables[svar].ndim; dim=C.variables[svar].shape; dimname=C.variables[svar].dimensions

            data=[]
            for i in arange(nt):
                k1=k1s[i]; k2=k2s[i]; rat=rats[i]

                #get variable values 
                if ('nSCHISM_hgrid_node' in dimname):
                    trii=array(C.variables[svar][i][bp.ip])
                elif ('nSCHISM_hgrid_face' in dimname): 
                    trii=array(C.variables[svar][i][bp.ie])
                else:
                    sys.exit('unknown variable format: {},{}'.format(svar,dim))

                #extend values in the bottom: dims[2] is nvrt
                if ('nSCHISM_vgrid_layers' in dimname):
                   sindp=arange(bp.nsta)
                   for nn in arange(3):
                       kbp=bp.kbp[:,nn]; btri=trii[sindp,nn,kbp]
                       for k in arange(vd.nvrt):
                           fp=k<kbp
                           trii[sindp[fp],nn,k]=btri[fp]

                #horizontal interp
                if ('nSCHISM_hgrid_node' in dimname):
                   if ndim==2: tri=(trii*bp.acor).sum(axis=1)
                   if ndim==3: tri=(trii*bp.acor[...,None]).sum(axis=1)
                   if ndim==4: tri=(trii*bp.acor[...,None,None]).sum(axis=1); rat=rat[:,None]
                else:
                   tri=trii

                #vertical interp
                if ('nSCHISM_vgrid_layers' in dimname):
                   datai=(tri[sindp,k1]*(1-rat)+tri[sindp,k2]*rat)
                else:
                   datai=tri
                data.append(datai)

            #save results
            exec('ds=[1,0,*arange(2,{}-1)]; Sii.{}.extend(array(data).transpose(ds))'.format(ndim,svar))
        for i in svars: exec('Si.{}.extend(Sii.{})'.format(i,i)) 

    #combine istack results
    for i in svars: exec('ds=[1,0,*arange(2,array(Si.{}).ndim)]; S.{}.extend(array(Si.{})[sinds].transpose(ds))'.format(i,i,i)) 
    dt00=time.time()-t00; print('finish reading stack={}; time={:0.2f}s'.format(istack,dt00)); sys.stdout.flush()
S.time=array(S.time); ['S.{}=array(S.{}).astype("float32")'.format(i,i) for i in svars]

#-----------------------------------------------------------------------------
#combine results from all ranks
#-----------------------------------------------------------------------------
if igather==1: save_npz('{}_{}'.format(sname,myrank),S)
comm.Barrier()
if igather==0: sdata=comm.gather(S,root=0)
if igather==1 and myrank==0: sdata=[loadz('{}_{}.npz'.format(sname,i)) for i in arange(nproc)]

if myrank==0:
   S=npz_data(); S.time=[]; S.bp=bp
   for i in svars: exec('S.{}=[]'.format(i))
   for i in arange(nproc):
       Si=sdata[i]; S.time.extend(Si.time)
       for m,svar in enumerate(svars): exec('S.{}.extend(Si.{})'.format(svar,svar))

   #save data        
   S.time=array(S.time); sind=argsort(S.time); S.time=S.time[sind]
   for i in svars: exec('ds=[1,0,*arange(2,array(S.{}).ndim)]; S.{}=array(S.{})[sind].transpose(ds)'.format(i,i,i)) 
   save_npz('{}'.format(sname),S)
   if igather==1: [os.remove('{}_{}.npz'.format(sname,i)) for i in arange(nproc)]

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
if qnode=='x5672' or qnode=='james':
   os._exit(0)
else:
   sys.exit(0)
