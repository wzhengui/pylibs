#!/usr/bin/env python3
'''
Extract SCHISM-FABM variable values at (x,y,z) from station.bp. 
  1). work for both uncombined and combined SCHISM outputs
  2). can extract multiple variables at the same time
  3). can work in interactive or batch mode 
  4). output in ACSII or *npz format 
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/sciclone/data10/wangzg/fabm_dev/RUN10a'     #run dir containing outputs
stacks=[1,73]
sname='./fabm_station' #name for results
svars=['cosine_MTZ1','cosine_MIDN']          #variable to be extracted
bpfile=run+'/station.bp'  #file name of station.bp
ifs=1                   #ifs=1: depth relative to surface; ifs=0: fixed depth (z coordiante)
fmt=0                   #fmt=0: output as *.npz format; fmt=1: output as ASCII

#optional
grid=run+'/grid.npz'  #saved grid info, to speed up; use hgrid.gr3 and vgrid.in if not exist
igather=1          #igather=1: save data on each rank,then combine; igather=0: use MPI  

#resource requst 
walltime='00:10:00'
#qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=2; ppn=12    #vortex, ppn=12
#qnode='femto'; nnode=2; ppn=12     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='skylake'; nnode=2; ppn=36   #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2    #viz3,haswell, ppn=24,or 28

jname='Rd_{}'.format(os.path.basename(run)) #job name
ibatch=1; scrout='screen.out'; bdir=os.path.abspath(os.path.curdir)
#-----------------------------------------------------------------------------
#on front node: 1). submit jobs first (qsub), 2) running parallel jobs (mpirun) 
#-----------------------------------------------------------------------------
if ibatch==0: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run locally
if os.getenv('job_on_node')==None:
   if os.getenv('param')==None: fmt=0; bcode=sys.argv[0]
   if os.getenv('param')!=None: fmt=1; bdir,bcode=os.getenv('param').split(); os.chdir(bdir)
   scode=get_hpc_command(bcode,bdir,jname,qnode,nnode,ppn,walltime,scrout,fmt=fmt)
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
if myrank==0:
   sdir=os.path.dirname(sname)
   if (not os.path.exists(sdir)) and sdir!='': os.mkdir(sdir)

#-----------------------------------------------------------------------------
#compute grid and bpfile information
#-----------------------------------------------------------------------------
#read grid information
t00=time.time()
if os.path.exists(grid):
   gd=loadz(grid).hgrid; vd=loadz(grid).vgrid
else:
   gd=read_schism_hgrid('{}/hgrid.gr3'.format(run))
   vd=read_schism_vgrid('{}/vgrid.in'.format(run))

#compute area coordinate for stations
bp=read_schism_bpfile(bpfile)
bp.ie,bp.ip,bp.acor=gd.compute_acor(c_[bp.x,bp.y]); #bp.ne,bp.np=gd.ne,gd.np
bp.dp=gd.dp[bp.ip]; bp.dp0=(bp.dp*bp.acor).sum(axis=1)
if vd.ivcor==1: bp.sigma=vd.sigma[bp.ip]; bp.kbp=vd.kbp[bp.ip]; vd.sigma=None

#check pts inside grid
sindn=nonzero(bp.ie==-1)[0]
if len(sindn)!=0: sys.exit('pts outside of domain: {}'.format(c_[bp.x[sindn],bp.y[sindn]]))
dt00=time.time()-t00; print('finish reading grid info: time={:0.2f}s, myrank={}'.format(dt00,myrank)); sys.stdout.flush()

#read subdomain info
t00=time.time()
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
dt00=time.time()-t00; print('finish reading subdomain info: time={:0.2f}s, myrank={}'.format(dt00,myrank)); sys.stdout.flush()

#-----------------------------------------------------------------------------
#extract data on each processor
#-----------------------------------------------------------------------------
nproc=min(nproc,len(isub))
istacks=[*arange(stacks[0],stacks[1]+1)]

#distribute jobs
for n,isubi in enumerate(isub):
    if n%nproc!=myrank: continue

    #open file
    P=ReadNC('{}/outputs/fabm_state_{:06}.nc'.format(run,isubi),1); sbp=sbps[n]
    ptime=array(P.variables['time'][:])

    #declear variables
    S=npz_data(); S.time=[]
    for i in svars: exec('S.{}=[]'.format(i))

    #read every variables
    for nn,istack in enumerate(istacks): 
        print('reading stack= {} on myrank = {}'.format(istack,myrank))
        C=ReadNC('{}/outputs/schout_{:04}_{}.nc'.format(run,isubi,istack),1)
        ctime=array(C.variables['time'][:]); nt=len(ctime)

        #extract elevation -> compute zcor -> vertical interploate
        eis=[]; k1s=[]; k2s=[]; rats=[]
        for i in arange(nt):
            eii=array(C.variables['elev'][i][sbp.ip]) if ('elev' in C.variables) else 0*sbp.dp
            ei=(eii*sbp.acor).sum(axis=1); eis.append(ei)
            if len(svars)==1 and svars[0]=='elev': continue

            #compute zcor
            zii=[]; kbpii=[]
            for k in arange(3):
                if vd.ivcor==1: ziii=vd.compute_zcor(sbp.dp[:,k],eii[:,k],sigma=sbp.sigma[:,k,:],kbp=sbp.kbp[:,k],method=1)
                if vd.ivcor==2: ziii,kbpiii=vd.compute_zcor(sbp.dp[:,k],eii[:,k],method=1,ifix=1); kbpii.append(kbpiii)
                zii.append(ziii)
            zi=(array(zii)*sbp.acor.T[...,None]).sum(axis=0).T
            if vd.ivcor==2: sbp.kbp=array(kbpii).T.astype('int')

            #station depth
            mzi=sbp.z.copy()
            if ifs==1: mzi=-mzi+ei

            #interpolation in the vertical
            k1=ones(sbp.nsta)*nan; k2=ones(sbp.nsta)*nan; rat=ones(sbp.nsta)*nan
            fp=mzi<=zi[0];  k1[fp]=0; k2[fp]=0; rat[fp]=0   #bottom
            fp=mzi>=zi[-1]; k1[fp]=(vd.nvrt-1); k2[fp]=(vd.nvrt-1); rat[fp]=1  #surface
            for k in arange(vd.nvrt-1):
                fp=(mzi>=zi[k])*(mzi<zi[k+1])
                k1[fp]=k; k2[fp]=k+1
                rat[fp]=(mzi[fp]-zi[k][fp])/(zi[k+1][fp]-zi[k][fp])
            if sum(isnan(r_[k1,k2,rat]))!=0: sys.exit('check vertical interpolation')
            k1s.append(k1); k2s.append(k2); rats.append(rat)
        eis=array(eis); k1s=array(k1s).astype('int'); k2s=array(k2s).astype('int'); rats=array(rats)
        if len(svars)==1 and svars[0]=='elev': Si.elev.extend(array(eis).T);  continue

        #for each variable at each time step
        Si=npz_data()
        for i in svars: exec('Si.{}=[]'.format(i))
        for i in arange(nt):
            it=nonzero(ptime==ctime[i])[0][0]; k1=k1s[i]; k2=k2s[i]; rat=rats[i]
            sindp=arange(sbp.nsta)
            for m,svar in enumerate(svars):
                exec('tri=array(P.variables["{}"][it,sbp.ie])'.format(svar,svar))
                trii=tri[sindp,k1]*(1-rat)+tri[sindp,k2]*rat
                exec('Si.{}.append(trii)'.format(svar))
        S.time.extend(ctime)
        for i in svars: exec('S.{}.extend(Si.{})'.format(i,i))

    #save results
    S.time=array(S.time)/86400
    for i in svars: exec('S.{}=array(S.{}).T'.format(i,i))

#-----------------------------------------------------------------------------
#combine results from all ranks
#-----------------------------------------------------------------------------
if igather==1 and myrank<nproc: save_npz('{}_{}'.format(sname,myrank),S)
comm.Barrier()
if igather==0: sdata=comm.gather(S,root=0)
if igather==1 and myrank==0: sdata=[loadz('{}_{}.npz'.format(sname,i)) for i in arange(nproc)]

if myrank==0: 
   S=npz_data(); S.time=[]; S.bp=bp
   for i in svars: exec('S.{}=[]'.format(i))
   for i in arange(nproc):
       Si=sdata[i]; S.time=Si.time
       for m,svar in enumerate(svars): exec('S.{}.extend(Si.{})'.format(svar,svar))

   #save data        
   for i in svars: exec('S.{}=array(S.{})[sinds]'.format(i,i))
   if fmt==0:
      save_npz('{}'.format(sname),S)
   else:
      #write out ASCII file
      for i in svars: exec('ds=[1,*arange(2,array(S.{}).ndim),0]; S.{}=array(S.{}).transpose(ds)'.format(i,i,i))
      fid=open('{}.dat'.format(sname),'w+')
      for i,ti in enumerate(S.time):
          datai=[]
          for svar in svars: exec('datai.extend(S.{}[{}].ravel())'.format(svar,i))
          fid.write(('{:12.6f}'+' {:10.6f}'*len(datai)+'\n').format(ti,*datai))
      fid.close()
   if igather==1: [os.remove('{}_{}.npz'.format(sname,i)) for i in arange(nproc)] #clean

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora'] else os._exit(0)
