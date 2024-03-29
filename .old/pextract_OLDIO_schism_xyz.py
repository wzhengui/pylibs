#!/usr/bin/env python3
'''
Extract SCHISM variable values at (x,y,z) from station.bp. 
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
run='/sciclone/data10/wangzg/fabm_dev/RUN12'     #run dir containing outputs
stacks=[1,146]                                    #stacks of schout_*.nc 
sname='RUN12/cosine'                             #name for results
svars=['elev','salt','temp','COS_1']        #SCHISM variables to be extracted
rvars=['elev','salt','temp','NO3']          #rename variable names 
bpfile='/sciclone/data10/wangzg/fabm_dev/RUN12/station.bp'  #file name of station.bp
icmb=0                  #icmb=0: work on uncombined; icmb=1: work on combined schout_*.nc
ifs=1                   #ifs=1: depth relative to surface; ifs=0: fixed depth (z coordiante) 
fmt=0                   #fmt=0: output as *.npz format; fmt=1: output as ASCII

#optional
grid='/sciclone/data10/wangzg/fabm_dev/RUN12/grid.npz'  #saved grid info, to speed up; use hgrid.gr3 and vgrid.in if not exist
igather=1          #igather=1: save data on each rank,then combine; igather=0: use MPI  

#resource requst 
walltime='00:10:00'
qnode='x5672'; nnode=2; ppn=8       #hurricane, ppn=8
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
#qnode='vortex'; nnode=2; ppn=12    #vortex, ppn=12
#qnode='femto'; nnode=2; ppn=12     #femto,ppn=32
#qnode='potomac'; nnode=4; ppn=8    #ches, ppn=12
#qnode='james'; nnode=5; ppn=20     #james, ppn=20
#qnode='frontera'; nnode=1; ppn=56  #frontera, ppn=56 (flex,normal)
#qnode='stampede2'; nnode=1; ppn=48 #stampede2, ppn=48 (skx-normal,skx-dev,normal,etc)

#additional information:  frontera,stampede2
qname='flex'                        #partition name
account='TG-OCE140024'              #stampede2: NOAA_CSDL_NWI,TG-OCE140024

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
if myrank==0: t0=time.time()

#-----------------------------------------------------------------------------
#do MPI work on each core
#-----------------------------------------------------------------------------
nproc=max(min(nproc,int(diff(stacks))),1)
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
if icmb==0:
   t00=time.time()
   subs=gd.read_prop('{}/outputs/global_to_local.prop'.format(run)).astype('int')[bp.ie]
   isub=unique(subs); sbps=[]; sindes=[]
   for i, isubi in enumerate(isub):
       sinde=nonzero(subs==isubi)[0] #elem index of stations 
     
       #build the iegl and ipgl 
       T=read_schism_local_to_global('{}/outputs/local_to_global_{}'.format(run,srank(isubi,run)))
       iegl=dict(zip(T.ielg,arange(T.ne))); ipgl=dict(zip(T.iplg,arange(T.np)))
     
       #compute subdomain ie,ip and acor,dp,z,sigma,kbp
       sbp=zdata(); #sbp.ne,sbp.np=T.ne,T.np
       sbp.ie=array([iegl[k] for k in bp.ie[sinde]])
       sbp.ip=array([[ipgl[k] for k in n ] for n in bp.ip[sinde]])
       sbp.acor=bp.acor[sinde]; sbp.dp=bp.dp[sinde]; sbp.z=bp.z[sinde]; sbp.nsta=len(sinde) 
       if vd.ivcor==1: sbp.sigma=bp.sigma[sinde]; sbp.kbp=bp.kbp[sinde]
       sbps.append(sbp); sindes.extend(sinde)
   sinds=argsort(array(sindes)) #indices to sort station order
   dt00=time.time()-t00; print('finish reading subdomain info: time={:0.2f}s, myrank={}'.format(dt00,myrank)); sys.stdout.flush()
else: 
   isub=[None]; sbps=[bp]; sinds=arange(bp.nsta)

#-----------------------------------------------------------------------------
#extract data on each processor
#-----------------------------------------------------------------------------
#distribute jobs
istacks=[i for i in arange(stacks[0],stacks[1]+1) if i%nproc==myrank]

#initilize data capsule
S=zdata(); S.time=[]; #S.bp=bp
for i in svars: exec('S.{}=[]'.format(i)) 

#extract (x,y,z) value for each stack and each subdomain
for n,istack in enumerate(istacks):
    t00=time.time(); Si=zdata()
    for m in svars: exec('Si.{}=[]'.format(m))
    for m,isubi in enumerate(isub):
        #open schout_*.nc
        if icmb==0: fname='{}/outputs/schout_{}_{}.nc'.format(run,srank(isubi,run),istack)
        if icmb==1: fname='{}/outputs/schout_{}.nc'.format(run,istack)
        if (not os.path.exists(fname)) and icmb==0: sys.exit('not exist: {}'.format(fname))
        C=ReadNC(fname,1); sbp=sbps[m]
        
        #read time
        mti=array(C.variables['time'][:])/86400; nt=len(mti); 
        if m==0: S.time.extend(mti)

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

        #compute (x,y,z) for each variables
        Sii=zdata()
        for mm, svar in enumerate(svars):
            exec('Sii.{}=[]'.format(svar))
            ndim=C.variables[svar].ndim; dim=C.variables[svar].shape; dimname=C.variables[svar].dimensions

            data=[]
            for i in arange(nt):
                k1=k1s[i]; k2=k2s[i]; rat=rats[i]

                #get variable values 
                if ('nSCHISM_hgrid_node' in dimname):
                    trii=array(C.variables[svar][i][sbp.ip])
                elif ('nSCHISM_hgrid_face' in dimname): 
                    trii=array(C.variables[svar][i][sbp.ie])
                else:
                    sys.exit('unknown variable format: {},{}'.format(svar,dim))

                #extend values in the bottom: dim[2] is nvrt
                if ('nSCHISM_vgrid_layers' in dimname):
                   sindp=arange(sbp.nsta)
                   if ('nSCHISM_hgrid_node' in dimname):
                      for nn in arange(3):
                          kbp=sbp.kbp[:,nn]; btri=trii[sindp,nn,kbp]
                          for k in arange(vd.nvrt):
                              fp=k<kbp
                              trii[sindp[fp],nn,k]=btri[fp]
                   elif ('nSCHISM_hgrid_face' in dimname):
                        kbe=sbp.kbp.max(axis=1); btri=trii[sindp,kbe]
                        for k in arange(vd.nvrt):
                            fp=k<kbe
                            trii[sindp[fp],k]=btri[fp]
                   else:
                        sys.exit('unknown variable format: {},{}'.format(svar,dim))

                #horizontal interp
                if ('nSCHISM_hgrid_node' in dimname):
                   if ndim==2: tri=(trii*sbp.acor).sum(axis=1)
                   if ndim==3: tri=(trii*sbp.acor[...,None]).sum(axis=1)
                   if ndim==4: tri=(trii*sbp.acor[...,None,None]).sum(axis=1); rat=rat[:,None]
                else:
                   tri=trii

                #vertical interp
                if ('nSCHISM_vgrid_layers' in dimname):
                   datai=(tri[sindp,k1]*(1-rat)+tri[sindp,k2]*rat)
                else:
                   datai=tri
                data.append(datai)

            #save result from each variables
            exec('ds=[1,0,*arange(2,{}-1)]; Sii.{}.extend(array(data).transpose(ds))'.format(ndim,svar))

        #save result form subdomain
        for i in svars: exec('Si.{}.extend(Sii.{})'.format(i,i)) 

    #combine istack results
    for i in svars: exec('ds=[1,0,*arange(2,array(Si.{}).ndim)]; S.{}.extend(array(Si.{})[sinds].transpose(ds))'.format(i,i,i)) 
    dt00=time.time()-t00; print('finish reading stack={}; time={:0.2f}s, myrank={}'.format(istack,dt00,myrank)); sys.stdout.flush()
S.time=array(S.time); ['S.{}=array(S.{}).astype("float32")'.format(i,i) for i in svars]

#-----------------------------------------------------------------------------
#combine results from all ranks
#-----------------------------------------------------------------------------
if igather==1 and myrank<nproc: savez('{}_{}'.format(sname,myrank),S)
comm.Barrier()
if igather==0: sdata=comm.gather(S,root=0)
if igather==1 and myrank==0: sdata=[loadz('{}_{}.npz'.format(sname,i)) for i in arange(nproc)]

if myrank==0:
   S=zdata(); S.time=[]; S.bp=bp
   for i in rvars: exec('S.{}=[]'.format(i))
   for i in arange(nproc):
       Si=sdata[i]; S.time.extend(Si.time)
       for m,[svar,rvar] in enumerate(zip(svars,rvars)): exec('S.{}.extend(Si.{})'.format(rvar,svar))

   #save data        
   S.time=array(S.time); sind=argsort(S.time); S.time=S.time[sind]
   for i in rvars: exec('ds=[1,0,*arange(2,array(S.{}).ndim)]; S.{}=array(S.{})[sind].transpose(ds)'.format(i,i,i)) 
   if fmt==0:
      savez('{}'.format(sname),S)
   else:
      #write out ASCII file
      for i in rvars: exec('ds=[1,*arange(2,array(S.{}).ndim),0]; S.{}=array(S.{}).transpose(ds)'.format(i,i,i)) 
      fid=open('{}.dat'.format(sname),'w+')
      for i,ti in enumerate(S.time):
          datai=[]
          for rvar in rvars: exec('datai.extend(S.{}[{}].ravel())'.format(rvar,i))
          fid.write(('{:12.6f}'+' {:10.6f}'*len(datai)+'\n').format(ti,*datai))
      fid.close()
   if igather==1: [os.remove('{}_{}.npz'.format(sname,i)) for i in arange(nproc)] #clean

#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
sys.exit(0) if qnode in ['bora'] else os._exit(0)
