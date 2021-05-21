#!/usr/bin/env python3
'''
Extract slabs of SCHISM variables 
  1). work for both uncombined and combined SCHISM outputs
  2). can extract multiple variables at the same time (only node-based variables)
  3). can work in interactive or batch mode 
  4). native schism level (>=1), or zcoor (<=0)
'''
from pylib import *
import time

#-----------------------------------------------------------------------------
#Input
#-----------------------------------------------------------------------------
run='/sciclone/data10/wangzg/DSP/RUN01c'     #run dir containing outputs
sname='slab'      #name for the resutls 
stacks=[7,9]        #stacks of schout_*.nc (1st and last)
svars=['salt','elev']      #variable to be extracted
levels=[1,10,100]        #schism level indices (1-nvrt: surface-bottom; (>nvrt): kbp level), or fix depthis ( z coordinate, values<=0)
nspool=12            #sub-sampling frequency within each stack (1 means all)
icmb=0              #icmb=0: work on uncombined; icmb=1: work on combined schout_*.nc
fmt=0               #fmt=0: one output file;   fmt=1: one output file for each stack

#optional
grid='./grid.npz'  #saved grid info, to speed up; use hgrid.gr3 and vgrid.in if not exist

#resource requst 
walltime='00:10:00'
qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
#qnode='bora'; nnode=2; ppn=20      #bora, ppn=20
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
dt=time.time()-t00; print('finish reading grid info: time={:0.2f}s, myrank={}'.format(dt,myrank)); sys.stdout.flush()

#read subdomain info
if icmb==0:
   t00=time.time(); nsub=getglob('{}/outputs/local_to_global_0000'.format(run)).nproc; isub=[]; subs=[]; iplg=[]
   for i in arange(nsub):
       if i%nproc!=myrank: continue  
       #build the ielg and ipgl tables and 
       S=read_schism_local_to_global('{}/outputs/local_to_global_{:04}'.format(run,i)); S.dp=gd.dp[S.iplg]
       if vd.ivcor==1: S.kbp=vd.kbp[S.iplg]; S.sigma=vd.sigma[S.iplg]
       if vd.ivcor==2: S.kbp=zeros(S.np).astype('int') #pure sigma
       isub.append(i); subs.append(S); iplg.extend(S.iplg)
   dt=time.time()-t00; print('finish reading subdomain info: time={:0.2f}s, myrank={}'.format(dt,myrank)); sys.stdout.flush()
else:
   S=npz_data(); S.np=gd.np; S.iplg=arange(gd.np).astype('int'); S.dp=gd.dp
   if vd.ivcor==1: S.kbp=vd.kbp; S.sigma=vd.sigma
   if vd.ivcor==2: S.kbp=zeros(gd.np).astype('int') #pure sigma
   isub=[0]; subs=[S]; iplg=S.iplg 
if vd.ivcor==1: vd.sigma=None
gd=None

#check level or z-coor
levels=array(levels)
if levels.min()>=1: imed=0 
if levels.max()<=0: imed=1
if (sum(levels>0)*sum(levels<0))!=0: sys.exit('check level type: {}'.format(levels))

#-----------------------------------------------------------------------------
#extract data on each processor
#-----------------------------------------------------------------------------
#distribute jobs
if icmb==0: istacks=[*arange(stacks[0],stacks[1]+1)]
if icmb==1: istacks=[i for i in arange(stacks[0],stacks[1]+1) if i%nproc==myrank]

#extract slab value for each stack and each subdomain
for n,istack in enumerate(istacks):
    t00=time.time(); S=npz_data(); S.iplg=array(iplg).astype('int'); S.ndim=[]
    for m,isubi in enumerate(isub):
        #open schout_*.nc
        if icmb==0: fname='{}/outputs/schout_{:04}_{}.nc'.format(run,isubi,istack)
        if icmb==1: fname='{}/outputs/schout_{}.nc'.format(run,istack)
        if (not os.path.exists(fname)) and icmb==0: sys.exit('not exist: {}'.format(fname))
        C=ReadNC(fname,1); sub=subs[m]

        #read time
        mti=array(C.variables['time'][:])/86400; nt=len(mti); 
        if m==0: S.time=mti[::nspool]

        #compute zcor, and do interpolation
        k1s=[]; rats=[]
        for i in arange(nt):
            if (i%nspool!=0) or imed==0: continue
            ei=array(C.variables['elev'][i]); k1si=[]; ratsi=[]
            if vd.ivcor==1: zi=vd.compute_zcor(sub.dp,ei,sigma=sub.sigma,kbp=sub.kbp,method=1)
            if vd.ivcor==2: zi=vd.compute_zcor(sub.dp,ei,method=1,ifix=1)

            #for interpolation
            for nn, level in enumerate(levels): 
                k1=ones(sub.np)*nan; rat=ones(sub.np)*nan; mzi=ones(sub.np)*level+ei 
                fpz=mzi>=zi[:,-1]; k1[fpz]=vd.nvrt-2;  rat[fpz]=1.0 #surface
                fpz=mzi<zi[:,0]; k1[fpz]=sub.kbp[fpz]; rat[fpz]=0.0 #bottom
                for k in arange(vd.nvrt-1):
                    fpz=(mzi>=zi[:,k])*(mzi<zi[:,k+1])
                    k1[fpz]=k; rat[fpz]=(mzi[fpz]-zi[fpz,k])/(zi[fpz,k+1]-zi[fpz,k]) 
                if sum(isnan(r_[k1,rat]))!=0: sys.exit('check vertical interpolation')
                k1si.append(k1); ratsi.append(rat)
            k1s.append(k1si); rats.append(ratsi)
        k1s=array(k1s).astype('int'); rats=array(rats)

        #compute value for each variables
        for mm, svar in enumerate(svars):
            dimname=C.variables[svar].dimensions; ivs=C.variables[svar].ivs
            if m==0: exec('S.{}=[]'.format(svar)); S.ndim.append(len(dimname))

            #extract slabs
            data=[]; irec=0
            for i in arange(nt):
                if i%nspool!=0: continue
                if 'nSCHISM_vgrid_layers' in dimname:  #3D
                   datai=[]
                   if imed==0: #for schism levels
                      for k, level in enumerate(levels):
                          if level>vd.nvrt:  dataii=array(C.variables[svar][i][arange(sub.np),sub.kbp])
                          if level<=vd.nvrt: dataii=array(C.variables[svar][i,:,vd.nvrt-int(level)])
                          datai.append(dataii)
                   else:
                      for k, level in enumerate(levels):
                          k1=k1s[irec,k]; k2=k1+1; rat=rats[irec,k]
                          dataii=array(C.variables[svar][i][arange(sub.np),k1]*(1-rat)+C.variables[svar][i][arange(sub.np),k2]*rat)
                          datai.append(dataii)
                   if ivs==1: data.append(array(datai).T)
                   if ivs==2: data.append(array(datai).transpose([1,0,2]))
                   irec=irec+1   
                else:  #2D 
                   data.append(array(C.variables[svar][i]))

            #remove missing data 
            data=array(data).transpose([1,0,*arange(2,len(dimname))]); fpn=data>1.e9; data[fpn]=nan
            exec('S.{}.extend(data)'.format(svar))

    save_npz('{}_{}_{}'.format(sname,istack,myrank),S)
    dt=time.time()-t00; print('finish reading stack={} on myrank={} time={:0.2f}s'.format(istack,myrank,dt)); sys.stdout.flush()
#-----------------------------------------------------------------------------
#combine results from all ranks
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: 
   if fmt==0:  C=npz_data(); C.time=[]; [exec('C.{}=[]'.format(m)) for m in svars]

   #combine result from subdomains
   for istack in arange(stacks[0],stacks[1]+1):
       S=npz_data(); [exec('S.{}=[]'.format(m)) for m in svars]; iplg=[]
       for n in arange(nproc):
           fname='{}_{}_{}.npz'.format(sname,istack,n)
           if not os.path.exists(fname): continue
           Si=loadz(fname); iplg.extend(Si.iplg); S.time=Si.time
           for m in svars: exec('S.{}.extend(Si.{})'.format(m,m))
       
       #sort
       tmp,sind=unique(array(iplg),return_index=True)
       for m,svar in enumerate(svars): 
           if Si.ndim[m]==2: exec('S.{}=array(S.{})[sind].transpose([1,0])'.format(svar,svar))
           if Si.ndim[m]==3: exec('S.{}=array(S.{})[sind].transpose([1,2,0])'.format(svar,svar))
           if Si.ndim[m]==4: exec('S.{}=array(S.{})[sind].transpose([1,2,0,3])'.format(svar,svar))
       S.levels=array(levels)
       if fmt==1: save_npz('{}_{}'.format(sname,istack),S) 

       #combine all stacks
       if fmt==0: 
          C.time.extend(S.time)
          for m in svars: exec('C.{}.extend(S.{})'.format(m,m))

   if fmt==0:
      for m in svars: exec('C.{}=array(C.{})'.format(m,m))
      C.time=array(C.time); C.levels=array(levels); save_npz(sname,C)

   #clean
   os.system('rm {}_*_*.npz'.format(sname))
#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
os._exit(0)
