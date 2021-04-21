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
ibatch=0           #ibatch=0: submit batch job;   ibatch=1: run script locally (interactive)

#resource requst 
walltime='1:00:00'
#qnode='bora'; nnode=1; ppn=5      #bora, ppn=20
#qnode='vortex'; nnode=10; ppn=12   #vortex, ppn=12
qnode='x5672'; nnode=2; ppn=8      #hurricane, ppn=8
#qnode='potomac'; nnode=1; ppn=2    #ches, ppn=12
#qnode='james'; nnode=1; ppn=2     #james, ppn=20
#qnode='femto'; nnode=1; ppn=2      #femto,ppn=32, not working yet
#qnode='skylake'; nnode=2; ppn=36    #viz3,skylake, ppn=36
#qnode='haswell'; nnode=2; ppn=2   #viz3,haswell, ppn=24,or 28

#-----------------------------------------------------------------------------
#pre-processing
#-----------------------------------------------------------------------------
nproc=nnode*ppn
bdir=os.path.abspath(os.path.curdir)
jname='Rd_{}'.format(os.path.basename(run)) #job name
scrout='screen.out'

if ibatch==1: os.environ['job_on_node']='1'; os.environ['bdir']=bdir #run local
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
       rcode="mpiexec -x job_on_node=1 -x bdir='{}' -n {} {} >& {}".format(bdir,nproc,bcode,scrout)
    elif qnode=='femto':
       pypath='/sciclone/home10/wangzg/bin/pylibs/Scripts/:/sciclone/home10/wangzg/bin/pylibs/Utility/'
       rcode="srun --export=job_on_node=1,bdir='{}',PYTHONPATH='{}' {} >& {}".format(bdir,pypath,bcode,scrout)
    elif qnode=='x5672' or qnode=='vortex' or qnode=='potomac' or qnode=='james':
       rcode="mvp2run -v -e job_on_node=1 -e bdir='{}' {} >& {}".format(bdir,bcode,scrout)
    elif qnode=='skylake' or qnode=='haswell':
       rcode="mpiexec --env job_on_node 1 --env bdir='{}' -np {} {} >& {}".format(bdir,nproc,bcode,scrout)
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
           exec('S.{}=array(S.{})[sind].transpose([1,0,*arange(2,{})])'.format(svar,svar,Si.ndim[m]))
       S.levels=array(levels)
       if fmt==1: save_npz('{}_{}'.format(sname,istack),S) 

       #combine all stacks
       if fmt==0: 
          C.time.extend(S.time)
          for m in svars: exec('C.{}.extend(S.{})'.format(m,m))

   if fmt==0:
      for m in svars: exec('C.{}=array(C.{})'.format(m,m))
      C.time=array(C.time); C.levles=array(levels); save_npz(sname,C)

   #clean
   os.system('rm {}_*_*.npz'.format(sname))
#-----------------------------------------------------------------------------
#finish MPI jobs
#-----------------------------------------------------------------------------
comm.Barrier()
if myrank==0: dt=time.time()-t0; print('total time used: {} s'.format(dt)); sys.stdout.flush()
if qnode=='x5672' or qnode=='james':
   os._exit(0)
else:
   sys.exit(0)
