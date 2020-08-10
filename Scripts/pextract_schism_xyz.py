#!/usr/bin/env python3
#python script used to extract slabs of schism results
from pylib import *
#import xarray as xr

#---------inputs------------------------------
run='run4ia'
stack=[2,3]
svars=['salt',['SED3D_1','SED3D_2','SED3D_3']]
snames=['salt','SED3D']
depths=[0,1]
qnode='haswell'
nproc=2

#-------flags---------------------------------
icmb=0  #0:normal read (parallel); 1: read when model running (serial);
        #2:read but not combine (parallel); 3: not read, only combine *npz (skip read)
ifs=0   #0: station depth relative to free surface; 1: fixed station depth

#------pre-processing-------------------------
bdir=os.path.abspath(os.path.curdir)
if (icmb==1)|(icmb==3): nproc=1

#-------on front node-------------------------
if os.getenv('param')==None:
   args=sys.argv
   #link results
   if not os.path.exists(run): os.mkdir(run)
   os.system('cd {}; ln -sf ../../{}/outputs/schout_?.nc ./;ln -sf ../../{}/outputs/schout_??.nc ./'.format(run,run,run))
   #submit job on node
   param=[bdir,args[0]]
   os.system('qsub {} -v param="{} {}", -N rd_{} -q {} -e Rd_outputs.e -o Rd_outputs.o -l nodes={}:ppn=1 -l walltime=100:00:00'.format(args[0],*param,run,qnode,nproc))
   #os.system('qsub {} -v param="{} {}", -N rd_{} -q {} -e Rd_outputs.e -o Rd_outputs.o -l procs={} -l walltime=100:00:00'.format(args[0],*param,run,qnode,nproc))
   sys.exit(0)

param=os.getenv('param').split();
param=[int(i) if i.isdigit() else i for i in param]
bdir=param[0]; fname0=param[1];

#submit jobs on each core
if os.getenv('job_on_node')==None:
   print("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname0))
   os.system("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname0))
   sys.exit()

#start to work on each core
os.chdir('{}/{}'.format(bdir,run))

#----get nproc and myrank--------
comm=MPI.COMM_WORLD
nproc=comm.Get_size()
myrank=comm.Get_rank()

t0=time.time()
#----distribute work--------------------------
stacks=arange(stack[0],stack[1]+1); istack=[];
for i in arange(len(stacks)):
    if i%nproc==myrank:
        istack.append(stacks[i])
istack=array(istack)
if(icmb==3): istack=[]

#read results
P=npz_data();
for istacki in istack:
    while ((icmb==1)*(os.path.exists('schout_{}.nc'.format(istacki+1)))): sleep(10)
    fname='schout_{}.nc'.format(istacki);
    S=npz_data();

    #read nc values
    C=Dataset(fname);
    mtime=array(C.variables['time'][:])/86400; S.time=mtime.astype('float32'); S.depth=array(depths).astype('float32')
    nt,np,nz=C.variables['zcor'].shape
    [exec('S.{}=[]'.format(sname)) for sname in snames]

    for i in arange(nt): 
        #compute matrix for vertical interpolation
        zcor=array(C.variables['zcor'][i,:,:]).T; #(nz,np)
        #treat invalid depths
        fp=abs(zcor)>1e10; zcor[fp]=-zcor[fp]
        srat=[]
        for m in arange(len(depths)):
            zs=-array(depths[m]).astype('float')
            if ifs==0: zs=zs+zcor[-1] 
            rat=zeros([np,nz])
            fp=zs<=zcor[0]; rat[fp,0]=1
            fp=zs>zcor[-1]; rat[fp,-1]=1
            for k in arange(1,nz):
                zi0=zcor[k-1]; zi=zcor[k]
                fp=(zs>zi0)*(zs<=zi)
                rati=(zs[fp]-zi0[fp])/(zi[fp]-zi0[fp])
                rat[fp,k]=rati
                rat[fp,k-1]=1-rati
            srat.append(rat)
        srat=array(srat)
        if sum(abs(srat.sum(axis=2)-1.0)>1e-6)!=0: sys.exit('wrong for srat: {},step={}'.format(fname,i))

        #read slices
        for n in arange(len(svars)):
            svari=svars[n]; sname=snames[n] 
            if i==0: print('reading {} (slab): {} '.format(fname,svari)); sys.stdout.flush()
            datai=[]
            for m in arange(len(depths)):
               if (svari=='elev')*(m!=0): continue
               #read raw data
               if type(svari)==list:
                  for k in arange(len(svari)):  
                      if k==0: 
                         exec("P.vi=C.variables['{}'][i]".format(svari[k]))
                      else:
                         exec("P.vi=P.vi+C.variables['{}'][i]".format(svari[k]))
               else:
                  exec("P.vi=C.variables['{}'][i]".format(svari))
               #extract
               if svari=='elev': 
                  dataii=P.vi 
               elif svari=='hvel':
                  dataii=(P.vi*tile(srat[m][:,:,None],[1,1,2])).sum(axis=1)
               else: 
                  dataii=(P.vi*srat[m]).sum(axis=1)
               datai.append(dataii)
            datai=array(datai)
            exec('S.{}.append(datai)'.format(sname)) 

    #save data
    for n in arange(len(svars)):
        svari=svars[n]; sname=snames[n] 
        if svari=='elev':
           exec("S.{}=squeeze(array(S.{})).astype('float32')".format(sname,sname))
        elif svari=='hvel': 
           exec("S.{}=array(S.{}).transpose([0,2,3,1]).astype('float32')".format(sname,sname))
        else:
           exec("S.{}=array(S.{}).transpose([0,2,1]).astype('float32')".format(sname,sname))

    #save data
    save_npz('{}_slab_{}'.format(run,istacki),S)

#colloect results
comm.Barrier()
if myrank==0:
    #wait all results
    while(True):
        iflag=len(stacks)
        for i in arange(len(stacks)):
            if os.path.exists('{}_slab_{}.npz'.format(run,stacks[i])): iflag=iflag-1
        if iflag==0: break
        if iflag!=0: time.sleep(1)

    #read result
    if icmb!=2:
       S=npz_data();
       for i in arange(len(stacks)):
           Si=loadz('{}_slab_{}.npz'.format(run,stacks[i]))

           if i==0:
              exec('S.time=Si.time; S.depth=Si.depth');
              for m in arange(len(snames)):
                  exec('S.{}=Si.{}'.format(snames[m],snames[m]))
           else:
              exec('S.time=r_[S.time,Si.time]');
              for m in arange(len(snames)):
                  exec('S.{}=r_[S.{},Si.{}]'.format(snames[m],snames[m],snames[m]))

       #save result
       save_npz('{}_slab.npz'.format(run),S)
       [os.system("rm {}_slab_{}.npz".format(run,i)) for i in stacks]

    #clean
    os.system("rm schout_*.nc ")
    dt=time.time()-t0
    print('finish reading {}: {}s'.format(run,dt)); sys.stdout.flush()
