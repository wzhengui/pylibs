#!/usr/bin/env python3
#python script used to extract time series at Station.bp
from pylib import *

#--------inputs--------------------------------
run='run4ia'
stack=[1,73]
#svars=['elev','temp','salt','hvel']
svars=[['SED3D_1','SED3D_2','SED3D_3']] #schism variable to be read
snames=['SED3D'] #name to be saved
stps=['Station.bp_COS'] #station.bp files
qnode='haswell' #'skylake' 
nproc=16

#-------flags----------------------------------
icmb_serial=0  #0:normal read; 1: read when model running; 2: only combine *npz (skip read)
ifs=0    #0: station depth relative to free surface; 1: fixed station depth

#------pre-processing-------------------------
bdir=os.path.abspath(os.path.curdir)
if icmb_serial>=1: nproc=1

#-------on front node-------------------------
if os.getenv('param')==None:
   args=sys.argv
   #link results
   if not os.path.exists(run): os.mkdir(run)
   os.system('cd {}; ln -sf ../../{}/outputs/schout_?.nc ./;ln -sf ../../{}/outputs/schout_??.nc ./'.format(run,run,run))
   os.system('cd {}; ln -sf ../Station/Station.bp_* ./'.format(run))
   os.system('cd {}; ln -sf ../../{}/hgrid.gr3 ./'.format(run,run))
   os.system('cd {}; ln -sf ../../{}/vgrid.in ./'.format(run,run))
   os.system('cd {}; ln -sf ../../{}/param.in ./'.format(run,run))
   os.system('cd {}; ln -sf ../../{}/cosine.in ./'.format(run,run))
   #submit job on node
   param=[bdir,args[0]]
   os.system('qsub {} -v param="{} {}", -N rd_{} -q {} -e Rd_outputs.e -o Rd_outputs.o -l nodes={}:ppn=1 -l walltime=100:00:00'.format(args[0],*param,run,qnode,nproc))
   sys.exit()

param=os.getenv('param').split();
param=[int(i) if i.isdigit() else i for i in param]
bdir=param[0]; fname0=param[1];

#submit jobs on each core
if os.getenv('job_on_node')==None:
   print("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname0))
   os.system("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname0))
   sys.exit()

#start working on each core
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
if (icmb_serial==2): istack=[] 

#compute area coordinates for stations
gd=read_schism_hgrid('hgrid.gr3'); gd.compute_ctr();
P=npz_data(); P.nstp=len(stps); P.SE=None; P.SW=None; P.SH=None
for m in arange(P.nstp):
    stp=stps[m]
    #read station.bp
    bp=read_schism_bpfile(stps[m]); bp.compute_acor(gd);

    #save information
    if stp=='Station.bp_Elev': 
       P.SE=bp
    elif (stp=='Station.bp_COS')|(stp=='Station.bp_WQ'):
       P.SW=bp
    elif stp=='Station.bp_Current': 
       P.SH=bp

#read results
for istacki in istack:
    while ((icmb_serial==1)*(~os.path.exists('schout_{}.nc'.format(istacki+1)))): sleep(10)

    fname='schout_{}.nc'.format(istacki);
    S=npz_data();

    #read nc values
    C=ReadNC(fname);
    mtime=array(C.variables['time'][:])/86400; S.mtime=mtime
    nt,tmp,nz=C.variables['zcor'].shape

    #compute ratio matrix for vertical interpolation
    for stp in ['SW','SH']: 
        exec('P.bp=P.{}'.format(stp))
        if P.bp==None: continue

        acor=tile(P.bp.acor[:,:,None],[1,1,nz]); #(nsta,3,nz)
        srat=[]
        for i in arange(nt):
           #compute veritcal coordinates for stations
           zcor0=array(C.variables['zcor'][i,:,:]);
           zcor=sum(zcor0[P.bp.ip,:]*acor,axis=1).T; #(nvrt,nsta) 
           #treat invalid depths 
           fp=abs(zcor)>1e10; zcor[fp]=-zcor[fp]
              
           zs=-P.bp.z #bp.z is positive 
           if ifs==0: zs=zs+zcor[-1]

           #compute index for vertical interpolation
           rat=zeros([P.bp.nsta,nz]);
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
        if sum(srat.sum(axis=2)!=1)!=0: sys.exit('wrong for srat: {}'.format(stp)) 
        if stp=='SW': wrat=srat 
        if stp=='SH': hrat=srat 
    #print('finsih reading {}: zcor'.format(fname)); sys.stdout.flush()

    #read variables
    for m in arange(len(svars)):
        svari=svars[m]; sname=snames[m]
        print('reading {}: {}'.format(fname,svari)); sys.stdout.flush()

        if type(svari)==list:
           #determine variable name
           # for i in arange(len(svari[0])):
           #     if svari[0][i]!=svari[1][i]:
           #        sind=i; break
           # sname=svari[0][:sind]
           # if sname.endswith('_'): sname=sname[:-1]

            #read variable
            for n in arange(len(svari)):
                vi=[]
                for i in arange(nt):
                    exec("P.vi=C.variables['{}'][i,:,:]".format(svari[n]))
                    vi.append(P.vi[P.SW.ip,:])
                vi=sum(array(vi)*tile(P.SW.acor[None,:,:,None],[nt,1,1,nz]),axis=2)
                if n==0:
                   datai=(vi*wrat).sum(axis=2)
                else:
                   datai=datai+(vi*wrat).sum(axis=2)
            
        elif svari=='elev':
            exec("P.vi=C.variables['{}'][:][:,P.SE.ip]".format(svari))
            datai=sum(array(P.vi)*tile(P.SE.acor[None,:,:],[nt,1,1]),axis=2)
        elif (svari in ['temp','salt'])|(svari.startswith('COS'))|(svari.startswith('SED3D')):
            vi=[]
            for i in arange(nt):
                exec("P.vi=C.variables['{}'][i,:,:]".format(svari))
                vi.append(P.vi[P.SW.ip,:])
            vi=sum(array(vi)*tile(P.SW.acor[None,:,:,None],[nt,1,1,nz]),axis=2)
            datai=(vi*wrat).sum(axis=2)
        elif svari=='hvel':
            vi=[]
            for i in arange(nt):
                exec("P.vi=C.variables['{}'][i,:,:,:]".format(svari)) #(np,nvrt,2)
                vi.append(P.vi[P.SH.ip,:,:])
            vi=sum(array(vi)*tile(P.SH.acor[None,:,:,None,None],[nt,1,1,nz,2]),axis=2)
            datai=(vi*tile(hrat[:,:,:,None],[1,1,1,2])).sum(axis=2)

        #save result
        exec('S.{}=datai'.format(sname))
    #save ith stack results
    save_npz('S_{}'.format(istacki),S)
    #print('finsih reading {}: variables'.format(fname)); sys.stdout.flush()

#collect results
comm.Barrier()
if myrank==0:
    #wait all results
    while(True):
        iflag=len(stacks)
        for i in arange(len(stacks)):
            if os.path.exists('S_{}.npz'.format(stacks[i])): iflag=iflag-1
        if iflag==0: break
        if iflag!=0: time.sleep(1)

    #read result
    S=npz_data();
    for i in arange(len(stacks)):
        Si=loadz('S_{}.npz'.format(stacks[i]))

        if i==0:
            exec('S.time=Si.mtime');
            for m in arange(len(snames)):
                exec('S.{}=Si.{}'.format(snames[m],snames[m]))
        else:
            exec('S.time=r_[S.time,Si.mtime]');
            for m in arange(len(snames)):
                exec('S.{}=r_[S.{},Si.{}]'.format(snames[m],snames[m],snames[m]))

    #save result
    S.SE=P.SE; S.SW=P.SW; S.SH=P.SH
    S.param=read_schism_param('param.in')
    save_npz('{}.npz'.format(run),S)

    #clean
    os.system("rm S_*.npz Station.bp_* *.in schout_*.nc hgrid.gr3")
    dt=time.time()-t0
    print('finish reading {}: {}s'.format(run,dt)); sys.stdout.flush()
    
