#!/usr/bin/env python3
from pylib import *
#import xarray as xr

#inputs
run='run2f'
stack=[1,73]
icmb_serial=0  # combine when running
svars=['elev','temp','salt','hvel']
#svars=['elev']
stps=['Station.bp_Elev','Station.bp_WQ','Station.bp_Current'] #station.bp files
qnode='haswell' #'skylake' #'haswell'
nproc=5

#------pre-processing-------------------------
bdir=os.path.abspath(os.path.curdir)
if icmb_serial>=1: nproc=1

#-------on front node-------------------------
if os.getenv('param')==None:
   args=sys.argv
   #link results
   if not os.path.exists(run): os.mkdir(run)
   os.system('cd {}; ln -sf ../../{}/outputs/schout_?.nc ./;ln -sf ../../{}/outputs/schout_??.nc ./'.format(run,run,run))
   os.system('cd {}; ln -sf ../../{}/hgrid.gr3 ./'.format(run,run))
   os.system('cd {}; ln -sf ../../{}/vgrid.in ./'.format(run,run))
   os.system('cd {}; ln -sf ../../{}/param.in ./'.format(run,run))
   os.system('cd {}; ln -sf ../../{}/cosine.in ./'.format(run,run))
   #os.system('cd {}; ln -sf ../Station/Station.bp_* ./'.format(run))
   #submit job on node
   param=[bdir,args[0]]
   os.system('qsub {} -v param="{} {}", -N rd_{} -q {} -e Rd_outputs.e -o Rd_outputs.o -l nodes={}:ppn=1 -l walltime=100:00:00'.format(args[0],*param,run,qnode,nproc))
   sys.exit(0)

param=os.getenv('param').split();
param=[int(i) if i.isdigit() else i for i in param]
bdir=param[0]; fname=param[1];

#submit jobs on each core
if os.getenv('job_on_node')==None:
   print("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname))
   os.system("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname))
   sys.exit()

#----get nproc and myrank--------
comm=MPI.COMM_WORLD
nproc=comm.Get_size()
myrank=comm.Get_rank()

#----distribute work--------------------------
stacks=arange(stack[0],stack[1]+1); istack=[];
for i in arange(len(stacks)):
    if i%nproc==myrank:
        istack.append(stacks[i])
istack=array(istack)

os.chdir('{}/{}'.format(bdir,run))
#read grid
gd=read_schism_hgrid('hgrid.gr3'); gd.compute_ctr();

#compute area coordinates
P=npz_data(); 
for m in arange(len(stps)):
    #read bp file
    bp=read_schism_bpfile(stps[m])

    #compute area coordinate
    ie=near_pts(c_[bp.x,bp.y],c_[gd.xctr,gd.yctr],N=100); ip0=gd.elnode[ie]; i34=gd.i34[ie]
    x1=gd.x[ip0][:,0]; x2=gd.x[ip0][:,1]; x3=gd.x[ip0][:,2]; x4=gd.x[ip0][:,3]; x=bp.x
    y1=gd.y[ip0][:,0]; y2=gd.y[ip0][:,1]; y3=gd.y[ip0][:,2]; y4=gd.y[ip0][:,3]; y=bp.y

    A1=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2
    A11=((x2-x)*(y3-y)-(x3-x)*(y2-y))/2
    A12=((x-x1)*(y3-y1)-(x3-x1)*(y-y1))/2
    A13=((x2-x1)*(y-y1)-(x-x1)*(y2-y1))/2
    fA1=(abs(A11)+abs(A12)+abs(A13)-abs(A1))/abs(A1);

    A2=((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
    A21=((x3-x)*(y4-y)-(x4-x)*(y3-y))/2
    A22=((x-x1)*(y4-y1)-(x4-x1)*(y-y1))/2
    A23=((x3-x1)*(y-y1)-(x-x1)*(y3-y1))/2
    fA2=(abs(A21)+abs(A22)+abs(A23)-abs(A2))/abs(A2);

    ip=[];acor=[];
    for i in arange(bp.nsta):
        if abs(fA1[i])<1e-5: #pt in 1st triange
            a1=A11[i]/A1[i]; a2=A12[i]/A1[i]; a3=1-a1-a2;
            if (a1<0)|(a2<0): sys.exit('check pt: {}, {}, {}, {}, {},{}'.format(i,i34[i],ie[i],A1[i],A11[i],A12[i]))
            ip.append(ip0[i,:3]); acor.append(array([a1,a2,a3]))
        else:
            a1=A21[i]/A2[i]; a2=A22[i]/A2[i]; a3=1-a1-a2;
            if (a1<0)|(a2<0)|(i34[i]==3)|(abs(fA2[i])>=1e-5): sys.exit('check pt: {}, {}, {}, {}, {},{}'.format(i,i34[i],ie[i],A1[i],A11[i],A12[i]))
            ip.append(ip0[i,array([0,2,3])]); acor.append(array([a1,a2,a3]))
    ip=array(ip); acor=array(acor)

    #save acor
    bp.ie=ie; bp.ip=ip; bp.acor=acor
    if stps[m]=='Station.bp_Elev': P.SE=bp
    if stps[m]=='Station.bp_WQ': P.SW=bp
    if stps[m]=='Station.bp_Current': P.SH=bp

#C=xr.open_dataset(fname)
#mtime=C.time.values.astype('datetime64[s]');
#zcor=C.zcor.values; nt,np,nz=zcor.shape

#read results
if (icmb_serial==2): istack=[] 
for istacki in istack:
    while ((icmb_serial==1)*(~os.path.exists('schout_{}.nc'.format(istacki+1)))): sleep(10)

    fname='schout_{}.nc'.format(istacki);
    S=npz_data();

    #read nc values
    C=Dataset(fname);
    mtime=array(C.variables['time'][:])/86400; S.mtime=mtime
    nt,tmp,nz=C.variables['zcor'].shape

    #compute acor
    acor_e=tile(P.SE.acor[None,:,:],[nt,1,1])
    acor_w=tile(P.SW.acor[None,:,:,None],[nt,1,1,nz])
    acor_h=tile(P.SH.acor[None,:,:,None],[nt,1,1,nz])

    print('reading {}: zcor'.format(fname)); sys.stdout.flush()
    #read zcor
    zw=[]; zh=[];
    for i in arange(nt):
        #print('zcor: {}'.format(i)); sys.stdout.flush()
        zcor0=array(C.variables['zcor'][i,:,:]); nz=zcor0.shape[-1]
        zwi=zcor0[P.SW.ip,:]; zw.append(zwi)
        zhi=zcor0[P.SH.ip,:]; zh.append(zhi)
    zw=sum(array(zw)*acor_w,axis=2);
    zh=sum(array(zh)*acor_h,axis=2);

    for m in arange(len(svars)):
        svari=svars[m]
        print('reading {}: {}'.format(fname,svari)); sys.stdout.flush()

        if svari=='elev':
            exec("P.vi=C.variables['{}'][:][:,P.SE.ip]".format(svari))
            datai=sum(array(P.vi)*acor_e,axis=2)
        elif (svari=='temp')|(svari=='salt'):
            bp=P.SW; vi=[]
            for i in arange(nt):
                #print('{}: {}'.format(svari,i)); sys.stdout.flush()
                exec("P.vi=C.variables['{}'][i,:,:]".format(svari))
                vi.append(P.vi[bp.ip,:])
            vi=sum(array(vi)*acor_w,axis=2)

            #interpolation
            datai=zeros([nt,bp.nsta])
            for i in arange(nt):
                for j in arange(bp.nsta):
                    zii=zw[i,j]
                    vii=vi[i,j]
                    bzi=zii[-1]-bp.z[j]
                    #limit 
                    fp=(zii<1e6)*(vii<1e6); zii=zii[fp]; vii=vii[fp]
                    bzi=min(max(zii.min(),bzi),zii.max())
                    fd=interpolate.interp1d(zii,vii,'linear')
                    datai[i,j]=fd(bzi)

            #interpolation,method 1
            #bti=arange(nt)*1e6; ti=tile(bti,[nz,1]).T; ds=ti.shape
            #datai=[];
            #for k in arange(P.SW.nsta):
            #    tii=reshape(ti,prod(ds));
            #    zii=reshape(zw[:,k,:],prod(ds))
            #    vii=reshape(vi[:,k,:],prod(ds))
            #    bzi=zw[:,k,-1]-P.SW.z[k]

            #    fp=(~isnan(zii))*(~isnan(vii))
            #    bvi=sp.interpolate.griddata(c_[tii[fp],zii[fp]],vii[fp],c_[bti,bzi],method='linear',fill_value=nan)
            #    bvin=sp.interpolate.griddata(c_[tii[fp],zii[fp]],vii[fp],c_[bti,bzi],'nearest')
            #    sys.exit()
            #    datai.append(bvi)
            #datai=array(datai).T
        elif svari=='hvel':
            bp=P.SH; vi=[]
            for i in arange(nt):
                #print('{}: {}'.format(svari,i)); sys.stdout.flush()
                exec("P.vi=C.variables['{}'][i,:,:,:]".format(svari))
                vi.append(P.vi[bp.ip,:,:])
            vi=sum(array(vi)*tile(acor_h[:,:,:,:,None],[2]),axis=2)

            #interpolation
            datai=zeros([nt,bp.nsta,2])
            for i in arange(nt):
                for j in arange(bp.nsta):
                    for k in arange(2):
                        zii=zh[i,j,:]
                        vii=vi[i,j,:,k]
                        bzi=zii[-1]-bp.z[j]
                        #limit 
                        fp=(zii<1e6)*(vii<1e6); zii=zii[fp]; vii=vii[fp]
                        bzi=min(max(zii.min(),bzi),zii.max())
                        fd=interpolate.interp1d(zii,vii,'linear')
                        datai[i,j,k]=fd(bzi)             

        #save result
        exec('S.{}=datai'.format(svari))
    #save ith stack results
    save_npz('S_{}'.format(istacki),S)
    print('reading {}: completed'.format(fname)); sys.stdout.flush()

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
            for m in arange(len(svars)):
                exec('S.{}=Si.{}'.format(svars[m],svars[m]))
        else:
            exec('S.time=r_[S.time,Si.mtime]');
            for m in arange(len(svars)):
                exec('S.{}=r_[S.{},Si.{}]'.format(svars[m],svars[m],svars[m]))

    #save result
    S.SE=P.SE; S.SW=P.SW; S.SH=P.SH
    S.param=read_schism_param('param.in')
    save_npz('{}.npz'.format(run),S)

    #clean
    os.system("rm S_*.npz Station.bp_* *.in schout_*.nc hgrid.gr3")
