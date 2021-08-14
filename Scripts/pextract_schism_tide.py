#!/usr/bin/env python3
#perform harmonic analysis on schism results
from pylib import *
close("all")

#------------------------------------------------------------------------------
#input
#------------------------------------------------------------------------------
run='run3'
stack=[1,12]
tide_name=['O1','K1','Q1','P1','M2','S2','K2','N2']
svar='elev'
qnode='haswell'
nproc=128
sdir='../{}'.format(run) #model run dir
tdir='/home/zwang/bin/Harmonic_Analysis' #where tidal_analyze exists
sname='{}_tide'.format(run)
include=None #'include.gr3'  #include=['include.gr3'/None]: only extract certian pts in "sdir/include.gr3"

#combine extracted results only
icmb_serial=0 #0:normal read; 1: only combine *npz (skip read)

#------------------------------------------------------------------------------
#-------on front node-------------------------
#------------------------------------------------------------------------------
#pre-processing
bdir=os.path.abspath(os.path.curdir)
sdir=os.path.abspath(sdir)
if icmb_serial>=1: nproc=1

if os.getenv('param')==None:
   args=sys.argv
   if not os.path.exists(run): os.mkdir(run)

   #write tidal_const.dat
   T=loadz('{}/tide_fac_const.npz'.format(tdir))
   fid=open('{}/tidal_const.dat'.format(run),'w+')
   fid.write('{}\n'.format(len(tide_name)))
   for m in arange(len(tide_name)):
       fp=T.name==tide_name[m]; freqi=squeeze(T.freq[fp])
       fid.write('{}\n {:e}\n'.format(tide_name[m],freqi))
   fid.close()

   #link results
   os.system('cd {}; ln -sf {}/outputs/schout_?.nc ./;ln -sf {}/outputs/schout_??.nc ./'.format(run,sdir,sdir))
   os.system('cd {}; ln -sf {}/hgrid.* ./'.format(run,sdir))
   if include is not None: os.system('cd {}; ln -sf {}/{} ./'.format(run,sdir,include))

   #submit job on node
   param=[bdir,args[0]]
   #os.system('qsub {} -v param="{} {}", -N HA_{} -q {} -e HA_outputs.e -o HA_outputs.o -l nodes={}:ppn=1 -l walltime=100:00:00'.format(args[0],*param,run,qnode,nproc))
   os.system('qsub {} -v param="{} {}", -N HA_{} -q {} -e HA_outputs.e -o HA_outputs.o -l procs={} -l walltime=100:00:00'.format(args[0],*param,run,qnode,nproc))
   sys.exit()

#------------------------------------------------------------------------------
#-------still on front node,but in batch mode---------
#------------------------------------------------------------------------------
param=os.getenv('param').split();
param=[int(i) if i.isdigit() else i for i in param]
bdir=param[0]; fname0=param[1];

#submit jobs on each core
if os.getenv('job_on_node')==None:
   print("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname0))
   os.system("cd {}; mpiexec -np {} --env job_on_node 1 {}>>screen.out".format(bdir,nproc,fname0))
   sys.exit()

#------------------------------------------------------------------------------
#-------still on computational node---------
#------------------------------------------------------------------------------
#start working on each core
os.chdir('{}/{}'.format(bdir,run))

#----get nproc and myrank--------
comm=MPI.COMM_WORLD
nproc=comm.Get_size()
myrank=comm.Get_rank()

t0=time.time()
#----distribute work--------------------------
if icmb_serial==0:
    #get all pts
    if include is None:
       if os.path.exists('hgrid.npz'): 
          gd=loadz('hgrid.npz').hgrid
       else:
          gd=read_schism_hgrid('hgrid.gr3')
    else:
       gd=read_schism_hgrid(include)
    ipts=nonzero(gd.dp!=0)[0]; npts=len(ipts)

    #distribute work
    for i in arange(nproc):
        if i!=myrank: continue
        if i!=(nproc-1):
           npt=int(npts/nproc); id1=npt*i; id2=npt*(i+1)
           ipt=ipts[id1:id2]
        else:
           npt=npts-int(npts/nproc)*(nproc-1); id1=npts-npt
           ipt=ipts[id1:]

    #read schism results
    ti=[]; yi=[]
    for istack in arange(stack[0],stack[1]+1):
        S=ReadNC('schout_{}.nc'.format(istack))
        tii=array(S.variables['time'][:])
        if svar=='elev':
           yii=array(S.variables[svar][:,ipt])
        else:
           sys.exit('coding according for {}'.format(svar))
        ti.extend(tii); yi.extend(yii)
        S.close()
    ti=array(ti); yi=array(yi)
    print('finish reading schout_[{}-{}].nc: myrank={}'.format(*stack,myrank)); sys.stdout.flush()

    #perform HA for each ipt
    HA=[]
    for i in arange(npt):
        fnamei='schism_HA_timeseries_{}_{}.dat'.format(myrank,i)
        snamei='schism_HA_components_{}_{}.dat'.format(myrank,i)
        fid=open(fnamei,'w+'); fid.writelines(['{} {}\n'.format(ii,kk) for ii,kk in zip(ti,yi[:,i])]); fid.close()
        os.system('{}/tidal_analyze {} tidal_const.dat {}'.format(tdir,fnamei,snamei))
        lines=array([line.strip().split() for line in open(snamei,'r').readlines()])
        tname=lines[:,0]; HA.extend(lines[:,1:][None,:].astype('float'))
        os.remove(fnamei); os.remove(snamei)
    tname=array(tname); HA=array(HA)

    #save for each proc
    S=zdata();
    S.inode=ipt; S.tidal_name=tname; S.amplitude=HA[:,:,0]; S.phase=HA[:,:,1]; S.note='e.g. gd.dp[inode]=S.amplitude[:,2]'
    savez('S_{}'.format(myrank),S)
    print('finish HA analysis: myrank={}'.format(myrank)); sys.stdout.flush()

#collect results
if myrank==0:
    if icmb_serial==0:
        #wait all results
        while(True):
            iflag=nproc
            for i in arange(nproc):
                if os.path.exists('S_{}.npz'.format(i)): iflag=iflag-1
            if iflag==0: break
            if iflag!=0: time.sleep(5)
    else:
        fnames=[i for i in os.listdir() if (i.startswith('S_')*i.endswith('.npz'))]
        nproc=len(fnames)

    #read and combine result
    S=zdata(); S.inode=[]; S.amplitude=[]; S.phase=[]
    for i in arange(nproc):
        if icmb_serial==0:
           Si=loadz('S_{}.npz'.format(i))
        elif icmb_serial==1:
           Si=loadz(fnames[i])
        else:
           sys.exit('fname wrong')

        S.inode.extend(Si.inode)
        S.amplitude.extend(Si.amplitude)
        S.phase.extend(Si.phase)
        S.tidal_name=Si.tidal_name; S.note=Si.note
    S.inode=array(S.inode); S.amplitude=array(S.amplitude); S.phase=array(S.phase)

    #save result
    savez('{}'.format(sname),S)

    #clean
    os.system("rm S_*.npz schout_?.nc schout_??.nc hgrid.* tidal_const.dat")
    if include is not None: os.system("rm {}".format(include))
    dt=time.time()-t0
    print('finish reading {}: {}s'.format(run,dt)); sys.stdout.flush()

#comm.Ibarrier()
sys.exit('done: myrank={}'.format(myrank))
#comm.Barrier()
