#!/usr/bin/env python3
from pylib import *
from PIL import Image
import multiprocessing as mp
from multiprocessing.shared_memory import SharedMemory as SM
mpl.use('agg'); close("all")

#--------------------------------------------------
#Inputs
#--------------------------------------------------
sdir='/sciclone/home/yinglong/DISKS/schism10/tsunami/Project_2024_NorthCoast_PTHA/RUN02c'  #run dir
xm=[-124.538,-123.283]; ym=[45.416,46.726]  #subset region
zlim=[-10,10]                               #elevation limit
stacks=[1,24]                               #output stacks to be plotted
nspool=5                                    #sub-sampling frequency within each stack (1 means all)
bpfile='station.bp'                         #station for plotitng time series; commment out if not need
nproc=32                                    #number of process
delay=200                                   #animation step (milli-seconds)
fname='animation_tsunami'                   #animation file name

#--------------------------------------------------
#read grid and output information
#--------------------------------------------------
istacks0=arange(stacks[0],stacks[1]+1); nstack=len(istacks0)
p=read(sdir+'/param.nml',3); dt=int(p.dt); nrec0=int(p.ihfskip/p.nspool); nrec=int(nrec0/nspool)
istacks=tile(istacks0,[nrec,1]).T.ravel(); irecs=tile(arange(nrec)*nspool,[nstack,1]).ravel()
its=arange(nstack*nrec); nt=len(its); mt=arange(nt)*dt*nspool

if not fexist('hgrid.npz'): read(sdir+'/hgrid.gr3').save('hgrid.npz')
gd=read('hgrid.npz')

#--------------------------------------------------
#fread: extract time series @station; fplot: plot
#--------------------------------------------------
def fread(it):
    irec=irecs[it]; istack=istacks[it]; fid=read('{}/outputs/out2d_{}.nc'.format(sdir,istack),1)
    mys[:,it]=sum(fid.elevation[irec][pip]*pacor,axis=1); fid.close()
    if it%10==0: print('done ts: {}/{}'.format(it,nt))

def fplot(its):
    for it in its:
        irec=irecs[it]; istack=istacks[it]; 
        fid=read('{}/outputs/out2d_{}.nc'.format(sdir,istack),1); elev=fid.elevation[irec]; fid.close()
        if npt==0:
           figure(figsize=[10,10])
           gd.plot(1,elev,clim=zlim,cb_aspect=60,cb_pad=0.01,levels=100, cmap='bwr',bnd=1)
           setp(gca(),xlim=xm,ylim=ym)
           title('Elevation (m): {:d} s'.format(mt[it]),fontsize=14,fontweight='bold')
        else:
           figure(figsize=[20,10])
           #map
           subplot(1,2,1)
           gd.plot(1,elev,clim=zlim,cb_aspect=60,cb_pad=0.01,levels=100, cmap='bwr',bnd=1)
           plot(*bp.xy.T,'g^',ms=12) #add bp pts
           for x,y,v in zip(bp.x,bp.y,bp.station): text(x,y,v,fontsize=12,color='g')
           setp(gca(),xlim=xm,ylim=ym)
           title('Elevation (m): {:d} s'.format(mt[it]),fontsize=14,fontweight='bold')

           #time series
           for m in arange(npt):
               subplot(max([npt,4]),2,2*m+2)
               my=mys[m]; y1=mys.min(); y2=mys.max(); dx=mt[-1]-mt[0]; dy=y2-y1; xmt=[mt[0],mt[-1]]; ymt=[y1-0.05*dy,y2+0.05*dy]
               plot(mt,0*mt,'k:',lw=0.1);plot(mt,my,'k-'); plot(mt[it],my[it],'ro',ms=10); plot([mt[it],mt[it]],ymt,'r')

               #note
               gca().xaxis.grid('on'); xshift=-0.06 if (mt[-1]-mt[it]<0.2*dx) else 0.02
               text(mt[it]+xshift*dx,mean(ymt),'{:d} s'.format(mt[it]),color='r',fontsize=12)
               if m==0: setp(gca(),xlim=xmt,ylim=ymt); xts=xticks()[0]; xls=[str(i) for i in xts]
               setp(gca(),xticks=xts,xticklabels=xls if m==(npt-1) else [],xlim=xmt,ylim=ymt)
               rtext(0.5,0.9,bp.station[m],fontsize=14,fontweight='bold',color='r')
        
        gcf().tight_layout()
        savefig('A_{:04}'.format(it),dpi=150); close()
        if it%10==0: print('done plot: {}/{}'.format(it,nt))

#--------------------------------------------------
#extract time series 
#--------------------------------------------------
if 'bpfile' in locals():
    bp=read(bpfile); npt=bp.npt; pip,pacor=gd.compute_acor(bp.xy)[1:]
    ds=[bp.npt,nt]; shm=SM(create=True,size=zeros(ds,'float32').nbytes); mys=np.ndarray(ds,'float32',buffer=shm.buf)
    for m in arange(int(nt/nproc)+1):
        i1=m*nproc; i2=min([(m+1)*nproc,nt]); nps=i2-i1+1
        if i1>=nt: continue
        pp=mp.Pool(nps); pp.map(fread,its[i1:i2])
    C=zdata(); C.time=mt; C.elev=mys; C.save('TS')
else:
   npt=0 

#--------------------------------------------------
#plot tsunami contour 
#--------------------------------------------------
nps=min([nt,nproc]); sindt=array_split(its,nps)
pp=mp.Pool(nps); pp.map(fplot,sindt)
shm.close(); shm.unlink()

#convert to gif
fns=[Image.open(i) for i in sorted(glob("A_*.png"))]
fns[0].save(fname+'.gif',save_all=True, append_images=fns[1:],duration=delay,loop=0)
for fn in fns: os.remove(fn)
