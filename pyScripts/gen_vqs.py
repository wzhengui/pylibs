#!/usr/bin/env
from pylib import *
close('all')

#------------------------------------------------------------------------------
#input
#------------------------------------------------------------------------------
grd='./hgrid.gr3'        #hgrid name (grid.npz, or hgrid.gr3)

#hsm: depth for each master grid;  nv: number for layers for each depth
hsm=array([0.5,1.0,1.5, *arange(2,11),*[11.5,14,17,21],*arange(25,51,5),*arange(60,81,10),95,125])
nhm=len(hsm); nv=2+arange(nhm)

#check transect
bname='./transect.bp'   #transect bpfile
#------------------------------------------------------------------------------
#compute master grid
#------------------------------------------------------------------------------
nvrt=nv[-1]; z_mas=ones([nhm,nvrt])*nan; eta=0.0
for m, [hsmi,nvi] in enumerate(zip(hsm,nv)):
    #strethcing funciton
    theta_b=0; theta_f=2.5; hc=10.5
    hc=min(hsmi,hc)

    for k in arange(nvi):
        sigma= k/(1-nvi)  #zi=-sigma #original sigma coordiante

        #compute zcoordinate
        cs=(1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)+theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
        z_mas[m,k]=eta*(1+sigma)+hc*sigma+(hsmi-hc)*cs

    #normalize z_mas
    z_mas[m]=-(z_mas[m]-z_mas[m,0])*hsmi/(z_mas[m,nvi-1]-z_mas[m,0])
s_mas=array([z_mas[i]/hsm[i] for i in arange(nhm)])

#check master grid
for i in arange(nhm-1):
    if min(z_mas[i,:nv[i]]-z_mas[i+1,:nv[i]])<0: print('check: master grid layer={}, hsm={}, nv={}'.format(i+1,hsm[i+1],nv[i+1]))

#plot master grid
figure(figsize=[10,5])
for i in arange(nhm): plot(i*ones(nvrt),z_mas[i],'k-',lw=0.3)
for k in arange(nvrt): plot(arange(nhm),z_mas.T[k],'k-',lw=0.3)
setp(gca(),xlim=[-0.5,nhm-0.5],ylim=[-hsm[-1],0.5])
gcf().tight_layout()
move_figure(gcf(),0,0)
savefig('Master_Grid',dpi=200)

# sys.exit()

#------------------------------------------------------------------------------
#compute vgrid
#------------------------------------------------------------------------------
#read hgrid
gd=loadz(grd).hgrid if grd.endswith('.npz') else read_schism_hgrid(grd)
fpz=gd.dp<hsm[0]; gd.dp[fpz]=hsm[0]

#find hsm index for all points
rat=ones(gd.np)*nan; nlayer=zeros(gd.np).astype('int');
ind1=zeros(gd.np).astype('int'); ind2=zeros(gd.np).astype('int')
for m, hsmi in enumerate(hsm):
    if m==0:
        fp=gd.dp<=hsm[m];
        ind1[fp]=0; ind2[fp]=0; rat[fp]=0; nlayer[fp]=nv[0]
    else:
        fp=(gd.dp>hsm[m-1])*(gd.dp<=hsm[m])
        ind1[fp]=m-1; ind2[fp]=m
        rat[fp]=(gd.dp[fp]-hsm[m-1])/(hsm[m]-hsm[m-1]); nlayer[fp]=nv[m]

znd=z_mas[ind1]*(1-rat[:,None])+z_mas[ind2]*rat[:,None]; #z coordinate
for i in arange(gd.np): znd[i,nlayer[i]-1]=-gd.dp[i]
snd=znd/gd.dp[:,None]; #sigma coordinate

#check vgrid
for i in arange(gd.np):
    for k in arange(nvrt-1):
        if znd[i,k]<=znd[i,k+1]:
            sys.exit('wrong vertical layers')

#write vgrid.in
fid=open('vgrid.in','w+')
fid.write('1  !average # of layers={:0.2f}\n{}\n'.format(mean(nlayer),nvrt))
for i in arange(gd.np):
    nlayeri=nlayer[i]; si=flipud(snd[i,:nlayeri])
    fstr='{:6}    {:2}   '+'{:10.6f}  '*nlayeri+'\n'
    fid.write(fstr.format(i+1,nvrt-nlayeri+1,*si))
fid.close()
print('Average number of layers is: {:0.2f}'.format(mean(nlayer)))

#------------------------------------------------------------------------------
#plot transect
#------------------------------------------------------------------------------
if os.path.exists(str(bname)):
    bp=read_schism_bpfile(str(bname))

    #compute dist
    dist=[0,]
    for i in arange(bp.nsta-1):
        disti=abs((bp.x[i+1]-bp.x[i])+1j*(bp.y[i+1]-bp.y[i]))+dist[i]
        dist.append(disti)
    dist=array(dist)

    #compute zcor
    sindp=near_pts(c_[bp.x,bp.y],c_[gd.x,gd.y]); zi=znd[sindp]
    for i in arange(bp.nsta): fpn=isnan(zi[i]); zi[i][fpn]=min(zi[i])


    #plot
    figure(figsize=[10,5])
    for k in arange(nvrt): plot(dist,zi[:,k],'k',lw=0.5)
    for i in arange(bp.nsta): plot(ones(nvrt)*dist[i],zi[i],'k',lw=0.5)
    setp(gca(),ylim=[zi.min()-1,0.5],xlim=[0,dist.max()])
    gcf().tight_layout()
    # move_figure(gcf(),0,0)
