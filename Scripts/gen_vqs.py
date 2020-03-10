from pylib import *

fname='./SF_1.gr3'; #grid name
vname='vgrid_1.in' #name for vgrid
bname='./transect1.bp'; #bp file

#------------------------------------------------------------------------------
#----Master Grid---------------------------------------------------------------
#------------------------------------------------------------------------------
hsm=[0.45, 0.9, 1.35, 1.9, 2.5,3.25,4,5,6,7,8.1,9.25,10.5];
[hsm.append(i) for i in arange(12,21,2)]
[hsm.append(i) for i in arange(22,41,3)]
[hsm.append(i) for i in arange(45,75,5)];
[hsm.append(i) for i in arange(80,130,10)]

hsm=array(hsm)        #depth for each master grid
nv=2+arange(len(hsm)) #number for layers

#--generate master grid-----------------
nvrt=nv[-1];
z_mas=ones([len(hsm),nvrt])*nan

eta=0.0
for m in arange(len(nv)):
    zsum=0.0
    for k in arange(nv[m]):
        sigma= k/(1-nv[m]) #original sigma coordiante
        zi=-sigma
        #strethcing funciton
        if m>=3:
            theta_b=0
            theta_f=1.25+hsm[m]*0.01
        else:
            theta_b=0
            theta_f=1
        cs=(1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)+theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
        z_mas[m,k]=eta*(1+sigma)+hsm[1]*sigma+(hsm[m]-hsm[1])*cs

    #normalize z_mas
    zsum=z_mas[m,nv[m]-1]-z_mas[m,0]
    z_mas[m,:]=-(z_mas[m,:]-z_mas[m,0])*hsm[m]/zsum

s_mas=z_mas.copy()
for m in arange(len(hsm)):
    s_mas[m,:]=s_mas[m,:]/hsm[m]

#------------------------------------------------------------------------------
#---compute zcoor--------------------------------------------------------------
#------------------------------------------------------------------------------
gd=read_schism_hgrid(fname);
#find hsm index for all points
rat=ones(gd.np)*nan; nlayer=rat.copy().astype('int');
ind1=rat.copy().astype('int'); ind2=rat.copy().astype('int')
for m in arange(len(hsm)):
    if m==0:
        fp=gd.dp<=hsm[m];
        ind1[fp]=0; ind2[fp]=0; rat[fp]=0; nlayer[fp]=nv[0]
    else:
        fp=(gd.dp>hsm[m-1])*(gd.dp<=hsm[m])
        ind1[fp]=m-1; ind2[fp]=m
        rat[fp]=(gd.dp[fp]-hsm[m-1])/(hsm[m]-hsm[m-1]); nlayer[fp]=nv[m]

znd=z_mas[ind1]*(1-rat[:,None])+z_mas[ind2]*rat[:,None]; #z coordinate
for i in arange(gd.np):
    znd[i,nlayer[i]-1]=-gd.dp[i]
snd=znd/gd.dp[:,None]; #sigma coordinate

#check
for i in arange(gd.np):
    for k in arange(nvrt-1):
        if znd[i,k]<=znd[i,k+1]:
            print('wrong vertical layers')
            sys.exit()

#----write vgrid.in-----------------
with open(vname,'w+') as fid:
    fid.write('   1 !average # of layers={:0.2f}\n'.format(sum(nlayer)/gd.np))
    fid.write('   {}\n'.format(nvrt))
    for i in arange(gd.np):
        nlayeri=nlayer[i]
        si=flipud(snd[i,:nlayeri]);
        fid.write('{:6}    {:2}   '.format(i+1,nvrt-nlayeri+1))
        formatstr='{:10.6f}  '*nlayeri
        fid.write(formatstr.format(*si))
        fid.write('\n')

print('Average number of layers is: {:0.2f}'.format(sum(nlayer)/gd.np));

#------------------------------------------------------------------------------
#---plot master grid-----------------------------------------------------------
#------------------------------------------------------------------------------
close('all')
figure(figsize=[10,5])

for i in arange(len(hsm)):
    xi=ones(nvrt)*i+1
    if hsm[i]<=3:
        plot(xi,z_mas[i,:],'r-',lw=0.3)
    elif hsm[i]<=30:
        plot(xi,z_mas[i,:],'k-',lw=0.3)
    elif hsm[i]<=40:
        plot(xi,z_mas[i,:],'k-',lw=0.3)
    elif hsm[i]<=60:
        plot(xi,z_mas[i,:],'k-',lw=0.3)
    elif hsm[i]<=100:
        plot(xi,z_mas[i,:],'k-',lw=0.3)
    else:
        plot(xi,z_mas[i,:],'k-',lw=0.3)
for i in arange(nvrt):
    xi=arange(len(hsm))+1;
    plot(xi,squeeze(z_mas[:,i]),'k-',lw=0.3)
setp(gca(),xlim=[0,len(hsm)],ylim=[-120,0])


gcf().tight_layout()
move_figure(gcf(),0,0)
savefig('Master_Grid',dpi=200)

#------------------------------------------------------------------------------
#---plot transect--------------------------------------------------------------
#------------------------------------------------------------------------------
bp=read_schism_bpfile(bname);

#-find nearest point------------------
P=gd.x+1j*gd.y;
b=bp.x+1j*bp.y;
P=reshape(P,[1,len(P)])
b=reshape(b,[len(b),1])
Dist=abs(P-b);
bind=[];L=[0]
for i in arange(bp.nsta):
    if i>0:
        L.append(L[i-1]+abs(b[i]-b[i-1]))
    disti=Dist[i];
    ind=nonzero(disti==min(disti))[0]
    bind.append(ind)
L=array(L); bind=array(bind)
zb=squeeze(znd[bind]);
#extend values-----------
for i in arange(bp.nsta):
    fp=isnan(zb[i]); ind=max(nonzero(~fp)[0]);
    zb[i,fp]=zb[i,ind]

#---plot transect---------
figure(figsize=[10,4])
for i in arange(bp.nsta):
    xi=ones(nvrt)*L[i];
    plot(xi,zb[i,:],'k-',lw=0.4)
for i in arange(nvrt):
    plot(L,zb[:,i],'k-',lw=0.4)

setp(gca(),ylim=[-120,0])
gcf().tight_layout()
move_figure(gcf(),0,0)
