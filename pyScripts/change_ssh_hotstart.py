#!/usr/bin/env python3
# This script replaces Sea Surface Height (SSH) in hotstart.nc with new SSH from another large-scale model
from pylib import *

########################################################
grd='./hgrid.gr3' #grd file
lmname='./dt_global_allsat_phy_l4_20160915.nc' #large-scale model for new SSH
sshname='adt' # name of SSH in the large-scale model
hname='./hotstart.nc' # old hotstart.nc
sname='./hotstart_AVISO_ssh.nc' # new hotstart.nc to be saved
########################################################

#------ Load data
gd=read_schism_hgrid('{}'.format(grd))
oldhot=ReadNC('{}'.format(hname))
C=ReadNC('{}'.format(lmname),1)
lxi0=gd.x; lyi0=gd.y; bxy=c_[lxi0,lyi0]
sx=array(C.variables['longitude'][:])
if sx.max() > 180: sx=sx-360 # only work for Atlantic ocean when using AVISO
sy=array(C.variables['latitude'][:]);

#------ Perform interpolation
#get interp index 
sxi,syi=meshgrid(sx,sy); sxy=c_[sxi.ravel(),syi.ravel()];
cvs=array(C.variables[sshname][0]); sindns=[]; sindps=[]
print('computing interpation index')
cv=cvs; ds=cv.shape; cv=cv.ravel()
fpn=abs(cv)>1e3; sindn=nonzero(fpn)[0]; sindr=nonzero(~fpn)[0]; sindp=sindr[near_pts(sxy[sindn],sxy[sindr])]
sindns.append(sindn); sindps.append(sindp)

#get interp index for pts
sx0=sx[:]; sy0=sy[:]; print('get new interp indices')
idx0=((lxi0[:,None]-sx0[None,:])>=0).sum(axis=1)-1; ratx0=(lxi0-sx0[idx0])/(sx0[idx0+1]-sx0[idx0])
idy0=((lyi0[:,None]-sy0[None,:])>=0).sum(axis=1)-1; raty0=(lyi0-sy0[idy0])/(sy0[idy0+1]-sy0[idy0])
exec("cv=array(C.variables['{}'][0])".format(sshname)); 
sindn,sindp=sindns[0],sindps[0]
cv=cv.ravel(); fpn=(abs(cv[sindn])>1e3)*(abs(cv[sindp])<1e3); cv[sindn]=cv[sindp]; fpn=abs(cv)>1e3 #init fix
if sum(fpn)!=0: fni=nonzero(fpn)[0]; fri=nonzero(~fpn)[0]; fpi=fri[near_pts(sxy[fni],sxy[fri])]; cv[fni]=cv[fpi] #final fix
cv=cv.reshape(ds)

#find parent pts
v0=array([cv[idy0,idx0],cv[idy0,idx0+1],cv[idy0+1,idx0],cv[idy0+1,idx0+1]])

#interp
v1=v0[0]*(1-ratx0)+v0[1]*ratx0;  v2=v0[2]*(1-ratx0)+v0[3]*ratx0
vi=v1*(1-raty0)+v2*raty0 
vi=array(vi) # interpolated SSH on nodes

#replace SSH with new SSH from large-scale model that you choose  
oldhot.eta2.val[:]=vi
oldhot.cumsum_eta.val[:]=vi

#save new hostart.nc
WriteNC('{}'.format(sname),oldhot)

print('--------done--------')

