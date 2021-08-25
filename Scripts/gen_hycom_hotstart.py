#!/usr/bin/env python3
#create hotstart condition based on hycom data
from pylib import *
close("all")

#------------------------------------------------------------------------------
#input
#------------------------------------------------------------------------------
StartT=datenum(2010,1,1)
grd='../grid.npz'
dir_hycom='Data_init'

#------------------------------------------------------------------------------
#interpolate hycom data to boundary
#------------------------------------------------------------------------------
#variables to be interpolated
svars=['water_temp','salinity']
mvars=['temp','salt']

#find hycom file
fnames=array([i for i in os.listdir(dir_hycom) if i.endswith('.nc')])
mti=array([datenum(*array(i.replace('.','_').split('_')[1:5]).astype('int')) for i in fnames])
fpt=nonzero(abs(mti-StartT)==min(abs(mti-StartT)))[0][0]; fname=fnames[fpt]

#read hgrid
gd=loadz(grd).hgrid; vd=loadz(grd).vgrid; gd.x,gd.y=gd.lon,gd.lat
ne,np,ns,nvrt=gd.ne,gd.np,gd.ns,vd.nvrt

#get node xyz
lxi=gd.x%360; lyi=gd.y; lzi0=abs(vd.compute_zcor(gd.dp)).T

#get hycom time, xyz
C=ReadNC('{}/{}'.format(dir_hycom,fname),1); #print(fname)
ctime=array(C.variables['time'])/24+datenum(2000,1,1); sx=array(C.variables['lon'][:])%360
sy=array(C.variables['lat'][:]); sz=array(C.variables['depth'][:])
fpz=lzi0>sz.max(); lzi0[fpz]=sz.max()-1e-6

#interp for ST
S=zdata(); [exec('S.{}=[]'.format(i)) for i in mvars]
for k in arange(nvrt):
    lzi=lzi0[k]; bxyz=c_[lxi,lyi,lzi]

    #get interp index
    idx=((lxi[:,None]-sx[None,:])>=0).sum(axis=1)-1; ratx=(lxi-sx[idx])/(sx[idx+1]-sx[idx])
    idy=((lyi[:,None]-sy[None,:])>=0).sum(axis=1)-1; raty=(lyi-sy[idy])/(sy[idy+1]-sy[idy])
    idz=((lzi[:,None]-sz[None,:])>=0).sum(axis=1)-1; ratz=(lzi-sz[idz])/(sz[idz+1]-sz[idz])

    #for each variable
    for m,svar in enumerate(svars):
        print(svar,k)
        mvar=mvars[m]
        exec("cv=array(C.variables['{}'][0])".format(svar))
        v0=array([cv[idz,idy,idx],cv[idz,idy,idx+1],cv[idz,idy+1,idx],cv[idz,idy+1,idx+1],
              cv[idz+1,idy,idx],cv[idz+1,idy,idx+1],cv[idz+1,idy+1,idx],cv[idz+1,idy+1,idx+1]])

        #remove nan pts
        for n in arange(8):
            fpn=abs(v0[n])>1e3
            v0[n,fpn]=sp.interpolate.griddata(bxyz[~fpn,:],v0[n,~fpn],bxyz[fpn,:],'nearest',rescale=True)

        v11=v0[0]*(1-ratx)+v0[1]*ratx;  v12=v0[2]*(1-ratx)+v0[3]*ratx; v1=v11*(1-raty)+v12*raty
        v21=v0[4]*(1-ratx)+v0[5]*ratx;  v22=v0[6]*(1-ratx)+v0[7]*ratx; v2=v21*(1-raty)+v22*raty
        vi=v1*(1-ratz)+v2*ratz

        #save
        exec('S.{}.append(vi)'.format(mvar))
[exec('S.{}=array(S.{})'.format(i,i)) for i in mvars]

#------------------------------------------------------------------------------
#creat netcdf
#------------------------------------------------------------------------------
nd=zdata()
nd.dimname=['node','elem','side','nVert','ntracers','one']; nd.dims=[np,ne,ns,nvrt,2,1]
nd.file_format='NETCDF4'

#--time step, time, and time series----
nd.vars=['time','iths','ifile','idry_e','idry_s','idry','eta2','we','tr_el',
         'tr_nd','tr_nd0','su2','sv2','q2','xl','dfv','dfh','dfq1','dfq2']
for i in nd.vars: exec('nd.{}=zdata(); nd.{}.attrs=[]'.format(i,i))

nd.time.dimname=('one',);nd.time.val=array(0.0) #time
nd.iths.dimname=('one',);nd.iths.val=array(0) #iths
nd.ifile.dimname=('one',);nd.ifile.val=array(1) #ifile

nd.idry_e.dimname=('elem',);nd.idry_e.val=zeros(ne).astype('int32') #idry_e
nd.idry_s.dimname=('side',);nd.idry_s.val=zeros(ns).astype('int32') #idry_s
nd.idry.dimname=('node',);nd.idry.val=zeros(np).astype('int32') #idry
nd.eta2.dimname=('node',);nd.eta2.val=zeros(np) #eta2

nd.we.dimname=('elem','nVert');nd.we.val=zeros([ne,nvrt]) #we

nd.su2.dimname=('side','nVert');nd.su2.val=zeros([ns,nvrt]) #su2
nd.sv2.dimname=('side','nVert');nd.sv2.val=zeros([ns,nvrt]) #sv2

nd.q2.dimname=('node','nVert');nd.q2.val=zeros([np,nvrt]) #q2
nd.xl.dimname=('node','nVert');nd.xl.val=zeros([np,nvrt]) #xl
nd.dfv.dimname=('node','nVert');nd.dfv.val=zeros([np,nvrt]) #dfv
nd.dfh.dimname=('node','nVert');nd.dfh.val=zeros([np,nvrt]) #dfh
nd.dfq1.dimname=('node','nVert');nd.dfq1.val=zeros([np,nvrt]) #dfq1
nd.dfq2.dimname=('node','nVert');nd.dfq2.val=zeros([np,nvrt]) #dfq2

tr_nd=r_[S.temp[None,...],S.salt[None,...]].T; tr_el=tr_nd[gd.elnode[:,:3]].mean(axis=1)
nd.tr_el.dimname=('elem','nVert','ntracers'); nd.tr_el.val=tr_el #tr_el
nd.tr_nd.dimname=('node','nVert','ntracers'); nd.tr_nd.val=tr_nd #tr_nd
nd.tr_nd0.dimname=('node','nVert','ntracers');nd.tr_nd0.val=tr_nd #tr_nd0

WriteNC('hotstart.nc',nd)

