#!/usr/bin/env python3
from pylib import *

#----inputs--------------
dir_hycom='./HYCOM/Data';
date_hycom=datenum(2005,1,1) #date for hycom file

Var=['water_temp','salinity'];
VarName=['temp','salt'];

t0=time.time();
#---read grid----
gd=read_schism_hgrid('hgrid.ll')

#index of pts for interpolation
node=arange(gd.np)

#---get z coordinate------
zcor=read_schism_vgrid('vgrid.in',gd,node=node,flag=1);
zcor=-zcor; fp=zcor<0; zcor[fp]=0; fp=zcor>5000; zcor[fp]=5000;
nvrt=zcor.shape[1]

#---interpolation pts----
loni=gd.x[node]; lati=gd.y[node];
bxy=c_[lati,loni]

lon2i=tile(loni,[nvrt,1]).T; lat2i=tile(lati,[nvrt,1]).T;
bxyz=c_[zcor.reshape(size(zcor)),lat2i.reshape(size(lat2i)),lon2i.reshape(size(lon2i))]

#--------interpolation --------------------------------------------------------
#interpolate is done only for salinity and temperature,but could modified to include elevation and velocity
Data=npz_data()
for i in arange(len(Var)):
    vari=Var[i]; varnamei=VarName[i];
    t1=num2date(date_hycom); t2=num2date(date_hycom+1-1/24/60)
    fname='{}/Hycom_{}_{}_{}.nc'.format(dir_hycom,varnamei,t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))

    if os.path.exists(fname):
        print(fname)

    C=ReadNC(fname,2)

    #get value
    exec('val=C.{}.val'.format(vari))
    val=array(val); fp=val<-29999;
    val[fp]=nan

    ti=datestr2num(C.time.time_origin)+array(C.time.val)/24

    cdepi=array(C.depth.val)
    clati=array(C.lat.val);
    cloni=array(C.lon.val); cloni=mod(cloni,360)-360

    #------define data region extracted----
    ind_dep=nonzero((cdepi<=zcor.max()+1000)*(cdepi>=zcor.min()-100))[0];
    ind_lat=nonzero((clati<=max(lati)+0.1)*(clati>=min(lati)-0.1))[0];
    ind_lon=nonzero((cloni<=max(loni)+0.1)*(cloni>=min(loni)-0.1))[0];
    i1_dep=ind_dep.min(); i2_dep=i1_dep+len(ind_dep)
    i1_lat=ind_lat.min(); i2_lat=i1_lat+len(ind_lat)
    i1_lon=ind_lon.min(); i2_lon=i1_lon+len(ind_lon)

    cdepi=cdepi[i1_dep:i2_dep]; cloni=cloni[i1_lon:i2_lon]; clati=clati[i1_lat:i2_lat];

    #----extract hycom data for the defined region
    valii=squeeze(val[0,i1_dep:i2_dep,i1_lat:i2_lat,i1_lon:i2_lon])

    #interpolation
    fd=sp.interpolate.RegularGridInterpolator((cdepi,clati,cloni),valii,fill_value=nan)
    fdn=sp.interpolate.RegularGridInterpolator((cdepi,clati,cloni),valii,'nearest')
    vi=fd(bxyz); vin=fdn(bxyz)

    #remove nan pts
    fp=isnan(vi);
    if sum(fp)!=0:
        #vi[fp]=vin[fp]
        vi[fp]=sp.interpolate.griddata(bxyz[~fp,:],vi[~fp],bxyz[fp,:],'nearest')

    fp=isnan(vi);
    if sum(fp)!=0: sys.exit()

    vi=vi.reshape(zcor.shape)
    exec('Data.{}=vi'.format(varnamei))

#---compute tr_nd and tr_el------
tr_nd=c_[Data.temp[:,:,None],Data.salt[:,:,None]].transpose([2,1,0])
tr_el=zeros([2,nvrt,gd.ne])
for i in arange(2):
    for k in arange(nvrt):
        gd.dp=squeeze(tr_nd[i,k,:]);
        gd.compute_ctr();
        tr_el[i,k,:]=gd.dpe;

#-----build hotstart.nc--------------------------------------------------------
nd=npz_data()
nd.dimname=['node','elem','side','nVert','ntracers','one']
nd.dims=[gd.np,gd.ne,gd.ns,nvrt,2,1]
nd.file_format='NETCDF4'

#--time step, time, and time series----
nd.vars=['time','iths','ifile','idry_e','idry_s','idry','eta2','we','tr_el','tr_nd',\
         'tr_nd0','su2','sv2','q2','xl','dfv','dfh','dfq1','dfq2']
ndi=npz_data(); ndi.attrs=[]; #template

nd.time=npz_data(); nd.time.dimname=('one',);nd.time.val=array(0.0) #time
nd.iths=npz_data(); nd.iths.dimname=('one',);nd.iths.val=array(0) #iths
nd.ifile=npz_data(); nd.ifile.dimname=('one',);nd.ifile.val=array(1) #ifile

nd.idry_e=npz_data(); nd.idry_e.dimname=('elem',);nd.idry_e.val=zeros(gd.ne).astype('int32') #idry_e
nd.idry_s=npz_data(); nd.idry_s.dimname=('side',);nd.idry_s.val=zeros(gd.ns).astype('int32') #idry_s
nd.idry=npz_data(); nd.idry.dimname=('node',);nd.idry.val=zeros(gd.np).astype('int32') #idry
nd.eta2=npz_data(); nd.eta2.dimname=('node',);nd.eta2.val=zeros(gd.np) #eta2

nd.we=npz_data(); nd.we.dimname=('nVert','elem');nd.we.val=zeros([nvrt,gd.ne]) #eta2

nd.tr_el=npz_data(); nd.tr_el.dimname=('ntracers','nVert','elem');nd.tr_el.val=tr_el #tr_el
nd.tr_nd=npz_data(); nd.tr_nd.dimname=('ntracers','nVert','node');nd.tr_nd.val=tr_nd #tr_nd
nd.tr_nd0=npz_data(); nd.tr_nd0.dimname=('ntracers','nVert','node');nd.tr_nd0.val=tr_nd #tr_nd0

nd.su2=npz_data(); nd.su2.dimname=('nVert','side');nd.su2.val=zeros([nvrt,gd.ns]) #su2
nd.sv2=npz_data(); nd.sv2.dimname=('nVert','side');nd.sv2.val=zeros([nvrt,gd.ns]) #sv2

nd.q2=npz_data(); nd.q2.dimname=('nVert','node');nd.q2.val=zeros([nvrt,gd.np]) #q2
nd.xl=npz_data(); nd.xl.dimname=('nVert','node');nd.xl.val=zeros([nvrt,gd.np]) #xl
nd.dfv=npz_data(); nd.dfv.dimname=('nVert','node');nd.dfv.val=zeros([nvrt,gd.np]) #dfv
nd.dfh=npz_data(); nd.dfh.dimname=('nVert','node');nd.dfh.val=zeros([nvrt,gd.np]) #dfh
nd.dfq1=npz_data(); nd.dfq1.dimname=('nVert','node');nd.dfq1.val=zeros([nvrt,gd.np]) #dfq1
nd.dfq2=npz_data(); nd.dfq2.dimname=('nVert','node');nd.dfq2.val=zeros([nvrt,gd.np]) #dfq2

for vari in nd.vars:
    exec('nd.{}.attrs=[]'.format(vari))

#--write ncfile----
ncname='hotstart.nc'
if os.path.exists(ncname): os.remove(ncname)
WriteNC(nd,ncname,med=2,order=1)

#---time used--------
print('time consumed for generating hotstart: {} s'.format(time.time()-t0));

sys.exit()
#----plot for check------------------------------------------------------------
for i in arange(12):
    subplot(3,4,i+1)
    gd.dp=vi[:,i*4]
    gd.compute_ctr()
    gd.plot_grid(plotz=1,ec=None,clim=[10,35]);
#    colorbar(gd.hc)
