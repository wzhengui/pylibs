#!/usr/bin/env python3
#generate SCHISM 3D and 2D boundary condition based on Hycom data
from pylib import *

#----inputs--------------
StartT=datenum(2005,1,1)
EndT=datenum(2005,1,6)

Var=['surf_el','water_temp','salinity',['water_u','water_v']];
VarName=['elev','temp','salt',['Ux','Uy']];
ncname=['elev2D.th.nc','TEM_3D.th.nc','SAL_3D.th.nc','uv3D.th.nc']

dir_hycom='HYCOM/Data';

#---read grid----
gd=read_schism_hgrid('../hgrid.ll')

#index of pts for interpolation
node=gd.iobn[0] #open boundary node index
#node=r_[gd.iobn[0],gd.iobn[1]] #open boundary node index

#---get z coordinate------
zcor=read_schism_vgrid('../vgrid.in',gd,node=node,flag=1);
zcor=-zcor; fp=zcor<0; zcor[fp]=0; fp=zcor>5000; zcor[fp]=5000;
nvrt=zcor.shape[1]

#---interpolation pts----
dt=1/4; Time=arange(StartT,EndT+dt,dt);
loni=gd.x[node]; lati=gd.y[node];
bxy=c_[lati,loni]

lon2i=tile(loni,[nvrt,1]).T; lat2i=tile(lati,[nvrt,1]).T;
bxyz=c_[zcor.reshape(size(zcor)),lat2i.reshape(size(lat2i)),lon2i.reshape(size(lon2i))]

#---read Hycom data and interpolate onto boundary nodes------------------
Dt=[];
for i in arange(len(Var)):
    vari=Var[i]; varnamei=VarName[i]; ncnamei=ncname[i]
    #----build npz file for 3Dth or 2Dth----------
    nd=npz_data()
    nd.dimname=['nOpenBndNodes', 'nLevels', 'nComponents', 'one', 'time']
    if varnamei=='elev':
        nd.dims=[len(node),1,1,1,len(Time)]
    elif varnamei=='salt' or varnamei=='temp':
        nd.dims=[len(node),nvrt,1,1,len(Time)]
    elif isinstance(varnamei,list):
        nd.dims=[len(node),nvrt,2,1,len(Time)]
    nd.file_format='NETCDF4'

    #--time step, time, and time series----
    nd.vars=['time_step', 'time', 'time_series'];
    nd.time_step=npz_data()
    nd.time_step.attrs=['long_name'];nd.time_step.long_name='time step in seconds';
    nd.time_step.dimname=('one',); nd.time_step.val=array(dt*86400).astype('float32');
    nd.time=npz_data()
    nd.time.attrs=['long_name'];nd.time.long_name='simulation time in seconds';
    nd.time.dimname=('time',); nd.time.val=(Time-Time[0])*86400;
    nd.time_series=npz_data()
    nd.time_series.attrs=[];
    nd.time_series.dimname=('nComponents','nLevels','nOpenBndNodes','time');

    t0=time.time();
    T0=[]; Data0=[];
    for ti in arange(StartT,EndT+1):
        t1=num2date(ti); t2=num2date(ti+1-1/24/60);

        if isinstance(vari,list):
            fname='Hycom_{}_{}_{}.nc'.format(varnamei[0],t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))
            fname2='Hycom_{}_{}_{}.nc'.format(varnamei[1],t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))
            if not os.path.exists(r'{}/{}'.format(dir_hycom,fname)): continue
            if not os.path.exists(r'{}/{}'.format(dir_hycom,fname2)): continue
            print(fname+'; '+fname2)
            C=ReadNC('{}/{}'.format(dir_hycom,fname),2); C2=ReadNC('{}/{}'.format(dir_hycom,fname2),2)

            #get value
            exec('val=C.{}.val'.format(vari[0]))
            exec('val2=C2.{}.val'.format(vari[1]))
            val=array(val); fp=val<-29999; val2=array(val2); fp2=val2<-29999;
            val[fp]=nan; val2[fp2]=nan;
        else:
            fname='Hycom_{}_{}_{}.nc'.format(varnamei,t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))
            if not os.path.exists(r'{}/{}'.format(dir_hycom,fname)): continue
            print(fname)
            C=ReadNC('{}/{}'.format(dir_hycom,fname),2)

            #get value
            exec('val=C.{}.val'.format(vari))
            val=array(val); fp=val<-29999;
            val[fp]=nan


        ti=datestr2num(C.time.time_origin)+array(C.time.val)/24
        cloni=array(C.lon.val); cloni=mod(cloni,360)-360
        clati=array(C.lat.val);

        #------define data region extracted
        ind_lon=nonzero((cloni<=max(loni)+0.1)*(cloni>=min(loni)-0.1))[0];
        ind_lat=nonzero((clati<=max(lati)+0.1)*(clati>=min(lati)-0.1))[0];
        i1_lon=ind_lon.min(); i2_lon=i1_lon+len(ind_lon)
        i1_lat=ind_lat.min(); i2_lat=i1_lat+len(ind_lat)

        cloni=cloni[i1_lon:i2_lon]; clati=clati[i1_lat:i2_lat]

        if varnamei=='elev':
            for m in arange(len(ti)):
                valii=squeeze(val[m,i1_lat:i2_lat,i1_lon:i2_lon])

                #interpolation
                fd=sp.interpolate.RegularGridInterpolator((clati,cloni),valii,fill_value=nan)
                vi=fd(bxy)

                #remove nan pts
                fp=isnan(vi);
                if sum(fp)!=0:
                    vi[fp]=sp.interpolate.griddata(bxy[~fp,:],vi[~fp],bxy[fp,:],'nearest')

                T0.append(ti[m]); Data0.append(vi);
        else:
            #------define data region extracted for depth
            cdepi=array(C.depth.val)
            ind_dep=nonzero((cdepi<=zcor.max()+1000)*(cdepi>=zcor.min()-100))[0];
            i1_dep=ind_dep.min(); i2_dep=i1_dep+len(ind_dep)
            cdepi=cdepi[i1_dep:i2_dep];

            for m in arange(len(ti)):
                T0.append(ti[m]);
                valii=squeeze(val[m,i1_dep:i2_dep,i1_lat:i2_lat,i1_lon:i2_lon])

                #interpolation
                fd=sp.interpolate.RegularGridInterpolator((cdepi,clati,cloni),valii,fill_value=nan)
                vi=fd(bxyz)

                #remove nan pts
                fp=isnan(vi);
                if sum(fp)!=0:
                    vi[fp]=sp.interpolate.griddata(bxyz[~fp,:],vi[~fp],bxyz[fp,:],'nearest')

                vi=vi.reshape(zcor.shape)

                #----if variable is velocity
                if isinstance(varnamei,list):
                    val2ii=squeeze(val2[m,i1_dep:i2_dep,i1_lat:i2_lat,i1_lon:i2_lon])

                    #interpolation
                    fd=sp.interpolate.RegularGridInterpolator((cdepi,clati,cloni),val2ii,fill_value=nan)
                    v2i=fd(bxyz)

                    #remove nan pts
                    fp=isnan(v2i);
                    if sum(fp)!=0:
                        v2i[fp]=sp.interpolate.griddata(bxyz[~fp,:],v2i[~fp],bxyz[fp,:],'nearest')

                    v2i=v2i.reshape(zcor.shape)
                    Data0.append(r_[expand_dims(vi,0),expand_dims(v2i,0)])

                else:
                    Data0.append(vi)

    T0=array(T0); Data0=array(Data0)

    #---check whether there is nan
    y0=Data0.reshape([size(Data0)]); fp=isnan(y0)
    if sum(fp)!=0:
        print('{} has NaN: check'.format(varnamei))
        sys.exit()

    #interpolation in time
    ds=Data0.shape; y0=Data0.reshape([ds[0],prod(ds[1:])]);
    #lpfilter
    if varnamei=='elev' or isinstance(varnamei,list):
        y0=lpfilt(y0,dt,0.9)
    fd=interpolate.interp1d(T0,y0,axis=0,fill_value='extrapolate');
    Data=reshape(fd(Time),[len(Time),*ds[1:]]);

    #----put data into ncfile
    if varnamei=='elev':
        Data=Data[:,:,None,None].transpose([3,2,1,0])
    elif varnamei=='salt' or varnamei=='temp':
        Data=Data[:,:,:,None].transpose([3,2,1,0])
    elif isinstance(varnamei,list):
        Data=Data.transpose([1,3,2,0])
    nd.time_series.val=Data.astype('float32');

    #--write ncfile----
    if os.path.exists(ncnamei): os.remove(ncnamei)
    WriteNC(nd,ncnamei,med=2,order=1)

    Dt.append(time.time()-t0);

#----print time consumed--------------------
for i in arange(len(Dt)):
    varnamei=VarName[i]
    print('reading {}: {}'.format(varnamei,Dt[i]))

#---dist---------------------------------
#bP=gd.x[gd.iobn[0]]+1j*gd.y[gd.iobn[0]];
#L=zeros(bP.shape);
#for i in arange(gd.nobn[0]-1):
#    L[i+1]=L[i]+abs(bP[i+1]-bP[i])

##----plot boundary grid---------------
#for i in arange(zcor.shape[0]):
#    xi=ones(zcor.shape[1])*L[i]
#    zi=zcor[i,:]
#    plot(xi,zi,'k-')
#
#for i in arange(zcor.shape[1]):
#    xi=L
#    zi=zcor[:,i]
#    plot(xi,zi,'k-')
