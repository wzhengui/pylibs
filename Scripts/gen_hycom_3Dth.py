#!/usr/bin/env python3
#create boundary condition based on hycom data
from pylib import *
close("all")

#------------------------------------------------------------------------------
#input
#------------------------------------------------------------------------------
StartT=datenum(2021,4,5); EndT=datenum(2021,4,10); dt=1/8
grd='../grid.npz'
dir_hycom='../../../HYCOM/Data'

#------------------------------------------------------------------------------
#interpolate hycom data to boundary
#------------------------------------------------------------------------------
mtime=arange(StartT,EndT+dt,dt); nt=len(mtime)

#variables for each files
snames=['elev2D.th.nc','TEM_3D.th.nc','SAL_3D.th.nc','uv3D.th.nc']
svars=['surf_el','water_temp','salinity',['water_u','water_v']]
mvars=['elev','temp','salt',['u','v']]

#find all hycom files
fnames=array([i for i in os.listdir(dir_hycom) if i.endswith('.nc')])
mti=array([datenum(*array(i.replace('.','_').split('_')[1:5]).astype('int')) for i in fnames])
fpt=(mti>=(StartT-1))*(mti<(EndT+1)); fnames=fnames[fpt]; mti=mti[fpt]
sind=argsort(mti); mti=mti[sind]; fnames=fnames[sind]

#read hgrid
gd=loadz(grd).hgrid; vd=loadz(grd).vgrid; gd.x,gd.y=gd.lon,gd.lat; nvrt=vd.nvrt

#get bnd node xyz
bind=gd.iobn[0]; nobn=gd.nobn[0]
lxi0=gd.x[bind]; lyi0=gd.y[bind]; bxy=c_[lxi0,lyi0] #for 2D
lxi=tile(lxi0,[nvrt,1]).T.ravel(); lyi=tile(lyi0,[nvrt,1]).T.ravel() #for 3D
lzi=abs(compute_zcor(vd.sigma,gd.dp[bind],ivcor=2,vd=vd)).ravel();  bxyz=c_[lxi,lyi,lzi]
#lzi=abs(compute_zcor(vd.sigma[bind],gd.dp[bind])).ravel();  bxyz=c_[lxi,lyi,lzi]
sx0,sy0,sz0=None,None,None

#for each variables
for n,sname in enumerate(snames):
    svar=svars[n]; mvar=mvars[n]
    if isinstance(svar,str): svar=[svar]; mvar=[mvar]

    #interp in space
    S=zdata(); S.time=[]
    [exec('S.{}=[]'.format(i)) for i in mvar]
    for m,fname in enumerate(fnames):
        C=ReadNC('{}/{}'.format(dir_hycom,fname),1); print(fname)
        ctime=array(C.variables['time'])/24+datenum(2000,1,1); sx=mod(array(C.variables['lon'][:])+360,360)-360
        sy=array(C.variables['lat'][:]); sz=array(C.variables['depth'][:])

        #get interp index
        if not array_equal(sx,sx0)*array_equal(sy,sy0)*array_equal(sz,sz0):
            sx0=sx[:]; sy0=sy[:]; sz0=sz[:]; print('get new interp indices: {}'.format(fname))
            idx0=((lxi0[:,None]-sx0[None,:])>=0).sum(axis=1)-1; ratx0=(lxi0-sx0[idx0])/(sx0[idx0+1]-sx0[idx0])
            idy0=((lyi0[:,None]-sy0[None,:])>=0).sum(axis=1)-1; raty0=(lyi0-sy0[idy0])/(sy0[idy0+1]-sy0[idy0])

            idx=((lxi[:,None]-sx0[None,:])>=0).sum(axis=1)-1; ratx=(lxi-sx0[idx])/(sx0[idx+1]-sx0[idx])
            idy=((lyi[:,None]-sy0[None,:])>=0).sum(axis=1)-1; raty=(lyi-sy0[idy])/(sy0[idy+1]-sy0[idy])
            idz=((lzi[:,None]-sz0[None,:])>=0).sum(axis=1)-1; ratz=(lzi-sz0[idz])/(sz0[idz+1]-sz0[idz])

        S.time.extend(ctime)
        for i, cti in enumerate(ctime):
            for k,svari in enumerate(svar):
                exec("cv=array(C.variables['{}'][{}])".format(svari,i)); mvari=mvar[k]

                #interp in space
                if mvari=='elev':
                    v0=array([cv[idy0,idx0],cv[idy0,idx0+1],cv[idy0+1,idx0],cv[idy0+1,idx0+1]])
                    #remove nan pts
                    for k in arange(4):
                        fpn=abs(v0[k])>1e3
                        v0[k,fpn]=sp.interpolate.griddata(bxy[~fpn,:],v0[k,~fpn],bxy[fpn,:],'nearest')

                    v1=v0[0]*(1-ratx0)+v0[1]*ratx0;  v2=v0[2]*(1-ratx0)+v0[3]*ratx0
                    vi=v1*(1-raty0)+v2*raty0

                else:
                    v0=array([cv[idz,idy,idx],cv[idz,idy,idx+1],cv[idz,idy+1,idx],cv[idz,idy+1,idx+1],
                              cv[idz+1,idy,idx],cv[idz+1,idy,idx+1],cv[idz+1,idy+1,idx],cv[idz+1,idy+1,idx+1]])

                    #remove nan pts
                    for k in arange(8):
                        fpn=abs(v0[k])>1e3
                        v0[k,fpn]=sp.interpolate.griddata(bxyz[~fpn,:],v0[k,~fpn],bxyz[fpn,:],'nearest',rescale=True)

                    v11=v0[0]*(1-ratx)+v0[1]*ratx;  v12=v0[2]*(1-ratx)+v0[3]*ratx; v1=v11*(1-raty)+v12*raty
                    v21=v0[4]*(1-ratx)+v0[5]*ratx;  v22=v0[6]*(1-ratx)+v0[7]*ratx; v2=v21*(1-raty)+v22*raty
                    vi=v1*(1-ratz)+v2*ratz

                #save data
                exec('S.{}.append(vi)'.format(mvari))
        C.close();
    S.time=array(S.time); [exec('S.{}=array(S.{})'.format(i,i)) for i in mvar]

    #interp in time
    for mvari in mvar:
        exec('vi=S.{}'.format(mvari))
        svi=interpolate.interp1d(S.time,vi,axis=0)(mtime)
        exec('S.{}=svi'.format(mvari))
    S.time=mtime

    #reshape the data, and save
    [exec('S.{}=S.{}.reshape([{},{},{}])'.format(i,i,nt,nobn,nvrt)) for i in mvar if i!='elev']
    #savez('hycom_{}'.format(mvar[0]),S) if len(mvar)==1 else savez('hycom_uv',S)

    #--------------------------------------------------------------------------
    #create netcdf
    #--------------------------------------------------------------------------
    nd=zdata(); nd.file_format='NETCDF4'

    #define dimensions
    nd.dimname=['nOpenBndNodes', 'nLevels', 'nComponents', 'one', 'time']
    if sname=='elev2D.th.nc':
        snvrt=1; ivs=1; vi=S.elev[...,None,None]
    elif sname=='uv3D.th.nc':
        snvrt=nvrt; ivs=2; vi=c_[S.u[...,None],S.v[...,None]]
    elif sname in ['TEM_3D.th.nc','SAL_3D.th.nc']:
        snvrt=nvrt; ivs=1; exec('vi=S.{}[...,None]'.format(mvar[0]))
    nd.dims=[nobn,snvrt,ivs,1,nt]

    # print(mvar,nd.dims,vi.shape)

    #--time step, time, and time series----
    nd.vars=['time_step', 'time', 'time_series']
    nd.time_step=zdata()
    nd.time_step.attrs=['long_name'];nd.time_step.long_name='time step in seconds'
    nd.time_step.dimname=('one',); nd.time_step.val=array(dt*86400).astype('float32')

    nd.time=zdata()
    nd.time.attrs=['long_name'];nd.time.long_name='simulation time in seconds'
    nd.time.dimname=('time',); nd.time.val=(S.time-S.time[0])*86400

    nd.time_series=zdata()
    nd.time_series.attrs=[]
    nd.time_series.dimname=('time','nOpenBndNodes','nLevels','nComponents')
    nd.time_series.val=vi

    WriteNC(sname,nd)

