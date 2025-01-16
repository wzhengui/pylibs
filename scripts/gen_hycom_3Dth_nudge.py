#!/usr/bin/env python3
#create 3D boundary condition/nudging files based on hycom data
from pylib import *
close("all")

#------------------------------------------------------------------------------
#input
#------------------------------------------------------------------------------
StartT=datenum(2018,1,1); EndT=datenum(2021,1,1)
grd='../grid.npz'
dir_hycom='../../../HYCOM/Data'

#parameters to each files
iflags=[1,0,0,1,1,1]                    #if iflag=0: skip generating file 
dts=[1/24, 1.0, 1.0, 1/24, 1.0, 1.0]    #time steps for each file (day)
iLP=[1,0,0,1,0,0];   fc=0.25            #iLP=1: remove tidal signal with cutoff frequency fc (day)
snames=['elev2D.th.nc','TEM_3D.th.nc','SAL_3D.th.nc','uv3D.th.nc','TEM_nu.nc','SAL_nu.nc']
svars=['surf_el','water_temp','salinity',['water_u','water_v'],'water_temp','salinity']
mvars=['elev','temp','salt',['u','v'],'temp','salt']
ifix=0  #ifix=0: fix hycom nan 1st, then interp;  ifix=1: interp 1st, then fixed nan
#------------------------------------------------------------------------------
#interpolate hycom data to boundary
#------------------------------------------------------------------------------
#find all hycom files
fnames=array([i for i in os.listdir(dir_hycom) if i.endswith('.nc')])
mti=array([datenum(*array(i.replace('.','_').split('_')[1:5]).astype('int')) for i in fnames])
fpt=(mti>=(StartT-1))*(mti<(EndT+1)); fnames=fnames[fpt]; mti=mti[fpt]
sind=argsort(mti); mti=mti[sind]; fnames=fnames[sind]

#read hgrid
gd=loadz(grd).hgrid; vd=loadz(grd).vgrid; gd.x,gd.y=gd.lon,gd.lat; nvrt=vd.nvrt

#for each variables
for n,[sname,svar,mvar,dt,iflag] in enumerate(zip(snames,svars,mvars,dts,iflags)):
    if isinstance(svar,str): svar=[svar]; mvar=[mvar]
    if iflag==0: continue

    #get bnd or nugding nodes
    if sname.endswith('_nu.nc'):
       if sname.startswith('TEM'): gdn=read_schism_hgrid('TEM_nudge.gr3')
       if sname.startswith('SAL'): gdn=read_schism_hgrid('SAL_nudge.gr3')
       bind=nonzero(gdn.dp!=0)[0]; nobn=len(bind)
    else:
       bind=gd.iobn[0]; nobn=gd.nobn[0]

    #compute node xyz
    lxi0=gd.x[bind]%360; lyi0=gd.y[bind]; bxy=c_[lxi0,lyi0] #for 2D
    lxi=tile(lxi0,[nvrt,1]).T.ravel(); lyi=tile(lyi0,[nvrt,1]).T.ravel() #for 3D
    if vd.ivcor==2:
        lzi=abs(compute_zcor(vd.sigma,gd.dp[bind],ivcor=2,vd=vd)).ravel()
    else:
        lzi=abs(compute_zcor(vd.sigma[bind],gd.dp[bind])).ravel();
    bxyz=c_[lxi,lyi,lzi]

    #interp in space
    S=zdata(); sdict=S.__dict__
    for i in ['time',*mvar]: sdict[i]=[]
    sx0,sy0,sz0=None,None,None #used for check whether HYCOM files have the same dimensions
    for m,fname in enumerate(fnames):
        C=ReadNC('{}/{}'.format(dir_hycom,fname),1); print(fname)
        ctime=array(C.variables['time'])/24+datestr2num(C.variables['time'].time_origin)
        sx=array(C.variables['lon'][:])%360; sy=array(C.variables['lat'][:]); sz=array(C.variables['depth'][:]); nz=len(sz)
        fpz=lzi>=sz.max(); lzi[fpz]=sz.max()-1e-6

        if not array_equal(sx,sx0)*array_equal(sy,sy0)*array_equal(sz,sz0):
            #get interp index for HYCOM data
            if ifix==0:
                sxi,syi=meshgrid(sx,sy); sxy=c_[sxi.ravel(),syi.ravel()];
                cvs=array(C.variables['water_temp'][0]); sindns=[]; sindps=[]
                for ii in arange(nz):
                    print('computing HYCOM interpation index: level={}/{}'.format(ii,nz))
                    cv=cvs[ii]; ds=cv.shape; cv=cv.ravel()
                    fpn=abs(cv)>1e3; sindn=nonzero(fpn)[0]; sindr=nonzero(~fpn)[0]
                    if len(sindr)!=0: 
                        sindp=sindr[near_pts(sxy[sindn],sxy[sindr])]
                    else:
                        sindp=array([])
                    sindns.append(sindn); sindps.append(sindp)

            #get interp index for pts
            sx0=sx[:]; sy0=sy[:]; sz0=sz[:]; print('get new interp indices: {}'.format(fname))
            idx0=((lxi0[:,None]-sx0[None,:])>=0).sum(axis=1)-1; ratx0=(lxi0-sx0[idx0])/(sx0[idx0+1]-sx0[idx0])
            idy0=((lyi0[:,None]-sy0[None,:])>=0).sum(axis=1)-1; raty0=(lyi0-sy0[idy0])/(sy0[idy0+1]-sy0[idy0])

            idx=((lxi[:,None]-sx0[None,:])>=0).sum(axis=1)-1; ratx=(lxi-sx0[idx])/(sx0[idx+1]-sx0[idx])
            idy=((lyi[:,None]-sy0[None,:])>=0).sum(axis=1)-1; raty=(lyi-sy0[idy])/(sy0[idy+1]-sy0[idy])
            idz=((lzi[:,None]-sz0[None,:])>=0).sum(axis=1)-1; ratz=(lzi-sz0[idz])/(sz0[idz+1]-sz0[idz])

        S.time.extend(ctime)
        for i, cti in enumerate(ctime):
            for k,[svari,mvari] in enumerate(zip(svar,mvar)):
                cv=array(C.variables[svari][i])
                if sum(abs(cv)<1e3)==0: sdict[mvari].append(sdict[mvari][-1]); continue #fix nan data at this time
                #interp in space
                if mvari=='elev':
                    #remove HYCOM nan pts
                    if ifix==0:
                        sindn,sindp=sindns[0],sindps[0]
                        cv=cv.ravel(); fpn=(abs(cv[sindn])>1e3)*(abs(cv[sindp])<1e3); cv[sindn]=cv[sindp]; fpn=abs(cv)>1e3 #init fix
                        if sum(fpn)!=0: fni=nonzero(fpn)[0]; fri=nonzero(~fpn)[0]; fpi=fri[near_pts(sxy[fni],sxy[fri])]; cv[fni]=cv[fpi] #final fix
                        cv=cv.reshape(ds)

                    v0=array([cv[idy0,idx0],cv[idy0,idx0+1],cv[idy0+1,idx0],cv[idy0+1,idx0+1]]) #find parent pts
                    for ii in arange(4): #remove nan in parent pts
                        if ifix==1: fpn=abs(v0[ii])>1e3; v0[ii,fpn]=sp.interpolate.griddata(bxy[~fpn,:],v0[ii,~fpn],bxy[fpn,:],'nearest')
                    v1=v0[0]*(1-ratx0)+v0[1]*ratx0; v2=v0[2]*(1-ratx0)+v0[3]*ratx0; vi=v1*(1-raty0)+v2*raty0 #interp
                else:
                    #remove HYCOM nan pts
                    if ifix==0:
                        for ii in arange(nz):
                            if sum(abs(cv[ii])<1e3)==0: cv[ii]=cv[ii-1] #fix nan data for whole level
                            sindn,sindp=sindns[ii],sindps[ii]
                            if len(sindp)!=0:
                               cvi=cv[ii].ravel(); fpn=(abs(cvi[sindn])>1e3)*(abs(cvi[sindp])<1e3); cvi[sindn]=cvi[sindp]; fpn=abs(cvi)>1e3 #init fix
                               if sum(fpn)!=0: fni=nonzero(fpn)[0]; fri=nonzero(~fpn)[0]; fpi=fri[near_pts(sxy[fni],sxy[fri])]; cvi[fni]=cvi[fpi] #final fix

                    v0=array([cv[idz,idy,idx],cv[idz,idy,idx+1],cv[idz,idy+1,idx],cv[idz,idy+1,idx+1],
                              cv[idz+1,idy,idx],cv[idz+1,idy,idx+1],cv[idz+1,idy+1,idx],cv[idz+1,idy+1,idx+1]]) #find parent pts
                    for ii in arange(8): #remove nan in parent pts
                        if ifix==1: fpn=abs(v0[ii])>1e3; v0[ii,fpn]=sp.interpolate.griddata(bxyz[~fpn,:],v0[ii,~fpn],bxyz[fpn,:],'nearest',rescale=True)
                    v11=v0[0]*(1-ratx)+v0[1]*ratx;  v12=v0[2]*(1-ratx)+v0[3]*ratx; v1=v11*(1-raty)+v12*raty
                    v21=v0[4]*(1-ratx)+v0[5]*ratx;  v22=v0[6]*(1-ratx)+v0[7]*ratx; v2=v21*(1-raty)+v22*raty
                    vi=v1*(1-ratz)+v2*ratz  #interp in space
                sdict[mvari].append(vi) #save data
        C.close();
    for i in ['time',*mvar]: sdict[i]=array(sdict[i])

    #interp in time
    mtime=arange(StartT,EndT+dt,dt); nt=len(mtime)
    for mvari in mvar:
        svi=interpolate.interp1d(S.time,sdict[mvari],axis=0)(mtime)
        if iLP[n]==1: svi=lpfilt(svi,dt,fc) #low-pass
        sdict[mvari]=svi
    S.time=mtime
    for i in setdiff1d(mvar,'elev'): sdict[i]=sdict[i].reshape([nt,nobn,nvrt]) #reshape the data 

    #--------------------------------------------------------------------------
    #create netcdf
    #--------------------------------------------------------------------------
    if sname.endswith('.th.nc'):
       #define dimensions
       dimname=['nOpenBndNodes', 'nLevels', 'nComponents', 'one', 'time']
       if sname=='elev2D.th.nc':
           dims=[nobn,1,1,1,nt]; vi=S.elev[...,None,None]
       elif sname=='uv3D.th.nc':
           dims=[nobn,nvrt,2,1,nt]; vi=c_[S.u[...,None],S.v[...,None]]
       elif sname in ['TEM_3D.th.nc','SAL_3D.th.nc']:
           dims=[nobn,nvrt,1,1,nt]; vi=sdict[mvar[0]][...,None]
       nd=zdata(); nd.dimname=dimname; nd.dims=dims

       #define variables
       z=zdata(); z.attrs=['long_name']; z.long_name='time step (sec)'; z.dimname=('one',); z.val=array(dt*86400); nd.time_step=z
       z=zdata(); z.attrs=['long_name']; z.long_name='time (sec)'; z.dimname=('time',); z.val=(S.time-S.time[0])*86400; nd.time=z
       z=zdata(); z.dimname=('time','nOpenBndNodes','nLevels','nComponents'); z.val=vi.astype('float32'); nd.time_series=z
    else:
       #define dimensions
       dimname=['time','node','nLevels','one']
       dims=[nt,nobn,nvrt,1]; vi=sdict[mvar[0]][...,None]
       nd=zdata(); nd.dimname=dimname; nd.dims=dims

       #nd.vars=['time', 'map_to_global_node', 'tracer_concentration']
       z=zdata(); z.dimname=('time',); z.val=(S.time-S.time[0])*86400; nd.time=z
       z=zdata(); z.dimname=('node',); z.val=bind+1; nd.map_to_global_node=z
       z=zdata(); z.dimname=('time','node','nLevels','one'); z.val=vi.astype('float32'); nd.tracer_concentration=z

    WriteNC(sname,nd)
