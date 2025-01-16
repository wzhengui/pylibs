#!/usr/bin/env python3
#build sflux directory based on NARR
from pylib import *
close("all")

#---------------------------------------------------------------------------
#inputs
#---------------------------------------------------------------------------
StartT=datenum(2009,12,31); EndT=datenum(2011,1,3)
xm=[-76.4,-76.25]; ym=[37.25,37.35]

sdir='/sciclone/data10/wangzg/narr' #narr source on sciclone
tdir='sflux'   #target dir
itag=1   #itag=[1 or 2],for sflux_air_itag.0691.nc

#---------------------------------------------------------------------------
#make sflux files
#---------------------------------------------------------------------------
if not os.path.exists(tdir): os.mkdir(tdir)
fid=open('{}/sflux_inputs.txt'.format(tdir),'w+'); fid.write('&sflux_inputs\n \n/'); fid.close()
mtime=arange(StartT,EndT+1); svars=['air','prc','rad']; ix1s,ix2s,iy1s,iy2s=[],[],[],[]
for irec,ti in enumerate(mtime):
    pti=num2date(ti); year,month,day=pti.year,pti.month,pti.day

    #link each file
    for m,svar in enumerate(svars):
        #NARR database
        fname='{}/{}_{:02}/narr_{}.{}_{:02d}_{:02d}.nc'.format(sdir,year,month,svar,year,month,day)
        bname='{}/sflux_{}_{}.{:04d}.nc'.format(tdir,svar,itag,irec+1)
        if day==1 and m==0: print('reading: {}'.format(fname))

        #compute indices of subdomaine
        if irec==0:
           C=ReadNC(fname,1);lon=array(C.variables['lon'][:]); lat=array(C.variables['lat'][:]); C.close(); 
           ix1,ix2,iy1,iy2=subdomain_index(lon,lat,xm,ym) #save indices
           ix1s.append(ix1); ix2s.append(ix2); iy1s.append(iy1); iy2s.append(iy2)

        #read sflux, and change dimension
        C=ReadNC(fname);  ix1,ix2,iy1,iy2=ix1s[m],ix2s[m],iy1s[m],iy2s[m]
        for n,dimname in enumerate(C.dimname):
            if dimname=='nx_grid': C.dims[n]=ix2-ix1
            if dimname=='ny_grid': C.dims[n]=iy2-iy1

        #change variables values
        for mvar in C.vars:
            if mvar in ['time']: continue
            if mvar in ['lon','lat']:
               exec('C.{}.val=C.{}.val[iy1:iy2,ix1:ix2]'.format(mvar,mvar))
            else:
               exec('C.{}.val=C.{}.val[:,iy1:iy2,ix1:ix2]'.format(mvar,mvar))
        WriteNC(bname,C) #write sflux
