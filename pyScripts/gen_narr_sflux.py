#!/usr/bin/env python3
from pylib import *
import urllib

#------------------------------------------------------------------------------------
#inputs: creat narr sflux database
#------------------------------------------------------------------------------------
years=[2021,]

#ftp info
bdir='ftp://ftp.cdc.noaa.gov/NARR/monolevel'
vars=['air','rad','prc']
svars=[('uwind','vwind','prmsl','stmp','spfh'),('dlwrf','dswrf'),('prate',)] #sflux variable
nvars=[('uwnd.10m','vwnd.10m','prmsl','air.2m','shum.2m'),('dlwrf','dswrf'),('prate',)] #narr variables

#download data, and precess data
for year in years:
    #create folder
    sdir='{}'.format(year);
    if not os.path.exists(sdir): os.mkdir(sdir)

    #pre-calculation
    S0=loadz('sflux_template.npz');
    # days=arange(datenum(year,1,1),datenum(year+1,1,1))

    #for each dataset
    for m in arange(len(vars)):
        vari=vars[m]; svari=svars[m]; nvari=nvars[m]
        #download data
        for nvarii in nvari:
            fname='{}.{}.nc'.format(nvarii,year)
            if os.path.exists(fname): continue
            url='{}/{}'.format(bdir,fname)
            print('downloading {}'.format(fname))
            urllib.request.urlretrieve(url,fname)

#        #read the data
        C=zdata()
        for svarii,nvarii in zip(svari,nvari):
            fname='{}.{}.nc'.format(nvarii,year)
            exec('C.{}=ReadNC("{}")'.format(svarii,fname))

        #processing data and write sflux
        exec('S=S0.{}'.format(vari))
        exec('time=datenum(1800,1,1)+array(C.{}.time.val)/24'.format(svari[0]))
        exec('xgrid=array(C.{}.x.val); ygrid=array(C.{}.y.val)'.format(svari[0],svari[0]))
        exec('lon=array(C.{}.lon.val); lat=array(C.{}.lat.val)'.format(svari[0],svari[0]))

        #get days
        days=unique(time.astype('int'))
        nt,nx,ny=[int(len(time)/len(days)),len(xgrid),len(ygrid)]
        for dayi in days:
            ti=num2date(dayi)
            fp=(time>=dayi)*(time<(dayi+1));

            #dims
            S.dims=[nx,ny,nt]
            #time, lon, lat
            S.time.base_date=array([ti.year,ti.month,ti.day,0]);
            S.time.units='days since {}'.format(ti.strftime('%Y-%m-%d'))
            S.time.dims=[nt]; S.time.val=time[fp]-dayi;
            S.lon.dims=[ny,nx]; S.lon.val=lon
            S.lat.dims=[ny,nx]; S.lat.val=lat
            #variables
            for svarii,nvarii in zip(svari,nvari):
                exec('S.{}.dims=[nt,ny,nx]'.format(svarii));
                exec('S.{}.val=C.{}.{}.val[fp,:,:]'.format(svarii,svarii,nvarii.split('.')[0]));

            #write narr files
            fname='narr_{}.{}.nc'.format(vari,ti.strftime('%Y_%m_%d'))
            print('writing {}'.format(fname))
            WriteNC('{}/{}'.format(sdir,fname),S)

    #move files
    if sys.platform.startswith('win'): continue
    for i in arange(1,13):
        subdir='{}_{:02}'.format(year,i);
        if not os.path.exists('{}/{}'.format(sdir,subdir)): os.mkdir('{}/{}'.format(sdir,subdir))
        os.system('cd {}; mv *{}*.nc {}'.format(sdir,subdir,subdir))
    os.system("ln -sf {}/{}_* ./ ".format(sdir,year))

#------------------------------------------------------------------------------
##--prepare template for sflux based on former sflux files
#------------------------------------------------------------------------------
#S=zdata();
#svars=['air','rad','prc']
#for svar in svars:
#    fname='sflux_{}_1.0001.nc'.format(svar)
#    Si=ReadNC(fname,2);
#    #clear variables
#    for vari in Si.vars:
#        exec('Si.{}.val=None'.format(vari))
#    exec('S.{}=Si'.format(svar));
#S.vars=svars;
#savez('sflux_template',S);
