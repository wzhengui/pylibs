#!/usr/bin/env python3
from pylib import *
import datetime
# This script download 
# Need to install motuclinet (Web Server) to use this script (https://github.com/clstoulouse/motu-client-python.git)
# Check user manual for detailed information of CMEMS data 
# analysis/forcast: https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-024.pdf
# reanlaysis: https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-030.pdf

###############################################
StartT,EndT=datenum(2018,12,31,fmt=1),datenum(2020,1,1,fmt=1) # time 
tres='daily' # hourly, daily, monthly. 
xl=[-163,-146] # longitude
yl=[16,26] #latitude

dl=[0, 8727.917] #depth
data='rana' # analysis:'GLOBAL_ANALYSISFORECAST_PHY_001_024-TDS'; reanlysis: GLOBAL_MULTIYEAR_PHY_001_030-TDS
myid='kpark' #user id
pswd='KyungminPark0616' #user password
sdir='./' #location to save
svrs=['so','thetao','uv','zos'] # so: water salinity, thetao: water temperature, uv: water velocities, zos: sea surface height

#######################################
#download CMEMS data
if data=='ana':
   serv='GLOBAL_ANALYSISFORECAST_PHY_001_024-TDS'
   print('Downloading analysis data from CMEMS')
   for svr in svrs:
       # Variable definition
       if svr=='so': prod='-so'; vr='--variable so'
       elif svr=='thetao': prod='-thetao'; vr='--variable thetao'
       elif svr=='uv': prod='-cur'; vr='--variable uo --variable vo';
       elif svr=='zos': prod=''; vr='--variable zos'
       else: print('Unrecognized variable --{}--.'.format(svr)); break

       # Temporal resolution definition
       if tres=='hourly':
          if svr=='zos': dt=datetime.timedelta(hours=1); pid='PT1H-m'
          else: dt=datetime.timedelta(hours=6); pid='PT6H-i'
       elif tres=='daily':
          dt=datetime.timedelta(hours=24); pid='P1D-m'
       elif tres=='monthly':
          dt=datetime.timedelta(days=31); pid='P1M-m'
       else: print('Unrecognized temporal resolution --{}--.'.format(tres)); break
       tlist=drange(StartT,EndT+dt,dt)
       dt=tlist[1]-tlist[0]
       # Download data
       for nn,ti in enumerate(tlist):
           if tres=='hourly' or tres=='daily':
              st=num2date(ti).strftime('%Y_%m_%d_%H')
           elif tres=='monthly':
              st=num2date(ti).strftime('%Y_%m_00_00')
           if os.path.isfile('cmems_{}_{}.nc'.format(svr,st)): print('Variable --{}-- on {} exsit'.format(svr,st)); continue
           furl='motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id {} --product-id cmems_mod_glo_phy{}_anfc_0.083deg_{} --longitude-min {} --longitude-max {} --latitude-min {} --latitude-max {} --date-min "{}" --date-max "{}" --depth-min {} --depth-max {} {} --out-dir {} --out-name cmems_{}_{}.nc --user {} --pwd {}'.format(serv,prod,pid,xl[0],xl[1],yl[0],yl[1],num2date(ti).strftime('%Y-%m-%d %H:%M:%S'),num2date(ti+dt*0.9).strftime('%Y-%m-%d %H:%M:%S'),dl[0],dl[1],vr,sdir,svr,st,myid,pswd)
           #print(furl)
           os.system(furl)

elif data=='rana':
       serv='GLOBAL_MULTIYEAR_PHY_001_030-TDS'
       print('Downloading reanalysis data from CMEMS')
       # Temporal resolution definition
       if tres=='daily':
          dt=datetime.timedelta(hours=24); pid='P1D-m'
       elif tres=='monthly':
          dt=datetime.timedelta(days=31); pid='P1M-m'
       else: print('Unrecognized temporal resolution --{}--.'.format(tres)); 
       tlist=drange(StartT,EndT+dt,dt)
       dt=tlist[1]-tlist[0]
       # Download data
       for nn,ti in enumerate(tlist):
           if tres=='daily':
              st=num2date(ti).strftime('%Y_%m_%d_%H')
           elif tres=='monthly':
              st=num2date(ti).strftime('%Y_%m_00_00')
           if os.path.isfile('cmems_{}.nc'.format(st)): print('Data on {} exsit'.format(st)); continue
           furl='motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id cmems_mod_glo_phy_my_0.083_{} --longitude-min {} --longitude-max {} --latitude-min {} --latitude-max {} --date-min "{}" --date-max "{}" --depth-min {} --depth-max {} --variable so --variable thetao --variable uo --variable vo --variable zos --out-dir {} --out-name cmems_{}.nc --user {} --pwd {}'.format(serv,pid,xl[0],xl[1],yl[0],yl[1],num2date(ti).strftime('%Y-%m-%d %H:%M:%S'),num2date(ti+dt*0.9).strftime('%Y-%m-%d %H:%M:%S'),dl[0],dl[1],sdir,st,myid,pswd)
           print(furl)
           os.system(furl)

else: print('Wrong data type. Choose ana or rana in data')

print('---------------done------------')
