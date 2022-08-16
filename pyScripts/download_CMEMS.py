#!/usr/bin/env python3
from pylib import *

# Need to install motuclinet (Web Server) to use this script (https://github.com/clstoulouse/motu-client-python.git)
# Check user manual for detailed information of CMEMS data 
# analysis/forcast: https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-024.pdf
# reanlaysis: https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-GLO-PUM-001-030.pdf

###############################################
StartT,EndT=datenum(2016,9,14),datenum(2016,9,16) # time period to download
xl=[-100,-50] # longitude
yl=[5,60] #latitude
dl=[0.494, 5727.917] #depth
myid='kpark' #user id
pswd='KyungminPark0616' #user password
sdir='./' #location to save
serv='GLOBAL_MULTIYEAR_PHY_001_030-TDS' # server name; read CMEMS manual above for detail 
prod='cmems_mod_glo_phy_my_0.083_P1D-m' # product name; read CMEMS manual above for detail
###############################################

#download CMEMS data
for ti in arange(StartT,EndT):

    #generating command for Motu    
    furl='motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id {} --product-id {} --longitude-min {} --longitude-max {} --latitude-min {} --latitude-max {} --date-min "{}" --date-max "{}" --depth-min {} --depth-max {} --variable thetao --variable so --variable uo --variable vo --variable zos --out-dir {} --out-name cmems_{}.nc --user {} --pwd {}'.format(serv,prod,xl[0],xl[1],yl[0],yl[1],num2date(ti).strftime('%Y-%m-%d %H:%M:%S'),num2date(ti+1).strftime('%Y-%m-%d %H:%M:%S'),dl[0],dl[1],sdir,num2date(ti+0.5).strftime('%Y_%m_%d_%H'),myid,pswd) 
    #call motu to download CMEMS data
    if os.path.exists('cmems_{}.nc'.format(num2date(ti+0.5).strftime('%Y_%m_%d_%H'))): print('cmems_{}.nc file exist'.format(num2date(ti+0.5).strftime('%Y_%m_%d_%H')));continue
    os.system(furl)      

print('----------Done-----------')
