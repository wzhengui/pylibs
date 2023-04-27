#!/usr/bin/env python3
# Dowload HF Radar-information--> https://secoora.org/high-frequency-radar/
from pylib import *
import pandas as pd

###############################################
StartT, EndT= '2015-1-1', '2015-1-2'
sdir='./' #location to save
rname='USEGC'# AKNS: Alaska North Slope, GAK: Gulf of Alaska, USHI: Hawaiian Islands, PRVI: Puerto Rico/Virgin Islands, USWC: U.S. West Coast, USEGC: U.S. East Coast and Gulf of Mexico
url='https://www.ncei.noaa.gov/data/oceans/ndbc/hfradar/rtv'
#######################################

tlist=pd.date_range(StartT, EndT, freq='MS').strftime("%Y-%m-%d")
tlist=datenum(tlist)
#download CMEMS data
for ti in tlist:
    furl='wget -N -r --no-parent -nH --reject "index.html*" -P {} {}/{}/{}/{}/'.format(sdir,url,num2date(ti).strftime('%Y'),num2date(ti).strftime('%Y%m'),rname)
    os.system(furl)
