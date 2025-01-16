#!/usr/bin/env python3
# generate flow.th using USGS river discharge (https://waterdata.usgs.gov/nwis/rt)

from pylib import *
import pandas as pd
import numpy as np

usgs=['USGS_02231254_RD.csv'] # USGS files (example --> https://nwis.waterservices.usgs.gov/nwis/iv/?sites=02231254&parameterCd=00060&startDT=2016-09-01T00:00:00.000-04:00&endDT=2016-11-01T23:59:59.999-04:00&siteStatus=all&format=rdb)
st=datenum(2016,9,8) # start time
et=datenum(2016,10,25) # end time
sname='flux.th' # save name
dt=900 # time step of flow.th (second)
pt=1 # check result 1:on


#generate flux.th
ntime=arange(0,(et-st)*86400,dt) # new time window
newset=ntime.copy() # new matrix to save data
for file in usgs:
    df = pd.read_csv(file,skiprows=26)
    df = df.drop([0])
    if df['tz_cd'][1]=='EDT': print('{} has EDT'.format(file));df['tz_cd'][1]='Etc/GMT+4'
    time=pd.to_datetime(df['datetime']).dt.tz_localize(df['tz_cd'][1]).dt.tz_convert('GMT')
    time= datenum(time.values.astype('str')).astype('float')
    rd=df[df.columns[4]].values.astype('float')*0.0283168

    #subset of time and data
    fpt=(time>=st)*(time<=et); time=time[fpt]; rd=rd[fpt]

    #interpolate the data to new time window and add it into new matrix
    time=(time-st)*86400; time,idx=unique(time,return_index=True); rd=rd[idx]
    nrd = -interpolate.interp1d(time, rd)(ntime)
    newset=column_stack((newset,nrd))


# save result
np.savetxt('{}'.format(sname),newset,fmt='%f')

#check result
if pt == 1:
    fs=loadtxt(sname)
    for nn in arange(shape(fs)[1]-1):
        plot(fs[:,0],fs[:,nn+1])
    xlabel('time (s)'); ylabel('River discharge (m^3/s)')
    show()
