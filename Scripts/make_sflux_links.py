#!/usr/bin/env python3
#link sflux files
from pylib import *

#input
#StartT=datenum(2009,2,10); EndT=datenum(2011,1,1)
StartT=datenum(2009,2,10); EndT=datenum(2009,3,5)
sdir='/ches/data10/yinglong/narr'  #source dir
tdir='sflux'                       #target dir
itag=1   #itag=[1 or 2],for sflux_air_itag.0691.nc 

#make links
if not os.path.exists(tdir): os.mkdir(tdir)
mtime=arange(StartT,EndT+1); svars=['air','prc','rad']
for irec,ti in enumerate(mtime):
    #get date 
    year=num2date(ti).year
    month=num2date(ti).month
    day=num2date(ti).day

    #link each file
    for m,svar in enumerate(svars):
        fname='{}/{}_{:02}/narr_{}.{}_{:02d}_{:02d}.nc'.format(sdir,year,month,svar,year,month,day)
        print(fname)
        os.system('cd {}; ln -sf {} sflux_{}_{}.{:04d}.nc'.format(tdir,fname,svar,itag,irec+1))

#write sflux_inputs.txt
fid=open('{}/sflux_inputs.txt'.format(tdir),'w+'); fid.write('&sflux_inputs\n/'); fid.close()
   

    

    

