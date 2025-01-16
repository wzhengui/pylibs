#!/usr/bin/env python3
#link sflux files
from pylib import *

#--------------------------------------------------------------------
#input
#--------------------------------------------------------------------
StartT=datenum(2009,2,10); EndT=datenum(2011,1,1)
tdir='sflux'                         #target dir
itag=1                               #itag=[1 or 2],for sflux_air_itag.0691.nc 
sdir='/sciclone/data10/wangzg/narr'  #narr source                

#--------------------------------------------------------------------
#make links
#--------------------------------------------------------------------
bdir=os.path.abspath(os.path.curdir); tdir=os.path.abspath(tdir)
if fexist(tdir): os.system('rm -rf {}'.format(tdir))
os.mkdir(tdir); os.chdir(tdir)
mtime=arange(StartT-2,EndT+2); svars=['air','prc','rad']
for irec,ti in enumerate(mtime):
    #link each file
    year=num2date(ti).year; month=num2date(ti).month; day=num2date(ti).day
    for m,svar in enumerate(svars):
        fname='{}/{}_{:02}/narr_{}.{}_{:02d}_{:02d}.nc'.format(sdir,year,month,svar,year,month,day)
        os.symlink(os.path.relpath(fname),'sflux_{}_{}.{:04d}.nc'.format(svar,itag,irec+1))
        if m==1 and day==1: print('    sflux: {:04d}-{:02d}-{:02d}'.format(year,month,day))
#write sflux_inputs.txt
fid=open('{}/sflux_inputs.txt'.format(tdir),'w+'); fid.write('&sflux_inputs\n \n/'); fid.close()
os.chdir(bdir)
