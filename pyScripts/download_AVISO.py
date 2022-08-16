#!/usr/bin/env python3
from ftplib import FTP
from pylib import *

###############################################
StartT=datenum('2016-9-1'); EndT=datenum('2016-9-4') # time period to download
host='ftp-access.aviso.altimetry.fr' # domain name 
user='kmpark19900616@gmail.com'; passwd = 'Y0PLwL' # user id and password
wdir='/duacs-experimental/dt-phy-grids/multiscale_interpolation_alti_drifters/version_01_00'# target directory to download
#ftp.dir() #list contents in current dir
###############################################

# connect to host and go to target directory
ftp = FTP(host)
ftp.login(user=user, passwd = passwd)
ftp.cwd(wdir)#cd to target directory

# organazing and sorting file names
fnames = ftp.nlst(); fnames=array(fnames)
mti=datenum(array([(i.replace('.','_').split('_')[5]) for i in fnames]))
fpt=(mti>=(StartT))*(mti<(EndT)); fnames=fnames[fpt]; mti=mti[fpt]
sind=argsort(mti); mti=mti[sind]; fnames=fnames[sind]

# dowload data
for fname in fnames:
    print('Downloading {}'.format(fname))
    with open(fname, "wb") as file:
        ftp.retrbinary(f"RETR {fname}", file.write)
print('--------------Done-------------')
