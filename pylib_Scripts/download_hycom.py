#!/usr/bin/env python3
#automatically download HYCOM data
from pylib import *
close("all")

#--------------------------------------------------------------------
#input
#--------------------------------------------------------------------
StartT,EndT=datenum(2009,12,31),datenum(2020,1,1)
xm=[-123.05,-122.35]; ym=[37.35,38.05]
sdir='Data'

#database
url_1='http://ncss.hycom.org/thredds/ncss'; url_2='https://ncss.hycom.org/thredds/ncss'
databases=[['GLBv0.08/expt_53.X',datestr2num('Jan-01-1994'),datestr2num('Dec-31-2015')],
           ['GLBv0.08/expt_56.3',datestr2num('Jan-01-2016'),datestr2num('Apr-30-2016')],
           ['GLBv0.08/expt_57.2',datestr2num('May-01-2016'),datestr2num('Jan-31-2017')],
           ['GLBv0.08/expt_92.8',datestr2num('Feb-01-2017'),datestr2num('May-31-2017')],
           ['GLBv0.08/expt_57.7',datestr2num('Jun-01-2017'),datestr2num('Sep-30-2017')],
           ['GLBv0.08/expt_92.9',datestr2num('Oct-01-2017'),datestr2num('Dec-31-2017')],
           ['GLBv0.08/expt_93.0',datestr2num('Jan-01-2018'),datestr2num('Feb-18-2020')],]

#download hycom data
if not os.path.exists(sdir): os.mkdir(sdir)
for ti in arange(StartT,EndT):
    #choose database
    bname=[i for i in databases if (ti>=i[1]) and (ti<=i[2])]
    bname=bname[0][0] if len(bname)==1 else sys.exit('database not found: {}'.format(num2date(ti)))

    #get fname
    if (ti>=datestr2num('Jan-01-1994'))*(ti<=datestr2num('Dec-31-2015')):
       url='{}/{}/data/{}?'.format(url_1,bname,num2date(ti).year) 
    else:
       url='{}/{}?'.format(url_2,bname) 

    #add variables
    for i in ['surf_el','salinity','water_temp','water_u','water_v']: url=url+'var={}&'.format(i)

    #add domain
    for i,k in zip(['south','north','west','east'],[*ym,*xm]): url=url+'{}={}&'.format(i,k)

    #horizontal stride and time
    for n in arange(0,24,3):
        fname='{}/hycom_{}_{:02}.nc'.format(sdir,num2date(ti).strftime('%Y_%m_%d'),n)
        furl=url+'horizStride=1&time={}T{:02}%3A00%3A00Z&vertCoord=&accept=netcdf4'.format(num2date(ti).strftime('%Y-%m-%d'),n)

        #download hycom data
        if os.path.exists(fname): continue
        try: 
           urllib.request.urlretrieve(furl,fname)
           print(fname)
        except: 
           pass
