#!/usr/bin/env python3
from pylib import *

#compile
os.system('gcc tidal_analysis.c nrutil.c -lm; mv a.out tidal_analyze');

#sys.exit()
#generate tidal consitituents. 
#tide_name=['O1','K1','Q1','P1','M2','S2','K2','N2'];
#C=loadz('tide_fac_const.npz');
#with open('tidal_const.dat','w+') as fid:
#    fid.write('{}\n'.format(len(tide_name)))
#    for m in arange(len(tide_name)):
#        fp=C.name==tide_name[m]; freqi=squeeze(C.freq[fp])
#        fid.write('{}\n {:e}\n'.format(tide_name[m],freqi))
#
##analyze tidal components
#os.system('./tidal_analyze Time_Series.txt tidal_const.dat t1.txt');
#os.system('./tidal_analyze Time_Series.txt tidal_const.dat.sample t1.txt');
