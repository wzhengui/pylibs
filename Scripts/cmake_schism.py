#!/usr/bin/env python3
'''
  auto compile schism model; work for cmake on sciclone/frontera
'''
from pylib import *

#-----------------------------------------------------------------------------------------------
#Inputs: 
#mods: OLDIO,PREC_EVAP,GOTM,HA,MARSH,SED2D,WWW,ICM,MICM,GEN,AGE,ECO,ICM,PH,COSINE,FIB,SED,FABM,ANALYSIS
#-----------------------------------------------------------------------------------------------
modules=['OLDIO','PREC_EVAP','AGE']

#directory of schism/fabm code
schism='~/schism'; fabm='~/fabm'

#-----------------------------------------------------------------------------------------------
#compile the code
#-----------------------------------------------------------------------------------------------
#write SCHSI.local.build
schism=schism.replace('~',os.path.expanduser("~")); fabm=fabm.replace('~',os.path.expanduser("~")); 
fname='{}/cmake/SCHISM.local.build'.format(schism)
rewrite(fname,replace=['ON','OFF'],include=['USE_','PREC_EVAP','OLDIO'])
rewrite(fname,replace=['OFF','ON'],include=[i if i in ['OLDIO','PREC_EVAP'] else 'USE_{}'.format(i) for i in modules])
rewrite(fname,replace=[],include=['FABM_BASE'])
rewrite(fname,append=['set( FABM_BASE {} CACHE STRING "Path to FABM base")\n'.format(fabm)])

#determine host
host=os.getenv('HOST').split('.')[0]
if host=='hurricane': host='whirlwind'

#compile 
if not fexist('{}/build'.format(schism)): os.mkdir('{}/build'.format(schism))
os.system('cd {}/build; rm -rf *; cmake -C ../cmake/SCHISM.local.build -C ../cmake/SCHISM.local.{} ../src; make -j8 pschism'.format(schism,host))

#put tag number
sname=os.listdir('{}/build/bin'.format(schism))[0]
irev=command_outputs('cd {}; git log'.format(schism)).stdout.split('\n')[0].split()[1][:8]
os.system('cp {}/build/bin/{} ./{}.{}'.format(schism,sname,sname,irev)); 
