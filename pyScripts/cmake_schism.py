#!/usr/bin/env python3
'''
  auto compile schism model; work for cmake on sciclone/frontera
'''
from pylib import *

#-----------------------------------------------------------------------------------------------
#Inputs: 
#mods: OLDIO,PREC_EVAP,GOTM,HA,MARSH,SED2D,WWW,ICM,MICM,GEN,AGE,ECO,ICM,PH,COSINE,FIB,SED,FABM,ANALYSIS
#-----------------------------------------------------------------------------------------------
modules=['OLDIO', 'ICM', 'PREC_EVAP'] 
#modules=['OLDIO','COSINE']

#directory of schism/fabm code
schism='~/schism'; fabm='~/fabm'

target=sys.argv[1] if len(sys.argv)>1 else 'pschism' #combine target
#-----------------------------------------------------------------------------------------------
#compile the code
#-----------------------------------------------------------------------------------------------
#write SCHSI.local.build
schism=schism.replace('~',os.path.expanduser("~")); fabm=fabm.replace('~',os.path.expanduser("~")); 
fname='{}/cmake/SCHISM.local.build'.format(schism) 

fid=open(fname,'r'); lines=fid.readlines(); fid.close() #save original file
rewrite(fname,replace=['ON','OFF'],include=['USE_','PREC_EVAP','OLDIO'])
rewrite(fname,replace=['OFF','ON'],include=[i if i in ['OLDIO','PREC_EVAP'] else 'USE_{}'.format(i) for i in modules])
rewrite(fname,replace=[],include=['FABM_BASE'])
rewrite(fname,append=['set( FABM_BASE {} CACHE STRING "Path to FABM base")\n'.format(fabm)])
if 'SED2D' not in modules: rewrite(fname,replace=['ON','OFF'],include=['SED2D'])

#determine host
host=os.getenv('HOST').split('.')[0]
if host=='hurricane': host='whirlwind'

try: 
   #compile 
   if not fexist('{}/build'.format(schism)): os.mkdir('{}/build'.format(schism))
   os.system('cd {}/build; rm -rf *; cmake -C ../cmake/SCHISM.local.build -C ../cmake/SCHISM.local.{} ../src; make -j8 {}'.format(schism,host,target))
   
   #put tag number
   sname=os.listdir('{}/build/bin'.format(schism))[0]
   irev=command_outputs('cd {}; git log'.format(schism)).stdout.split('\n')[0].split()[1][:8]
   os.system('cp {}/build/bin/{} ./{}.{}'.format(schism,sname,sname,irev)); 

   #write original file
   fid=open(fname,'w+'); fid.writelines(lines); fid.close() 
except: 
   #write original file
   fid=open(fname,'w+'); fid.writelines(lines); fid.close()
