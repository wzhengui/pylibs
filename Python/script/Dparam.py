#!/usr/bin/env python3
#show difference between SCHISM param.in files
# usage: Dparam.py param.in_1 param.in_2
#        Dparam.py -d param.in_1 param.in_2 (only show common varibles all files have)

#import os, sys
#from read_schism_file import read_schism_param
from pylib import *

flagd=0;
if len(sys.argv)==1:
   files=['param.in']
else:
   files=[];
   for vi in sys.argv[1:]:
      if vi[0]!='-':
         files.append(vi)
      elif vi=='-d':
         flagd=1

Par=[]; Key=[]; Val=[];
for i in range(len(files)):
   Pi=read_schism_param(files[i]);
   Par.append(Pi)    
   Key.append(Pi.keys())    
   Val.append(Pi.values())    
   if i==0:
      AKey=set(Key[i])
   else:
      [AKey.add(ki) for ki in Key[i]]

print("{:20s}:    ".format('Parameters')+', '.join("{:10s}".format(vi) for vi in files))
AKey=sorted(list(AKey)); 
for ki in AKey:
    v=[p[ki] if (ki in k) else 'N/A' for k, p in zip(Key,Par)]  
    if len(set(v))!=1:
       if flagd==1 and 'N/A' in v:
          pass
       else:
          print("{:20s}:    ".format(ki)+', '.join("{:10s}".format(vi) for vi in v))
