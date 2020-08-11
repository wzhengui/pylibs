#!/usr/bin/env python3
'''
  compare package difference for different python verions
  usage: DPackage pip3.5 pip3.7 (only shows missing pkgs) 
         DPackage -a pip3.5 pip3.7 (shows all information) 
'''
from pylib import *

#read inputs
pv=sys.argv[1:]

sflag=True
if pv[0]=='-a': 
   sflag=False
   pv=pv[1:]

#parse all the packages for each python version
S=[]; svars=[]
for pvi in pv:
    lines=command_outputs('{} list'.format(pvi)).stdout.split('\n')
    Si=dict(); 
    for line in lines:
        sline=line.split()

        #check whether it is a package line
        if len(sline)!=2: continue
        if sline[1].split('.')[0][0]=='(': 
           if not sline[1].split('.')[0][1].isdigit(): continue
        else:
           if not sline[1].split('.')[0].isdigit(): continue
        
        #add package
        Si[sline[0]]=sline[1]

    #save 
    S.append(Si)
    svars.extend(Si.keys())

#compare different python packages 
svars=unique(array(svars))
sformat='{:25} '+'{:25} '*len(pv)

print(sformat.format('Package:',*pv))
for svar in svars:
    vi=[]; sflag2=True
    for m in arange(len(pv)):
        if svar in S[m].keys():
           vi.append(S[m][svar])
        else:
           vi.append('')
           sflag2=False
    if sflag2*sflag: continue
    print(sformat.format(svar,*vi))
