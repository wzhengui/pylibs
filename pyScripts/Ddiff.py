#!/usr/bin/env python3
'''
 compare directories and show differences
'''
import os,sys
from numpy import sort,unique,array,zeros
import subprocess

#get dirnames
if len(sys.argv)<3: 
   sys.exit('usage: Ddiff.py dir1 dir2 ...')
else:
   sdirs=sys.argv[1:]
ndir=len(sdirs); slens=[max(len(i),10) for i in sdirs]

#get fnames for each directory
fnames=[]; [fnames.extend(os.listdir(i)) for i in sdirs]
fnames=sort(unique(array([i for i in fnames if not i.startswith('.')])))

#compare each file
print("Files below are different(0:doesn't exist; same number: files are the same; directory is skipped")
fstr="{:10s} "*len(sdirs); print(fstr.format(*sdirs))
for fname in fnames:
    #print(fname)
    #get real path of each files
    rnames=[]
    for n,sdir in enumerate(sdirs): 
        ffname='{}/{}'.format(sdir,fname)
        if os.path.exists(ffname):
           rnames.append(os.path.realpath(ffname))
        else:
           rnames.append(None)

    #compare each file
    snames=[]; fns=zeros(ndir).astype('int')
    for n,rname in enumerate(rnames): 
        if rname is None: continue 
        if len(snames)==0:  #1st existing file
           snames.append(rname); fns[n]=len(snames)
        else: #find a new file 
           fn=None
           #compare with previouse files
           for m,sname in enumerate(snames):
               if rname==sname:  fn=m+1; break
               if os.path.isdir(rname): continue
               #use diff to compare two actual files
               code='diff -q {} {}'.format(rname,sname)
               p=subprocess.Popen(code,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True) 
               stdout,stderr=p.communicate()
               if stdout!=None: stdout=stdout.decode('utf-8')
               if 'differ' not in stdout: fn=m+1; break
           
           #save file status 
           if fn is None:
              snames.append(rname); fns[n]=len(snames)
           else:
              fns[n]=fn 
   
    #print status
    fns=array(fns);  fstr=''
    if len(unique(fns))==1: continue  
    for i,[fn,slen] in enumerate(zip(fns,slens)):
        fstr=fstr+'    {}{}'.format(fn,' '*(slen-4))
        if i==(ndir-1): fstr=fstr+': {}'.format(fname)
    print(fstr)
