#!/usr/bin/env python3
#check runtime for each run

#-----input------------------------------------------------------------
runs=['run1a','run2f','run3','run4','run4a','run5a']

#--compute runtime----------------------------------------------------
import os,sys,re,subprocess,time,datetime
from get_schism_param import read_schism_param
for run in runs:
    fname='{}/mirror.out'.format(run);
    if not os.path.exists(fname): continue

    #find start time
    with open(fname) as fid:
        line=fid.readline();    
    
    R=re.findall('(\d+), (\d+).(\d+)',line)[0]
    yyyy=int(R[0][:4]); mm=int(R[0][4:6]); dd=int(R[0][6:]); 
    HH=int(R[1][:2]); MM=int(R[1][2:4]); SS=int(R[1][4:]); mSS=int(R[2])*1000
    t0=datetime.datetime(yyyy,mm,dd,HH,MM,SS,mSS)

    #get lines in the end
    code='tail -n 60 {}'.format(fname)
    p=subprocess.Popen(code,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=True)
    out,err=p.communicate()
    if out!=None: out=out.decode('utf-8')  
    lines=out.split('\n');
    if lines[-1]=='': lines.pop()  
    lines.reverse()
    
    #find end time
    line=lines[0]
    if 'Run completed successfully' in line:
        R=re.findall('(\d+), (\d+).(\d+)',line)[0]
        yyyy=int(R[0][:4]); mm=int(R[0][4:6]); dd=int(R[0][6:]); 
        HH=int(R[1][:2]); MM=int(R[1][2:4]); SS=int(R[1][4:]);mSS=int(R[2])*1000
        t1=datetime.datetime(yyyy,mm,dd,HH,MM,SS,mSS)
    else:
        t1=datetime.datetime.fromtimestamp(os.path.getmtime(fname))

    #find time step 
    P=read_schism_param('{}/param.in'.format(run)); dt=float(P['dt'])
    
    #find number of time step  completed
    nstep=None
    for line in lines:
        if re.match('TIME STEP=',line): 
            nstep=int(re.findall('(\d+); ',line)[0])
            break
    if nstep==None: continue

    ds=(t1-t0).total_seconds() 

    #print results
    print('{}: {:.2f} days finished in {:.1f} hours (or {:.0f} minutes); need {:.1f} hours for 365 days'.format(run,dt*nstep/86400,ds/3600,ds/60,ds*365*24/dt/nstep))

    
