#!/usr/bin/env python3
#check runtime for each run
from pylib import *
import subprocess

#-----input------------------------------------------------------------
runs=['RUN01a','RUN01b']

#--compute runtime----------------------------------------------------
for run in runs:
    fname='{}/mirror.out'.format(run)
    if not os.path.exists(fname): fname='{}/outputs/mirror.out'.format(run);
    if not os.path.exists(fname): continue

    #find start time
    with open(fname) as fid:
        line=fid.readline();    
    
    R=re.findall('(\d+), (\d+).(\d+)',line)[0]
    yyyy=int(R[0][:4]); mm=int(R[0][4:6]); dd=int(R[0][6:]); 
    HH=int(R[1][:2]); MM=int(R[1][2:4]); SS=int(R[1][4:]); mSS=int(R[2])*1000
    t0=datetime.datetime(yyyy,mm,dd,HH,MM,SS,mSS)

    #get the start time step
    fid=open(fname,'r'); nstep0=0
    while (True): 
        line=fid.readline()
        if len(line)==0: sys.exit()
        if line.strip().startswith('TIME STEP='): 
           nstep0=float(line.split(';')[0].split('=')[1])
           break
      
    fid.close()

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
    if os.path.exists('{}/param.in'.format(run)):
       P=read_schism_param('{}/param.in'.format(run)); dt=float(P['dt']); rnday=float(P['rnday'])
    else:
       P=read_schism_param('{}/param.nml'.format(run)); dt=float(P['dt']); rnday=float(P['rnday'])
    
    #find number of time step  completed
    nstep1=None
    for line in lines:
        if re.match('TIME STEP=',line): 
            nstep1=int(re.findall('(\d+); ',line)[0])
            break
    if nstep1==None: continue

    #output values
    nstep=nstep1-nstep0
    ds=(t1-t0).total_seconds()
    RTR=(dt*nstep/86400)/(ds/86400)

    nday=rnday; nday0=(dt*nstep0)/86400; nday1=(dt*nstep1)/86400 
    time_all=nday*24/RTR
    time_left=(nday-nday1)*24/RTR
    time_365=365*24/RTR

    #print results
    if nday0*86400<2*dt: 
       print('{}: RTR={:.1f}; {:.2f} days finished in {:.1f} hrs ({:.0f} min); ({:.1f}, {:.1f}, {:.1f}) hrs needed for ({:.1f}, {:.1f}, {:.0f}) days'.format \
             (run,RTR,nday1,ds/3600,ds/60, time_left,time_all,time_365,(nday-nday1),nday, 365))
    else:
       print('{}: RTR={:.1f}; {:.2f} ({:.1f}) days finished in {:.1f} hrs ({:.0f} min); ({:.1f}, {:.1f}, {:.1f}) hrs needed for ({:.1f}, {:.1f}, {:.0f}) days'.format \
             (run,RTR,(nday1-nday0),nday1,ds/3600,ds/60, time_left,time_all,time_365,(nday-nday1),nday, 365))

    
sys.exit()

    
