#!/usr/bin/env python3
'''
  this script is convert source_sink.in, msource.th, vsource, vsink to source.nc
'''
from pylib import *
close("all")

#-----------------------------------------------------------
#inputs
#-----------------------------------------------------------
sdir='../RUN01b'   #directory of source_sink.in and *.th  

#-----------------------------------------------------------
#read and save source_sink information
#-----------------------------------------------------------
#save netcdf data
C=zdata(); C.vars=[]; C.file_format='NETCDF4'
nsources,nsinks,ntracers,ntm,ntv,nts=0,0,0,0,0,0; dtm,dtv,dts=1e6*ones(3)

#source_sink.in
lines=[i.strip() for i in open('{}/source_sink.in'.format(sdir),'r').readlines()]
nsources=int(lines[0].split()[0]); isource=array([int(i) for i in lines[1:(nsources+1)]])
nsinks=int(lines[2+nsources].split()[0]); isink=array([int(i) for i in lines[(nsources+3):(nsources+3+nsinks)]])

#vsource.th and msource.th
if nsources>0: 
   #read data
   fdata=loadtxt('{}/msource.th'.format(sdir)).T; mti=fdata[0]; msource=fdata[1:]
   fdata=loadtxt('{}/vsource.th'.format(sdir)).T; vti=fdata[0]; vsource=fdata[1:]
   ntm=len(mti); ntv=len(vti); dtm=floor(diff(mti)[0]); dtv=floor(diff(vti)[0]); ntracers=int(len(msource)/len(vsource))
   vsource=vsource.T;  msource=msource.reshape([ntracers,nsources,len(mti)]).transpose([2,0,1])
   
   #save data
   C.vars.extend(['source_elem','vsource','msource'])
   vi=zdata(); vi.dimname=('nsources',); vi.val=isource; C.source_elem=vi
   vi=zdata(); vi.dimname=('time_vsource','nsources'); vi.val=vsource; C.vsource=vi
   vi=zdata(); vi.dimname=('time_msource','ntracers','nsources'); vi.val=msource; C.msource=vi
   
   #not needed
   #C.vars.extend(['time_msource','time_vsource',])
   #vi=zdata(); vi.dimname=('time_msource',); vi.val=mti; C.time_msource=vi
   #vi=zdata(); vi.dimname=('time_vsource',); vi.val=vti; C.time_vsource=vi

#vsink.th
if nsinks>0: 
   #read data
   fdata=loadtxt('{}/vsink.th'.format(sdir)).T; sti=fdata[0]; vsink=fdata[1:].T
   nts=len(sti); dts=floor(diff(sti)[0])

   #save data
   C.vars.extend(['sink_elem','vsink'])
   vi=zdata(); vi.dimname=('nsinks',); vi.val=isink; C.sink_elem=vi
   vi=zdata(); vi.dimname=('time_vsink','nsinks',); vi.val=vsink; C.vsink=vi

   #not needed
   #C.vars.extend(['time_vsink',])
   #vi=zdata(); vi.dimname=('time_vsink',); vi.val=sti; C.time_vsink=vi

#assign dimension value
C.dimname=['nsources','nsinks','ntracers','time_msource','time_vsource','time_vsink','one']
C.dims=[nsources,nsinks,ntracers,ntm,ntv,nts,1]

#add time step 
C.vars.extend(['time_step_vsource','time_step_msource','time_step_vsink'])
vi=zdata(); vi.dimname=('one',); vi.val=dtm; C.time_step_msource=vi
vi=zdata(); vi.dimname=('one',); vi.val=dtv; C.time_step_vsource=vi
vi=zdata(); vi.dimname=('one',); vi.val=dts; C.time_step_vsink=vi

#save as netcdf
WriteNC('source.nc',C)
