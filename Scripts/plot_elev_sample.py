#!/usr/bin/env python3
from pylib import *
close("all")

#input
StartT=datenum(2018,8,24) #plot start time
MT1=datenum(2018,6,17)+1; MT2=MT1+118
# runs=['elev.Coast_0.RUN1_ZG','elev.Coast_0.RUN1c_ZG']
runs=['elev.Coast_0.RUN1c_ZG','elev.Coast_1.RUN1c_ZG']
# runs=['elev.Coast_0.RUN1_ZG',]

#read station.bp
fid=open('stations.txt','r'); lines=fid.readlines(); fid.close()
stations=array([i.strip().split(',')[0] for i in lines ])
stanames=array([i.strip().split(',')[1] for i in lines ])
ss=dict(zip(stations,stanames))
bp=read_schism_bpfile('Coast_0.bp')
bp.staname=array([ss[i] for i in bp.station])

#read model results
Model=[]
for m, run in enumerate(runs):
    S=loadtxt(run); mti=S[:,0]+datenum(2018,6,17); S=S[:,1:]
    Model.append(S)
    # S1=loadtxt(run); mti=S1[:,0]+datenum(2018,6,17); S1=S1[:,1:]

#read obs data
C=loadz('./noaa_elev_navd.npz')
C1=loadz('./noaa_elev_msl_2010_2020.npz')

#-----------------------plot---------------------------------------------------
MAE=[]; Corr=[]; MSL_shift=[]
for m in arange(5):
    i1=m*16; i2=(m+1)*16; iflag=0
    figure(figsize=[18,9.4])
    for i in arange(i1,i2):
        iflag=iflag+1;
        station=int(bp.station[i])

        subplot(8,2,iflag)

        #get obs
        # fp=(C.station==station)*(C.time>=(StartT-120))*(C.time<=(StartT+120))
        fp=(C.station==station)*(C.time>=MT1)*(C.time<=MT2); oti=C.time[fp]; oyi=C.elev[fp]
        if sum(fp)==0:
            fp=(C1.station==station)*(C1.time>=MT1)*(C1.time<=MT2); oti=C1.time[fp]; oyi=C1.elev[fp]
        plot(oti-StartT,oyi,'r')

        #plot
        colors='bg'
        for n, run in enumerate(runs):
            myi=Model[n][:,i]
            plot(mti-StartT,myi,color=colors[n])
        # setp(gca(),xticks=[*arange(0,120,10)],xlim=[-70,50])
        setp(gca(),xticks=[*arange(0,120,10)],xlim=[10,15])
        if iflag<=14: setp(gca(),xticklabels=[])
        gca().xaxis.grid('on')
        text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'{}: {}'.format(bp.station[i],bp.staname[i]))
        if iflag==1: legend(['Obs','RUN1'])

        #do statistics

        if sum(oti)!=0:
            MAEi=[]; Corri=[]; MSL_shifti=[]
            for n, run in enumerate(runs):
                myi=Model[n][:,i]
                myii=interpolate.interp1d(mti,myi)(oti)
                st=get_stat(myii,oyi); msl_shift=mean(myii)-mean(oyi)
                if n==0: text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.6*diff(ylim()),'MAE={:6.3f}, Corr={:6.3f},MSL_shift={:6.3f}'.format(st.MAE,st.R,msl_shift))
                MAEi.append(st.MAE); Corri.append(st.R); MSL_shifti.append(msl_shift)
            MAE.append(MAEi); Corr.append(Corri); MSL_shift.append(MSL_shifti)
    gcf().tight_layout()

MAE=array(MAE); Corr=array(Corr); MSL_shift=array(MSL_shift)

print('MAE: ', MAE.mean(axis=0))
print('Corr: ', Corr.mean(axis=0))
print('MSL_shift: ', MSL_shift.mean(axis=0))
# savefig('elev_cmp_RUN1_RUN1b_{}'.format(m))


















#------------------------------------------------------------------------------
# plot old
#------------------------------------------------------------------------------

# #---plot run1 and obs------------------------------------------

# StartT=datenum(2018,6,17)

# #read sta
# fid=open('stations.txt','r'); lines=fid.readlines(); fid.close()
# stations=array([i.strip().split(',')[0] for i in lines ])
# stanames=array([i.strip().split(',')[1] for i in lines ])
# ss=dict(zip(stations,stanames))
# bp=read_schism_bpfile('Coast_0.bp')
# bp.staname=array([ss[i] for i in bp.station])


# S1=loadtxt('elev.Coast_0.RUN1_ZG'); t1=S1[:,0]+StartT-datenum(2018,8,24); S1=S1[:,1:]
# C=loadz('./noaa_elev_navd.npz')
# # C=loadz('./noaa_elev_msl_2010_2020.npz')

# for m in arange(5):
#     i1=m*16; i2=(m+1)*16; iflag=0
#     figure(figsize=[18,9.4])
#     for i in arange(i1,i2):
#         iflag=iflag+1;
#         station=int(bp.station[i])

#         subplot(8,2,iflag)
#         y1=S1[:,i]; y2=S2[:,i]
#         #get model restls
#         fp=(C.station==station)*(C.time>=StartT)*(C.time<=(StartT+120))
#         oti=C.time[fp]-datenum(2018,8,24); oyi=C.elev[fp]
#         # plot(t1,y1,'b',t1,y2,'g')
#         plot(oti,oyi,'r',t1,y1,'b')
#         # setp(gca(),xticks=[20,21,22,23,24,25],xlim=[20,25])
#         setp(gca(),xticks=[*arange(10,20)],xlim=[10,15])
#         if iflag<=14: setp(gca(),xticklabels=[])
#         gca().xaxis.grid('on')
#         text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'{}: {}'.format(bp.station[i],bp.staname[i]))
#         if iflag==1: legend(['RUN1','RUN1b'])

#     gcf().tight_layout()
#     # savefig('elev_cmp_RUN1_RUN1b_{}'.format(m))

#------------------------------------------------------------------------------
#comp run1 and run1b
#------------------------------------------------------------------------------

# #read sta
# fid=open('stations.txt','r'); lines=fid.readlines(); fid.close()
# stations=array([i.strip().split(',')[0] for i in lines ])
# stanames=array([i.strip().split(',')[1] for i in lines ])
# ss=dict(zip(stations,stanames))
# bp=read_schism_bpfile('Coast_0.bp')
# bp.staname=array([ss[i] for i in bp.station])

# S1=loadtxt('elev.Coast_0.RUN1_ZG'); t1=S[:,0]; S1=S1[:,1:]
# S2=loadtxt('elev.Coast_0.RUN1b_ZG'); t2=S[:,0]; S2=S2[:,1:]


# for m in arange(5):
#     i1=m*16; i2=(m+1)*16; iflag=0
#     figure(figsize=[18,9.4])
#     for i in arange(i1,i2):
#         iflag=iflag+1;

#         subplot(8,2,iflag)
#         y1=S1[:,i]; y2=S2[:,i]
#         plot(t1,y1,'b',t1,y2,'g')
#         setp(gca(),xticks=[20,21,22,23,24,25],xlim=[20,25])
#         if iflag<=14: setp(gca(),xticklabels=[])
#         gca().xaxis.grid('on')
#         text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'{}: {}'.format(bp.station[i],bp.staname[i]))
#         if iflag==1: legend(['RUN1','RUN1b'])

#     gcf().tight_layout()
#     savefig('elev_cmp_RUN1_RUN1b_{}'.format(m))