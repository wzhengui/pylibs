#!/usr/bin/env python3
from pylib import *
close("all")
#------------------------------------------------------------------------------
#input
# ------------------------------------------------------------------------------
StartT=datenum(2018,7,13); StartT_model=datenum(2018,6,17) #plot start time, model start time
# StartT=datenum(2017,8,20); StartT_model=datenum(2017,8,4) #plot start time, model start time

# runs=['elev.Coast_6b.RUN1k_ZG','elev.Coast_6b.RUN03c_ZG']; tags=['RUN1k','RUN03c']; bpfile='Coast_6b.bp'
runs=['elev.Coast_6a.RUN1j_ZG','elev.Coast_6a.RUN03c_ZG']; tags=['RUN1j','RUN03c']; bpfile='Coast_6a.bp'
# runs=['elev.Coast_6b.RUN04a_ZG']; tags=['RUN04a']; bpfile='Coast_6b.bp'
# runs=['elev.TB_5_Original.RUN04a_ZG','elev.TB_5.RUN04a_ZG']; tags=['RUN04a_Original','RUN04a']; bpfile='TB_5.bp'
# runs=['elev.TB_5_Original.RUN03c_ZG','elev.TB_5.RUN03c_ZG','elev.TB_5.RUN1k_ZG']; tags=['RUN03c_Original','RUN03c','RUN1k']; bpfile='TB_5.bp'

#do Harmonic Analysis
iHA=1

#figure save names.
sname=None
# sname='NWM_elev_Coast_6b_RUN1k_RUN03c_selected' # name for savefig
# sname='NWM_elev_Coast_6b_RUN1k_RUN03c' # name for savefig
# sname='NWM_elev_Coast_6a_RUN1j_RUN03c_BerginPoint' # name for savefig
# sname='NWM_elev_Coast_6b_RUN04a' # name for savefig
# sname='Test_TB_4_RUN1k' # name for savefig
sname='NWM_figure_tmp'

#linestypes for different runs
colors='kgbc'; lstyle=['-','-','-']; Markers=['*','.','None']
# colors='gbc'; lstyle=['-','-','-']; Markers=['None','.','None']

#msl_to_navd shift
# msl_shift_station=array([8760721]); msl_shift=[-0.5]
msl_shift_station=None; msl_shift=None

#shift all data
# all_shift=[0,-0.2059,  0,0,0,0]
all_shift=[0,0, 0,0,0,0]

# ym=[0,1]
ym=None

#stations to plots
stations=None
# stations=[110,*range(112,115),116,*range(28,32),*range(33,40)]


#------------------------------------------------------------------------------
#read station data
#------------------------------------------------------------------------------
#read station.bp
fid=open('./stations.txt'); staname=dict(zip(*array([i.strip().split(',') for i in fid.readlines()]).T)); fid.close()
bp=read_schism_bpfile(bpfile);  bp.staname=array([staname[i] for i in bp.station])

#for subset of stations
if stations is not None:
    staid=array(stations).astype('int')-1
    bp.nsta=len(stations); bp.x=bp.x[staid]; bp.y=bp.y[staid]; bp.z=bp.z[staid]; bp.station=bp.station[staid]; bp.staname=bp.staname[staid]
else:
    staid=arange(bp.nsta).astype('int')
    stations=arange(bp.nsta)+1

#read model results
Model=[]
for m, run in enumerate(runs):
    Si=npz_data(); Data=loadtxt(run); Si.time=Data[:,0]+StartT_model; Si.elev=Data[:,1:][:,staid]
    Model.append(Si)

#read obs data
#C=loadz('./noaa_elev_navd.npz')
#C1=loadz('./noaa_elev_msl_2010_2020.npz')
C=loadz('/sciclone/data10/wangzg/Database/elev/noaa_elev_navd.npz')
C1=loadz('/sciclone/data10/wangzg/Database/elev/noaa_elev_msl_2010_2020.npz')

#-----------------------plot---------------------------------------------------
MAE=[]; Corr=[]; ME=[]; Amp_model=[]; datums=[]; Pha_model=[];  Amp_obs=[]; Pha_obs=[]
for m in arange(10):
    i1=m*16; i2=(m+1)*16; iflag=0
    if bp.nsta<=i1: continue
    figure(figsize=[19,9.4])
    for i in arange(i1,min(i2,bp.nsta)):
        iflag=iflag+1;
        station=int(bp.station[i])

        subplot(8,2,iflag)

        #plot obs
        lobs='None'
        fp=(C.station==station)*(C.time>(StartT_model-10))*(C.time<(StartT_model+130)); oti=C.time[fp]; oyi=C.elev[fp]; lobs='navd'
        if sum(fp)==0:
            fp=(C1.station==station)*(C1.time>(StartT_model-10))*(C1.time<(StartT_model+130)); oti=C1.time[fp]; oyi=C1.elev[fp]; lobs='msl'

        #add nan data between oti
        if len(oti)>100:
            ts=find_continuous_sections(oti,1.0); eoti=array([i[-1]+1/24 for i in ts.sections]); eoyi=ones(len(eoti))*nan
            oti=r_[oti,eoti]; oyi=r_[oyi,eoyi]; sind=argsort(oti); oti=oti[sind]; oyi=oyi[sind]

        plot(oti-StartT,oyi,'r')

        #plot model
        for n, run in enumerate(runs):
            mti=Model[n].time; myi=Model[n].elev[:,i]
            if all_shift[n]!=0: myi=myi+all_shift[n]

            #add msl_to_navd shift
            if msl_shift_station is not None:
                if station in msl_shift_station:
                    sid=nonzero(station==msl_shift_station)[0][0]
                    myi=myi+msl_shift[sid]
                    if n==0: lobs=lobs+', shift={}m'.format(msl_shift[sid])

            plot(mti-StartT,myi,color=colors[n],linestyle=lstyle[n],Marker=Markers[n],ms=3,alpha=0.85)

        if ym is not None: setp(gca(),ylim=ym)

        #note
        xts=[*arange(0,120,1)]; xls=[num2date(StartT).strftime('%m/%d/%Y'),*xts[1:]]
        setp(gca(),xticks=xts,xticklabels=xls,xlim=[0,10])
        if (iflag<=14) and (i<(bp.nsta-2)): setp(gca(),xticklabels=[])
        gca().xaxis.grid('on')
        text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'{}, {}: {}, {}'.format(stations[i],bp.station[i],bp.staname[i],lobs))
        if iflag==1: legend(['Obs',*tags])

        #do statistics
        MAEi=[]; Corri=[]; MEi=[]
        Amp_modeli=[]; Pha_modeli=[];  Amp_obsi=[]; Pha_obsi=[]

        fpn=~isnan(oyi); oti=oti[fpn]; oyi=oyi[fpn]
        fpt=(oti>=Model[-1].time.min())*(oti<=Model[-1].time.max()); otii=oti[fpt]
        if len(otii)!=0:
            for n, run in enumerate(runs):
                #get data pairs
                mti=Model[n].time; myi=Model[n].elev[:,i];
                if all_shift[n]!=0: myi=myi+all_shift[n]
                fpt=(oti>=mti.min())*(oti<=mti.max()); otii=oti[fpt]; oyii=oyi[fpt]
                myii=interpolate.interp1d(mti,myi)(otii)

                #statistics
                st=get_stat(myii,oyii); MEii=mean(myii)-mean(oyii)
                if n==(len(runs)-1): text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.6*diff(ylim()),'MAE={:6.3f}, Corr={:6.3f},ME={:6.3f}'.format(st.MAE,st.R,MEii))
                MAEi.append(st.MAE); Corri.append(st.R); MEi.append(MEii)

                #HA
                if iHA==1:
                    ts=find_continuous_sections(otii,1)
                    if diff(ts.section_max)>15:
                        fp=(otii>=ts.section_max[0])*(otii<=ts.section_max[1])
                        HAm=harmonic_analysis(myii[fp],1/24,StartT=datenum(2018,6,17))
                        HAo=harmonic_analysis(oyii[fp],1/24,StartT=datenum(2018,6,17))
                    else:
                        HAm=npz_data(); HAm.amplitude=ones(9)*nan; HAm.phase=ones(9)*nan
                        HAo=npz_data(); HAo.amplitude=ones(9)*nan; HAo.phase=ones(9)*nan

                    #save
                    Amp_modeli.append(HAm.amplitude); Pha_modeli.append(HAm.phase)
                    Amp_obsi.append(HAo.amplitude); Pha_obsi.append(HAo.phase)

            MAE.append(MAEi); Corr.append(Corri); ME.append(MEi)
        else:
            Amp_modeli=ones([len(runs),9])*nan; Pha_modeli=ones([len(runs),9])*nan
            Amp_obsi=ones([len(runs),9])*nan; Pha_obsi=ones([len(runs),9])*nan
            MAE.append([*(ones(len(runs))*nan)]); Corr.append([*(ones(len(runs))*nan)]); ME.append([*(ones(len(runs))*nan)])

        Amp_model.append(Amp_modeli); Pha_model.append(Pha_modeli)
        Amp_obs.append(Amp_obsi); Pha_obs.append(Pha_obsi); datums.append(lobs)

    gcf().tight_layout()
    # savefig('NWM_elev_Coast_4_RUN1c_RUN1e_{}'.format(m),dpi=300)
    if sname!=None: savefig('{}_{}'.format(sname,m),dpi=300)

MAE=array(MAE); Corr=array(Corr); ME=array(ME)
if iHA==1: Amp_model=array(Amp_model); Pha_model=array(Pha_model); Amp_obs=array(Amp_obs); Pha_obs=array(Pha_obs); tidal_name=HAm.tidal_name

print('MAE: ', MAE.mean(axis=0))
print('Corr: ', Corr.mean(axis=0))
print('ME: ', ME.mean(axis=0))

#----plot HA------------
if iHA==0: sys.exit()
# figure(figsize=[19,9.4])
figure(figsize=[40,9.4])
nsta=bp.nsta; xi=arange(nsta)+1

markers=['*','^','s','+']
subplot(4,1,1)
plot(xi,Amp_obs[:,0,5],'r.',ms=8,alpha=0.75)
for m, run in enumerate(runs):
    plot(xi,Amp_model[:,m,5],'.',color=colors[m],marker=markers[m],ms=6,alpha=0.5)
setp(gca(),xlim=[0.5,nsta+0.5],xticks=[*(arange(nsta)+1)])
gca().xaxis.grid('on')
legend(['obs',*tags])
text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'Amp.: M2',fontsize=14,fontweight='bold')

subplot(4,1,2)
plot(xi,Pha_obs[:,0,5],'r.',ms=8,alpha=0.75)
for m, run in enumerate(runs):
    plot(xi,Pha_model[:,m,5],'.',color=colors[m],marker=markers[m],ms=6,alpha=0.5)
setp(gca(),yticks=[-pi,-pi/2,0,pi/2,pi],yticklabels=[*arange(-180,200,90)],ylim=[-pi,pi])
setp(gca(),xlim=[0.5,nsta+0.5],xticks=[*(arange(nsta)+1)])
gca().xaxis.grid('on')
text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'Phase: M2',fontsize=14,fontweight='bold')


subplot(4,1,3)
plot(xi,Amp_obs[:,0,2],'r.',ms=8,alpha=0.75)
for m, run in enumerate(runs):
    plot(xi,Amp_model[:,m,2],'.',color=colors[m],marker=markers[m],ms=6,alpha=0.5)
setp(gca(),xlim=[0.5,nsta+0.5],xticks=[*(arange(nsta)+1)])
gca().xaxis.grid('on')
legend(['obs',*tags])
text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'Amp.: K1',fontsize=14,fontweight='bold')

subplot(4,1,4)
plot(xi,Pha_obs[:,0,2],'r.',ms=8,alpha=0.75)
for m, run in enumerate(runs):
    plot(xi,Pha_model[:,m,2],'.',color=colors[m],marker=markers[m],ms=6,alpha=0.5)
setp(gca(),yticks=[-pi,-pi/2,0,pi/2,pi],yticklabels=[*arange(-180,200,90)],ylim=[-pi,pi])
setp(gca(),xlim=[0.5,nsta+0.5],xticks=[*(arange(nsta)+1)],xticklabels=bp.station)
# xticks([*bp.station],rotation='vertical')
xticks(rotation=90)
gca().xaxis.grid('on')
text(xlim()[0]+0.1*diff(xlim()),ylim()[0]+0.8*diff(ylim()),'Phase: K1',fontsize=14,fontweight='bold')


gcf().tight_layout()
if sname is not None: savefig('HA_{}'.format(sname))
