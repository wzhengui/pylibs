#!/usr/bin/env python3
#generate bctides.in's tide harmonics
#specify inputs in the bottom 
from pylib import *

def get_tide_nodal(tide_name,StartT,nday):
    tdir='tide_fac_improved'
    #compile the code
    os.system('cd {}; ifort -o tide_fac_improved tf_main.f90 tf_selfe.f90'.format(tdir))

    #write input
    with open('{}/tide.in'.format(tdir),'w+') as fid:
        fid.write('{}\n'.format(nday))
        fid.write('{} {} {} {}\n'.format(*flipud(StartT)))
        fid.write('0\n');

    #run the code
    os.system('cd {}; ./tide_fac_improved <tide.in'.format(tdir))

    #read outputs
    with open('{}/tide_fac.out'.format(tdir),'r') as fid:
        lines=fid.readlines()

    tnodal=[]; tear=[]
    for ti in tide_name:
        for line in lines:
            tline=line.strip()
            R=re.match('{}'.format(ti),tline)
            if R!=None:
                num=re.findall(r"[-+]?\d*\.\d+|\d+", tline[2:])
                tnodal.append(num[0])
                tear.append(num[1])

    tnodal=array(tnodal); tear=array(tear)
    return tnodal,tear

def get_tide_amp_freq(tide_name):
    T=loadz('./tide_fac_const/tide_fac_const.npz');

    amp=[]; freq=[];
    for ti in tide_name:
        ind=nonzero(T.name==ti)[0];
        amp.append(T.amp[ind])
        freq.append(T.freq[ind])

    amp=squeeze(array(amp)); freq=squeeze(array(freq))
    return amp, freq

def copy_inputs(bdir):

    #link database
    files=['eastward_velocity','northward_velocity','fes2014a_loadtide','fes2014b_elevations','fes2014b_elevations_extrapolated']
    for fi in files:
       os.system("ln -sf {}/{}".format(bdir,fi))

    #copy files
    files=['gen_harm_FES.m','README']
    for fi in files:
       os.system("cp {}/{} ./".format(bdir,fi))

    #copy tide dir
    files=['tide_fac_improved','tide_fac_const']
    for fi in files:
       os.system("cp -r {}/{} ./".format(bdir,fi))

if __name__=="__main__":
    tide_name=['O1','K1','Q1','P1','M2','S2','K2','N2'];#check order in gen_harm_FES.m
    StartT=[2005,1,1,0]; #year,month,day,hour
    nday=367; #number of days
    ibnds=[1];  #order of boundary segments

    #---setup inputs-------------------
    bdir='/home/zwang/FES2014';
    copy_inputs(bdir);

    #----get tide amplitude and frequcy--
    tamp,tfreq=get_tide_amp_freq(tide_name)

    #---get nodal factor
    tnodal,tear=get_tide_nodal(tide_name,StartT,nday)

    #write 1st open boundary information
    gd0=read_schism_hgrid('../hgrid.gr3')
    gd=read_schism_hgrid_ll('../hgrid.ll',gd0)

    #write bctides file
    with open('bctides.in','w+') as fid:
        fid.write('{:02}/{:02}/{:4} {:02}:00:00 GMT\n'.format(StartT[1],StartT[2],StartT[0],StartT[3]))
        fid.write('{:d} 40.0 !ntip\n'.format(len(tide_name)))
        for i in arange(len(tide_name)):
            fid.write('{}\n'.format(tide_name[i]))
            fid.write('{} {:<.6f}  {:<.9e}  {}  {}\n'.format(tide_name[i][1],tamp[i],tfreq[i],tnodal[i],tear[i]))

        fid.write('{:d} !nbfr\n'.format(len(tide_name)))
        for i in arange(len(tide_name)):
            fid.write('{}\n'.format(tide_name[i]))
            fid.write('  {:<.9e}  {}  {}\n'.format(tfreq[i],tnodal[i],tear[i]))

        fid.write('{} !nope\n'.format(gd.nob))
        
        #write each open boundary information
        for ibnd in ibnds:        
            fid.write('{} 5 5 4 4 !ocean\n'.format(gd.nobn[ibnd-1]))
        
            #generate boundary information
            bnodes=gd.iobn[ibnd-1].astype('int')
            lxi=gd.x[bnodes]; lyi=gd.y[bnodes]; lzi=gd.dp[bnodes]
            with open('open.ll','w+') as fid2:
                for i in arange(len(bnodes)):
                    fid2.write("{:10d}  {:.7e}  {:.7e}  {:.7e}\n".format(bnodes[i]+1,lxi[i],lyi[i],lzi[i]))
    
            #generate ap.out
            os.system('matlab -nodisplay <gen_harm_FES.m')
            with open('ap.out','r') as fid3:
                lines_ap=fid3.readlines();
            
            for line in lines_ap:
                fid.write(line)
            fid.write('1.0 !TEM nudge\n1.0 !SAL nudge\n')
        
        #write river boundary information
        #ibnds_river=[2,3,4,5]
        #rivers=['Coyote','San Joaquin','Sacramento','Napa']
        #for m in arange(len(rivers)):
        #    ibnd=ibnds_river[m]
        #    fid.write('{} 0 1 1 1 !{}\n'.format(gd.nobn[ibnd-1],rivers[m]))
        #    fid.write('1.0 !TEM nudge\n1.0 !SAL nudge\n')
        

          

