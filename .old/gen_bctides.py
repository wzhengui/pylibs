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

def get_tide_name(fname='gen_harm_FES.m',keyword='const'):
    #read in fname
    fid=open(fname,'r'); lines=fid.readlines(); fid.close()

    #find the line starting with keyword
    for line in lines:
        if line.strip().startswith(keyword):
            sline=line.strip()
            break

    #parse tidal names
    i1=sline.find('{')+1; i2=sline.find('}')
    tide_name=[i.replace("'",'').upper() for i in sline[i1:i2].split(',')]

    return tide_name

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
    #input
    StartT=[2018,6,17,0]; #year,month,day,hour
    nday=121; #number of days
    ibnds=[1];  #order of open boundary segments (starting from 1)
    iZ0=1; zconst=0.05 #iZ0==1 add Z0 constant (value=zconst)

    #---setup inputs-------------------
    bdir='/sciclone/data10/wangzg/FES2014';
    copy_inputs(bdir);

    #----get tide names from gen_harm_FES.m
    tide_name=get_tide_name(fname='gen_harm_FES.m',keyword='const')

    #----get tide amplitude and frequcy--
    tamp,tfreq=get_tide_amp_freq(tide_name)

    #---get nodal factor
    tnodal,tear=get_tide_nodal(tide_name,StartT,nday)

    #write 1st open boundary information
    gd0=read_schism_hgrid('../hgrid.gr3')
    gd=read_schism_hgrid_ll('../hgrid.ll',gd0)

    #write bctides file
    with open('bctides.in','w+') as fid:
        fid.write('!{:02}/{:02}/{:4} {:02}:00:00 UTC\n'.format(StartT[1],StartT[2],StartT[0],StartT[3]))
        fid.write(' {:d}  50.000 !number of earth tidal potential, cut-off depth for applying tidal potential\n'.format(len(tide_name)))
        for i in arange(len(tide_name)):
            fid.write('{}\n'.format(tide_name[i]))
            fid.write('{} {:<.6f}  {:<.9e}  {}  {}\n'.format(tide_name[i][1],tamp[i],tfreq[i],tnodal[i],tear[i]))

        fid.write('{:d} !nbfr\n'.format(int(len(tide_name)+iZ0)))
        if iZ0==1: fid.write('Z0\n0. 1. 0.\n') 
        for i in arange(len(tide_name)):
            fid.write('{}\n'.format(tide_name[i]))
            fid.write('  {:<.9e}  {}  {}\n'.format(tfreq[i],tnodal[i],tear[i]))

        fid.write('{} !nope\n'.format(gd.nob))

        #write each open boundary information
        for ibnd in ibnds:
            fid.write('{} 3 3 0 0 !ocean\n'.format(gd.nobn[ibnd-1]))

            #generate boundary information
            bnodes=gd.iobn[ibnd-1].astype('int');
            lxi=gd.x[bnodes]; lyi=gd.y[bnodes]; lzi=gd.dp[bnodes];
            with open('open.ll','w+') as fid2:
                for i in arange(len(bnodes)):
                    fid2.write("{:10d}  {:.7e}  {:.7e}  {:.7e}\n".format(bnodes[i]+1,lxi[i],lyi[i],lzi[i]))

            #generate ap.out
            os.system('matlab -nodisplay <gen_harm_FES.m')
            with open('ap.out','r') as fid3:
                lines_ap=fid3.readlines();

            #add Z0
            if iZ0==1: fid.write('Z0\n'+'{} 0\n'.format(zconst)*gd.nobn[ibnd-1])   

            for i,line in enumerate(lines_ap):
                fid.write(line)

                #add Z0
                if (i+1)==(1+gd.nobn[ibnd-1])*len(tide_name) and iZ0==1:
                    fid.write('Z0\n'+'0 0 0 0\n'*gd.nobn[ibnd-1])   
                
            #fid.write('1.0 !TEM nudge\n1.0 !SAL nudge\n')

        #write river boundary information
        #ibnds_river=[2,3,4,5]
        #rivers=['Coyote','San Joaquin','Sacramento','Napa']
        #for m in arange(len(rivers)):
        #    ibnd=ibnds_river[m]
        #    fid.write('{} 0 1 1 1 !{}\n'.format(gd.nobn[ibnd-1],rivers[m]))
        #    fid.write('1.0 !TEM nudge\n1.0 !SAL nudge\n')




