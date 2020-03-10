#!/usr/bin/evn python3
from pylib import *

#-------misc-------------------------------------------------------------------
def close_data_loop(xi):
    #close data loop by add xi[0] in the end
    if len(xi)<3:
        vi=xi;
    else:
        if xi.ndim==1:
            if xi[0]!=xi[-1]:
                vi=r_[xi,xi[0]];
            else:
                vi=xi
        elif xi.ndim==2:
            if xi[0][0]!=xi[-1][0]:
                vi=r_[xi,xi[0][None,:]]
            else:
                vi=xi
    return vi

#-----find all continous sections of a time series
def find_continuous_sections(xi,dx):
    #analyze time series of xi; when a gap is larger than dx, it is a break;
    dxi=diff(xi); ind=nonzero(dxi>dx)[0];

    #divide xi into sub-sections
    sections=[];
    if len(ind)==0:
        sii=[xi[0],xi[-1]]; sections.append(sii); section_max=sii;
    else:
        i1=0; i2=ind[0]; sii=[xi[i1],xi[i2]]; sections.append(sii); section_max=sii;
        for i in arange(len(ind)-1):
            i1=ind[i]+1; i2=ind[i+1]
            sii=[xi[i1],xi[i2]]; sections.append(sii);
            if diff(section_max)<diff(sii): section_max=sii
        i1=ind[-1]+1; sii=[xi[i1],xi[-1]]; sections.append(sii)
        if diff(section_max)<diff(sii): section_max=sii

    #save results
    S=npz_data();
    S.sections=array(sections);
    S.section_max=array(section_max)

    return S

#-------str2num----------------------------------------------------------------
def str2num(line,*args):
    num=str2num_process(line,*args)
    if isinstance(num[0],float):
        num=num.astype('float64')
    else:
        num=[s.astype('float64') for s in num]
        num=array(num)
    return num

@np.vectorize
def str2num_process(line,*args):
    if len(args)>0:
        if len(args)>1:
            for i in range(len(args)-1):
                line=line.replace(arg)
        line=line.replace(args[0],',')
    else:
        line=line.replace(';',',').replace(' ',',')
    linei=[s for s in line.split(',') if s]
    fc=np.vectorize(lambda x: np.float64(x))
    return fc(linei).astype('object')

@np.vectorize
def remove_tail(line):
    li=line.rstrip();
    ind=li.find('!');
    if ind!=-1:
        li=li[:ind]
    ind=li.find('=');
    if ind!=-1:
        li=li[:ind]
    return li

#-------date_proc--------------------------------------------------------------
def datenum_0(*args):
    if len(args)==1:
        args=args[0];

    args=array(args)
    args=args.astype('int')
    return datetime.datetime(*args)

def datenum(*args,doy=0):
    args=array(args)
    e1=args[0]

    if hasattr(e1, "__len__"):
        if not hasattr(e1[0],"__len__"):
            f=datenum_0(*e1)
        else:
            f=apply_along_axis(datenum_0,1,e1)
    else:
        f=datenum_0(*args)
    if doy==0:
        return date2num(f)
    else:
        return f

def get_xtick(xi=None,sformat='%Y',method=0):
    #return time ticks for plot
    #method: 0=outputs ticks of xi; xi should be datenum
    #        1=yearly interal; 2=monthly interval; 3=daily interaly (xi is list of year if not None)
    #sformat:
    #   %d=01;     %-d=1                           (day of month)
    #   %b=Jan;    %B=January;   %m=01;     %-m=1  (month)
    #   %y=01;     %-y=1;        %Y=2000           (year)
    #   %H=09;     %-H=9;        %I=[00,12] %-I=1  (hour)
    #   %M=09;     %-M=9;                          (minute)
    #   %S=09;     %-S=9                           (second)
    #
    #   %a=MON;    %A=Monday;    %w=[0,6]          (week)
    #   %U=[00,53];    %W=[00,53];                 (week of year)
    #   %j=045;    %-j=45;                         (day of year)
    #   %p=[AM,PM]                                 (moring, or afternoon)
    #
    #   %c='Fri Jan 25 04:05:02 2008'
    #   %x='01/25/08';
    #   %X='04:05:02'

    #-----determine ti --------------------------------
    if method==0:
        if xi==None:
            ti=array([datenum(i,1,1) for i in arange(1950,2050)])
        else:
            ti=array(xi)
    elif method==1:
        if xi==None: xi=[1950,2050]
        ti=[]
        [ti.append(datenum(yeari,1,1)) for yeari in arange(xi[0],xi[-1])];
        ti=array(ti)
    elif method==2:
        if xi==None: xi=[1950,2050]
        ti=[];
        [[ti.append(datenum(yeari,j+1,1)) for j in arange(12)] for yeari in arange(xi[0],xi[-1])];
        ti=array(ti)
    elif method==3:
        if xi==None: xi=[1950,2050]
        ti=[];
        [ti.append(i) for i in arange(datenum(xi[0],1,1),datenum(xi[-1]+1,1,1))]
        ti=array(ti)

    #---------------------------------------------------
    xtick=ti;
    xticklabel=array([num2date(tii).strftime(sformat) for tii in ti]);

    return xtick,xticklabel

#-------loadz------------------------------------------------------------------
class npz_data(object):
    def __init__(self):
        pass

def save_npz(fname,C):
    #npz_vars=[ npz_vari.split(':')[0] for npz_vari in C.VINFO ];
    npz_vars=list(C.__dict__.keys())
    if 'VINFO' in npz_vars: npz_vars.remove('VINFO')

    save_str='savez_compressed("{}" '.format(fname);
    for vari in npz_vars:
        save_str=save_str+',{}=C.{}'.format(vari,vari)
    save_str=save_str+')';
    #print(save_str)
    exec(save_str)

def loadz(fname,med=1):
    #med=1: return class format; med=2:return dict format
    data0=load(fname,allow_pickle=True)
    keys0=data0.keys()

    if med==1:
        zdata=npz_data();
    else:
        zdata2={}

    VINFO=[]
    for keyi in keys0:
        datai=data0[keyi];
        if datai.dtype==dtype('O'): datai=datai[()]
        if med==1:
            exec('zdata.'+keyi+'=datai')
        else:
            zdata2[keyi]=datai

        #gather information about datai
        vinfo=keyi+": "+type(datai).__name__
        if isinstance(datai,list):
            vinfo=vinfo+'('+str(len(datai))+'), '
        elif isinstance(datai,np.ndarray):
            vinfo=vinfo+str(datai.shape)+', dtype='+str(datai.dtype)
        VINFO.append(vinfo)
    VINFO=array(VINFO)
    zdata.VINFO=VINFO

    return zdata if med==1 else zdata2

#-------mfft-------------------------------------------------------------------
def mfft(xi,dt):
    #input
    #xi: time series
    #dt: interval
    #
    #output
    #perid[period],afx[amplitude],pfx[phase]
    N=xi.size;
    fx=fft(xi);
    afx=abs(fx[1:N//2])*2.0/N;
    pfx=angle(fx[1:N//2]);
    period=dt*N/arange(1,N//2);
    return period,afx,pfx

#def compute_cofficient(myi,oyi):
#    #compute different evaluation coefficients
#    N=len(myi)
#    mmyi=myi.mean(); moyi=oyi.mean();
#    emyi=myi-mmyi; eoyi=oyi-moyi; e=myi-oyi
#    stdm=std(emyi); stdo=std(eoyi)
#
#    SS_tot=sum((oyi-moyi)**2)
#    SS_reg=sum((myi-moyi)**2);
#    SS_res=sum(e**2);
#
#    #evaluation coefficient
#    ms=1-SS_res/sum((abs(myi-moyi)+abs(oyi-moyi))**2) #model skill
#    r=mean((myi-mmyi)*(oyi-moyi))/stdm/stdo #correlation coefficient
#    r2=1-SS_res/SS_tot  #R2
#    rmse=sqrt(sum(e**2)/N) #RMSE
#    mae=mean(abs(e))         #MAE
#    me=mean(e)               #ME
#
#    return [ms,r,r2,rmse,mae,me]

def command_outputs(code,shell=True):
    import subprocess
    p=subprocess.Popen(code,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=shell)
    stdout,stderr=p.communicate()
    if stdout!=None: stdout=stdout.decode('utf-8')
    if stderr!=None: stderr=stderr.decode('utf-8')

    out=npz_data()
    out.stdout=stdout
    out.stderr=stderr
    return out

#def near_pts(pts,pts0):
#    #pts[n,2]: xy of points
#    #pts0[n,2]: xy of points
#    #return index of pts0(x0,y0) that pts(x,y) is nearest
#    x=pts[:,0]; y=pts[:,1]
#    x0=pts0[:,0]; y0=pts0[:,1]
#    dist=(x[None,:]-x0[:,None])**2+(y[None,:]-y0[:,None])**2
#
#    ind=[];
#    for i in arange(x.shape[0]):
#        disti=dist[:,i];
#        ind.append(nonzero(disti==min(disti))[0][0])
#    ind=array(ind);
#
#    return ind

def near_pts(pts,pts0,method=0,N=100):
    #pts[n,2]: xy of points
    #pts0[n,2]: xy of points
    #return index of pts0(x0,y0) that pts(x,y) is nearest
    #method=0: quick method by subgroups (N); method=1: slower methods
    if method==0:
        p=pts[:,0]+1j*pts[:,1]
        p0=pts0[:,0]+1j*pts0[:,1]

        N=min(N,len(p)-1);
        #--divide pts into subgroups based on dist, each group size is about N---------
        ps0=[]; ps=[]; ds=[]; inds=[]
        i0=0; mval=1e50; pflag=p==None
        while(True):
            ps0.append(p[i0]);
            dist=abs(p-p[i0]); dist[pflag]=mval;
            ind=argsort(dist); sdist=dist[ind]; mds=sdist[N]; i0=ind[N]
            fp=dist<mds; psi=p[fp]; dsi=max(dist[fp]); pflag[fp]=True;
            ps.append(psi); ds.append(dsi); inds.append(nonzero(fp)[0])
            if mds==mval: break

        #---find nearest pts for each subgroup-----------------------------------------
        inds0=[]
        for m in arange(len(ps)):
            dist=abs(p0-ps0[m]); dsi=ds[m]; psi=ps[m]
            #---find radius around ps0[m]
            dsm=dsi
            while(True):
                fp=dist<=dsm;
                if sum(fp)>0:
                    break;
                else:
                    dsm=dsm+dsi;
            #-----subgroup pts of p0---------------------------------
            fp=dist<=(dsm+2*dsi); ind0=nonzero(fp)[0]; p0i=p0[ind0];
            psii=psi[:,None]; p0ii=p0i[None,:]
            dist=abs(psii-p0ii); indi=dist.argmin(axis=1);
            inds0.append(ind0[indi]);

        #-arrange index-----------------
        pind=array([]).astype('int'); pind0=array([]).astype('int')
        for m in arange(len(inds)):
            pind=r_[pind,inds[m]]
            pind0=r_[pind0,inds0[m]]

        ind=argsort(pind);
        sind=pind0[ind]
    else:
        n=pts.shape[0]; n0=pts0.shape[0]
        N=max(min(1e7//n0,n),1e2)
        print('total pts: {}'.format(n));

        i0=int(0); i1=int(N);
        while(True):
            print('processing pts: {}-{}'.format(i0,i1));
            x=pts[i0:i1,0]; y=pts[i0:i1,1]
            x0=pts0[:,0]; y0=pts0[:,1]
            dist=(x[None,:]-x0[:,None])**2+(y[None,:]-y0[:,None])**2

            indi=[];
            for i in arange(x.shape[0]):
                disti=dist[:,i];
                indi.append(nonzero(disti==min(disti))[0][0])

            if i0==0:
                ind=array(indi);
            else:
                ind=r_[ind,squeeze(array(indi))];
            #next step
            i0=int(i0+N); i1=int(min(i1+N,n))
            if i0>=n: break
        sind=ind

    return sind

def inside_polygon(pts,px,py):
    #pts[n,2]: xy of points
    #px[nploy,x]: x coordiations of polygons
    #py[nploy,y]: x coordiations of polygons
    #return index of polygons that pts resides in
    npy=px.shape[0]; nv=px.shape[1];
    ind=[];
    for i in arange(pts.shape[0]):
        pxi=pts[i,0]; pyi=pts[i,1]; fi=ones(npy)
        for m in arange(nv):
            xi=c_[ones(npy)*pxi,px[:,m],px[:,mod(m+1,nv)]]
            yi=c_[ones(npy)*pyi,py[:,m],py[:,mod(m+1,nv)]]
            area=signa(xi,yi)
            fp=area<=0; fi[fp]=0;
        indi=nonzero(fi)
        ind.append(indi)

    ind=squeeze(array(ind))

    return ind

def signa(x,y):
    #compute signed area
    if x.ndim==1:
        area=((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]))/2
    elif x.ndim==2:
        area=((x[:,0]-x[:,2])*(y[:,1]-y[:,2])-(x[:,1]-x[:,2])*(y[:,0]-y[:,2]))/2
    area=squeeze(area)
    return area

def mdivide(A,B):
    #perform matrix division B/A
    if A.ndim==1: A=A[None,:]
    if B.ndim==1: B=B[None,:]
    A2=inv(A@A.T);
    B2=B@A.T
    return squeeze(B2@A2)

def lpfilt(data,delta_t,cutoff_f):
    #import gc
    #low pass filter for 1D (data[time]) or 2D (data[time,array]) array;
    ds=data.shape
    mn=data.mean(axis=0)
    data=data-mn
    P=fft(data,axis=0)
    #data=None; gc.collect()
    N=ds[0];
    filt=ones(N)
    k=int(floor(cutoff_f*N*delta_t))
    filt[k]=0.715
    filt[k+1]=0.24
    filt[k+2]=0.024
    filt[k+3:N-(k+4)]=0.0
    filt[N-(k+4)]=0.024
    filt[N-(k+3)]=0.24
    filt[N-(k+2)]=0.715

    #expand filter dimensions
    fstr=',None'*(len(ds)-1); s=npz_data()
    exec('s.filt=filt[:{}]'.format(fstr))

    #remove high freqs
    P=P*s.filt

    #lp results
    fdata=real(ifft(P,axis=0))+mn

    return fdata

def reload(gfunc):
    #reload modules,mainly used for coding debug 
    #usage: reload(globals())
    import inspect,imp
    mods=['mylib','schism_file','mpas_file']
    for modi in mods:
        imp.reload(sys.modules[modi])
        #get all module functions
        fs=[];afs=inspect.getmembers(sys.modules[modi],inspect.isfunction);
        for fsi in afs: 
            if inspect.getmodule(fsi[1]).__name__!=modi: continue
            if fsi[0] not in gfunc.keys(): continue 
            fs.append(fsi)  
        #refresh module functions
        for fsi in fs:
            if gfunc[fsi[0]]!=fsi[1]: gfunc[fsi[0]]=fsi[1] 
    return

def wipe(n=50):
    print('\n'*n)
    return

#def clear_globals():
#    allvar = [var for var in globals() if var[0] != "_"]
#    for var in allvar:
#       #global var
#       #del var;
#       #del globals()[var]
#       #exec('del global()['+var+']')
#       exec('global '+var)
#       exec('del '+var)

def smooth(xi,N):
    #xi: time series
    #N: window size (if N is even, then N=N+1)
    if mod(N,2)==0: N=N+1
    nz=int((N-1)/2)

    X=xi; SN=ones(len(xi))*N
    for i in arange(nz):
        nmove=i+1; ii=-i-1;
        xext=zeros(nmove)
        X=X+r_[xext,xi[:ii]]
        for j in range(nmove):
            SN[j]=SN[j]-1

    for i in arange(nz):
        nmove=i+1; ii=i+1
        xext=zeros(nmove)
        X=X+r_[xi[ii:],xext]
        for j in arange(nmove):
            jj=-j-1
            SN[jj]=SN[jj]-1
    SX=np.divide(X,SN)
    return SX;

def DaytimeLength(Lat,Doy):
    #calculate daytime length
    #D=DaytimeLength(Lat,Doy)
    #Lat: latitude, Doy: (1-365), sunrise=12-D/2, sunset=12+D/2
    P=arcsin(0.39795*cos(0.2163108 + 2*arctan(0.9671396*tan(0.00860*(Doy-186)))));
    T=(sin(0.8333*pi/180)+sin(Lat*pi/180)*sin(P))/cos(Lat*pi/180)/cos(P);
    D=24-(24/pi)*arccos(T);
    return D

def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)

def proj(fname0,format0,proj0,fname1,format1,proj1):
    #tranfrom projection of files: proj(fname0,format0,proj0,fname1,format1,proj1)
    #fname: file name
    #format: 0: SCHISM gr3 file; 1: SCHISM bp file; 2: xyz file; 3: xyz file with line number
    #proj: projection name (e.g. 'epsg:26918', 'epsg:4326')

    from read_schism_file import read_schism_hgrid
    #read file
    if format0==0:
        gd=read_schism_hgrid(fname0)
        np=gd.np; x=gd.x; y=gd.y; z=gd.dp
    elif format0==1:
        gd=read_schism_bpfile(fname0)
        np=gd.nsta; x=gd.x; y=gd.y; z=gd.z
    elif format0==2:
        gd=loadtxt(fname0)
        np=gd.shape[0]; x=gd[:,0]; y=gd[:,1]; z=gd[:,2]
    elif format0==3:
        gd=loadtxt(fname0)
        np=gd.shape[0]; x=gd[:,1]; y=gd[:,2]; z=gd[:,3]
    else:
        sys.exit('unknown format of {}'.format(fname0))

    #transform coordinate
    x1,y1=transform(Proj(init=proj0),Proj(init=proj1),x,y);

    #write file
    if format1==0:
        if format0!=0: sys.exit('{} should have gr3 format'.format(fname0))
        gd.x=x1; gd.y=y1;
        gd.write_hgrid(fname1,Info='!coordinate transformed: {}=>{}'.format(proj0,proj1))
    elif format1==1:
        if format0==1:
            gd.x=x1; gd.y=y1
        else:
            gd=schism_bpfile()
            gd.nsta=np; gd.x=x1; gd.y=y1; gd.z=z;
        gd.write_bpfile(fname1)
    elif format1==2 or format1==3:
        with open(fname1,'w+') as fid:
            for i in arange(np):
                if format1==2: fid.write('{} {} {}\n'.format(x1[i],y1[i],z[i]))
                if format1==3: fid.write('{} {} {} {}\n'.format(i+1,x1[i],y1[i],z[i]))
    else:
        sys.exit('unknown format of {}'.format(fname1))

def get_prj_file(prjname='epsg:4326',method=0):
    #return project name or entire database (dict)
    #
    #--online method-----
    # function to generate .prj file information using spatialreference.org
    #def getWKT_PRJ (epsg_code):
    #     import urllib
    #     # access projection information
    #     wkt = urllib.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code))
    #     # remove spaces between charachters
    #     remove_spaces = wkt.read().replace(" ","")
    #     # place all the text on one line
    #     output = remove_spaces.replace("\n", "")
    #     return output

    bdir=os.path.dirname(__file__)

    if method==0:
        S=loadz('{}/prj.npz'.format(bdir))
        return S.prj[prjname]
    elif method==1:
        S=loadz('{}/prj.npz'.format(bdir))
        return S.prj
    elif method==-1:
        #dir of *.prj files
        prj_dir=r'D:\Work\Database\projection\prj_files'

        #---processing all prj files-------------------------------------------
        fnames=os.listdir(prj_dir)

        prj=dict()
        for fname in fnames:
            R=re.match('epsg.(\d+).prj',fname);
            if not R: continue
            prj_num=int(R.groups()[0])
            with open('{}/{}'.format(prj_dir,fname),'r') as fid:
                line=fid.readline().strip();
                prj['epsg:{}'.format(prj_num)]=line
        #save prj file as a database
        S=npz_data();
        S.prj=prj
        save_npz('{}/prj'.format(bdir),S)

        return 0

#----------------convert_matfile_format----------------------------------------
def convert_matfile_format(file):
    #input for a directory or a matfile
    # r'C:\Users\Zhengui\Desktop\Observation2\DWR\SFBay_DWRData_SSI.mat'];
    #file=r'D:\OneDrive\Python\tem.mat';
    fname=[];
    if os.path.isdir(file):
        sfile=os.listdir(file)
        for sfilei in sfile:
            ename=sfilei.rstrip().split('.')[-1]
            if ename=='mat':
                fname.append(file+os.sep+sfilei)
    else:
        fname=[file];

    fid=open('log_matfile_convert.txt','w+')
    #convert mat format from v7.3 to v7
    cdir=os.getcwd();
    os.chdir('D:\OneDrive\Python')
    import matlab.engine
    eng = matlab.engine.start_matlab()

    for fn in fname:
        print('converting matfile: '+fn)
        dname=os.path.dirname(fn)
        bname=os.path.basename(fn).split('.')[0]
        fnv7=dname+os.sep+bname+'_v7'
        fnz=dname+os.sep+bname
        eflag=eng.convert_matfile_format(fn,fnv7,nargout=1)
        if eflag!=0:
            print('convert flag is %d: %s\n' % (eflag, fn));
            fid.write('convert flag is %d: %s\n' % (eflag,fn))
            continue
        convert_matfile(fnz,fnv7)
        os.remove(fnv7+'.mat')
    fid.close()
    os.chdir(cdir)

#convert mat to npz
def convert_matfile(fnz,fnv7):
    fc=np.vectorize(lambda x: x[0])
    C=sp.io.loadmat(fnv7+'.mat')
    vn=C.keys();

    iflag=0;Y={};
    for vni in vn:
        if vni[:2]=='__':
            continue
        Ci=C[vni];
        if issubdtype(Ci.dtype,np.number):
            Yi=Ci.copy();
        else:
            Yi=fc(Ci)
        if vni=='Doy' or vni=='doy':
            Yi=Yi-366;
        Y[vni]=Yi
    savez_compressed(fnz,**Y)

#---------------------shpfile--------------------------------------------------
def read_shapefile_data(fname):
    import shapefile as shp
    with shp.Reader(fname) as C:
        #----read shapefile----------------
        S=npz_data();
        S.nrec=C.numRecords
        S.type=C.shapeTypeName

        #----read pts----------------------------------------------------------
        #works for pts and polygon, may not work for other geomerty (need update in this case)
        S.xy=[];
        for i in arange(S.nrec):
            xyi=array(C.shape(i).points);
            parti=array(C.shape(i).parts,dtype='int');
            #insert nan for delimiter
            #to get original index: ind=nonzero(isnan(xyi[:,0]))[0]-arange(len(parti));
            S.xy.append(insert(xyi,parti,nan,axis=0))
        S.xy=squeeze(array(S.xy))

        #---read attributes----------------------------------------------------
        S.attname=array([C.fields[i][0] for i in arange(1,len(C.fields))]);
        stype=array([type(C.record()[m]) for m in S.attname])
        svalue=array(C.records(),dtype='O');
        S.attvalue=array(zeros(len(S.attname))).astype('O')
        for m in arange(len(S.attname)):
            S.attvalue[m]=svalue[:,m].astype(stype[m])
        S.atttype=stype

        #read prj file if exist---
        bdir=os.path.dirname(os.path.abspath(fname));
        bname=os.path.basename(fname).split('.')[0]
        prjname='{}/{}.prj'.format(bdir,bname)
        if os.path.exists(prjname):
            with open(prjname,'r') as fid:
                S.prj=fid.readline().strip()

    return S

def write_shapefile_data(fname,S,float_len=18,float_decimal=8):
    import shapefile as shp
    #---get nrec-----
    if S.type=='POINT':
        if S.xy.dtype==dtype('O'):
            print('S.xy has a dtype="O" for POINT'); sys.exit()
        else:
            nrec=S.xy.shape[0];
    elif S.type=='POLYLINE' or S.type=='POLYGON':
        if S.xy.dtype==dtype('O'):
            nrec=len(S.xy)
        else:
            nrec=1;

    #---check nrec
    if hasattr(S,'nrec'):
        if nrec!=S.nrec:
            print('nrec inconsistent')
            sys.exit()

    #---write shapefile---------
    with shp.Writer(fname) as W:
        W.autoBalance=1;
        #define attributes
        if hasattr(S,'attname'):
            if S.attvalue.dtype==dtype('O'):
                stype=[type(S.attvalue[m][0]) for m in arange(len(S.attname))]
            else:
                if S.attvalue.ndim==1:
                    stype=[type(S.attvalue[0])]
                elif S.attvalue.ndim==2:
                    stype=[type(S.attvalue[m][0]) for m in arange(len(S.attname))]

            for m in arange(len(stype)):
                if stype[m] in [np.int,np.int8,np.int16,np.int32,np.int64]:
                    W.field(S.attname[m],'N')
                elif stype[m] in [np.float,np.float16,np.float32,np.float64]:
                    W.field(S.attname[m],'F',float_len,float_decimal)
                elif stype[m] in [np.str0,np.str,np.str_,np.string_]:
                    W.field(S.attname[m],'C',100)
                else:
                    print('attribute type not included: {}'.format(stype[m]))
                    sys.exit()
        else:
            W.field('field','C')
            W.record()

        #put values
        for i in arange(nrec):
            if S.type=='POINT': #point, W.multipoint(S.xy) is multiple pts features
                vali=S.xy[i]
                W.point(*vali)
            elif S.type=='POLYLINE':
                if S.xy.dtype==dtype('O'):
                    vali=S.xy[i]
                else:
                    vali=S.xy
                #reorganize the shape of vali
                valii=delete_shapefile_nan(vali,0)
                W.line(valii)
            elif S.type=='POLYGON':
                if S.xy.dtype==dtype('O'):
                    vali=S.xy[i]
                else:
                    vali=S.xy
                #reorganize the shape of vali
                valii=delete_shapefile_nan(vali,1)
                W.poly(valii)

            #add attribute
            if hasattr(S,'attname'):
                if S.attvalue.dtype==dtype('O'):
                    atti=[S.attvalue[m][i] for m in arange(len(stype))]
                else:
                    if S.attvalue.ndim==1:
                        atti=[S.attvalue[i]]
                    elif S.attvalue.ndim==2:
                        atti=[S.attvalue[m][i] for m in arange(len(stype))]
                W.record(*atti)

        #----write projection------------
        bname=os.path.basename(fname).split('.')[0]
        bdir=os.path.dirname(os.path.abspath(fname));
        if hasattr(S,'prj'):
            with open('{}/{}.prj'.format(bdir,bname),'w+') as fid:
                fid.write(S.prj)

def delete_shapefile_nan(xi,iloop=0):
    #----delete nan (head and tail), and get ind for the rest
    if xi.ndim==1:
        i1=0; i2=xi.shape[0]
        if isnan(xi[0]): i1=1
        if isnan(xi[-1]): i2=i2-1
        yi=xi[i1:i2]; ind=nonzero(isnan(yi))[0]
    elif xi.ndim==2:
        i1=0; i2=xi.shape[0]
        if isnan(xi[0,0]): i1=1
        if isnan(xi[-1,0]): i2=i2-1
        yi=xi[i1:i2]; ind=nonzero(isnan(yi[:,0]))[0]

    #------reorganize-----------
    if len(ind)==0:
        #close the geomety
        if iloop==1: yi=close_data_loop(yi)

        vi=[yi];
    else:
        vi=[];
        yii=yi[:ind[0]];
        if iloop==1: yii=close_data_loop(yii)
        vi.append(yii)
        for m in arange(len(ind)-1):
            i1=ind[m]+1; i2=ind[m+1];
            yii=yi[i1:i2];
            if iloop==1: yii=close_data_loop(yii);
            vi.append(yii)
        yii=yi[(ind[-1]+1):];
        if iloop==1: yii=close_data_loop(yii)
        vi.append(yii)

    return vi

#---------------------netcdf---------------------------------------------------
def ReadNC(fname,med=1,order=0):
    #ReadNC(fname,med=1)
    #if med=1: return netcdf.Dateset(fname)
    #if med=2: reorgnaized Dataset similar to npz_data
    #order=1: only works for med=2
    #order=0: variable dimension order read not changed for python format
    #order=1: variable dimension order read reversed follwoing in matlab/fortran format
    C=Dataset(fname);

    if med==1:
        return C
    else:

        F=npz_data();
        F.file_format=C.file_format

        #read dims
        ncdims=[i for i in C.dimensions]; 
        F.dimname=ncdims; F.dims=[]; F.dim_unlimited=[]
        for i in ncdims:
            F.dims.append(C.dimensions[i].size)
            F.dim_unlimited.append(C.dimensions[i].isunlimited())
        
        #read attrbutes
        ncattrs=C.ncattrs(); F.attrs=ncattrs
        for i in ncattrs:
            exec('F.{}=C.getncattr("{}")'.format(i,i))
           
        ncvars=[i for i in C.variables]; F.vars=ncvars
        #read variables
        for i in ncvars:
            fi=npz_data();
            dimi=C.variables[i].dimensions;
            fi.dimname=dimi
            fi.dims=[C.dimensions[j].size for j in dimi]
            fi.val=C.variables[i][:]
            fi.attrs=C.variables[i].ncattrs()
            for j in C.variables[i].ncattrs():
                ncattri=C.variables[i].getncattr(j);
                exec('fi.{}=ncattri'.format(j))

            if order==1:
                fi.dimname=list(flipud(fi.dimname))
                fi.dims=list(flipud(fi.dims))
                nm=flipud(arange(ndim(fi.val)));
                fi.val=fi.val.transpose(nm)

            exec('F.{}=fi'.format(i))

        return F

def WriteNC(C,fname,med=1,order=0):
    #WriteNC(C,fname,med=1)
    #C is data source
    #if med=1, C has netcdf.Dataset format
    #if med=2, C has different format
    #order=0: variable dimension order written not changed for python format
    #order=1: variable dimension order written reversed follwoing in matlab/fortran format
    if med==1:
        #----write NC files-------------
        fid=Dataset(fname,'w',format=C.file_format); #C.file_format
        fid.setncattr('file_format',C.file_format)

        #set attrs
        ncattrs=C.ncattrs()
        for i in ncattrs:
            exec("fid.setncattr('{}',C.{})".format(i,i))

        #set dim
        ncdims=[i for i in C.dimensions]
        for i in ncdims:
            if C.dimensions[i].isunlimited() is True:
               fid.createDimension(i,None)
            else:
               fid.createDimension(i,C.dimensions[i].size)

        #set variable
        ncvars=[i for i in C.variables]
        if order==0:
            for vari in ncvars:
                vid=fid.createVariable(vari,C.variables[vari].dtype,C.variables[vari].dimensions)
                for attri in C.variables[vari].ncattrs():
                    vid.setncattr(attri,C.variables[vari].getncattr(attri))
                fid.variables[vari][:]=C.variables[vari][:]
        elif order==1:
            for vari in ncvars:
                vid=fid.createVariable(vari,C.variables[vari].dtype,flipud(C.variables[vari].dimensions))
                for attri in C.variables[vari].ncattrs():
                    vid.setncattr(attri,C.variables[vari].getncattr(attri))
                nm=flipud(arange(ndim(C.variables[vari][:])));
                fid.variables[vari][:]=C.variables[vari][:].transpose(nm)

        fid.close()
    else:
        #----write NC files-------------
        fid=Dataset(fname,'w',format=C.file_format); #C.file_format
       
        #set attrs
        fid.setncattr('file_format',C.file_format)
        if hasattr(C,'attrs'):
           for i in C.attrs:
               exec("fid.setncattr('{}',C.{})".format(i,i))
        
        #set dimension
        for i in range(len(C.dims)):
            if hasattr(C,'dim_unlimited'): 
               dim_flag=C.dim_unlimited[i]
            else:
               dim_flag=False
     
            if dim_flag is True:
               fid.createDimension(C.dimname[i],None)
            else:
               fid.createDimension(C.dimname[i],C.dims[i])

        #set variable 
        if order==0:
            for vari in C.vars:
                vi=eval('C.{}'.format(vari));
                vid=fid.createVariable(vari,vi.val.dtype,vi.dimname)
                for j in vi.attrs:
                    attri=eval('vi.{}'.format(j))
                    vid.setncattr(j,attri)
                fid.variables[vari][:]=vi.val
        elif order==1:
            for vari in C.vars:
                vi=eval('C.{}'.format(vari));
                vid=fid.createVariable(vari,vi.val.dtype,flipud(vi.dimname))
                for j in vi.attrs:
                    attri=eval('vi.{}'.format(j))
                    vid.setncattr(j,attri)
                if ndim(vi.val)>=2:
                    nm=flipud(arange(ndim(vi.val)));
                    fid.variables[vari][:]=vi.val.transpose(nm)
                else:
                    fid.variables[vari][:]=vi.val
        fid.close()

def harmonic_fit(oti,oyi,dt,mti,tidal_const=0):
    #use harmonic analysis to fit time series: used to fill data gap, only for tidal parts.
    #[oti,oyi,dt]: continunous time series
    #mti: time to be interpolated
    #tidal_const=[0,1]: choose tidal const options.
    #         or can be tidal_const=array(['O1','K1','Q1','P1','M2','S2',])
    import platform

    #choose tidal consts
    if tidal_const==0:
        tnames=array(['O1','K1','Q1','P1','M2','S2','K2','N2']);
    elif tidal_const==1:
        tnames=array(['O1','K1','Q1','P1','M2','S2','K2','N2','M3','M4','M6','M7','M8','M10'])
    elif tidal_const==2:
        tnames=array(['SA','SSA','MM','MSF','O1','K1','Q1','P1','M2','S2','K2','N2','M3','M4','M6','M7','M8','M10','N4','S4','S6'])


    #get frequeny consituents
    tfreqs=[]
    sname='.temporary_tidal_consituents_for_HA.dat'
    if platform.system()=='Windows':
        P=loadz(r'D:\Work\Harmonic_Analysis\tide_fac_const.npz')  #directories where database may exist
    elif platform.system()=='Linux':
        P=loadz(r'~/bin/Harmonic_Analysis/tide_fac_const.npz')
    else:
        sys.exit('unknow system')

    for i in arange(len(tnames)):
        ind=nonzero(P.name==tnames[i])[0]
        freqi=P.freq[ind]
        tfreqs.append(freqi)
    tfreqs=squeeze(array(tfreqs))

    with open(sname,'w+') as fid:
        fid.write('{}\n'.format(len(tnames)));
        for i in arange(len(tnames)):
            fid.write('{}\n'.format(tnames[i]));
            fid.write('{}\n'.format(tfreqs[i]));

    #do harmnoic analysis
    H=harmonic_analysis(oyi,dt,tidal_const=sname); os.remove(sname)

    #fit
    mti=mti-oti[0]; myi=mti*0+H.amplitude[0];
    for k in arange(len(tnames)):
        Ai=H.amplitude[k+1]
        Pi=H.phase[k+1]
        Fi=tfreqs[k]*86400
        myi=myi+Ai*cos(Fi*mti-Pi)

    return myi

def harmonic_analysis(yi,dt,StartT=0,executable=None,tidal_const=None):
    #harmonic analyze time series
    #yi: time series; dt: time step (day); StartT: starting time of time series (day)
    #can specified paths of executable and tidal_const separately

    import subprocess
    if executable is None:
        import platform
        #check OS type, and locate the executable
        if platform.system()=='Windows':
            bdirs=[r'D:\Work\Harmonic_Analysis'] #directories where exectuable may exist
            for bdir in bdirs:
                if os.path.exists('{}/tidal_analyze.exe'.format(bdir)): break
            executable='{}/tidal_analyze.exe'.format(bdir)
            if tidal_const is None: tidal_const='{}/tidal_const.dat'.format(bdir)
        elif platform.system()=='Linux':
            bdirs=[r'~/bin/Harmonic_Analysis']
            for bdir in bdirs:
                bdir_ex=os.path.expanduser(bdir)
                if os.path.exists('{}/tidal_analyze'.format(bdir_ex)): break
            executable='{}/tidal_analyze'.format(bdir_ex)
            if tidal_const is None: tidal_const='{}/tidal_const.dat'.format(bdir_ex)
        else:
            print('Operating System unknow: {}'.format(platform.system()))

    #----write time series, HA, and return results--------------------
    fname='temporary_time_series_for_HA.dat'
    sname='temporary_tidal_consituents_for_HA.dat'
    with open(fname,'w+') as fid:
        for i in arange(len(yi)):
            fid.write('{:12.1f} {:12.7f}\n'.format((StartT+i*dt)*86400,yi[i]))
    #print([executable,fname,tidal_const,sname])
    subprocess.call([executable,fname,tidal_const,sname]) #HA
    fid=open(sname,'r'); lines=fid.readlines(); fid.close()
    S=npz_data();
    S.tidal_name=[]; S.amplitude=[]; S.phase=[]
    for line in lines:
        linei=line.split()
        S.tidal_name.append(linei[0])
        S.amplitude.append(linei[1])
        S.phase.append(linei[2])
    S.tidal_name=array(S.tidal_name)
    S.amplitude=array(S.amplitude).astype('float')
    S.phase=array(S.phase).astype('float')

    #clean temporaray files--------
    os.remove(fname); os.remove(sname)

    return S

def get_stat(xi_model,xi_obs):
    #compute statistics between two time series
    #x1, x2 must have the same dimension
    #x1: model; x2: obs
    #import matlab.engine
    #eng=matlab.engine.start_matlab()
    x1=xi_model; x2=xi_obs
    mx1=mean(x1); mx2=mean(x2)

    S=npz_data(); dx=x1-x2; std1=std(x1); std2=std(x2)
    #---save data
    S.R=corrcoef(x1,x2)[0,1] #R
    S.ME=mean(dx)
    S.MAE=mean(abs(dx))
    S.RMSD=sqrt((dx**2).mean())
    S.std=std(dx)
    S.ms=1-sum(dx**2)/sum((abs(x1-mx2)+abs(x2-mx2))**2)
    a,pvalue=sp.stats.pearsonr(x1,x2)
    S.pvalue=pvalue
    S.std1=std1; S.std2=std2
    S.taylor=array([sqrt(mean((x1-x1.mean())**2))/S.std2,sqrt(((x1-x1.mean()-x2+x2.mean())**2).mean())/S.std2,S.R])

    return S

def get_hycom(Time,xyz,vind,hdir='./HYCOM'):
    #extract Hycom time series at stations
    #ti: time seires; xyz=c_[loni,lati,depi];
    #vind: list of index for variables to be extracted. [0,1,2,3] for ['elev','temp','salt','uv']
    #hdir: directory for hycom data

    #variable names
    Var=['surf_el','water_temp','salinity',['water_u','water_v']];
    VarName=['elev','temp','salt',['Ux','Uy']];

    #time
    StartT=floor(Time.min())
    EndT=ceil(Time.max())

    #---interpolation pts----
    if xyz.ndim==1: xyz=xyz[None,:]
    loni=xyz[:,0]; lati=xyz[:,1]
    bxy=c_[lati,loni];
    if xyz.shape[1]==3:
       depi=xyz[:,2]
       bxyz=c_[depi,lati,loni];

    #-----save data----------------
    S=npz_data();
    S.x=xyz[:,0]; S.y=xyz[:,1]
    if xyz.shape[1]==3: S.z=xyz[:,2]

    #---read Hycom data and interpolate onto boundary nodes------------------
    p=npz_data();
    for i in vind: #arange(len(Var)):
        vari=Var[i]; varnamei=VarName[i];

        t0=time.time();
        T0=[]; Data0=[];
        for ti in arange(StartT,EndT+1):
            t1=num2date(ti); t2=num2date(ti+1-1/24/60);

            if isinstance(vari,list):
                fname='Hycom_{}_{}_{}.nc'.format(varnamei[0],t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))
                fname2='Hycom_{}_{}_{}.nc'.format(varnamei[1],t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))
                if not os.path.exists(r'{}/{}'.format(hdir,fname)): continue
                if not os.path.exists(r'{}/{}'.format(hdir,fname2)): continue
                print(fname+'; '+fname2)
                C=ReadNC('{}/{}'.format(hdir,fname),2); C2=ReadNC('{}/{}'.format(hdir,fname2),2)

                #get value
                exec('p.val=C.{}.val'.format(vari[0]))
                exec('p.val2=C2.{}.val'.format(vari[1]))
                p.val=array(p.val); fp=p.val<-29999; p.val2=array(p.val2); fp2=p.val2<-29999;
                p.val[fp]=nan; p.val2[fp2]=nan;
            else:
                fname='Hycom_{}_{}_{}.nc'.format(varnamei,t1.strftime('%Y%m%dZ%H%M00'),t2.strftime('%Y%m%dZ%H%M00'))
                if not os.path.exists(r'{}/{}'.format(hdir,fname)): continue
                print(fname)
                C=ReadNC('{}/{}'.format(hdir,fname),2)

                #get value
                exec('p.val=C.{}.val'.format(vari))
                p.val=array(p.val); fp=p.val<-29999;
                p.val[fp]=nan

            ti=datestr2num(C.time.time_origin)+array(C.time.val)/24
            cloni=array(C.lon.val); cloni=mod(cloni,360)-360
            clati=array(C.lat.val);

            #------define data region extracted
            ind_lon=nonzero((cloni<=max(loni)+0.1)*(cloni>=min(loni)-0.1))[0];
            ind_lat=nonzero((clati<=max(lati)+0.1)*(clati>=min(lati)-0.1))[0];
            i1_lon=ind_lon.min(); i2_lon=i1_lon+len(ind_lon)
            i1_lat=ind_lat.min(); i2_lat=i1_lat+len(ind_lat)

            cloni=cloni[i1_lon:i2_lon]; clati=clati[i1_lat:i2_lat]

            if varnamei=='elev':
                for m in arange(len(ti)):
                    valii=squeeze(p.val[m,i1_lat:i2_lat,i1_lon:i2_lon])

                    #interpolation
                    fd=sp.interpolate.RegularGridInterpolator((clati,cloni),valii,fill_value=nan)
                    vi=fd(bxy)

                    #remove nan pts
                    fp=isnan(vi);
                    if sum(fp)!=0:
                        vi[fp]=sp.interpolate.griddata(bxy[~fp,:],vi[~fp],bxy[fp,:],'nearest')

                    T0.append(ti[m]); Data0.append(vi);
            else:
                #------define data region extracted for depth
                cdepi=array(C.depth.val)
                ind_dep=nonzero((cdepi<=depi.max()+1000)*(cdepi>=depi.min()-100))[0];
                i1_dep=ind_dep.min(); i2_dep=i1_dep+len(ind_dep)
                cdepi=cdepi[i1_dep:i2_dep];

                for m in arange(len(ti)):
                    T0.append(ti[m]);
                    valii=squeeze(p.val[m,i1_dep:i2_dep,i1_lat:i2_lat,i1_lon:i2_lon])

                    #interpolation
                    fd=sp.interpolate.RegularGridInterpolator((cdepi,clati,cloni),valii,fill_value=nan)
                    vi=fd(bxyz)

                    #remove nan pts
                    fp=isnan(vi);
                    if sum(fp)!=0:
                        vi[fp]=sp.interpolate.griddata(bxyz[~fp,:],vi[~fp],bxyz[fp,:],'nearest')

                    #----if variable is velocity
                    if isinstance(varnamei,list):
                        val2ii=squeeze(p.val2[m,i1_dep:i2_dep,i1_lat:i2_lat,i1_lon:i2_lon])

                        #interpolation
                        fd=sp.interpolate.RegularGridInterpolator((cdepi,clati,cloni),val2ii,fill_value=nan)
                        v2i=fd(bxyz)

                        #remove nan pts
                        fp=isnan(v2i);
                        if sum(fp)!=0:
                            v2i[fp]=sp.interpolate.griddata(bxyz[~fp,:],v2i[~fp],bxyz[fp,:],'nearest')

                        Data0.append(r_[expand_dims(vi,0),expand_dims(v2i,0)])

                    else:
                        Data0.append(vi)

        #interpolate in time, and save data
        T0=array(T0); Data0=array(Data0).T; S.time=Time;
        if Data0.ndim==2:
            #interpolate
            datai=[];
            for k in arange(Data0.shape[0]):
                fd=interpolate.interp1d(T0,Data0[k]);
                datai.append(fd(Time));
            datai=array(datai)
            exec('S.{}=datai'.format(varnamei))
        else:
            #interpolate
            datai=[]; data2i=[];
            for k in arange(Data0.shape[0]):
                fd=interpolate.interp1d(T0,Data0[k,0]);
                fd2=interpolate.interp1d(T0,Data0[k,1]);
                datai.append(fd(Time)); data2i.append(fd2(Time))
            datai=array(datai); data2i=array(data2i)
            exec('S.{}=datai'.format(varnamei[0]))
            exec('S.{}=data2i'.format(varnamei[1]))

    return S

#------------------------------------------------------------------------------

if __name__=="__main__":
    pass

##---------get_hycom------------------------------------------------------------
#    #----inputs--------------
#    StartT=datenum(2014,1,1)
#    EndT=datenum(2014,1,4)
#    ti=arange(StartT,EndT,1/24);
#
#    #--convert coordinate if not lon&lat---
#    bp=read_schism_bpfile('Station.bp_WQ')
#    loni,lati=transform(Proj(init='epsg:26919'),Proj(init='epsg:4326'),bp.x,bp.y);
#    xyz=c_[loni,lati,zeros(loni.shape)]
#    xyz=array([-65.91084725,  42.33315317,   0.])
#
#    vind=[0,1,2,3]; #variable index for extract, vind=['elev','temp','salt','uv']
#
#    S=get_hycom(ti,xyz,vind,hdir='./Data');

##-------------HA --------------------------------------------------------------
#    data=loadtxt(r'D:\Work\Harmonic_Analysis\TS.txt');
#    S=harmonic_analysis(data[:,1],1/24,tidal_const=r'D:\Work\Harmonic_Analysis\tmp.dat')
#    S=harmonic_analysis(data[:,1],1/24)

#---------------------projection-----------------------------------------------
#    fname0=r'D:\Work\E3SM\ChesBay\Grid\ChesBay_1a.gr3';
#    fname1=r'D:\Work\E3SM\ChesBay\Grid\ChesBay_1a.ll';
#    proj(fname0,0,'epsg:26918',fname1,3,'epsg:26918')
#    proj(fname0,0,'epsg:26918',fname1,0,'epsg:4326')
#    proj(fname0,3,'epsg:26918',fname1,2,'epsg:4326')

#------------------------------------------------------------------------------
#    clear_globals()
#------------------------------------------------------------------------------
#    xi=arange(100)*0.1;
#    yi=np.random.rand(100)+sin(xi)
#    y1=smooth(yi,3); y2=smooth(yi,5); y3=smooth(yi,7)
#    plot(xi,yi,'r-',xi,y1,'b-',xi,y2,'g-',xi,y3,'k-')
#    show()

#------------------str2num-----------------------------------------------------
#        x='3.5, 4, 5; 5  6.5, 78';
#    x='3.5 4 5 5 6.5   78'
#        xi=str2num(x);
#    x=['3 4 5','4 3 6 8']
#    x=['3 4 5','4 3 6']
#    xi=str2num(x)
#    print(xi)

#     #--test files--
#     fname=r'C:\Users\Zhengui\Desktop\Python\learn\hgrid.gr3'
#     with open(fname,'r') as fid:
#         lines=fid.readlines()
#     line=lines[24535:24545]
#     rline=remove_tail(line)
#     print(line)
#     print(rline)

#-------date_proc---------------------------------------------------------------
#    n1=(2006.,1,1)
#    n2=array([2006,2,1]);
#    n3=array([[2006,3,2,1,1],[2006,3,3,3,0]]);
#
#    f1=datenum(*n1);
#    f2=datenum(n2);
#    f3=datenum(n3,doy=1);

#-------mfft--------------------------------------------------------------------
#    plt.close('all')
#    T=10; dt=0.01; N=T/dt;
#    x=linspace(0.0, T, N);
#    y=4*cos(2.0*pi*(x-0.3)/0.5)+2*cos(2.0*pi*(x-0.4)/1.0)+4*cos(2.0*pi*(x-0.5)/2.0)
#    f,a,p=mfft(y,dt)
#
#    subplot(2,1,1)
#    plot(x,y,'k-')
#    subplot(2,1,2)
#    plot(f,a,'k.',ms=20)
#    setp(gca(),'xlim',[0,5])

#---------------------shpfile--------------------------------------------------

#---------------------netcdf---------------------------------------------------
#    # read NC files
#    C=Dataset('sflux_air_1.002.nc')
#
#    ncdims=[i for i in C.dimensions]
#    ncvars=[i for i in C.variables]
#
#
#    [print("{}".format(i)) for i in ncdims]
#    [print("{}".format(i)) for i in ncvars]
#
#    #----write NC files-------------
#    fid=Dataset('test.nc','w',format='NETCDF3_CLASSIC'); #C.file_format
#
#    for dimi in ncdims:
#        fid.createDimension(dimi,C.dimensions[dimi].size)
#
#    for vari in ncvars:
#        vid=fid.createVariable(vari,C.variables[vari].dtype,C.variables[vari].dimensions)
#        for attri in C.variables[vari].ncattrs():
#           vid.setncattr(attri,C.variables[vari].getncattr(attri))
#        fid.variables[vari][:]=C.variables[vari][:]
#    fid.close()
#
#    ## check results
#    F=Dataset('test.nc');

#----------------convert_matfile_format----------------------------------------
#    fname=r'C:\Users\Zhengui\Desktop\Observation2\USGS\SFBay_USGSData_MAL.mat'
#    fname=r'C:\Users\Zhengui\Desktop\convert_matfile\tem.mat'
#    cmat.convert_matfile_format(fname)
