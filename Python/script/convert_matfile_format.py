#!/usr/bin/evn python3
from pylib import *

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


if __name__ == "__main__":
    pass
#    fname=r'C:\Users\Zhengui\Desktop\Observation2\USGS\SFBay_USGSData_MAL.mat'
#    fname=r'C:\Users\Zhengui\Desktop\convert_matfile\tem.mat'
#    cmat.convert_matfile_format(fname)
