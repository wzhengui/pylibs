#!/usr/bin/evn python3
from pylib import *

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
    data0=load(fname)
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


if __name__=='__main__':
    pass
    fname='D:\Work\SFBay\Observation\CMON\DWR\SFBay_DWRData_Turb.npz'
    T=loadz(fname)
    wipe()

