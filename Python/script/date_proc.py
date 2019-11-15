#!/usr/bin/evn python3
from pylib import *

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



if __name__=='__main__':
    pass
#    n1=(2006.,1,1)
#    n2=array([2006,2,1]);
#    n3=array([[2006,3,2,1,1],[2006,3,3,3,0]]);
#
#    f1=datenum(*n1);
#    f2=datenum(n2);
#    f3=datenum(n3,doy=1);

