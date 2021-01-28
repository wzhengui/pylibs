#!/usr/bin/evn python3
from pylib import *

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

if __name__=="__main__":
    pass
#        x='3.5, 4, 5; 5  6.5, 78';
#    x='3.5 4 5 5 6.5   78'
#        xi=str2num(x);
#    x=['3 4 5','4 3 6 8']
#    x=['3 4 5','4 3 6']
#    xi=str2num(x)
#    print(xi)

##-----test files----------
#     fname=r'C:\Users\Zhengui\Desktop\Python\learn\hgrid.gr3'
#     with open(fname,'r') as fid:
#         lines=fid.readlines()
#     line=lines[24535:24545]
#     rline=remove_tail(line)
#     print(line)
#     print(rline)
