#!/usr/bin/env python3
import sys

def read_schism_param(fname,*args):
  with open(fname,'r') as fid:
    lines=fid.readlines()

  param={}
  for line in lines:
    line=line.strip()
    if len(line)==0 or line[0]=='!': continue
    ind=line.find('!');
    if(ind!=-1): line=line[0:ind];
    ind=line.find('=');
    keyi=line[:ind].strip();
    vali=line[(ind+1):].strip();
    param[keyi]=vali
    if((len(args)>0) and (args[0]==1)):
       if vali.lstrip('-').replace('.','',1).isdigit(): param[keyi]=float(vali)
       #try:
       #   param[keyi]=float(vali)
       #except:
       #   pass
       
  return param;

def write_schism_param(fname,param):
    pkeys=sorted(param.keys())
    with open(fname,'w+') as fid:
        #[fid.write('{:10}= {:}\n'.format(x,y)) for x,y in zip(param.keys(),param.values())];    
        #[fid.write('{:10}= {:}\n'.format(i,param[i])) for i in pkeys];    
        for i in range(len(pkeys)):
           fid.write('{:10}= {:}\n'.format(pkeys[i],param[pkeys[i]]))    
    

if __name__=="__main__":
  if len(sys.argv)<2:
    print("not enough arguments !!\nExample: get_schism_param param.in ihot dt\n");
    sys.exit()

  fname=sys.argv[1];
  var=sys.argv[2:];
  param=read_schism_param(fname)

  #print params
  if len(var)!=0:
    for keyi in var:
      if param.get(keyi)==None:
        print('{} : not exist'.format(keyi))
      else:
        print('{} : {}'.format(keyi,param[keyi]))
  else:
    for keyi in param:
      print('{} : {}'.format(keyi,param[keyi]))
