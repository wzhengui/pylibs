#!/usr/bin/env python3
'''
  script to copy setup 
'''
from pylib import *

#input
base='RUN02m'
runs=['RUN02n']

#output dir
sdir='/sciclone/scr10/wangzg/DSP';

#copy runs
for run in runs:
    if os.path.exists(run): os.system('rm -r {}'.format(run))
    os.system("mkdir {}; cd {}; cp -a ../{}/* ./; rm *.out".format(run,run,base))
    print('copy {} => {}'.format(base,run))
    outdir='{}/{}/outputs'.format(sdir,run)
    os.system("cd {}; rm outputs; mkdir -p {}; ln -sf {} ./".format(run,outdir,outdir))
