#!/usr/bin/env python3
#usage:
#   backup.py: compare Python directory with backup directory bk
#   backup.py 1: backup files into bk
#   backup.py 2: restore files from bk
#   backup.py 1 viz3: rsync files from local to viz3
#   backup.py 2 viz3: rsync files from viz3 to local

import os, sys

args=sys.argv

if len(args)==1: #only compare
    print('-'*100)
    os.system('Ddiff bk ./; Ddiff bk/script script')
    print('\n'+'-'*100+'\n#usage:')
    print('#   backup.py: compare Python directory with backup directory bk')
    print('#   backup.py 1: backup files into bk')
    print('#   backup.py 2: restore files from bk')
    print('#   backup.py 1 viz3: rsync files from local to viz3')
    print('#   backup.py 2 viz3: rsync files from viz3 to local')
    print('-'*100)
elif len(args)==2: #backup or restore in local')
    flag=int(args[1])
    if flag==1:
        os.system('cp *.* ./bk; cp script/*.* ./bk/script/')
    elif flag==2: #restore
        #backup files in case
        bdir='.backup'
        if not os.path.exists(bdir): os.mkdir(bdir)
        os.system('cp *.* {}; cp -r script {}'.format(bdir,bdir))

        #restore
        os.system('cp ./bk/*.* ./; cp ./bk/script/*.* ./script/')
    else:
        print('unknow flag: {}'.format(flag))
elif len(args)==3: # rsync with viz3
    flag=int(args[1]); rname=args[2]
    if rname=='viz3':
        viz3=os.environ['viz3']
        if flag==1: #rsync from local to viz3
            os.system("rsync -ra ./*.* {}/bin/Python/; cd script; rsync -ra ./*.py {}/bin/Python/script".format(viz3,viz3));
        elif flag==2: #rsync from viz3 to local
            #backup files in case
            bdir='.backup'
            if not os.path.exists(bdir): os.mkdir(bdir)
            os.system('cp *.* {}; cp -r script {}'.format(bdir,bdir))

            #rsync
            os.system("rsync -ra {}/bin/Python/*.* ./; cd script; rsync -ra {}/bin/Python/script/*.py ./".format(viz3,viz3));
