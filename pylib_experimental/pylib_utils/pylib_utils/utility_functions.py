#!/usr/bin/evn python3
from numpy import array, ones, zeros, arange, tile, r_, nan, isnan, unique, where, \
                    concatenate, append, delete, savetxt, loadtxt, linspace, \
                    meshgrid, arctan2, pi, cos, sin, sqrt, exp, log, log10, \
                    sum, mean, std, min, max, abs, round, floor, ceil, \
                    squeeze, transpose, flipud, fliplr, rot90, roll, copy, \
                    int8, int16, int32, int64, float16, float32, float64, \
                    setdiff1d, load, angle, argsort, nonzero,\
                    mod, c_, argmax, real
from numpy.linalg import inv
from scipy.fftpack import fft, ifft

import numpy as np
import pickle
import sys
import os
from copy import deepcopy as dcopy
import scipy as sp
import matplotlib as mpl
from netCDF4 import Dataset
from pyproj import Transformer

def close_data_loop(xi):
    '''
    constructe data loop along the first dimension.
    if xi[0,...]!=xi[-1,...], then,add xi[0,...] in the end
    '''
    if array_equal(xi[0,...].ravel(),xi[-1,...].ravel()):
        vi=xi
    else:
        vi=r_[xi,xi[0,...][None,...]]
    return vi

#-------loadz------------------------------------------------------------------
def get_VINFO(data):
    '''
    collect information about object's attributes
    '''
    atts=[]; sdict=data.__dict__; skeys=sdict.keys(); fnc=0
    if ('dimname' in skeys) and ('dims' in skeys) and ('file_format' in skeys): fnc=1 #netcdf
    stypes=[int,int8,int16,int32,int64, float,float16,float32,float64]
    snames=['int','int8','int16','int32','int64','float','float16','float32','float64']
    for i in skeys:
        try:
            vi=sdict[i]; dt=type(vi); dta=''
            if (fnc==1) and (dt is zdata): vi=vi.val; dt=np.ndarray #netcdf file
            #get data information
            if dt is list:
               nd=': list({},)'.format(len(vi))
               if (fnc==1) and (i in ['dimname','dims']): nd=': list{}'.format(vi) #netcdf file
            elif dt is dict:
               nd=': dict({},)'.format(len(vi))
            elif dt is str:
               nd=': "{}", string'.format(vi[:30])
            elif dt is np.ndarray:
               nd=': array{}'.format(vi.shape); dta=str(vi.dtype)
               if (fnc==1) and (i in ['dimname','dims']): nd=': array({})'.format(vi) #netcdf file
               if vi.size==1 and (vi.dtype in stypes): nd=': {}, array({})'.format(squeeze(vi),1)
            elif dt in stypes:
               nd=': {}, {} '.format(vi,snames[stypes.index(dt)])
            else:
               nd=': {}'.format(type(vi))

            #output
            ms=min(6,max([len(k) for k in skeys])); fs1='{:'+str(ms)+'s}{}'; fs2=fs1+', {}'
            fstr=fs2.format(i,nd,dta) if dta!='' else fs1.format(i,nd)
            atts.append(fstr.strip())
        except:
            pass
    return atts

class zdata:
    '''
    self-defined data structure by Zhengui Wang.  Attributes are used to store data
    '''
    def __init__(self):
        pass

    @property
    def VINFO(self):
        return get_VINFO(self)

def savez(fname,data,fmt=0):
    '''
    save data as self-defined python format
       fmt=0: save data as *.npz (small filesize, reads slower)
       fmt=1: save data as *.pkl (large filesize, reads faster)
    if fname endswith *.npz or *.pkl, then fmt is reset to match fname
    '''

    #determine format
    if fname.endswith('.npz'): fmt=0; fname=fname[:-4]
    if fname.endswith('.pkl'): fmt=1; fname=fname[:-4]
    if fmt==1: fname=fname+'.pkl'

    #save data
    if fmt==0:
       #get all attribute
       svars=list(data.__dict__.keys())
       if 'VINFO' in svars: svars.remove('VINFO')

       #check whether there are functions. If yes, change function to string
       rvars=[]
       for svar in svars:
           if hasattr(data.__dict__[svar], '__call__'):
              import cloudpickle
              try:
                 data.__dict__[svar]=cloudpickle.dumps(data.__dict__[svar])
              except:
                 print('function {} not saved'.format(svar))
                 rvars.append(svar)
       svars=setdiff1d(svars,rvars)

       #check variable types for list,string,int and float
       lvars=[]; tvars=[]; ivars=[]; fvars=[]
       for svar in svars:
           if isinstance(data.__dict__[svar],int): ivars.append(svar)
           if isinstance(data.__dict__[svar],float): fvars.append(svar)
           if isinstance(data.__dict__[svar],str): tvars.append(svar)
           if isinstance(data.__dict__[svar],list):
              lvars.append(svar)
              data.__dict__[svar]=array(data.__dict__[svar],dtype='O')

       #constrcut save_string
       save_str='savez_compressed("{}" '.format(fname)
       for svar in svars: save_str=save_str+',{}=data.{}'.format(svar,svar)
       save_str=save_str+',_list_variables=lvars,_str_variables=tvars,_int_variables=ivars,_float_variables=fvars)'
       exec(save_str)
    elif fmt==1:
       fid=open(fname,'wb'); pickle.dump(data,fid,pickle.HIGHEST_PROTOCOL); fid.close()

def loadz(fname,svars=None):
    '''
    load self-defined data "fname.npz" or "fname.pkl"
         svars: list of variables to be read
    '''

    if fname.endswith('.npz'):
       #get data info
       data0=load(fname,allow_pickle=True)
       keys0=data0.keys() if svars is None else svars
       ivars=list(data0['_int_variables']) if ('_int_variables' in keys0) else []
       fvars=list(data0['_float_variables']) if ('_float_variables' in keys0) else []
       tvars=list(data0['_str_variables']) if ('_str_variables' in keys0) else []
       lvars=list(data0['_list_variables']) if ('_list_variables' in keys0) else []

       #extract data
       vdata=zdata()
       for keyi in keys0:
           if keyi in ['_list_variables','_str_variables','_int_variables','_float_variables']: continue

           #get variable
           datai=data0[keyi]
           if datai.dtype==dtype('O'): datai=datai[()] #restore object
           if keyi in ivars: datai=int(datai)          #restore int variable
           if keyi in fvars: datai=float(datai)        #restore int variable
           if keyi in tvars: datai=str(datai)          #restore str variable
           if keyi in lvars: datai=datai.tolist()      #restore list variable

           #if value is a function
           if 'cloudpickle.cloudpickle' in str(datai):
              import pickle
              try:
                 datai=pickle.loads(datai)
              except:
                 continue
           vdata.__dict__[keyi]=datai
    elif fname.endswith('.pkl'):
       import pickle
       vdata=zdata(); fid=open(fname,'rb')
       data=pickle.load(fid)
       vdata.__dict__=dcopy(data).__dict__.copy()
       fid.close()
    else:
       sys.exit('unknown format: {}'.format(fname))
    return vdata

def least_square_fit(X,Y):
    '''
    perform least square fit
    usage: CC,fit_data=least_square_fit(X,Y)
        X: data maxtrix (npt,n)
        Y: data to be fitted (npt)
    where CC(n) is coefficient, fit_data is data fitted
    '''
    if X.shape[0]!=len(Y): X=X.T
    CC=inv(X.T@X)@(X.T@Y); fy=X@CC

    return [CC,fy]

#-------mfft-------------------------------------------------------------------
def mfft(xi,dt):
    '''
    Perform FFT for a time series, with a time interval specified

    usage: period,afx,pfx=mfft(xi,dt)
    input:
       xi: time series
       dt: time interval

    output:
       period[period],afx[amplitude],pfx[phase]
    '''
    N=xi.size;
    fx=fft(xi);
    afx=abs(fx[1:N//2])*2.0/N;
    pfx=angle(fx[1:N//2]);
    period=dt*N/arange(1,N//2);
    return period,afx,pfx

def command_outputs(code,shell=True):
    '''
    Capture the command output from the system

    usage:
         S=command_outputs('ls')
         print(S.stderr)
         print(S.stdout) # normally this is the results
    '''
    import subprocess
    p=subprocess.Popen(code,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=shell)
    stdout,stderr=p.communicate()
    if stdout!=None: stdout=stdout.decode('utf-8')
    if stderr!=None: stderr=stderr.decode('utf-8')

    out=zdata()
    out.stdout=stdout
    out.stderr=stderr
    return out

def near_pts(pts,pts0,method=0,N=100):
    '''
    return index of pts0 that pts is nearest
       usage: sind=near_pts(pts,pts0)
       pts0: c_[x0,y0];  pts: c_[x,y]
       algorithm: using sp.spatial.cKDTree (default)

    old methods: method=1 and  method=2
       usage: sind=near_pts(pts,pts0,method=1[2],N=100)
       pts[n,2]: xy of points
       pts0[n,2]: xy of points

       method=1: quick method by subgroups (N);
       method=2: slower methods
    '''

    if method==0:
        sind=sp.spatial.cKDTree(pts0).query(pts)[1]
    elif method==1:
        p=pts[:,0]+(1j)*pts[:,1]
        p0=pts0[:,0]+(1j)*pts0[:,1]

        # N=min(N,len(p)-1);
        # #--divide pts into subgroups based on dist, each group size is about N---------
        # ps0=[]; ps=[]; ds=[]; inds=[]
        # i0=0; mval=1e50; pflag=p==None
        # while(True):
        #     ps0.append(p[i0]);
        #     dist=abs(p-p[i0]); dist[pflag]=mval;
        #     ind=argsort(dist); sdist=dist[ind]; mds=sdist[N]; i0=ind[N]
        #     fp=dist<mds; psi=p[fp]; dsi=max(dist[fp]); pflag[fp]=True;
        #     ps.append(psi); ds.append(dsi); inds.append(nonzero(fp)[0])
        #     if mds==mval: break

        N=min(N,len(p));
        #--divide pts into subgroups based on dist, each group size is about N---------
        ps0=[]; ps=[]; ds=[]; inds=[]; inum=arange(len(p))
        while(True):
            if len(inum)==0: break
            dist=abs(p[inum]-p[inum[0]]); sind=argsort(dist); inum=inum[sind]; dist=dist[sind]; sN=min(N,len(inum))
            ps0.append(p[inum[0]]); ps.append(p[inum[:sN]]); ds.append(dist[sN-1]); inds.append(inum[:sN])
            inum=inum[sN:]

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
    elif method==2:
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

def inside_polygon(pts,px,py,fmt=0,method=0):
    '''
    check whether points are inside polygons

    usage: sind=inside_polygon(pts,px,py):
       pts[npt,2]: xy of points
       px[npt] or px[npt,nploy]: x coordiations of polygons
       py[npt] or py[npt,nploy]: y coordiations of polygons
       (npt is number of points, nploy is number of polygons)

       fmt=0: return flags "index[npt,nploy]" for whether points are inside polygons (1 means Yes, 0 means No)
       fmt=1: only return the indices "index[npt]" of polygons that pts resides in
             (if a point is inside multiple polygons, only one indice is returned; -1 mean pt outside of all Polygons)

       method=0: use mpl.path.Path
       method=1: use ray method explicitly

    note: For method=1, the computation time is proportional to npt**2 of the polygons. If the geometry
          of polygons are too complex, dividing them to subregions will increase efficiency.
    '''
    #----use ray method-----------------------------
    #check dimension
    if px.ndim==1:
       px=px[:,None]; py=py[:,None]

    #get dimensions
    npt=pts.shape[0]; nv,npy=px.shape

    if nv==3 and fmt==1:
       #for triangles, and only return indices of polygons that points resides in
       px1=px.min(axis=0); px2=px.max(axis=0); py1=py.min(axis=0); py2=py.max(axis=0)

       sind=[];
       for i in arange(npt):
           pxi=pts[i,0]; pyi=pts[i,1]
           sindp=nonzero((pxi>=px1)*(pxi<=px2)*(pyi>=py1)*(pyi<=py2))[0]; npy=len(sindp)
           if npy==0:
               sind.append(-1)
           else:
               isum=ones(npy)
               for m in arange(nv):
                   xi=c_[ones(npy)*pxi,px[m,sindp],px[mod(m+1,nv),sindp]]
                   yi=c_[ones(npy)*pyi,py[m,sindp],py[mod(m+1,nv),sindp]]
                   area=signa(xi,yi)
                   fp=area<0; isum[fp]=0;
               sindi=nonzero(isum!=0)[0]

               if len(sindi)==0:
                   sind.append(-1)
               else:
                   sind.append(sindp[sindi[0]])
       sind=array(sind)

    else:
        if method==0:
            sind=[]
            for m in arange(npy):
                sindi=mpl.path.Path(c_[px[:,m],py[:,m]]).contains_points(pts)
                sind.append(sindi)
            sind=array(sind).T+0  #convert logical to int
        elif method==1  :
            #using ray method explicitly
            sind=ones([npt,npy])
            x1=pts[:,0][:,None]; y1=pts[:,1][:,None]
            for m in arange(nv):
                x2=px[m,:][None,:]; y2=py[m,:][None,:]; isum=zeros([npt,npy])
                # sign_x1_x2=sign(x1-x2)
                for n in arange(1,nv-1):
                    x3=px[(n+m)%nv,:][None,:]; y3=py[(n+m)%nv,:][None,:]
                    x4=px[(n+m+1)%nv,:][None,:]; y4=py[(n+m+1)%nv,:][None,:]

                    #formulation for a ray to intersect with a line
                    fp1=((y1-y3)*(x2-x1)+(x3-x1)*(y2-y1))*((y1-y4)*(x2-x1)+(x4-x1)*(y2-y1))<=0 #intersection inside line P3P4
                    fp2=((y2-y1)*(x4-x3)-(y4-y3)*(x2-x1))*((y4-y3)*(x1-x3)+(y3-y1)*(x4-x3))<=0 #P1, P2 are in the same side of P3P4

                    fp12=fp1*fp2
                    isum[fp12]=isum[fp12]+1
                fp=((isum%2)==0)|((x1==x2)*(y1==y2))
                sind[fp]=0

        #change format
        if fmt==1:
            sindm=argmax(sind,axis=1)
            sindm[sind[arange(npt),sindm]==0]=-1
            sind=sindm
        elif fmt==0 and npy==1:
            sind=sind[:,0]
    return sind

def mdist(xy1,xy2,fmt=0,outfmt=0):
    '''
      compute distance or (minimum distance) between geometry (pts,lines) to geomtry(pts,lines,polygons)
      format of xy:
        points: c_[x,y]
        lines:  c_[x,y], or c_[x1,y1,x2,y2]
        polygon: c_[x,y];  if polygon is not closed, (x[0],y[0]) will be added in the end
      fmt: choose geometry types
          0: pts and pts;      1: pts and lines;       2: lines and lines
          3: pts and polygon;  4: lines and polygon;   5: polygon and polygon
      outfmt=0: return an array of minimum distance;  outfmt=1: return an matrix of distance
    '''
    def _reshape_lines(lxy):
        '''
          if lxy=c_[x1,y1,x2,y2]: each row is a line (x1,y1) -> (x2,y2)
          if lxy=c_[px,py]: lines are obtained with (px[i],px[i]) -> (px[i+1],px[i+1])
        '''
        if lxy.shape[1]==2:
            lxy=c_[lxy[:-1],lxy[1:]]
        elif lxy.shape[1]!=4:
            sys.exit('unknown shape of lines: lxy={}'.format(lxy.shape))
        return lxy

    def pts_pts(xy,xy0,outfmt=0):
        '''
        find the minimum distance of c_[x,y] to a set of points c_[x0,y0]
          outfmt=0: minimum distance of pts to all pts (return an array)
          outfmt=1: minimum distance of pts to each pts (return a matrix)
        '''
        if outfmt==0:
            pid=near_pts(xy,xy0); dist=abs((xy[:,0]+1j*xy[:,1])-(xy0[pid,0]+1j*xy0[pid,1]))
        else:
            dist=abs((xy[:,0]+1j*xy[:,1])[:,None]-(xy0[pid,0]+1j*xy0[pid,1])[None,:])
        return dist

    def pts_lines(xy,lxy,outfmt=0):
        '''
        find the minimum distance of c_[x,y] to a set of lines lxy
          1). lxy=c_[x1,y1,x2,y2]: each row is a line (x1,y1) -> (x2,y2)
          2). lxy=c_[px,py]: lines are obtained with (px[i],px[i]) -> (px[i+1],px[i+1])
          outfmt=0: minimum distance of pts to all lines (return an array)
          outfmt=1: minimum distance of pts to each lines (return a matrix)
        '''

        #reshape pts and lines
        x,y=xy.T[:,:,None]; x1,y1,x2,y2=_reshape_lines(lxy).T[:,None,:]

        dist=nan*ones([x.size,x1.size])
        #compute the foot of a perpendicular, and distance between pts
        k=-((x1-x)*(x2-x1)+(y1-y)*(y2-y1))/((x1-x2)**2+(y1-y2)**2); xn=k*(x2-x1)+x1; yn=k*(y2-y1)+y1
        fpn=(k>=0)*(k<=1); dist[fpn]=abs((x+1j*y)-(xn+1j*yn))[fpn] #normal line
        dist[~fpn]=array(r_[abs((x+1j*y)-(x1+1j*y1))[None,...], abs((x+1j*y)-(x2+1j*y2))[None,...]]).min(axis=0)[~fpn] #pt-pt dist
        if outfmt==0: dist=dist.min(axis=1)
        return dist

    def lines_lines(lxy,lxy0,outfmt=0):
        '''
        find the minimum distance of lines lxy to a set of lines lxy0
          1). lxy or lxy0=c_[x1,y1,x2,y2]: each row is a line (x1,y1) -> (x2,y2)
          2). lxy or lxy0=c_[px,py]: lines are obtained with (px[i],px[i]) -> (px[i+1],px[i+1])
        outfmt=0: minimum distance of lines lxy to all lines lxy0 (return an array)
        outfmt=1: minimum distance of lines lxy to each lines in lxy0 (return a matrix)
        Note: if lines are intersecting, the minimum distance is zero.
        '''
        #reshape lines
        x1,y1,x2,y2=_reshape_lines(lxy).T[:,:,None]; x3,y3,x4,y4=_reshape_lines(lxy0).T[:,None,:]

        #check whether intersects
        fc1=(((x1-x3)*(y4-y3)+(x3-x4)*(y1-y3))*((x2-x3)*(y4-y3)+(x3-x4)*(y2-y3)))<=0
        fc2=(((x3-x1)*(y2-y1)+(x1-x2)*(y3-y1))*((x4-x1)*(y2-y1)+(x1-x2)*(y4-y1)))<=0
        fpn=(fc1*fc2); x1,y1,x2,y2=x1[:,0],y1[:,0],x2[:,0],y2[:,0]; x3,y3,x4,y4=x3[0],y3[0],x4[0],y4[0]

        if outfmt==1:
            dist=nan*ones([x1.size,x3.size]); dist[fpn]=0

            print('ZG', x1.shape,x3.shape)
            if sum(~fpn)!=0:
                #check distance between pts to lines
                dist[~fpn]=array([pts_lines(c_[x1,y1],c_[x3,y3,x4,y4],1),pts_lines(c_[x2,y2],c_[x3,y3,x4,y4],1),
                   pts_lines(c_[x3,y3],c_[x1,y1,x2,y2],1).T,pts_lines(c_[x4,y4],c_[x1,y1,x2,y2],1).T]).min(axis=0)[~fpn]
        else:
            dist=nan*ones(x1.size); fpn=fpn.sum(axis=1)>0; dist[fpn]=0
            x1,y1,x2,y2=x1[~fpn],y1[~fpn],x2[~fpn],y2[~fpn]
            if sum(~fpn)!=0:
                dist[~fpn]=array([pts_lines(c_[x1,y1],c_[x3,y3,x4,y4]),pts_lines(c_[x2,y2],c_[x3,y3,x4,y4]),
                   pts_lines(c_[x3,y3],c_[x1,y1,x2,y2],1).min(axis=0),pts_lines(c_[x4,y4],c_[x1,y1,x2,y2],1).min(axis=0)]).min(axis=0)
        return dist

    def pts_polygon(xy,pxy,outfmt=0):
        '''
        find the minimum distance of c_[x,y] to a polygon pxy
        Note: If pts are inside the polygon, the distance (minimum distance) is zero
        '''
        pxy=close_data_loop(pxy); fpn=inside_polygon(xy, pxy[:,0],pxy[:,1])==1

        if outfmt==0:
            dist=nan*ones(len(xy)); dist[fpn]=0
            if sum(~fpn)!=0: dist[~fpn]=pts_lines(xy[~fpn],pxy)
        else:
            dist=nan*ones([len(xy),len(pxy)-1]); dist[fpn]=0
            if sum(~fpn)!=0: dist[~fpn]=pts_lines(xy[~fpn],pxy,1)
        return dist

    def lines_polygon(lxy,pxy,outfmt=0):
        '''
        find the minimum distance of lines lxy to a polygon pxy
        Note: if pts of lines are inside the polygon, the minimum distance is zero.
        '''
        x1,y1,x2,y2=_reshape_lines(lxy).T;  pxy=close_data_loop(pxy)

        #check pts inside the polygon
        fpn1=inside_polygon(c_[x1,y1],pxy[:,0],pxy[:,1])==1
        fpn2=inside_polygon(c_[x2,y2],pxy[:,0],pxy[:,1])==1
        fpn=fpn1|fpn2

        dist=nan*ones(x1.size) if outfmt==0 else nan*ones([x1.size,len(pxy)-1])
        dist[fpn]=0; x1,y1,x2,y2=x1[~fpn],y1[~fpn],x2[~fpn],y2[~fpn]
        dist[~fpn]=lines_lines(c_[x1,y1,x2,y2],pxy,outfmt)
        return dist

    if fmt==0: dist=pts_pts(xy1,xy2,outfmt)
    if fmt==1: dist=pts_lines(xy1,xy2,outfmt)
    if fmt==2: dist=lines_lines(xy1,xy2,outfmt)
    if fmt==3: dist=pts_polygon(xy1,xy2,outfmt)
    if fmt==4: dist=lines_polygon(xy1,xy2,outfmt)
    return dist

def signa(x,y):
    '''
        compute signed area for triangles along the last dimension (x[...,0:3],y[...,0:3])
    '''
    if x.ndim==1:
        area=((x[0]-x[2])*(y[1]-y[2])-(x[1]-x[2])*(y[0]-y[2]))/2
    elif x.ndim==2:
        # area=((x[:,0]-x[:,2])*(y[:,1]-y[:,2])-(x[:,1]-x[:,2])*(y[:,0]-y[:,2]))/2
        area=((x[...,0]-x[...,2])*(y[...,1]-y[...,2])-(x[...,1]-x[...,2])*(y[...,0]-y[...,2]))/2
    area=np.squeeze(area)
    return area

def mdivide(A,B):
    '''
        perform matrix division B/A
    '''
    if A.ndim==1: A=A[None,:]
    if B.ndim==1: B=B[None,:]
    A2=inv(A@A.T);
    B2=B@A.T
    return squeeze(B2@A2)

def bpfilt(data,delta_t,band_f):
    '''
    band pass filter for 1D (data[time]) or nD (data[time,...]) array along the first dimension
       delta_t: time interval
       band_f: frequency band (eg. band_f=[12,48])
    '''

    #data dimension
    ds=data.shape; N=ds[0]

    #remove mean, do fft
    data=data-data.mean(axis=0)[None,...]
    fdata=fft(data,axis=0)

    #design filter
    filt=ones(N)
    k1=int(floor(band_f[0]*N*delta_t))-1; k1=max(min(k1,N//2),1)
    k2=int(ceil(band_f[1]*N*delta_t))+1;  k2=max(min(k2,N//2),1)
    filt[:k1]=0; filt[-k1:]=0; filt[k2:N-k2]=0
    filt=(ones([*ones(data.ndim).astype('int')])*filt).T

    #remove low and high freqs
    fdata=fdata*filt

    #bp results
    bpdata=real(ifft(fdata,axis=0))

    return bpdata

def lpfilt(data,delta_t,cutoff_f):
    '''
    low pass filter for 1D (data[time]) or nD (data[time,...]) array along the first dimension

    Note: there is no phase shift for this LP-filter
    '''
    #import gc #discard

    ds=data.shape

    #fft original data
    mdata=data.mean(axis=0)[None,...]
    #print(data.shape,mdata.shape)
    data=data-mdata
    fdata=fft(data,axis=0)

    #desgin filter
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

    #expand dimension of filt
    filt=(ones([*ones(data.ndim).astype('int')])*filt).T

    #remove high freqs
    fdata=fdata*filt

    #lp results
    lfdata=real(ifft(fdata,axis=0))+mdata

    return lfdata

def smooth(xi,N):
    '''
    smooth average (on the 1st dimension):
       xi[time,...]: time series
       N: window size (if N is even, then N=N+1)
    '''

    #window span
    if mod(N,2)==0: N=N+1
    nz=int((N-1)/2)

    #moving averaging
    X=xi.copy(); SN=N*ones(xi.shape)
    for nmove in arange(1,nz+1):
        #sum in the beginning and end
        X[nmove:,...]=X[nmove:,...]+xi[:-nmove,...]
        X[:-nmove,...]=X[:-nmove,...]+xi[nmove:,...]

        #count
        SN[:nmove,...]=SN[:nmove,...]-1
        SN[-nmove:]=SN[-nmove:]-1
    SX=X/SN
    return SX

def daytime_length(lat,doy):
    '''
    calculate daytime length based on latitutde and day_of_year
    lat: latitude, doy: (1-365), sunrise=12-daytimelength/2, sunset=12+daytimelength/2
    '''
    P=arcsin(0.39795*cos(0.2163108 + 2*arctan(0.9671396*tan(0.00860*(doy-186)))));
    T=(sin(0.8333*pi/180)+sin(lat*pi/180)*sin(P))/cos(lat*pi/180)/cos(P);
    dt=24-(24/pi)*arccos(T);
    return dt

def move_figure(x=0,y=0,f=None):
    '''
    Move figure (f=gcf()) to upper left corner to pixel (x, y)
    e.g. move_figure(0,0,gcf())
    '''
    if f is None: f=gcf()
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)

def proj(fname0=None,fmt0=None,prj0=None,fname1=None,fmt1=None,prj1=None,x=None,y=None,lon0=None,lat0=None,order0=0,order1=0):
    '''
    tranfrom projection of files: proj(fname0,fmt0,prj0,fname1,fmt1,prj1,x,y,lon0,lat0,order0,order1)
       fname: file name
       fmt: 0: SCHISM gr3 file; 1: SCHISM bp file; 2: xyz file; 3: xyz file with line number
       prj: projection name (e.g. 'epsg:26918', 'epsg:4326','cpp'), or projection string (e.g. prj0=get_prj_file('epsg:4326'))
       lon0,lat0: center for transformation between lon&lat and x&y; if lon0=None or lat0=None, then x0=mean(lon), and y0=mean(lat)
       order=0 for projected coordinate; order=1 for lat&lon (if prj='epsg:4326', order=1 is applied automatically')

    tranform data directly: px,py=proj(prj0='epsg:26918',prj1='epsg:4326',x=px,y=py)

    function used:
          lat,lon=Transformer.from_crs('epsg:26918','epsg:4326').transform(x,y)
          x,y=Transformer.from_crs('epsg:4326','epsg:26918').transform(lat,lon)
    #x1,y1=transform(Proj(proj0),Proj(proj1),x,y); #not used anymore
    '''

    from pylib_essentials.schism_file import read_schism_hgrid,read_schism_bpfile,schism_bpfile

    #read file
    if fmt0==0:
        gd=read_schism_hgrid(fname0)
        np=gd.np; x=gd.x; y=gd.y; z=gd.dp
    elif fmt0==1:
        gd=read_schism_bpfile(fname0)
        np=gd.nsta; x=gd.x; y=gd.y; z=gd.z
    elif fmt0==2:
        gd=loadtxt(fname0)
        np=gd.shape[0]; x=gd[:,0]; y=gd[:,1]; z=gd[:,2]
    elif fmt0==3:
        gd=loadtxt(fname0)
        np=gd.shape[0]; x=gd[:,1]; y=gd[:,2]; z=gd[:,3]
    else:
        if x is None or y is None: sys.exit('unknown format of input files: {}, {}'.format(fname0,fname1))

    #order of data
    prj0=prj0.lower(); prj1=prj1.lower(); icpp=0
    if prj0=='epsg:4326': order0=1
    if prj1=='epsg:4326': order1=1
    if 'cpp' in [prj0,prj1]: icpp=1; rearth=6378206.4

    #transform coordinate
    if icpp==1 and prj0=='cpp':
        if prj1!='epsg:4326': sys.exit('projection wrong: prj0={}, prj1={}'.format(prj0,prj1))
        if (lon0 is None) or (lat0 is None): sys.exit('need lon0 and lat0 for cpp=>ll transform')
        x1=lon0+x*180/(pi*rearth*cos(lat0*pi/180)); y1=y*180/(pi*rearth)
    elif icpp==1 and prj1=='cpp':
        if prj0!='epsg:4326': sys.exit('projection wrong: prj0={}, prj1={}'.format(prj0,prj1))
        if lon0 is None: lon0=mean(x)
        if lat0 is None: lat0=mean(y)
        x1=rearth*(x-lon0)*(pi/180)*cos(lat0*pi/180); y1=rearth*y*pi/180
    else:
        if order0==1: x,y=y,x
        fpn=~(isnan(x)|isnan(y)); x1=arange(len(x))*nan; y1=arange(len(y))*nan
        x1[fpn],y1[fpn]=Transformer.from_crs(prj0,prj1).transform(x[fpn],y[fpn])
        if order1==1: x1,y1=y1,x1
        if (sum(isnan(x1[fpn]))!=0) | (sum(isnan(y1[fpn]))!=0): sys.exit('nan found in tranformation: x1,y1') #check nan

    # if (sum(isnan(x))!=0) and (sum(isnan(y))!=0): sys.exit('nan found in x,y') #check nan
    #x1,y1=Transformer.from_crs(prj0,prj1).transform(x,y)
    #x1,y1=transform(Proj(init=proj0),Proj(init=proj1),x,y);
    #x1,y1=transform(Proj(proj0),Proj(proj1),x,y);

    #write file
    if fmt1==0:
        if fmt0!=0: sys.exit('{} should have gr3 format'.format(fname0))
        gd.x=x1; gd.y=y1;
        gd.write_hgrid(fname1,Info='!coordinate transformed: {}=>{}'.format(prj0,prj1))
    elif fmt1==1:
        if fmt0==1:
            gd.x=x1; gd.y=y1
        else:
            gd=schism_bpfile()
            gd.note='coordinate transformed: {}=>{}'.format(prj0,prj1)
            gd.nsta=np; gd.x=x1; gd.y=y1; gd.z=z;
        gd.write_bpfile(fname1)
    elif fmt1==2 or fmt1==3:
        with open(fname1,'w+') as fid:
            for i in arange(np):
                if fmt1==2: fid.write('{} {} {}\n'.format(x1[i],y1[i],z[i]))
                if fmt1==3: fid.write('{} {} {} {}\n'.format(i+1,x1[i],y1[i],z[i]))
    else:
       return [x1,y1]

def proj_pts(x,y,prj1='epsg:4326',prj2='epsg:26918'):
    '''
    convert projection of points from prj1 to prj2
      x,y: coordinate of pts
      prj1: name of original projection
      prj2: name of target projection
    '''
    px,py=proj(prj0=prj1,prj1=prj2,x=x,y=y)
    return [px,py]

def get_prj_file(prjname='epsg:4326',fmt=0,prj_dir=r'D:\Work\Database\projection\prj_files'):
    '''
    return projection name or entire database (dict)
        fmt=0: get one projection
        fmt=1: return dict of projection database
        fmt=-1: process *.prj files from prj_dir

    #-------online method-----------------
    #function to generate .prj file information using spatialreference.org
    #def getWKT_PRJ (epsg_code):
    #     import urllib
    #     # access projection information
    #     wkt = urllib.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code))
    #     # remove spaces between charachters
    #     remove_spaces = wkt.read().replace(" ","")
    #     # place all the text on one line
    #     output = remove_spaces.replace("\n", "")
    #     return output
    '''

    #get location of database
    bdir=os.path.dirname(__file__)

    if fmt==0:
        S=loadz('{}/../pyScripts/prj.npz'.format(bdir))
        return S.prj[prjname]
    elif fmt==1:
        S=loadz('{}/../pyScripts/prj.npz'.format(bdir))
        return S.prj
    elif fmt==-1:
        #dir of *.prj files
        #prj_dir=r'D:\Work\Database\projection\prj_files'

        #---processing all prj files-------------------------------------------
        fnames=os.listdir(prj_dir)

        prj=dict()
        for fname in fnames:
            import re
            R=re.match('epsg.(\d+).prj',fname);
            if not R: continue
            prj_num=int(R.groups()[0])
            with open('{}/{}'.format(prj_dir,fname),'r') as fid:
                line=fid.readline().strip();
                prj['epsg:{}'.format(prj_num)]=line

        # #save prj file as a database
        # S=zdata();
        # S.prj=prj
        # savez('{}/prj'.format(bdir),S)
        return prj
    else:
        sys.exit('unknow fmt')

#----------------convert_matfile_format----------------------------------------
#def convert_matfile_format(file):
#    '''
#    convert a matlab data to self-defined python format
#    file: input, a directory or a matfile
#    '''
#    # eg. file="C:\Users\Zhengui\Desktop\Observation2\DWR\SFBay_DWRData_SSI.mat"
#    # file="D:\OneDrive\Python\tem.mat"
#
#    fname=[]
#    if os.path.isdir(file):
#        sfile=os.listdir(file)
#        for sfilei in sfile:
#            ename=sfilei.rstrip().split('.')[-1]
#            if ename=='mat':
#                fname.append(file+os.sep+sfilei)
#    else:
#        fname=[file];
#
#    fid=open('log_matfile_convert.txt','w+')
#    #convert mat format from v7.3 to v7
#    cdir=os.getcwd();
#    os.chdir('D:\OneDrive\Python')
#    import matlab.engine
#    eng = matlab.engine.start_matlab()
#
#    for fn in fname:
#        print('converting matfile: '+fn)
#        dname=os.path.dirname(fn)
#        bname=os.path.basename(fn).split('.')[0]
#        fnv7=dname+os.sep+bname+'_v7'
#        fnz=dname+os.sep+bname
#        eflag=eng.convert_matfile_format(fn,fnv7,nargout=1)
#        if eflag!=0:
#            print('convert flag is %d: %s\n' % (eflag, fn));
#            fid.write('convert flag is %d: %s\n' % (eflag,fn))
#            continue
#        convert_matfile(fnz,fnv7)
#        os.remove(fnv7+'.mat')
#    fid.close()
#    os.chdir(cdir)
#

#convert between MATLAB matfile and zdata
def npz2mat(npz_data,fname):
    '''
      convert Python *.npz file/data to Matlab *mat file
      Exmaples
        1. npz2mat('test.npz','test.mat')
        2. npz2mat(C,'test.mat') #C=loadz('test.npz')
    '''
    from scipy import io

    if isinstance(npz_data,str): npz_data=loadz(npz_data)
    sp.io.savemat(fname,npz_data.__dict__)

def convert_matfile(matfile,fname=None):
    '''
    Convert Matlab *.mat file to Python *.npz file (alias: mat2npz)
    Inputs:
      matfile: name of matlab file ('name.mat' or 'name')
      fname: name of *.npz file to be saved  ('name.npz' or 'name')

    Examples: (mat2npz is alias to convert_matfile)
      1. mat2npz('A.mat','A.npz')
    '''
    from scipy import io

    #check name
    if matfile.endswith('.mat'): matfile=matfile[:-4]

    #read matfile and convert
    C=sp.io.loadmat(matfile+'.mat',simplify_cells=True); S=zdata()
    for x in C.keys():
       if x.startswith('__') or x=='VINFO':continue

       #for list of strings
       if isinstance(C[x],np.ndarray):
          if C[x].size!=1 and ('<U' in str(C[x].dtype)):
             C[x]=array([i.strip() for i in C[x]])
       S.__dict__[x]=C[x]
    if fname is not None: savez(fname,S)
    return S

def cindex(index,shape):
    '''
    convert array index: same as unravel_index and ravel_multi_index
      1) cindex(id, ds):  convert flat index (1D) to indices (nD)
      2) cindex(c_[id1,id2,...], ds):  convert flat indices (nD) to index (1D)
    '''
    index=array(index)
    if index.ndim!=1 and index.shape[0]!=len(shape): index=index.T

    if index.ndim==1:
       cid=[*unravel_index(index,shape)]
    else:
       cid=ravel_multi_index(index,shape)
    return cid

def EOF(data,npc=8,scale=0,center=False,**args):
    '''
    EOF analysis
        data(time,...): data with 1st dimension as time, other dimes for space
        npc: number of principal components (PCs) returned
        scale=0: normalized by eigenvalue; scale=1: normalized by the value of time series of each PC.
        Note: the mean value is not removed (center=False) in the analysis at default.
    Outputs: (PC,CC,VC,solver)
        PC: principal components
        CC: coefficients of each PC
        VC: variation of each PC
        solver: EOF analysis solver, and all results can be derived from it (e.g., solver.reconstructedField(8))
    '''
    from eofs.standard import Eof
    solver=Eof(data,center=center,**args)
    PC=solver.eofs(neofs=npc,eofscaling=2)
    CC=solver.pcs(npcs=npc,pcscaling=1).T
    VC=solver.varianceFraction(npc)

    #normalize
    for m,cc in enumerate(CC):
        if cc.mean()<0: CC[m],PC[m]=-CC[m],-PC[m]
        if scale==1: rat=abs(cc).mean(); CC[m], PC[m]=CC[m]/rat,PC[m]*rat
    return PC, CC, VC, solver

def get_stat(xi_model,xi_obs,fmt=0):
    '''
    compute statistics between two time series
    x1, x2 must have the same dimension
    x1: model; x2: obs
    fmt=1: compute pvalue using scipy.stats.pearsonr

    #import matlab.engine
    #eng=matlab.engine.start_matlab()
    '''

    x1=xi_model; x2=xi_obs
    mx1=mean(x1); mx2=mean(x2)

    S=zdata(); dx=x1-x2; std1=std(x1); std2=std(x2)
    #---save data
    S.R=corrcoef(x1,x2)[0,1] #R
    S.ME=mean(dx)
    S.MAE=mean(abs(dx))
    S.RMSD=sqrt((dx**2).mean())
    S.std=std(dx)
    S.ms=1-sum(dx**2)/sum((abs(x1-mx2)+abs(x2-mx2))**2)
    if fmt==1: a,S.pvalue=sp.stats.pearsonr(x1,x2)
    S.std1=std1; S.std2=std2
    S.taylor=array([sqrt(mean((x1-x1.mean())**2))/S.std2,sqrt(((x1-x1.mean()-x2+x2.mean())**2).mean())/S.std2,S.R])

    return S

#---------------------shpfile--------------------------------------------------
def read_shapefile_data(fname):
    '''
    read shapefile (*.shp) and return its content

    note:  works for pts and polygon only, may not work for other geomerties (need update in these cases)
    '''

    import shapefile as shp
    with shp.Reader(fname) as C:
        #----read shapefile----------------
        S=zdata();
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
        S.xy=squeeze(array(S.xy,dtype='O'))

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

def write_shapefile_data(fname,data,prj='epsg:4326',float_len=18,float_decimal=8):
    '''
    write shapefiles
      fname: file name
      data: data in zdata() format
      prj: projection

    example: 
       S=zdata(); S.type=='POINT'; S.prj='epsg:4326'
       S.xy=c_[slon[:],slat[:]]
       S.attname=['station']; S.attvalue=station[:]
       write_shp('test.shp',S)

    note: only works for geometry: POINT, POLYLINE, POLYGON
    '''

    import shapefile as shp
    S=data
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
    else:
        sys.exit('unknow type')
    if hasattr(S,'nrec'):
        if nrec!=S.nrec: sys.exit('nrec inconsistent')

    #---write shapefile---------
    with shp.Writer(fname) as W:
        W.autoBalance=1;
        #define attributes
        if hasattr(S,'attname'):
            if not hasattr(S.attvalue,'dtype'): S.attvalue=array(S.attvalue)
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
                valii=[[[*k] for k in valii[0]]]
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
        bdir=os.path.dirname(os.path.abspath(fname))
        prjs=S.prj if hasattr(S,'prj') else prj
        if len(prjs)<20: prjs=get_prj_file(prjs)
        fid=open('{}/{}.prj'.format(bdir,bname),'w+'); fid.write(prjs); fid.close()

def delete_shapefile_nan(xi,iloop=0):
    '''
    delete nan (head and tail), and get ind for the rest
    '''
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
def ReadNC(fname,fmt=0,mode='r',order=0):
    '''
    read netcdf files, and return its values and attributes
        fname: file name

        fmt=0: reorgnaized Dataset with format of zdata
        fmt=1: return netcdf.Dateset(fname)
        fmt=2: reorgnaized Dataset (*npz format), ignore attributes

        mode='r+': change variable in situ. other choices are: 'r','w','a'

        order=1: only works for med=2; change dimension order
        order=0: variable dimension order read not changed for python format
        order=1: variable dimension order read reversed follwoing in matlab/fortran format
    '''

    #get file handle
    C=Dataset(fname,mode=mode);

    if fmt==1:
        return C
    elif fmt in [0,2]:
        F=zdata(); F.file_format=C.file_format

        #read dims
        ncdims=[i for i in C.dimensions];
        F.dimname=ncdims; F.dims=[]; F.dim_unlimited=[]
        for i in ncdims:
            F.dims.append(C.dimensions[i].size)
            F.dim_unlimited.append(C.dimensions[i].isunlimited())

        #read attrbutes
        ncattrs=C.ncattrs(); F.attrs=ncattrs
        for i in ncattrs: F.__dict__[i]=C.getncattr(i)

        ncvars=[*C.variables]; F.vars=array(ncvars)
        #read variables
        for i in ncvars:
            if fmt==0:
               vi=zdata()
               dimi=C.variables[i].dimensions
               vi.dimname=dimi
               vi.dims=[C.dimensions[j].size for j in dimi]
               vi.val=C.variables[i][:]
               vi.attrs=C.variables[i].ncattrs()
               for j in C.variables[i].ncattrs():
                   ncattri=C.variables[i].getncattr(j)
                   vi.__dict__[j]=ncattri

               if order==1:
                   vi.dimname=list(flipud(vi.dimname))
                   vi.dims=list(flipud(vi.dims))
                   nm=flipud(arange(ndim(vi.val)))
                   vi.val=vi.val.transpose(nm)
            elif fmt==2:
               vi=array(C.variables[i][:])
               if order==1: vi=vi.transpose(flip(arange(vi.ndim)))
            F.__dict__[i]=vi

        C.close()
        return F
    else:
        sys.exit('wrong fmt')

def WriteNC(fname,data,fmt=0,order=0):
    '''
    write zdata to netcdf file
        fname: file name
        data:  soure data (can be zdata or capsule of zdatas)
        fmt=0, data has zdata() format
        fmt=1, data has netcdf.Dataset format
        order=0: variable dimension order written not changed for python format
        order=1: variable dimension order written reversed follwoing in matlab/fortran format
    '''

    if fmt==0: #write NC files for zdata
        #--------------------------------------------------------------------
        # pre-processing zdata format
        #--------------------------------------------------------------------
        if not fname.endswith('.nc'): fname=fname+'.nc'
        if not hasattr(data,'file_format'): data.file_format='NETCDF4'
        S=data.__dict__
        #check variables
        if not hasattr(data,'vars'):
           svars=[i for i,vi in S.items() if isinstance(vi,zdata) or isinstance(vi,np.ndarray)] 
        else:
           svars=data.vars
        for svar in svars:  #change np.array to zdata
            if isinstance(S[svar],np.ndarray): fvi=zdata(); fvi.val=S[svar]; S[svar]=fvi

        #check global dimensions
        if not hasattr(data,'dims'):  
           dims=[]; [dims.extend(S[svar].val.shape) for svar in svars]; dims=unique(array(dims)) 
           dimname=['dim_{}'.format(i) for i in dims]
        else:
           dims,dimname=array(data.dims),data.dimname
        if (not hasattr(data,'dim_unlimited')) or (len(data.dim_unlimited)!=len(data.dims)):
           dim_unlimited=[False for i in dims]
        else: 
           dim_unlimited=data.dim_unlimited

        #check variable dimensions and attributes
        for svar in svars:
            if not hasattr(S[svar],'dimname'): S[svar].dimname=[dimname[nonzero(dims==i)[0][0]] for i in S[svar].val.shape]
            if not hasattr(S[svar],'attrs'): S[svar].attrs=[i for i in S[svar].__dict__ if i!='val']

        #check global attribute
        if not hasattr(data,'attrs'):  
           attrs=[i for i in S if i not in svars]
        else:
           attrs=data.attrs

        #write zdata as netcdf file
        fid=Dataset(fname,'w',format=data.file_format)   #open file
        fid.setncattr('file_format',data.file_format)    #set file_format
        for i in attrs: fid.setncattr(i,S[i])            #set attrs
        for ds,dn,dF in zip(dims,dimname,dim_unlimited): #set dimension
            fid.createDimension(dn,None) if (dF is True) else fid.createDimension(dn,ds)
        for svar in svars:  #set variable
            vi=S[svar]; nm=vi.val.ndim; vdm=vi.dimname
            vid=fid.createVariable(svar,vi.val.dtype,vdm) if order==0 else fid.createVariable(svar,vi.val.dtype,vdm[::-1])
            if hasattr(vi,'attrs'):
               [vid.setncattr(i,vi.__dict__[i]) if i!='_FillValue' else vid.setncattr('_fillvalue',vi.__dict__[i]) for i in vi.attrs]
            fid.variables[svar][:]=vi.val if order==0 else vi.val.T
        fid.close()
    elif fmt==1: #write netcdf.Dateset as netcdf file
        C=data; cdict=C.__dict__; cdim=C.dimensions; cvar=C.variables
        fid=Dataset(fname,'w',format=C.file_format)     #open file
        fid.setncattr('file_format',C.file_format)      #set file_format
        for i in C.ncattrs(): fid.setncattr(i,cdict[i]) #set attrs
        for i in cdim: fid.createDimension(i,None) if cdim[i].isunlimited() else fid.createDimension(i,cdim[i].size) #set dims
        for i in cvar: #set vars
            vdim=cvar[i].dimensions
            vid=fid.createVariable(i,cvar[i].dtype,vdim) if order==0 else fid.createVariable(i,cvar[i].dtype,vdim[::-1])
            for k in cvar[i].ncattrs(): vid.setncattr(k,cvar[i].getncattr(k))
            fid.variables[i][:]=cvar[i][:] if order==0 else cvar[i][:].T
        fid.close()

def read(fname,*args0,**args):
    '''
    generic function in read files with standard suffixs
    suffix supported:  npz, gr3, ll, ic, vgrid.in, bp, reg, prop, xlsx, nc, shp, nml, asc, tif, tiff
    '''
    from .schism_file import (read_schism_hgrid, read_schism_vgrid, read_schism_bpfile, read_schism_reg,
                              read_schism_prop, read_schism_param)

    #determine read function
    F=None
    if fname.endswith('.asc') or fname.endswith('.tif') or fname.endswith('.tiff'): F=read_dem
    if fname.endswith('.gr3') or fname.endswith('.ll') or fname.endswith('.ic'): F=read_schism_hgrid
    if fname.endswith('vgrid.in'):  F=read_schism_vgrid
    if fname.endswith('.npz'):  F=loadz
    if fname.endswith('.bp'):   F=read_schism_bpfile
    if fname.endswith('.reg'):  F=read_schism_reg
    if fname.endswith('.prop'): F=read_schism_prop
    if fname.endswith('.xlsx'): F=read_excel
    if fname.endswith('.nc'):   F=ReadNC
    if fname.endswith('.shp'):  F=read_shapefile_data
    if fname.endswith('.nml'):  F=read_schism_param
    if F is None: sys.exit('unknown type of file: '+fname)

    return F(fname,*args0,**args)

def harmonic_fit(oti,oyi,dt,mti,tidal_names=0, **args):
    '''
    construct time series using harmonics: used to fill data gap, or extrapolation
       [oti,oyi,dt]: observed (time, value, time interval)
       mti: time to interpolate/extrapolate
       tidal_names=0/1/2: tidal names options; tidal_names=['O1','K1','Q1','P1','M2','S2',]: list of tidal consituent names
    '''
    #choose tidal consts
    if not hasattr(tidal_names,'__len__'):
       if tidal_names==0:
          tidal_names=array(['O1','K1','Q1','P1','M2','S2','K2','N2'])
       elif tidal_names==1:
          tidal_names=array(['O1','K1','Q1','P1','M2','S2','K2','N2','M3','M4','M6','M7','M8','M10'])
       elif tidal_names==2:
          tidal_names=array(['SA','SSA','MM','MSF','O1','K1','Q1','P1','M2','S2','K2','N2','M3','M4','M6','M7','M8','M10','N4','S4','S6'])

    #do harmnoic analysis
    H=harmonic_analysis(oyi,dt,tidal_names=tidal_names,**args)

    #fit data
    myi=zeros(mti.shape)+H.amplitude[0]
    for k,tname in enumerate(tidal_names): myi=myi+H.amplitude[k+1]*cos(H.freq[k+1]*(mti-oti[0])*86400-H.phase[k+1])

    return myi

def harmonic_analysis(data,dt,t0=0,tidal_names=None,code=None,tname=None,fname=None,sname=None):
    '''
    harmonic analyze time series
        data: time series
        dt: time step (day)
        t0: starting time of time series (day); normally use datenum (e.g. datenum(2010,0,0))
        tidal_names: 1) path of tidal_const.dat, or 2) names of tidal consituents
        code: path of executable (tidal_analyze, or tidal_analyze.exe) for HA analysis
        [tname,fname,sname]: temporary names for tidal_const, time series, and HA results
    '''
    import subprocess

    #names of temporary file
    if tname is None: tname='.temporary_tidal_const_for_HA.dat'
    if fname is None: fname='.temporary_time_series_for_HA.dat'
    if sname is None: sname='.temporary_tidal_consituents_for_HA.dat'

    #check OS type, and locate the executable
    sdir='{}/../pyScripts/Harmonic_Analysis'.format(os.path.dirname(__file__))
    if code is None:
        import platform
        #directories where exectuable may exist
        bdirs=[sdir, r'D:\Work\Harmonic_Analysis',r'~/bin/Harmonic_Analysis']
        if platform.system().lower()=='windows':
           code=['{}/tidal_analyze.exe'.format(i) for i in bdirs if os.path.exists('{}/tidal_analyze.exe'.format(i))][0]
        elif platform.system().lower()=='linux':
           code=['{}/tidal_analyze'.format(i) for i in bdirs if os.path.exists('{}/tidal_analyze'.format(i))][0]
        else:
            print('Operating System unknow: {}'.format(platform.system()))
        if code is None: sys.exit('exectuable "tidal_analyze" was not found')

    #locate or write tidal_const.dat
    if tidal_names is None:
       tidal_names='{}/tidal_const.dat'.format(sdir)
    else:
       if not isinstance(tidal_names,str): #names of tidal consituents
          C=loadz('{}/tide_fac_const.npz'.format(sdir)); tdict=dict(zip(C.name,C.freq))
          fid=open(tname,'w+'); fid.write('{}\n'.format(len(tidal_names)))
          for i in tidal_names: fid.write('{}\n{}\n'.format(i.upper(),tdict[i.upper()]))
          fid.close(); tidal_names=tname

    #write time series, HA, and return results
    fid=open(fname,'w+')
    for i,datai in enumerate(data): fid.write('{:12.1f} {:12.7f}\n'.format((t0+i*dt)*86400,datai))
    fid.close(); subprocess.call([code,fname,tidal_names,sname]) #HA

    #save results
    S=zdata(); S.tidal_name,S.amplitude,S.phase=array([i.strip().split() for i in open(sname,'r').readlines()]).T
    S.amplitude=array(S.amplitude).astype('float'); S.phase=array(S.phase).astype('float')
    S.freq=r_[0,array([i.strip() for i in open(tidal_names,'r').readlines()[2::2]]).astype('float')]

    #clean temporaray files--------
    os.remove(fname); os.remove(sname)
    if os.path.exists(tname): os.remove(tname)
    return S

if __name__=="__main__":
    pass
