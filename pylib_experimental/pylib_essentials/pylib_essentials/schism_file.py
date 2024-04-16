#/usr/bin/env python3
import os
import sys
from numpy import array, array_equal, argsort, arange, c_, r_, diff,\
  sort, unique, zeros, ones, setdiff1d, fliplr, flipud, tile, nonzero,\
  nan, isnan, loadtxt, savetxt, load, linspace, meshgrid, concatenate,\
  squeeze, mod, angle, pi, sign, sqrt, maximum
from matplotlib.pyplot import figure, show, savefig, tricontour,\
    tricontourf, tripcolor, colorbar, gca, gcf, plot, axis, xlim, ylim,\
    triplot, title, xlabel, ylabel, gca, gcf, setp, getp, close
import matplotlib.cm as cm
import matplotlib as mpl

from pylib_essentials.utility_functions import loadz, savez, zdata, signa, \
    get_VINFO, proj, WriteNC, inside_polygon, near_pts

# for experimental features
from pathlib import Path
import numpy as np
import pickle
import pandas as pd
import copy
import unittest


class schism_grid:
    def __init__(self, fname=None):
        '''
        Initialize to empty instance if fname is not provided;
        otherwise, read from three supported file format
        '''
        if fname is None:
            pass
        elif fname.endswith('gr3') or fname.endswith('.ll'):
            self.read_hgrid(fname)
        elif fname.endswith('.pkl') or fname.endswith('.npz'):
            self.__dict__=loadz(fname).hgrid.__dict__.copy()
        else:
            raise Exception('hgrid file format {} not recognized'.format(fname))
        self.source_file = fname

    @property
    def INFO(self):
        return get_INFO(self)
    @property
    def VINFO(self):
        return get_INFO(self)

    @property
    def z(self):
        return self.dp

    @property
    def zctr(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return self.dpe

    @property
    def zcj(self):
        if not hasattr(self,'zcj'): self.compute_side(2)
        return self.dps

    def plot(self,ax=None,fmt=0,value=None,ec=None,fc=None,lw=0.1,levels=None,shading='gouraud',xy=0,
             ticks=None,xlim=None,ylim=None,clim=None,extend='both',method=0,cb=True,cb_aspect=30,cb_pad=0.02,**args):
        '''
        plot grid with default color value (grid depth)
        fmt=0: plot grid only; fmt=1: plot filled contours
        fmt=2: plot contour lines at levels; colors and linewidths can be provided for each contour
        value: color value size(np,or ne)
        ec: color of grid line;  fc: element color; lw: grid line width
        levels=100: number of colors for depths; levels=array([v1,v2,...]): depths for plot
        ticks=[v1,v2,...]: colorbar ticks; ticks=10: number of ticks
        clim=[vmin,vmax]: value range for plot/colorbar
        method=0: using tricontourf/tricontouf; method=1: using tripcolor
        cb=False: not add colorbar
        cb_aspect: adjust colorbar width
        shading: only used for method=1, and value.size=gd.np
        xy=0: plot with gd.x,gd.y;  xy=1: use gd.lon,gd.lat;  xy=c_[x,y]: use provided xy coordinates
        '''

        if ec is None: ec='None'
        if fc is None: fc='None'
        if levels is None: levels=51
        if ax is None: ax=gca()
        x,y=[self.x,self.y] if xy==0 else [self.lon,self.lat] if xy==1 else xy.T
        fp3=self.i34==3; fp4=~fp3; vm=clim

        if fmt in [1,2]: #plot contours
           trs=r_[self.elnode[:,:3],c_[self.elnode[fp4,0],self.elnode[fp4,2:]]]
           if value is None: value=self.dp
           if vm is None: fpn=~isnan(value); vm=[min(value[fpn]),max(value[fpn])]
           if vm[0]==vm[1] or (vm[1]-vm[0])/(abs(vm[0])+abs(vm[1]))<1e-10: vm[1]=vm[1]+max((vm[1]-vm[0])*1e-10,1e-10)

           #plot
           if fmt==1 and method==1:  #tripcolor
              if value.size==self.np: hg=tripcolor(x,y,trs,value,vmin=vm[0],vmax=vm[1],shading=shading,**args)
              if value.size==self.ne: hg=tripcolor(x,y,trs,facecolors=r_[value,value[fp4]],vmin=vm[0],vmax=vm[1],**args)
              if value.size==self.ne+sum(fp4) and sum(fp4)!=0: hg=tripcolor(x,y,trs,facecolors=value,vmin=vm[0],vmax=vm[1],**args)
           else:  #contourf or contour
              if sum(isnan(value))!=0: trs=trs[~isnan(value[trs].sum(axis=1))] #set mask
              if value.size==self.ne: value=self.interp_elem_to_node(value=value) #elem value to node value
              if not hasattr(levels,'__len__'): levels=linspace(*vm,int(levels)) #detemine levels
              if fmt==1: hg=tricontourf(x,y,trs,value,levels=levels,vmin=vm[0],vmax=vm[1],extend=extend,**args)
              if fmt==2: hg=tricontour(x,y,trs,value,levels=levels, vmin=vm[0],vmax=vm[1],extend=extend,**args)

           #add colobar
           cm.ScalarMappable.set_clim(hg,vmin=vm[0],vmax=vm[1])
           if cb==True:
              hc=colorbar(hg,aspect=cb_aspect,pad=cb_pad); self.hc=hc
              if ticks is not None:
                 if not hasattr(ticks,'__len__'):
                    hc.set_ticks(linspace(*vm,int(ticks)))
                 else:
                    hc.set_ticks(ticks)

        if (fmt==0)|(ec!='None'): #plot grid
           if ec=='None': ec=['k','k']
           if isinstance(ec,str): ec=[ec,ec]
           if not hasattr(lw,'__len__'): lw=[lw,lw*0.75]
           iqd=self.elnode[fp4]; iqd=c_[iqd,iqd[:,0],tile(0,len(iqd))].ravel()
           x3,y3=x[iqd],y[iqd]; x3[5::6]=nan; y3[5::6]=nan
           hg0=[triplot(x,y,self.elnode[fp3,:3],lw=lw[0],color=ec[0],**args), plot(x3,y3,lw=lw[1],color=ec[1],**args)]

        hg=hg0 if fmt==0 else hg if ec=='None' else [*hg0,hg]; self.hg=hg
        if xlim is not None: setp(ax,xlim=xlim)
        if ylim is not None: setp(ax,ylim=ylim)
        if mpl.get_backend().lower() in ['qt5agg','qtagg']:
           acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
           if 'bp' not in ats: self.bp=schism_bpfile()
           if 'reg' not in ats: self.reg=schism_bpfile(fmt=1)
           if 'query' not in ats: self.query_pt()
           if 'bnd' not in ats: self.create_bnd()
           if 'node' not in ats: self.show_node()
        return hg

        #-------------------------------------------------
        #for reference: old grid plot method
        #-------------------------------------------------
        #elif gridplot==1:
        #   if not hasattr(self,'isidenode'): self.compute_side(fmt=1)
        #   x3=c_[self.x[self.isidenode],nan*zeros(self.ns)].ravel(); x4=[]
        #   y3=c_[self.y[self.isidenode],nan*zeros(self.ns)].ravel(); y4=[]
        #else:
        #   #tri
        #   tri=self.elnode[fp3,:3]; tri=c_[tri,tri[:,0]]
        #   x3=self.x[tri]; y3=self.y[tri]
        #   x3=c_[x3,ones([sum(fp3),1])*nan]; x3=reshape(x3,x3.size)
        #   y3=c_[y3,ones([sum(fp3),1])*nan]; y3=reshape(y3,y3.size)
        #   #quad
        #   quad=self.elnode[fp4,:]; quad=c_[quad,quad[:,0]]
        #   x4=self.x[quad]; y4=self.y[quad]
        #   x4=c_[x4,ones([sum(fp4),1])*nan]; x4=reshape(x4,x4.size)
        #   y4=c_[y4,ones([sum(fp4),1])*nan]; y4=reshape(y4,y4.size)

        #-------------------------------------------------
        #for reference: old method in plot contourf using PolyCollection
        #-------------------------------------------------
        #elif method==1:
        #   #creat polygon
        #   xy4=c_[self.x,self.y][self.elnode];
        #   xy4=array([s[0:-1,:] if (i34==3 and len(s)==4) else s for s,i34 in zip(xy4,self.i34)])

        #   #elem value
        #   if value is None:
        #      if not hasattr(self,'dpe'): self.compute_ctr()
        #      value=self.dpe
        #   else:
        #      if len(value)==self.np:
        #         value=self.interp_node_to_elem(value=value)
        #      elif len(value)!=self.ne:
        #         sys.exit('value has wrong size: {}'.format(value.shape))

        #   # apply mask
        #   if mask is not None: xy4=xy4[mask]; value=value[mask]
        #      #ind=nonzero(mask)[0]
        #      #xy4=xy4[ind];
        #      #value=value[ind]

        #   #get clim
        #   if clim is None: clim=[min(value),max(value)]

        #   #plot
        #   if fmt==0:
        #      hg=mpl.collections.PolyCollection(xy4,lw=lw,edgecolor=ec,facecolor=fc,antialiased=False,**args)
        #   else:
        #      hg=mpl.collections.PolyCollection(xy4,lw=lw,edgecolor=ec,array=value,clim=clim,antialiased=False,**args)

        #      #add colorbar
        #      if cb:
        #         hc=colorbar(hg); self.hc=hc;
        #         if ticks is not None: hc.set_ticks(ticks)
        #         hc.set_clim(clim);

        #   #add to figure
        #   ax.add_collection(hg)
        #   ax.autoscale_view()
        #   self.hg=hg; #show(block=False)
        #   if xlim is not None: setp(ax,xlim=xlim)
        #   if ylim is not None: setp(ax,ylim=ylim)
        #   return hg
    def plot_grid(self,**args):
        '''
        alias for gd.plot()
        '''
        return self.plot(**args)

    def plot_bnd(self,c='k',lw=0.5,ax=None,xy=0,**args):
        '''
          plot schims grid boundary
          xy=0: plot with gd.x,gd.y;  xy=1: use gd.lon,gd.lat;  xy=c_[x,y]: use provided xy coordinates
          gd.plot_bnd(): plot bnd
          gd.plot_bnd(c='rb'): open bnd in red,land bnd in blue
        '''
        if ax!=None: sca(ax)
        x,y=[self.x,self.y] if xy==0 else [self.lon,self.lat] if xy==1 else xy.T
        if not hasattr(self,'nob'): self.compute_bnd()

        #get indices for bnds
        sindo=[]
        for i in arange(self.nob):
            sindo=r_[sindo,-1,self.iobn[i]]
        sindo=array(sindo).astype('int'); fpn=sindo==-1
        bx1=x[sindo]; by1=y[sindo]
        bx1[fpn]=nan; by1[fpn]=nan

        sindl=[]
        for i in arange(self.nlb):
            if self.island[i]==0:
               sindl=r_[sindl,-1,self.ilbn[i]]
            else:
               sindl=r_[sindl,-1,self.ilbn[i],self.ilbn[i][0]]
        sindl=array(sindl).astype('int'); fpn=sindl==-1
        bx2=x[sindl]; by2=y[sindl]
        bx2[fpn]=nan; by2[fpn]=nan

        if len(c)==1:
           hb=plot(r_[bx1,nan,bx2],r_[by1,nan,by2],c,lw=lw,**args); self.hb=hb
        else:
          hb1=plot(bx1,by1,c[0],lw=lw,**args); hb2=plot(bx2,by2,c[-1],lw=lw,**args); self.hb=[hb1,hb2]
        if mpl.get_backend().lower() in ['qt5agg','qtagg']:
           acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
           if 'bp' not in ats: self.bp=schism_bpfile()
           if 'reg' not in ats: self.reg=schism_bpfile(fmt=1)
           if 'query' not in ats: self.query_pt()
           if 'bnd' not in ats: self.create_bnd()
           if 'node' not in ats: self.show_node()
        return self.hb

    def read_hgrid(self,fname,*args):
        #attribute tracking the file originally read, mainly used for savez and save_pkl
        self.source_file = fname

        fid=open(fname,'r'); lines=fid.readlines(); fid.close()

        #read ne and np; lx,ly and dp
        self.ne,self.np=[*array(lines[1].split()[0:2]).astype('int')]
        self.x,self.y,self.dp=array([i.split()[1:4] for i in lines[2:(2+self.np)]]).astype('float').T
        if len(lines)<(2+self.np+self.ne): return

        #read elnode and i34
        fdata=[i.strip().split() for i in lines[(2+self.np):(2+self.np+self.ne)]]
        fdata=array([i if len(i)==6 else [*i,'-1'] for i in fdata]).astype('int')
        self.i34=fdata[:,1]; self.elnode=fdata[:,2:]-1; fdata=None

        #compute ns
        self.compute_side()
        if len(lines)<(4+self.np+self.ne): return

        #read open bnd info
        n=2+self.np+self.ne; self.nob=int(lines[n].strip().split()[0]); n=n+2; self.nobn=[]; self.iobn=[]
        for i in arange(self.nob):
            self.nobn.append(int(lines[n].strip().split()[0]))
            self.iobn.append(array([int(lines[n+1+k].strip().split()[0])-1 for k in arange(self.nobn[-1])]))
            n=n+1+self.nobn[-1]
        self.nobn=array(self.nobn); self.iobn=array(self.iobn,dtype='O')
        if len(self.iobn)==1: self.iobn=self.iobn.astype('int')

        #read land bnd info
        self.nlb=int(lines[n].strip().split()[0]); n=n+2; self.nlbn=[]; self.ilbn=[]; self.island=[]
        for i in arange(self.nlb):
            sline=lines[n].split('=')[0].split(); self.nlbn.append(int(sline[0])); ibtype=0
            self.ilbn.append(array([int(lines[n+1+k].strip().split()[0])-1 for k in arange(self.nlbn[-1])]))
            n=n+1+self.nlbn[-1]

            #add bnd type info
            if len(sline)==2: ibtype=int(sline[1])
            if self.ilbn[-1][0]==self.ilbn[-1][-1]: ibtype=1
            self.island.append(ibtype)
        self.island=array(self.island); self.nlbn=array(self.nlbn); self.ilbn=array(self.ilbn,dtype='O');
        if len(self.ilbn)==1: self.ilbn=self.ilbn.astype('int')

    def read_prop(self,fname):
        '''
        read schism prop, and return the values
        '''
        evi=read_schism_prop(fname)
        if len(evi)!=self.ne: sys.exit("check dimension: ne={}, prop={}".format(self.ne,len(evi)))
        return evi

    def interp_node_to_elem(self,value=None):
        '''
        interpolate node values to element values (works for multi-dimensional data)
            default is self.dp => self.dpe
        '''
        #interpolate
        dp=self.dp if (value is None) else value
        dms=[*dp.shape]; ip=dms.index(self.np); idm=arange(len(dms)); dms[ip]=self.ne
        if len(dms)>1:  idm[0],idm[ip]=ip,0; dms[0],dms[ip]=dms[ip],dms[0] #put dim=np 1st for multi-dimensional data
        fp3=self.i34==3; fp4=~fp3; dp=dp.transpose(idm); dpe=zeros(dms)
        dpe[fp3]=dp[self.elnode[fp3,:3]].mean(axis=1); dpe[fp4]=dp[self.elnode[fp4]].mean(axis=1)
        return dpe.transpose(idm)

    def interp_elem_to_node(self,value=None,fmt=0,p=1):
        '''
        interpolate element values to nodes
        if value not given, dpe is used
        fmt=0: element area weighted avarage in nodal ball
        fmt=1: inverse distance (power=p)
        fmt=2: maximum of surrounding nodal values
        fmt=3: minimum of surrounding nodal values
        fmt=4: simple avarage in nodal ball
        '''
        #element values
        if not hasattr(self,'nne'): self.compute_nne()
        if (value is None) and (not hasattr(self,'dpe')): self.compute_ctr()
        v0=self.dpe if (value is None) else value

        #interpolation
        vs=v0[self.ine]
        if fmt==0:
           if not hasattr(self,'area'): self.compute_area()
           fpn=self.ine!=-1; a=self.area[self.ine]
           ta=sum(a*fpn,axis=1); tv=sum(a*vs*fpn,axis=1)
           fpz=ta!=0; v=zeros(self.np)*nan; v[fpz]=tv[fpz]/ta[fpz]
        if fmt==1:
              dist=abs((self.xctr[self.ine]+1j*self.yctr[self.ine])-(self.x+1j*self.y)[:,None])
              w=1/(dist**p); w[self.ine==-1]=0; tw=w.sum(axis=1); v=(w*vs).sum(axis=1)/tw
        if fmt==2: vs[self.ine==-1]=v0.min()-1; v=vs.max(axis=1)
        if fmt==3: vs[self.ine==-1]=v0.max()+1; v=vs.min(axis=1)
        if fmt==4:
           w=self.ine!=-1; tw=w.sum(axis=1)
           if sum(isnan(v0))!=0:
              vs[~w]=0; v=vs.sum(axis=1)/tw
           else:
              v=(w*vs).sum(axis=1)/tw
        return v

    def compute_all(self,fmt=0):
        '''
        compute all geometry information of hgrid by invoking:
           functions: compute_ctr(),compute_area(),compute_side(fmt=2),compute_nne(),compute_ic3()
           attrs: (xctr,yctr,dpe),(area),(ns,isidenode,isdel,xcj,ycj,dps,distj),(nne,mnei,indel,ine),(ic3,elside)
           fmt=0: skip if already attrs are already available; fmt=1: recompute all attrs
        '''
        if (not hasattr(self,'dpe'))  or fmt==1: self.compute_ctr()
        if (not hasattr(self,'area')) or fmt==1: self.compute_area()
        if (not hasattr(self,'dps'))  or fmt==1: self.compute_side(fmt=2)
        if (not hasattr(self,'ine'))  or fmt==1: self.compute_nne(fmt=1)
        if (not hasattr(self,'ic3'))  or fmt==1: self.compute_ic3()

    def compute_ctr(self):
        '''
        compute element center information: (xctr,yctr,dpe)
        '''
        if not hasattr(self,'xctr'):
           fp3=self.i34==3; fp4=~fp3; self.xctr,self.yctr,self.dpe=zeros([3,self.ne])
           self.xctr[fp3],self.yctr[fp3],self.dpe[fp3]=c_[self.x,self.y,self.dp][self.elnode[fp3,:3]].mean(axis=1).T
           self.xctr[fp4],self.yctr[fp4],self.dpe[fp4]=c_[self.x,self.y,self.dp][self.elnode[fp4]].mean(axis=1).T
        return self.dpe

    def compute_area(self):
        x1,x2,x3,x4=self.x[self.elnode].T; y1,y2,y3,y4=self.y[self.elnode].T
        fp=self.elnode[:,-1]<0; x4[fp]=x1[fp]; y4[fp]=y1[fp]
        self.area=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)+(x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
        return self.area

    def compute_gradient(self,fmt=0,value=None,outfmt=0,cpp=None,lon0=None,lat0=None):
        '''
        Compute gradient on each element first, then convert to node value
          fmt: interploation method (see details in interp_elem_to_node())
              (0: simple avarage; fmt=1: inverse distance; fmt=2: maximum of surrounding nodal values)
          value=None: gd.dp is used;    value: array of [np,] or [ne]
          outfmt=0: return (dpdx,dpdy,dpdxy); outfmt=1: return (dpdx,dpdy,dpdxy,dpedx,dpedy,dpedxy)
          cpp=1: CPP projection. It project decimal degress (lon/lat) to CPP coordinate (https://www.xmswiki.com/wiki/CPP_Coordinate_System).
                 If lon0 and lat0 are not defined, they will be the center of the grid
        '''
        if not hasattr(self,'area'): self.compute_area()
        if not hasattr(self,'dpe'): self.compute_ctr()

        if cpp==1:
            if lon0==None or lat0==None: lon0=(self.x.min()+self.x.max())/2; lat0=(self.y.min()+self.y.max())/2
            x,y=proj_pts(self.x,self.y,'epsg:4326','cpp',lon0=lon0,lat0=lat0)
        else:
            x=self.x; y=self.y

        #get node value
        v0=self.dp if value is None else value
        if len(v0)==self.ne: v0=self.interp_elem_to_node(value=v0)

        #get pts
        x1,x2,x3,x4=x[self.elnode].T; y1,y2,y3,y4=y[self.elnode].T; v1,v2,v3,v4=v0[self.elnode].T
        fp=self.elnode[:,-1]<0; fpn=~fp; x4[fp],y4[fp],v4[fp]=x1[fp],y1[fp],v1[fp]
        a1=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2; a2=((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2

        #compute gradients
        dpedx=(v1*(y2-y3)+v2*(y3-y1)+v3*(y1-y2))/(2*a1)
        dpedy=((x3-x2)*v1+(x1-x3)*v2+(x2-x1)*v3)/(2*a1)

        #for quads
        dpedx2=(v1[fpn]*(y3[fpn]-y4[fpn])+v3[fpn]*(y4[fpn]-y1[fpn])+v4[fpn]*(y1[fpn]-y3[fpn]))/(2*a2[fpn])
        dpedy2=((x4[fpn]-x3[fpn])*v1[fpn]+(x1[fpn]-x4[fpn])*v3[fpn]+(x3[fpn]-x1[fpn])*v4[fpn])/(2*a2[fpn])
        dpedx[fpn]=(dpedx[fpn]+dpedx2)/2;  dpedy[fpn]=(dpedy[fpn]+dpedy2)/2

        #interp to node
        dpedxy=sqrt(dpedx**2+dpedy**2)
        dpdx=self.interp_elem_to_node(value=dpedx,fmt=fmt)
        dpdy=self.interp_elem_to_node(value=dpedy,fmt=fmt)
        dpdxy=self.interp_elem_to_node(value=dpedxy,fmt=fmt)

        if outfmt==0:
           return dpdx,dpdy,dpdxy
        else:
           return dpdx,dpdy,dpdxy,dpedx,dpedy,dpedxy

    def compute_side(self,fmt=0):
        '''
        compute side information of schism's hgrid
        fmt=0: compute ns (# of sides) only
        fmt=1: compute (ns,isidenode,isdel)
        fmt=2: compute (ns,isidenode,isdel), (xcj,ycj,dps,distj), and (nns,ins)
        '''

        #collect sides
        fp3=self.i34==3; self.elnode[fp3,-1]=self.elnode[fp3,0]; sis=[]; sie=[]
        for i in arange(4):
            sis.append(c_[self.elnode[:,(i+1)%4],self.elnode[:,(i+2)%4]]); sie.append(arange(self.ne))
        sie=array(sie).T.ravel(); sis=array(sis).transpose([1,0,2]).reshape([len(sie),2])
        fpn=diff(sis,axis=1)[:,0]!=0; sis=sis[fpn]; sie=sie[fpn]; self.elnode[fp3,-1]=-2

        #sort sides
        usis=sort(sis,axis=1).T; usis,sind,sindr=unique(usis[0]+1j*usis[1],return_index=True,return_inverse=True)
        self.ns=len(sind)

        if fmt==0:
           return self.ns
        elif fmt in [1,2]:
           #build isidenode
           sinda=argsort(sind); sinds=sind[sinda]; self.isidenode=sis[sinds]

           #build isdel
           se1=sie[sinds]; se2=-ones(self.ns).astype('int')
           sindl=setdiff1d(arange(len(sie)),sind); se2[sindr[sindl]]=sie[sindl]; se2=se2[sinda]
           self.isdel=c_[se1,se2]; fps=(se1>se2)*(se2!=-1); self.isdel[fps]=fliplr(self.isdel[fps])

           #compute xcj,ycj and dps
           if fmt==2:
              self.xcj,self.ycj,self.dps=c_[self.x,self.y,self.dp][self.isidenode].mean(axis=1).T
              self.distj=abs(diff(self.x[self.isidenode],axis=1)+1j*diff(self.y[self.isidenode],axis=1))[:,0]

              inode=self.isidenode.ravel(); iside=tile(arange(self.ns),2) #node-side table
              sind=argsort(inode); inode=inode[sind]; iside=iside[sind]
              self.nns=unique(inode,return_counts=True)[1]; self.ins=-ones([self.np,self.nns.max()]).astype('int'); n=0
              for i in arange(self.np): self.ins[i,:self.nns[i]]=iside[n:(n+self.nns[i])]; n=n+self.nns[i]

           return self.ns,self.isidenode,self.isdel

    def compute_bnd(self,bxy=None,nb_max=500000,method=0):
        '''
        compute boundary information. If bxy is provided, define open/land boundries

        bxy: endpoint coordinates of open boundaries. Examples:
            1). bxy=[x1,x2,y1,y2]  #only one open boundary
            2). bxy=[[x1,x2,y1,y2],[x1,x2,y1,y2],...]  #multiple open boundaries
            3). bxy="bpfile", with paired build points sequentially defining each boundaries

        nb_max: maximum boundary nodes allowed for each bnd segments
        method=0: latest method in computing bnd; method=1: old method
        '''
        print('computing grid boundaries')
        if not hasattr(self,'isdel') or not hasattr(self,'isidenode'): self.compute_side(fmt=1)

        if method==0:
           #find boundary side and element, and sort it
           fpn=self.isdel[:,-1]==-1; be=self.isdel[fpn][:,0]; n1,n2=self.isidenode[fpn].T; i34=self.i34[be]; nbs=len(be)
           for i in arange(4): fp=(n1==self.elnode[be,(i+1)%i34])*(n2==self.elnode[be,i%i34]); n1[fp],n2[fp]=n2[fp],n1[fp]
           sinds=dict(zip(n1,[[i] for i in n2]))

           #---------------------------------------------------------------------------
           #this part is used to dealing with illegal mesh
           #---------------------------------------------------------------------------
           #find the nodes connecting to more than 2 sides,and build the side dict
           ibnx,nbnx=unique(n1,return_counts=True); ibnx=ibnx[nbnx>1]
           for i in ibnx:
               #find all nodes in the nodal ball, and compute the angles
               if not hasattr(self,'xctr'): self.compute_ctr()
               ips=unique(self.elnode[self.indel[i]]); ips=setdiff1d(ips[1:] if ips[0]<0 else ips,i); nps=len(ips)
               A=angle((self.x[ips]-self.x[i])+1j*(self.y[ips]-self.y[i])); iA=argsort(A); A,ips=A[iA],ips[iA]
               cA=angle((self.xctr[self.indel[i]]-self.x[i])+1j*(self.yctr[self.indel[i]]-self.y[i])) #angle for element center

               #get all boundary nodes
               ib1=nonzero(n1==i)[0]; ip1=n2[ib1];  ib2=nonzero(n2==i)[0]; ip2=n1[ib2]; sinds[i]=[]
               for n in arange(nps):
                   i1=ips[n]; i2=ips[(n+1)%nps]; A1=A[n]; A2=A[(n+1)%nps]
                   if A2<A1: A2=A2%(2*pi); cA=cA%(2*pi)
                   if sum((cA>A1)*(cA<A2))!=0: continue #element center between i1 and i2
                   if (i1 in ip1) and (i2 in ip2):  #find a pair
                       sinds[i].append([i2,i,i1])
                   elif (i2 in ip1) and (i1 in ip2):
                       sinds[i].append([i1,i,i2])
                   else:
                       sys.exit('wrong in connect side pairs')

           #reconnect sides, but skip the nodes connecting to more than 2 sides
           while True:#append sides to their downstream sides
               iflag=0
               for i in ibnx:
                   for n, side in enumerate(sinds[i]):
                       if side is None: continue
                       if side[-1] in ibnx:
                           iflag=1
                           for m,dside in enumerate(sinds[side[-1]]):
                               if dside is None: continue
                               nsd=dside.index(side[-1])
                               if array_equal(side[-(nsd+1):],dside[:(nsd+1)]):
                                   sinds[i][n]=[*side,*dside[(nsd+1):]]
                                   sinds[side[-1]][m]=None; break
               if iflag==0: break

           while True: #append all sides to their upstream ordinary sides
               iflag=0
               for i in ibnx:
                   for n, side in enumerate(sinds[i]):
                       if side is None: continue
                       sinds[i][n]=None; iflag=1
                       if side[0] in ibnx: #append to upstream side
                           for m,uside in enumerate(sinds[side[0]]):
                               if uside is None: continue
                               isd=uside.index(side[0]); nsd=len(uside)-isd
                               if array_equal(uside[isd:],side[:nsd]):
                                   sinds[side[0]][m]=[*uside,*side[nsd:]]
                                   sinds[i][n]=None; break
                       else: #append to other side
                           sinds[side[0]]=side[1:]
               if iflag==0: break
           for i in ibnx: del sinds[i]
           #---------------------------------------------------------------------------

           #search boundary from other nodes
           ids=list(sinds.keys()); ns=len(ids); sindf=dict(zip(ids,arange(ns))); ifb=ones(ns).astype('int'); nb=0; nbn=[]; ibn=[]
           while(sum(ifb)!=0):
               #start points
               ib=nonzero(ifb==1)[0][0]; id0=ids[ib]; id=sinds[id0][-1]; ifb[sindf[id0]]=0; ibni=[id0,*sinds[id0]]; iloop=0

               #search each segments
               while True:
                   idp=id; id=sinds[id][-1]; ifb[sindf[idp]]=0
                   if(id==id0): ibni.extend(sinds[idp][:-1]); break
                   ibni.extend(sinds[idp]); iloop=iloop+1
                   if iloop>nb_max: print('grid boundary search failed: check your grid'); return array(ibni)
               nb=nb+1; nbn.append(len(ibni)); ibn.append(array(ibni))
        else:
            #find boundary side and element
            fpn=self.isdel[:,-1]==-1;  isn=self.isidenode[fpn]; be=self.isdel[fpn][:,0]; nbs=len(be)

            #sort isn
            i2=ones(nbs).astype('int'); fp3=nonzero(self.i34[be]==3)[0]; fp4=nonzero(self.i34[be]==4)[0]
            for i in arange(4):
                if i==3:
                    i1=self.elnode[be[fp4],3]; i2=self.elnode[be[fp4],0]
                    fp=(isn[fp4,0]==i2)*(isn[fp4,1]==i1); isn[fp4[fp]]=fliplr(isn[fp4[fp]])
                else:
                    i1=self.elnode[be,i]; i2[fp3]=self.elnode[be[fp3],(i+1)%3]; i2[fp4]=self.elnode[be[fp4],i+1]
                    fp=(isn[:,0]==i2)*(isn[:,1]==i1); isn[fp]=fliplr(isn[fp])

            #compute all boundaries
            ifb=ones(nbs).astype('int'); nb=0; nbn=[]; ibn=[]
            sinds=dict(zip(isn[:,0],arange(nbs))) #dict for sides

            #find the nodes connecting to more than 2 sides
            ibnx,nbnx=unique(isn[:,0],return_counts=True); idx=[]; nbx=0
            for i in ibnx[nbnx>1]: k=nonzero(isn[:,0]==i)[0]; idx.extend(k); nbx=nbx+len(k)

            #search boundary from other nodes
            while(sum(ifb)!=0):
                #start points
                if nbx>0:
                    nbx=nbx-1; id0,id=isn[idx[nbx]]; ifb[idx[nbx]]=0
                else:
                    id0=isn[nonzero(ifb==1)[0][0],0]; id=isn[sinds[id0],1]; ifb[sinds[id0]]=0
                ibni=[id0,id]; ifb[sinds[id]]=0; iloop=0

                #search each segments
                while True:
                    id=isn[sinds[id],1]; ifb[sinds[id]]=0
                    if(id==id0): break
                    ibni.append(id); iloop=iloop+1
                    if iloop>nb_max: print('grid boundary search failed: check your grid'); return array(ibni)
                nb=nb+1; nbn.append(len(ibni)); ibn.append(array(ibni))

        #sort bnd
        nbn=array(nbn); ibn=array(ibn,dtype='O'); fps=flipud(argsort(nbn)); nbn,ibn=nbn[fps],ibn[fps]
        if ibn.shape[0]==1: ibn=ibn.astype('int')

        #find the outline
        island=ones(nb).astype('int'); bid=[]
        for i in arange(nb):
            px=self.x[ibn[i].astype('int')]; i0=nonzero(px==px.min())[0][0]
            sid=ibn[i][array([(i0-1)%nbn[i],i0,(i0+1)%nbn[i]])].astype('int')
            if signa(self.x[sid],self.y[sid])>0: island[i]=0; bid.append(i)

        #put outline bnd ahead
        if len(bid)!=0:
           bid=array(bid); nbid=setdiff1d(arange(nb),bid); island=array([*island[bid],*island[nbid]])
           nbn=array([*nbn[bid],*nbn[nbid]]); ibn=array([*ibn[bid],*ibn[nbid]],dtype='O')
           if ibn.shape[0]==1: ibn=ibn.astype('int')

        #save boundary information
        if not hasattr(self,'bndinfo'): self.bndinfo=zdata()
        ip=[]; sind=[]; S=self.bndinfo
        for m,ibni in enumerate(ibn): ip.extend(ibni); sind.extend(tile(m,len(ibni)))
        ip=array(ip); sind=array(sind); S.sind=sind; S.ip=ip; S.island=island
        S.nb=nb; S.nbn=nbn; S.ibn=ibn; S.x=self.x[ip]; S.y=self.y[ip]

        #add to grid bnd info
        if (not hasattr(self,'nob')) and (bxy is None):
           self.nob=0; self.nobn=array([]); self.iobn=array([[]])
           self.nlb=S.nb+1; self.nlbn=array([2,*S.nbn]); self.island=array([0,*S.island])
           self.ilbn=array([array([S.ibn[0][-1],S.ibn[0][0]]).astype('int'), *S.ibn],dtype='O')

        #----------------------------------------------------------------------------
        #define open/land/island boundaries
        #----------------------------------------------------------------------------
        if bxy is not None:
           if isinstance(bxy,str): #bxy is a bpfile
              bp=read_schism_bpfile(bxy); bxy=[]
              if bp.nsta%2!=0: sys.exit('even number of points are required')
              for i in arange(int(bp.nsta/2)): bxy.append([bp.x[2*i],bp.x[2*i+1],bp.y[2*i],bp.y[2*i+1]])

           S=self.bndinfo; bf=array([zeros(i) for i in S.nbn],dtype='O')
           pxy=array([*bxy]) if array([*bxy]).ndim==2 else array([*bxy])[None,:]
           gxy=c_[self.x[S.ip],self.y[S.ip]]; p1=near_pts(pxy[:,0::2],gxy); p2=near_pts(pxy[:,1::2],gxy)
           sindb,sindb2=S.sind[p1],S.sind[p2]; p1,p2=S.ip[p1],S.ip[p2];
           if not array_equal(sindb,sindb2): sys.exit('piont pairs are not on the same boundary')

           #parse each segment for open boundaries
           self.nob=len(pxy); self.nobn=[]; self.iobn=[]; self.nlb=0; self.nlbn=[]; self.ilbn=[]; self.island=[]
           for m in sindb[sort(unique(sindb, return_index=True)[1])]:
               sind=nonzero(sindb==m)[0];  sids=[]; bf=ones(S.nbn[m])
               if len(sind)==0: continue
               for n,sindi in enumerate(sind): # add open bnd
                   id1=nonzero(S.ibn[m]==p1[sindi])[0][0]; id2=nonzero(S.ibn[m]==p2[sindi])[0][0]; sids.extend([id1,id2])
                   itype=(((id1<id2) and (S.island[m]==0)) or ((id1>id2) and (S.island[m]==1)))
                   sid=arange(id1,id2+1) if (id1<id2) else r_[arange(id1,S.nbn[m]),arange(id2+1)]
                   self.nobn.append(sid.size); self.iobn.append(S.ibn[m][sid]); bf[sid]=0
               for n,id1 in enumerate(sids):  #add land boundaries
                   id2=sids[(n+1)%len(sids)]
                   itype=(((id1<id2) and (S.island[m]==0)) or ((id1>id2) and (S.island[m]==1)))
                   sid=arange(id1,id2+1) if (id1<id2) else r_[arange(id1,S.nbn[m]),arange(id2+1)]
                   if sum(bf[sid])!=0:
                      self.nlb=self.nlb+1; self.nlbn.append(sid.size); self.ilbn.append(S.ibn[m][sid]); self.island.append(0)
           for m in arange(S.nb): #add remaining land bnd
               sind=nonzero(sindb==m)[0]
               if len(sind)!=0: continue
               self.nlb=self.nlb+1; self.nlbn.append(S.nbn[m]); self.ilbn.append(S.ibn[m]); self.island.append(S.island[m])
           self.nobn=array(self.nobn); self.iobn=array(self.iobn,dtype='O')
           self.nlbn=array(self.nlbn); self.ilbn=array(self.ilbn,dtype='O'); self.island=array(self.island)

    def compute_node_ball(self,**args):
        '''
        alias for compute_nne()
        '''
        return self.compute_nne(**args)

    def compute_nne(self,fmt=0):
        '''
        compute nodal ball information: nne,mnei,indel,ine,inp
        where:
             nne:   number of elements in nodal ball
             mnei:  maximum number of elements in nodal ball
             indel: indices for each nodal ball
             ine:   indices for each nodal ball, but in maxtrix " shape=[np,max(nne)]"
             nnp:   number of nodes in nodal ball (exclude itself)
             mnpi:  maximum number of nodes in nodal ball
             inp:   node indices for each nodal ball
             indnd: node indices for each nodal ball

        fmt=0: not compute inp; fmt=1: compute inp
        '''

        #get index of all node and elements
        elem=tile(arange(self.ne),[4,1]).T.ravel(); node=self.elnode.ravel()
        fpn=node!=-2;      elem,node=elem[fpn],node[fpn]
        fpn=argsort(node); elem,node=elem[fpn],node[fpn]

        #compute nne,ine,indel
        unode,sind,self.nne=unique(node,return_index=True,return_counts=True); self.mnei=self.nne.max()
        self.ine=-ones([self.np,self.mnei]).astype('int')
        for i in arange(self.mnei): fpe=self.nne>i; sinde=sind[fpe]+i; self.ine[fpe,i]=elem[sinde]
        self.indel=array([array(i[:k]) for i,k in zip(self.ine,self.nne)],dtype='O')

        #compute inp
        if fmt==1:
           #self.inp=array([setdiff1d(self.elnode[self.indel[i]].ravel(),[i,-2]) for i in arange(self.np)],dtype='O')
           inp=self.elnode[self.ine].reshape([self.np,4*self.mnei])
           inp=[set(i[:4*k]) for m,[k,i] in enumerate(zip(self.nne,inp))]
           [i.remove(-2) for i in inp if (-2 in i)]; [i.remove(k) for k,i in enumerate(inp)]
           self.inp=array([array([*i]) for i in inp],dtype='O')
           self.nnp=array([len(i) for i in self.inp]); self.mnpi=self.nnp.max(); self.indnd=-ones([self.np,self.mnpi]).astype('int')
           for i,[k,n] in enumerate(zip(self.nnp,self.inp)): self.indnd[i,:k]=n

        return self.nne

    def compute_ic3(self):
        '''
        compute element-to-side table, where
             elside: element-to-side table
        '''

        #get index for all elements and sides
        if not hasattr(self,'isdel'): self.compute_side(fmt=1)
        side=tile(arange(self.ns),[2,1]).T.ravel(); elem=self.isdel.ravel()
        fpn=elem!=-1; side,elem=side[fpn],elem[fpn]
        fpn=argsort(elem); side,elem=side[fpn],elem[fpn]

        #build elside
        uelem,sind=unique(elem,return_index=True); self.elside=-ones([self.ne,4]).astype('int'); m34=self.i34.max()
        for i in arange(m34):
            fps=nonzero(self.i34>i)[0]; i34=self.i34[fps]; sinds=sind[fps]+i
            sd=side[sinds];  n1,n2=self.isidenode[sd].T
            for k in arange(m34): #sort order of sides
                id1,id2=(k+1)%i34,(k+2)%i34
                fpk=((self.elnode[fps,id1]==n1)*(self.elnode[fps,id2]==n2))|((self.elnode[fps,id1]==n2)*(self.elnode[fps,id2]==n1))
                self.elside[fps[fpk],k]=sd[fpk]
        self.elside[self.i34==3,-1]=-1

        #build ic3
        self.ic3=-ones([self.ne,4]).astype('int'); ie=arange(self.ne)
        for i in arange(m34):
            es=self.isdel[self.elside[:,i]]; fp=es[:,0]==ie
            self.ic3[fp,i]=es[fp,1]; self.ic3[~fp,i]=es[~fp,0]
        self.ic3[self.elside==-1]=-1
        return self.ic3,self.elside

    def compute_acor(self,pxy,fmt=0):
        '''
        compute acor coodinate for points pxy[npt,2]

        usage: ie,ip,acor=compute_acor(c_[xi,yi]), where xi and yi are array of coordinates
        outputs: ie[npt],ip[npt,3],acor[npt,3]
               ie:  the element number
               ip:  the nodal indices of the ie
               acor: the area coordinate
               fmt=0: (default) faster method by searching the neighbors of elements and nodes
               fmt=1: slower method using point-wise comparison

               Note: for interpolation of few pts on a large grid, fmt=1 can be faster than fmt=0
        '''

        npt=len(pxy); pip=-ones([npt,3]).astype('int'); pacor=zeros([npt,3])
        if fmt==0:
           pie=-ones(npt).astype('int'); sindp=arange(npt)
           #search element ctr
           if not hasattr(self,'xctr'): self.compute_ctr()
           #if hasattr(self,'bndinfo'): sindp=sindp[self.inside_grid(pxy)==1]
           sinde=near_pts(pxy[sindp],c_[self.xctr,self.yctr]); fps,sip,sacor=self.inside_elem(pxy[sindp],sinde)
           if len(fps)!=0: pie[sindp[fps]]=sinde[fps]; pip[sindp[fps]]=sip; pacor[sindp[fps]]=sacor

           #search the direct neighbors of element ctr
           fp=pie[sindp]==-1; sindp,sinde=sindp[fp],sinde[fp]
           if len(sindp)!=0:
               if not hasattr(self,'ic3'): self.compute_ic3()
               for i in arange(self.i34.max()):
                   ie=self.ic3[sinde,i]; fps,sip,sacor=self.inside_elem(pxy[sindp],ie)
                   if len(fps)!=0: pie[sindp[fps]]=ie[fps]; pip[sindp[fps]]=sip; pacor[sindp[fps]]=sacor

                   #update sindp
                   fp=pie[sindp]==-1; sindp,sinde=sindp[fp],sinde[fp]
                   if len(sindp)==0: break

           #search elements inside node ball
           if len(sindp)!=0:
               if not hasattr(self,'ine'): self.compute_nne()
               sindn=near_pts(pxy[sindp],c_[self.x,self.y]); pip[sindp]=sindn[:,None]; pacor[sindp,0]=1
               for i in arange(self.mnei):
                   ie=self.ine[sindn,i]; fps,sip,sacor=self.inside_elem(pxy[sindp],ie)
                   if len(fps)!=0: pie[sindp[fps]]=ie[fps]; pip[sindp[fps]]=sip; pacor[sindp[fps]]=sacor

                   #update sindp
                   if i<(self.mnei-1): fp=(pie[sindp]==-1)*(self.ine[sindn,i+1]!=-1); sindp,sindn=sindp[fp],sindn[fp]
                   if len(sindp)==0: break

           #use point-wise method for the remain pts
           sindp=nonzero(pie==-1)[0]; sindp=sindp[self.inside_grid(pxy[sindp])==1]
           if len(sindp)!=0: pie[sindp],pip[sindp],pacor[sindp]=self.compute_acor(pxy[sindp],fmt=1)

        elif fmt==1:
            #check 1st triangle
            sindn=self.elnode.T[:3]; pie=inside_polygon(pxy,self.x[sindn],self.y[sindn],fmt=1)
            fps=pie!=-1; pip[fps]=sindn.T[pie[fps]]

            #check 2nd triangle
            sind4=nonzero(self.i34==4)[0]; sind2=nonzero(~fps)[0]
            if len(sind2)!=0 and len(sind4)!=0:
               sindn=self.elnode[sind4].T[array([0,2,3])]; pie2=inside_polygon(pxy[sind2],self.x[sindn],self.y[sindn],fmt=1)
               fps=pie2!=-1; pie[sind2[fps]]=sind4[pie2[fps]]; pip[sind2[fps]]=sindn.T[pie2[fps]]

            #compute acor
            fpn=pie!=-1
            if sum(fpn)!=0:
               x1,x2,x3=self.x[pip[fpn]].T; y1,y2,y3=self.y[pip[fpn]].T; x,y=pxy[fpn].T
               A1=signa(c_[x,x2,x3],c_[y,y2,y3]); A2=signa(c_[x1,x,x3],c_[y1,y,y3])
               A=signa(c_[x1,x2,x3],c_[y1,y2,y3]); pacor[fpn]=c_[A1/A,A2/A,1-(A1+A2)/A]
            if sum(~fpn)!=0:
               sindn=near_pts(pxy[~fpn],c_[self.x,self.y]); pip[~fpn]=sindn[:,None]; pacor[~fpn,0]=1

        return pie,pip,pacor

    def compute_kb(self,kbp,fmt=0):
        '''
        fmt=0: compute bottom element indices
        fmt=1: compute bottom side indices
        '''
        if fmt==0: kb=kbp[self.elnode]; kb[self.i34==3,-1]=-1; kb=kb.max(axis=1)
        if fmt==1:
           if not hasattr(self,'isidenode'): self.compute_side(fmt=1)
           kb=kbp[self.isidenode].max(axis=1)
        return kb

    def compute_angle(self):
        '''
        compute internal angles
        '''
        ie=arange(self.ne); angles=[]
        for i in arange(4):
            i1=(i-1)%self.i34; i2=i%self.i34; i3=(i+1)%self.i34
            id=self.elnode[ie[:,None],c_[i1,i2,i3]]; x1,x2,x3=self.x[id].T; y1,y2,y3=self.y[id].T
            a=(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))%(2*pi); angles.append(a*180/pi)
        self.angles=array(angles).T
        return self.angles

    def interp(self,pxy,value=None,fmt=0):
        '''
        interpolate to get value at pxy
          pxy: c_[x,y]
          value=None: gd.dp is used; value: array of [np,] or [ne,]
          fmt=0: (default) faster method by searching the neighbors of elements and nodes
          fmt=1: slower method using point-wise comparison

          Note: for interpolation of few pts on a large grid, fmt=1 can be faster than fmt=0
        '''
        #get node value
        vi=self.dp if value is None else value
        if len(vi)==self.ne: vi=self.interp_elem_to_node(value=vi)

        #interp
        pip,pacor=self.compute_acor(pxy,fmt=fmt)[1:]
        return (vi[pip]*pacor).sum(axis=1)

    def scatter_to_grid(self,fmt=0,reg_in=1,reg_out=1,**args):
        '''
        construct a new schism grid based on:
           fmt=0:  node
           fmt=1:  element center
           fmt=2:  side center

        reg_in=1:  delete elements inside islands;      reg_in=0: not delete
        reg_out=1: delete elements outside grid domain; reg_out=0: not delete

        '''
        #get regions
        if reg_in==1 or reg_out==1:
           if not hasattr(self,'bndinfo'): self.compute_bnd()
           if reg_in==1:  reg_in=[c_[self.x[i],self.y[i]] for i in self.bndinfo.ibn[1:]]
           if reg_out==1: reg_out=[c_[self.x[i],self.y[i]] for i in self.bndinfo.ibn[:1]]
        if reg_in==0:  reg_in=None
        if reg_out==0: reg_out=None

        #get scatter
        if fmt==0: xyz=c_[self.x,self.y,self.z]
        if fmt==1:
           if not hasattr(self,'dpe'): self.compute_ctr()
           xyz=c_[self.xctr,self.yctr,self.dpe]
        if fmt==2:
           if not hasattr(self,'dps'): self.compute_side(fmt=2)
           xyz=c_[self.xcj,self.ycj,self.dps]

        #build grid
        gdn=scatter_to_schism_grid(xyz,reg_out=reg_out,reg_in=reg_in,**args)
        return gdn

    def smooth(self,dist,value=None,fmt=0,ms=1e8):
        '''
        smooth field by averaging values within radius of dist
          dist: averaging radius
          value=None: gd.dp is used; value: array of [np,] or [ne,]
          fmt=0: return smoothed nodal values; fmt=1: return smoothed elem. values
          ms: matrix size; reduce ms to avoid "out of memory" error
        '''
        from scipy.spatial.distance import cdist

        if not hasattr(self,'dpe'): self.compute_ctr()
        if not hasattr(self,'area'): self.compute_area()
        if value is None: value=self.dpe
        if value.size==self.np: value=self.interp_node_to_elem(value)

        #smooth unstructure data
        ne,x,y=self.ne,self.xctr,self.yctr; exy=x+1j*y; pxy=x+1j*y
        w=self.area; v=value*w; sinde=arange(ne); pv=zeros(ne); sindp=arange(ne)
        while sindp.size!=0:
            #print(sindp.size,self.ne,sindp.size/self.ne)
            fpc=abs(pxy-pxy[0])<=dist; pid=sindp[fpc] #node inside circle of dc
            eid=sinde[nonzero(abs(exy-pxy[0])<=(2*dist))[0]]; #elem. inside circle of 2*dc
            if eid.size==1: pv[pid]=value[eid]

            #find dist maxtrix, get weight and value, and assign average value at node
            ds=pid.size*eid.size; nloop=int(ceil(ds/ms)); dsb=int(ceil(pid.size/nloop))
            for n in arange(nloop):
                pidi=pid[arange((dsb*n),min(dsb*(n+1),pid.size)).astype('int')]
                pdist=cdist(c_[x[pidi],y[pidi]],c_[x[eid],y[eid]]); fpd=pdist>dist
                ew=tile(w[eid],[pidi.size,1]); ev=tile(v[eid],[pidi.size,1]); ew[fpd]=0; ev[fpd]=0
                pv[pidi]=ev.sum(axis=1)/ew.sum(axis=1)
            sindp,pxy=sindp[~fpc],pxy[~fpc] #update remaining nodes
        if fmt==1:
           return pv
        else:
           return self.interp_elem_to_node(pv)

    def compute_zcor(self,vgrid,eta=None):
        '''
        compute zcor @ each nodes; need inputs:
           vgrid: read_schism_vgrid('vgrid.in')
           eta (optional): surface elevation (dim=[np] or [1])
        '''
        if eta is None:
           zcor=self.zcor if hasattr(self,'zcor') else vgrid.compute_zcor(self.dp)
           if not hasattr(self,'zcor'): self.zcor=zcor #save zcor
        else:
           zcor=vgrid.compute_zcor(self.dp,eta=eta)
        return zcor

    def compute_volume(self,vgrid,fmt=0,value=None,eta=None):
        '''
        compute volume or volume*value for schism grid
        Inputs:
           vgrid: read_schism_vgrid('vgrid.in')
           eta   (optional): surface elevation (dim=[np] or [1])
           value (optional): tracer concentration with dimension size (ne, nvrt) or (np, nvrt).
        Outputs:
          fmt=0: each prism   (ne,nvrt-1)
          fmt=1: each column  (ne,)
          fmt=2: total volumn (1,)
        '''
        #compute prism volume
        if not hasattr(self,'area'): self.compute_area()
        zcor=self.compute_zcor(vgrid,eta=eta); ze=zeros([self.ne,vgrid.nvrt]); fp3=self.i34==3; fp4=~fp3
        ze[fp3]=zcor[self.elnode[fp3,:3]].mean(axis=1); ze[fp4]=zcor[self.elnode[fp4]].mean(axis=1)
        v=(ze[:,1:]-ze[:,:-1])*self.area[:,None]

        #multiplee value
        if value is not None:
           if value.shape[0]==self.np: value=self.interp_node_to_elem(value); value[:,1:]=(value[:,:-1]+value[:,1:])/2
           v=v*value[:,1:]

        return v if fmt==0 else v.sum(axis=1) if fmt==1 else v.sum()

    def compute_contour(self,levels,value=None):
        '''
        compute contour lines
        '''
        value=self.dp if value is None else value
        if value.size==self.ne: value=self.interp_elem_to_node(value)

        #plot contour and extract the lines
        fp4=self.i34==4; trs=r_[self.elnode[:,:3],c_[self.elnode[fp4,0],self.elnode[fp4,2:]]]
        hf=figure(); hf.set_visible(False)
        P=tricontour(self.x,self.y,trs,value,levels=levels); close(hf); cxy=[]
        for k in arange(len(P.collections)):
            p=P.collections[k].get_paths()
            for i in arange(len(p)):
                xii,yii=p[i].vertices.T
                xi=r_[xii,NaN] if i==0 else r_[xi,xii,NaN]
                yi=r_[yii,NaN] if i==0 else r_[yi,yii,NaN]
            cxy.append(c_[xi,yi])
        return array(cxy) if len(cxy)==1 else array(cxy,dtype='O')

    def write(self, fname=None,**args):
        '''
        generic function to save grid as different format
        format: (*.gr3, *.ll, *.ic, *.prop, *.npz, *.pkl, *.shp, *.bnd), (exceptions are in *.gr3 format)
        examples:
           1). gd.write('grid.npz')  or gd.write('grid')
           2). gd.write('grid.pkl')
           3). gd.write('hgrid.gr3') or gd.write('hgrid.ll') or gd.write('temp.ic',value=5)
        '''
        #check fname, and output as *.npz or *.pkl format
        if fname is None: fname = '{}.npz'.format(os.path.splitext(self.source_file)[0]) #check fname
        if fname.endswith('.npz') or fname.endswith('.pkl'): savez(fname,self,**args); return

        #outputs grid as other format
        F=None
        if fname.endswith('.prop'): F=self.write_prop
        if fname.endswith('.bnd'):  F=self.write_bnd
        if fname.endswith('.shp'):  F=self.write_shp
        if fname.endswith('.gr3') or fname.endswith('.ll') or fname.endswith('.ic') or (F is None): F=self.write_hgrid
        F(fname,**args)

    def save(self, fname=None,**args):
        '''
        alias to self.write
        '''
        self.write(fname,**args)

    def write_hgrid(self,fname,value=None,fmt=0,outfmt='{:<.8f}',elnode=1,bndfile=None,Info=None):
        '''
        write *.gr3 file
            fname: file name
            value: depth value to be outputted
                   value=const: uniform value in space
                   value=dp[np]: specify depth value
                   value=None:  grid's default depth self.dp is used
            fmt=0: not output grid boundary info.; fmt=1: output grid boundary info.
            outfmt: depth format
            elnode=1: output grid connectivity; elnode=0: not output grid connectivity
            bndfile=filepath:  if bndfile is not None, append it at the end of file
            Info: annotation of the gr3 file
        '''

        #get depth value
        if value is None:
           dp=self.dp
        else:
           if hasattr(value,"__len__"):
              dp=value
           else:
              dp=ones(self.np)*value
        if fmt==1: elnode=1

        #write *gr3
        with open(fname,'w+') as fid:
            fid.write('!grd info:{}\n'.format(Info))
            fid.write('{} {}\n'.format(self.ne,self.np))
            lineformat='{:<d} {:<.8f} {:<.8f} '+outfmt+'\n'
            for i in arange(self.np):
                fid.write(lineformat.format(i+1,self.x[i],self.y[i],dp[i]))
            if elnode!=0:
                for i in arange(self.ne):
                    if self.i34[i]==3: fid.write('{:<d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))
                    if self.i34[i]==4: fid.write('{:<d} {:d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))

            #write bnd information
            if fmt==1 and bndfile is None:
               if not hasattr(self,'nob'): self.compute_bnd()
               self.write_bnd(fid=fid)
            if bndfile is not None: fid.writelines(open(bndfile,'r').readlines())

    def write_bnd(self,fname='grd.bnd',fid=None):
        '''
        write grid's boundary information
            fname: name of boundary information
            fid: file handle
        '''
        if not hasattr(self,'nob'): self.compute_bnd()
        bid=open(fname,'w+') if fid is None else fid

        #open bnd
        bid.write('{} = Number of open boundaries\n'.format(self.nob))
        bid.write('{} = Total number of open boundary nodes\n'.format(int(sum(self.nobn))))
        for i in arange(self.nob):
            bid.write('{} = Number of nodes for open boundary {}\n'.format(self.nobn[i],i+1))
            bid.writelines(['{}\n'.format(k+1) for k in self.iobn[i]])

        #land bnd
        bid.write('{} = number of land boundaries\n'.format(self.nlb))
        bid.write('{} = Total number of land boundary nodes\n'.format(int(sum(self.nlbn)))); nln=int(sum(self.island==0))
        for i in arange(self.nlb):
            if self.island[i]==0:
               bid.write('{} {} = Number of nodes for land boundary {}\n'.format(self.nlbn[i],self.island[i],i+1))
            else:
               bid.write('{} {} = Number of nodes for island boundary {}\n'.format(self.nlbn[i],self.island[i],i+1-nln))
            bid.writelines(['{}\n'.format(k+1) for k in self.ilbn[i]])
        if fid is not None: bid.close()

    def write_prop(self,fname='schism.prop',value=None,fmt='{:8.5f}'):
        '''
        write schism prop file.
            fname: file name
            value: prop value;
                   1). if value=None, self.dpe is outputed.
                   2). value=const
                   3). value=array[gd.ne]
            fmt:   output format of prop value
        '''

        #get prop value
        if value is None:
           if not hasattr(self,'dpe'): self.compute_ctr()
           pvi=self.dpe.copy()
        else:
           if hasattr(value,"__len__"):
              pvi=value
           else:
              pvi=ones(self.ne)*value
        if 'd' in fmt: pvi=pvi.astype('int')

        #prepare values
        fstr=('{:d} '+fmt+' \n')*int(self.ne)
        fval=array([range(1,self.ne+1),pvi],dtype='O').T

        #write prop value
        fid=open(fname,'w+'); fid.writelines(fstr.format(*fval.ravel())); fid.close()

    def grd2sms(self,fname='hgrid.2dm'):
        '''
          convert grid to *.2dm format and save
        '''

        lines=[]; lines.append('MESH2D\n')
        for i in arange(self.ne):
            if self.i34[i]==3: lines.append('E3T {} {} {} {} 1\n'.format(i+1,*(self.elnode[i,:3]+1)))
            if self.i34[i]==4: lines.append('E4Q {} {} {} {} {} 1\n'.format(i+1,*(self.elnode[i]+1)))
        for i in arange(self.np):
            lines.append('ND {} {:.8f} {:.8f} {:.8f}\n'.format(i+1,self.x[i],self.y[i],self.dp[i]))
        fid=open(fname,'w+'); fid.writelines(lines); fid.close()

    def split_quads(self,angle_ratio=None,angle_min=None,angle_max=None,fname=None):
        '''
        1). check the quality of quads, violations can be:
            a. angle_min/angle_max <= angle_ratio
            b. internal angle <= angle_min,
            c. internal angle >= angle_max,
            if no parameters are provided, angle_ratio=0.5 will be use
        2). fname: output a new grid "fname" if fname!=None
        '''
        if not hasattr(self,'index_bad_quad'): self.check_quads(angle_ratio,angle_min,angle_max)

        #compute (angle_max-angle_min) in splitted triangle
        qind=self.index_bad_quad; x=self.x[self.elnode[qind,:]]; y=self.y[self.elnode[qind,:]]

        #compute difference between internal angles
        for i in arange(4):
            nid=r_[(i-1)%4,i,(i+1)%4]; x1,x2,x3=x[:,nid].T; y1,y2,y3=y[:,nid].T
            a1=(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))%(2*pi)
            a2=(angle((x2-x3)+1j*(y2-y3))-angle((x1-x3)+1j*(y1-y3)))%(2*pi)
            a3=(angle((x3-x1)+1j*(y3-y1))-angle((x2-x1)+1j*(y2-y1)))%(2*pi)
            a=c_[a1,a2,a3]*180/pi; da=a.max(axis=1)-a.min(axis=1) #compute amax-amin
            A=da if i==0 else c_[A,da]
        flag=sign(A[:,0]+A[:,2]-A[:,1]-A[:,3]) #flag for spliting quads

        #split quads
        nea=len(self.index_bad_quad); elnode=tile(-2,[nea,4]); nds=self.elnode[qind]; self.i34[qind]=3
        fp=flag>=0; ie=qind[fp]; self.elnode[ie,:3]=nds[fp,:3]; self.elnode[ie,-1]=-2; elnode[fp,:3]=nds[fp][:,array([2,3,0])]
        fp=flag<0;  ie=qind[fp]; self.elnode[ie,:3]=nds[fp][:,array([1,2,3])]; self.elnode[ie,-1]=-2; elnode[fp,:3]=nds[fp][:,array([3,0,1])]
        self.ne=self.ne+nea; self.i34=r_[self.i34,tile(3,nea)]; self.elnode=r_[self.elnode,elnode]

        #write new grids
        if fname is not None: self.write_hgrid(fname)

    def check_quads(self,angle_ratio=None,angle_min=None,angle_max=None,fname='bad_quad.bp'):
        '''
        1). check the quality of quads, violations can be:
            a. angle_min/angle_max <= angle_ratio
            b. internal angle <= angle_min,
            c. internal angle >= angle_max,
            if no parameters are provided, angle_ratio=0.5 will be used
        2). the locations of bad quads are saved in file "fname"
        '''

        if angle_ratio is None and angle_min is None and angle_max is None: angle_ratio=0.5
        if not hasattr(self,'angles'): self.compute_angle()
        if not hasattr(self,'xctr'): self.compute_ctr()

        #check violation
        qind=nonzero(self.i34==4)[0]; A=self.angles[qind]; fp=qind==-999
        if angle_ratio is not None: fp=fp|((A.min(axis=1)/A.max(axis=1))<angle_ratio)
        if angle_min is not None: fp=fp|(A.min(axis=1)<=angle_min)
        if angle_max is not None: fp=fp|(A.max(axis=1)>=angle_max)
        self.index_bad_quad=qind[nonzero(fp)[0]]

        #output bad_quad location as bp file
        sbp=schism_bpfile(); sbp.x,sbp.y=c_[self.xctr,self.yctr][self.index_bad_quad].T; sbp.save(fname)
        return self.index_bad_quad

    def plot_bad_quads(self,color='r',ms=12,*args):
        #plot grid with bad quads
        if not hasattr(self,'index_bad_quad'): self.check_quads()
        if not hasattr(self,'xctr'): self.compute_ctr()

        qxi=self.xctr[self.index_bad_quad]; qyi=self.yctr[self.index_bad_quad]
        self.plot()
        plot(qxi,qyi,'.',color=color,ms=ms,*args)
        #show(block=False)
        pass

    def split_quads_wwm(self,fname='hgrid_WWM.gr3'):
        '''
        split quads for WWM model, and output "hgrid_WWM.gr3"
        '''
        sid=nonzero(self.i34==4)[0]; ne,nea=self.ne,len(sid); self.elnode=resize(self.elnode,[ne+nea,4])
        self.elnode[ne:,:3]=self.elnode[sid][:,array([0,2,3])] #add new elements
        self.ne=ne+nea; self.elnode[sid,-1]=-2; self.i34=tile(3,ne+nea)
        self.save(fname)

    def proj(self,prj0,prj1='epsg:4326',fmt=0,x=None,y=None,lon0=None,lat0=None):
        '''
        transform the projection of schism grid's coordinates
        Inputs:
            prj0: projection name of schism grid
            prj1: target projection name; default is 'epsg:4326'
            fmt=0: change gd.x,gd.y to transform xy; fmt=1: only return transformed xy
            x,y: values of grid coordiantes; default is (gd.x, gd.y)
            lon0,lat0: lon&lat of cpp projection center; needed only if 'cpp' in [prj0,prj1]
                       if ("ll"=>"cpp") and (lon0 or lat0 is not provided): lon0=mean(x); lat0=mean(y)
        '''
        if (x is None) or (y is None): x=self.x; y=self.y
        x1,y2=proj(prj0=prj0,prj1=prj1,x=x,y=y,lon0=lon0,lat0=lat0)
        if fmt==0: self.x,self.y=x1,y2
        return [x1,y2]

    def check_skew_elems(self,threshold=None,angle_min=None,fname=None):
        '''
        check schism grid's skewness based on criteria:
            a. xmgredit length_ratio > threshold
            b. internal angle , angle_min
            Note: if no parameters, threshold=30 is used

        Inputs:
            threshold: xmgredit length_ratio of longest_side_length/quivalent_element_radius
            angle_min: minimum angle
            fname!=None: output skew_element bpfile
        Outputs:
            element indices of skew elements

        Examples:
            1). gd.check_skew_elems(), or gd.check_skew_elems(threshold=30)
            2). gd.check_skew_elems(15,10), or gd.check_skew_elems(threshold=15,angle_min=10)
        '''
        if threshold is None and angle_min is None: threshold=30
        if not hasattr(self,'dpe'): self.compute_ctr()

        if threshold is not None: #check length ratio
           if not hasattr(self,'distj'):  self.compute_side(fmt=2)
           if not hasattr(self,'elside'): self.compute_ic3()
           if not hasattr(self,'area'):   self.compute_area()
           distj=self.distj[self.elside]; distj[self.elside==-1]=0  #side length
           srat=distj.max(axis=1)/sqrt(maximum(self.area/pi,sys.float_info.epsilon)) #ratio
           sindw=nonzero(srat>threshold)[0]

        if angle_min is not None: #check minimum angle
           if not hasattr(self,'angles'): self.compute_angle()
           a=self.angles; a[a<0]=999; sindw=nonzero(a.min(axis=1)<=angle_min)[0]

        #combine and save
        sindw=array(sindw)
        if fname is not None: C=schism_bpfile(); C.x,C.y,C.z=self.xctr[sindw],self.yctr[sindw],self.dpe[sindw]; C.save(fname)
        return sindw

    def inside_elem(self,pxy,ie):
        '''
        check whether pts are inside elements, then compute area coordinates for pts in elements
           pxy: c_[x,y]
           ie: array of element indices corresponding to each pt
        '''
        sind=[]; pip=[]; pacor=[]
        for i in arange(self.i34.max()-2):
            #get pts and element info
            if i==0: ip=self.elnode[ie,:3]; x1,x2,x3=self.x[ip].T; y1,y2,y3=self.y[ip].T; xi,yi=pxy.T
            if i==1:
                fpr=(~fps)*(self.i34[ie]==4); sindr=nonzero(fpr)[0];
                ip=self.elnode[ie[fpr]][:,array([0,2,3])]; x1,x2,x3=self.x[ip].T; y1,y2,y3=self.y[ip].T; xi,yi=pxy[fpr].T

            #compute area coordinates
            A0=signa(c_[x1,x2,x3],c_[y1,y2,y3]); A1=signa(c_[xi,x2,x3],c_[yi,y2,y3])
            A2=signa(c_[x1,xi,x3],c_[y1,yi,y3]); A3=signa(c_[x1,x2,xi],c_[y1,y2,yi])
            fps=(A1>=0)*(A2>=0)*(A3>=0); ac1=A1[fps]/A0[fps]; ac2=A2[fps]/A0[fps]
            if not isinstance(fps,np.ndarray): fps=array([fps])

            #get index of pts
            if i==0: sind.extend(nonzero(fps)[0])
            if i==1: sind.extend(sindr[fps])
            pip.extend(ip[fps]); pacor.extend(c_[ac1,ac2,1-ac1-ac2])
        return array(sind),array(pip),array(pacor)

    def inside_grid(self,pxy):
        '''
        check whether pts are inside grid
        usage:
            sind=gd.inside_grid(pxy)
            sind=0: outside; sind=1: inside
        '''
        npt=len(pxy); sind=zeros(npt).astype('int')
        if not hasattr(self,'bndinfo'): self.compute_bnd()
        if not hasattr(self.bndinfo,'nb'): self.compute_bnd()
        for i in arange(self.bndinfo.nb): #inside outline
            if self.bndinfo.island[i]==1: continue
            fpb=self.bndinfo.ibn[i].astype('int'); fp=inside_polygon(pxy,self.x[fpb],self.y[fpb])==1; sind[fp]=1
        for i in arange(self.bndinfo.nb): #inside island
            if self.bndinfo.island[i]==0: continue
            fpb=self.bndinfo.ibn[i].astype('int'); fp=inside_polygon(pxy,self.x[fpb],self.y[fpb])==1; sind[fp]=0
        return sind

    def subset(self,xy):
        '''
        return a new grid for elements inside a polygon defined by xy
        xy:  subdomin region (c_[x,y], or reg file)
        '''
        return get_schism_grid_subdomain(self,xy)

    def write_shp(self,fname,fmt=0,prj='epsg:4326'):
        '''
        generic function to write grid elem/node/bnd as shapefile
        fmt=0: elem (default);   fmt=1: node;    fmt=2: bnd
        '''
        if fmt==0: self.write_shapefile_element(fname,prj)
        if fmt==1: self.write_shapefile_node(fname,prj)
        if fmt==2: self.write_shapefile_bnd(fname,prj)

    def write_shapefile_bnd(self,fname,prj='epsg:4326'):
        self.shp_bnd=zdata()
        self.shp_bnd.type='POLYLINE'; xy=array([[],[]]).T
        for i in arange(self.nob):
            ind=self.iobn[i]
            xyi=c_[self.x[ind],self.y[ind]];
            xyi=insert(xyi,0,nan,axis=0);
            xy=r_[xy,xyi]
        for i in arange(self.nlb):
            ind=self.ilbn[i]
            xyi=c_[self.x[ind],self.y[ind]];
            if self.island[i]==1: xyi=close_data_loop(xyi)
            xyi=insert(xyi,0,nan,axis=0)
            xy=r_[xy,xyi]
        self.shp_bnd.xy=xy
        self.shp_bnd.prj=get_prj_file(prj)
        write_shapefile_data(fname,self.shp_bnd)

    def write_shapefile_node(self,fname,prj='epsg:4326'):
        self.shp_node=zdata()
        self.shp_node.type='POINT'
        self.shp_node.xy=c_[self.x,self.y]
        self.shp_node.attname=['id_node']
        self.shp_node.attvalue=arange(self.np)+1;
        self.shp_node.prj=get_prj_file(prj)
        write_shapefile_data(fname,self.shp_node)

    def write_shapefile_element(self,fname,prj='epsg:4326'):
        self.shp_elem=zdata()
        self.shp_elem.type='POLYGON'
        elnode=self.elnode; fp=elnode[:,-1]<0; elnode[fp,-1]=elnode[fp,0]
        elnode=fliplr(elnode)
        for i in arange(4):
            xyi=c_[self.x[elnode[:,i]],self.y[elnode[:,i]]]
            if i==0:
                xy=xyi[:,:,None]
            else:
                xy=c_[xy,xyi[:,:,None]]
        xy=transpose(xy,[0,2,1]);
        self.shp_elem.xy=zeros(self.ne).astype('O')
        for i in arange(self.ne):
            self.shp_elem.xy[i]=xy[i]

        self.shp_elem.attname=['id_elem']
        self.shp_elem.attvalue=arange(self.ne)+1;
        self.shp_elem.prj=get_prj_file(prj)
        write_shapefile_data(fname,self.shp_elem)

    def create_bnd(self):
        '''
        create open and land boundaries for grid
        '''
        def connect_actions():
            self.cidbnd=gcf().canvas.mpl_connect('button_press_event', onclick)
            if not hasattr(S,'nb'): self.compute_bnd()
            if not hasattr(S,'hb0'): S.hb0=[plot(self.x[r_[i,i[0]]],self.y[r_[i,i[0]]],'b',lw=0.5) for i in S.ibn]
            acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]); ac=acs[nonzero(ats=='Pan')[0][0]]
            if not ac.isChecked(): ac.trigger()
            gcf().canvas.draw()

        def onclick(sp):
            dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
            #double click
            if dlk==1 and btn==1: add_pt(bx,by)
            if dlk==1 and btn==3: remove_pt(bx,by)
            if dlk==1 and btn==2:
               if S.npt%2==1: print('open boundary needs to be defined'); return #don't allow finish

               #add a new land bnd to the end of the segment
               bid=S.bid[-1]; ip=S.pt[-1]
               #if S.nlb<S.nob:
               if S.ibn[bid][0]!=ip:
                  bid=S.bid[-1]; pid=nonzero(S.ibn[bid]==S.pt[-1])[0][0]
                  S.nlb=S.nlb+1; ibni=r_[S.ibn[bid][pid:],S.ibn[bid][0]]; S.ilbn.append(ibni)
                  hlb=plot(self.x[ibni],self.y[ibni],'g-'); S.hb.append(hlb); S.ihb.append(0)

               #save boundary information
               self.nob=S.nob; self.iobn=array(S.iobn,dtype='O'); self.nobn=array([len(i) for i in self.iobn])
               sid=setdiff1d(unique(S.sind),unique(array(S.bid)))
               self.nlb=S.nlb+len(sid); self.ilbn=array([*S.ilbn,*[S.ibn[i] for i in sid]],dtype='O')
               self.nlbn=array([len(i) for i in self.ilbn]); self.island=r_[tile(0,S.nlb),tile(1,len(sid))]

               #finish
               gcf().canvas.mpl_disconnect(self.cidbnd)
               acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]); ac=acs[nonzero(ats=='Pan')[0][0]]
               if ac.isChecked(): ac.trigger()
               gcf().canvas.draw()

        def add_pt(x,y):
            distp=squeeze(abs((S.x-x)+1j*(S.y-y))); sid=nonzero(distp==distp.min())[0][0]
            ip=S.ip[sid]; bid=S.sind[sid]; pid=nonzero(S.ibn[bid]==ip)[0][0]
            if S.npt!=0:
               bid0=S.bid[-1]; pid0=nonzero(S.ibn[bid0]==S.pt[-1])[0][0]
               if S.npt%2==1 and bid!=bid0: return  #two pts are not on the same boundary
               if S.npt%2==1 and (pid0>=pid and pid!=0): return  #the 2nd pt is ahead of the 1st pt (except open bnd to 1st node)
            if bid not in S.bid: S.ibn[bid]=r_[S.ibn[bid][pid:],S.ibn[bid][:pid]] #reorder boundary points

            #new bnd pt
            S.pt.append(ip); S.bid.append(bid); S.npt=S.npt+1
            hp=plot(self.x[ip],self.y[ip],'ro'); S.hp.append(hp)

            #new open bnd
            if S.npt%2==0:
               S.nob=S.nob+1; ibni=r_[S.ibn[bid][pid0:],S.ibn[bid][0]] if pid==0 else S.ibn[bid][pid0:(pid+1)]
               hob=plot(self.x[ibni],self.y[ibni],'r-'); S.hb.append(hob); S.ihb.append(1); S.iobn.append(ibni)

            #new land bnd
            if S.npt>2 and S.npt%2==1 and bid0==bid and pid!=pid0:
               S.nlb=S.nlb+1; ibni=S.ibn[bid][pid0:(pid+1)]; S.ilbn.append(ibni)
               hlb=plot(self.x[ibni],self.y[ibni],'g-'); S.hb.append(hlb); S.ihb.append(0)

            #add a new land bnd to the end of the segment
            if S.npt>=2 and bid0!=bid and pid0!=0:
               S.nlb=S.nlb+1; ibni=r_[S.ibn[bid0][pid0:],S.ibn[bid0][0]]; S.ilbn.append(ibni)
               hlb=plot(self.x[r_[ibni,S.ibn[bid0][0]]],self.y[r_[ibni,S.ibn[bid0][0]]],'g-'); S.hb.append(hlb); S.ihb.append(0)
            gcf().canvas.draw()

        def remove_pt(x,y):
            if S.npt==0: return
            bid=S.bid[-1]; pid=nonzero(S.ibn[bid]==S.pt[-1])[0][0]

            #remove 1st land bnd
            if len(S.hb)!=0:
               if S.ihb[-1]==0 and S.ilbn[-1][-1]!=S.pt[-1]: S.hb[-1][0].remove(); S.hb.pop(); S.ihb.pop(); S.ilbn.pop(); S.nlb=S.nlb-1

            #remove bnd pt
            S.hp[-1][0].remove(); S.hp.pop(); S.pt.pop(); S.bid.pop(); S.npt=S.npt-1

            #remove open bnd
            if len(S.hb)!=0:
               if S.npt%2==1 and S.ihb[-1]==1: S.hb[-1][0].remove(); S.hb.pop(); S.ihb.pop(); S.nob=S.nob-1; S.iobn.pop()

            #remove land bnd
            if len(S.hb)!=0:
               if S.ihb[-1]==0 and S.ilbn[-1][-1]!=S.pt[-1]: S.hb[-1][0].remove(); S.hb.pop(); S.ihb.pop(); S.ilbn.pop(); S.nlb=S.nlb-1
            gcf().canvas.draw()

        #add bnd icon
        if mpl._pylab_helpers.Gcf.get_active() is None: self.plot()
        acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
        abn=acs[nonzero(ats=='bnd')[0][0]] if 'bnd' in ats else gcf().canvas.toolbar.addAction('bnd')

        #add bndinfo capsule
        if not hasattr(self,'bndinfo'): self.bndinfo=zdata()
        S=self.bndinfo; S.hp=[]; S.hb=[]; S.ihb=[]; S.nob=0; S.iobn=[]; S.nlb=0; S.ilbn=[]; S.npt=0; S.pt=[]; S.bid=[]

        #connect to actions
        abn.triggered.connect(connect_actions)
        gcf().canvas.draw()

    def query_pt(self):
        '''
        add function for querying depth
        '''
        def connect_actions():
            self.cidquery=gcf().canvas.mpl_connect('button_press_event', onclick)

        def onclick(sp):
            dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
            if dlk==0 and btn==1:
               acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]);ac=acs[nonzero(ats=='bp')[0][0]]
               if hasattr(ac,'bp'):
                  if ac.bp.nsta==0: return
                  distp=squeeze(abs((ac.bp.x-bx)+1j*(ac.bp.y-by))); sid=nonzero(distp==distp.min())[0][0]
                  pie,pip,pacor=self.compute_acor(c_[ac.bp.x[sid],ac.bp.y[sid]]); pzi=(self.dp[pip]*pacor).sum()
                  print('query: bp depth= {}'.format(pzi))
            elif dlk==0 and btn==3:
               pie,pip,pacor=self.compute_acor(c_[bx,by]); pzi=(self.dp[pip]*pacor).sum()
               print('query: depth= {}'.format(pzi))
            elif dlk==0 and btn==2:
               gcf().canvas.mpl_disconnect(self.cidquery)

        acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
        abp=acs[nonzero(ats=='query')[0][0]] if 'query' in ats else gcf().canvas.toolbar.addAction('query')
        #if not abp.isCheckable(): abp.setCheckable(True)
        abp.triggered.connect(connect_actions)

    def show_node(self):
        '''
        show node/element number
        '''
        def _show_node():
            if not hasattr(self,'dpe'): self.compute_ctr()
            if len(self.hts)==0:
               xm=xlim(); ym=ylim()
               sind=nonzero((self.x>=xm[0])*(self.x<=xm[1])*(self.y>=ym[0])*(self.y<=ym[1]))[0]
               for i in sind: ht=text(self.x[i],self.y[i],'{}'.format(i+1),fontsize=6); self.hts.append(ht)
               sind=nonzero((self.xctr>=xm[0])*(self.xctr<=xm[1])*(self.yctr>=ym[0])*(self.yctr<=ym[1]))[0]
               for i in sind: ht=text(self.xctr[i],self.yctr[i],'{}'.format(i+1),fontsize=6); self.hts.append(ht)
            else:
               for i in arange(len(self.hts)): self.hts.pop().remove()
            gcf().canvas.draw()
        acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]); self.hts=[]
        abp=acs[nonzero(ats=='node')[0][0]] if 'node' in ats else gcf().canvas.toolbar.addAction('node')
        abp.triggered.connect(_show_node)

class schism_bpfile:
    def __init__(self,x=None,y=None,z=None,station=None,fmt=0):
        self.nsta=0; self.x=array([]); self.y=array([]); self.z=array([])
        self.station=[]; self.hp=[]; self.ht=[]; self.fmt=fmt
        if x is not None: self.x=x
        if y is not None: self.y=y
        if z is not None: self.z=z
        if station is not None: self.station=station
        self.check(); self.edit()

    @property
    def INFO(self):
        return get_INFO(self)
    @property
    def VINFO(self):
        return get_INFO(self)

    def read_reg(self,*args0,**args):
        self.read(fname,*args0,**args)

    def read_bpfile(self,*args0,**args):
        self.read(fname,*args0,**args)

    def read(self,fname,delimiter='!',fmt=None):
        '''
        read SCHISM bpfile or reg file
        fmt=0: bpfile;  fmt=1: regfile
        '''
        #pre-proc
        lines=[i.strip() for i in open(fname,'r').readlines()]
        self.fmt=1 if len(lines[2].split())==2 else 0 #check file format
        if (fmt is not None) and self.fmt!=fmt: sys.exit('fmt is wrong; fmt is not need anymore!')

        if self.fmt==0: #bpfile
           self.nsta=int(lines[1].split()[0])
           snames=[i.split(delimiter)[-1] if (delimiter in i) else ''  for i in lines[2:]] #read stations
           self.station=[] if (len(unique(snames))==1 and snames[0]=='') else array(snames,dtype='U')
           self.x,self.y,self.z=array([i.split()[1:4] for i in lines[2:(2+self.nsta)]]).astype('float').T
        else: #regfile
            self.nsta=int(lines[2].split()[0])
            self.x,self.y=array([i.split()[:2] for i in lines[3:(3+self.nsta)]]).astype('float').T
        self.check()

    def write(self,fname,**args):
        '''
        generic fun in saving file in different format (*.bp, *.reg, *.shp)
        when other format is provided, output as *.bp
        '''
        F=None
        if fname.endswith('.npz') or fname.endswith('.pkl'): savez(fname,self,**args); return
        if fname.endswith('.reg'): F=self.write_reg
        if fname.endswith('.shp'): F=self.write_shapefile; fname=fname[:-4]
        if fname.endswith('.bp') or (F is None):  F=self.write_bpfile
        F(fname,**args)

    def save(self,fname,**args):
        '''
        alias to generic function self.write
        '''
        self.write(fname,**args)

    def write_reg(self,fname):
        self.write_bpfile(fname,fmt=1)

    def write_bpfile(self,fname,fmt=0):
        '''
        fmt=0: write ACE/gredit *.bp file
        fmt=1: write ACE/gredit *.reg file
        '''

        self.check(); fid=open(fname,'w+')
        #write header
        if hasattr(self,'note'): fid.write('ACE/gredit: {}'.format(self.note))
        if fmt==0: fid.write('bpfile in ACE/gredit format\n{}\n'.format(self.nsta))
        if fmt==1: fid.write('Region in ACE/gredit format\n1\n{} 1\n'.format(self.nsta))

        #write pts
        for i in arange(self.nsta):
            if fmt==0: fid.write('{:<d} {:<.8f} {:<.8f} {:<.8f} !{}\n'.format(i+1,self.x[i],self.y[i],self.z[i],self.station[i]))
            if fmt==1: fid.write('{:<.8f} {:<.8f}\n'.format(self.x[i],self.y[i]))
        fid.close()

    def check(self):
        '''
        fill up missing field in bpfile based on self.x
        '''
        npt=len(self.x); self.nsta=npt; self.npt=npt; nz0=len(self.z); nsta0=len(self.station)
        if nz0<npt: self.z=r_[array(self.z),zeros(npt-nz0)]
        if nsta0<npt: self.station=r_[array(self.station,dtype='U'),arange(nsta0+1,npt+1).astype('U')]

    def get_unique_pts(self,fmt=0):
        '''
        compute unique pts
          fmt=0: unique xy; fmt=1: unique xyz
        '''
        sind=sort(unique(c_[self.x,self.y] if fmt==0 else c_[self.x,self.y,self.z],axis=0,return_index=True)[1])
        self.ux,self.uy,self.uz,self.ustation=self.x[sind],self.y[sind],self.z[sind],self.station[sind]
        self.unsta=self.unpt=len(sind)
        return [self.ux,self.uy,self.uz,self.ustation]

    def proj(self,prj0,prj1='epsg:4326',fmt=0,lon0=None,lat0=None):
        '''
        transform the projection of bpfile points' coordinates
        Inputs:
            prj0: projection name of bpfile
            prj1: target projection name; default is 'epsg:4326'
            fmt=0: change bp.x,bp.y to transform xy; fmt=1: only return transformed xy
            lon0,lat0: lon&lat of cpp projection center; needed only if 'cpp' in [prj0,prj1]
                       if ("ll"=>"cpp") and (lon0 or lat0 is not provided): lon0=mean(x); lat0=mean(y)
        '''
        self.check(); px,py=proj(prj0=prj0,prj1=prj1,x=self.x,y=self.y,lon0=lon0,lat0=lat0)
        if fmt==0: self.x,self.y=px,py
        return [px,py]

    def inside(self,xy,fmt=0):
        '''
        check whether pts c_[x,y] are inside the polygon of bp points.
        fmt=0: return indices of pts inside; fmt=1: return boolean flag
        fmt=2: return indices of pts outside region; fmt=3: return boolean flag outside
        '''

        fp=inside_polygon(xy,self.x,self.y)==1
        return nonzero(fp)[0] if fmt==0 else fp if fmt==1 else nonzero(~fp)[0] if fmt==2 else ~fp
    def outside(self,xy,fmt=2):
        return self.inside(xy,fmt=fmt)

    def write_shapefile(self,fname,prj='epsg:4326'):
        self.shp_bp=zdata()
        self.shp_bp.type='POINT'
        self.shp_bp.xy=c_[self.x,self.y]
        self.shp_bp.prj=get_prj_file(prj)

        if hasattr(self,'station'):
            self.shp_bp.attname=['station']
            self.shp_bp.attvalue=self.station
        write_shapefile_data(fname,self.shp_bp)

    def plot_station(self,**args):
        return self.plot(**args)

    def plot(self,ax=None,color=None,ls=None,label=None,marker='.',markersize=6,fmt=0,connect_mpl=1,**args):
        '''
        plot points on current figure
          fmt=0: plot all points; fmt=1: plot unique points
        '''
        #pre-processing
        if connect_mpl==1: self.edit()
        if color is None: color='r'
        if ls is None: ls='None' if self.fmt==0 else '-'
        if label is None: label=1
        if ax is None: ax=gca()

        #get pts
        self.check(); npt,xs,ys,zs,stations=self.npt,self.x,self.y,self.z,self.station
        if fmt==1: xs,ys,zs,stations=self.get_unique_pts(); npt=len(xs)

        #label
        if label==1:
           [i.remove() for i in self.ht]; self.ht=[]
           for xi,yi,ti in zip(xs,ys,stations): hti=text(xi,yi,ti,color=color,fontsize=10); self.ht.append(hti)

        #plot
        if self.fmt==1: #for region file
           if hasattr(self,'hp2'): self.hp2.remove(); delattr(self,'hp2')
           if npt>=3: [self.hp2]=fill(xs,ys,facecolor='m',alpha=0.15)
           if npt!=0: xs=r_[xs,xs[0]]; ys=r_[ys,ys[0]] #close polygon
        if len(self.hp)!=0: self.hp[0].set_xdata(xs); self.hp[0].set_ydata(ys)
        if len(self.hp)==0: self.hp=plot(xs,ys,marker=marker,markersize=markersize,color=color,linestyle=ls,**args)
        gcf().canvas.draw()

    def compute_acor(self,gd):
        #compute areal coordinates, and gd is the schism grid
        self.ie,self.ip,self.acor=gd.compute_acor(c_[self.x,self.y])
        return self.ie,self.ip,self.acor

    def disconnect_edit(self):
        if hasattr(self,'cidpress'): gcf().canvas.mpl_disconnect(self.cidpress)
        if hasattr(self,'cidmove'):  gcf().canvas.mpl_disconnect(self.cidmove)
        acs=gcf().canvas.toolbar.actions(); ats=[i.iconText() for i in acs]; ap=acs[ats.index('Pan')]
        if ap.isChecked(): ap.trigger()
        gcf().canvas.draw()

    def edit(self):
        def connect_actions():
            self.cidmove=gcf().canvas.mpl_connect('motion_notify_event', onmove)
            self.cidpress=gcf().canvas.mpl_connect('button_press_event', onclick)
            if self.nsta!=0 and len(self.hp)==0: self.plot_station()
            acs=gcf().canvas.toolbar.actions(); ats=[i.iconText() for i in acs]
            ap=acs[ats.index('Pan')]; ab=acs[ats.index('reg' if self.fmt==0 else 'bp')]
            if not ap.isChecked(): ap.trigger()
            if hasattr(ab,'bp'): ab.bp.disconnect_edit()
            gcf().canvas.draw()
            print('double click: left=add pt, right=remove pt; middle=finish edit')
            print('single click: middle=move pt')

        def onmove(sp):
            if sp.button is not None:
               dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
               if dlk==0 and btn==2: move_pt(bx,by)

        def onclick(sp):
            dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
            #double click
            if dlk==1 and btn==1: add_pt(bx,by)
            if dlk==1 and btn==3: remove_pt(bx,by)
            if dlk==1 and btn==2: self.disconnect_edit()
            if dlk==1 and btn in [1,3]: self.station=(arange(len(self.x))+1).astype('U'); self.plot(connect_mpl=0)

        def add_pt(x,y):
            self.x=r_[self.x,x]; self.y=r_[self.y,y]; self.plot(connect_mpl=0)

        def remove_pt(x,y):
            if self.nsta==0: return
            distp=squeeze(abs((self.x-x)+1j*(self.y-y))); sid=nonzero(distp==distp.min())[0][0]
            self.x=r_[self.x[:sid],self.x[(sid+1):]]; self.y=r_[self.y[:sid],self.y[(sid+1):]]; self.plot()

        def move_pt(xi,yi):
            distp=squeeze(abs((self.x-xi)+1j*(self.y-yi))); sid=nonzero(distp==distp.min())[0][0]
            self.x[sid]=xi; self.y[sid]=yi; self.plot()

        if mpl._pylab_helpers.Gcf.get_active() is not None:
            acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]); at='bp' if self.fmt==0 else 'reg'
            abp=acs[nonzero(ats==at)[0][0]] if at in ats else gcf().canvas.toolbar.addAction(at)
            #if not abp.isCheckable(): abp.setCheckable(True)

            #disconnect and clean previous bpfile
            if hasattr(abp,at):
               if self is not abp.bp:
                  nhp=len(abp.bp.hp)
                  for i in arange(nhp):
                      abp.bp.hp[-1][0].remove(); abp.bp.ht[-1].remove()
                      del abp.bp.hp[-1],abp.bp.ht[-1]
               abp.triggered.disconnect()

            #connect to new object
            abp.triggered.connect(connect_actions); abp.bp=self
            gcf().canvas.draw()

def read_schism_hgrid(fname):
    gd=schism_grid(); gd.read_hgrid(fname)
    return gd

def read_schism_bpfile(fname,fmt=0):
    '''
    read schism *bp (fmt=0) or *.reg (fmt=1) file created by ACE/gredit
    '''
    bp=schism_bpfile(); bp.read(fname); bp.edit()
    return bp

def read_schism_reg(fname):
    '''
    read schism *bp (fmt=0) or *.reg (fmt=1) file created by ACE/gredit
    '''
    bp=schism_bpfile(); bp.read(fname); bp.edit()
    return bp

def read_schism_prop(fname):
    '''
    read schism *.prop file (element based), and return the values
    '''
    try:
      pdata=loadtxt(fname)
    except:
      pdata=array([i.strip().split() for i in open(fname,'r').readlines() if len(i.split())==2]).astype('float')
    pvalue=pdata[:,1] if pdata.ndim==2 else pdata[None,:][:,1]
    return pvalue

def save_schism_grid(fmt=0,path='.',method=0):
    '''
    save schism grid information in *.npz format (hgrid.gr3,hgrid.ll,vgrid.in}
       path:  directory of grids
       fmt=0: save as grid.npz (include hgrid and vgrid); fmt=1: save as hgrid.npz and vgrid.npz
       method=1: save hgrid full geometry information
    '''
    #read grids
    grd=path+'/hgrid.gr3'; grd0=path+'/hgrid.ll'; vrd='{}/vgrid.in'.format(path)
    gd=None; gd0=None; vd=None
    if os.path.exists(grd):  gd=read_schism_hgrid(grd)
    if os.path.exists(grd0): gd0=read_schism_hgrid(grd0)
    if os.path.exists(vrd):  vd=read_schism_vgrid(vrd)
    if (gd is None) and (gd0 is not None): gd=gd0
    if (gd is not None) and (gd0 is not None): gd.lon,gd.lat=gd0.x,gd0.y
    if (gd is not None) and method==1: gd.compute_all()

    #save grid
    if fmt==1:
       if gd is not None: gd.save('hgrid.npz')
       if vd is not None: vd.save('vgrid.npz')
    else:
       S=zdata()
       if gd is not None: S.hgrid=gd
       if vd is not None: S.vgrid=vd
       S.save('grid.npz')

class schism_vgrid:
    def __init__(self):
        pass

    @property
    def INFO(self):
        return get_INFO(self)
    @property
    def VINFO(self):
        return get_INFO(self)

    def read_vgrid(self,fname):
        #read schism vgrid
        fid=open(fname,'r'); lines=fid.readlines(); fid.close()

        self.ivcor=int(lines[0].strip().split()[0]); self.nvrt=int(lines[1].strip().split()[0])
        if self.ivcor==1:
            #read vgrid info
            lines=lines[2:]; sline=array(lines[0].split()).astype('float')
            if sline.min()<0: #old format
               self.kbp=array([int(i.split()[1])-1 for i in lines]); self.np=len(self.kbp)
               self.sigma=-ones([self.np,self.nvrt])
               for i,line in enumerate(lines):
                   self.sigma[i,self.kbp[i]:]=array(line.strip().split()[2:]).astype('float')
            else:
              sline=sline.astype('int'); self.kbp=sline-1; self.np=len(sline)
              self.sigma=array([i.split()[1:] for i in lines[1:]]).T.astype('float')
              fpm=self.sigma<-1; self.sigma[fpm]=-1
        elif self.ivcor==2:
            self.kz,self.h_s=lines[1].strip().split()[1:3]; self.kz=int(self.kz); self.h_s=float(self.h_s)

            #read z grid
            self.ztot=[]; irec=2
            for i in arange(self.kz):
                irec=irec+1
                self.ztot.append(lines[irec].strip().split()[1])
            self.ztot=array(self.ztot).astype('float')

            #read s grid
            self.sigma=[]; irec=irec+2
            self.nsig=self.nvrt-self.kz+1
            self.h_c,self.theta_b,self.theta_f=array(lines[irec].strip().split()[:3]).astype('float')
            for i in arange(self.nsig):
                irec=irec+1
                self.sigma.append(lines[irec].strip().split()[1])
            self.sigma=array(self.sigma).astype('float')
        return self.sigma

    def compute_zcor(self,dp,eta=0,fmt=0,method=0,sigma=None,kbp=None,ifix=0):
        '''
        compute schism zcor (ivcor=1 and 2)
            dp:  depth at nodes (dim=[np] or [1])
            eta: surface elevation (dim=[np] or [1])
            fmt: output format of zcor
                 fmt=0: bottom depths byeond kbp are extended
                 fmt=1: bottom depths byeond kbp are nan
            method=1 and ivcor=1: used for computing zcor for subset of nodes (need sigma,kbp)
            method=1 and ivcor=2: return zcor and kbp
            ifix=1 and ivcor=2: using traditional sigma in shallow if error raise
        '''
        if self.ivcor==1:
           if method==0: return compute_zcor(self.sigma,dp,eta=eta,fmt=fmt,kbp=self.kbp)
           if method==1: return compute_zcor(sigma,dp,eta=eta,fmt=fmt,kbp=kbp)
        elif self.ivcor==2:
           zcor,kbp=compute_zcor(self.sigma,dp,eta=eta,fmt=fmt,ivcor=2,vd=self,method=1,ifix=ifix)
           if method==0: return zcor
           if method==1: return [zcor,kbp]
    def save(self,fname=None,fmt=0,**args):
        '''
        alias to write_vgrid
        '''
        self.write_vgrid(fname,fmt,**args)

    def write_vgrid(self,fname='vgrid.in',fmt=0,**args):
        '''
        write schism vertical grid
            fmt=0: write vgrid.in in latest format of ivcor=1 (one line per lelvel)
            fmt=1: write vgrid.in in old format of ivcor=1    (one line per node)
        '''
        if fname.endswith('.npz') or fname.endswith('.pkl'): savez(fname,self,**args); return
        if self.ivcor==1:
           nvrt,np,kbp,sigma=self.nvrt,self.np,self.kbp.copy(),self.sigma.copy()
           fid=open(fname,'w+'); fid.write('1    !average # of layers={}\n{}  \n'.format(mean(nvrt-kbp),nvrt))
           if fmt==0:
              for i in arange(np): sigma[i,:kbp[i]]=-9
              fstr='    '+' {:10d}'*np+'\n'; kbp=kbp+1; fid.write(fstr.format(*kbp))
              fstr='{:8d}'+' {:10.6f}'*np+'\n'; sigma=sigma.T
              [fid.write(fstr.format(i+1,*k)) for i,k in enumerate(sigma)]
           elif fmt==1:
              for i,[k,sigma] in enumerate(zip(kbp,sigma)):
                  fstr='{:9d} {:3d}'+' {:11.6f}'*(nvrt-k)+'\n'
                  fid.write(fstr.format(i+1,k+1,*sigma[k:]))
           fid.close()
        elif self.ivcor==2:
           fid=open(fname,'w+'); fid.write('2  !ivcor\n')
           fid.write('{} {} {} !nvrt, kz, h_s \nZ levels\n'.format(self.nvrt,self.kz,self.h_s))
           for k,zlevel in enumerate(self.ztot): fid.write('{} {}\n'.format(k+1,zlevel))
           fid.write('S levels\n{} {} {} !h_c, theta_b, theta_f\n'.format(self.h_c,self.theta_b,self.theta_f))
           for k,slevel in enumerate(self.sigma): fid.write('{} {:9.6f}\n'.format(k+1,slevel))
           fid.close()
        else:
           sys.exit('unknow ivcor={}'.format(self.ivcor))

def read_schism_vgrid(fname):
    '''
    read schism vgrid information
    '''
    vd=schism_vgrid(); vd.read_vgrid(fname)
    return vd

def compute_zcor(sigma,dp,eta=0,fmt=0,kbp=None,ivcor=1,vd=None,method=0,ifix=0):
    '''
    compute schism zcor (ivcor=1 and 2)
        sigma: sigma cooridinate (dim=[np,nvrt])
        dp: depth at nodes (dim=[np] or [1])
        eta: surface elevation (dim=[np] or [1])
        fmt: output format of zcor
            fmt=0: bottom depths byeond kbp are extended
            fmt=1: bottom depths byeond kbp are nan
        kbp: index of bottom layer (not necessary, just to speed up if provided for ivcor=1)
        method=1 and ivcor=2: return zcor and kbp
        ifix=1 and ivcor=2: using traditional sigma in shallow if error raise

    usage:
      for ivcor=1: compute_zcor(sigma,dp,eta)
      for ivcor=2: compute_zcor(sigma,dp,eta,vd=vd,ivcor=2)
    '''

    if ivcor==1:
        np=sigma.shape[0]
        if not hasattr(dp,'__len__'):  dp=ones(np)*dp
        if not hasattr(eta,'__len__'): eta=ones(np)*eta

        #get kbp
        if kbp is None:
            kbp=array([nonzero(abs(i+1)<1e-10)[0][-1] for i in sigma])

        #thickness of water column
        hw=dp+eta; hw[hw<0]=0

        #add elevation
        zcor=hw[:,None]*sigma+eta[:,None]
        fpz=hw<0; zcor[fpz]=-dp[fpz][:,None]

        #change format
        if fmt==1:
            for i in arange(np):
                zcor[i,:kbp[i]]=nan
        return zcor
    elif ivcor==2:
        #arange data
        np,dp=[1,array([dp])] if not hasattr(dp,'__len__') else [len(dp),dp]
        eta=ones(np)*eta if not hasattr(eta,'__len__') else eta

        #sz parameter
        nvrt,sigma,theta_b,theta_f,h_c,h_s,kz,ztot=vd.nvrt,vd.sigma,vd.theta_b,vd.theta_f,vd.h_c,vd.h_s,vd.kz,vd.ztot
        cs=(1-theta_b)*sinh(theta_f*sigma)/sinh(theta_f)+ theta_b*(tanh(theta_f*(sigma+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)

        #for sigma layer
        zcor=zeros([nvrt,np]); hmod=dp.copy(); hmod[hmod>h_s]=h_s
        for k in arange(kz-1,nvrt):
            kin=k-kz+1
            #hmod<=h_c
            fps=hmod<=h_c; htot=hmod[fps]+eta[fps]; htot[htot<0]=0
            zcor[k,fps]=sigma[kin]*(htot)+eta[fps]

            #hmod>h_c
            sindc=nonzero(~fps)[0]; fpc=eta[sindc]<=(-h_c-(hmod[sindc]-h_c)*theta_f/sinh(theta_f)); sindf=sindc[fpc]; sinds=sindc[~fpc]
            if sum(fpc)>0 and ifix==0: sys.exit('Pls choose a larger h_c: {}'.format(h_c))
            if sum(fpc)>0 and ifix==1: zcor[k:,sindf]=(eta[sindf]+hmod[sindf])*sigma[kin]+eta[sindf]
            zcor[k,sinds]=eta[sinds]*(1+sigma[kin])+h_c*sigma[kin]+cs[kin]*(hmod[sinds]-h_c)
        zcor[:(kz-1),:]=nan if fmt==1 else zcor[kz-1,:][None,:] #extend

        #for z layer
        kbp=zeros(np).astype('int'); fp=dp>h_s; sind=nonzero(fp)[0]; kbp[~fp]=kz-1
        for k in arange(kz-1):
            fpz=(-dp[sind]>=ztot[k])*(-dp[sind]<ztot[k+1]); sindz=sind[fpz]; kbp[sindz]=k
            zcor[:(k+1),sindz]=-dp[sindz][None,:]; zcor[(k+1):(kz-1),sindz]=ztot[(k+1):(kz-1),None]
            if fmt==1: zcor[:k,sindz]=nan
        zcor=zcor.T; vd.kbp=kbp
        if method==0: return zcor
        if method==1: return [zcor,kbp]

def create_schism_vgrid(fname='vgrid.in',ivcor=2,nvrt=10,zlevels=-1.e6,h_c=10,theta_b=0.5,theta_f=1.0):
    '''
    create schism vertical grid:
        fname: name of the grid
        nvrt: number of vertical layers
        ivcor=2: SZ grid
              zlevels: 1). Z levels or 2). single number for h_s
              h_c, theta_b, theta_f:  strench constants for sigma grid
        ivcor=1: LCS^2 grid
    '''
    vd=schism_vgrid(); vd.ivcor,vd.nvrt=ivcor,nvrt
    if ivcor==2:
        if hasattr(zlevels,'__len__'):
           vd.kz,vd.ztot,vd.h_s=len(zlevels),zlevels,-zlevels[-1]
        else:
           vd.kz,vd.ztot,vd.h_s=1,[zlevels],-zlevels
        vd.h_c,vd.theta_b,vd.theta_f=h_c,theta_b,theta_f
        vd.sigma=linspace(-1,0,nvrt+1-vd.kz)
        vd.write_vgrid(fname)
    else:
        sys.exit('ivcor=1 option not available yet')

def read_schism_local_to_global(fname):
    '''
    read schism partition information
    '''
    lines=open(fname,'r').readlines()[2:]

    #get ne, np, ns, i34,elnode,
    S=zdata()
    ne=int(lines[0].strip()); np=int(lines[ne+1].strip()); ns=int(lines[ne+np+2].strip())
    S.ielg=array([i.strip().split() for i in lines[1:(ne+1)]])[:,1].astype('int')-1
    S.iplg=array([i.strip().split() for i in lines[(ne+2):(ne+np+2)]])[:,1].astype('int')-1
    S.islg=array([i.strip().split() for i in lines[(ne+np+3):(ne+np+ns+3)]])[:,1].astype('int')-1

    #find line for np,ne
    for i in arange(ne+np+ns+3,len(lines)):
        sline=lines[i].strip().split()
        if len(sline)!=2: continue
        if int(float(sline[0]))==np and int(float(sline[1]))==ne: nd=i; break;

    slines=array([i.strip().split() if len(i.split())==5 else [*i.strip().split(),'-1'] for i in lines[(nd+np+1):(nd+np+ne+1)]]).astype('int')
    i34=slines[:,0].astype('int'); elnode=slines[:,1:].astype('int')-1

    S.ne,S.np,S.ns,S.i34,S.elnode=ne,np,ns,i34,elnode
    return S

def read_schism_param(fname,fmt=0):
    '''
    read schism parameters from param.nml/param.in/cosine.in
      fmt=0: return dictionary with all field as string
      fmt=1: return dictionary with all field as float if possible
      fmt=2: return zdata with all field as string attributes
      fmt=3: return zdata with all field as float attributes if possible
    '''

    #read all lines first
    fid=open(fname,'r'); lines=[i.strip() for i in fid.readlines()]; fid.close()
    lines=[i for i in lines if ('=' in i) and (i!='') and (i[0]!='!') and (i[0]!='&')]

    #parse each line
    P={}
    for line in lines:
      if '!' in line: line=line[:line.find('!')]
      keyi,vali=line.split('='); keyi=keyi.strip(); vali=vali.strip()
      if fmt in [1,3]:  #convert string to float
         try:
            vali=[float(i) if (('.' in i) or ('e' in i) or ('E' in i)) else int(i) for i in vali.replace(',',' ').replace(';',' ').split()]
            if len(vali)==1: vali=vali[0]
         except:
            pass
      P[keyi]=vali
    param=P

    #change output format as zdata
    if fmt in [2,3]:
       S=zdata()
       for keyi,vali in P.items(): S.__dict__[keyi]=vali
       param=S

    return param

def change_schism_param(fname,param=None,value=None,source=None,note_delimiter='!'):
    '''
    change parameter values
      fname: the name of parameter file (param.nml,param.in, ...)
      param: parameter name to be changed
      value: new parameter value are to be assigned
      source: either a reference parameter file, or and object of "read_schism_param"
    '''
    def _newline(line,param,value):
       eid=line.find('='); nid=line.find(note_delimiter)
       ns1=len(line)-len(line.lstrip())
       ns2=len(line[:eid])-len(line[:eid].rstrip())
       ns3=len(line[:nid])-len(line[:nid].rstrip())
       sline=' '*ns1+param+' '*ns2+'= '+str(value)
       sline=sline+' '*ns3+note_delimiter+line[(nid+1):] if nid!=-1 else sline+'\n'
       return sline

    #read fname information
    fid=open(fname,'r'); slines=fid.readlines(); fid.close()

    #change parameter value based on refernece parameter file
    if source is not None:
       P=read_schism_param(fname,1); lines=slines[:]; slines=[]
       S=read_schism_param(source,1) if isinstance(fname,str) else source
       for line in lines:
           eid=line.find('='); pname=line[:eid].strip(); sline=line
           if (eid!=-1) and (pname in S):
              if P[pname]!=S[pname]:
                 svalue=' '.join([str(i) for i in S[pname]]) if isinstance(S[pname],list) else str(S[pname])
                 sline=_newline(line,pname,svalue)
           slines.append(sline)

    #change parameter value based on arguments
    if param is not None:
       lines=slines[:]; slines=[]
       for line in lines:
          if line[:max(line.find('='),0)].strip()==param:
             slines.append(_newline(line,param,value))
          else:
             slines.append(line)

    #write parameter file
    fid=open(fname,'w+'); fid.writelines(slines); fid.close()

def write_schism_param(fname,param):
    pkeys=sorted(param.keys())
    with open(fname,'w+') as fid:
        for i in range(len(pkeys)):
           fid.write('{:10}= {:}\n'.format(pkeys[i],param[pkeys[i]]))

def sms2grd(sms,grd=None):
    '''
      1). read SMS *2dm grid, and return grid object
      2). if grd!=None, save grid as *gr3 format
          e.g. gd=sms2gr3('hgrid.2dm','hgrid.gr3')
    '''

    #read 2dm file
    fid=open(sms,'r'); lines=fid.readlines(); fid.close()

    #for traingle and quads elements
    E3=array([ [*i.strip().split()[1:-1],'-1'] for i in lines if i.startswith('E3T')]).astype('int')
    E4=array([i.strip().split()[1:-1] for i in lines if i.startswith('E4Q')]).astype('int')
    E34=array([*E3,*E4]); sind=argsort(E34[:,0]); E34=E34[sind]

    #for nodes
    ND=array([i.strip().split()[1:] for i in lines if i.startswith('ND')]).astype('float')
    sind=argsort(ND[:,0]); ND=ND[sind]

    #save grid information
    gd=schism_grid(); gd.ne=E34.shape[0]; gd.np=ND.shape[0]
    gd.elnode=E34[:,1:]-1; gd.x,gd.y,gd.dp=ND[:,1:].T
    gd.i34=4*ones(gd.ne).astype('int'); fp3=E34[:,-1]==-1; gd.i34[fp3]=3

    if grd is not None: gd.write_hgrid(grd)
    return gd

def grd2sms(grd,sms):
    '''
      convert schism hgrid to SMS *.2dm format
      usage:
           1). grd2sms('hgrid.gr3','hgrid.2dm')
           2). grd2sms(gd,'hgrid.2dm'), or gd.grd2sms('hgrid.2dm')
           note  gd=read_schism_hgrid('hgrid.gr3')
    '''

    #read grid
    if isinstance(grd,str):
       gd=read_schism_hgrid(grd)
    elif isinstance(grd,schism_grid):
       gd=grd
    else:
       sys.exit('unknow format of grd: {}'.format(grd))

    #save grid save *2dm format
    gd.grd2sms(sms)


# ---------------------experimental functions and classes-------------------------------------

def cread_schism_hgrid(fname):
    '''
    A wrapper function to read SCHISM hgrid.gr3 file using c++ binding,
    then copy the data to pylib's grid object.
    Usefull when the grid is large and the python reading is slow.
    '''
    import hgrid_pybind

    hgrid_obj = hgrid_pybind.HGrid(fname, True)  # read c++ HGrid object from file and optionally get side information

    gd=schism_grid()  # initialize empty pylib's grid object
    gd.source_file = str(fname)

    # copy the member variables from the c++ HGrid object to pylib's grid object
    gd.np = hgrid_obj.np
    gd.ne = hgrid_obj.ne
    gd.x = hgrid_obj.x
    gd.y = hgrid_obj.y
    gd.dp = hgrid_obj.z
    gd.elnode = hgrid_obj.elements
    gd.i34 = hgrid_obj.i34
    gd.iobn = np.array(hgrid_obj.openBoundaryNodes, dtype=object)
    gd.nobn = np.array(hgrid_obj.nobn)
    gd.nob = np.array(hgrid_obj.nob)
    gd.ilbn = np.array(hgrid_obj.landBoundaryNodes, dtype=object)
    gd.nlbn = np.array(hgrid_obj.nlbn)
    gd.nlb = np.array(hgrid_obj.nlb)
    gd.island = np.array(hgrid_obj.island)
    gd.ns = hgrid_obj.ns

    return gd

def read_schism_vgrid_cached(vg_filename, overwrite_cache=False):
    vg_cache_fname = os.path.splitext(vg_filename)[0] + '.pkl'
    if not overwrite_cache:
        try:
            with open(vg_cache_fname, 'rb') as handle:
                vg = pickle.load(handle)
        except Exception as e:
            print(f'{e}\nError reading cache file {vg_cache_fname}.\nReading from original file.')
            # read original file
            vg = read_schism_vgrid(vg_filename)
            with open(vg_cache_fname, 'wb') as handle:
                pickle.dump(vg, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        vg = read_schism_vgrid(vg_filename)
        with open(vg_cache_fname, 'wb') as handle:
            pickle.dump(vg, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return vg

def combine_dataframes(A, B, weights=[1.0, 1.0]):
    import numbers

    AB = copy.deepcopy(A)

    # warn if B's time period does not contain A's time
    if B.index[0] > A.index[0] or B.index[-1] < A.index[-1]:
        print('Warning: B\'s time period does not contain A\'s time period')
        print('In the interpolated B, NaN will be filled for the time period that is not covered by B')
        print('and the NaN in the interpolated B will be treated as 0 when summing up')

    # Interpolate B to A's index
    B_interpolated = B.reindex(A.index).interpolate(method='time')

    # Find the columns that are in B but not in A
    new_columns = B_interpolated.columns.difference(A.columns)

    # append these columns to A
    AB = pd.concat([A, B_interpolated[new_columns]], axis=1)

    # Find the common columns in A and B
    common_columns = A.columns.intersection(B_interpolated.columns)

    # Sum the values from both dataframes for the common columns
    if isinstance(weights[0], numbers.Number):
        AB[common_columns] = weights[0] * A[common_columns] + weights[1] * B_interpolated[common_columns]
    elif isinstance(weights[0], TimeHistory):
        vs_A = weights[0].df[common_columns]
        vs_B = weights[1].df[common_columns]
        vs_B = vs_B.reindex(vs_A.index).interpolate(method='time')
        weights_A = vs_A / (vs_A + vs_B)
        weights_B = vs_B / (vs_A + vs_B)
        AB[common_columns] = weights_A * A[common_columns] + weights_B * B_interpolated[common_columns]
    else:
        raise ValueError('weights must be a list of two numbers or two TimeHistory objects')

    return AB

class TimeHistory:
    """Class for handling SCHISM's *.th file format.
    A *.th file is a simple ascii file containing time series data,
    the 1st column is time in seconds or days,
    the 2nd to last column is data, so each column is a time series.
    However, the *.th file itself lacks information about the time units and start time
    and this class is to bind those information in a dataframe.
    The class also sets the time column as the index of the dataframe
    and assign each data column with a meaningful name (e.g., station name)
    to facilitate data queries and processing.
    """
    @property
    def time(self):
        # original time in *th's format
        seconds = (self.df.index - self.df.index[0]).total_seconds().values
        t = seconds / self.sec_per_time_unit
        return t

    @property
    def delta_t(self): # time step
        return self.time[1] - self.time[0]

    @property
    def n_time(self):
        return self.df.shape[0]

    @property
    def n_station(self):
        return self.df.shape[1]

    @property
    def stations(self):
        return self.df.columns.values.tolist()

    @property
    def data(self):
        # original data excluding time
        return self.df.values


    def __init__(self, start_time_str="2000-01-01 00:00:00", data_array=None, columns=None, th_unit='seconds'):
        """
        Initialize a TimeHistory object from a data array,
        assuming the first column is time, and the rest are data.
        Note that you need to provide the start time and the time unit.
        """

        # list of main attributes
        self.df = pd.DataFrame(data_array)
        self.meaningful_columns = False  # some functions only work when column names are meaningful,
        self.th_unit = th_unit
                                         # e.g., station names or element ids, but not datetime + column numbers
        self.sec_per_time_unit = None

        # set time unit
        unit_dict = {'seconds': 1, 'minutes': 60, 'hours': 3600, 'days': 86400, 'weeks': 604800, 'years': 31536000}
        self.sec_per_time_unit = unit_dict[th_unit]  # only related to the time column of a *.th file

        # set column names, which usually are datetime + station ids
        if type(columns) is list:  # from a user-specified list
            if len(columns) == self.df.shape[1]:
                self.df.columns = [str(x) for x in columns]
                self.meaningful_columns = True
            elif len(columns) == self.df.shape[1]-1:
                print('number of column labels does not match the data array, assuming the first column of the data is time')
                self.df.columns = [str(x) for x in ['datetime'] + columns]
                self.meaningful_columns = True
            else:
                raise Exception('number of columns does not match')
        elif columns is None:  # first col is time and the rest are column numbers
            self.df.columns = ['datetime'] + [str(x) for x in range(1, self.df.shape[1])]
            self.meaningful_columns = False
        else:
            raise Exception('unknown columns type')

        # Lay some ground rules
        # force the first column's name to "datetime"
        self.df.rename(columns={0: 'datetime'}, inplace=True)
        second_series = self.df['datetime'].values * self.sec_per_time_unit
        # convert the first column to pandas DatetimeIndex
        time_stamps = pd.to_timedelta(second_series, unit='s') + pd.to_datetime(start_time_str)
        self.df['datetime'] = time_stamps
        # set datetime as index
        self.df.set_index('datetime', inplace=True)
        # force column names to be string
        self.df.columns = self.df.columns.map(str)

    @classmethod
    def from_file(cls, file_name, start_time_str="2000-01-01 00:00:00", th_unit='seconds', columns=None):
        """
        Initialize from a file.
        Note that the *.th file doen't have information about the time units and start time,
        and columns names, so you need to provide them.
        """
        data = np.loadtxt(file_name)
        return cls(data_array=data, th_unit=th_unit, start_time_str=start_time_str, columns=columns)

    def __getitem__(self, selector):
        """Subset the TimeHistory object by column names"""

        # parse row and col selectors from selector
        if type(selector) is str:
            selector = [selector]
        elif isinstance(selector, np.ndarray):
            if len(selector.shape) == 1:  # 1D array of column names
                selector = selector.astype(str).tolist()
            else:
                raise IndexError("Column names must be a 1D array")

        if type(selector) is list and all(isinstance(x, str) for x in selector):  # subset by column names
            column_names = selector
            subset_data = np.array(self.df[column_names])
            return TimeHistory(
                start_time_str=self.df.index[0],
                data_array=np.c_[self.time, subset_data],
                columns=column_names,
                th_unit=self.th_unit
            )
        elif isinstance(selector, tuple):  # subset by row and column indices; column index does not include time
            if len(selector) != 2:
                raise IndexError("Only 2D indexing is supported")
            row_idx, col_idx = selector
            subset_data = self.data[row_idx, col_idx]
            subset_time_str = self.df.index[row_idx]
            subset_time = self.time[row_idx]
        elif isinstance(selector, slice):  # subset by row index only, i.e., by time
            subset_df = self.df.loc[selector]
            subset_data = subset_df.values
            subset_time_str = subset_df.index
            subset_time = np.array(subset_df.index - subset_df.index[0]).astype('timedelta64[s]').astype(float) / self.sec_per_time_unit
            col_idx = range(len(self.df.columns))

            # if isinstance(row_idx, slice) or isinstance(col_idx, slice):
            #     # Handle slices here
            #     rows = self.data[row_idx] if isinstance(row_idx, int) else [self.data[i] for i in range(*row_idx.indices(len(self.data)))]
            #     if isinstance(col_idx, int):
            #         return [row[col_idx] for row in rows]
            #     else:
            #         return [row[col_idx] for row in rows for col_idx in range(*col_idx.indices(len(row)))]
            # else:
            #     # Regular integer indexing
            #     data = self.data[row_idx][col_idx]

            return TimeHistory(
                start_time_str=subset_time_str[0],
                data_array=np.c_[subset_time, subset_data],
                columns=self.df.columns[col_idx].tolist()
            )
        else:
            raise IndexError("Unknown type of index")

    def __add__(self, other, weights=[1.0, 1.0]):
        """
        Add two TimeHistory objects together
        Interpolate other to self's time stamps;
        other attributes also inherit from self.
        When combining columns with the same name,
        the default weights are 1.0 for both self and other,
        meaning that the values are summed up.
        if you want to average the values, set weights=[0.5, 0.5]
        You can also specify two TimeHistory objects as the surrogates of the weights,
        e.g., when combining two sets of msource, you can set weights=[vsource1, vsource2],
        i.e., a weighted average based on the volume of the two sources,
        """

        A = copy.deepcopy(self)
        B = other

        A.df = combine_dataframes(A.df, B.df, weights=weights)

        return A

    def __eq__(self, other) -> bool:
        """Check if two TimeHistory objects are equal"""
        for att in ['sec_per_time_unit']:
            if getattr(self, att) != getattr(other, att):
                print(f'{att} is not equal')
                return False

        # dataframes, strict test
        # return self.df.equals(other.df)

        # test if the data are equal within a tolerance
        if np.any(self.df.columns != other.df.columns):
            print('column labels are not equal')
            return False
        return np.allclose(self.df, other.df, rtol=0.001, atol=0.0001)

    def writer(self, file_name, np_savetxt_args={'fmt':'%.4f', 'delimiter':' ', 'newline':'\n'}):
        # assemble data array in *.th format and write to file
        np.savetxt(file_name, np.c_[self.time, self.data], **np_savetxt_args)


class SourceSinkIn():
    def __init__(self, ele_groups=[[], []]):
        self.ele_groups = ele_groups  # 0: source; 1: sink

    @property
    def n_group(self):
        return len(self.ele_groups)

    @property
    def np_group(self):
        return [len(x) for x in self.ele_groups]

    @property
    def ip_group(self):
        return [np.array(x) for x in self.ele_groups]

    @property
    def n_source(self):
        return self.np_group[0]

    @property
    def n_sink(self):
        return self.np_group[1]

    @classmethod
    def from_file(cls, filename):
        ele_groups = [[], []]
        with open(filename, 'r') as file:
            for k in range(0, 2):  # 0: source; 1: sink
                num_points = int(file.readline().strip().split()[0])
                for _ in range(num_points):
                    ele_groups[k].append(int(file.readline().strip().split()[0]))
                # blank line between groups
                if k == 0:
                    file.readline()
        source_sink_in = cls(ele_groups=ele_groups)
        source_sink_in.print_info()

        return source_sink_in

    def print_info(self):
        print(f"nsource: {self.n_source}")
        if self.n_source > 0:
            print(f"first and last ele: {self.ele_groups[0][0]}, {self.ele_groups[0][-1]}")

        print(f"nsink: {self.n_sink}")
        if self.n_sink > 0:
            print(f"first and last ele: {self.ele_groups[1][0]}, {self.ele_groups[1][-1]}")

    def writer(self, filename=None):
        with open(filename, 'w') as fout:
            for k in range(0, self.n_group):
                print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                fout.write(f"{self.np_group[k]}\n")
                for i in range(0, self.np_group[k]):
                    fout.write(f"{self.ele_groups[k][i]}\n")
                fout.write("\n")  # empty line

    def __eq__(self, other) -> bool:
        for g1, g2 in zip(self.ele_groups, other.ele_groups):
            if g1.tolist() != g2.tolist():
                return False
        return True


class source_sink:
    """
    Class for handling all source/sink inputs:

    source_sink.in,
    vsource.th,
    vsink.th,
    msource.th

    These files are always used together,
    so a class is created to bind and process them,
    for example, adding or comparing two source_sink objects

    In addition, there are no complete time info in the *.th files,
    so the TimeHistory class is also used to bind the time info and
    provide additional functions.
    """
    @property
    def source_eles(self):  # element index starts from 1
        return self.source_sink_in.ele_groups[0]

    @property
    def sink_eles(self):  # element index starts from 1
        return self.source_sink_in.ele_groups[1]

    @property
    def nsource(self):
        return self.source_sink_in.n_source

    @property
    def nsink(self):
        return self.source_sink_in.n_sink

    def __init__(self, vsource:TimeHistory, vsink:TimeHistory, msource:list):
        """initialize from TimeHistory objects,
        vsource: TimeHistory object for volume source
        vsink: TimeHistory object for volume sink
        msource: list of TimeHistory objects for mass source

        Note that msource is different from SCHISM's native msource.th
        because saving all tracers in one array is not convenient for
        subsequent processing; instead, each tracer is saved in a separate
        TimeHistory object in a list
        """

        # list of main attributes
        self.vsource = vsource
        self.vsink = vsink
        self.msource = msource
        self.ntracers = None
        self.source_sink_in = None

        # if vsource, vsink, and msources are properly set,
        # then ntracers and source_sink_in can be decided without additional inputs
        if vsource is not None:
            self.ntracers = len(msource)
            source_eles = vsource.df.columns.astype(int).values  # index starts from 1
        else:
            self.ntracers = 0
            source_eles = []

        if vsink is not None:
            sink_eles = vsink.df.columns.astype(int).values
        else:
            sink_eles = []

        self.source_sink_in = SourceSinkIn(ele_groups=[source_eles, sink_eles])

        self.sanity_check()

    @classmethod
    def dummy(cls, start_time_str='2000-01-01 00:00:00', timestamps=[0.0, 86400.0*365*100],
                 source_eles=[], sink_eles=[], ntracers=2):
        """create a dummy source_sink object for testing purpose"""

        # initialize a set of source/sink files from scratch
        nt = len(timestamps)
        nsources = len(source_eles)
        nsinks = len(sink_eles)
        vsource = None
        vsink = None
        msource = [None] * ntracers

        if nsources > 0:
            vsource = TimeHistory(start_time_str=start_time_str,
                                  data_array=np.c_[np.array(timestamps), np.zeros([nt, nsources])],
                                  columns=source_eles)
            # dummy temperature, set to -9999, i.e., ambient temperature
            msource[0] = TimeHistory(start_time_str=start_time_str,
                                     data_array=np.c_[np.array(timestamps), -9999*np.ones([nt, nsources])],
                                     columns=source_eles)
            # dummy salinity, set to 0
            msource[1] = TimeHistory(start_time_str=start_time_str,
                                     data_array=np.c_[np.array(timestamps), np.zeros([nt, nsources])],
                                     columns=source_eles)

        if nsinks > 0:
            vsink = TimeHistory(start_time_str=start_time_str,
                                data_array=np.c_[np.array(timestamps), np.zeros([nt, nsinks])],
                                columns=sink_eles)

        return cls(vsource, vsink, msource)

    @classmethod
    def from_files(cls, source_dir, start_time_str='2000-01-01 00:00:00'):
        '''
        Initialize from existing source/sink files under the source_dir.
        Note that these files don't have start time information,
        and you need to provide it.
        '''
        vsource = None
        msource = None
        vsink = None

        # ele_groups are defined in source_sink.in
        source_sink_in = SourceSinkIn.from_file(filename=f"{source_dir}/source_sink.in")

        if source_sink_in.n_source > 0:
            print('reading vsource\n')
            vsource = TimeHistory.from_file(
                f"{source_dir}/vsource.th",
                start_time_str=start_time_str,
                columns=source_sink_in.ele_groups[0],  # source_eles, index starts from 1
            )

            print('reading msource\n')
            msource_total = loadtxt(f"{source_dir}/msource.th")
            msource_t = msource_total[:, 0]
            ntracers = (msource_total.shape[1] - 1) / vsource.n_station
            if (int(ntracers) != ntracers):
                raise ValueError("Number of tracers must be an integer, vsource and msource don't match")
            else:
                ntracers = int(ntracers)

            # Split msource_total into separate TimeHistory objects,
            # one for each tracer. This facilitates subsequent processing.
            msource = [None] * ntracers
            for i in range(ntracers):
                idx1 = i * vsource.n_station + 1
                idx2 = (i + 1) * vsource.n_station + 1
                msource[i] = TimeHistory(
                    data_array=np.c_[msource_t, msource_total[:, idx1:idx2]],
                    start_time_str=start_time_str,
                    columns=source_sink_in.ele_groups[0]
                )

        if source_sink_in.n_sink > 0:
            print('reading vsink\n')
            vsink = TimeHistory.from_file(
                f"{source_dir}/vsink.th",
                start_time_str=start_time_str,
                columns=source_sink_in.ele_groups[1]
            )

        source_sink = cls(vsource, vsink, msource)
        source_sink.sanity_check()
        return source_sink

    def subset_by_time(self, start_time_str, end_time_str):
        '''
        Subset source/sink files by time.
        '''
        time_slice = slice(start_time_str, end_time_str)
        if self.vsource is not None:
            subset_vsource = self.vsource[time_slice]
            subset_msource = [x[time_slice] for x in self.msource]
        else:
            subset_vsource = None
            subset_msource = None

        if self.vsink is not None:
            subset_vsink = self.vsink[time_slice]
        else:
            subset_vsink = None

        return source_sink(subset_vsource, subset_vsink, subset_msource)

    def subset_by_idx(self, source_idx, sink_idx):
        '''
        Subset source/sink files by index.
        '''
        if self.vsource is not None:
            subset_vsource = self.vsource[:, source_idx]
            subset_msource = [x[:, source_idx] for x in self.msource]
        else:
            subset_vsource = None
            subset_msource = None

        if self.vsink is not None:
            subset_vsink = self.vsink[:, sink_idx]
        else:
            subset_vsink = None

        return source_sink(subset_vsource, subset_vsink, subset_msource)

    def subset_by_ele(self, source_eles=[], sink_eles=[]):
        '''subset source/sink files by element ids (index starts from 1)'''
        if self.vsource is not None:
            if source_eles == []:  # no subsetting
                subset_vsource = self.vsource
                subset_msource = self.msource
            else:
                subset_vsource = self.vsource[source_eles]
                subset_msource = [x[source_eles] for x in self.msource]
        else:
            subset_vsource = None
            subset_msource = None

        if self.vsink is not None:
            if sink_eles == []:
                subset_vsink = self.vsink  # no subsetting
            else:
                subset_vsink = self.vsink[sink_eles]
        else:
            subset_vsink = None

        return source_sink(subset_vsource, subset_vsink, subset_msource)


    def clip_by_polygons(self, hgrid, polygons_xy=[]):
        '''
        Select source/sink elements by polygons.
        An hgrid of schism_grid type is required to get element coordinates.
        polygons: a list of 2D np arrays of (x, y) coordinates of each polygon
        '''
        hgrid.compute_ctr()

        # select source and sink elements
        inside_source = np.zeros(self.nsource, dtype=bool)
        inside_sink = np.zeros(self.nsink, dtype=bool)
        for polygon_xy in polygons_xy:
            ele_xy = np.c_[hgrid.xctr[self.source_eles-1], hgrid.yctr[self.source_eles-1]]
            inside_source += inside_polygon(ele_xy, polygon_xy[:, 0], polygon_xy[:, 1]).astype(bool)

            ele_xy = np.c_[hgrid.xctr[self.sink_eles-1], hgrid.yctr[self.sink_eles-1]]
            inside_sink += inside_polygon(ele_xy, polygon_xy[:, 0], polygon_xy[:, 1]).astype(bool)

        inside_ss = self.subset_by_idx(inside_source, inside_sink)
        outside_ss = self.subset_by_idx(~inside_source, ~inside_sink)

        return inside_ss, outside_ss

    def writer(self, output_dir):
        '''
        Write source/sink files to the output_dir.
        '''
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        self.source_sink_in.writer(f"{output_dir}/source_sink.in")
        if self.vsource is not None:
            self.vsource.writer(f"{output_dir}/vsource.th")

            msource_total = self.msource[0].time
            for i in range(self.ntracers):
                msource_total = np.c_[msource_total, self.msource[i].df.values]
            np.savetxt(f"{output_dir}/msource.th", msource_total)

        if self.vsink is not None:
            self.vsink.writer(f"{output_dir}/vsink.th")

        # additional outputs in *.nc format
        self.nc_writer(output_dir=output_dir)

    def diag_writer(self, hgrid, output_dir):
        '''writer for diagnostic files'''
        hgrid.compute_ctr()

        if self.vsource:
            np.savetxt(
                f'{output_dir}/sources.xyz',
                np.c_[hgrid.xctr[self.source_eles-1], hgrid.yctr[self.source_eles-1], self.vsource.df.mean().values]
            )
        if self.vsink:
            np.savetxt(
                f'{output_dir}/sinks.xyz',
                np.c_[hgrid.xctr[self.sink_eles-1], hgrid.yctr[self.sink_eles-1], self.vsink.df.mean().values]
            )

    def nc_writer(self, output_dir=None):
        if output_dir is None:
            raise Exception("output_dir is required.")
        os.makedirs(output_dir, exist_ok=True)

        # create netcdf data
        C=zdata(); C.vars=[]; C.file_format='NETCDF4'

        # create dummy source/sink if they are empty
        if self.nsource == 0:
            dummy_ss = source_sink.dummy(start_time_str=self.vsink.df.index[0], timestamps=self.vsink.time, source_eles=['1'], sink_eles=['1'], ntracers=2)
            vsource = dummy_ss.vsource
            msource = dummy_ss.msource
            vsink = self.vsink
            nsource = 1
            nsink = self.nsink
            ntracers = 2
        else:
            dummy_ss = source_sink.dummy(start_time_str=self.vsource.df.index[0], timestamps=self.vsource.time, source_eles=['1'], sink_eles=['1'], ntracers=self.ntracers)
            vsource = self.vsource
            msource = self.msource
            vsink = dummy_ss.vsink
            nsource = self.nsource
            nsink = 1
            ntracers = self.ntracers

        C.dimname=['nsources', 'nsinks', 'ntracers', 'time_msource','time_vsource','time_vsink','one']
        C.dims=[nsource, nsink, ntracers, msource[0].n_time, vsource.n_time, vsink.n_time, 1]

        C.vars.extend(['source_elem','vsource','msource'])
        vi=zdata(); vi.dimname=('nsources',); vi.val=vsource.df.columns.values.astype(int); C.source_elem=vi
        vi=zdata(); vi.dimname=('time_vsource','nsources'); vi.val=vsource.data; C.vsource=vi
        msource_data = np.stack([x.data for x in msource], axis=1)  # cast into a 3D array of shape (nt, ntracers, nsources)
        vi=zdata(); vi.dimname=('time_msource','ntracers','nsources'); vi.val=msource_data; C.msource=vi

        C.vars.extend(['sink_elem','vsink'])
        vi=zdata(); vi.dimname=('nsinks',); vi.val=vsink.df.columns.values.astype(int); C.sink_elem=vi
        vi=zdata(); vi.dimname=('time_vsink','nsinks',); vi.val=vsink.data; C.vsink=vi

        C.vars.extend(['time_step_vsource','time_step_msource','time_step_vsink'])
        vi=zdata(); vi.dimname=('one',); vi.val=vsource.delta_t; C.time_step_vsource=vi
        vi=zdata(); vi.dimname=('one',); vi.val=msource[0].delta_t; C.time_step_msource=vi
        vi=zdata(); vi.dimname=('one',); vi.val=vsink.delta_t; C.time_step_vsink=vi

        WriteNC(f'{output_dir}/source.nc', C)

    def sanity_check(self):
        # check consistency of source_sink_in and vsource/vsink/msource
        if self.vsource is None and self.vsink is None:
            raise ValueError("vsource and vsink cannot be both None")

        if self.vsource is not None:
            if len(self.msource) != self.ntracers:
                raise Exception('inconsistent number of tracers')
            if self.nsource != self.msource[0].n_station:
                raise Exception('inconsistent number of msource stations')
            if self.nsource != self.vsource.n_station:
                raise Exception('inconsistent number of vsource stations')
            if np.min(self.vsource.df.values, axis=None) < 0:
                raise Exception('vsource must be non-negative')

        if self.vsink is not None:
            if self.nsink != self.vsink.n_station:
                raise Exception('inconsistent number of sink stations')
            if np.max(self.vsink.df.values, axis=None) > 0:
                raise Exception('vsink must be non-positive')

    def __add__(self, other):
        '''
        Add source/sink other to source/sink self,
        retaining self's time stamps.
        '''
        A = self
        B = other

        # sanity check
        if A.nsource == 0 and B.nsource == 0 and A.nsink == 0 and B.nsink == 0:
            raise Exception('both source and sink are empty')

        # most cases are trivial unless both A and B have source
        if A.nsource == 0 and B.nsource == 0:  # neither has source
            vsource = None
            msource = None
        elif A.nsource == 0:  # only B has source
            vsource = B.vsource
            msource = B.msource
        elif B.nsource == 0:  # only A has source
            vsource = A.vsource
            msource = A.msource
        else:  # both have source
            vsource = A.vsource + B.vsource  # using TimeHistory.__add__
            msource = [None] * A.ntracers
            for i in range(A.ntracers):  # also using TimeHistory.__add__, but with weights
                msource[i] = A.msource[i].__add__(B.msource[i], weights=[A.vsource, B.vsource])

        # most cases are trivial unless both A and B have sink
        if A.nsink == 0 and B.nsink == 0:  # neither has sink
            vsink = None
        elif A.nsink == 0:  # only B has sink
            vsink = B.vsink
        elif B.nsink == 0:  # only A has sink
            vsink = A.vsink
        else:  # both have sink
            vsink = A.vsink + B.vsink  # using TimeHistory.__add__

        return type(self)(vsource=vsource, vsink=vsink, msource=msource)

    def __eq__(self, other):
        for att in ['source_sink_in', 'vsource', 'vsink']:
            if getattr(self, att) != getattr(other, att):
                print(f'{att} not equal')
                return False
        for i, [ms_A, ms_B] in enumerate(zip(self.msource, other.msource)):
            if ms_A != ms_B:
                print('msource {i} not equal')
                return False
        return True

# ---------------------------- unit test ----------------------------
class test_add_source_sink(unittest.TestCase):
    def test_add_source_sink(self):
        print('*************** test_add_source_sink ****************')
        # read from prepared sample files
        A = source_sink.from_files('/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Source_sink/Test_data/source_sink_sample1')
        B = source_sink.from_files('/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Source_sink/Test_data/source_sink_sample2')
        AB = source_sink.from_files('/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Source_sink/Test_data/source_sink_sample12')

        A_add_B = A + B

        self.assertEqual(A_add_B, AB)

class test_combine_dataframes(unittest.TestCase):
    def setUp(self):
        # create some random time series for test
        index_A = pd.date_range(start='2023-01-01', end='2023-12-31', freq='D')
        index_B = pd.date_range(start='2023-01-01', end='2023-10-31', freq='3D')  # lower frequency
        self.df_A = pd.DataFrame(np.random.rand(len(index_A), 3), index=index_A, columns=['A', 'B', 'C'])
        self.df_B = pd.DataFrame(np.random.rand(len(index_B), 3), index=index_B, columns=['B', 'C', 'D'])

    def test_combine_dataframes(self):
        print('*************** test combine dataframes ****************')
        df_combined = combine_dataframes(self.df_A, self.df_B)

        print("First few rows of DataFrame A:")
        print(self.df_A.head())
        print("First few rows of DataFrame B:")
        print(self.df_B.head())
        print("First few rows of Combined DataFrame:")
        print(df_combined.head())

        # Check the shape of the resulting dataframe
        self.assertEqual(df_combined.shape[0], self.df_A.shape[0])  # check row number
        self.assertEqual(df_combined.shape[1], 4)  # check column number

        # Check the column names of the resulting dataframe
        self.assertTrue(all(np.isin(['A', 'B', 'C', 'D'], df_combined.columns)))  # check column names

        # Check that the values of the resulting dataframe are correct
        common_columns = self.df_A.columns.intersection(self.df_B.columns)
        for col in common_columns:
            self.assertTrue(all(np.isclose(df_combined[col].loc[self.df_B.index], self.df_A[col].loc[self.df_B.index] + self.df_B[col], atol=1e-5)))

if __name__ == "__main__":

    unittest.main()
