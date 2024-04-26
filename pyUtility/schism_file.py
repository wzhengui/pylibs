#/usr/bin/env python3
from pylib import *

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
           self.data=value

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

    def compute_area(self,fmt=0):
        '''
        compute element area
        fmt=0: return element area
        fmt=1: return element area, 1st triangle area, and 2nd triangle area
        '''
        x1,x2,x3,x4=self.x[self.elnode].T; y1,y2,y3,y4=self.y[self.elnode].T
        #fp=self.elnode[:,-1]<0; x4[fp]=x1[fp]; y4[fp]=y1[fp]
        fp=self.i34==3; x4[fp]=x1[fp]; y4[fp]=y1[fp]
        a1=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2; a2=((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2;  self.area=a1+a2
        return self.area if fmt==0 else [self.area, a1, a2]

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

        #compute gradients
        x1,x2,x3,x4=x[self.elnode].T; y1,y2,y3,y4=y[self.elnode].T; v1,v2,v3,v4=v0[self.elnode].T
        a1, a2=self.compute_area(fmt=1)[1:]
        dpedx=(v1*(y2-y3)+v2*(y3-y1)+v3*(y1-y2))/(2*a1); dpedy=((x3-x2)*v1+(x1-x3)*v2+(x2-x1)*v3)/(2*a1)

        #for quads, average 2 triangles
        fpn=self.i34==4
        if sum(fpn)!=0:
           dpedx[fpn]=(dpedx[fpn]+(v1*(y3-y4)+v3*(y4-y1)+v4*(y1-y3))[fpn]/(2*a2[fpn]))/2
           dpedy[fpn]=(dpedy[fpn]+((x4-x3)*v1+(x1-x4)*v3+(x3-x1)*v4)[fpn]/(2*a2[fpn]))/2

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
        if not hasattr(self,'isdel') or not hasattr(self,'isidenode'): self.compute_side(fmt=1)
        if self.ns>1e6: print('computing grid boundaries')

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
               A=signa(c_[x1,x2,x3],c_[y1,y2,y3]); fp=(A1+A2-A)>0; A[fp]=A1[fp]+A2[fp]
               pacor[fpn]=c_[A1/A,A2/A,1-(A1+A2)/A]
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

        vi=self.dp if value is None else value; npt=len(vi) #get value
        pie,pip,pacor=self.compute_acor(pxy,fmt=fmt)        #get interp coeff

        #interp
        if npt==self.np:
           return (vi[pip]*pacor).sum(axis=1)
        elif npt==self.ne:
           return vi[pie]
        else:
           sys.exit('unknown data size: {}'.format(npt))

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

    def compute_CFL(self,dt=120, u=0):
        '''
        compute CFL number
        dt: time step (second);  u: flow velocity (m/s)
        '''
        if not hasattr(self,'area'): self.compute_area()
        if not hasattr(self,'dpe'): self.compute_ctr()
        self.cfl=0.5*dt*(abs(u)+sqrt(9.81*self.dpe))/sqrt(self.area/pi); fp4=self.i34==4
        self.cfl[fp4]=self.cfl[fp4]*sqrt(2)
        return self.cfl

    def compute_curl(self,u,v):
        '''
        compute curl of vector filed (u, v):  curl=dv/dx-du/dy
        note: the rotation rate of eddy flow field is: omiga=0.5*curl (radian/second)
        '''
        #pre-proc
        if u.shape[0]==self.ne: u=self.interp_elem_to_node(u); v=self.interp_elem_to_node(v)
        x1,x2,x3,x4=self.x[self.elnode].T; y1,y2,y3,y4=self.y[self.elnode].T
        u1,u2,u3,u4=u[self.elnode].T; v1,v2,v3,v4=v[self.elnode].T;  a1, a2=self.compute_area(fmt=1)[1:]

        #compute gradient
        dvdx=(v1*(y2-y3)+v2*(y3-y1)+v3*(y1-y2))/(2*a1)
        dudy=((x3-x2)*u1+(x1-x3)*u2+(x2-x1)*u3)/(2*a1)

        #for quads, average 2 triangles
        fpn=self.i34==4
        if sum(fpn)!=0:
           dvdx[fpn]=(dvdx[fpn]+(v1*(y3-y4)+v3*(y4-y1)+v4*(y1-y3))[fpn]/(2*a2[fpn]))/2
           dudy[fpn]=(dudy[fpn]+((x4-x3)*u1+(x1-x4)*u3+(x3-x1)*u4)[fpn]/(2*a2[fpn]))/2
        return dvdx-dudy

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
            ac3=1-ac1-ac2; fpn=ac3<0; ac2[fpn]=1-ac1[fpn]; ac3[fpn]=0
            pip.extend(ip[fps]); pacor.extend(c_[ac1,ac2,ac3])
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
            dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata; pz=0
            if dlk==0 and btn==3:
               acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]);ac=acs[nonzero(ats=='bp')[0][0]]
               if hasattr(ac,'bp'):
                  if ac.bp.nsta==0: return
                  distp=squeeze(abs((ac.bp.x-bx)+1j*(ac.bp.y-by))); sid=nonzero(distp==distp.min())[0][0]
                  px=ac.bp.x[sid]; py=ac.bp.y[sid]; pz=1
            elif dlk==0 and btn==1:
               px=bx; py=by; pz=1
            elif dlk==0 and btn==2:
               self.hqt.remove(); self.hqp.remove(); delattr(self,'hqt'); delattr(self,'hqp'); gcf().canvas.draw()
               gcf().canvas.mpl_disconnect(self.cidquery)

            #annotate text
            if dlk==0 and (btn in [1,3]) and pz==1:
               pz=self.interp(c_[px,py],value=self.data if hasattr(self,'data') else self.dp)
               if not hasattr(self,'hqp'):
                  self.hqp=plot(px,py,'g^',ms=6,alpha=1)[0]
                  self.hqt=text(px,py,'',color='orangered',fontsize=12,bbox=dict(facecolor='w',alpha=0.75))
               self.hqt.set_x(px); self.hqt.set_y(py); self.hqt.set_text(' {}'.format(pz[0]))
               self.hqp.set_xdata(px); self.hqp.set_ydata(py); gcf().canvas.draw()

        acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
        abp=acs[nonzero(ats=='query')[0][0]] if 'query' in ats else gcf().canvas.toolbar.addAction('query')
        if hasattr(self,'hqt'): delattr(self,'hqt')
        if hasattr(self,'hqp'): delattr(self,'hqp')
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

def interp_schism_3d(gd,vd,pxy,pz,values,pind=None,zind=None,fmt=0):
    '''
    3D interpolation for multiple variables; interplation only for the variables with a dimension of ne or np
        (gd,  vd):    (hgrid,vgrid) from save_schism_grid
        (pxy, pz): (c_[x,y],   z[:,:])  or (hgrid,vgrid) from save_schism_grid
        values: list of np.array, or np.array
        pind: indice of ne/np dimension for each variable; use -1 to skip interpolation of the variables
        zind: indice of nvrt dimension for each variable;  use -1 to skip the vertical dimension
        fmt=0: higher efficency for many pts; fmt=1: higher efficency for fewer pts
    '''
    #get xyz of interpolation pts
    if hasattr(pxy,'dp') and hasattr(pz,'ivcor'):
        pz=-pz.compute_zcor(pxy.dp); pxy=c_[pxy.x,pxy.y]
    elif isinstance(pxy,np.ndarray) and isinstance(pz,np.ndarray):
        if pz.ndim==1: pz=pz[:,None]
    else:
        sys.exit('(pxy,pz) should be either (schism_grid, schism_vgrid) or (np.array, np.array)')
    npt=len(pxy); nz=pz.shape[1]

    #compute zcor
    if not hasattr(vd,'ivcor'): sys.exit('vgrid should a "schism_vgrid" object')
    nvrt=vd.nvrt; zcor=-vd.compute_zcor(gd.dp)

    #compute area coordinate
    pie,pip,pacor=gd.compute_acor(pxy,fmt=fmt); zcor=(zcor[pip]*pacor[...,None]).sum(axis=1)

    #limit pz
    zm=tile(zcor[:,-1],[nz,1]).T; fpz=pz<zm; pz[fpz]=zm[fpz]
    zm=tile(zcor[:,0],[nz,1]).T;  fpz=pz>zm; pz[fpz]=zm[fpz]; zm=None

    #get variables' dimension information
    ilst,nvar,values=[0,1,[values,]] if isinstance(values,np.ndarray) else [1,len(values),values]
    dims=[array(i.shape) for i in values]; ndims=[len(i) for i in dims]
    if pind is None: pind=nan*ones(nvar)
    if zind is None: zind=nan*ones(nvar)
    pind=array(pind).astype('int'); zind=array(zind).astype('int')

    #modify variables dimension, interp from elem to node, interp horizontally,
    tind=[]; w0=None; pvalues=[]
    for i,[value,ds,ndim] in enumerate(zip(values,dims,ndims)):
        sindp=nonzero(ds==gd.np)[0]; sinde=nonzero(ds==gd.ne)[0]; sindk=nonzero(ds==nvrt)[0]

        #get index of ne or np
        if pind[i]!=-1:
            if len(sindp)==1: pind[i]=sindp[0]
            if len(sindp)==0 and len(sinde)==1: pind[i]=sinde[0]
            if len(sindp)==0 and len(sinde)==0: pind[i]=-1
            if len(sindp)>1 and len(sinde)>1:  sys.exit('need "pind" for np/ne in {}th variable: dims={}, np={}, ne={}'.format(i+1,ds,gd.ne,gd.np))

        #get index for nvrt
        if zind[i]!=-1:
            if len(sindk)==1: zind[i]=sindk[0]
            if len(sindk)==0: zind[i]=-1
            if len(sindk)>1: sys.exit('need "zind" for nvrt in {}th variable: dims={}, nvrt={} '.format(i+1,ds,nvrt))

        #transpose the dimensions
        if pind[i]!=-1 and zind[i]!=-1:
            sind=[pind[i],zind[i],*setdiff1d(arange(ndim),[pind[i],zind[i]])]
        elif pind[i]!=-1 and zind[i]==-1:
            sind=[pind[i],*setdiff1d(arange(ndim),[pind[i]])]
        else:
            sind=None
        tind.append(sind)
        if sind is not None: values[i]=value.transpose(sind); dims[i]=dims[i][sind]

        #compute weight function, and interpolate from elem to node
        if pind[i]!=-1 and dims[i][0]==gd.ne:
            if not hasattr(gd,'nne'): gd.compute_nne()
            if w0 is None: w0=gd.ine!=-1; tw0=w0.sum(axis=1)
            w=w0.copy(); tw=tw0.copy()
            for n in arange(ndim-1): w=expand_dims(w,axis=2); tw=expand_dims(tw,axis=1)
            values[i]=(w*values[i][gd.ine]).sum(axis=1)/tw

        #create variables for pxyz
        dsn=dims[i].copy()
        if pind[i]!=-1: dsn[0]=npt
        if zind[i]!=-1: dsn[1]=nz
        pvalue=zeros(dsn)

        #interp horizontal
        if pind[i]!=-1:
            acor=pacor.copy()
            for n in arange(ndim-1): acor=expand_dims(acor,axis=2)
            values[i]=(values[i][pip]*acor).sum(axis=1)
            if zind[i]==-1: pvalue=values[i]
        pvalues.append(pvalue)

    #interp in vertical
    for k, pzi in enumerate(pz.T):
        for n in arange(nvrt):
            z1,z2=zcor[:,min(nvrt-1,n+1)],zcor[:,n]
            if n==(nvrt-1):
                fpz=pzi==z1; ratz=zeros(sum(fpz))
            else:
                fpz=(pzi>z1)*(pzi<=z2); ratz=(pzi[fpz]-z1[fpz])/(z2[fpz]-z1[fpz])
            if sum(fpz)==0: continue

            for i,[value,ds,ndim,ip,iz] in enumerate(zip(values,dims,ndims,pind,zind)):
                if ip==-1 or iz==-1: continue
                v1,v2=value[fpz,min(nvrt-1,n+1)],value[fpz,n]; rat=ratz.copy()
                for m in arange(ndim-2):rat=expand_dims(rat,axis=1)
                pvalues[i][fpz,k]=v1*(1-rat)+v2*rat

    #restore dimension order
    for i,[pvalue,sind] in enumerate(zip(pvalues,tind)):
        if sind is None: continue
        sinds=argsort(sind); pvalues[i]=pvalues[i].transpose(sinds)
    if ilst==0: pvalues=pvalues[0]

    return pvalues

def getglob(dirpath='.',fmt=0):
    '''
    get global information about schism run (ne,ns,np,nvrt,nproc,ntracers,ntrs)
    dirpath: run directory or outputs directory
    fmt=0: default is 0; fmt(!=0) are for eariler schism versions
    '''

    rstr,bdir=srank(0,dirpath=dirpath,fmt=1)
    fname='{}/local_to_global_{}'.format(bdir,rstr) #local_to_global_0000 or local_to_global_000000

    #get info
    S=zdata()
    S.info=array(open(fname,'r').readline().strip().split()).astype('int')
    if fmt==0:
       S.ns,S.ne,S.np,S.nvrt,S.nproc,S.ntracers=S.info[:6]
       S.ntrs=S.info[6:]
    else:
       sys.exit('fmt unknown')
    return S

def srank(rank=0,dirpath='.',fmt=0):
    '''
    return string of schism rank number ('0032', or '000032')
    dirpath: run directory, or RUN*/outputs
    fmt=0: return rank string; fmt=1: return rank string and the location dir.
    '''
    bdir=None;str_rank=''

    #old format with 4 digits
    if os.path.exists('{}/local_to_global_0000'.format(dirpath)): bdir=os.path.abspath(dirpath); str_rank='{:04}'.format(rank)
    if os.path.exists('{}/outputs/local_to_global_0000'.format(dirpath)): bdir=os.path.abspath('{}/outputs/'.format(dirpath)); str_rank='{:04}'.format(rank)

    #new format with 6 digits
    if os.path.exists('{}/local_to_global_000000'.format(dirpath)): bdir=os.path.abspath(dirpath); str_rank='{:06}'.format(rank)
    if os.path.exists('{}/outputs/local_to_global_000000'.format(dirpath)): bdir=os.path.abspath('{}/outputs/'.format(dirpath)); str_rank='{:06}'.format(rank)

    if fmt==0:
       return str_rank
    elif fmt==1:
       return [str_rank,bdir]

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

def zcor_to_schism_grid(zcor,x=None,value=None):
    '''
    convert z-coordinate transect to a schism grid
    Inputs:
        zcor[npt,nvrt]: z-coordinates for each points
        x[npt]: distances of each pionts from starting point
        value[npt,nvrt]: value assicated for each point @(x,y)
    '''

    #get x,y,z
    ys=zcor.T; nvrt,npt=ys.shape
    xs=tile(arange(npt) if x is None else x,[nvrt,1]).astype('float')
    vs=zeros([nvrt,npt]) if value is None else value.T

    #create schism grid
    gd=schism_grid(); ip=arange(ys.size).reshape([nvrt,npt]);  #original point index
    for k in arange(nvrt-1)[::-1]: fp=ys[k]==ys[k+1]; ip[k,fp]=ip[k+1,fp]  #remove repeated points
    elnode=array([ip[:-1,:-1],ip[:-1,1:],ip[1:,1:],ip[1:,:-1]]).reshape([4,(nvrt-1)*(npt-1)]).T #all quads
    y=ys.ravel()[elnode]; fp=(y[:,0]==y[:,3])*(y[:,1]==y[:,2]); elnode=elnode[~fp]; p=elnode #valid quads
    sindp,sindv=unique(elnode,return_inverse=True) #get unique points
    gd.x,gd.y,gd.dp=xs.ravel()[sindp],ys.ravel()[sindp],vs.ravel()[sindp]
    elnode=arange(sindp.size)[sindv].reshape(elnode.shape) #renumber node index
    p=elnode; p[p[:,3]==p[:,0],3]=-2; p[p[:,2]==p[:,1],2]=-2; fp=p[:,2]==-2; p[fp]=p[fp][:,array([3,0,1,2])] #for triangle
    gd.elnode=elnode; gd.np,gd.ne=len(gd.x),len(gd.elnode); gd.i34=sum(gd.elnode!=-2,axis=1)
    gd.compute_area(); fp3=(gd.i34==3)*(gd.area<0); fp4=(gd.i34==4)*(gd.area<0)
    gd.elnode[fp3,:3]=gd.elnode[fp3,2::-1]; gd.elnode[fp4,:]=gd.elnode[fp4,::-1]; gd.compute_area()

    return gd

def compute_schism_volume(hgrid,vgrid,fmt=0,value=None,eta=None):
    '''
    compute volume or volume*value for schism grid; see help from gd.compute_volume
    '''
    return hgrid.compute_volume(vgrid,fmt,value,eta)

def scatter_to_schism_grid(xyz,angle_min=None,area_max=None,side_min=None,side_max=None,reg_in=None,reg_out=None):
    '''
    create schism grid from scatter pts
        xyz: c_[x,y] or c_[x,y,z]
        angle_min: remove element with internal_angle < angle_min
        area_max:  remove element with element_area   > area_max
        side_min:  remove element with side_length    < side_min
        side_max:  remove element with side_length    > side_max
        reg_in:    ACE/xgredit region file. remove elements in region if reg_in is provided
        reg_out:   ACE/xgredit region file. remove elements outside of region if reg_out is provided
    '''

    #get xyz
    x,y=xyz.T[:2]; np=len(x)
    z=xyz[:,2] if xyz.shape[1]>=3 else zeros(np)

    #triangulate scatter
    gd=schism_grid(); gd.x,gd.y,gd.dp=x,y,z; gd.np=len(gd.x)
    #tr=sp.spatial.Delaunay(c_[x,y]); gd.elnode=tr.simplices #method 1
    tr=mpl.tri.Triangulation(x,y); gd.elnode=unique(tr.triangles,axis=0) #method 2
    gd.ne=len(gd.elnode); gd.elnode=resize(gd.elnode,[gd.ne,4],-2); gd.i34=tile(3,gd.ne)

    #clean mesh
    gd=delete_schism_grid_element(gd,angle_min=angle_min,area_max=area_max,side_min=side_min,side_max=side_max,reg_in=reg_in,reg_out=reg_out)
    return gd

def delete_schism_grid_element(gd,angle_min=5,area_max=None,side_min=None,side_max=None,reg_in=None,reg_out=None,method=0):
    '''
    delete schism grid's elements
        grd: schism_grid object
        angle_min: remove element with internal_angle < angle_min
        area_max:  remove element with element_area   > area_max
        side_min:  remove element with side_length    < side_min
        side_max:  remove element with side_length    > side_max
        reg_in:    ACE/xgredit region file. remove elements in region if reg_in is provided
        reg_out:   ACE/xgredit region file. remove elements outside of region if reg_out is provided
        method=0: use side_max for dangling pts; method=1: use angle_min for dangling pts
    '''

    fpe=tile(True,gd.ne)
    if angle_min is not None: #check min angle
        if not hasattr(gd,'angles'): gd.compute_angle()
        fpe=fpe*(gd.angles.min(axis=1)>=angle_min)

    if area_max is not None: #check maximum area
        if not hasattr(gd,'area'): gd.compute_area()
        fpe=fpe*(gd.area<=area_max)

    if (side_min is not None) or (side_max is not None): #check side length
        slen=[]; ie=arange(gd.ne)
        for n in arange(4):
            i1=gd.elnode[ie,n%gd.i34]; i2=gd.elnode[ie,(n+1)%gd.i34]
            slen.append(abs((gd.x[i2]-gd.x[i1])+1j*(gd.y[i2]-gd.y[i1])))
        slen=array(slen).T
        if side_min is not None: fpe=fpe*(slen.min(axis=1)>=side_min)
        if side_max is not None: fpe=fpe*(slen.max(axis=1)<=side_max)

    for m in arange(2): #check region inside and outside
        regs=reg_in if m==0 else reg_out
        if regs is not None:
            if not hasattr(gd,'xctr'): gd.compute_ctr()
            if isinstance(regs,str) or array(regs[0]).ndim==1: regs=[regs]
            for n,rxy in enumerate(regs):
                if isinstance(rxy,str): bp=read(rxy); rxy=c_[bp.x,bp.y]
                fpe=fpe*(inside_polygon(c_[gd.xctr,gd.yctr], rxy[:,0],rxy[:,1])==(0 if m==0 else 1))

    #add back one element if nodal is zero
    ip=unique(gd.elnode[fpe])[1:]
    if ip.size<gd.np:
       if not hasattr(gd, 'ine'): gd.compute_nne()
       if not hasattr(gd,'angles'): gd.compute_angle()
       sindp=setdiff1d(arange(gd.np),ip); sinde=gd.ine[sindp]; fpn=sinde==-1
       A=gd.angles[sinde]; dA=A.max(axis=2)-A.min(axis=2); dA[fpn]=1e3
       ie=sinde[arange(len(sinde)),dA.argmin(axis=1)]; fpe[ie]=True

    #delete elements and update grid
    gd.ne=sum(fpe); gd.i34=gd.i34[fpe]; gd.elnode=gd.elnode[fpe]
    if hasattr(gd,'angles'): gd.angles=gd.angles[fpe]
    if hasattr(gd,'area'): gd.angles=gd.area[fpe]
    if hasattr(gd,'xctr'): gd.xctr=gd.xctr[fpe]; gd.yctr=gd.yctr[fpe]; gd.dpe=gd.dpe[fpe]
    gd.compute_nne()
    return gd

def combine_icm_output(rundir='.',sname='icm.nc',fmt=0,outfmt=0):
    '''
    combine schism-icm station outputs
      rundir: run directory
      sname: save name for combined output
      fmt=1: copy subdomain output file, and then read
      outfmt=0: float32;  outfmt=1: float64
    '''
    from time import time as gettime
    from shutil import copyfile
    from netCDF4 import Dataset

    outdir=rundir+'/outputs/'
    bp=read_schism_bpfile(rundir+'/istation.in')
    fid=Dataset(outdir+sname,'w',format='NETCDF4'); fvar=fid.variables; fdim=fid.dimensions

    #get fnames information
    fnames=glob(outdir+'icm_*.nc'); mts=[]
    if fmt==1: [copyfile(i,i+'.copy') for i in fnames]; fnames=[i+'.copy' for i in fnames] #copy output files
    for i in fnames: C=ReadNC(i,1); mts.append(array(C.variables['time'])); C.close() #time length
    t1=min([i.min() for i in mts]); t2=max([i[abs(i)<1e15].max() for i in mts]); dt=min([i[1]-i[0] for i in mts])
    mti=arange(t1,t2+dt,dt); nt=mti.size #; sindt=[intersect1d(i,mti,return_indices=True)[1:] for i in mts]

    #combine station output
    for n,fname in enumerate(fnames):
        C=ReadNC(fname,1); cvar=C.variables; cdim=C.dimensions; t0=gettime()
        sind=array(cvar['istation'][:])-1; itc,itf=intersect1d(mts[n],mti,return_indices=True)[1:]
        if n==0:
           #def dim
           for dn in cdim:
               if dn=='nstation':
                  fid.createDimension(dn,bp.nsta)
               elif dn=='time':
                  fid.createDimension(dn,nt)
               else:
                  fid.createDimension(dn,cdim[dn].size)

           #def variables
           for i,cn in enumerate(cvar):
               cdn=[*cvar[cn].dimensions]
               cdn=(cdn[1:] if cdn[0]=='dim_01' else [cdn[1],cdn[0],cdn[2]]) if len(cdn)==3 else cdn
               if len(cdn)>=2 and outfmt==0:
                  fid.createVariable(cn,np.float32,cdn,fill_value=-99999)
               else:
                  fid.createVariable(cn,cvar[cn].dtype,cdn,fill_value=-99999)
           fid.createVariable('station',str,['nstation'],fill_value=False)
           fvar['time'][:]=array(mti)/86400 #set time

        #set variables
        for i,cn in enumerate(cvar):
            cdn=[*cvar[cn].dimensions]; cds=[cdim[k].size for k in cdn]
            if cdn[0]=='nstation': fvar[cn][sind]=array(cvar[cn][:])
            if len(cds)==3 and cds[0]==1: fvar[cn][sind,itf]=array(cvar[cn][0])[:,itc]
            if len(cds)==3 and cds[0]!=1: fvar[cn][sind,:,itf]=array(cvar[cn])[...,itc].transpose([1,0,2])

        fvar['station'][:]=bp.station
        C.close(); dt=gettime()-t0
        if dt>1: print('finish combining {}/{}: {}, {:0.1f}'.format(n+1,len(fnames),fname,dt))
        if fmt==1: os.remove(fname)
    fid.close()

def combine_schism_hotstart(outdir='.',fmt=0,irec=None):
    '''
    combine schism hotstart
       outdir: schism output directory
       fmt=0: combine the last hotstart; fmt=1: combine all hotstart, skip combined; fmt=2: re-combine all
       irec: step number of hotstrat (fmt is disabled when irec is set)
    '''
    from netCDF4 import Dataset

    #get all hotstart records
    irecs=sort([int(i[:-3].split('_')[-1]) for i in os.listdir(outdir) if i.startswith('hotstart_000000_')])
    if fmt==0: irecs=irecs[-1:]
    if irec is not None: irecs=[irec]

    #get subdomain information
    S=get_schism_output_info(outdir,fmt=2); dname=['nResident_node','nResident_elem','nResident_side']
    dnn=['node','elem','side']; dsn=[S.npg,S.neg,S.nsg] #global dim

    fnames=[]
    for m,irec in enumerate(irecs):
        fname='hotstart.nc_{}'.format(irec); fnames.append(fname)
        if fmt==1 and os.path.exists(outdir+'/'+fname): continue
        fid=Dataset(outdir+'/'+fname,'w',format='NETCDF4');  #open file
        print('combining '+fname)
        for n in arange(S.nproc):  #assemble subdomain information
            C=ReadNC('{}/hotstart_{:06}_{}.nc'.format(outdir,n,irec),1)
            cvar=C.variables; cdim=C.dimensions

            #def dimension and variables
            if n==0:
               dn=[*cdim]; ds=[cdim[i].size for i in dn]                    #dim. name and size
               for dni in dname: k=dn.index(dni); dn[k],ds[k]=dnn[k],dsn[k] #change dimension
               for dni,dsi in zip(dn,ds): fid.createDimension(dni,dsi)      #def dims
               sids=tile(-1,len(cvar)) #save index for node/elem/side
               for i,cn in enumerate(cvar): #def variables
                   cdn=[*cvar[cn].dimensions]
                   if cdn[0] in dname: k=dname.index(cdn[0]); cdn[0]=dnn[k]; sids[i]=k
                   fid.createVariable(cn,cvar[cn].dtype,cdn,fill_value=False)

            #set variables
            for i,cn in enumerate(cvar):
               if sids[i]==-1:
                  if n==0: fid.variables[cn][:]=array(cvar[cn][:])
               else:
                  sind=S.ips[n] if sids[i]==0 else (S.ies[n] if sids[i]==1 else S.iss[n])
                  fid.variables[cn][sind]=cvar[cn][:]
            C.close()
        fid.renameVariable('it','iths'); fid.close()
    return fnames

def get_schism_grid_subdomain(grd,xy):
   '''
   compute indice for sub-domian
     grd: SCHISM grid ("grid.npz", "hgrid.gr3", or schism_grid object)
     xy:  subdomin region (c_[x,y], or reg file)
   '''
   #read grid, and region
   gd=(loadz(grd).hgrid if grd.endswith('npz') else read_schism_hgrid(grd)) if isinstance(grd,str) else grd
   gd.compute_ctr(); gd.compute_side(1)
   if isinstance(xy,str): bp=read_schism_reg(xy); xy=c_[bp.x,bp.y]

   #get subdomin information for node,elem. and elnode
   ie=nonzero(inside_polygon(c_[gd.xctr,gd.yctr],xy[:,0],xy[:,1])==1)[0]
   el0=gd.elnode[ie]; ip=unique(el0.ravel())[1:]; ne,np=len(ie),len(ip)
   pid=zeros(gd.np); pid[ip]=arange(np); el=pid[el0].astype('int'); el[el0==-2]=-2

   #create new grid for sub-domain, and compute index for side
   gn=schism_grid(); gn.ne,gn.np,gn.x,gn.y,gn.dp,gn.elnode,gn.i34=ne,np,gd.x[ip],gd.y[ip],gd.dp[ip],el,gd.i34[ie]
   gn.compute_ctr(); gn.compute_side(1); sn=sort(ip[gn.isidenode],axis=1); s0=sort(gd.isidenode,axis=1)
   iA,iB=intersect1d(sn[:,0]+1j*sn[:,1],s0[:,0]+1j*s0[:,1],return_indices=True)[1:]; isd=iB[argsort(iA)]
   gn.sindp,gn.sinde,gn.sinds=ip,ie,isd #save indices of node,elem. and side for subset
   return gn

def get_schism_output_subset(fname,sname,xy=None,grd=None):
   '''
   compute subset of SCHIMS outputs
     fname: original schism outputs (*.nc)
     sname: subset of schism outputs (*.nc)
     xy:    subdomin region (c_[x,y], or reg file)
     grd:   schism grid. a): old grid with xy; b): results from get_schism_grid_subdomain(grd,xy)
   '''
   from netCDF4 import Dataset

   #get subset index
   if (grd is None) and (xy is None): sys.exit('both grd and xy are None')
   if grd is None:
      C=ReadNC(fname,1); cvar=C.variables; gd=schism_grid()
      gd.x=array(cvar['SCHISM_hgrid_node_x'][:]); gd.y=array(cvar['SCHISM_hgrid_node_y'][:])
      gd.dp=array(cvar['depth'][:]); gd.elnode=array(cvar['SCHISM_hgrid_face_nodes'][:])-1; C.close()
      gd.np,gd.ne=len(gd.dp),len(gd.elnode); gd.i34=4*ones(gd.ne).astype('int'); gd.i34[gd.elnode[:,-1]<0]=3
   else:
      gd=(loadz(grd).hgrid if grd.endswith('npz') else read_schism_hgrid(grd)) if isinstance(grd,str) else grd
   if (not hasattr(gd,'sindp')) or (not hasattr(gd,'sinde')) or (not hasattr(gd,'sinds')): gd=get_schism_grid_subdomain(gd,xy)

   #get subset data
   sd={'nSCHISM_hgrid_node':gd.np,'nSCHISM_hgrid_face':gd.ne,'nSCHISM_hgrid_edge':gd.ns}
   sv={'nSCHISM_hgrid_node':gd.sindp,'nSCHISM_hgrid_face':gd.sinde,'nSCHISM_hgrid_edge':gd.sinds}
   def _subset(dname,ds,fmt=0):
       #fmt=0: change dim;  fmt=1: get output subset
       if dname in sd: ds=sd[dname] if fmt==0 else ds[sv[dname]]
       return ds

   #create subset file
   C=ReadNC(fname,1); cdict=C.__dict__; cdim=C.dimensions; cvar=C.variables
   fid=Dataset(sname,'w',format=C.file_format)     #open file
   fid.setncattr('file_format',C.file_format)      #set file_format
   for i in C.ncattrs(): fid.setncattr(i,cdict[i]) #set attrs
   for i in cdim: fid.createDimension(i,None) if cdim[i].isunlimited() else \
                  fid.createDimension(i,_subset(i,cdim[i].size)) #set dims
   for i in cvar: #set vars
       vd=cvar[i].dimensions #;  print(i,vd)
       vid=fid.createVariable(i,cvar[i].dtype,vd,fill_value=True)
       [vid.setncattr(k,cvar[i].getncattr(k)) for k in cvar[i].ncattrs() if (k not in ['_FillValue'])]
       if i=='SCHISM_hgrid_face_nodes':
          fid.variables[i][:]=gd.elnode+1
       elif i=='SCHISM_hgrid_edge_nodes':
          fid.variables[i][:]=gd.isidenode+1
       elif i=='time':
          fid.variables[i][:]=cvar[i][:]
       else:
          if vd[0]=='time' and len(vd)>1:
             for n,k in enumerate(cvar[i][:]): fid.variables[i][n]=_subset(vd[1],k,1)
          else:
             fid.variables[i][:]=_subset(vd[0],cvar[i][:],1)
   fid.close(); C.close()
   return gd

def read_schism_output(run,varname,xyz,stacks=None,ifs=0,nspool=1,sname=None,fname=None,hgrid=None,vgrid=None,fmt=0,mdt=None,extend=0,prj=None,zcor=0,sgrid=None):
    '''
    extract time series of SCHISM results @xyz or transects @xy (works for scribe IO and combined oldIO)
       run:     run directory where (grid.npz or hgrid.gr3) and outputs are located
       varname: variables to be extracted; accept shortname(s) or fullname(s) (elev, hvel, horizontalVelX, NO3, ICM_NO3, etc. )
       xyz:     c_[x,y,z], or bpfile, or c_[x,y]
       fmt:     (optional) 0: read time series @xyz;     1: read transect @xy
       stacks:  (optional) output stacks to be extract; all avaiable stacks will be extracted if not specified
       ifs=0:   (optional) extract results @xyz refers to free surface (default); ifs=1: refer to fixed levels
       nspool:  (optional) sub-sampling frequency within each stack (npsool=1 means all records)
       mdt:     (optional) time window (day) for averaging output
       sname:   (optional) variable name for save
       fname:   (optional) save the results as fname.npz
       hgrid:   (optional) hgrid=read_schism_hgrid('hgrid.gr3'); hgrid.compute_all(); used to speed up
       vgrid:   (optional) vgrid=read_schism_vgrid('vgrid.in'); used to speed up
       extend:  (optional) 0: extend bottom value beyond;  1: assign nan for value beyond bottom
       prj:     (optional) used to tranform xy (e.g. prj=['epsg:26918','epsg:4326'])
       zcor:    (optional) 0: zcoordinate computed from elev;  1: read from outputs
       sgrid:   (optional) side-based grid used to extract side values (FEM method); see gd.scatter_to_grid(fmt=2)
    '''

    #get schism outputs information
    bdir=run+'/outputs'; modules,outfmt,dstacks,dvars,dvars_2d =get_schism_output_info(bdir,1)

    #variables
    if isinstance(varname,str): varname=[varname]
    if isinstance(sname,str): sname=[sname]

    #read grid
    fexist=os.path.exists; fgz=run+'/grid.npz'; fgd=run+'/hgrid.gr3'; fvd=run+'/vgrid.in'
    gd=hgrid if (hgrid is not None) else loadz(fgz,'hgrid') if fexist(fgz) else read_schism_hgrid(fgd)
    if zcor==0: vd=vgrid if (vgrid is not None) else loadz(fgz,'vgrid') if fexist(fgz) else read_schism_vgrid(fvd)

    #read station coordinates (xyz)
    if isinstance(xyz,str): bp=read_schism_bpfile(xyz); xyz=c_[bp.x,bp.y,bp.z]
    if array(xyz).ndim==1: xyz=array(xyz)[None,:]
    lx,ly=xyz.T[:2]; npt=len(lx); lz=xyz.T[2] if xyz.shape[1]==3 else zeros(npt)
    if prj is not None: lx,ly=proj_pts(lx,ly,prj[0],prj[1])
    pie,pip,pacor=gd.compute_acor(c_[lx,ly]); pip,pacor=pip.T,pacor.T; P=zdata()
    def _sindex(sgrid,P): #for output@side
        if sgrid is None:
           if not hasattr(gd,'nns'): gd.compute_side(fmt=2)
           if not hasattr(P,'nns'): P.nns,P.ins=gd.nns[pip],gd.ins[pip]; P.ds=P.ins.shape; P.fp=P.ins!=0
        else:
           if isinstance(sgrid,int): sgrid=gd.scatter_to_grid(fmt=2)
           if not hasattr(P,'pie'): P.pie,P.pip,P.pacor=sgrid.compute_acor(c_[lx,ly])

    #extract time series@xyz or transect@xy
    mtime=[]; mdata=[[] for i in varname]
    stacks=dstacks if (stacks is None) else [*array(stacks).ravel()] #check outputs stacks
    for istack in stacks:
        print('reading stack: {}'.format(istack))
        C0=ReadNC('{}/{}_{}.nc'.format(bdir,'out2d' if outfmt==0 else 'schout',istack),1)
        if fmt==0 and zcor==1: #open zcoordinates channel
            if outfmt==0:
              Z0=ReadNC('{}/zCoordinates_{}.nc'.format(bdir,istack),1); Z=Z0.variables['zCoordinates']
            else:
              Z=C0.variables['zcor']
        mti0=array(C0.variables['time'])/86400; nt=len(mti0); mti=mti0[::nspool]; nrec=len(mti); mtime.extend(mti); ntime=len(mti); zii=None
        if istack==stacks[0]: nvrt=C0.dimensions['nSCHISM_vgrid_layers'].size

        for m,varnamei in enumerate(varname):
            vs=[]; svars=get_schism_var_info(varnamei,modules,fmt=outfmt) #get variable information
            for n,[vari,svar] in enumerate(svars):
                if svar in dvars_2d: #2D
                    C=C0.variables[svar]; np=C.shape[1]
                    if np==gd.ns: _sindex(sgrid,P)
                    if np==gd.np: vi=array([(array([C[i,j] for j in pip])*pacor).sum(axis=0) for i in arange(0,nt,nspool)])
                    if np==gd.ne: vi=array([C[i,pie] for i in arange(0,nt,nspool)])
                    if (np==gd.ns) and (sgrid is None): vi=array([(sum(array(C[i,P.ins.ravel()]).reshape(P.ds)*P.fp,axis=2)*pacor/P.nns).sum(axis=0) for i in arange(0,nt,nspool)])
                    if (np==gd.ns) and (sgrid is not None): vi=array([sum(array(C[i][P.pip]*P.pacor),axis=1) for i in arange(0,nt,nspool)])
                    vs.append(vi)
                else: #3D
                    #read zcoor,and extend zcoor downward
                    if (zii is None) and fmt==0:
                       if n==0:
                          if zcor==0:
                             eta=array([C0['elevation'][i][pip] for i in arange(0,nt,nspool)]).ravel(); sindp=tile(pip,[ntime,1,1]).ravel()
                             sigma=vd.sigma[sindp] if vd.ivcor==1 else 0; zii0=compute_zcor(sigma,gd.dp[sindp],eta,ivcor=vd.ivcor,vd=vd)
                             zii=(zii0.reshape([ntime,3,npt,nvrt])*pacor[None,...,None]).sum(axis=1).transpose([2,0,1])
                          else:
                             zii=(array([array(Z[i])[pip] for i in arange(0,nt,nspool)])*pacor[None,...,None]).sum(axis=1).transpose([2,0,1])
                       for k in arange(nvrt-1): z1=zii[nvrt-k-2]; z2=zii[nvrt-k-1]; z1[abs(z1)>1e8]=z2[abs(z1)>1e8]
                       if ifs==0: zii=zii-zii[-1][None,...]

                    #read data for the whole vertical
                    C1=ReadNC('{}/{}_{}.nc'.format(bdir,svar,istack),1) if outfmt==0 else C0
                    C=C1.variables[svar]; np=C.shape[1]; nd=C.ndim
                    if np==gd.ns: _sindex(sgrid,P)
                    if np==gd.np and nd==3: vii=(array([array(C[i])[pip] for i in arange(0,nt,nspool)])*pacor[None,...,None]).sum(axis=1)
                    if np==gd.np and nd==4: vii=(array([array(C[i])[pip] for i in arange(0,nt,nspool)])*pacor[None,...,None,None]).sum(axis=1)
                    if np==gd.ne: vii=array([C[i][pie] for i in arange(0,nt,nspool)])
                    if (np==gd.ns) and (sgrid is None): vii=array([(sum(C[i][P.ins]*P.fp[...,None],axis=2)*pacor[...,None]/P.nns[...,None]).sum(axis=0) for i in arange(0,nt,nspool)])
                    if (np==gd.ns) and (sgrid is not None): vii=array([sum(C[i][P.pip]*P.pacor[...,None],axis=1) for i in arange(0,nt,nspool)])
                    vii=vii.transpose([2,0,1]) if nd==3 else vii.tranpose([2,0,1,3]) #from (nt,npt,nvrt) to (nvrt,nt,npt)
                    if extend==0:
                       for k in arange(nvrt-1)[::-1]: z1=vii[k]; z2=vii[k+1]; fpn=abs(z1)>1e8; z1[fpn]=z2[fpn] #extend value at bottom
                    else:
                       fpn=abs(vii)>1e8; vii[fpn]=nan
                    if outfmt==0: C1.close()

                    #interp in the vertical
                    if fmt==0: #time series
                       vii=interp_vertical(vii,zii,-lz[None,None,:]); vs.append(vii)
                       if nd==4: vs=vs[0].transpose([2,0,1])
                    else: #transect
                       if nd==3: vs.append(vii)
                       if nd==4: vs.append(vii[...,0]); vs.append(vii[...,1])
            mdata[m].extend(array(vs).transpose([1,2,0])) if (fmt==0 or (svar in dvars_2d))*(array(vs).ndim==3) else mdata[m].extend(array(vs).transpose([2,3,1,0]))
        C0.close()
        if fmt==0 and outfmt==0 and zcor==1: Z0.close()

    #save data
    S=zdata(); sdict=S.__dict__; mt=array(mtime); mdata=[array(i) for i in mdata]
    for m,k in enumerate(varname if sname is None else sname):
        vi=mdata[m] if (mdt is None) else array([mdata[m][(mt>=i)*(mt<(i+mdt))].mean(axis=0) for i in arange(mt[0],mt[-1],mdt)]) #average
        sdict[k]=squeeze(vi.transpose([1,0,2] if vi.ndim==3 else [1,0,2,3]))
    S.time=mt if (mdt is None) else array([mt[(mt>=i)*(mt<(i+mdt))].mean() for i in arange(mt[0],mt[-1],mdt)]) #average
    if fname is not None: S.save(fname)
    return S

def read_schism_slab(run,varname,levels,stacks=None,nspool=1,mdt=None,sname=None,fname=None,reg=None):
    '''
    extract slabs of SCHISM results (works for scribe IO and combined oldIO)
       run:     run directory where (grid.npz or hgrid.gr3,vgrid.in) and outputs are located
       varname: variables to be extracted; accept shortname(s) or fullname(s) (elev, hvel, horizontalVelX, NO3, ICM_NO3, etc. )
       levels:  schism level indices (1-nvrt: surface-bottom; (>nvrt): kbp level; "all": all layers)
       stacks:  (optional) output stacks to be extract; all avaiable stacks will be extracted if not specified
       nspool:  (optional) sub-sampling frequency within each stack (npsool=1 means all records)
       mdt:     (optional) time window (day) for averaging output
       sname:   (optional) variable name for save
       fname:   (optional) save the results as fname.npz
       reg:     (optional) subsetting reslts inside region (*.reg, or *.bp, or gd_subgrid)
    '''
    #proc
    fgz=run+'/grid.npz'; fgd=run+'/hgrid.gr3'; fvd=run+'/vgrid.in'; fexist=os.path.exists; P=zdata()
    bdir=run+'/outputs'; modules,outfmt,dstacks,dvars,dvars_2d =get_schism_output_info(bdir,1)
    gd=read(fgz,'hgrid') if fexist(fgz) else read(fgd); np,ne,ns=gd.np,gd.ne,gd.ns
    if isinstance(varname,str): varname=[varname]
    if sname is None: sname=varname
    if stacks is None: stacks=dstacks
    if outfmt==1: sys.exit('OLDIO not supported yet')
    if reg is not None: #build subgrid
       sgd=gd.subset(reg) if isinstance(reg,str) else reg; sindp,sinde,sinds=sgd.sindp,sgd.sinde,sgd.sinds

    #read output
    mtime=[]; mdata=[[] for i in varname]
    for istack in [*unique(stacks).ravel()]:
        if outfmt==0:
           C0=ReadNC('{}/out2d_{}.nc'.format(bdir,istack),1); nvrt=C0.dimensions['nSCHISM_vgrid_layers'].size
           mt=array(C0.variables['time'][:])/86400; nt=len(mt); mtime.extend(mt[::nspool])
           zs=arange(nvrt,0,-1) if levels=='all' else array(levels)

        for ivar, varnamei in enumerate(varname):
            svars=get_schism_var_info(varnamei,modules,fmt=outfmt)
            for m,[vari,svar] in enumerate(svars):
                C=C0 if (svar in dvars_2d) else read('{}/{}_{}.nc'.format(bdir,svar,istack),1); cvar=C.variables[svar]; vi=[]
                if svar in dvars_2d:  #2D
                   vi=array([array(cvar[i]).astype('float32') for i in arange(nt) if i%nspool==0])
                   if reg is not None: npt=cvar.shape[1]; vi=vi[:,sindp if npt==np else sinde if npt==ne else sinds] #subset
                else:   #3D
                   for n,k in enumerate(zs):
                       if k>nvrt:
                          if not hasattr(P,'kbp'): vd=read(fgz,'vgrid') if fexist(fgz) else read(fvd); P.sindp,P.kbp=arange(np),vd.kbp
                          vii=array([array(cvar[i][P.sindp,P.kbp]).astype('float32') for i in arange(nt) if i%nspool==0])
                       else:
                          vii=array([array(cvar[i,:,nvrt-k]).astype('float32') for i in arange(nt) if i%nspool==0])
                       if reg is not None: npt=cvar.shape[1]; vii=vii[:,sindp if npt==np else sinde if npt==ne else sinds] #subset
                       vi.append(vii)
                   vi=vi[0] if len(zs)==1 else array(vi).transpose([1,0,2]); C.close()
                vs=vi if m==0 else c_[vs[...,None],vi[...,None]]
            mdata[ivar].extend(vs)
        C0.close()

    #save data
    S=zdata(); sdict=S.__dict__; S.zs=array(zs); mt=array(mtime); mdata=[array(i) for i in mdata]
    for m, k in enumerate(sname):
        sdict[k]=mdata[m] if (mdt is None) else array([mdata[m][(mt>=i)*(mt<(i+mdt))].mean(axis=0) for i in arange(mt[0],mt[-1],mdt)])
    sdict['time']=mt if (mdt is None) else array([mt[(mt>=i)*(mt<(i+mdt))].mean() for i in arange(mt[0],mt[-1],mdt)])
    if reg is not None: S.gd=sgd
    if fname is not None: S.save(fname)
    return S

def get_schism_output_info(run,fmt=0):
    '''
    get info. about SCHISM ["modules", "output format", "stacks", "output variables", "output 2D variables"]
    Inputs:
        run:  SCHISM outputs directory
        fmt=0: remove schsim default variables related to grid; fmt=1: not remove

    Outputs:
        fmt=0/1: [modules,outfmt,stacks,svars,svars_2d]
        fmt=2: information about domain-decomposition
        fmt=3: gather run information
        fmt=4: return reconstructed SCHISM hgrid
        (Note: fmt=3/4 works for new schism outputs)
    '''

    if fmt in [0,1]:
       #default variables to be excluded
       dvars=('out2d','schout','time','minimum_depth','SCHISM_hgrid','crs','SCHISM_hgrid_node_x','SCHISM_hgrid_node_y',
              'depth','bottom_index_node','SCHISM_hgrid_face_x','SCHISM_hgrid_face_y','SCHISM_hgrid_edge_x',
              'SCHISM_hgrid_edge_y','SCHISM_hgrid_face_nodes','SCHISM_hgrid_edge_nodes','dryFlagNode','Cs',
              'coordinate_system_flag','dry_value_flag','edge_bottom_index','ele_bottom_index','node_bottom_index',
              'sigma','sigma_h_c','sigma_maxdepth','sigma_theta_b','sigma_theta_f','wetdry_elem','wetdry_node',
              'wetdry_side','dryFlagElement', 'dryFlagSide')
       if fmt==1: dvars=('out2d','schout','time')

       #get SCHISM variables
       ovars=unique([i[:i.rfind('_')] for i in os.listdir(run) if (i.endswith('.nc') \
             and (not i.startswith('hotstart_')) and i.rfind('_')!=-1)])
       outfmt=0 if ('out2d' in ovars) else 1 #output format: newIO or oldIO

       #get additional SCHISM variables
       if outfmt==0:
          stacks=unique([int(i[:-3].split('_')[-1]) for i in glob(run+'/out2d_*.nc')])
          C=ReadNC(run+'/out2d_{}.nc'.format(stacks[0]),1); svars_2d=setdiff1d(array([*C.variables]),dvars); C.close()
          svars=setdiff1d(r_[ovars,svars_2d],dvars)
       else:
          stacks=unique([int(i[:-3].split('_')[-1]) for i in glob(run+'/schout_*.nc')])
          fname=[i for i in os.listdir(run) if (i.startswith('schout_') and i.endswith('_{}.nc'.format(stacks[0])))][0]
          C=ReadNC(run+'/{}'.format(fname),1); svars=setdiff1d(array([*C.variables]),dvars)
          svars_2d=array([i for i in svars if C.variables[i].ndim==2]); C.close()

       #get SCHISM modules
       M=get_schism_var_info(fmt=outfmt).__dict__
       modules=array([k for k in M if len(intersect1d([*M[k].values()],svars))!=0])
       return [modules,outfmt,stacks,svars,svars_2d]
    elif fmt==2:
       nes=[]; nps=[]; nss=[]; ies=[]; ips=[]; iss=[]
       nsg,neg,npg,nvrt,nproc,ntr=[int(i) for i in open(run+'/local_to_global_000000','r').readline().split()[:6]]
       for n in arange(nproc):
           lines=open('{}/local_to_global_{:06}'.format(run,n),'r').readlines()[2:]
           ne=int(lines[0].strip()); np=int(lines[ne+1].strip()); ns=int(lines[ne+np+2].strip())
           iei=array([i.split()[1] for i in lines[1:(ne+1)]]).astype('int')-1
           ipi=array([i.split()[1] for i in lines[(ne+2):(ne+np+2)]]).astype('int')-1
           isi=array([i.split()[1] for i in lines[(ne+np+3):(ne+np+ns+3)]]).astype('int')-1
           nes.append(ne); nps.append(np); nss.append(ns); ies.append(iei); ips.append(ipi); iss.append(isi)
       S=zdata(); S.nsg,S.neg,S.npg,S.nvrt,S.nproc,S.ntr=nsg,neg,npg,nvrt,nproc,ntr
       S.nes,S.nps,S.nss,S.ies,S.ips,S.iss=array(nes),array(nps),array(nss),ies,ips,iss
       return S
    elif fmt==3:
       fns=array(glob(run+'/out2d_*.nc')); ifs=array([int(i.replace('_','.').split('.')[-2]) for i in fns]); idx=argsort(ifs); ifs,fns=ifs[idx],fns[idx] #files
       f=ReadNC(fns[0],1); v=f.variables['time']; StartT=datenum(*[int(i.split('.')[0]) for i in v.base_date.split()]); dt=(v[1]-v[0])/86400; ns=len(v); f.close() #StartT
       mis=array([arange(ns)+(i-1)*ns for i in ifs]).ravel(); mts=dt*mis+StartT; EndT=mts[-1] #time
       irec=mod(mis,ns); istack=array([tile(i,ns) for i in ifs]).ravel(); f=ReadNC(fns[-1],1); nm=len(f.variables['time'])-ns; f.close() #stack and record
       if nm!=0: mts,istack,irec=mts[:nm],istack[:nm],irec[:nm]
       S=zdata(); S.outdir,S.StartT, S.EndT, S.dt, S.nrec, S.ifs, S.fns, S.istack, S.irec, S.mts=run,StartT,EndT,dt,ns,ifs,fns,istack,irec,mts
       return S
    elif fmt==4:
       sinfo=get_schism_output_info(run,3); f=ReadNC(sinfo.fns[0],1); cvar=f.variables; gd=schism_grid()
       gd.x=array(cvar['SCHISM_hgrid_node_x']); gd.y=array(cvar['SCHISM_hgrid_node_y']); gd.dp=array(cvar['depth']); gd.elnode=array(cvar['SCHISM_hgrid_face_nodes'])-1
       gd.np,gd.ne=gd.dp.size,len(gd.elnode); gd.elnode[(gd.elnode<0)|(abs(gd.elnode)>1e8)]=-2; gd.i34=sum(gd.elnode!=-2,axis=1); gd.ns=cvar['SCHISM_hgrid_edge_x'].size; f.close()
       return gd

def get_schism_var_info(svar=None,modules=None,fmt=0):
    '''
      usage:
         svar: variable name (either short and long names)
         modules: modules that svar may belong to
         fmt=0: scribe IO;   fmt=1: OLDIO (schout_*.nc)

         1. get_schism_var_info(): return all module information
         2. varname,fullname=get_schism_var_info('elev')
            or varname,fullname=get_schism_var_info('elev',['Hydro',]): get SCHISM output variable information
    '''

    #dicts for SCHISM output variables
    oHydro={'elev':'elev','apres':'air_pressure','atemp':'air_temperature',
           'shum':'specific_humidity','srad':'solar_radiation','heat_sen':'sensible_flux',
           'heat_lat':'latent_heat','rad_up':'upward_longwave','rad_down':'downward_longwave',
           'heat_tot':'total_heat_flux','evap':'evaporation','prec':'precipitation',
           'tau_bot':'bottom_stress','wind':'wind_speed','tau_wind':'wind_stress','dahv':'dahv',
           'zvel':'vertical_velocity','temp':'temp','salt':'salt','rho':'water_density',
           'diff':'diffusivity','vis':'viscosity','TKE':'TKE','ML':'mixing_length',
           'zcor':'zCoordinates','hvel':'hvel', 'hvel_side':'hvel_side','zvel_elem':'wvel_elem',
           'temp_elem':'temp_elem','salt_elem':'salt_elem','pres':'pressure_gradient'}

    Hydro={'elev':'elevation','apres':'airPressure','atemp':'airTemperature',
           'shum':'specificHumidity','srad':'solarRadiation','heat_sen':'sensibleHeat',
           'heat_lat':'latentHeat','rad_up':'upwardLongwave','rad_down':'downwardLongwave',
           'heat_tot':'totalHeat','evap':'evaporationRate','prec':'precipitationRate',
           'tau_bot_x':'bottomStressX','tau_bot_y':'bottomStressY','wind_x':'windSpeedX','wind_y':'windSpeedY',
           'tau_wind_x':'windStressX','tau_wind_y':'windStressY','dahv_x':'depthAverageVelX','dahv_y':'depthAverageVelY',
           'zvel':'verticalVelocity','temp':'temperature','salt':'salinity','rho':'waterDensity',
           'diff':'diffusivity','vis':'viscosity','TKE':'turbulentKineticEner','ML':'mixingLength',
           'zcor':'zCoordinates','hvel_x':'horizontalVelX','hvel_y':'horizontalVelY',
           'hvel_side_x':'horizontalSideVelX','hvel_side_y':'horizontalSideVelY',
           'zvel_elem':'verticalVelAtElement','temp_elem':'temperatureAtElement','salt_elem':'salinityAtElement',
           'pres_x':'barotropicPresGradX','pres_y':'barotropicPresGradY'}

    ICM={'PB1':'ICM_PB1','PB2':'ICM_PB2','PB3':'ICM_PB3',  #core: phytoplankton
        'RPOC':'ICM_RPOC','LPOC':'ICM_LPOC','DOC':'ICM_DOC', #core: carbon
        'RPON':'ICM_RPON','LPON':'ICM_LPON','DON':'ICM_DON','NH4':'ICM_NH4','NO3':'ICM_NO3', #core: nitrogen
        'RPOP':'ICM_RPOP','LPOP':'ICM_LPOP','DOP':'ICM_DOP','PO4':'ICM_PO4', #core: phosphorus
        'COD':'ICM_COD','DO':'ICM_DOX', #core: COD and DO
        'SU':'ICM_SU','SA':'ICM_SA', #silica
        'ZB1':'ICM_ZB1','ZB2':'ICM_ZB2', #zooplankton
        'TIC':'ICM_TIC','ALK':'ICM_ALK','CA':'ICM_CA','CACO3':'ICM_CACO3', #pH
        'sleaf':'ICM_sleaf','sstem':'ICM_sstem','sroot':'ICM_sroot', #sav: leaf/stem/root
        'stleaf':'ICM_stleaf','ststem':'ICM_ststem','stroot':'ICM_stroot','sht':'ICM_sht', #sav: total leaf/stem/root, height
        'vtleaf1':'ICM_vtleaf1','vtleaf2':'ICM_vtleaf2','vtleaf3':'ICM_vtleaf3', #veg: leaf
        'vtstem1':'ICM_vtstem1','vtstem2':'ICM_vtstem2','vtstem3':'ICM_vtstem3', #veg: stem
        'vtroot1':'ICM_vtroot1','vtroot2':'ICM_vtroot2','vtroot3':'ICM_vtroot3', #veg: root
        'vht1':'ICM_vht1','vht2':'ICM_vht2','vht3':'ICM_vht3',                   #veg: vht
        'bPOC1':'ICM_bPOC1','bPOC2':'ICM_bPOC2','bPOC3':'ICM_bPOC3','bPON1':'ICM_bPON1', #SFM
        'bPON2':'ICM_bPON2','bPON3':'ICM_bPON3','bPOP1':'ICM_bPOP1','bPOP2':'ICM_bPOP2', #SFM
        'bPOP3':'ICM_bPOP3','bNH4':'ICM_bNH4','bNO3':'ICM_bNO3','bPO4':'ICM_bPO4','bH2S':'ICM_bH2S', #SFM
        'bCH4':'ICM_bCH4','bPOS':'ICM_bPOS','bSA':'ICM_bSA','bstc':'ICM_bstc','bSTR':'ICM_bSTR', #SFM
        'bThp':'ICM_bThp','bTox':'ICM_bTox','SOD':'ICM_SOD','JNH4':'ICM_JNH4','JNO3':'ICM_JNO3', #SFM
        'JPO4':'ICM_JPO4','JSA':'ICM_JSA','JCOD':'ICM_JCOD','BA':'ICM_BA'} #SFM

    COS={'NO3':'COS_NO3','SiO4':'COS_SiO4','NH4':'COS_NH4','S1':'COS_S1','S2':'COS_S2','Z1':'COS_Z1',
         'Z2':'COS_Z2','DN':'COS_DN','DSi':'COS_DSi','PO4':'COS_PO4','DO':'COS_DOX','CO2':'COS_CO2'}

    SED={'sed_dp':'sedBedThickness','sed_str':'sedBedStress','sed_rough':'sedBedRoughness',
         'sed_por':'sedPorocity','sed_eflux':'sedErosionalFlux','sed_dflux':'sedDepositionalFlux',
         'sed_frac1':'sedBedFraction_1','sed_frac2':'sedBedFraction_2','sed_frac3':'sedBedFraction_3','sed_frac4':'sedBedFraction_4',
         'sed_conc1':'sedConcentration_1','sed_conc2':'sedConcentration_2','sed_conc3':'sedConcentration_3','sed_conc4':'sedConcentration_4',
         'sed_tconc':'totalSuspendedLoad'}
    WWM={'wwm_WVHT':'sigWaveHeight','wwm_DPD':'dominantDirection','wwm_Tp':'peakPeriod'}

    if fmt==1: Hydro=oHydro
    mdict={'Hydro':Hydro,'ICM':ICM,'COS':COS,'SED':SED,'WWM':WWM}; S=zdata(); C=[]
    if (modules is not None) and isinstance(modules,str): modules=[modules,]
    if svar is None: #get all module info
       for n,m in mdict.items(): S.__dict__[n]=m
       return S
    else: #get variable info
       for n,m in mdict.items():
           if (modules is not None) and (n not in modules): continue
           for key,value in m.items():
               if svar in [key,value]: C.append((key,value))
               if key.endswith('_x') or key.endswith('_y'):
                  if svar in [key[:-2],value[:-2]]: C.append((key,value))
       if len(C)==0: C=[(svar,svar)]
       return C

def convert_schism_source(run='.',fname='source.nc'):
    '''
    convert schism source_sink format from ASCII (source_sink.in, vsource.th, vsink.th, msource.th) to source.nc
    '''
    sdir=os.path.abspath(run)+os.path.sep
    #open data cap
    C=zdata(); C.vars=[]; C.file_format='NETCDF4'
    nsources,nsinks,ntracers,ntm,ntv,nts=0,0,0,0,0,0; dtm,dtv,dts=1e6*ones(3)
    
    #source_sink.in
    lines=[i.strip() for i in open(sdir+'source_sink.in','r').readlines()]
    nsources=int(lines[0].split()[0]); isource=array([int(i.split()[0]) for i in lines[1:(nsources+1)]])
    nsinks=int(lines[2+nsources].split()[0]); isink=array([int(i.split()[0]) for i in lines[(nsources+3):(nsources+3+nsinks)]])
    
    #vsource.th and msource.th
    if nsources>0:
       #read data
       fdata=loadtxt(sdir+'msource.th').T; mti=fdata[0]; msource=fdata[1:]
       fdata=loadtxt(sdir+'vsource.th').T; vti=fdata[0]; vsource=fdata[1:]
       ntm=len(mti); ntv=len(vti); dtm=floor(diff(mti)[0]); dtv=floor(diff(vti)[0]); ntracers=int(len(msource)/len(vsource))
       vsource=vsource.T;  msource=msource.reshape([ntracers,nsources,len(mti)]).transpose([2,0,1])
    
       #save data
       C.vars.extend(['source_elem','vsource','msource'])
       vi=zdata(); vi.dimname=('nsources',); vi.val=isource; C.source_elem=vi
       vi=zdata(); vi.dimname=('time_vsource','nsources'); vi.val=vsource; C.vsource=vi
       vi=zdata(); vi.dimname=('time_msource','ntracers','nsources'); vi.val=msource; C.msource=vi
    
    #vsink.th
    if nsinks>0:
       #read data
       fdata=loadtxt('{}/vsink.th'.format(sdir)).T; sti=fdata[0]; vsink=fdata[1:].T
       nts=len(sti); dts=floor(diff(sti)[0])
    
       #save data
       C.vars.extend(['sink_elem','vsink'])
       vi=zdata(); vi.dimname=('nsinks',); vi.val=isink; C.sink_elem=vi
       vi=zdata(); vi.dimname=('time_vsink','nsinks',); vi.val=vsink; C.vsink=vi
    
    #assign dimension value
    C.dimname=['nsources','nsinks','ntracers','time_msource','time_vsource','time_vsink','one']
    C.dims=[nsources,nsinks,ntracers,ntm,ntv,nts,1]
    
    #add time step 
    C.vars.extend(['time_step_vsource','time_step_msource','time_step_vsink'])
    vi=zdata(); vi.dimname=('one',); vi.val=dtm; C.time_step_msource=vi
    vi=zdata(); vi.dimname=('one',); vi.val=dtv; C.time_step_vsource=vi
    vi=zdata(); vi.dimname=('one',); vi.val=dts; C.time_step_vsink=vi
    
    #save as netcdf
    WriteNC(sdir+fname,C)

class schism_view:
    def __init__(self, run='.'):
        #note: p is a capsule including all information about a figure
        self.figs=[]; self.fns=[]; self._nf=0; self.sbp=[] #list of figure objects
        self.run_info(run)
        self.window, self.wp=self.init_window(); self.window.title('SCHSIM Visualization : '+self.run+' (Author: Z. WANG)')
        self.hold='off'  #animation
        self.play='off'  #animation
        self.curve_method=0 #the method in extracting time series (0: nearest, 1: interpolation)
        self.window.mainloop()
        #todo: 1). for cases that out2d*.nc not exist

    def init_plot(self,fmt=0):
        w=self.wp; fn=w.fn.get()
        if fn=='add': #new figure
            if fmt==1: return #return if entry from control window
            p=zdata(); self.get_param(p); p.xm,p.ym=self.xm,self.ym; p.hf=figure(figsize=p.figsize,num=self.nf())
            cid=len(self.fns); self.figs.append(p); self.fns.append('{}: {}'.format(len(self.fns)+1,p.var))
        else: #old figure
            cid=self.fns.index(fn); p=self.figs[cid]
            if not fignum_exists(p.hf.number): p.hf=figure(figsize=p.figsize,num=p.hf.number) #restore closed figure
            if fmt==0: #modify current figure
                self.get_param(p); self.fns[cid]=self.fns[cid][:3]+p.var; p.hf.clear(); p.bm=None
            elif fmt==1: #restore old figure setting
                self.update_panel('old',p)
        w._fn['values']=['add',*self.fns]; w.fn.set(self.fns[cid])
        self.fig=p; figure(p.hf.number) #bring figure to front
        return p

    def schism_plot(self,fmt=0):
        if not hasattr(self,'hgrid'): print('wait: still reading grid'); return
        if self.play=='on' and fmt==1: self.play='off'; self.hold='off'; return
        if self.hold=='on': return #add hold to avoid freeze when user press too frequently
        w=self.wp; gd=self.hgrid
        p=self.init_plot(0) if fmt==0 else self.init_plot(1)
        if p is None: return
        if fmt==2: p.it=max(p.it-p.ns,0)
        if fmt==3: p.it=min(p.it+p.ns,len(self.irec)-1)
        if fmt==4: p.it=0
        if fmt==5: it=len(self.irec)-1; p.it=it; p.it2=it; self.update_panel('it2',p)

        #plot figure and save the backgroud
        self.hold='on'
        if fmt==0:
           p.hp=[]; p.hg=[]; p.hb=[]; p.hv=[]; anim=True if p.med==0 else False
           if p.var!='none':
               self.get_data(p); v=self.data
               if p.med==0: p.hp=[gd.plot(fmt=1,method=1,value=v,clim=p.vm,ticks=11,animated=True,cmap='jet',zorder=1,cb_aspect=50)]
               if p.med==1: p.hp=[gd.plot(fmt=1,method=0,value=v,clim=p.vm,ticks=11,cmap='jet',zorder=1,cb_aspect=50)]
           if p.vvar!='none': u,v=self.get_vdata(p); p.hv=[quiver(p.vx,p.vy,u,v,animated=anim,scale=1.0/p.zoom,scale_units='inches',width=0.001,zorder=3)]
           if p.vvar!='none': quiverkey(p.hv[0], X=0.92, Y=1.01, U=1, label='1.0 m/s',color='r', labelpos='E',zorder=4)
           if p.grid==1: hg=gd.plot(animated=anim,zorder=2); p.hg=[hg[0][0],*hg[1]]
           if p.bnd==1: p.hb=gd.plot_bnd(lw=0.5,alpha=0.5,animated=anim)
           p.ht=title('{}, layer={}, {}'.format(p.var,p.layer,self.mls[p.it]),animated=anim)

           #add pts for time series
           m=20; n=p.npt; x=array([*p.px,*tile(0.0,m-n)]); y=array([*p.py,*tile(nan,m-n)])
           fpn=nonzero(~((x[:n]>p.xm[0])*(x[:n]<p.xm[1])*(y[:n]>p.ym[0])*(y[:n]<p.ym[1])))[0]; x[fpn]=0.0; y[fpn]=nan
           p.hpt=plot(x,y,'r.',ms=6,alpha=0.75,animated=anim)
           for i in arange(m):
               [xi,yi,k]=[x[i],y[i],str(i+1)] if (i<n and (i not in fpn)) else [0,0,'']
               p.hpt.append(text(xi,yi,k,color='r',animated=anim))
           setp(gca(),xlim=p.xm,ylim=p.ym); gcf().tight_layout(); p.ax=gca(); pause(0.05)

           #associcate with actions
           p.hf.canvas.mpl_connect("draw_event", self.update_panel)
           p.hf.canvas.mpl_connect("button_press_event", self.onclick)
           p.hf.canvas.mpl_connect("button_press_event", self.query)
           #p.hf.canvas.mpl_connect('motion_notify_event', self.onmove) #todo: this fun is not ready yet, as it cause screen freeze
           if p.med==0: p.bm=blit_manager([p.ht,*p.hp,*p.hg,*p.hb,*p.hv,*p.hpt],p.hf); p.bm.update()
           self.update_panel('it',p)
           if hasattr(p,'qxy'): self.query(0)

        #animation
        if fmt!=0 and (p.var not in ['depth','none',*self.gr3] or p.vvar!='none'):
            if fmt==1: w.player['text']='stop'; self.window.update(); self.play='on'; it0=p.it; its=arange(it0+p.ns,p.it2,p.ns)
            if fmt in [2,3,4,5]: its=[p.it]; self.play='on'
            if p.anim!=None: savefig('.{}_{:06}'.format(p.anim,p.it)) #savefig for animation
            for p.it in its:
                if self.play=='off': break
                if p.var not in ['depth','none']: # contourf
                    self.get_data(p); v=self.data; self.query(1)
                    if p.med==0:
                        p.hp[0].set_array(v if v.size==gd.np else r_[v,v[self.fp4]])
                    else:
                        for i in arange(len(p.ax.collections)): p.ax.collections.pop()
                        gd.plot(ax=p.ax,fmt=1,value=v,clim=p.vm,ticks=11,cmap='jet',cb=False,zorder=1)
                if p.vvar!='none':  #vector
                   u,v=self.get_vdata(p)
                   if p.med==0: p.hv[0].set_UVC(u,v)
                   if p.med==1: p.hv=[quiver(p.vx,p.vy,u,v,scale=1/p.zoom,scale_units='inches',width=0.001,zorder=3)]
                p.ht.set_text('{}, layer={}, {}'.format(p.var,p.layer,self.mls[p.it]))
                self.update_panel('it',p); self.window.update()
                if p.med==0: p.bm.update()
                if p.anim!=None: savefig('.{}_{:06}'.format(p.anim,p.it)) #save fig for animation
                if p.med==1: pause(0.1)
                if hasattr(p,'pause'): pause(max(p.pause,0.0001))
                if self.play=='off': break
                if fmt in [2,3,4,5]: self.play='off'
            if fmt==1: w.player['text']='play'; self.window.update()
            if p.anim!=None:
               from PIL import Image
               ims=['.{}_{:06}.png'.format(p.anim,i) for i in  [it0,*its]]; fms=[Image.open(i) for i in ims]; adt=max(p.pause*1e3,50) if hasattr(p,'pause') else 200
               fms[0].save(p.anim+'.gif',format='GIF', append_images=fms[1:], save_all=True, duration=adt, loop=0)
               [os.remove(i) for i in ims]
        self.hold='off'

    def plotts(self):
        import threading
        #function to extract data
        def get_tsdata(ts,x,y,svar,layer,ik1,ik2):
            w.curve['text']='wait'; ts.x=x; ts.y=y; ts.var=svar; ts.layer=layer; ts.ik1=ik1; ts.ik2=ik2; ts.mys=[]; nt=0
            for ik in arange(ik1,ik2+1):
                fname='{}/out2d_{}.nc'.format(self.outputs,ik) if svar in self.vars_2d else '{}/{}_{}.nc'.format(self.outputs,svar,ik)
                C=self.fid(fname); nt0,npt=C.variables[svar].shape[:2]; t00=time.time()
                if ik==ik1 and self.curve_method==0: sindp=near_pts(c_[x,y],c_[gd.x,gd.y]) if npt==gd.np else near_pts(c_[x,y],c_[gd.xctr,gd.yctr]) #compute index
                if ik==ik1 and self.curve_method==1: pie,pip,pacor=gd.compute_acor(c_[x,y],fmt=1); sindp=pip.ravel() if npt==gd.np else pie #compute index for interp
                if svar in self.vars_2d:
                    data=array(C.variables[svar][:,sindp])
                else:
                    ks=(self.kbp[sindp] if npt==gd.np else self.kbe[sindp]) if layer=='bottom' else (-tile(1 if layer=='surface' else int(layer),sindp.size))
                    data=array([C.variables[svar][:,i,k] for i,k in zip(sindp,ks)]).T
                if npt==gd.np and self.curve_method==1: data=sum(reshape(data,[nt0,*pip.shape])*pacor[None,...],axis=2)
                ts.mys.extend(data); nt=nt+nt0; print('extracting {} from {}: {:0.2f}'.format(svar,fname,time.time()-t00))
            ts.mys=array(ts.mys).T; ts.mt=array(self.mts[it1:(it1+nt)]); ts.mls=array(self.mls[it1:(it1+nt)]); p.ts=ts
            print('done in extracting'); w.curve['text']='curve'

        def update_xts(event):
            if event!=0 and type(event)!=mpl.backend_bases.DrawEvent: return
            t1,t2=xlim(); dt1=abs(mt-t1); dt2=abs(mt-t2); i1=nonzero(dt1==dt1.min())[0][0]; i2=nonzero(dt2==dt2.min())[0][0]
            ns=max(int(floor((i2-i1+1)/5)),1); mti=mt[i1:i2:ns]; mlsi=mls[i1:i2:ns]
            if hasattr(self,'StartT'): mlsi=[i[:10]+'\n'+i[11:] for i in mlsi]
            s.ax.set_xticks(mti); s.ax.set_xticklabels(mlsi)

        #prepare info. about time sereis
        if not hasattr(self,'fig'): return
        p=self.fig; w=self.wp; gd=self.hgrid; gd.compute_ctr()
        svar,layer=p.var,p.layer; x=array(p.px); y=array(p.py)
        if svar=='depth' or len(x)==0: return
        ik1=self.istack[p.it]; ik2=self.istack[p.it2-1]; it1=self.istack.index(ik1); it2=len(self.istack)-self.istack[::-1].index(ik2)
        fpc=(array_equal(x,p.ts.x) and array_equal(y,p.ts.y) and svar==p.ts.var and layer==p.ts.layer and ik1>=p.ts.ik1 and ik2<=p.ts.ik2) if hasattr(p,'ts') else False
        if not fpc: ts=zdata(); threading.Thread(target=get_tsdata,args=(ts,x,y,svar,layer,ik1,ik2)).start(); return

        #plot time series
        s=p.ts; mt=s.mt; mls=s.mls; s.hf=figure(figsize=[6.5,3.5],num=self.nf()); cs='rgbkcmy'; ss=['-',':','--']; lstr=[]
        for n,my in enumerate(s.mys):
            plot(mt,my,color=cs[n%7],ls=ss[int(n/7)]); lstr.append('pt_{}'.format(n+1))
        s.ax=gca(); ym=ylim(); plot(mt,zeros(mt.size),'k:',lw=0.3,alpha=0.5)
        setp(s.ax,xticks=mt[:2],xticklabels=mls[:2],xlim=[mt.min(),mt.max()],ylim=ym); s.ax.xaxis.grid('on')
        title('{}, layer={}, ({}, {})'.format(s.var,s.layer,mls[0][:10],mls[-1][:10]))
        legend(lstr); s.hf.tight_layout(); show(block=False)
        update_xts(0); s.hf.canvas.mpl_connect("draw_event", update_xts)

    def plotsc(self):
        print('profile function not available yet'); return

    def query(self,sp=None):
        if not (hasattr(self,'fig') and hasattr(self.fig,'hpt') and hasattr(self.fig,'hf')): return
        p=self.fig; qr=self.wp.query; hf=p.hf; gd=self.hgrid
        hp=p.hpt[0]; ht=p.hpt[-1]; x,y=hp.get_xdata(),hp.get_ydata()

        if sp is None: #set query state
           import tkinter as tk
           qr.config(relief=tk.RAISED) if qr['relief'].lower()=='sunken' else qr.config(relief=tk.SUNKEN)
           if hasattr(p,'qxy'): delattr(p,'qxy')
           x[-1]=nan; y[-1]=nan; hp.set_xdata(x); hp.set_ydata(y); ht.set_text(''); hf.canvas.draw()
        elif qr['relief'].lower()=='sunken':
            if sp not in [0,1]: #save query pt info
               dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
               if not (dlk==0 and btn==2): return
               p.qxy=[bx,by,*gd.compute_acor(c_[bx,by])]
            #update query
            bx,by,pie,pip,pacor=p.qxy
            x[-1]=bx; y[-1]=by; hp.set_xdata(x); hp.set_ydata(y); ht.set_x(bx); ht.set_y(by)
            pv=self.data[pie[0]] if self.data.size==gd.ne else sum(self.data[pip]*pacor,axis=1)[0]
            ht.set_text('{:0.6f}'.format(pv) if pv>1e-6 else '{:0.6e}'.format(pv))
            if sp!=1: ht.set_backgroundcolor([1,1,1,0.5]); ht.set_color('k'); hf.canvas.draw()

    def onclick(self,sp):
        if sp.button is None: return
        if not fignum_exists(self.fig.hf.number): return
        p=self.fig; dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
        if dlk==1 and btn==1: #animation on and off
            if self.play=='on':
                self.play='off'
            else:
                self.schism_plot(1)
        elif dlk==1 and btn==3: #add pts
            hp=p.hpt[0]; x,y=hp.get_xdata(),hp.get_ydata(); x[p.npt]=bx; y[p.npt]=by; hp.set_xdata(x); hp.set_ydata(y)
            ht=p.hpt[p.npt+1]; ht.set_x(bx); ht.set_y(by); ht.set_text(str(p.npt+1))
            p.px.append(bx); p.py.append(by); p.npt=p.npt+1
        elif dlk==1 and btn==2: #remove pts
            if p.npt==0: return
            p.npt=p.npt-1; p.px.pop(); p.py.pop()
            hp=p.hpt[0]; y=hp.get_ydata(); y[p.npt]=nan; hp.set_ydata(y); p.hpt[p.npt+1].set_text('')
        if p.med==0: p.bm.update()
        if p.med==1: p.hf.canvas.draw()

    def onmove(self,sp): #move pts
        if sp.button is None: return
        p=self.fig; dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
        if dlk==0 and btn==2 and p.npt!=0:
            x=array(p.px); y=array(p.py); distp=squeeze(abs((x-bx)+1j*(y-by))); sid=nonzero(distp==distp.min())[0][0]
            hp=p.hpt[0]; x,y=hp.get_xdata(),hp.get_ydata(); x[sid]=bx; y[sid]=by; hp.set_xdata(x); hp.set_ydata(y)
            p.hpt[sid+1].set_x(bx); p.hpt[sid+1].set_y(by); p.px[sid]=bx; p.py[sid]=by
        if p.med==0: p.bm.update(); pause(0.001)
        if p.med==1: p.hf.canvas.draw()

    def get_data(self,p):  #slab data
        svar,layer,istack,irec=p.var,p.layer,self.istack[p.it],self.irec[p.it]; gd=self.hgrid
        if p.var=='depth': self.data=gd.dp; return
        if p.var in self.gr3: self.data=read_schism_hgrid(self.runpath+os.sep+p.var).z; return
        C=self.fid('{}/out2d_{}.nc'.format(self.outputs,istack) if svar in self.vars_2d else '{}/{}_{}.nc'.format(self.outputs,svar,istack))
        if svar in self.vars_2d:
            data=array(C.variables[svar][irec])
        else:
            if layer=='bottom':
                npt=C.variables[svar].shape[1]
                if npt==gd.np: data=array(C.variables[svar][irec][arange(gd.np),self.kbp])
                if npt==gd.ne: data=array(C.variables[svar][irec][arange(gd.ne),self.kbe])
            else:
                layer=1 if layer=='surface' else int(layer); data=array(C.variables[svar][irec,:,-layer])
        data[abs(data)>1e20]=nan
        self.data=data

    def get_vdata(self,p): #get vector data
        svar,layer,istack,irec=p.vvar,p.layer,self.istack[p.it],self.irec[p.it]; gd=self.hgrid
        if svar+'X' in self.vars_2d:
            C=self.fid('{}/out2d_{}.nc'.format(self.outputs,istack)); npt=C.variables[svar+'X'].shape[1]; fmt=0
        else:
            C=[self.fid('{}/{}X_{}.nc'.format(self.outputs,svar,istack)),self.fid('{}/{}Y_{}.nc'.format(self.outputs,svar,istack))]; npt=C[0].variables[svar+'X'].shape[1]; fmt=1
        if npt==gd.ne and (not hasattr(gd,'xctr')): gd.compute_ctr()
        if npt==gd.ns and (not hasattr(gd,'xcj')):  gd.compute_side(fmt=2)
        if npt==gd.ns and (not hasattr(self,'kbs')): self.kbs=gd.compute_kb(self.kbp,fmt=1)
        if p.sindv is None: #get xy coordinate and indices
           x,y=[gd.x,gd.y] if npt==gd.np else [gd.xctr,gd.yctr] if npt==gd.ne else [gd.xcj,gd.ycj]
           p.sindv=nonzero((x>=p.xm[0])*(x<=p.xm[1])*(y>=p.ym[0])*(y<=p.ym[1]))[0]; p.vx=x[p.sindv]; p.vy=y[p.sindv]

        sind=p.sindv #read record of vector variables
        if fmt==0:
            u=array(C.variables[svar+'X'][irec][sind]); v=array(C.variables[svar+'Y'][irec][sind])
        else:
            if layer=='bottom':
                kb=(self.kbp if npt==gd.np else self.kbe if npt==gd.ne else self.kbs)[sind]
                u=array(C[0].variables[svar+'X'][irec][sind,kb]); v=array(C[1].variables[svar+'Y'][irec][sind,kb])
            else:
                layer=1 if layer=='surface' else int(layer); u,v=array(C[0].variables[svar+'X'][irec,:,-layer][sind]),array(C[1].variables[svar+'Y'][irec,:,-layer][sind])
        return [u,v]

    def update_panel(self,event,p=None):  #update control panel
        w=self.wp
        if event in ['time','stack','julian','it','it2']: #set time variables
            mls=self.mls if w.time.get()=='time' else self.julian if w.time.get()=='julian' else self.stacks
            if event=='it':
                if w.time.get()=='stack': mls=self.istack
                w.StartT.set(mls[p.it]); p.StartT=mls[p.it]
            elif event=='it2':
                w.EndT.set(mls[-1]); p.EndT=mls[-1]
            else:
                w.StartT.set(mls[0]); w.EndT.set(mls[-1]); w._StartT['values']=mls; w._EndT['values']=mls; w._StartT['width']=6; w._EndT['width']=6
                if event=='time': w._StartT['width']=18; w._EndT['width']=18
        elif event=='vm': #reset data limit
            if not hasattr(self,'hgrid'): return
            p=self.get_param()
            if p.var!='none': self.get_data(p); w.vmin.set(self.data.min()); w.vmax.set(self.data.max())
        elif event=='old': #reset all panel variables from a previous figure
            w.var.set(p.var); w.fn.set(p.fn); w.layer.set(p.layer); w.time.set(p.time)
            w.StartT.set(p.StartT); w.EndT.set(p.EndT); w.vmin.set(p.vm[0]); w.vmax.set(p.vm[1])
            w.xmin.set(p.xm[0]); w.xmax.set(p.xm[1]); w.ymin.set(p.ym[0]); w.ymax.set(p.ym[1])
            w.ns.set(p.ns); w.grid.set(p.grid); w.bnd.set(p.bnd); w.vvar.set(p.vvar); w.zoom.set(p.zoom); w.med.set(p.med)
        elif type(event)==mpl.backend_bases.DrawEvent:
            if not fignum_exists(self.fig.hf.number): return
            p=self.fig; ax=p.ax; xm=ax.get_xlim(); ym=ax.get_ylim(); p.xm=[*xm]; p.ym=[*ym]
            w=self.wp; w.xmin.set(xm[0]); w.xmax.set(xm[1]); w.ymin.set(ym[0]); w.ymax.set(ym[1])
            p.figsize=[p.hf.get_figwidth(),p.hf.get_figheight()]
        self.window.update()

    def get_param(self,p=None):
        if p is None: p=zdata()
        w=self.wp; p.var=w.var.get(); p.fn=w.fn.get(); p.layer=w.layer.get(); p.time=w.time.get()
        p.StartT=w.StartT.get(); p.EndT=w.EndT.get(); p.vm=[w.vmin.get(),w.vmax.get()]
        p.xm=[w.xmin.get(),w.xmax.get()]; p.ym=[w.ymin.get(),w.ymax.get()]; p.med=w.med.get()
        p.ns=w.ns.get(); p.grid=w.grid.get(); p.bnd=w.bnd.get(); p.vvar=w.vvar.get(); p.zoom=w.zoom.get()
        p.anim=None; p.sindv=None; p.figsize=[7.2,5.5]
        if not hasattr(p,'npt'): p.npt=0; p.px=[]; p.py=[]

        #get time index
        if p.time=='time':
            p.it=self.mls.index(p.StartT); p.it2=self.mls.index(p.EndT)+1
        elif p.time=='julian':
            p.it=self.julian.index(float(p.StartT)); p.it2=self.julian.index(float(p.EndT))+1
        else: #stacks
            p.it=self.istack.index(int(p.StartT)); p.it2=len(self.istack)-self.istack[::-1].index(int(p.EndT))
        return p

    def fid(self,fname): #output chanenl
        if not hasattr(self,'_fid'): self._fid={}
        if fname not in self._fid: self._fid[fname]=ReadNC(fname,1)
        return self._fid[fname]

    def nf(self): #get a new figure number
        self._nf=self._nf+1
        return self._nf

    def window_exit(self):
        if hasattr(self,'_fid'):
            for i in self._fid: self._fid[i].close()
        for p in self.sbp:
            try:
               p.kill()
            except:
               pass
        if os.path.exists(os.curdir+os.sep+'.schism_check'): os.remove(os.curdir+os.sep+'.schism_check')
        if os.path.exists(os.curdir+os.sep+'.schism_compare'): os.remove(os.curdir+os.sep+'.schism_compare')
        close('all'); self.window.destroy()

    @property
    def INFO(self):
        return get_INFO(self)
    @property
    def VINFO(self):
        return get_INFO(self)

    def run_info(self,run):
        import threading,time

        #check output
        print('reading grid and output info.')
        fns=glob(run+'/out2d_*.nc'); fns2=glob(run+'/outputs/out2d_*.nc'); self.gr3=[os.path.basename(i) for i in glob(run+'/*.gr3')]; iout=1
        self.outputs,fnames=[run,fns] if len(fns)!=0 else [run+'/outputs',fns2]; self.runpath=os.path.abspath(run); self.run=os.path.basename(self.runpath)
        if len(fnames)!=0:
           [iks,self.vars,self.vars_2d]=get_schism_output_info(self.outputs)[2:]
           iks=sort(iks); ik0=iks[0]; C=self.fid('{}/out2d_{}.nc'.format(self.outputs,ik0)); cvar=C.variables; cdim=C.dimensions
           self.vvars=[i[:-1] for i in self.vars if (i[-1]=='X') and (i[:-1]+'Y' in self.vars)]
        else:
           self.vars=[]; self.svars_2d=[]; self.vvars=[]; iout=0; cvar=None
           print('schism outputs dir. not found')

        #read grid and param
        grd=run+'/grid.npz'; gr3=run+'/hgrid.gr3'; vrd=run+'/vgrid.in'; par=run+'/param.nml'
        if os.path.exists(par): p=read_schism_param(par,3); self.param=p; self.StartT=datenum(p.start_year,p.start_month,p.start_day,p.start_hour)
        def _read_grid():
           gd=loadz(grd).hgrid if os.path.exists(grd) else read_schism_hgrid(gr3) if os.path.exists(gr3) else None
           vd=loadz(grd).vgrid if os.path.exists(grd) else read_schism_vgrid(vrd) if os.path.exists(vrd) else None
           if gd==None: #create hgrid
              if cvar==None: return
              gd=get_schism_output_info(self.outputs,4)
           self.hgrid=gd; self.xm=[gd.x.min(),gd.x.max()]; self.ym=[gd.y.min(),gd.y.max()]; self.vm=[gd.dp.min(),gd.dp.max()]; self.fp3=nonzero(gd.i34==3)[0]; self.fp4=nonzero(gd.i34==4)[0]
           self.kbp, self.nvrt=[vd.kbp, vd.nvrt] if vd!=None else [array(cvar['bottom_index_node']), cdim['nSCHISM_vgrid_layers'].size]; self.kbe=gd.compute_kb(self.kbp)
           while not hasattr(self,'wp'): time.sleep(0.01)
           w=self.wp; w._layer['values']=['surface','bottom',*arange(2,self.nvrt+1)]; print('schismview ready')
           w.vmin.set(self.vm[0]); w.vmax.set(self.vm[1]); w.xmin.set(self.xm[0]); w.xmax.set(self.xm[1]); w.ymin.set(self.ym[0]); w.ymax.set(self.ym[1])
        self.nvrt=2; self.xm=[0,1]; self.ym=[0,1]; self.vm=[0,1]
        threading.Thread(target=_read_grid).start()

        #read available time
        if iout==0: self.stacks=[0]; self.julian=[0]; self.istack=[0]; self.irec=[0]; self.mls=['0']; return
        self.stacks=[]; self.julian=[]; self.istack=[]; self.irec=[]
        ti=array(cvar['time'])/86400; nt=len(ti); t0=ti[0]  #assume all stacks have the same number of records
        if (not hasattr(self,'StartT')) and hasattr(C.variables['time'],'base_date'):
           self.StartT=datenum(*[int(float(i)) for i in C.variables['time'].base_date.split()])
        if len(iks)>1 and nt==1:
            C1=self.fid('{}/out2d_{}.nc'.format(self.outputs,iks[1])); ti1=array(C1.variables['time'])/86400; dt=ti1-ti[0]
        else:
            dt=0 if nt==0 else ti[1]-ti[0]
        for ik in iks:
            fn='{}{}out2d_{}.nc'.format(self.outputs,os.path.sep,ik); nt1=nt
            if fn not in fnames: continue
            if ik==iks[-1]:
               try:
                   C=self.fid(fn); nt1=C.dimensions['time'].size
                   if nt1==0: continue
               except:
                   pass
            self.julian.extend(t0+(ik-ik0)*nt*dt+arange(nt1)*dt); self.irec.extend(arange(nt1))
            self.stacks.append(ik); self.istack.extend(tile(ik,nt1))

        self.mts=self.julian[:]; self.mls=['{}'.format(i) for i in self.mts]
        if hasattr(self,'StartT'): #get time_string
           def _set_time():
               self.mls=[i.strftime('%Y-%m-%d, %H:%M') for i in num2date(self.mts)]
               while not hasattr(self,'wp'): time.sleep(0.01)
               self.wp._StartT['values']=self.mls; self.wp._EndT['values']=self.mls
           self.mts=[*(array(self.mts)+self.StartT)]
           self.mls[0]=num2date(self.mts[0]).strftime('%Y-%m-%d, %H:%M')
           self.mls[-1]=num2date(self.mts[-1]).strftime('%Y-%m-%d, %H:%M')
           threading.Thread(target=_set_time).start()

    def cmd_window(self):
        import tkinter as tk
        from tkinter import ttk
        cw=tk.Toplevel(self.window); cw.geometry("400x200"); cw.title('command input')
        cw.rowconfigure(0,minsize=150, weight=1); cw.columnconfigure(0,minsize=2, weight=1)
        txt=tk.Text(master=cw,width=150,height=14); txt.grid(row=0,column=0,pady=2,padx=2,sticky='nsew')
        rbn=ttk.Button(cw, text= "run",command=lambda: self.cmd_exec(txt.get('1.0',tk.END))); rbn.grid(row=1,column=0,padx=10)
        cw.update(); xm=max(txt.winfo_width(),rbn.winfo_width()); ym=txt.winfo_height()+rbn.winfo_height()+12
        if hasattr(self,'cmd'): txt.insert('1.0',self.cmd)
        cw.geometry('{}x{}'.format(xm,ym)); cw.update()

    def cmd_exec(self,cmd):
        self.cmd=cmd
        if hasattr(self,'fig'): #get figure handles
           p=self.fig; fig, hf, ax, hp, hv, hg, hb, hpt = p, p.hf, p.ax, p.hp, p.hv, p.hg, p.hb, p.hpt
           hp=None if len(hp)==0 else hp[0]; hv=None if len(hv)==0 else hv[0]; hb=None if len(hb)==0 else hb[0]
        for i in cmd.strip().split('\n'): #run command
            if i=='': continue
            try:
               print('run: '+i); exec(i)
            except:
               print('fail: '+i)
        if hasattr(self,'fig'): self.fig.hf.canvas.draw()
        self.window.update()

    def anim_window(self):
        import tkinter as tk
        from tkinter import ttk
        if not hasattr(self,'fig'): return
        p=self.fig; p._anim=tk.StringVar(self.window)
        cw=tk.Toplevel(self.window); cw.geometry("300x30"); cw.title('animation')
        fm=tk.Frame(master=cw); fm.grid(row=0,column=0)
        ttk.Label(master=fm,text='name').grid(row=0,column=0)
        ttk.Entry(master=fm,textvariable=p._anim,width=20).grid(row=0,column=1,pady=5,padx=3)
        rbn=ttk.Button(fm, text= "save",command=self.anim_exec,width=6); rbn.grid(row=0,column=2)
        cw.update(); cw.geometry('{}x{}'.format(fm.winfo_width()+12,rbn.winfo_height()+5)); cw.update()

    def show_node(self):
        p=self.fig; gd=self.hgrid
        if hasattr(p,'hns'):
           for i in arange(len(p.hns)): p.hns.pop().remove()
           delattr(p,'hns'); p.hf.canvas.draw()
        else:
           xm=xlim(); ym=ylim(); p.hns=[]
           if not hasattr(gd,'xctr'): gd.compute_ctr()
           sind=nonzero((gd.x>=xm[0])*(gd.x<=xm[1])*(gd.y>=ym[0])*(gd.y<=ym[1]))[0]
           for i in sind: ht=text(gd.x[i],gd.y[i],'{}'.format(i+1),fontsize=6,zorder=3); p.hns.append(ht)
           sind=nonzero((gd.xctr>=xm[0])*(gd.xctr<=xm[1])*(gd.yctr>=ym[0])*(gd.yctr<=ym[1]))[0]
           for i in sind: ht=text(gd.xctr[i],gd.yctr[i],'{}'.format(i+1),fontsize=6,zorder=3); p.hns.append(ht)
           p.hf.canvas.draw()
        return

    def anim_exec(self):
        p=self.fig; anim=self.fig._anim.get()
        p.anim=anim[:-4] if anim.endswith('.gif') else anim if anim.strip()!='' else None
        self.play='off'; self.schism_plot(1); p.anim=None

    def schism_compare(self):
        from tkinter import filedialog
        import subprocess
        crun=filedialog.askdirectory(initialdir = self.runpath,title = "choose run to compare")
        cfile='.schism_compare'; fid=open(cfile,'w+'); os.chmod(cfile,0o777)
        fid.write("#!/usr/bin/env python3\nfrom pylib import *\nmpl.use('TkAgg')\nschism_view('{}')".format(crun)); fid.close()
        p=subprocess.Popen(os.curdir+os.sep+cfile); self.sbp.append(p)

    def schism_check(self):
        import subprocess
        cfile='.schism_check'; fid=open(cfile,'w+'); os.chmod(cfile,0o777)
        fid.write("#!/usr/bin/env python3\nfrom pylib import *\nmpl.use('TkAgg')\nschism_check('{}')".format(self.runpath)); fid.close()
        p=subprocess.Popen(os.curdir+os.sep+cfile); self.sbp.append(p)

    def init_window(self):
        #open an window
        import tkinter as tk
        from tkinter import ttk

        wd=tk.Tk(); fms=[]
        wd.title("SCHSIM Visualization")
        wd.rowconfigure([0,1,2], minsize=5, weight=1)
        wd.columnconfigure(0, minsize=20, weight=1)
        w=zdata() #capsule used to store parameters

        #frame 1
        fm=ttk.Frame(master=wd); fm.grid(row=0,column=0,sticky='NW',pady=10); fms.append(fm)

        #variable
        w.var=tk.StringVar(wd); w.var.set('depth')
        ttk.Label(master=fm,text='  variable').grid(row=0,column=0,sticky='W',pady=4)
        svar=ttk.Combobox(fm,textvariable=w.var,values=['none','depth',*self.vars,*self.gr3],width=15,); svar.grid(row=0,column=1)
        svar.bind("<<ComboboxSelected>>",lambda x: self.update_panel('vm'))

        #figure
        sfm=ttk.Frame(master=fm); sfm.grid(row=0,column=2,sticky='E',padx=5)
        ttk.Label(master=sfm,text='        figure').grid(row=0,column=0,sticky='E',padx=2)
        w.fn=tk.StringVar(wd); w.fn.set('add')
        w._fn=ttk.Combobox(sfm,textvariable=w.fn,values=['add'],width=10); w._fn.grid(row=0,column=1)
        w._fn.bind("<<ComboboxSelected>>",lambda x: self.init_plot(1))

        #layer
        w.layer=tk.StringVar(wd); w.layer.set('surface')
        ttk.Label(master=fm,text='  layer').grid(row=1,column=0,sticky='W',pady=4)
        w._layer=ttk.Combobox(fm,textvariable=w.layer,values=['surface','bottom',*arange(2,self.nvrt+1)],width=15); w._layer.grid(row=1,column=1)

        #grid, bnd, method
        sfm2=ttk.Frame(master=fm); sfm2.grid(row=1,column=2)
        w.grid=tk.IntVar(wd); w.grid.set(0); w.bnd=tk.IntVar(wd); w.bnd.set(0); w.med=tk.IntVar(wd); w.med.set(0)
        tk.Checkbutton(master=sfm2,text='grid',variable=w.grid,onvalue=1,offvalue=0).grid(row=0,column=0)
        tk.Checkbutton(master=sfm2,text='bnd',variable=w.bnd,onvalue=1,offvalue=0).grid(row=0,column=1,sticky='W')
        tk.Checkbutton(master=sfm2,text='ctr',variable=w.med,onvalue=1,offvalue=0).grid(row=0,column=2,sticky='W')

        #time
        w.time=tk.StringVar(wd); w.StartT=tk.StringVar(wd); w.EndT=tk.StringVar(wd); w.mls=self.mls; w.StartT.set(self.mls[0]); w.EndT.set(self.mls[-1])
        ttk.OptionMenu(fm,w.time,'time','time','stack','julian',command=self.update_panel).grid(row=2,column=0,sticky='W',pady=4)
        w._StartT=ttk.Combobox(master=fm,textvariable=w.StartT,values=self.mls,width=18); w._StartT.grid(row=2,column=1,padx=0,sticky='W')
        w._EndT=ttk.Combobox(master=fm,textvariable=w.EndT,values=self.mls,width=18); w._EndT.grid(row=2,column=2,sticky='W',padx=1)

        #limit
        w.vmin=tk.DoubleVar(wd); w.vmax=tk.DoubleVar(wd); w.vmin.set(self.vm[0]); w.vmax.set(self.vm[1])
        ttk.Label(master=fm,text='  limit').grid(row=3,column=0,sticky='W')
        ttk.Entry(fm,textvariable=w.vmin,width=10).grid(row=3,column=1,sticky='W',padx=2)
        ttk.Entry(fm,textvariable=w.vmax,width=10).grid(row=3,column=2,sticky='W')

        #frame2: vector, time_series
        fm=ttk.Frame(master=wd); fm.grid(row=1,column=0,sticky='NW'); fms.append(fm)

        #vector
        w.vvar=tk.StringVar(wd); w.zoom=tk.DoubleVar(wd); w.vvar.set('none'); w.zoom.set(1.0)
        fm3=ttk.Frame(master=fm); fm3.grid(row=0,column=0,sticky='W')
        ttk.Label(master=fm3,text='  vector ').grid(row=1,column=0,sticky='W',pady=4)
        vvar=ttk.Combobox(fm3,textvariable=w.vvar,values=['none',*self.vvars],width=14,); vvar.grid(row=1,column=1)
        ttk.Entry(fm3,textvariable=w.zoom,width=5).grid(row=1,column=2)

        #time series
        fm0=ttk.Frame(master=fm); fm0.grid(row=0,column=1)
        w.curve=ttk.Button(master=fm0,text='curve',command=self.plotts,width=5); w.curve.grid(row=0,column=1)
        w.query=tk.Button(master=fm0,text='query',bg='gray88',command=self.query,width=5); w.query.grid(row=0,column=2,padx=1,pady=2)
        #tk.Button(master=fm0,text='profile',bg='darkgray',command=self.plotsc,width=7).grid(row=0,column=3,padx=1,pady=2)

        #xlim, ylim
        w.xmin=tk.DoubleVar(wd); w.xmax=tk.DoubleVar(wd); w.xmin.set(self.xm[0]); w.xmax.set(self.xm[1])
        w.ymin=tk.DoubleVar(wd); w.ymax=tk.DoubleVar(wd); w.ymin.set(self.ym[0]); w.ymax.set(self.ym[1])
        fm1=ttk.Frame(master=fm); fm1.grid(row=2,column=0); fm2=ttk.Frame(master=fm); fm2.grid(row=2,column=1,padx=10)
        ttk.Label(fm1,text='  xlim',width=6).grid(row=0,column=0,sticky='E')
        ttk.Entry(fm1,textvariable=w.xmin,width=11).grid(row=0,column=1,sticky='W')
        ttk.Entry(fm1,textvariable=w.xmax,width=11).grid(row=0,column=2,sticky='W',padx=2)
        ttk.Label(fm2,text='ylim',width=4).grid(row=0,column=0,sticky='E')
        ttk.Entry(fm2,textvariable=w.ymin,width=11).grid(row=0,column=1,sticky='W')
        ttk.Entry(fm2,textvariable=w.ymax,width=11).grid(row=0,column=2,sticky='W',padx=2)

        #frame 3: control
        fm=ttk.Frame(master=wd); fm.grid(row=2,column=0,sticky='W',pady=2); fms.append(fm)
        sfm0=ttk.Frame(master=fm); sfm0.pack(side=tk.LEFT)
        ttk.Button(master=sfm0,text='exit',command=self.window_exit,width=5).pack(side=tk.LEFT)
        mbar=ttk.Menubutton(sfm0,text='option',width=6); mbar.pack(side=tk.LEFT)
        menu=tk.Menu(mbar,tearoff=0)
        menu.add_command(label="command", command=self.cmd_window)
        menu.add_command(label="save animation", command=self.anim_window)
        menu.add_command(label="show node/element", command=self.show_node)
        menu.add_command(label="schismcheck", command=self.schism_check)
        menu.add_command(label="compare", command=self.schism_compare)
        mbar['menu']=menu; mbar['direction']='below'

        sfm=ttk.Frame(master=fm); sfm.pack(side=tk.LEFT); w.ns=tk.IntVar(wd); w.ns.set(1)
        L1=ttk.Label(master=sfm,text=''); L1.grid(row=0,column=0,sticky='W')
        ttk.Button(master=sfm,text='|<',width=3, command=lambda: self.schism_plot(4)).grid(row=0,column=1,sticky='W',padx=1)
        ttk.Button(master=sfm,text='<',width=3, command=lambda: self.schism_plot(2)).grid(row=0,column=2,sticky='W',padx=1)
        w.player=ttk.Button(master=sfm,text='play',width=5,command=lambda: self.schism_plot(1)); w.player.grid(row=0,column=3,sticky='W',padx=1)
        ttk.Button(master=sfm,text='>',width=3,command=lambda: self.schism_plot(3)).grid(row=0,column=4,sticky='W',padx=1)
        ttk.Button(master=sfm,text='>|',width=3,command=lambda: self.schism_plot(5)).grid(row=0,column=5,sticky='W',padx=1)
        ttk.Label(master=sfm,text='  skip:',width=6).grid(row=0,column=6,sticky='E')
        ttk.Entry(sfm,textvariable=w.ns,width=3).grid(row=0,column=7,sticky='W')
        L2=ttk.Label(master=sfm,text=''); L2.grid(row=0,column=8,sticky='W')
        ttk.Button(master=fm,text='draw',width=4,command=lambda: self.schism_plot(0)).pack(side=tk.RIGHT,padx=1)

        #resize window
        wd.protocol("WM_DELETE_WINDOW",self.window_exit)
        wd.geometry('600x180'); wd.update(); xm=max([i.winfo_width() for i in fms])
        for i in arange(0,100,3):
            L1['text']=' '*i; L2['text']=' '*i; wd.update(); xs=fms[-1].winfo_width()
            if xs>xm: wd.geometry('{}x210'.format(xs)); wd.update(); break
        return wd,w

class schism_check(zdata):
   def __init__(self,run='.'):
       self.params,self.figs,self.fids,self.fmts={},{},{},{}

       self.run_info(run) #get run information
       self.init_window(); self.update_panel()
       self.window.mainloop()

   def update_panel(self,option=0):
       '''
       option=0: When input is changed
       option=1: When variable of input is changed
       '''
       import tkinter as tk
       from tkinter import ttk

       #get input fname and type
       wd,fm,params,fmts,ap=self.window,self.frame,self.params,self.fmts,self.ap
       fname=self.input.get(); self.fname=fname
       if fname not in params: p=zdata(); p.init=0; params[fname]=p
       if fname not in fmts:
           if fname.split('.')[-1] in ['gr3','ll','ic','prop']: fmts[fname]=0
           if fname.endswith('D.th.nc') or fname.endswith('_nu.nc'): fmts[fname]=1
           if fname.startswith('hotstart.nc') or fname=='ICM_param.nc': fmts[fname]=2
           if fname=='source.nc': fmts[fname]=3
           if fname=='source_input':
              sname='.source.nc'; fmts[fname]=3
              if not os.path.exists(self.run+sname): print('convert schism source_sink format: '+sname); convert_schism_source(self.run,sname)
           if fname=='*.th': fmts[fname]=4
           self.fmt=fmts[fname]; self.read_input_info()
       self.fmt=fmts[fname]; p=params[fname]

       for widget in fm.winfo_children(): widget.destroy() #clean plot option frame
       self.info.config(text=p.info); ap.dns=[]
       #design plot option frame
       if self.fmt==0: #update gr3 file parameters and panel
          if p.init==0:
             p.grid=tk.IntVar(wd); p.bnd=tk.IntVar(wd); p.ctr=tk.IntVar(wd); p.grid.set(0); p.bnd.set(0); p.ctr.set(1)
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.sflux=tk.StringVar(wd); p.vmin.set(0); p.vmax.set(0); p.sflux.set('None')

          #update plot options
          self.var.set('z'); self.vars['values']='z'; p.dnames=[]; p.dvars=[]
          sfm1=ttk.Frame(master=fm); sfm1.grid(row=0,column=0,sticky='W',pady=5)
          ttk.Label(master=sfm1,text='  ').grid(row=0,column=0,sticky='W')
          tk.Checkbutton(master=sfm1,text='grid',variable=p.grid,onvalue=1,offvalue=0).grid(row=0,column=1)
          tk.Checkbutton(master=sfm1,text='bnd',variable=p.bnd,onvalue=1,offvalue=0).grid(row=0,column=2,sticky='W')
          tk.Checkbutton(master=sfm1,text='contour',variable=p.ctr,onvalue=1,offvalue=0).grid(row=0,column=3,sticky='W')
          if self.sflux is not None:
             ttk.Label(master=sfm1,text=' sflux_grid').grid(row=0,column=4,sticky='W')
             ttk.Combobox(sfm1,textvariable=p.sflux,values=['None',*self.sflux],width=17).grid(row=0,column=5)

          #add limit panel
          sfm=ttk.Frame(master=fm); sfm.grid(row=1,column=0,sticky='W')
          ttk.Label(master=sfm,text='  limit').grid(row=0,column=0,sticky='W')
          ttk.Entry(sfm,textvariable=p.vmin,width=10).grid(row=0,column=1,sticky='W',padx=2)
          ttk.Entry(sfm,textvariable=p.vmax,width=10).grid(row=0,column=2,sticky='W')
       elif self.fmt==1: # *D.th.nc or _nu.nc
          if p.init==0:
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.vmin.set(0); p.vmax.set(0)
             p.transpose=tk.IntVar(wd); p.transpose.set(0)
             if fname.endswith('_nu.nc') or fname=='uv3D.th.nc':
                p.ctr=tk.IntVar(wd); p.grid=tk.IntVar(wd); p.bnd=tk.IntVar(wd); p.scale=tk.DoubleVar(wd)
                p.ctr.set(0); p.grid.set(0); p.bnd.set(0); p.scale.set(1.0)
             p.dvars=[tk.StringVar(wd) for i in p.dims]; p.xs=[[*arange(i)] for i in p.dims]

          #update plot options
          self.var.set(p.var); self.vars['values']=p.var
          sfm1=ttk.Frame(master=fm); sfm1.grid(row=0,column=0,sticky='W',pady=5)
          for n,[dn,ds,dvar] in enumerate(zip(p.dnames,p.dims,p.dvars)):
              if p.init==0:
                  if ds==1 or n==3:
                      dvar.set(0)
                  else:
                      dvar.set('all') if n==0 else dvar.set('all') if n==1 else dvar.set(0)#; dvar.set(p.xs[n][-1])
              vs=p.xs[n] if ds==1 else ['all','mean','min','max','sum',*p.xs[n]]
              dw=7 if n==2 else 2 if n==3 else 5
              sfm11=ttk.Frame(master=sfm1); sfm11.grid(row=0,column=n,sticky='W',pady=5)
              ttk.Label(master=sfm11,text='  '+dn).grid(row=0,column=0,sticky='W')
              mm=ttk.Combobox(sfm11,textvariable=dvar,values=vs,width=dw,); mm.grid(row=0,column=1,sticky='W')
              mm.bind("<<ComboboxSelected>>",lambda x: self.set_anim()); ap.dns.append(mm)

          #add limit panel
          sfm=ttk.Frame(master=fm); sfm.grid(row=1,column=0,sticky='W')
          ttk.Label(master=sfm,text='  limit').grid(row=0,column=0,sticky='W')
          ttk.Entry(sfm,textvariable=p.vmin,width=10).grid(row=0,column=1,sticky='W',padx=2)
          ttk.Entry(sfm,textvariable=p.vmax,width=10).grid(row=0,column=2,sticky='W')
          tk.Checkbutton(master=sfm,text='transpose',variable=p.transpose,onvalue=1,offvalue=0).grid(row=0,column=3,sticky='W')

          #add gd.plot for *_nu.nc
          if fname.endswith('_nu.nc') or fname=='uv3D.th.nc':
             sfm=ttk.Frame(master=fm); sfm.grid(row=2,column=0,sticky='W')
             ttk.Label(master=sfm,text='  ').grid(row=0,column=0,sticky='W')
             tk.Checkbutton(master=sfm,text='grid',variable=p.grid,onvalue=1,offvalue=0).grid(row=0,column=1,sticky='W')
             tk.Checkbutton(master=sfm,text='bnd',variable=p.bnd,onvalue=1,offvalue=0).grid(row=0,column=2,sticky='W')
             if fname.endswith('_nu.nc'):
                tk.Checkbutton(master=sfm,text='dispaly on mesh',variable=p.ctr,onvalue=1,offvalue=0).grid(row=0,column=3,sticky='W')
             else:
                ttk.Label(master=sfm,text='   zoom').grid(row=0,column=3,sticky='W')
                ttk.Entry(sfm,textvariable=p.scale,width=5).grid(row=0,column=4,sticky='W')
       elif self.fmt in [2,3]:  #hotstart.nc or source.nc
          if p.init==1 and option==1: #save parameter
             p0=zdata(); p0.dvars=[i.get() for i in p.dvars]; p0.vmin=p.vmin.get(); p0.vmax=p.vmax.get(); p0.scale=p.scale.get()
             p0.transpose=p.transpose.get(); p0.grid=p.grid.get(); p0.bnd=p.bnd.get(); p0.dims=[*p.dims]
          if p.init==0:
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.scale=tk.DoubleVar(wd); p.vmin.set(0); p.vmax.set(0); p.scale.set(1)
             p.transpose=tk.IntVar(wd); p.grid=tk.IntVar(wd); p.bnd=tk.IntVar(wd); p.transpose.set(0); p.grid.set(0); p.bnd.set(0)
             if self.fmt==3: p.sctr=tk.IntVar(wd); p.srat=tk.DoubleVar(wd); p.ctr=tk.IntVar(wd); p.id=tk.IntVar(wd); p.sctr.set(0); p.srat.set(1); p.ctr.set(0); p.id.set(0)
          if option==0: self.var.set(p.var); self.vars['values']=p.vars

          #update panel
          fid=self.fids[fname]; p.var=self.var.get(); p.cvar=fid[p.var]
          p.dims=p.cvar.shape; p.dnames=[*p.cvar.dimensions]
          if p.var in['su2','sv2']: p.dims=[*p.dims,1]; p.dnames=[*p.dnames,'uv']
          p.info='  dim={}'.format(p.dims); self.info.config(text=p.info)
          p.dvars=[tk.StringVar(wd) for i in p.dims];
          p.xs=[*p.dims] if self.fmt==2 else [[*arange(i)] for i in p.dims]

          sfm1=ttk.Frame(master=fm); sfm1.grid(row=0,column=0,sticky='W',pady=5)
          for n,[dn,ds,dvar] in enumerate(zip(p.dnames,p.dims,p.dvars)):
              if ds==1:
                 dvar.set(0); vs=[0]; dw=2
                 if dn=='uv': vs=['all','mean','min','max','sum',0]; dw=5
              elif self.fmt==3 and p.var=='source_elem':
                 vs=['all']; dvar.set('all'); dw=5
              elif dn in ['node', 'elem', 'side']:
                 vs=['all','mean','min','max','sum']; dvar.set('all'); dw=5
              else:
                 vs=['all','mean','min','max','sum',*arange(ds)]; dvar.set(0); dw=5
                 if dn=='nsources': dvar.set('all')
                 if dn=='nsinks': dvar.set('all')
              if dn.startswith('time_'): dn='time'
              if dn=='ntracers': dn='tracer'
              if dn=='nVert': dn='layer'
              if dn=='nsources': dn='source'
              if dn=='nsinks': dn='sink'
              if fname=='ICM_param.nc':
                 if not hasattr(self,'hgrid'): self.read_hgrid()
                 if ds==self.hgrid.ne or ds==self.hgrid.np:
                    dn='node/elem'; dvar.set('all')
                 else:
                    dn='d{}'.format(n)
              sfm11=ttk.Frame(master=sfm1); sfm11.grid(row=0,column=n,sticky='W',pady=5)
              ttk.Label(master=sfm11,text='  '+dn).grid(row=0,column=0,sticky='W'); p.dnames[n]=dn
              mm=ttk.Combobox(sfm11,textvariable=dvar,values=vs,width=dw,); mm.grid(row=0,column=1,sticky='W')
              mm.bind("<<ComboboxSelected>>",lambda x: self.set_anim()); ap.dns.append(mm)

          if self.fmt==2 and p.var in ['su2','sv2']:
             ttk.Label(master=sfm1,text='  zoom').grid(row=0,column=len(p.dims)+1,sticky='W')
             ttk.Entry(sfm1,textvariable=p.scale,width=5).grid(row=0,column=len(p.dims)+2,sticky='W',padx=2)

          #add limit panel
          sfm=ttk.Frame(master=fm); sfm.grid(row=1,column=0,sticky='W')
          ttk.Label(master=sfm,text='  limit').grid(row=0,column=0,sticky='W')
          ttk.Entry(sfm,textvariable=p.vmin,width=10).grid(row=0,column=1,sticky='W',padx=2)
          ttk.Entry(sfm,textvariable=p.vmax,width=10).grid(row=0,column=2,sticky='W')
          tk.Checkbutton(master=sfm,text='transpose',variable=p.transpose,onvalue=1,offvalue=0).grid(row=0,column=3,sticky='W')
          if self.fmt==2:
             tk.Checkbutton(master=sfm,text='grid',variable=p.grid,onvalue=1,offvalue=0).grid(row=0,column=4)
             tk.Checkbutton(master=sfm,text='bnd',variable=p.bnd,onvalue=1,offvalue=0).grid(row=0,column=5,sticky='W')

          #panel for source.nc
          if self.fmt==3:
             sfm=ttk.Frame(master=fm); sfm.grid(row=2,column=0,sticky='W',pady=5)
             tk.Checkbutton(master=sfm,text='scatter',variable=p.sctr,onvalue=1,offvalue=0).grid(row=0,column=0,sticky='W')
             ttk.Label(master=sfm,text='  zoom').grid(row=0,column=1,sticky='W')
             ttk.Entry(sfm,textvariable=p.srat,width=10).grid(row=0,column=2,sticky='W',padx=2)
             tk.Checkbutton(master=sfm,text='grid',variable=p.grid,onvalue=1,offvalue=0).grid(row=0,column=3)
             tk.Checkbutton(master=sfm,text='bnd',variable=p.bnd,onvalue=1,offvalue=0).grid(row=0,column=4,sticky='W')
             tk.Checkbutton(master=sfm,text='color',variable=p.ctr,onvalue=1,offvalue=0).grid(row=0,column=5,sticky='W')
             tk.Checkbutton(master=sfm,text='id',variable=p.id,onvalue=1,offvalue=0).grid(row=0,column=6,sticky='W')
             wd.geometry('400x185')

          #restore parameters if dims are the same
          if p.init==1 and option==1 and array_equal(array(p0.dims),array(p.dims)):
             [i.set(k) for i,k in zip(p.dvars,p0.dvars)]; p.vmin.set(p0.vmin); p.vmax.set(p0.vmax); p.scale.set(p0.scale)
             p.transpose.set(p0.transpose); p.grid.set(p0.grid); p.bnd.set(p0.bnd)
       elif self.fmt==4: #*.th file
          if p.init==0:
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.transpose=tk.IntVar(wd); p.vmin.set(0); p.vmax.set(0); p.transpose.set(0)
             p.dvars=[tk.StringVar(wd), tk.StringVar(wd)]; self.var.set(''); p.axs={'':[[],[]]}; p.var=''
          if option==1:
             self.convert_thfile(); p.dvars[0].set('all'); p.dvars[1].set(0)
             p.var=self.var.get(); p.dims=[*p.cvar.shape]; self.info.config(text=' dim={}'.format(p.dims))
          else:
             self.var.set(p.var)
          p.dnames=['time','dims']; p.xs=p.axs[p.var]

          #option panel
          self.vars['values']=self.thfiles
          sfm=ttk.Frame(master=fm); sfm.grid(row=0,column=0,sticky='W',pady=5)
          for m in arange(2):
              ttk.Label(master=sfm,text='  time' if m==0 else '  dims').grid(row=0,column=2*m,sticky='W')
              xs=['all','mean','min','max','sum',*arange(len(p.xs[m]))]
              mm=ttk.Combobox(sfm,textvariable=p.dvars[m],values=xs,width=7,); mm.grid(row=0,column=2*m+1,sticky='W')
              mm.bind("<<ComboboxSelected>>",lambda x: self.set_anim()); ap.dns.append(mm)

          #add limit panel
          sfm=ttk.Frame(master=fm); sfm.grid(row=1,column=0,sticky='W')
          ttk.Label(master=sfm,text='  limit').grid(row=0,column=0,sticky='W')
          ttk.Entry(sfm,textvariable=p.vmin,width=10).grid(row=0,column=1,sticky='W',padx=2)
          ttk.Entry(sfm,textvariable=p.vmax,width=10).grid(row=0,column=2,sticky='W')
          tk.Checkbutton(master=sfm,text='transpose',variable=p.transpose,onvalue=1,offvalue=0).grid(row=0,column=3,sticky='W')

       p.init=1; wd.update(); self.set_anim()
       self.init_plot(fmt=1)

   def set_anim(self):
       #set available dimensions for animation
       ap,p=self.ap,self.params[self.fname]
       if not hasattr(p,'dnames'): return
       dnames=p.dnames; dvars=p.dvars; avars=[]
       for dname,dvar,dn in zip(dnames,dvars,ap.dns):
           if dvar.get() in ['all','mean','min','max','sum'] or len(dn['values'])==1: continue
           avars.append(dname)
       ap.vars['values']=avars; ap.var.set(avars[0] if len(avars)!=0 else 'none')

   def anim(self,fmt=0):
       #animation for along one dimension
       ap,p=self.ap,self.params[self.fname]; dn=ap.var.get()
       if fmt==4: ap.stop=1; return
       if dn=='none': return

       #determine loop arange
       sid=p.dnames.index(dn); dvar=p.dvars[sid]; ds=ap.dns[sid]['values']  #find loop dimension
       nrec=len(ds); iv0=ds.index('0'); iv=ds.index(dvar.get()); ns=ap.skip.get() #get dim info.
       ap.stop=0 if fmt in [3,5] else 1  #used to stop the animation
       if fmt==1: ts=[ds[iv0],]   #first record
       if fmt==7: ts=[ds[-1], ]   #last  record
       if fmt==2: ts=[ds[max(iv-ns,iv0)],]  #back one
       if fmt==6: ts=[ds[min(iv+ns,nrec-1)],] #forward one
       if fmt==3: ts=ds[iv0:(iv+1)][::-1]  #loop backward
       if fmt==5: ts=ds[iv:-1]    #loop backward

       #loop plot in loop
       for ti in ts[::ns]:
           dvar.set(ti); self.plot(); pause(0.1)
           if ap.stop==1: break

   def init_plot(self,fmt=0):
       def _fig_close(hf): #mark closed figure
           hf.closed=True
       def _fig_resize(args): #update figsize
           hf,p=args; p.figsize=[hf.get_figwidth(),hf.get_figheight()]
       def _fig_xy(p): #update figure xm, and ym
           if hasattr(p,'hp'): p.xm=p.hp.get_xlim(); p.ym=p.hp.get_ylim()

       fname=self.fname; p=self.params[fname]; istat=0
       if fname in self.figs: #restore old figure
           hf=self.figs[fname]
           if hf.closed: #open a new figure canvas if original figure is closed
              if fmt==1: return
              hf=figure(figsize=p.figsize,num=hf.number); hf.closed=False; self.figs[fname]=hf; istat=1
           else:
              figure(hf.number)
       elif fmt==0 and (fname not in self.figs): #new figure
           hf=figure(len(self.figs)); self.figs[fname]=hf; hf.closed=False; _fig_resize([hf,p]); istat=1

       #define figure actions
       if istat==1:
           hf.canvas.mpl_connect("close_event",lambda x: _fig_close(hf))
           hf.canvas.mpl_connect("resize_event",lambda x: _fig_resize([hf,p]))
           hf.canvas.mpl_connect("draw_event",lambda x: _fig_xy(p))
       if fmt==0: return hf

   def plot(self):
       #plot for each input file
       if self.var.get()=='': return
       fname,run,params,fids=self.fname,self.run,self.params,self.fids; p=params[fname]
       self.read_input(); hf=self.init_plot(); hf.clf(); p.curve_1d=False
       if hasattr(p,'data') and (p.data is None):  close(hf); self.ap.stop=1; return

       #save axis limit for reset function
       def slimit(x,y,v=None):
           p.xm0=[array(x).min(),array(x).max()]; p.ym0=[array(y).min(),array(y).max()]
           if v is not None: p.vm0=[v.min(),v.max()]

       if self.fmt==0:  #gr3 files
          gd=self.hgrid; p=self.params[fname]; data=fids[fname]; pfmt=0; fxy=0
          if p.sflux.get()!='None':
             if not hasattr(gd,'lon'): gd0=read_schism_hgrid(self.run+'hgrid.ll'); gd.lon=gd0.x; gd.lat=gd0.y
             sid=read(self.run+'sflux'+os.sep+p.sflux.get()+'.nc',1); sx=array(sid['lon'][:]); sy=array(sid['lat'][:]); sid.close()
             for i,k in zip(sx,sy): plot(i,k,'-',color='orange',lw=0.5,alpha=1,zorder=1); fxy=1
             for i,k in zip(sx.T,sy.T): plot(i,k,'-',color='orange',lw=0.5,alpha=1,zorder=1)
             if abs(gd.lon-gd.x).mean()>1: pfmt=-1; slimit(gd.lon,gd.lat,data)
          if p.ctr.get()==1:  gd.plot(fmt=1,value=data,clim=[p.vmin.get(),p.vmax.get()],ticks=11,cmap='jet',method=1,xy=fxy,zorder=0)
          if p.grid.get()==1: gd.plot(xy=fxy,zorder=2)
          if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=1,xy=fxy,zorder=2)
          p.hp=gca(); slimit(gd.x,gd.y,data)
       elif ((self.fmt==2 and (p.var in ['su2','sv2'])) or (self.fmt==1 and fname=='uv3D.th.nc')) and p.data.ndim==2 and p.data.shape[1]==2: #vector (su2,sv2) in hotstart, or uv3d.th.nc
          if not hasattr(self,'hgrid'): self.read_hgrid()
          gd=self.hgrid
          if p.grid.get()==1: gd.plot()
          if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=1)
          if self.fmt==1:
             bid=hstack(gd.iobn); hv=quiver(gd.x[bid],gd.y[bid],p.data[:,0],p.data[:,1],scale=1/p.scale.get(),width=0.001,scale_units='inches'); p.hp=gca()
          elif self.fmt==2:
             if not hasattr(self.hgrid,'xcj'): self.hgrid.compute_side(fmt=2)
             if hasattr(p,'xm') and hasattr(p,'ym'):
                fpm=(gd.xcj>=p.xm[0])*(gd.xcj<=p.xm[1])*(gd.ycj>=p.ym[0])*(gd.ycj<=p.ym[1])
             else:
                fpm=~isnan(p.data[:,0])
             hv=quiver(gd.xcj[fpm],gd.ycj[fpm],p.data[fpm,0],p.data[fpm,1],scale=1/p.scale.get(),width=0.001,scale_units='inches'); p.hp=gca()
          quiverkey(hv,X=0.92, Y=1.01, U=1, label='1.0 m/s',color='r', labelpos='E',zorder=4)
          pfmt=7; slimit(gd.x,gd.y,p.data)
       elif self.fmt==3 and p.sctr.get()==1 and p.data.ndim==1 and p.data.size==p.dims[-1]: #source.nc
            vm=[p.vmin.get(),p.vmax.get()]
            if not hasattr(self,'hgrid'): self.read_hgrid()
            if not hasattr(self.hgrid,'dpe'): self.hgrid.compute_ctr()
            if not hasattr(p,'isc') and ('source_elem' in self.fids[fname]): p.isc=array(self.fids[fname]['source_elem'][:])-1
            if not hasattr(p,'isk') and ('sink_elem' in self.fids[fname]): p.isk=array(self.fids[fname]['sink_elem'][:])-1

            #prepare scatter data
            gd=self.hgrid; srat=p.srat.get()
            sind=p.isc if p.var in ['source_elem','msource','vsource'] else p.isk; xi,yi=gd.xctr[sind],gd.yctr[sind]; eid=arange(len(sind))
            data=ones(p.data.shape) if p.var in ['source_elem','sink_elem'] else p.data.copy()
            fpn=data!=-9999; xi,yi,data,eid=xi[fpn],yi[fpn],data[fpn],eid[fpn] #remove -9999 values
            if data.max()<=0: data=-data #plot negative values (vsink)
            fpn=data>0; xi,yi,data,eid=xi[fpn],yi[fpn],data[fpn],eid[fpn] #only keep data>0
            if data.size==0: print('no valid points found!'); close(hf); return

            #plot and label
            if p.ctr.get()==0:
               hg=scatter(xi,yi,s=data*srat,c='r')
            else:
               hg=scatter(xi,yi,s=srat*10,c=data)
            p.hp=gca(); slimit(gd.x,gd.y,data); pfmt=2
            if p.grid.get()==1: gd.plot()
            if p.bnd.get()==1:  gd.plot_bnd(c='k',lw=0.3)
            if p.id.get()==1: #plot source id
               if hasattr(p,'fmt')  and hasattr(p,'xm') and p.fmt==pfmt: fp=(xi>=p.xm[0])*(xi<=p.xm[1])*(yi>=p.ym[0])*(yi<=p.ym[1]); xi,yi,eid=xi[fp],yi[fp],eid[fp]
               for xii,yii,eidi in zip(xi,yi,eid): text(xii,yii,'{}'.format(eidi),fontsize=7)

            #legend
            if p.ctr.get()==0:
               v1,v2=data.min(),data.max();  m1,m2=int(log10(v1)),int(log10(v2))
               m1=max(0,m1) if m2>=0 else m2; ms=[i for i in arange(m1,m2+1) if (10.0**i>=v1) and (10.0**i<=v2)]
               if len(ms)==0:
                  hl=legend(*hg.legend_elements("sizes", num=[(v1+v2)/2])) #legend
               else:
                  hl=legend(*hg.legend_elements("sizes", num=[srat*10.0**i for i in ms])) #legend
                  for i,m in enumerate(ms): hl.texts[i].set_text('$10^{'+str(m)+'}$')  #set legend value
            else:
               cm.ScalarMappable.set_clim(hg,vmin=vm[0],vmax=vm[1])
               hc=colorbar(fraction=0.05,aspect=50,spacing='proportional',pad=0.02); hc.set_ticks(linspace(*vm,11)); hc.ax.set_ylim(vm)
       elif self.fmt in [1,2,3,4]: # bnd, nudge, hotstart, source.nc
          vm=[p.vmin.get(),p.vmax.get()]
          if p.data.ndim==1:
              i0=p.ax[0]; xi,xn=p.xs[i0], p.dnames[i0]
              if not hasattr(self,'hgrid'):self.read_hgrid()
              if fname.endswith('_nu.nc') and xn=='node' and p.ctr.get()==1:
                  gd=self.hgrid
                  if not hasattr(p,'sindn'): p.sindn=array(read(self.run+os.path.sep+fname,1).variables['map_to_global_node'])-1
                  vi=zeros(gd.np); vi[p.sindn]=p.data
                  gd.plot(fmt=1,value=vi,clim=[p.vmin.get(),p.vmax.get()],ticks=11,cmap='jet',method=1); p.hp=gca()
                  if p.grid.get()==1: gd.plot()
                  if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=1)
                  pfmt=6; slimit(gd.x,gd.y,p.data)
              elif self.fmt==2 and (xn in ['node', 'elem','node/elem','dim_{}'.format(self.hgrid.np),'dim_{}'.format(self.hgrid.ne)]): #schism grid plot
                  gd=self.hgrid
                  gd.plot(fmt=1,value=p.data,clim=[p.vmin.get(),p.vmax.get()],ticks=11,cmap='jet',method=1); p.hp=gca()
                  if p.grid.get()==1: gd.plot()
                  if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=1)
                  pfmt=3; slimit(gd.x,gd.y,p.data)
              else: #1D line
                  if self.fmt==2: xi=arange(xi)
                  plot(xi,p.data,'k-'); p.hp=gca()
                  xlabel(xn); pfmt=1; slimit(xi,p.data); p.curve_1d=True

                  #set ylim
                  y1,y2=p.vmin.get(),p.vmax.get()
                  if y1!=y2 and ~isnan(y1+y2):
                     p.ym=[y1,y2]
                  else:
                     p.ym=[p.data.min(),p.data.max()]
                  #if hasattr(p,'ym') and (p.ym[0]>=p.data.max() or p.ym[1]<=p.data.min()): p.ym=[p.data.min(),p.data.max()]
          else: #2D plots
              i1,i2=p.ax[:2]; xi,yi=p.xs[i1],p.xs[i2]; xn,yn=p.dnames[i1],p.dnames[i2]
              if self.fmt==2: xi=arange(xi); yi=arange(yi)
              if p.transpose.get()==1:
                 hg=contourf(yi,xi,p.data,vmin=vm[0],vmax=vm[1],levels=50)
                 xlabel(yn); ylabel(xn); pfmt=5; slimit(yi,xi,p.data)
              else:
                 hg=contourf(xi,yi,p.data.T,vmin=vm[0],vmax=vm[1],levels=50)
                 xlabel(xn); ylabel(yn); pfmt=4; slimit(xi,yi,p.data)
              cm.ScalarMappable.set_clim(hg,vmin=vm[0],vmax=vm[1]); p.hp=gca()
              hc=colorbar(fraction=0.05,aspect=50,spacing='proportional',pad=0.02); hc.set_ticks(linspace(*vm,11)); hc.ax.set_ylim(vm)

       #add title
       from matplotlib import rc
       hstr='{}'.format(fname) if len(self.vars['values'])==1 else '{}: {}'.format(fname,p.var); dstrs=[]
       for m,[dn,dvar] in enumerate(zip(p.dnames,p.dvars)):
           if dvar.get() in ['all','mean','min','max','sum']:
              dstrs.append(dn)
           else:
              dstrs.append('{}={}'.format(dn,dvar.get()))
       title('{} \n ({})'.format(hstr,', '.join(dstrs)),fontsize=12) if len(dstrs)!=0 else title(hstr,fontsize=12)

       self.reset_limit(1,fmt=pfmt)
       self.figs[fname]=hf; hf.tight_layout(); hf.canvas.draw(); hf.show()

   def reset_limit(self,method=0,fmt=0):
       fname=self.fname; p=self.params[fname]
       if method==1: #modify axes range
          if hasattr(p,'fmt') and p.fmt==fmt:
             if hasattr(p,'xm'): setp(p.hp,xlim=p.xm)
             if hasattr(p,'ym'): y1,y2=p.ym; y2=y1+max(y2-y1,1e-10); setp(p.hp,ylim=[y1,y2])
          else:
             p.fmt=int(fmt)
             if hasattr(p,'xm'): delattr(p,'xm')
             if hasattr(p,'ym'): delattr(p,'ym')
       elif method==0: #reset
            if hasattr(p,'xm0'): p.xm=p.xm0[:]
            if hasattr(p,'ym0'): p.ym=p.ym0[:]
            if hasattr(p,'vm0'): p.vmin.set(p.vm0[0]); p.vmax.set(p.vm0[1]); self.window.update()
            if hasattr(p,'ym0') and p.curve_1d: p.vmin.set(p.ym0[0]); p.vmax.set(p.ym0[1]); self.window.update()
            self.plot()

   def read_input(self):
       wd,fname=self.window,self.fname; fids,p=self.fids,self.params[fname]

       if self.fmt==0: #gr3 and prop files
          if fname not in fids:
              if fname=='hgrid.gr3':
                 run=self.run; fn1=run+'hgrid.npz'; fn2=run+'grid.npz'; fn3=run+'hgrid.gr3'; fexist=os.path.exists
                 gd=read(fn1) if fexist(fn1) else read(fn2).hgrid if fexist(fn2) else read(fn3); self.hgrid=gd; fids[fname]=gd.z
              if fname!='hgrid.gr3' and fname.split('.')[-1] in ['gr3','ll','ic']: gd=read(fname); fids[fname]=gd.z
              if fname.endswith('.prop'): fids[fname]=read(fname); self.read_hgrid()
              if not hasattr(self,'hgrid'): self.hgrid=gd
              z=fids[fname]; p.info='  np={}, ne={}, [{}, {}]'.format(self.hgrid.np,self.hgrid.ne,'{:15f}'.format(z.min()).strip(),'{:15f}'.format(z.max()).strip())
              if p.vmin.get()==0 and p.vmax.get()==0: p.vmin.set(z.min()); p.vmax.set(z.max())
              self.info.config(text=p.info); wd.update()
       elif self.fmt in [1,2,3,4]: #for bnd, nudge files, or hotstart
          #get dimension index, read data slice
          cvar=fids[fname] if self.fmt==1 else p.cvar
          dns=array([i.get() for i in p.dvars])
          if not hasattr(p,'dns0'): #return if parameters are the same as old ones
             p.dns0=dns; p.cvar0=cvar
          else:
             if cvar==p.cvar0 and array_equal(dns,p.dns0):
                return
             else:
                p.dns0=dns; p.cvar0=cvar
          if sum(dns=='all')>=3 or sum(dns=='all')==0 : print("can't plot scalar or 3D data"); return

          #read data
          dind=','.join([':' if (i in ['all','mean','min','max','sum']) else str(i) for i in dns])
          if p.var in ['su2','sv2']: dind=dind[:-2]
          if p.var in ['su2','sv2'] and dns[-1]!='0': #deal with vector in hotstart
             fid=fids[fname]; exec('p.data=array(concatenate((fid["su2"][{}][...,None],fid["sv2"][{}][...,None]),axis=-1))'.format(dind,dind))
          else:
             exec('p.data=array(cvar[{}])'.format(dind))
          if (p.vmin.get()==0 and p.vmax.get()==0) or (hasattr(p,'itr') and p.itr!=p.dvars[-1].get()) or isnan(p.vmin.get()+p.vmax.get()):
              p.vmin.set(p.data.min()); p.vmax.set(p.data.max())
          if self.fmt==1: p.itr=p.dvars[-1].get()
          p.info='  dim={}, [{}, {}]'.format(p.dims,'{:15f}'.format(p.data.min()).strip(),'{:15f}'.format(p.data.max()).strip())

          #sanity check
          errmsg=' found in "{}, {}[{}]" !!'.format(fname,p.var,dind)
          if sum(isnan(p.data))!=0: print('NaN value'+errmsg); p.data=None; return #check nan values
          if fname in ['source.nc','source_input']: #check source.nc
             if p.var=='msource' and sum(p.data[p.data!=-9999]<0)!=0: print('negative mass concentration'+errmsg); p.data=None; return
             if p.var=='vsource' and sum(p.data<0)!=0: print('negative volume source'+errmsg); p.data=None; return
             if p.var=='vsink' and sum(p.data>0)!=0: print('positive volume sink'+errmsg); p.data=None; return

          #perform operation, and get axis
          p.ax=[]; isht=0
          for n, dn in enumerate(dns):
              if dn=='all': p.ax.append(n)
              if dn in ['mean','min','max','sum']: exec('p.data=p.data.{}(axis={})'.format(dn,n-isht))
              if dn!='all': isht=isht+1
          #if p.data.ndim==2 and p.transpose.get()==1: p.ax=p.ax[::-1]; p.data=p.data.T
          if 'sum' in dns: p.info='  dim={}, [{}, {}]'.format(p.dims,'{:15f}'.format(p.data.min()).strip(),'{:15f}'.format(p.data.max()).strip())
          if not hasattr(p,'ax0'):
             p.ax0=p.ax
          else:
             if (not array_equal(array(p.ax0),array(p.ax))) and hasattr(p,'xm'): delattr(p,'xm'); delattr(p,'ym');
             p.ax0=p.ax
          self.info.config(text=p.info); wd.update()

   def read_input_info(self):
       '''
       read input information
       '''
       fname,run,params,fids=self.fname,self.run,self.params,self.fids; p=params[fname]

       #read input information
       if self.fmt==0: #gr3, prop files
          if hasattr(self,'hgrid'):
             p.info='  np={}, ne={}'.format(self.hgrid.np,self.hgrid.ne)
          elif fname.split('.')[-1] in ['gr3','ll','ic']:
             fid=open(run+fname); fid.readline(); ne,np=fid.readline().strip().split()[:2]
             p.info='  np={}, ne={}'.format(np,ne); fid.close()
          else:
             p.info=' '
       elif self.fmt==1: #  *D.th.nc, *_nu.nc files
           C=read(run+fname,1)
           p.var='time_series' if 'time_series' in [*C.variables] else 'tracer_concentration'
           fid=C.variables[p.var]; fids[fname]=fid
           p.dnames=['time','node','layer','tracer']; p.dims=fid.shape
           p.info='  dim={}'.format(p.dims)
       elif self.fmt==2: #hostart.nc
           cvar=read(run+fname,1).variables
           p.vars=[*cvar]; fids[fname]=cvar
           p.var='tr_el' if fname.startswith('hotstart') else p.vars[0] 
           p.info='  dim={}'.format(cvar[p.var].shape)
       elif self.fmt==3: #source.nc
           cvar=read(run+fname,1).variables if fname=='source.nc' else read(run+'.source.nc',1).variables
           p.vars=[i for i in cvar if (i in ['source_elem','sink_elem','vsource','vsink','msource'])]; fids[fname]=cvar; p.var='msource'
           p.info='  dim={}'.format(cvar[p.var].shape)
       elif self.fmt==4: #*.th
           if fname not in fids:
              sname='.th_{}.nc'.format(len([i for i in os.listdir(self.run) if i.startswith('.th_') and i.endswith('.nc')]))
              zdata().save(self.run+sname)
              fids[fname]=read(self.run+sname,1); p.info='convert *.th to '+sname

   def read_hgrid(self):
       if hasattr(self,'hgrid'): return self.hgrid
       run=self.run; fn1=run+'hgrid.npz'; fn2=run+'grid.npz'; fn3=run+'hgrid.gr3'; fexist=os.path.exists
       fnames=[i for i in self.fnames if (i.split('.')[-1] in ['gr3','ll','ic'])]
       if not (fexist(fn1) or fexist(fn2) or fexist(fn3) or len(fnames)!=0): sys.exit('can"t find hgrid.gr3')
       self.hgrid=read(fn1) if fexist(fn1) else read(fn2).hgrid if fexist(fn2) else read(fn3) if fexist(fn3) else read(fnames[0])
       return self.hgrid

   def convert_thfile(self):
       fname=self.fname; p,fid=self.params[fname],self.fids[fname]; p.var=self.var.get()
       if p.var in fid.variables: p.cvar=fid.variables[p.var]; return

       #add *.th to th.nc
       sname=fid.filepath(); fid.close() #save file information
       data=loadtxt(self.run+p.var); ti=data[:,0]/86400; data=data[:,1:] #read *.th
       ds=data.shape; p.axs[p.var]=[ti,arange(ds[1])]; #save dimension info
       fid=ReadNC(sname,1,mode='r+'); dims=[fid.dimensions[i].size for i in fid.dimensions]
       [fid.createDimension('dim_{}'.format(i),i) for i in ds if i not in dims] #add new dimension if necessar
       dnames=[*fid.dimensions]; dims=[fid.dimensions[i].size for i in dnames]; dn=[dnames[dims.index(i)] for i in ds] #get diminfo
       fid.createVariable(p.var,data.dtype,dn,zlib=True); fid.variables[p.var][:]=data; fid.close() #add variables
       fid=read(sname,1); p.cvar=fid.variables[p.var]; self.fids[fname]=fid

   def run_info(self,run):
       '''
       collect schism input information
       '''
       self.run=os.path.abspath(run)+os.path.sep; fexist=os.path.exists
       snames=os.listdir(self.run); fnames=[]
       [fnames.append(i) for i in snames if i=='hgrid.gr3']         #hgrid.gr3
       [fnames.append(i) for i in snames if i=='hotstart.nc']       #hotstart.nc
       [fnames.append(i) for i in snames if i=='ICM_param.nc']      #ICM_param.nc
       [fnames.append(i) for i in snames if i.endswith('D.th.nc')]  #3D bnd
       [fnames.append(i) for i in snames if i.endswith('_nu.nc')]   #3D nudge
       [fnames.append(i) for i in snames if i=='source.nc']         #source.nc
       [fnames.append('source_input') for i in snames if i=='source_sink.in']           #source_input
       [snames.remove(i) for i in snames if i in ['source_sink.in','vsource.th','vsink.th','msource.th']]  #remove source_input
       snames=[i for i in snames if not i.startswith('vsource.')]                                        #remove vsource.th
       fnames.extend(unique(['*.th' for i in snames if i.endswith('.th')]))          #*.th
       [fnames.append(i) for i in snames if i.endswith('.ll')]                       #hgrid.gr3
       [fnames.append(i) for i in snames if i.endswith('.gr3') and (i not in fnames)]#gr3
       [fnames.append(i) for i in snames if i.endswith('.prop')]                     #prop
       mc=[i for i in snames if i.endswith('.ic') and ('hvar' in i)]                 #ic files
       [fnames.append(i) for i in snames if i.endswith('.ic') and (i not in mc)]; fnames.extend(mc)   #ic
       [fnames.append(i) for i in snames if i.endswith('.nc') and (i not in ['.source.nc',*fnames]) and (not i.startswith('.th_'))]  #other nc files
       self.fnames=fnames; self.thfiles=[i for i in snames if i.endswith('.th')]     #*.th
       self.sflux=[i[:-3] for i in os.listdir(self.run+'sflux') if i.endswith('nc')] if fexist(self.run+'sflux') else None
       # self.StartT=0 #update this later

   def cmd_window(self):
        import tkinter as tk
        from tkinter import ttk
        cw=tk.Toplevel(self.window); cw.geometry("400x200"); cw.title('command input')
        cw.rowconfigure(0,minsize=150, weight=1); cw.columnconfigure(0,minsize=2, weight=1)
        txt=tk.Text(master=cw,width=150,height=14); txt.grid(row=0,column=0,pady=2,padx=2,sticky='nsew')
        rbn=ttk.Button(cw, text= "run",command=lambda: self.cmd_exec(txt.get('1.0',tk.END))); rbn.grid(row=1,column=0,padx=10)
        cw.update(); xm=max(txt.winfo_width(),rbn.winfo_width()); ym=txt.winfo_height()+rbn.winfo_height()+12
        if hasattr(self,'cmd'): txt.insert('1.0',self.cmd)
        cw.geometry('{}x{}'.format(xm,ym)); cw.update()
        print('control vars: [fname,p,fid,hf,ax,gd]')

   def cmd_exec(self,cmd):
        self.cmd=cmd
        fname=self.fname; p,fid=self.params[fname],self.fids[fname]
        if fname in self.figs: hf=self.figs[fname]; ax=gca()
        if hasattr(self,'hgrid'): gd=self.hgrid

        #run command
        for i in cmd.strip().split('\n'):
            if i=='': continue
            try:
               print('run: '+i); exec(i)
            except:
               print('fail: '+i)
        print('\n')
        if fname in self.figs: hf.canvas.draw()
        self.window.update()

   def window_exit(self):
       for fn in self.figs: close(self.figs[fn])
       if os.path.exists(self.run+'.source.nc'): os.remove(self.run+'.source.nc')
       if '*.th' in self.fids: os.remove(self.fids['*.th'].filepath())
       self.window.destroy()

   def init_window(self):
       import tkinter as tk
       from tkinter import ttk

       #init
       wd=tk.Tk(); self.window=wd
       wd.title("SCHISM check (Author: Z. WANG)")
       wd.rowconfigure([0,1,2], minsize=5, weight=3)
       wd.columnconfigure(0, minsize=50, weight=1)

       #Input file
       fm=ttk.Frame(master=wd);   fm.grid(row=0,column=0,sticky='NW',pady=0)
       sfm=ttk.Frame(master=fm);  sfm.grid(row=0,column=0,sticky='W',pady=3)
       self.input=tk.StringVar(wd);  self.input.set(self.fnames[0])
       ttk.Label(master=sfm,text='  Input ').grid(row=0,column=0,sticky='W',pady=0)
       svar=ttk.Combobox(sfm,textvariable=self.input,values=self.fnames,width=13,); svar.grid(row=0,column=1)
       svar.bind("<<ComboboxSelected>>",lambda x: self.update_panel())

       ttk.Label(master=sfm,text='     Variable ').grid(row=0,column=2,sticky='W',pady=0); self.var=tk.StringVar(wd)
       vars=ttk.Combobox(sfm,textvariable=self.var,values='',width=16,); vars.grid(row=0,column=3,sticky='E')
       vars.bind("<<ComboboxSelected>>",lambda x: self.update_panel(1)); self.vars=vars

       #input info
       sfm=ttk.Frame(master=fm);  sfm.grid(row=2,column=0,sticky='W',pady=0)
       self.info=ttk.Label(master=sfm,text='   '); self.info.grid(row=0,column=0,sticky='W',pady=0)

       #action options for each input files
       fm=ttk.Frame(master=wd); fm.grid(row=1,column=0,sticky='W',pady=0,padx=0); self.frame=fm

       #animation
       ap=zdata(); ap.var=tk.StringVar(wd); ap.skip=tk.IntVar(wd); ap.var.set('none'); ap.skip.set(1); self.ap=ap
       sfm=ttk.Frame(master=wd); sfm.grid(row=2,column=0,sticky='SW',pady=1,padx=1)
       ttk.Label(master=sfm,text='anim').grid(row=0,column=0,sticky='W')
       ap.vars=ttk.Combobox(sfm,textvariable=ap.var,values=['none',],width=5,); ap.vars.grid(row=0,column=1)
       ttk.Button(master=sfm,text='|<',    width=3,command=lambda: self.anim(1)).grid(row=0,column=2,sticky='W',padx=1)
       ttk.Button(master=sfm,text='<',     width=3,command=lambda: self.anim(2)).grid(row=0,column=3,sticky='W',padx=1)
       ttk.Button(master=sfm,text='\u25C0',width=3,command=lambda: self.anim(3)).grid(row=0,column=4,sticky='W',padx=1)
       ttk.Button(master=sfm,text='\u2551',width=3,command=lambda: self.anim(4)).grid(row=0,column=5,sticky='W',padx=1)
       ttk.Button(master=sfm,text='\u25B6',width=3,command=lambda: self.anim(5)).grid(row=0,column=6,sticky='W',padx=1)
       ttk.Button(master=sfm,text='>',     width=3,command=lambda: self.anim(6)).grid(row=0,column=7,sticky='W',padx=1)
       ttk.Button(master=sfm,text='>|',    width=3,command=lambda: self.anim(7)).grid(row=0,column=8,sticky='W',padx=1)
       ttk.Label(master=sfm,text=' skip:',width=6).grid(row=0,column=9,sticky='E')
       ttk.Entry(sfm,textvariable=ap.skip,width=3).grid(row=0,column=10,sticky='W')

       #draw and exit
       fm=ttk.Frame(master=wd); fm.grid(row=3,column=0,sticky='nesw',pady=1,padx=1)
       ttk.Button(master=fm,text='draw',width=6,command=lambda: self.plot()).pack(side=tk.RIGHT)
       ttk.Button(master=fm,text='exit',command=self.window_exit,width=5).pack(side=tk.LEFT,padx=1)
       ttk.Button(master=fm,text='command',width=8,command=self.cmd_window).pack(side=tk.LEFT,padx=1)
       ttk.Button(master=fm,text='reset', command=self.reset_limit,width=6).pack(side=tk.LEFT,padx=1)

       #resize window
       wd.protocol("WM_DELETE_WINDOW",self.window_exit)
       wd.geometry('400x175'); wd.update()

if __name__=="__main__":
    pass
