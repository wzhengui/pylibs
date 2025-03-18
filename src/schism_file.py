#/usr/bin/env python3
#Copyright 2021 Zhengui WANG
#Apache License, Version 2.0; http://www.apache.org/licenses/LICENSE-2.0

from pylib import *

class schism_grid(zdata):
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

    #------------------------------------------------------------------------------------------------
    #alias for grid properties
    #------------------------------------------------------------------------------------------------
    @property
    def backend(self):
        bkn=mpl.get_backend().lower(); backend=1 if (bkn in ['qt5agg','qtagg']) else 2 if (bkn in ['tkagg',]) else 0
        return backend

    #node alias
    @property
    def z(self):
        return self.dp
    @z.setter
    def z(self,value):
        self.dp=value
    @property
    def xy(self):
        return c_[self.x,self.y]
    @property
    def xyz(self):
        return c_[self.x,self.y,self.z]
    @property
    def cxy(self):
        return self.x+1j*self.y
    @property
    def lxy(self):
        return c_[self.lon,self.lat] if hasattr(self,'lon') else self.xy
    @property
    def lxyz(self):
        return c_[self.lon,self.lat,self.z] if hasattr(self,'lon') else self.xyz
    @property
    def xm(self):
        return [self.x.min(),self.x.max()]
    @property
    def ym(self):
        return [self.y.min(),self.y.max()]
    @property
    def zm(self):
        return [self.z.min(),self.z.max()]

    #element alias
    @property
    def xe(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return self.xctr
    @property
    def ye(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return self.yctr
    @property
    def ze(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return self.dpe
    @property
    def zctr(self):
        return self.ze
    @property
    def exy(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return c_[self.xctr,self.yctr]
    @property
    def exyz(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return c_[self.xctr,self.yctr,self.zctr]
    @property
    def ecxy(self):
        if not hasattr(self,'dpe'): self.compute_ctr()
        return self.xctr+1j*self.yctr
    @property
    def fp3(self):
        return self.i34==3
    @property
    def fp4(self):
        return self.i34==4

    #side alias
    @property
    def xs(self):
        if not hasattr(self,'dps'): self.compute_side(2)
        return self.xcj
    @property
    def ys(self):
        if not hasattr(self,'dps'): self.compute_side(2)
        return self.ycj
    @property
    def zs(self):
        if not hasattr(self,'dps'): self.compute_side(2)
        return self.dps
    @property
    def zcj(self):
        return self.zs
    @property
    def sxy(self):
        if not hasattr(self,'xcj'): self.compute_side(2)
        return c_[self.xcj,self.ycj]
    @property
    def sxyz(self):
        if not hasattr(self,'xcj'): self.compute_side(2)
        return c_[self.xcj,self.ycj,self.zcj]
    @property
    def scxy(self):
        if not hasattr(self,'xcj'): self.compute_side(2)
        return self.xcj+1j*self.ycj

    #wrap-around element
    @property
    def wrap(self):
        if not hasattr(self,'_wrap'):
           xm,ym=self.xm,self.ym; self._wrap=0
           if xm[1]-xm[0]<=360 and xm[0]>=-360 and xm[1]<=360 and ym[1]-ym[0]<=180 and ym[0]>=-90 and ym[0]<=90: self._wrap=1
        return self._wrap
    @wrap.setter
    def wrap(self,value):
        self._wrap=value
    #------------------------------------------------------------------------------------------------
    #alias for grid properties
    #------------------------------------------------------------------------------------------------

    def plot(self,fmt=0,value=None,ec=None,fc=None,lw=0.1,levels=None,shading='gouraud',xy=0,ticks=None,xlim=None,
             ylim=None,clim=None,extend='both',method=0,cb=True,cb_aspect=30,cb_pad=0.02,ax=None,mask=None,bnd=0,
             cmap='jet',wrap=None,dx_wrap=270,actions=True,**args):
        '''
        plot grid with default color value (grid depth)
        fmt=0: plot grid only; fmt=1: plot filled contours; fmt=3: plot bnd
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
        mask: change value to nan for plotting. 1). mask=number: fp=value==mask;  2). mask=[v1,v2]: fp=(value<v1)|(value>v2), 3). mask='<v1'; fp=value<v1 
        bnd=1: plots grid boundary; bnd=0: not plot
        wrap=1: plot wrap-around element seperately; wrap=0: not
        dx_wrap: the criteria for wrap-around element (Xmax-Xmin>dx_wrap)
        '''
        if (wrap is None): wrap=self.wrap
        if wrap==1: self.check_wrap_elem(dx_wrap=dx_wrap) #check whether to deal wrap-around elem.
        if hasattr(self,'pid') and (value is not None) and value.ndim==2: value=value[self.pid] #for schism transect
        if ec is None: ec='None'
        if fc is None: fc='None'
        if levels is None: levels=51
        if ax is None: ax=gca()
        if (mask is None) and ('nodata' in args): mask=args['nodata'] #for back compatiblity 
        x,y=[self.x,self.y] if xy==0 else [self.lon,self.lat] if xy==1 else xy.T
        fp3,fp4=self.fp3,self.fp4; vm=clim

        if fmt in [1,2]: #plot contours
           trs=r_[self.elnode[:,:3],c_[self.elnode[fp4,0],self.elnode[fp4,2:]]]
           if value is None: value=self.dp
           if mask is not None:
              if isnumber(mask): mask=float(mask)
              if isinstance(mask,list):
                 fpnd=(value<mask[0])|(value>mask[1]); mask=value[fpnd]
              elif isinstance(mask,str):
                 if mask[0]=='=' and mask[1]!='=': mask='='+mask
                 s=zdata(); exec('s.fp=value{}'.format(mask)); fpnd=s.fp
              else:
                 fpnd=value==mask
              value_nd=value[fpnd]; value[fpnd]=nan
           if vm is None: fpn=~isnan(value); vm=[min(value[fpn]),max(value[fpn])]
           if vm[0]==vm[1] or (vm[1]-vm[0])/(abs(vm[0])+abs(vm[1]))<1e-10: vm[1]=vm[1]+max([(vm[1]-vm[0])*1e-10,1e-10])

           #plot
           value0=value; npt=value.size #save original data
           if wrap==1: #move the location of wrap-around elements
              ie=pindex(self.fpg); ix=c_[ie,ie,ie].ravel(); iy=tile([0,1,2],len(ie))
              idg=cindex(c_[c_[ie,ie,ie].ravel(),tile([0,1,2],len(ie))],trs.shape); ip=trs.ravel()[idg]
              fp=x[ip]-x.min()>dx_wrap; trs.ravel()[idg[fp]]=arange(x.size,x.size+sum(fp)); x=r_[x,x[ip[fp]]-360]; y=r_[y,y[ip[fp]]]
              self.sindg=ip[fp]; value=r_[value,value[ip[fp]]] if npt==self.np else value
           if fmt==1 and method==1:  #tripcolor
              if npt==self.np: hg=tripcolor(x,y,trs,value,vmin=vm[0],vmax=vm[1],shading=shading,cmap=cmap,**args)
              if npt==self.ne: hg=tripcolor(x,y,trs,facecolors=r_[value,value[fp4]],vmin=vm[0],vmax=vm[1],cmap=cmap,**args)
              if npt==self.ne+sum(fp4) and sum(fp4)!=0: hg=tripcolor(x,y,trs,facecolors=value,vmin=vm[0],vmax=vm[1],cmap=cmap,**args)
           else:  #contourf or contour
              if npt==self.ne: value=self.interp_elem_to_node(value=value) #elem value to node value
              if sum(isnan(value))!=0: trs=trs[~isnan(value[trs].sum(axis=1))] #set mask
              if not hasattr(levels,'__len__'): levels=linspace(*vm,int(levels)) #detemine levels
              if fmt==1: hg=tricontourf(x,y,trs,value,levels=levels,vmin=vm[0],vmax=vm[1],extend=extend,cmap=cmap,**args)
              if fmt==2: hg=tricontour(x,y,trs,value,levels=levels, vmin=vm[0],vmax=vm[1],extend=extend,cmap=cmap,**args)
           if mask is not None: value0[fpnd]=value_nd
           self.data=value0

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
           #if not hasattr(lw,'__len__'): lw=[lw,lw*0.75]
           hg0=plot(*self.lines(wrap=wrap,dx_wrap=dx_wrap).T,lw=lw,color=ec[0],**args)
        if fmt==3 or bnd!=0: hb=self.plot_bnd(lw=lw)
        hg=hg0 if fmt==0 else hb if fmt==3 else hg if ec=='None' else [*hg0,hg]; self.hg=hg
        if xlim is not None: setp(ax,xlim=xlim)
        if ylim is not None: setp(ax,ylim=ylim)
        self.add_actions()
        return hg

        #-------------------------------------------------
        #for reference: old grid plot method
        #-------------------------------------------------
        #if gridplot==0:
        #   x4,y4=r_[self.xy[close_data_loop(self.elnode[fp4].T)],nan*ones([1,sum(fp4),2])].T; x4=x4.ravel(); y4=y4.ravel()
        #   if sum(fp3)!=0:
        #      hg0=[triplot(x,y,self.elnode[fp3,:3],lw=lw[0],color=ec[0],**args), plot(x4,y4,lw=lw[1],color=ec[1],**args)]
        #   else:
        #      hg0=[plot(x4,y4,lw=lw[1],color=ec[1],**args),]
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
        #      #ind=pindex(mask)
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

    def plot_bnd(self,c='k',lw=0.5,ax=None,xy=0,wrap=None,dx_wrap=270,**args):
        '''
          plot schims grid boundary
          xy=0: plot with gd.x,gd.y;  xy=1: use gd.lon,gd.lat;  xy=c_[x,y]: use provided xy coordinates
          gd.plot_bnd(): plot bnd
          gd.plot_bnd(c='rb'): open bnd in red,land bnd in blue
          wrap=1: plot wrap-around element seperately; wrap=0: not
          dx_wrap: the criteria for wrap-around element (Xmax-Xmin>dx_wrap)
        '''

        if (wrap is None): wrap=self.wrap
        if wrap==1: self.check_wrap_elem(dx_wrap=dx_wrap) #check whether to deal wrap-around elem.

        if ax!=None: sca(ax)
        x,y=[self.x,self.y] if xy==0 else [self.lon,self.lat] if xy==1 else xy.T
        if not hasattr(self,'nob'): self.compute_bnd()
        xy1,xy2=self.lines(1,wrap=wrap,dx_wrap=dx_wrap)
        if len(c)==1:
           hb=plot(*r_[xy1,xy2].T,c,lw=lw,**args); self.hb=hb
        else:
          hb1=plot(*xy1.T,c[0],lw=lw,**args); hb2=plot(*xy2.T,c[-1],lw=lw,**args); self.hb=[hb1,hb2]
        self.add_actions()
        return self.hb

    def lines(self,fmt=0,xy=0,wrap=0,dx_wrap=270):
        '''
        return lines in format of c_[x,y] for plotting purpose
          fmt=0: grid lines for triangles and quadlaterals
          fmt=1: lines for open and land bounaries 
          xy=0: gd.x,gd.y;  xy=1: use gd.lon,gd.lat xy=c_[x,y]: use provided xy coordinates
          wrap=0: exclude wrap-around elements; wrap=0: include wrap-around elements
        '''
        xy=self.xy if xy==0 else c_[self.lon,self.lat] if xy==1 else xy
        if fmt==0:
           fp3,fp4,fp3g,fp4g=[self.fp3,self.fp4,0,0] if wrap==0 else [self.fp3*(~self.fpg),self.fp4*(~self.fpg),self.fp3*self.fpg,self.fp4*self.fpg]
           xy3=r_[xy[close_data_loop(self.elnode[fp3,:3].T)],nan*ones([1,sum(fp3),2])].transpose([1,0,2]).reshape([5*sum(fp3),2])
           xy4=r_[xy[close_data_loop(self.elnode[fp4].T)],nan*ones([1,sum(fp4),2])].transpose([1,0,2]).reshape([6*sum(fp4),2])
           if wrap==1: #for wrap-around elem.
              xy3g=r_[xy[close_data_loop(self.elnode[fp3g,:3].T)],nan*ones([1,sum(fp3g),2])].transpose([1,0,2]).reshape([5*sum(fp3g),2])
              xy4g=r_[xy[close_data_loop(self.elnode[fp4g].T)],nan*ones([1,sum(fp4g),2])].transpose([1,0,2]).reshape([6*sum(fp4g),2])
              fp=((xy3g[:,0]-self.xm[0])>dx_wrap)*(~isnan(xy3g[:,0])); xy3g[fp,0]=xy3g[fp,0]-360
              fp=((xy4g[:,0]-self.xm[0])>dx_wrap)*(~isnan(xy4g[:,0])); xy4g[fp,0]=xy4g[fp,0]-360
           return r_[xy3,xy4] if wrap==0 else r_[xy3,xy4,xy3g,xy4g]
        elif fmt==1:
           xy1=[]; [xy1.extend(r_[nan*ones([1,2]),xy[i]]) for i in self.iobn if len(i)!=0] #open
           xy2=[]; [xy2.extend(r_[nan*ones([1,2]),xy[i if k==0 else r_[i,i[0]]]]) for i,k in zip(self.ilbn,self.island) if len(i)!=0]
           xy1=zeros([0,2]) if len(xy1)==0 else array(xy1); xy2=zeros([0,2]) if len(xy2)==0 else array(xy2)
           if wrap==1:
              x1,y1=xy1.T; x2,y2=xy2.T; n1=sum(abs(diff(x1))>dx_wrap); n2=sum(abs(diff(x2))>dx_wrap)
              x1=list(x1); y1=list(y1); x2=list(x2); y2=list(y2)
              for k in arange(n1): i=pindex(abs(diff(x1))>dx_wrap)[0]+1; f=sign(x1[i]-x1[i-1]); \
                                     x1.insert(i,x1[i]-f*360); y1.insert(i,y1[i]); x1.insert(i+1,nan); y1.insert(i+1,nan)
              for k in arange(n2): i=pindex(abs(diff(x2))>dx_wrap)[0]+1; f=sign(x2[i]-x2[i-1]); \
                                     x2.insert(i,x2[i]-f*360); y2.insert(i,y2[i]); x2.insert(i+1,nan); y2.insert(i+1,nan)
              xy1=c_[x1,y1]; xy2=c_[x2,y2]

           return xy1,xy2

    def add_actions(self):
        self.toolbar=gcf().canvas.toolbar
        if not hasattr(self,'backend'): self.get_backend()
        if self.backend==0: return
        if self.backend==1: ats=array([i.iconText() for i in self.toolbar.actions()])
        if self.backend==2: ats=[i[0] for i in self.toolbar.toolitems]
        if 'bp' not in ats: self.bp=schism_bpfile()
        if 'reg' not in ats: self.reg=schism_bpfile(fmt=1)
        if 'query' not in ats: self.query_pt()
        if 'bnd' not in ats: self.create_bnd()
        if 'node' not in ats: self.show_node_elem(fmt=0)
        if 'elem' not in ats: self.show_node_elem(fmt=1)
        if 'save' not in ats: self.save_item()

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
        if len(self.iobn)==1 or self.iobn.ndim==2: self.iobn=self.iobn.astype('int')

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
        if len(self.ilbn)==1 or self.ilbn.ndim==2: self.ilbn=self.ilbn.astype('int')

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
        fp3,fp4=self.fp3,self.fp4; dp=dp.transpose(idm); dpe=zeros(dms)
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
        if ((value is None) or fmt==1) and (not hasattr(self,'dpe')): self.compute_ctr()
        v0=self.dpe if (value is None) else value

        #interpolation
        vs=v0[self.ine]; fpn=(self.ine==-1)|isnan(vs); fpv=~fpn; vs[fpn]=0; v=zeros(self.np)*nan
        if fmt==0:
           if not hasattr(self,'area'): self.compute_area()
           a=self.area[self.ine]; ta=sum(a*fpv,axis=1); tv=sum(a*vs,axis=1)
           fp=ta!=0; v[fp]=tv[fp]/ta[fp]
        if fmt==1:
              dist=abs((self.xctr[self.ine]+1j*self.yctr[self.ine])-(self.x+1j*self.y)[:,None])
              w=1/(dist**p); w[fpn]=0; tw=w.sum(axis=1); fp=tw!=0; v[fp]=(w*vs).sum(axis=1)[fp]/tw[fp]
        if fmt==2: fp=sum(fpv,axis=1)!=0; vs[fpn]=-np.Inf; v[fp]=vs.max(axis=1)[fp]
        if fmt==3: fp=sum(fpv,axis=1)!=0; vs[fpn]=np.Inf;  v[fp]=vs.min(axis=1)[fp]
        if fmt==4:
           w=fpv; tw=w.sum(axis=1); fp=tw!=0; v[fp]=vs.sum(axis=1)[fp]/tw[fp]
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
           fp3,fp4=self.fp3,self.fp4; self.xctr,self.yctr,self.dpe=zeros([3,self.ne])
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
               if not hasattr(self,'indel'): self.compute_nne()
               ips=unique(self.elnode[self.indel[i]]); ips=setdiff1d(ips[1:] if ips[0]<0 else ips,i); nps=len(ips)
               if self.i34[self.indel[i]].max()==4: ips=array([*n1[n2==i],*n2[n1==i]]); nps=len(ips) #check quads
               A=angle((self.x[ips]-self.x[i])+1j*(self.y[ips]-self.y[i])); iA=argsort(A); A,ips=A[iA],ips[iA]
               cA=angle((self.xctr[self.indel[i]]-self.x[i])+1j*(self.yctr[self.indel[i]]-self.y[i])) #angle for element center

               #get all boundary nodes
               ib1=pindex(n1,i); ip1=n2[ib1];  ib2=pindex(n2,i); ip2=n1[ib2]; sinds[i]=[]
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
               ib=pindex(ifb,1)[0]; id0=ids[ib]; id=sinds[id0][-1]; ifb[sindf[id0]]=0; ibni=[id0,*sinds[id0]]; iloop=0

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
            i2=ones(nbs).astype('int'); fp3=pindex(self.i34[be],3); fp4=pindex(self.i34[be],4)
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
            for i in ibnx[nbnx>1]: k=pindex(isn[:,0],i); idx.extend(k); nbx=nbx+len(k)

            #search boundary from other nodes
            while(sum(ifb)!=0):
                #start points
                if nbx>0:
                    nbx=nbx-1; id0,id=isn[idx[nbx]]; ifb[idx[nbx]]=0
                else:
                    id0=isn[pindex(ifb,1)[0],0]; id=isn[sinds[id0],1]; ifb[sinds[id0]]=0
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
            px=self.x[ibn[i].astype('int')]; i0=pindex(px,px.min())[0]
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
               sind=pindex(sindb,m);  sids=[]; bf=ones(S.nbn[m])
               if len(sind)==0: continue
               for n,sindi in enumerate(sind): # add open bnd
                   id1=pindex(S.ibn[m],p1[sindi])[0]; id2=pindex(S.ibn[m],p2[sindi])[0]; sids.extend([id1,id2])
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
               sind=pindex(sindb,m)
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

    def compute_nee(self):
        '''
        compute element ball:
          mnee:  maximum number of elements in element ball (exclude itself)
          nee:   number of elements in element ball
          ielel: indices for each element ball
          iee:   indices for each element ball, but in maxtrix " shape=[ne,max(nee)]"
        '''
        if not hasattr(self,'nne'): self.compute_nne()
        nee=zeros(self.ne).astype('int'); iee=arange(self.ne).astype('int')[:,None]
        for m in arange(4):
            ies=self.ine[self.elnode[:,m]].T; ies[:,self.elnode[:,m]==-2]=-1
            for n in arange(self.mnei):
                ie=ies[n]; mnee=iee.shape[1]; fp=ie!=-1
                for k in arange(mnee): fp=fp*(ie!=iee[:,k])
                if sum(fp)==0: continue
                nee[fp]=nee[fp]+1 #add new elem
                if nee.max()>=mnee: iee=c_[iee,-ones(self.ne).astype('int')[:,None]]
                iee[pindex(fp),nee[fp]]=ie[fp]
        self.nee=nee; self.iee=iee[:,1:]; self.mnee=self.iee.shape[1]; self.ielel=array([k[:i] for i,k in zip(nee,self.iee)],dtype='O')
        return self.nee

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
            fps=pindex(self.i34>i); i34=self.i34[fps]; sinds=sind[fps]+i
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

    def compute_acor(self,pxy,fmt=0,out=1):
        '''
        compute acor coodinate for points pxy[npt,2]

        usage: ie,ip,acor=compute_acor(c_[xi,yi]), where xi and yi are array of coordinates
        outputs: ie[npt],ip[npt,3],acor[npt,3]
               ie:  the element number
               ip:  the nodal indices of the ie
               acor: the area coordinate
               fmt=0: (default) faster method by searching the neighbors of elements and nodes
               fmt=1: slower method using point-wise comparison
               out=0: use nearest node if points outside grid domain; out=1: interp using nearest boundary element

               Note: for interpolation of few pts on a large grid, fmt=1 can be faster than fmt=0
        '''

        npt=len(pxy); pip=-ones([npt,3]).astype('int'); pacor=zeros([npt,3])
        if fmt==0:
           pie=-ones(npt).astype('int'); sindp=arange(npt)
           #search element ctr
           if not hasattr(self,'xctr'): self.compute_ctr()
           #if hasattr(self,'bndinfo'): sindp=sindp[self.inside_grid(pxy)==1]
           sinde=near_pts(pxy[sindp],self.exy); fps,sip,sacor=self.inside_elem(pxy[sindp],sinde)
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
           sindp=pindex(pie,-1); sindp=sindp[self.inside_grid(pxy[sindp])==1]
           if len(sindp)!=0: pie[sindp],pip[sindp],pacor[sindp]=self.compute_acor(pxy[sindp],fmt=1)

           #interp on nearest side for points outside of domain
           if out==1:
              sindp=pindex(pie,-1)
              if sum(sindp)!=0:
                 if not hasattr(self,'xcj'): self.compute_side(fmt=2)
                 sinds=pindex(self.isdel[:,1],-1); pie[sindp]=self.isdel[sinds[near_pts(pxy[sindp],c_[self.xcj[sinds],self.ycj[sinds]])],0]
                 fps,sip,sacor=self.inside_elem(pxy[sindp],pie[sindp],out=1); pip[sindp[fps]]=sip; pacor[sindp[fps]]=sacor

        elif fmt==1:
            #check 1st triangle
            sindn=self.elnode.T[:3]; pie=inside_polygon(pxy,self.x[sindn],self.y[sindn],fmt=1)
            fps=pie!=-1; pip[fps]=sindn.T[pie[fps]]

            #check 2nd triangle
            sind4=pindex(self.i34,4); sind2=nindex(fps)
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
        if isinstance(kbp,int) or len(kbp)==1: return kbp
        #if fmt==0: kb=kbp[self.elnode]; kb[self.i34==3,-1]=-1; kb=kb.max(axis=1)
        if fmt==0: kb=kbp[self.elnode]; kb[self.i34==3,-1]=10000; kb=kb.min(axis=1)+1
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

    def interp(self,pxy,value=None,fmt=0,out=1,p=1):
        '''
        interpolate to get value at pxy
          pxy: 1). c_[x,y];   2).if pxy is 1D array(ne,np), interp_elem_to_node or interp_node_to_elem is used
          value=None: gd.dp is used; value: array of [np,] or [ne,]
          fmt=0: (default) faster method by searching the neighbors of elements and nodes
          fmt=1: slower method using point-wise comparison
          out=0: use nearest node if points outside grid domain; out=1: interp using nearest element
          p: only used for interp_node_to_elem when fmt==1

          Note: for interpolation of few pts on a large grid, fmt=1 can be faster than fmt=0
        '''
        pxy=array(pxy); ndim=pxy.ndim; np=len(pxy)
        if ndim==1 and np==self.np:
           return self.interp_node_to_elem(pxy)
        elif ndim==1 and np==self.ne:
           return self.interp_elem_to_node(pxy,fmt=fmt,p=p)
        else: #interp to pts
           vi=self.dp if value is None else value; npt=len(vi)  #get value
           pie,pip,pacor=self.compute_acor(pxy,fmt=fmt,out=out) #get interp coeff
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
           xyz=c_[self.exy,self.dpe]
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
            eid=sinde[pindex(abs(exy-pxy[0])<=(2*dist))] #elem. inside circle of 2*dc
            if eid.size==1: pv[pid]=value[eid]

            #find dist maxtrix, get weight and value, and assign average value at node
            ds=pid.size*eid.size; nloop=int(ceil(ds/ms)); dsb=int(ceil(pid.size/nloop))
            for n in arange(nloop):
                pidi=pid[arange((dsb*n),min([dsb*(n+1),pid.size])).astype('int')]
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
        zcor=self.compute_zcor(vgrid,eta=eta); ze=zeros([self.ne,vgrid.nvrt]); fp3,fp4=self.fp3,self.fp4
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
        self.cfl=0.5*dt*(abs(u)+sqrt(9.81*self.dpe))/sqrt(self.area/pi)
        self.cfl[self.fp4]=self.cfl[self.fp4]*sqrt(2)
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
        trs=r_[self.elnode[:,:3],c_[self.elnode[self.fp4,0],self.elnode[self.fp4,2:]]]
        hf=figure(); hf.set_visible(False)
        P=tricontour(self.x,self.y,trs,value,levels=levels); close(hf); cxy=[]
        for k in arange(len(P.collections)):
            p=P.collections[k].get_paths()
            for i in arange(len(p)):
                xii,yii=p[i].vertices.T
                xi=r_[xii,nan] if i==0 else r_[xi,xii,nan]
                yi=r_[yii,nan] if i==0 else r_[yi,yii,nan]
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
        qind=pindex(self.i34,4); A=self.angles[qind]; fp=qind==-999
        if angle_ratio is not None: fp=fp|((A.min(axis=1)/A.max(axis=1))<angle_ratio)
        if angle_min is not None: fp=fp|(A.min(axis=1)<=angle_min)
        if angle_max is not None: fp=fp|(A.max(axis=1)>=angle_max)
        self.index_bad_quad=qind[pindex(fp)]

        #output bad_quad location as bp file
        sbp=schism_bpfile(); sbp.x,sbp.y=self.exy[self.index_bad_quad].T; sbp.save(fname)
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
        sid=pindex(self.i34,4); ne,nea=self.ne,len(sid); self.elnode=resize(self.elnode,[ne+nea,4])
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
           sindw=pindex(srat>threshold)

        if angle_min is not None: #check minimum angle
           if not hasattr(self,'angles'): self.compute_angle()
           a=self.angles; a[a<0]=999; sindw=pindex(a.min(axis=1)<=angle_min)

        #combine and save
        sindw=array(sindw)
        if fname is not None: C=schism_bpfile(); C.x,C.y,C.z=self.xctr[sindw],self.yctr[sindw],self.dpe[sindw]; C.save(fname)
        return sindw

    def check_wrap_elem(self,fmt=0,dx_wrap=270):
        '''
        return indice from wrap-around elements with element width>dx_wrap
        fmt=0: return element indices; fmt=1: return boolean value; fmt=2: return both
        '''
        if not hasattr(self,'fpg'): x=self.x[self.elnode]; fp=self.elnode[:,-1]==-2; x[fp,-1]=x[fp,0]; self.fpg=(x.max(1)-x.min(1))>dx_wrap
        return pindex(self.fpg) if fmt==0 else self.fpg if fmt==1 else [pindex(self.fpg),self.fpg]

    def inside_elem(self,pxy,ie,out=0):
        '''
        check whether pts are inside elements, then compute area coordinates for pts in elements
           pxy: c_[x,y]
           ie: array of element indices corresponding to each pt
           out=0: omit points outside of elements: out=1: include points even if they are outside of elements
        '''
        sind=[]; pip=[]; pacor=[]; nloop=self.i34[ie].max()-2; fp4=self.i34[ie]==4
        for i in arange(nloop):
            #get pts and element info
            sindp=arange(len(pxy)) if i==0 else pindex((~fps)*fp4 if out==0 else fp4)
            if len(sindp)==0: continue
            ip=self.elnode[ie[sindp]][:,array([0,1,2] if i==0 else [0,2,3])]
            x1,x2,x3=self.x[ip].T; y1,y2,y3=self.y[ip].T; xi,yi=pxy[sindp].T

            #compute area coordinates
            A1=signa(c_[xi,x2,x3],c_[yi,y2,y3]); A2=signa(c_[x1,xi,x3],c_[y1,yi,y3]); A3=signa(c_[x1,x2,xi],c_[y1,y2,yi])
            if out==0:
               A0=signa(c_[x1,x2,x3],c_[y1,y2,y3]); fps=(A1>=0)*(A2>=0)*(A3>=0)
            else: #include points outside
               A1=abs(A1); A2=abs(A2); A3=abs(A3); A0=A1+A2+A3; fps=A0>=0
               if i==0: Am=c_[A1,A2,A3][fp4].min(axis=1)
            ac1=A1[fps]/A0[fps]; ac2=A2[fps]/A0[fps]; ac3=1-ac1-ac2
            fpn=ac3<0; ac2[fpn]=1-ac1[fpn]; ac3[fpn]=0
            #if not isinstance(fps,np.ndarray): fps=array([fps])

            #save results
            if i==1 and out==1:
               pip=array(pip); pacor=array(pacor)
               fpn=c_[A1,A2,A3].min(axis=1)<Am; pip[sindp[fpn]]=ip[fpn]; pacor[sindp[fpn]]=array(c_[ac1,ac2,ac3][fpn])
            else:
               sind.extend(sindp[fps]); pip.extend(ip[fps]); pacor.extend(c_[ac1,ac2,ac3])
        fps=argsort(sind)
        return array(sind)[fps],array(pip)[fps],array(pacor)[fps]

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

    def subgrid(self,nsub=2,reg=None):
        '''
        return a subgrid with elements refined
           nsub: number of segment each side divided; each element will be divided to nsub*nsub parts
           reg:  only refine element inside region;  1). string; 2). xy of region.  3). array of element indices
        '''
        if not hasattr(self,'dpe'): self.compute_ctr()
        if not hasattr(self,'elside'): self.compute_ic3()
        if reg is None:
           sinde=arange(self.ne)
        elif isinstance(reg,str):
           sinde=read(reg).inside(self.exy)
        elif isinstance(reg,np.ndarray):
           if reg.ndim==1: sinde=reg
           if reg.ndim==2: sinde=pindex(inside_polygon(self.exy,*rxy.T),1)

        nsb=nsub; nsp=nsub+1; xn=[]; yn=[]; zn=[]; iep=[]; npa=self.np #additional nodes
        #get the side information, and then divide them into nsb sections
        sinds=unique(self.elside[sinde]); sinds=sinds[sinds!=-1] #triangles, and all sides
        ne,ns=len(sinde),len(sinds); sdict=-ones(self.ns,'int'); sdict[sinds]=arange(ns) #side dict
        ds=[nsp,ns]; xs,ys,zs=zeros(ds),zeros(ds),zeros(ds); isd=zeros(ds,'int') #the sub-node information
        sindp=self.isidenode[sinds].T; x1,x2=self.x[sindp]; y1,y2=self.y[sindp]; z1,z2=self.z[sindp]
        xs[0]=x1; xs[-1]=x2; ys[0]=y1; ys[-1]=y2; isd[0]=sindp[0]; isd[-1]=sindp[1]
        for i in arange(1,nsb): xs[i]=((nsb-i)*x1+i*x2)/nsb; ys[i]=((nsb-i)*y1+i*y2)/nsb; zs[i]=((nsb-i)*z1+i*z2)/nsb
        isd[1:nsb]=arange(ns*(nsb-1),dtype='int').reshape([nsb-1,ns])+npa;
        xn.extend(xs[1:nsb].ravel()); yn.extend(ys[1:nsb].ravel()); zn.extend(zs[1:nsb].ravel());  npa=npa+ns*(nsb-1)

        #form matrix for each element to store (x,y,index)
        ds=[ne,nsp,nsp]; inode=zeros(ds,'int'); xb=zeros(ds); yb=zeros(ds); i34=self.i34[sinde]; fp3=i34==3; fp4=~fp3; ipp=[]
        for i in arange(4): #create outer part of matrix
            sid=self.elside[sinde,i%i34]; m=sdict[sid]; xi,yi,isi=xs[:,m],ys[:,m],isd[:,m]
            fp=self.elnode[sinde,(i+1)%i34]==self.isidenode[sid,1]; xi[:,fp]=xi[::-1,fp]; yi[:,fp]=yi[::-1,fp]; isi[:,fp]=isi[::-1,fp]
            for k in arange(nsb):
                if i<=2: #for triangles
                    i1,i2=[k,nsb-k] if i==0 else [nsb-k,0] if i==1 else [0,k]; fp=fp3
                    inode[fp,i1,i2]=isi[k,fp]; xb[fp,i1,i2]=xi[k,fp]; yb[fp,i1,i2]=yi[k,fp]
                i1,i2=[k,nsb] if i==0 else [nsb,nsb-k] if i==1 else [nsb-k,0] if i==2 else [0,k]; fp=fp4
                inode[fp,i1,i2]=isi[k,fp]; xb[fp,i1,i2]=xi[k,fp]; yb[fp,i1,i2]=yi[k,fp]
        for i in arange(1,nsb):   #fill the maxtrix
            for k in arange(2):
                [fp,dx]=[fp3,i] if k==0 else [fp4,0]; ix=arange(1,nsb-dx); ixa=ix[None,:]; nx=len(ix)
                xbi=(xb[fp,i,0][:,None]*(nsb-dx-ixa)+xb[fp,i,nsb-dx][:,None]*ixa)/(nsb-dx)
                ybi=(yb[fp,i,0][:,None]*(nsb-dx-ixa)+yb[fp,i,nsb-dx][:,None]*ixa)/(nsb-dx)
                inodei=arange(sum(fp)*nx,dtype='int').reshape([sum(fp),nx])+npa
                for n,ixi in enumerate(ix): xb[fp,i,ixi]=xbi[:,n]; yb[fp,i,ixi]=ybi[:,n]; inode[fp,i,ixi]=inodei[:,n]
                xn.extend(xbi.ravel()); yn.extend(ybi.ravel()); ipp.extend(tile(sinde[fp],[nx,1]).T.ravel()); npa=npa+sum(fp)*nx
        if len(ipp)!=0: #get the new node depth
           npp=len(ipp); fps,sip,sacor=self.inside_elem(c_[array(xn)[-npp:],array(yn)[-npp:]],array(ipp))
           zn.extend((self.dp[sip]*sacor).sum(axis=1)[argsort(fps)])

        #collect new element
        elnode=-2*ones([ne,nsb*nsb,4],'int'); ie=0; iepi=[]
        for i in arange(nsb): #for triangles
            m=arange(nsb-i,dtype='int'); n=i*ones(len(m),'int')
            i1=c_[r_[n,n[1:]],r_[n,n[1:]+1],r_[n+1,n[1:]+1]]; i2=c_[r_[m,m[1:]], r_[m+1,m[1:]],r_[m,m[:-1]]]; inodei=inode[:,i1,i2][fp3]
            for n in arange(len(i1)): elnode[:sum(fp3),ie,:3]=inodei[:,n]; ie=ie+1; iepi.append(sinde[fp3])
        iep.extend(array(iepi).T.ravel())
        #for quads
        elnode[sum(fp3):]=c_[inode[fp4,:nsb,:nsb][...,None],inode[fp4,:nsb,1:nsp][...,None],inode[fp4,1:nsp,1:nsp][...,None],inode[fp4,1:nsp,:nsb][...,None]].reshape([sum(fp4),nsb*nsb,4])
        iep.extend(tile(sinde[fp4],[nsb*nsb,1]).T.ravel())

        #prepare for new grid
        xn=r_[self.x,xn]; yn=r_[self.y,yn]; zn=r_[self.z,zn]; elnode=elnode.reshape([prod(elnode.shape[:2]),4])
        tinde=setdiff1d(self.isdel[sinds].ravel(),sinde); tinde=tinde[tinde!=-1] #elements in transition zone
        fps=self.elside[tinde]; tinds=sdict[fps]; tinds[fps==-1]=-1 #sides in transition zone

        #for transition zone, find triangles (rinde,rinds) with only 1 side divided
        fpe=(self.i34[tinde]==3)*(sum(tinds!=-1,axis=1)==1); rinde=tinde[fpe]
        rinds=sort(tinds[fpe],axis=1)[:,-1]; isn=isd[:,rinds]; in0=zeros(len(rinde),'int')
        for i in arange(3): #find all nodes
            n1=self.elnode[rinde,(i+1)%3]; n2=self.elnode[rinde,(i+2)%3]
            fp1=(n1==isn[0])*(n2==isn[-1]); fp2=(n1==isn[-1])*(n2==isn[0]); fp=fp1|fp2
            in0[fp]=self.elnode[rinde[fp],i]; isn[:,fp2]=flipud(flipud(isn[:,fp2]))
        for i in arange(nsb): elnode=r_[elnode,c_[in0,isn[i:(i+2),:].T,-2*ones([len(rinde),1],'int')]]; iep.extend(rinde)  #add divided triangles

        #for transition zone, with multipe side divided
        fp=((self.i34[tinde]==3)*(sum(tinds!=-1,axis=1)!=1))|(self.i34[tinde]==4); pinde=tinde[fp]; pinds=tinds[fp]; trs=[]
        for ie,ps in zip(pinde,pinds):
            isn=isd[:,ps[ps!=-1]]; ips=r_[self.elnode[ie,:self.i34[ie]],isn[1:-1].ravel()]
            T=mpl.tri.Triangulation(xn[ips],yn[ips]); tr=ips[unique(T.triangles,axis=0)]; flag=ones(len(tr))
            for isni in isn.T: #remove invalid triangles
                sd=dict(zip(ips,zeros(len(ips),'int')))
                for i in isni: sd[i]=1
                flag[array([sd[i] for i in tr.ravel()]).reshape(tr.shape).sum(axis=1)==3]=0
            fp=flag==1; trs.extend(array(tr)[fp]); iep.extend(tile(ie,sum(fp)))
        elnode=r_[elnode,resize(array(trs,'int'),[len(trs),4],-2)]

        #construct a new grid
        gdn=schism_grid(); gdn.np,gdn.x,gdn.y,gdn.dp=npa,xn,yn,zn; dinde=setdiff1d(arange(self.ne),r_[sinde,tinde])
        gdn.elnode=r_[self.elnode[dinde],elnode]; gdn.ne=len(gdn.elnode); gdn.i34=sum(gdn.elnode!=-2,axis=1); gdn.iep=r_[dinde,iep]
        return gdn

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
            if self.backend==1:
               ac=[i for i in self.toolbar.actions() if i.iconText()=='Pan'][0]
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
                  bid=S.bid[-1]; pid=pindex(S.ibn[bid],S.pt[-1])[0]
                  S.nlb=S.nlb+1; ibni=r_[S.ibn[bid][pid:],S.ibn[bid][0]]; S.ilbn.append(ibni)
                  hlb=plot(self.x[ibni],self.y[ibni],'g-'); S.hb.append(hlb); S.ihb.append(0)

               #save boundary information
               self.nob=S.nob; self.iobn=array(S.iobn,dtype='O'); self.nobn=array([len(i) for i in self.iobn])
               sid=setdiff1d(unique(S.sind),unique(array(S.bid)))
               self.nlb=S.nlb+len(sid); self.ilbn=array([*S.ilbn,*[S.ibn[i] for i in sid]],dtype='O')
               self.nlbn=array([len(i) for i in self.ilbn]); self.island=r_[tile(0,S.nlb),tile(1,len(sid))]

               #finish
               gcf().canvas.mpl_disconnect(self.cidbnd)
               if self.backend==1:
                  ac=[i for i in self.toolbar.actions() if i.iconText()=='Pan'][0]
                  if ac.isChecked(): ac.trigger()
               gcf().canvas.draw()

        def add_pt(x,y):
            distp=squeeze(abs((S.x-x)+1j*(S.y-y))); sid=pindex(distp,distp.min())[0]
            ip=S.ip[sid]; bid=S.sind[sid]; pid=pindex(S.ibn[bid],ip)[0]
            if S.npt!=0:
               bid0=S.bid[-1]; pid0=pindex(S.ibn[bid0],S.pt[-1])[0]
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
            bid=S.bid[-1]; pid=pindex(S.ibn[bid],S.pt[-1])[0]

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

        #add bndinfo capsule
        if not hasattr(self,'bndinfo'): self.bndinfo=zdata()
        S=self.bndinfo; S.hp=[]; S.hb=[]; S.ihb=[]; S.nob=0; S.iobn=[]; S.nlb=0; S.ilbn=[]; S.npt=0; S.pt=[]; S.bid=[]

        #add bnd icon
        if mpl._pylab_helpers.Gcf.get_active() is None: self.plot()
        if self.backend==1:
           acs=self.toolbar.actions(); ats=array([i.iconText() for i in acs])
           abn=acs[pindex(ats,'bnd')[0]] if 'bnd' in ats else self.toolbar.addAction('bnd')
           abn.triggered.connect(connect_actions) #connect to actions
        elif self.backend==2:
            acs=gcf().canvas.toolbar.toolitems; ats=[i[0] for i in acs]
            if 'bnd' in ats:
               ac=acs[ats.index('bnd')][1]
            else:
               import tkinter as tk
               ac=tk.Button(gcf().canvas.toolbar,text='bnd', command=connect_actions); ac.pack(side=tk.LEFT)
               gcf().canvas.toolbar.toolitems=(*gcf().canvas.toolbar.toolitems,('bnd',ac))
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
               #acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]);ac=acs[pindex(ats,'bp')[0]]
               if self.backend==1: ac=[i for i in self.toolbar.actions() if i.iconText()=='bp'][0]
               if self.backend==2: ac=[i[1] for i in self.toolbar.toolitems if i[0]=='bp'][0]
               if hasattr(ac,'bp'):
                  if ac.bp.nsta==0: return
                  distp=squeeze(abs((ac.bp.x-bx)+1j*(ac.bp.y-by))); sid=pindex(distp,distp.min())[0]
                  px=ac.bp.x[sid]; py=ac.bp.y[sid]; pz=1
            elif dlk==0 and btn==1:
               px=bx; py=by; pz=1
            elif dlk==0 and btn==2:
               if not hasattr(self,'hqt'): return
               self.hqt.remove(); self.hqp.remove(); delattr(self,'hqt'); delattr(self,'hqp'); gcf().canvas.draw()
               gcf().canvas.mpl_disconnect(self.cidquery)

            #annotate text
            if dlk==0 and (btn in [1,3]) and pz==1:
               pz=self.interp(c_[px,py],value=self.data if hasattr(self,'data') else self.dp)
               if not hasattr(self,'hqp'):
                  self.hqp=plot(px,py,'r^',ms=6,alpha=1)[0]
                  self.hqt=text(px,py,'',color='orangered',fontsize=12,bbox=dict(facecolor='w',alpha=0.5))
               self.hqt.set_x(px); self.hqt.set_y(py); self.hqt.set_text(' {}'.format(pz[0]))
               self.hqp.set_xdata([px]); self.hqp.set_ydata([py]); gcf().canvas.draw()

        if hasattr(self,'hqt'): delattr(self,'hqt')
        if hasattr(self,'hqp'): delattr(self,'hqp')
        if self.backend==1:
           acs=self.toolbar.actions(); ats=array([i.iconText() for i in acs])
           abp=acs[pindex(ats,'query')[0]] if 'query' in ats else self.toolbar.addAction('query')
           abp.triggered.connect(connect_actions)
        elif self.backend==2:
            acs=gcf().canvas.toolbar.toolitems; ats=[i[0] for i in acs]
            if 'query' in ats:
               ac=acs[ats.index('query')][1]
            else:
               import tkinter as tk
               ac=tk.Button(gcf().canvas.toolbar,text='query', command=connect_actions); ac.pack(side=tk.LEFT)
               gcf().canvas.toolbar.toolitems=(*gcf().canvas.toolbar.toolitems,('query',ac))

    def show_node_elem(self,fmt):
        '''
        show node/element number
        '''
        def _show_node_elem():
            xm=xlim(); ym=ylim(); hts=self.hts_p if fmt==0 else self.hts_e
            if len(hts)==0:
               if fmt==0:
                  sind=pindex((self.x>=xm[0])*(self.x<=xm[1])*(self.y>=ym[0])*(self.y<=ym[1]))
                  for i in sind: ht=text(self.x[i],self.y[i],'{}'.format(i+1),fontsize=6); hts.append(ht)
               else:
                  sind=pindex((self.xe>=xm[0])*(self.xe<=xm[1])*(self.ye>=ym[0])*(self.ye<=ym[1]))
                  for i in sind: ht=text(self.xctr[i],self.yctr[i],'{}'.format(i+1),fontsize=6); hts.append(ht)
            else:
               for i in arange(len(hts)): hts.pop().remove()
            gcf().canvas.draw()

        at='node' if fmt==0 else 'elem'
        if self.backend==1:
           acs=self.toolbar.actions(); ats=array([i.iconText() for i in acs])
           abp=acs[pindex(ats,at)[0]] if at in ats else gcf().canvas.toolbar.addAction(at)
           abp.triggered.connect(_show_node_elem)
        else:
           acs=self.toolbar.toolitems; ats=[i[0] for i in acs]
           if at in ats:
              ac=acs[ats.index(at)][1]
           else:
              import tkinter as tk
              ac=tk.Button(gcf().canvas.toolbar,text=at, command=_show_node_elem); ac.pack(side=tk.LEFT)
              gcf().canvas.toolbar.toolitems=(*gcf().canvas.toolbar.toolitems,(at,ac))
        if fmt==0: self.hts_p=[]
        if fmt==1: self.hts_e=[]

    def save_item(self):
        '''
        add method to save bp/reg file
        '''
        def _save_item(gd,hf):
            import tkinter as tk
            from tkinter import ttk,filedialog

            def chdirpath():
                sdir.set(filedialog.askdirectory(initialdir=sdir.get(), title = "choose save dir")); resize()
            def resize():
                ns=len(sdir.get()); pdir.config(width=max(50,ns+3)); wd.geometry('{}x120'.format(max(450,int(8.5*ns)))); wd.update()
            def savefile():
                ftype=file.get(); fn=sdir.get()+os.path.sep+fname.get()
                if ftype==0: gd.bp.save(fn)
                if ftype==1: gd.reg.write_reg(fn)
                if ftype==2: gd.save(fn,fmt=1)
                if ftype==3: gd.write_bnd(fn)
                if ftype==4: savefig(fn,hf)

            wd=tk.Tk(); wd.title('save schism files')
            file=tk.IntVar(wd); sdir=tk.StringVar(wd); fname=tk.StringVar(wd)

            fm=ttk.Frame(master=wd); fm.grid(row=0,column=0,sticky='NW',pady=6)
            ttk.Radiobutton(fm,text='bp',   variable=file, value=0).grid(row=0,column=0,sticky='W',padx=2)
            ttk.Radiobutton(fm,text='reg',  variable=file, value=1).grid(row=0,column=1,sticky='W',padx=2)
            ttk.Radiobutton(fm,text='hgrid',variable=file, value=2).grid(row=0,column=2,sticky='W',padx=2)
            ttk.Radiobutton(fm,text='bnd',  variable=file, value=3).grid(row=0,column=3,sticky='W',padx=2)
            ttk.Radiobutton(fm,text='fig',  variable=file, value=4).grid(row=0,column=4,sticky='W',padx=2)

            fm=ttk.Frame(master=wd); fm.grid(row=1,column=0,sticky='NW',pady=6)
            ttk.Label(fm,text='dir',width=3).grid(row=0,column=0,sticky='W',padx=1)
            pdir=ttk.Entry(fm,textvariable=sdir,width=50); pdir.grid(row=0,column=1,sticky='W')
            ttk.Button(fm,text='..',command=chdirpath,width=2).grid(row=0,column=2,sticky='W',padx=2)
            sdir.set(os.path.abspath(os.curdir))

            fm=ttk.Frame(master=wd); fm.grid(row=2,column=0,sticky='NW',pady=4)
            ttk.Label(fm,text='filename',width=8).grid(row=0,column=0,sticky='W',padx=1)
            ttk.Entry(fm,textvariable=fname,width=20).grid(row=0,column=1,sticky='W')
            ttk.Button(fm,text='save',width=4,command=savefile).grid(row=0,column=2,sticky='W',padx=3)
            resize(); wd.mainloop()

        def init_save():
            #import multiprocessing as mp
            #p=mp.Process(target=_save_item,args=[self,gcf()]); p.start()
            import threading
            threading.Thread(target=_save_item,args=(self,gcf())).start()

        if self.backend==1:
           acs=self.toolbar.actions(); ats=array([i.iconText() for i in acs]); self.hts=[]
           abp=acs[pindex(ats,'save')[0]] if 'save' in ats else gcf().canvas.toolbar.addAction('save')
           abp.triggered.connect(init_save)
        else:
           acs=self.toolbar.toolitems; ats=[i[0] for i in acs]; self.hts=[]
           if 'save' in ats:
              ac=acs[ats.index('save')][1]
           else:
              import tkinter as tk
              ac=tk.Button(gcf().canvas.toolbar,text='save', command=init_save); ac.pack(side=tk.LEFT)
              gcf().canvas.toolbar.toolitems=(*gcf().canvas.toolbar.toolitems,('save',ac))

class schism_bpfile(zdata):
    def __init__(self,x=None,y=None,z=None,station=None,fmt=0):
        self.nsta=0; self.x=array([]); self.y=array([]); self.z=array([])
        self.station=[]; self.hp=[]; self.ht=[]; self.fmt=fmt
        if x is not None: self.x=x
        if y is not None: self.y=y
        if z is not None: self.z=z
        if station is not None: self.station=station
        bkn=mpl.get_backend().lower(); self.backend=1 if (bkn in ['qt5agg','qtagg']) else 2 if (bkn in ['tkagg',]) else 0
        self.check(); self.edit()

    @property
    def xyz(self):
        return c_[self.x,self.y,self.z]

    @property
    def xy(self):
        return c_[self.x,self.y]

    @property
    def cxy(self):
        return self.x+1j*self.y

    @property
    def dist(self):
        if not hasattr(self,'_dist'): self._dist=r_[0,cumsum(abs(diff(self.cxy)))]
        return self._dist

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

    def inside(self,xy,fmt=0,prj=None):
        '''
        check whether pts c_[x,y] are inside the polygon of bp points.
        fmt=0: return indices of pts inside; fmt=1: return boolean flag
        fmt=2: return indices of pts outside region; fmt=3: return boolean flag outside
        prj=['prj1','prj2']: convert xy from prj1 to prj2
        '''
        if prj is not None: xy=c_[proj_pts(*xy.T,*prj)].T
        fp=inside_polygon(xy,self.x,self.y)==1
        return pindex(fp) if fmt==0 else fp if fmt==1 else nindex(fp) if fmt==2 else ~fp
    def outside(self,xy,fmt=2):
        return self.inside(xy,fmt=fmt)

    def index(self,station,z=None):
        '''
        return station index of indices
        1). for single station, return all the indice of the station
        2). for a list of stations, return indice of stations (return 1st index if multiple indice exist)
        '''
        if isinstance(station,str):
           sid=pindex(self.station==station if (z is None) else (self.station==station)*(self.z==z))
           return None if len(sid)==0 else sid[0] if len(sid)==1 else sid
        else:
           return array([pindex(self.station==i)[0] if (i in self.station) else None for i in station])

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
        self.toolbar=gcf().canvas.toolbar
        gcf().canvas.draw()

    def compute_acor(self,gd):
        #compute areal coordinates, and gd is the schism grid
        self.ie,self.ip,self.acor=gd.compute_acor(c_[self.x,self.y])
        return self.ie,self.ip,self.acor

    def disconnect_edit(self):
        if hasattr(self,'cidpress'): gcf().canvas.mpl_disconnect(self.cidpress)
        if hasattr(self,'cidmove'):  gcf().canvas.mpl_disconnect(self.cidmove)
        if self.backend==1:
           acs=gcf().canvas.toolbar.actions(); ats=[i.iconText() for i in acs]; ap=acs[ats.index('Pan')]
           if ap.isChecked(): ap.trigger()
        gcf().canvas.draw()

    def edit(self):
        def connect_actions():
            self.cidmove=gcf().canvas.mpl_connect('motion_notify_event', onmove)
            self.cidpress=gcf().canvas.mpl_connect('button_press_event', onclick)
            if self.nsta!=0 and len(self.hp)==0: self.plot_station()
            #stop other bp/reg edit mode
            at='reg' if self.fmt==0 else 'bp'
            if self.backend==1:
               ap=[i for i in gcf().canvas.toolbar.actions() if i.iconText()=='Pan'][0]
               ab=[i for i in gcf().canvas.toolbar.actions() if i.iconText()==at][0]
               if not ap.isChecked(): ap.trigger()
            elif self.backend==2:
               ab=[i[1] for i in gcf().canvas.toolbar.toolitems if i[0]==at][0]
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
            distp=squeeze(abs((self.x-x)+1j*(self.y-y))); sid=pindex(distp,distp.min())[0]
            self.x=r_[self.x[:sid],self.x[(sid+1):]]; self.y=r_[self.y[:sid],self.y[(sid+1):]]; self.plot()

        def move_pt(xi,yi):
            distp=squeeze(abs((self.x-xi)+1j*(self.y-yi))); sid=pindex(distp,distp.min())[0]
            self.x[sid]=xi; self.y[sid]=yi; self.plot()

        if mpl._pylab_helpers.Gcf.get_active() is not None:
            at='bp' if self.fmt==0 else 'reg'
            if self.backend==1:
               acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
               abp=acs[pindex(ats,at)[0]] if at in ats else gcf().canvas.toolbar.addAction(at)
            elif self.backend==2:
               acs=gcf().canvas.toolbar.toolitems; ats=[i[0] for i in acs]
               if at in ats:
                  abp=acs[ats.index(at)][1]
               else:
                  import tkinter as tk
                  abp=tk.Button(gcf().canvas.toolbar,text =at, command=connect_actions); abp.pack(side=tk.LEFT)
                  gcf().canvas.toolbar.toolitems=(*gcf().canvas.toolbar.toolitems,(at,abp))
            if self.backend!=0:
               #disconnect and clean previous bpfile
               if hasattr(abp,at):
                  if self is not abp.bp:
                     for i in abp.bp.ht: i.remove()
                     for i in abp.bp.hp: i.remove()
                     abp.bp.ht=[]; abp.bp.hp=[]
                  if self.backend==1: abp.triggered.disconnect()

               #connect to new object
               if self.backend==1: abp.triggered.connect(connect_actions)
               abp.bp=self
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

def read_schism_th(fname,StartT=0,fmt=0):
    '''
    read schism *th file
    fmt=0: convert time (second) to day; fmt=1: not convert
    StartT=0: add starting date for the time
    '''
    C=zdata(); f=loadtxt(fname); C.time,C.data=f[:,0],f[:,1:]
    if fmt==0: C.time=C.time/86400
    C.time=C.time+StartT
    return C

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

class schism_vgrid(zdata):
    def __init__(self):
        pass

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
            self.sigma=array(self.sigma).astype('float'); self.kbp=0
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
    def save(self,fname='vgrid.in',fmt=0,**args):
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

def read_schism_grid(source='.',fmt=0):
    '''
    smartly read schism hgrid and vgrid from source
    1). schism rundir
    2). grid.npz or hgrid.npz
    3). hgrid.gr3, hgrid.ll, *.gr3, *.ic

    fmt=0: hgrid;  fmt=1: vgrid;  fmt=2: hgrid and vgrid
    '''
    if nargout()==2: fmt=2
    def _endswith(fn,svars): #check str ending
        return max([1 if fn.endswith(i) else 0 for i in svars])==1

    if fmt in [0,2]: #read hgrid
       if source.endswith('.npz'): #1). npz format
          svars=read(source,'vars')
          gd=read(source,'hgrid') if ('hgrid' in svars) else read(source) if {'x','y','elnode'}.issubset(svars) else None
          if gd is None: sys.exit('unknown *.npz format for schism_hgrid: '+source)
       elif _endswith(source,['.gr3','.ll','.ic']): #grid file
            gd=read(source)
       elif isinstance(source,str): #rundir
            fns=glob(os.path.abspath(source)+os.sep+'*'); fnames=[]
            [fnames.append(i) for i in fns if i.endswith('grid.npz')]; [fnames.append(i) for i in fns if i.endswith('hgrid.gr3')]
            [fnames.append(i) for i in fns if i.endswith('hgrid.ll')]; [fnames.append(i) for i in fns if i.endswith('hgrid.ll')]
            [fnames.append(i) for i in fns if i.endswith('.gr3')]; [fnames.append(i) for i in fns if i.endswith('.ic')]
            gd=read_schism_grid(fnames[0])
    
    if fmt in [1,2]: #read vgrid
       if source.endswith('.npz'): #1). npz format
          svars=read(source,'vars')
          vd=read(source,'vgrid') if ('vgrid' in svars) else read(source) if {'ivcor','nvrt'}.issubset(svars) else None
          if vd is None: sys.exit('unknown *.npz format for schism_vgrid: '+source)
       elif source.endswith('.in'): #grid file
            vd=read(source)
       elif isinstance(source,str): #rundir
            fns=glob(os.path.abspath(source)+os.sep+'*'); fnames=[]
            [fnames.append(i) for i in fns if i.endswith('grid.npz')]; [fnames.append(i) for i in fns if i.endswith('vgrid.in')]
            vd=read_schism_grid(fnames[0],1)
    return gd if fmt==0 else vd if fmt==1 else [gd,vd] 

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
            kbp=array([pindex(abs(i+1)<1e-10)[-1] for i in sigma])

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
            sindc=nindex(fps); fpc=eta[sindc]<=(-h_c-(hmod[sindc]-h_c)*theta_f/sinh(theta_f)); sindf=sindc[fpc]; sinds=sindc[~fpc]
            if sum(fpc)>0 and ifix==0: sys.exit('Pls choose a larger h_c: {}'.format(h_c))
            if sum(fpc)>0 and ifix==1: zcor[k:,sindf]=(eta[sindf]+hmod[sindf])*sigma[kin]+eta[sindf]
            zcor[k,sinds]=eta[sinds]*(1+sigma[kin])+h_c*sigma[kin]+cs[kin]*(hmod[sinds]-h_c)
        zcor[:(kz-1),:]=nan if fmt==1 else zcor[kz-1,:][None,:] #extend

        #for z layer
        kbp=zeros(np).astype('int'); fp=dp>h_s; sind=pindex(fp); kbp[~fp]=kz-1
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
        sindp=pindex(ds,gd.np); sinde=pindex(ds,gd.ne); sindk=pindex(ds,nvrt)

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
            z1,z2=zcor[:,min([nvrt-1,n+1])],zcor[:,n]
            if n==(nvrt-1):
                fpz=pzi==z1; ratz=zeros(sum(fpz))
            else:
                fpz=(pzi>z1)*(pzi<=z2); ratz=(pzi[fpz]-z1[fpz])/(z2[fpz]-z1[fpz])
            if sum(fpz)==0: continue

            for i,[value,ds,ndim,ip,iz] in enumerate(zip(values,dims,ndims,pind,zind)):
                if ip==-1 or iz==-1: continue
                v1,v2=value[fpz,min([nvrt-1,n+1])],value[fpz,n]; rat=ratz.copy()
                for m in arange(ndim-2):rat=expand_dims(rat,axis=1)
                pvalues[i][fpz,k]=v1*(1-rat)+v2*rat

    #restore dimension order
    for i,[pvalue,sind] in enumerate(zip(pvalues,tind)):
        if sind is None: continue
        sinds=argsort(sind); pvalues[i]=pvalues[i].transpose(sinds)
    if ilst==0: pvalues=pvalues[0]

    return pvalues

def check_schism_ihot(dirpath='.'):
    '''
    check schism run setup with ihot=2, and abort when outputs/flux.out missing
    '''
    bdir=os.path.abspath(dirpath); fname=bdir+'/param.nml'; odir=bdir+'/outputs'
    if not (os.path.exists(fname) and os.path.exists(odir)): return
    ihot=read_schism_param(fname,3).ihot
    if ihot==2 and (not os.path.exists(odir+'/flux.out')): sys.exit('wrong setup: flux.out missig for ihot={}'.format(ihot))

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
          if line[:max([line.find('='),0])].strip()==param:
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
    elif 'schism_file.schism_grid' in str(type(grd)):
       gd=grd
    else:
       sys.exit('unknow format of grd: {}'.format(grd))

    #save grid save *2dm format
    gd.grd2sms(sms)

class schism_transect(schism_grid):
    def __init__(self,bpfile=None,grid=None,hgrid=None,vgrid=None,eta=None,zcor=None,dist=None):
        '''
        create schism grid from a transect with inputs
          bpfile (optional): build points (npt,)
          grid   (optional): run dir, or grid.npz
          hgrid  (optional): hgrid.gr3 
          vgrid  (optional): vgrid.in
          eta    (optional): surface elevation (npt,)
          zcor   (optional): directly use z-coordinate (npt,nvrt)
          dist   (optional): distance from starting points (npt,)

        usage: 1). gdm=schism_transect('james_transect.bp','RUN10a')
               2). gdm=schism_transect('james_transect.bp','grid.npz')
               3). gdm=schism_transect('james_transect.bp',hgrid='hgrid.gr3',vgrid='vgrid.in')
               4). gdm=schism_transect(zcor=zcor) or gdm=schism_transect(zcor=zcor,dist=dist)
        '''

        #read bpfile and grid
        bp,grd,gd,vd=bpfile,grid,hgrid,vgrid
        bp=bp if ('schism_bpfile' in str(type(bp))) else read(bp) if isinstance(bp,str) else None
        gd=gd if ('schism_grid'  in str(type(gd))) else read_schism_grid(gd)   if (gd is not None) else read_schism_grid(grd)   if (grd is not None) else None
        vd=vd if ('schism_vgrid' in str(type(vd))) else read_schism_grid(vd,1) if (vd is not None) else read_schism_grid(grd,1) if (grd is not None) else None
        if (bp is not None) or (zcor is not None): self.create_transect(bp,gd,vd,eta,zcor,dist)

    def create_transect(self,bp=None,gd=None,vd=None,eta=None,zcor=None,dist=None):
        #get sigma coordiante
        self._eta=0 if (eta is None) else eta
        if (bp!=None)*(gd!=None)*(vd!=None):
           pie,pip,pacor=gd.compute_acor(bp.xy); self.z0=(gd.z[pip]*pacor).sum(axis=1); dist=bp.dist
           self.sigma=(vd.sigma[pip]*pacor[...,None]).sum(axis=1); zcor=compute_zcor(self.sigma,self.z0,self.eta)
        elif zcor is not None:
           npt,nvrt=zcor.shape; dist=arange(npt) if (dist is None) else dist
           self.z0=-zcor[:,0]; self._eta=zcor[:,-1]; self.sigma=(zcor-zcor[:,-1][:,None])/(self.eta+self.z0)[:,None]
        else:
           sys.exit('wrong inputs for creating SCHISM transect')
        
        #create profile holder
        npt,nvrt=zcor.shape; ns=npt*nvrt; ip=arange(ns).reshape([npt,nvrt]); xs=tile
        for k in arange(nvrt-1)[::-1]: fp=zcor[:,k]==zcor[:,k+1]; ip[fp,k]=ip[fp,k+1]  #remove repeated points
        p=unique(array([ip[:-1,1:],ip[:-1,:-1],ip[1:,:-1],ip[1:,1:]]).reshape([4,(nvrt-1)*(npt-1)]),axis=1) #p=elnode,e=elnode
        fpn=(p[0]==p[1])*(p[2]==p[3]); p=p[:,~fpn] #remove invalid elem.
        sindp,sindv=unique(p,return_inverse=True); e=arange(sindp.size)[sindv].reshape(p.shape) #get unique points, and renumber elem.
        e[3,e[2]==e[3]]=-2; fp=e[0]==e[1]; e[:,fp]=r_[e[array([0,2,3])][:,fp],-2*ones([1,sum(fp)],'int')]; e=e.T #for triangle
        
        #create grid
        self.pid=tile(False,[npt,nvrt]); self.pid.ravel()[sindp]=True
        self.x=tile(dist,[nvrt,1]).T[self.pid]; self.y=zcor[self.pid]; self.elnode=e
        self.np=len(sindp); self.ne=len(e); self.i34=sum(e!=-2,axis=1); self.dp=zeros(self.np)

    @property
    def z(self):
        return self.dp
    @z.setter
    def z(self,value):
        if isinstance(value,np.ndarray):
           self.dp=value if value.ndim==1 else value[self.pid] if value.ndim==2 else (ones(self.np)*value).ravel()[:self.np]
        else:
           self.dp=(ones(self.np)*value).ravel()[:self.np]
    @property
    def eta(self):
        return self._eta
    @eta.setter
    def eta(self,value):
        self.y=compute_zcor(self.sigma,self.z0,value)[self.pid]; self._eta=value

def zcor_to_schism_grid(zcor,x=None,value=None):
    '''
    convert z-coordinate transect to a schism grid
    Inputs:
        zcor[npt,nvrt]: z-coordinates for each points
        x[npt]: distances of each pionts from starting point
        value[npt,nvrt]: value assicated for each point @(x,y)
    '''
    gd=schism_transect(zcor=zcor,dist=x); gd.z=gd.z if (value is None) else value[gd.pid]
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
                fpe=fpe*(inside_polygon(gd.exy,rxy[:,0],rxy[:,1])==(0 if m==0 else 1))

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
   ie=pindex(inside_polygon(gd.exy,xy[:,0],xy[:,1]),1)
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
           if not hasattr(P,'pie'): P.pie,P.pip,P.pacor=sgrid.compute_acor(c_[lx,ly],out=0)

    #extract time series@xyz or transect@xy
    mtime=[]; mdata=[[] for i in varname]
    stacks=dstacks if (stacks is None) else [*array(stacks).ravel()] #check outputs stacks
    for istack in stacks:
        print('reading stack: {}'.format(istack)); sys.stdout.flush()
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
                    vii=vii.transpose([2,0,1]) if nd==3 else vii.transpose([2,0,1,3]) #from (nt,npt,nvrt) to (nvrt,nt,npt)
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
       levels:  schism level indices (1-nvrt: surface-bottom; (>nvrt): kbp level; "all": all layers (note: vertical dim is reversed))
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
                C=C0 if (svar in dvars_2d) else read('{}/{}_{}.nc'.format(bdir,svar,istack),1); cvar=C.variables[svar]; vi=[]; npt=cvar.shape[1]
                if svar in dvars_2d:  #2D
                   vi=array([array(cvar[i]).astype('float32') for i in arange(nt) if i%nspool==0])
                   if reg is not None: vi=vi[:,sindp if npt==np else sinde if npt==ne else sinds] #subset
                else:   #3D
                   if levels=='all':
                      vi=array([flipud(array(cvar[i]).astype('float32').T) for i in arange(nt) if i%nspool==0])
                   else:
                      for n,k in enumerate(zs):
                          if k>nvrt:
                             if (not hasattr(P,'kbp')) and npt==np: vd=read(fgz,'vgrid') if fexist(fgz) else read(fvd); P.sindp,P.kbp=arange(np),vd.kbp
                             if (not hasattr(P,'kbe')) and npt==ne: vd=read(fgz,'vgrid') if fexist(fgz) else read(fvd); P.sinde,P.kbe=arange(ne),gd.compute_kb(vd.kbp)
                             sindi,kbi=[P.sindp,P.kbp] if npt==np else [P.sinde,P.kbe]
                             vii=array([array(cvar[i][sindi,kbi]).astype('float32') for i in arange(nt) if i%nspool==0])
                          else:
                             vii=array([array(cvar[i,:,nvrt-k]).astype('float32') for i in arange(nt) if i%nspool==0])
                          if reg is not None: vii=vii[:,sindp if npt==np else sinde if npt==ne else sinds] #subset
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
           'zcor':'zcor','hvel':'hvel', 'hvel_side':'hvel_side','zvel_elem':'wvel_elem',
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

    ICM={'PB1':'ICM_PB1','PB2':'ICM_PB2','PB3':'ICM_PB3','CHLA':'ICM_CHLA',  #core: phytoplankton
        'RPOC':'ICM_RPOC','LPOC':'ICM_LPOC','DOC':'ICM_DOC', #core: carbon
        'RPON':'ICM_RPON','LPON':'ICM_LPON','DON':'ICM_DON','NH4':'ICM_NH4','NO3':'ICM_NO3', #core: nitrogen
        'RPOP':'ICM_RPOP','LPOP':'ICM_LPOP','DOP':'ICM_DOP','PO4':'ICM_PO4', #core: phosphorus
        'COD':'ICM_COD','DO':'ICM_DOX', #core: COD and DO
        'SRPOC':'ICM_SRPOC','SRPON':'ICM_SRPON','SRPOP':'ICM_SRPOP','PIP':'ICM_PIP', #SRM module
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

class schism_view(zdata):
    def __init__(self, run='.',scaling=None):
        #note: p is a capsule including all information about a figure
        self.figs=[]; self.fns=[]; self._nf=0; self.sbp=[] #list of figure objects
        self.run_info(run)
        self.window, self.wp=self.init_window(scaling); self.window.title('SCHSIM Visualization : '+self.run+' (Author: Z. WANG)')
        self.cmap='jet'
        self.hold='off'  #animation
        self.play='off'  #animation
        self.curve_method=0 #the method in extracting time series (0: nearest, 1: interpolation)
        self.itp=0 #transect plot
        self.window.mainloop()
        #todo: 1). for cases that out2d*.nc not exist

    def init_plot(self,fmt=0):
        w=self.wp; fn=w.fn.get()
        if fn=='add': #new figure
            if fmt==1: return #return if entry from control window
            p=zdata(); self.get_param(p); #p.xm,p.ym=[self.xm,self.ym] if p.map=='none' else [self.xml,self.yml]
            p.hf=figure(figsize=p.figsize,num=self.nf())
            cid=len(self.fns); self.figs.append(p); self.fns.append('{}: {}'.format(len(self.fns)+1,p.var))
        else: #old figure
            cid=self.fns.index(fn); p=self.figs[cid]
            if not fignum_exists(p.hf.number): p.hf=figure(figsize=p.figsize,num=p.hf.number) #restore closed figure
            if fmt==0: #modify current figure
                self.get_param(p); self.fns[cid]=self.fns[cid][:3]+p.var; p.hf.clear(); p.bm=None
            elif fmt==1: #restore old figure setting
                self.update_panel('old',p)
        w._fn['values']=['add',*self.fns]; w.fn.set(self.fns[cid])
        self.fig=p; p.itp=self.itp; figure(p.hf.number) #bring figure to front
        return p

    def schism_plot(self,fmt=0):
        if self.wp.var.get() not in self.pvars: print(self.wp.var.get()+' not exist'); return
        if not hasattr(self,'hgrid'): print('wait: still reading grid'); return
        if self.play=='on' and fmt==1: self.play='off'; self.hold='off'; return
        if self.hold=='on': return #add hold to avoid freeze when user press too frequently
        p=self.init_plot(0) if fmt==0 else self.init_plot(1)
        w=self.wp; gd=self.hgrid if self.itp==0 else p.td
        if p is None: return
        if fmt==2: p.it=max([p.it-p.ns,0])
        if fmt==3: p.it=min([p.it+p.ns,len(self.irec)-1])
        if fmt==4: p.it=0
        if fmt==5: it=len(self.irec)-1; p.it=it; p.it2=it; self.update_panel('it2',p)

        #plot figure and save the backgroud
        self.hold='on'; mask=None if p._nan=='none' else p.nan
        mls=self.mls if w.time.get()=='time' else self.julian if w.time.get()=='julian' else self.stacks
        if fmt==0:
           p.hp=[]; p.hg=[]; p.hb=[]; p.hv=[]; anim=True if p.med==0 else False
           if p.var!='none':
               self.get_data(p); v=self.data
               if p.med==0: p.hp=[gd.plot(fmt=1,method=1,value=v,clim=p.vm,mask=mask,ticks=11,animated=True,cmap=self.cmap,zorder=1,cb_aspect=50)]
               if p.med==1: p.hp=[gd.plot(fmt=1,method=0,value=v,clim=p.vm,mask=mask,ticks=11,cmap=self.cmap,zorder=1,cb_aspect=50)]
           if p.vvar!='none': u,v=self.get_vdata(p); p.hv=[quiver(p.vx,p.vy,u,v,animated=anim,scale=1.0/p.zoom,scale_units='inches',width=0.001,zorder=3)]
           if p.vvar!='none': quiverkey(p.hv[0], X=0.92, Y=1.01, U=1, label='1.0 m/s',color='r', labelpos='E',zorder=4)
           if p.grid==1: p.hg=gd.plot(animated=anim,zorder=2)
           if p.bnd==1: p.hb=gd.plot_bnd(lw=0.5,alpha=0.5,animated=anim)
           if self.itp==1 and p.bnd==1: p.hb.extend(plot(gd.xm,[0,0],'k:',lw=1,zorder=3))
           if p.map!='none': self.add_map()
           p.ht=title('{}, layer={}, {}'.format(p.var,p.layer,mls[p.it]),animated=anim)

           #add pts for time series
           m=20; n=p.npt; x=array([*p.px,*tile(0.0,m-n)]); y=array([*p.py,*tile(nan,m-n)])
           fpn=nindex((x[:n]>p.xm[0])*(x[:n]<p.xm[1])*(y[:n]>p.ym[0])*(y[:n]<p.ym[1])); x[fpn]=0.0; y[fpn]=nan
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
                    self.get_data(p); v=self.data; self.query(1); s=zdata()
                    if mask is not None:
                       if mask[0]=='=' and mask[1]!='=': mask='='+mask
                       if isnumber(mask): fpnd=v==float(mask); v_nd=v[fpnd]; v[fpnd]=nan
                       if not isnumber(mask): exec('s.fp=v'+mask); fpnd=s.fp; v_nd=v[fpnd]; v[fpnd]=nan
                    if p.med==0:
                        if self.itp==1: p.hp[0]._triangulation.y=gd.y #update transect grid
                        if self.itp==1 and p.grid==1: p.hg[0].set_ydata(gd.lines(0)[:,1])
                        if self.itp==1 and p.bnd==1: p.hb[0].set_ydata(r_[gd.lines(1)][:,1])
                        p.hp[0].set_array(r_[v,v[gd.sindg]] if (v.size==gd.np and gd.wrap==1) else v if v.size==gd.np else r_[v,v[gd.fp4]])
                    else:
                        [i.remove() for i in p.ax.collections] 
                        gd.plot(ax=p.ax,fmt=1,value=v,clim=p.vm,mask=mask,ticks=11,cmap=self.cmap,cb=False,zorder=1)
                    if mask is not None: v[fpnd]=v_nd
                if p.vvar!='none':  #vector
                   u,v=self.get_vdata(p)
                   if p.med==0: p.hv[0].set_UVC(u,v)
                   if p.med==1: p.hv=[quiver(p.vx,p.vy,u,v,scale=1/p.zoom,scale_units='inches',width=0.001,zorder=3)]
                p.ht.set_text('{}, layer={}, {}'.format(p.var,p.layer,mls[p.it]))
                self.update_panel('it',p); self.window.update()
                if p.med==0: p.bm.update()
                if p.anim!=None: savefig('.{}_{:06}'.format(p.anim,p.it)) #save fig for animation
                if p.med==1: pause(0.1)
                if hasattr(p,'pause'): pause(max([p.pause,0.0001]))
                if self.play=='off': break
                if fmt in [2,3,4,5]: self.play='off'
            if fmt==1: w.player['text']='play'; self.window.update()
            if p.anim!=None:
               from PIL import Image
               ims=['.{}_{:06}.png'.format(p.anim,i) for i in  [it0,*its]]; fms=[Image.open(i) for i in ims]; adt=max([p.pause*1e3,50]) if hasattr(p,'pause') else 200
               fms[0].save(p.anim+'.gif',format='GIF', append_images=fms[1:], save_all=True, duration=adt, loop=0)
               [os.remove(i) for i in ims]
        self.hold='off'

    def plotts(self):
        import threading
        #function to extract data
        def get_tsdata(ts,x,y,svar,layer,ik1,ik2):
            w.curve['text']='wait'; ts.x=x; ts.y=y; ts.var=svar; ts.layer=layer; ts.ik1=ik1; ts.ik2=ik2; ts.mys=[]; nt=0
            for ik in arange(ik1,ik2+1):
                def _ts(outputs):
                   fname='{}/out2d_{}.nc'.format(outputs,ik) if svar in self.vars_2d else '{}/{}_{}.nc'.format(outputs,svar,ik)
                   C=self.fid(fname); nt0,npt=C.variables[svar].shape[:2]
                   if ik==ik1 and self.curve_method==0: t.sindp=near_pts(c_[x,y],gd.xy) if npt==gd.np else near_pts(c_[x,y],gd.exy) #compute index
                   if ik==ik1 and self.curve_method==1: pie,pip,pacor=gd.compute_acor(c_[x,y],fmt=1); t.sindp=pip.ravel() if npt==gd.np else pie #compute index for interp
                   if svar in self.vars_2d:
                       data=array(C.variables[svar][:,t.sindp])
                   else:
                       ks=(self.kbp[t.sindp] if npt==gd.np else self.kbe[t.sindp]) if layer=='bottom' else (-tile(1 if layer=='surface' else int(layer),t.sindp.size))
                       data=array([C.variables[svar][:,i,k] for i,k in zip(t.sindp,ks)]).T
                   if npt==gd.np and self.curve_method==1: data=sum(reshape(data,[nt0,*pip.shape])*pacor[None,...],axis=2)
                   return data,nt0,fname
                t00=time.time(); data,nt0,fname=_ts(self.outputs)
                if p.cmp==1: data=data-_ts(p.run0+os.path.sep+'outputs')[0]
                ts.mys.extend(data); nt=nt+nt0; print('extracting {} from {}: {:0.2f}'.format(svar,fname,time.time()-t00))
            ts.mys=array(ts.mys).T; ts.mt=array(self.mts[it1:(it1+nt)]); ts.mls=array(self.mls[it1:(it1+nt)]); p.ts=ts
            print('done in extracting'); w.curve['text']='curve'

        def update_xts(event):
            if event!=0 and type(event)!=mpl.backend_bases.DrawEvent: return
            t1,t2=xlim(); dt1=abs(mt-t1); dt2=abs(mt-t2); i1=pindex(dt1,dt1.min())[0]; i2=pindex(dt2,dt2.min())[0]
            ns=max([int(floor((i2-i1+1)/5)),1]); mti=mt[i1:i2:ns]; mlsi=mls[i1:i2:ns]
            if hasattr(self,'StartT'): mlsi=[i[:10]+'\n'+i[11:] for i in mlsi]
            s.ax.set_xticks(mti); s.ax.set_xticklabels(mlsi)

        #prepare info. about time sereis
        if not hasattr(self,'fig'): return
        p=self.fig; w=self.wp; gd=self.hgrid; gd.compute_ctr()
        svar,layer=p.var,p.layer; x=array(p.px); y=array(p.py)
        if svar=='depth' or len(x)==0: return
        ik1=self.istack[p.it]; ik2=self.istack[p.it2-1]; it1=self.istack.index(ik1); it2=len(self.istack)-self.istack[::-1].index(ik2); t=zdata()
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

    def profile(self,sp=None,itp=None):
        import tkinter as tk
        w=self.wp; tp=w.tp
        if itp==None: #for transect
           if not hasattr(self,'td'): self.update_transect()
           td=self.td; gd=self.hgrid
           if tp['relief'].lower()=='sunken':
              tp.config(relief=tk.RAISED,bg='gray88'); self.itp=0; gd=self.hgrid
           else:
              tp.config(relief=tk.SUNKEN,bg='grey'); self.itp=1; gd=self.td
           xm,ym=gd.xm,gd.ym; w.xmin.set(xm[0]); w.xmax.set(xm[1]); w.ymin.set(ym[0]); w.ymax.set(ym[1])
        else:
           tp.config(relief=tk.RAISED,bg='gray88' if itp==0 else 'grey'); self.itp=itp

    def add_map(self):
        if hasattr(self,'mp'):
           mp=self.mp; sd={'Image':'World_Imagery','Topo':'World_Topo_Map','Street':'World_Street_Map'}; p=self.fig
           mp.llcrnrlon=p.xm[0]; mp.urcrnrlon=p.xm[1]; mp.llcrnrlat=p.ym[0]; mp.urcrnrlat=p.ym[1]
           mp.arcgisimage(xpixels=1500,service=sd[p.map],cachedir=self.runpath)
        else:
           from mpl_toolkits.basemap import Basemap
           xm,ym=self.xml,self.yml; self.mp=Basemap(llcrnrlon=xm[0],urcrnrlon=xm[1],llcrnrlat=ym[0],urcrnrlat=ym[1],epsg=4326); self.add_map()

    def query(self,sp=None):
        if not (hasattr(self,'fig') and hasattr(self.fig,'hpt') and hasattr(self.fig,'hf')): return
        p=self.fig; qr=self.wp.query; hf=p.hf; gd=self.hgrid
        hp=p.hpt[0]; ht=p.hpt[-1]; x,y=hp.get_xdata(),hp.get_ydata()

        if sp is None: #set query state
           import tkinter as tk
           qr.config(relief=tk.RAISED,bg='gray88') if qr['relief'].lower()=='sunken' else qr.config(relief=tk.SUNKEN,bg='grey')
           if hasattr(p,'qxy'): delattr(p,'qxy')
           x[-1]=nan; y[-1]=nan; hp.set_xdata(x); hp.set_ydata(y); ht.set_text(''); hf.canvas.draw()
        elif qr['relief'].lower()=='sunken':
            if sp not in [0,1]: #save query pt info
               dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
               if not (dlk==0 and btn==2): return
               p.qxy=[bx,by,*gd.compute_acor(c_[bx,by])]
            #update query
            if not hasattr(p,'qxy'): return
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
            x=array(p.px); y=array(p.py); distp=squeeze(abs((x-bx)+1j*(y-by))); sid=pindex(distp,distp.min())[0]
            hp=p.hpt[0]; x,y=hp.get_xdata(),hp.get_ydata(); x[sid]=bx; y[sid]=by; hp.set_xdata(x); hp.set_ydata(y)
            p.hpt[sid+1].set_x(bx); p.hpt[sid+1].set_y(by); p.px[sid]=bx; p.py[sid]=by
        if p.med==0: p.bm.update(); pause(0.001)
        if p.med==1: p.hf.canvas.draw()

    def get_data(self,p):  #slab data
        svar,layer,istack,irec,pdict=p.var,p.layer,self.istack[p.it],self.irec[p.it],self.pdict; gd=self.hgrid; idry=self.wp.dry.get(); itp=0
        if p.var=='depth': self.data=gd.dp; return
        if p.var in self.gr3: self.data=read_schism_hgrid(self.runpath+os.sep+p.var).z; return
        if self.itp==1 and hasattr(p,'td'): itp=1; tp,td=p.td.tp,p.td
        def _data0(output,var): #get original value of 'var' from outputs dir.
           C0=self.fid('{}/out2d_{}.nc'.format(output,istack))
           if var in self.vars_2d:
               if itp==0: #slab data
                  data=array(C0.variables[var][irec])
               elif itp==1: #transect data
                  p.td.eta=array(C0.variables['elevation'][irec][tp.ip]*tp.acor).sum(axis=1)
                  data=array(C0.variables[var][irec][tp.ip]*tp.acor).sum(axis=1)
           else:
               C=self.fid('{}/{}_{}.nc'.format(output,var,istack))
               if itp==0: #slab
                 if layer=='bottom':
                     npt=C.variables[var].shape[1]
                     if npt==gd.np: data=array(C.variables[var][irec][arange(gd.np),self.kbp])
                     if npt==gd.ne: data=array(C.variables[var][irec][arange(gd.ne),self.kbe])
                 else:
                     ilayer=1 if layer=='surface' else int(layer); data=array(C.variables[var][irec,:,-ilayer])
               elif itp==1: #transect
                 npt=C.variables[var].shape[1]
                 p.td.eta=array(C0.variables['elevation'][irec][tp.ip]*tp.acor).sum(axis=1)
                 if npt==gd.np: data=array(C.variables[var][irec][tp.ip]*tp.acor[...,None]).sum(axis=1)
                 if npt==gd.ne: data=array(C.variables[var][irec][tp.ie])
           if itp==0: data[abs(data)>1e20]=nan
           if itp==1 and data.ndim==2:
              for k in arange(data.shape[1]-1)[::-1]: fp=abs(data[:,k])>1e20; data[fp,k]=data[fp,k+1]
              data=data[td.pid]
           if idry==1: #add wetting and drying mask
              npt=len(data);  dvar='dryFlagNode' if npt==gd.np else 'dryFlagElement' if npt==gd.ne else 'dryFlagSide'
              fmask=pindex(array(C0.variables[dvar][irec]),1); data[fmask]=nan
           return data
        def _data1(output,var):  #get defined variable
           s=zdata(); [s.attr(i,_data0(self.outputs,i) if pdict[i]==0 else _data1(self.outputs,i)) for i in pdict[var].vars]
           exec('s.data={}'.format(pdict[var].exp)); return s.data

        self.data=_data0(self.outputs,svar) if pdict[svar]==0 else _data1(self.outputs,svar)
        if p.cmp==1: odir=p.run0+os.path.sep+'outputs'; self.data=self.data-_data0(odir,svar) if pdict[svar]==0 else _data1(odir,svar)

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
           x,y=(gd.xy if npt==gd.np else gd.exy if npt==gd.ne else gd.sxy).T
           p.sindv=pindex((x>=p.xm[0])*(x<=p.xm[1])*(y>=p.ym[0])*(y<=p.ym[1])); p.vx=x[p.sindv]; p.vy=y[p.sindv]

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
                w.StartT.set(mls[0]); w.EndT.set(mls[-1]); w._StartT['values']=mls; w._EndT['values']=mls; #w._StartT['width']=6; w._EndT['width']=6
                if event=='time': w._StartT['width']=20; w._EndT['width']=20
        elif event=='vm': #reset data limit
            if not hasattr(self,'hgrid'): return
            p=self.get_param(); s=self.pdict[p.var]
            if p.var!='none': self.get_data(p); w.vmin.set(self.data.min()); w.vmax.set(self.data.max())
            if s!=0: print('\ndefine var: {}\nevaluation: {}={}'.format(s.equation,s.name,s.exp))
            w._nan.set('none'); self.update_panel('mask')
        elif event=='mask': #set mask 
            w.sfm.grid(row=0,column=2,sticky='W') if w._nan.get()=='mask' else w.sfm.grid_forget()
        elif event=='old': #reset all panel variables from a previous figure
            w.var.set(p.var); w.fn.set(p.fn); w.layer.set(p.layer); w.time.set(p.time)
            w.StartT.set(p.StartT); w.EndT.set(p.EndT); w.vmin.set(p.vm[0]); w.vmax.set(p.vm[1])
            w._nan.set(p._nan); w.nan.set(p.nan); self.update_panel('mask'); w.cmp.set(p.cmp)
            w.xmin.set(p.xm[0]); w.xmax.set(p.xm[1]); w.ymin.set(p.ym[0]); w.ymax.set(p.ym[1]); w.map.set(p.map)
            w.ns.set(p.ns); w.grid.set(p.grid); w.bnd.set(p.bnd); w.vvar.set(p.vvar); w.zoom.set(p.zoom); w.med.set(p.med)
            self.profile(itp=p.itp)
            if p.run0 is not None: self.run0=p.run0
        elif type(event)==mpl.backend_bases.DrawEvent:
            if not fignum_exists(self.fig.hf.number): return
            p=self.fig; ax=p.ax; xm=ax.get_xlim(); ym=ax.get_ylim(); p.xm=[*xm]; p.ym=[*ym]
            w=self.wp; w.xmin.set(xm[0]); w.xmax.set(xm[1]); w.ymin.set(ym[0]); w.ymax.set(ym[1])
            p.figsize=[p.hf.get_figwidth(),p.hf.get_figheight()]
        self.window.update()

    def get_param(self,p=None):
        if p is None: p=zdata()
        w=self.wp; p.var=w.var.get(); p.fn=w.fn.get(); p.layer=w.layer.get(); p.time=w.time.get()
        p.StartT=w.StartT.get(); p.EndT=w.EndT.get(); p.vm=[w.vmin.get(),w.vmax.get()]; p._nan=w._nan.get(); p.nan=w.nan.get()
        p.xm=[w.xmin.get(),w.xmax.get()]; p.ym=[w.ymin.get(),w.ymax.get()]; p.med=w.med.get(); p.map=w.map.get()
        p.ns=w.ns.get(); p.grid=w.grid.get(); p.bnd=w.bnd.get(); p.vvar=w.vvar.get(); p.zoom=w.zoom.get()
        p.anim=None; p.sindv=None; p.figsize=[7.2,5.5]
        p.cmp=w.cmp.get(); p.run0=self.run0 if hasattr(self,'run0') else None
        if not hasattr(p,'npt'): p.npt=0; p.px=[]; p.py=[]
        if self.itp==1 and not hasattr(p,'td'): p.td=self.td

        #get time index
        if p.time=='time':
            p.it =self.mls.index(p.StartT) if p.StartT in self.mls else p.it if hasattr(p,'it') else 0
            p.it2=self.mls.index(p.EndT)+1 if p.EndT in self.mls else p.it2 if hasattr(p,'it2') else len(self.mls)
        elif p.time=='julian':
            p.it=self.julian.index(float(p.StartT)); p.it2=self.julian.index(float(p.EndT))+1
        else: #stacks
            p.it=self.istack.index(int(p.StartT)); p.it2=len(self.istack)-self.istack[::-1].index(int(p.EndT))
        return p

    def update_xy(self):
        w=self.wp; gd=self.hgrid; mp=w.map.get()
        if gd.ics==1:
           gd.x,gd.y,xm,ym=[gd.x0,gd.y0,self.xm,self.ym] if mp=='none' else [gd.lon,gd.lat,self.xml,self.yml]
           if hasattr(self,'fig'):
              p=self.fig
              if int(p.map=='none') +int(mp=='none')==1: w.xmin.set(xm[0]); w.xmax.set(xm[1]); w.ymin.set(ym[0]); w.ymax.set(ym[1]); p.xm=xm; p.ym=ym
              p.map=mp
           else:
              w.xmin.set(xm[0]); w.xmax.set(xm[1]); w.ymin.set(ym[0]); w.ymax.set(ym[1])

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
        for i in glob(self.runpath+os.sep+'World_*.png'): os.remove(i)
        close('all'); self.window.destroy()

    def run_info(self,run):
        import threading,time

        #check output
        print('reading grid and output info.')
        fns=glob(run+'/out2d_*.nc'); fns2=glob(run+'/outputs/out2d_*.nc'); self.gr3=[os.path.basename(i) for i in glob(run+'/*.gr3')]; iout=1
        self.outputs,fnames=[run,fns] if len(fns)!=0 else [run+'/outputs',fns2]; self.runpath=os.path.abspath(run); self.run=os.path.basename(self.runpath)
        if len(fnames)!=0:
           [iks,self.vars,self.vars_2d]=get_schism_output_info(self.outputs)[2:]; iks=sort(iks)
           try:
             self.fid('{}/out2d_{}.nc'.format(self.outputs,iks[-1]))
           except:
             iks=iks[:-1]
           ik0=iks[0]; C=self.fid('{}/out2d_{}.nc'.format(self.outputs,ik0)); cvar=C.variables; cdim=C.dimensions
           self.vvars=[i[:-1] for i in self.vars if (i[-1]=='X') and (i[:-1]+'Y' in self.vars)]
        else:
           self.vars=[]; self.svars_2d=[]; self.vvars=[]; iout=0; cvar=None
           print('schism outputs dir. not found')
        self.pvars=['none','depth',*self.vars,*self.gr3]; self.pdict=dict(zip(self.pvars,zeros(len(self.pvars))))

        #read grid and param
        grd=run+'/grid.npz'; gr3=run+'/hgrid.gr3'; grl=run+'/hgrid.ll'; vrd=run+'/vgrid.in'; par=run+'/param.nml'; fexist=os.path.exists; fpath=os.path.abspath
        if os.path.exists(par): p=read_schism_param(par,3); self.param=p; self.StartT=datenum(p.start_year,p.start_month,p.start_day,p.start_hour)
        def _read_grid():
           def _read_hgrid():
               gd0=read_schism_hgrid(gr3); gd2=read_schism_hgrid(grl) if (fexist(grl) and fpath(gr3)!=fpath(grl)) else gd0; gd0.lon,gd0.lat=gd2.x,gd2.y
               return gd0
           gd=loadz(grd).hgrid if fexist(grd) else _read_hgrid() if fexist(gr3) else None
           vd=loadz(grd).vgrid if fexist(grd) else read_schism_vgrid(vrd) if fexist(vrd) else None
           if gd==None: #create hgrid
              if cvar==None: return
              gd=get_schism_output_info(self.outputs,4)
           if not hasattr(gd,'x0'): gd.x0,gd.y0=[gd.x,gd.y]
           if not hasattr(gd,'lon'): gd.lon,gd.lat=[gd.x,gd.y]
           gd.ics=1 if abs(gd.x-gd.lon).max()>1e-5 else 2
           self.hgrid=gd; self.xm=gd.xm; self.ym=gd.ym; self.xml=[gd.lon.min(),gd.lon.max()]; self.yml=[gd.lat.min(),gd.lat.max()]
           self.vm=gd.zm; self.fp3=pindex(gd.i34,3)
           self.kbp, self.nvrt=[vd.kbp, vd.nvrt] if vd!=None else [array(cvar['bottom_index_node']), cdim['nSCHISM_vgrid_layers'].size]; self.kbe=gd.compute_kb(self.kbp)
           while not hasattr(self,'wp'): time.sleep(0.01)
           w=self.wp; w._layer['values']=['surface','bottom',*arange(2,self.nvrt+1)]; print('schismview ready')
           w.vmin.set(self.vm[0]); w.vmax.set(self.vm[1]); w.xmin.set(self.xm[0]); w.xmax.set(self.xm[1]); w.ymin.set(self.ym[0]); w.ymax.set(self.ym[1])
           if gd.lon.min()<-180 or gd.lon.max()>180 or gd.lat.min()<-90 or gd.lat.max()>90: w._map['values']=['none',]
           try:
             from mpl_toolkits.basemap import Basemap
           except:
             w._map['values']=['none',]
        self.nvrt=2; self.xm=[0,1]; self.ym=[0,1]; self.vm=[0,1]
        threading.Thread(target=_read_grid).start()
        if iout==0: self.stacks=[0]; self.julian=[0]; self.istack=[0]; self.irec=[0]; self.mls=['0']; return
        self.run_info_time()

    def run_info_time(self,fmt=0):
        '''
        read time information of outputs
        '''
        iks=sort(get_schism_output_info(self.outputs)[2]); fnames=glob(self.outputs+'/out2d_*.nc'); ik0=iks[0]
        C=self.fid('{}/out2d_{}.nc'.format(self.outputs,ik0)); cvar=C.variables

        #read available time
        self.stacks=[]; self.julian=[]; self.istack=[]; self.irec=[]
        ti=array(cvar['time'])/86400; nt=len(ti); t0=ti[0]  #assume all stacks have the same number of records
        if (not hasattr(self,'StartT')) and hasattr(C.variables['time'],'base_date'):
           self.StartT=datenum(*[int(float(i)) for i in C.variables['time'].base_date.split()])
        if len(iks)>1 and nt==1:
            C1=self.fid('{}/out2d_{}.nc'.format(self.outputs,iks[1])); ti1=array(C1.variables['time'])/86400; dt=ti1-ti[0]
        else:
            dt=0 if nt in [0,1] else ti[1]-ti[0]
        for ik in iks:
            fn='{}{}out2d_{}.nc'.format(self.outputs,os.path.sep,ik); nt1=nt
            if fn not in fnames: continue
            if ik==iks[-1]:
               try:
                   C=self.fid(fn); nt1=C.dimensions['time'].size
                   if nt1==0: continue
               except:
                   continue #last stack not readable
            self.julian.extend(t0+(ik-ik0)*nt*dt+arange(nt1)*dt); self.irec.extend(arange(nt1))
            self.stacks.append(ik); self.istack.extend(tile(ik,nt1))

        self.mts=self.julian[:]; self.mls=['{}'.format(i) for i in self.mts]
        if hasattr(self,'StartT'): #get time_string
           def _set_time():
               self.mls=[i.strftime('%Y-%m-%d,%H:%M:%S') for i in num2date(self.mts)]
               while not hasattr(self,'wp'): time.sleep(0.01)
               self.wp._StartT['values']=self.mls; self.wp._EndT['values']=self.mls
           self.mts=[*(array(self.mts)+self.StartT)]
           self.mls[0]=num2date(self.mts[0]).strftime('%Y-%m-%d,%H:%M:%S')
           self.mls[-1]=num2date(self.mts[-1]).strftime('%Y-%m-%d,%H:%M:%S')
           threading.Thread(target=_set_time).start()

        if fmt==1: #update panel and parameter
           p=self.fig if hasattr(self,'fig') else None; mls=self.mls; self.wp.EndT.set(mls[-1])
           if p!=None: p.EndT=mls[-1]; p.it2=self.mls.index(mls[-1])+1

    def cmd_window(self,fmt,scaling):
        import tkinter as tk
        from tkinter import ttk, filedialog,font
        def cmd_from_file():
            fname=filedialog.askopenfilename(initialdir=self.runpath, title = "choose file")
            if len(fname)==0 or fname.strip()=='': return
            fid=open(fname,'r'); txt.insert('1.0',''.join([*fid.readlines(),'\n'])); fid.close()
        if hasattr(self,'fig'): p=self.fig; hf,ax=p.hf,p.ax
        T0='control vars: [self,p,hf,ax]'; T1='define new variable: rho=999.8+0.8*salt-0.2*temp'; print(T0 if fmt==0 else T1)
        cw=tk.Toplevel(self.window); cw.geometry("400x300"); cw.title('command input' if fmt==0 else 'define variables')
        cw.rowconfigure(0,minsize=150, weight=1); cw.columnconfigure(0,minsize=2, weight=1)
        txt=tk.Text(master=cw,width=150,height=14); txt.grid(row=0,column=0,pady=2,padx=2,sticky='nsew')
        sbar=tk.Scrollbar(cw, orient="vertical", command=txt.yview); sbar.grid(row=0,column=1,sticky="ns"); txt.config(yscrollcommand=sbar.set)
        sfm=ttk.Frame(master=cw); sfm.grid(row=1,column=0,padx=10)
        if fmt==0:
           rbn=ttk.Button(sfm, text= "run",command=lambda: self.cmd_exec(txt.get('1.0',tk.END)))
        else:
           rbn=ttk.Button(sfm, text= "add",command=lambda: self.add_var(txt.get('1.0',tk.END)))
        fbn=ttk.Button(sfm, text= "file",command=cmd_from_file); dbn=ttk.Button(sfm, text= "clean",command=lambda: txt.delete("1.0",tk.END))
        fbn.grid(row=0,column=0,padx=10); dbn.grid(row=0,column=1,padx=10); rbn.grid(row=0,column=2,padx=10)
        cw.update(); xm=max([txt.winfo_width(),rbn.winfo_width()]); ym=txt.winfo_height()+rbn.winfo_height()+12
        if fmt==0: txt.insert('1.0',self.cmd0) if hasattr(self,'cmd0') else txt.insert('1.0','#'+T0)
        if fmt==1: txt.insert('1.0',self.cmd1) if hasattr(self,'cmd1') else txt.insert('1.0','#'+T1)
        if scaling is not None: #scaling
           srat=float(scaling); fs=int(10*srat); xm=int(xm*srat); ym=int(ym*srat)
           DF=font.Font(family="Arial",size=fs); ttk.Style().configure("Custom.TButton", font=("Arial", fs))
           txt.config(font=DF); rbn.config(style="Custom.TButton"); fbn.config(style="Custom.TButton"); dbn.config(style="Custom.TButton")
        cw.geometry('{}x{}'.format(xm,ym)); cw.update()

    def cmd_exec(self,cmd):
        self.cmd0=cmd
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

    def add_var(self,cmd):
        self.cmd1=cmd; lines=[i.strip() for i in cmd.split('\n') if (not i.startswith('#')) and len(i.strip())!=0]
        modules=get_schism_output_info(self.outputs)[0]
        for line in lines:  #parse each formulation
            try:
              pind=[i for i,k in enumerate(line) if k in "=+-*/()% '"]; s=[] #s contain all items
              s.append(line[:pind[0]]) # first var should be the new variable to be add
              for i in arange(len(pind)-1):  #find all variables and numbers
                  s.append(line[pind[i]])
                  if pind[i]==pind[i+1]-1: continue
                  s.append(line[(pind[i]+1):pind[i+1]])
              s.append(line[pind[-1]])
              if pind[-1]!=len(s)-1: s.append(line[(pind[-1]+1):])

              #add variables information
              C=zdata(); C.name=s[0]; C.exp=s[2:]; C.equation=line; C.vars=[]; iflag=0
              for i,svar in enumerate(C.exp):
                  if (svar in "=+-*/()% '") or (svar[0] in '0123456789') or svar.lower() in ['sqrt','none']: continue
                  vs=get_schism_var_info(svar,modules)
                  if len(vs)!=1 or vs[0][1] not in self.pvars: print('check variable in {}: '.format(line), vs); iflag=1
                  C.exp[i]='s.'+vs[0][1]; C.vars.append(vs[0][1])
              C.exp=''.join(C.exp) #evaluation expression
              if iflag==1: continue
              if (C.exp.lower() in ["''", 'none']) and C.name in self.pvars: #remove variables
                 vid=self.pvars.index(C.name); self.pvars=[*self.pvars[:vid],C.name,*self.pvars[(vid+1):]]; self.wp.svar.config(value=self.pvars)
                 self.pdict.pop(C.name); print('remove variable: '+C.name); continue
              if C.name not in self.pvars: #add variables
                 self.pvars=[*self.pvars[:2],C.name,*self.pvars[2:]]; self.wp.svar.config(value=self.pvars); print('add variable: '+C.name)
              else:
                 print('update variable: '+C.name)
              self.pdict[C.name]=C
            except:
              print('invalid formulation: '+line)

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

    def update_transect(self):
        import tkinter as tk
        from tkinter import filedialog
        fname = filedialog.askopenfilename(title='select SCHISM transect file (e.g. *bp)')
        tp=read_schism_bpfile(fname); tp.ie,tp.ip,tp.acor=tp.compute_acor(self.hgrid)
        td=schism_transect(tp,hgrid=self.hgrid,vgrid=self.runpath); td.tp=tp; self.td=td
        figure(0); self.hgrid.plot_bnd(); plot(tp.x,tp.y,'r.-'); show(block=False)
        if hasattr(self,'fig'): figure(self.fig.hf.number)

    def reset_limit(self):
        p,w=self.fig,self.wp; gd=self.hgrid if self.itp==0 else p.td
        w.xmin.set(gd.xm[0]); w.xmax.set(gd.xm[1]); w.ymin.set(gd.ym[0]); w.ymax.set(gd.ym[1])
        fpn=~isnan(self.data); w.vmin.set(self.data[fpn].min()); w.vmax.set(self.data[fpn].max())
        self.schism_plot(0)

    def show_node(self):
        p=self.fig; gd=self.hgrid
        if hasattr(p,'hns'):
           for i in arange(len(p.hns)): p.hns.pop().remove()
           delattr(p,'hns'); p.hf.canvas.draw()
        else:
           xm=xlim(); ym=ylim(); p.hns=[]
           if not hasattr(gd,'xctr'): gd.compute_ctr()
           sind=pindex((gd.x>=xm[0])*(gd.x<=xm[1])*(gd.y>=ym[0])*(gd.y<=ym[1]))
           for i in sind: ht=text(gd.x[i],gd.y[i],'{}'.format(i+1),fontsize=6,zorder=3); p.hns.append(ht)
           sind=pindex((gd.xctr>=xm[0])*(gd.xctr<=xm[1])*(gd.yctr>=ym[0])*(gd.yctr<=ym[1]))
           for i in sind: ht=text(gd.xctr[i],gd.yctr[i],'{}'.format(i+1),fontsize=6,zorder=3); p.hns.append(ht)
           p.hf.canvas.draw()
        return

    def anim_exec(self):
        p=self.fig; anim=self.fig._anim.get()
        p.anim=anim[:-4] if anim.endswith('.gif') else anim if anim.strip()!='' else None
        self.play='off'; self.schism_plot(1); p.anim=None

    def schism_instance(self,event):
        from tkinter import filedialog
        import subprocess

        #get directory
        if event=='compare':
           icmp=self.wp.cmp.get(); sdir=self.run0 if hasattr(self,'run0') else self.runpath
           if icmp==0: return
        else:
           sdir=self.runpath
        crun=filedialog.askdirectory(initialdir=sdir, title = "choose run for "+ event)

        if not (isinstance(crun,str) and os.path.exists(crun)): return
        if event in ['schism_check','schism_view']: #open another window for schismview or schismcheck
           p=subprocess.Popen('python3 -c "from pylib import *; {}(\'{}\')"'.format(event,crun),shell=True);  self.sbp.append(p)
        else: #add base run
           self.run0=os.path.abspath(crun)

    def init_window(self,scaling=None):
        #open an window
        import tkinter as tk
        from tkinter import ttk,font

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
        svar=ttk.Combobox(fm,textvariable=w.var,values=self.pvars,width=15,); svar.grid(row=0,column=1)
        svar.bind("<<ComboboxSelected>>",lambda x: self.update_panel('vm')); w.svar=svar

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
        w.map=tk.StringVar(wd); w.map.set('none')
        tk.Checkbutton(master=sfm2,text='grid',variable=w.grid,onvalue=1,offvalue=0).grid(row=0,column=0)
        tk.Checkbutton(master=sfm2,text='bnd',variable=w.bnd,onvalue=1,offvalue=0).grid(row=0,column=1,sticky='W')
        tk.Checkbutton(master=sfm2,text='ctr',variable=w.med,onvalue=1,offvalue=0).grid(row=0,column=2,sticky='W')
        ttk.Label(master=sfm2,text=', map').grid(row=0,column=3,sticky='W')
        w._map=ttk.Combobox(master=sfm2,textvariable=w.map,values=['none','Image','Topo','Street'],width=6); w._map.grid(row=0,column=4,sticky='W')
        w._map.bind("<<ComboboxSelected>>",lambda x: self.update_xy())

        #time
        w.time=tk.StringVar(wd); w.StartT=tk.StringVar(wd); w.EndT=tk.StringVar(wd); w.mls=self.mls; w.StartT.set(self.mls[0]); w.EndT.set(self.mls[-1])
        ttk.OptionMenu(fm,w.time,'time','time','stack','julian',command=self.update_panel).grid(row=2,column=0,sticky='W',pady=4)
        w._StartT=ttk.Combobox(master=fm,textvariable=w.StartT,values=self.mls,width=20); w._StartT.grid(row=2,column=1,padx=0,sticky='W')
        w._EndT=ttk.Combobox(master=fm,textvariable=w.EndT,values=self.mls,width=20); w._EndT.grid(row=2,column=2,sticky='W',padx=1)

        #limit
        sfm3=ttk.Frame(master=fm); sfm3.grid(row=3,column=2,sticky='W')
        w.vmin=tk.DoubleVar(wd); w.vmax=tk.DoubleVar(wd); w.vmin.set(self.vm[0]); w.vmax.set(self.vm[1])
        ttk.Label(master=fm,text='  limit').grid(row=3,column=0,sticky='W')
        ttk.Entry(fm,textvariable=w.vmin,width=10).grid(row=3,column=1,sticky='W',padx=2)
        ttk.Entry(sfm3,textvariable=w.vmax,width=10).grid(row=0,column=0,sticky='W',padx=0)
        #w._nan=tk.StringVar(wd); w.nan=tk.DoubleVar(wd); w._nan.set('none'); w.nan.set(0)
        w._nan=tk.StringVar(wd); w.nan=tk.StringVar(wd); w._nan.set('none'); w.nan.set('0')
        ttk.OptionMenu(sfm3,w._nan,'','none','mask',command=lambda x: self.update_panel('mask')).grid(row=0,column=1,padx=5,sticky='W')
        sfm4=ttk.Frame(master=sfm3); ttk.Entry(sfm4,textvariable=w.nan,width=10).grid(row=0,column=0,sticky='W'); w.sfm=sfm4

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
        w.tp=tk.Button(master=fm0,text='profile',bg='gray88',command=self.profile,width=5); w.tp.grid(row=0,column=2,padx=1,pady=2)
        w.query=tk.Button(master=fm0,text='query',bg='gray88',command=self.query,width=5); w.query.grid(row=0,column=3,padx=1,pady=2)

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
        menu=tk.Menu(mbar,tearoff=0); w.cmp=tk.IntVar(wd); w.dry=tk.IntVar(wd)
        menu.add_command(label="command", command=lambda: self.cmd_window(0,scaling))
        menu.add_command(label="reset",   command=self.reset_limit)
        menu.add_command(label="add_variable", command=lambda: self.cmd_window(1,scaling))
        menu.add_command(label="update outputs",   command=lambda: self.run_info_time(1))
        menu.add_command(label="save animation", command=self.anim_window)
        menu.add_checkbutton(label="wetting/drying",onvalue=1,offvalue=0,variable=w.dry)
        menu.add_command(label="show node/element", command=self.show_node)
        menu.add_command(label="schism_check", command=lambda: self.schism_instance('schism_check'))
        menu.add_command(label="schism_view", command=lambda: self.schism_instance('schism_view'))
        menu.add_checkbutton(label="schism_compare",onvalue=1,offvalue=0,variable=w.cmp,command=lambda: self.schism_instance('compare'))
        menu.add_command(label="schism_transect", command=self.update_transect)
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
        wd.geometry('600x180'); wd.update(); ym=210

        #scaling all elements
        if scaling is not None:
           def get_elem(fm):
               xs=[]; [xs.extend(get_elem(i)) if isinstance(i,ttk.Frame) else xs.append(i) for i in fm.winfo_children()]
               return xs
           xs=get_elem(wd); srat=float(scaling); fs=int(srat*10); ym=int(srat*ym) #scaling factor
           iw,ih=wd.geometry().split('+')[0].split('x'); wd.geometry('{}x{}'.format(int(srat*int(iw)),int(srat*int(ih))))
           #DF=font.Font(family="TkDefaultFont",size=fs); ttk.Style().configure("Custom.TButton", font=("TkDefaultFont", fs))
           #wd.option_add("*TCombobox*Listbox.font", ("TkDefaultFont", fs)) #font size in dropdown list
           DF=font.Font(family="Arial",size=fs); ttk.Style().configure("Custom.TButton", font=("Arial", fs))
           wd.option_add("*TCombobox*Listbox.font", ("Arial", fs)) #font size in dropdown list
           [i.config(font=DF) for i in xs if isinstance(i,ttk.Label) or isinstance(i,tk.Button) or isinstance(i,ttk.Entry) or isinstance(i,tk.Checkbutton) or isinstance(i,ttk.Combobox)]
           [i.config(style="Custom.TButton") for i in xs if isinstance(i,ttk.Button) or isinstance(i,ttk.OptionMenu) or isinstance(i,ttk.Menubutton)]
           [i['menu'].config(font=DF) for i in xs if isinstance(i,ttk.OptionMenu)]
           menu.config(font=DF); wd.update() #menu is treated seperately

        #adjust last frame
        xm=max([i.winfo_width() for i in fms])
        for i in arange(0,100):
            L1['text']=' '*i; L2['text']=' '*i; wd.update(); xs=fms[-1].winfo_width()
            if xs>xm: wd.geometry('{}x{}'.format(xs,ym)); wd.update(); break
        return wd,w

class schism_check(zdata):
   def __init__(self,run='.',scaling=None):
       self.params,self.figs,self.fids,self.fmts,self.scaling={},{},{},{},scaling

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
           if fname.startswith('hotstart.nc'): fmts[fname]=2
           if fname=='source.nc': fmts[fname]=3
           if fname=='source_input':
              sname='.source.nc'; fmts[fname]=3
              if not os.path.exists(self.run+sname): print('convert schism source_sink format: '+sname); convert_schism_source(self.run,sname)
           if fname=='*.th': fmts[fname]=4
           if fname not in fmts:
              if fname.endswith('.nc'):
                 fmts[fname]=2
              else:
                 sys.exit('unknown format: {}'.format(fname))
           self.fmt=fmts[fname]; self.read_input_info()
       self.fmt=fmts[fname]; p=params[fname]

       for widget in fm.winfo_children(): widget.destroy() #clean plot option frame
       self.info.config(text=p.info); ap.dns=[]
       #design plot option frame
       if self.fmt==0: #update gr3 file parameters and panel
          if p.init==0:
             p.grid=tk.IntVar(wd); p.bnd=tk.IntVar(wd); p.ctr=tk.IntVar(wd); p.grid.set(0); p.bnd.set(0); p.ctr.set(1)
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.sflux=tk.StringVar(wd); p.vmin.set(0); p.vmax.set(0); p.sflux.set('None')
             p._nan=tk.StringVar(wd); p.nan=tk.DoubleVar(wd); p._nan.set('none'); p.nan.set(0)

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
          ttk.OptionMenu(sfm,p._nan,'','none','mask',command=lambda x: self.update_panel('mask')).grid(row=0,column=3,padx=5,sticky='W')
          sfm2=ttk.Frame(master=sfm); ttk.Entry(sfm2,textvariable=p.nan,width=5).grid(row=0,column=0,sticky='W'); p.sfm=sfm2
          if option=='mask': p.sfm.grid(row=0,column=4,sticky='W') if p._nan.get()=='mask' else p.sfm.grid_forget()
       elif self.fmt==1: # *D.th.nc or _nu.nc
          if p.init==0:
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.vmin.set(0); p.vmax.set(0)
             p.transpose=tk.IntVar(wd); p.transpose.set(0); p._nan=tk.StringVar(wd); p.nan=tk.DoubleVar(wd); p._nan.set('none'); p.nan.set(0)
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
          if fname.endswith('_nu.nc'):
             ttk.OptionMenu(sfm,p._nan,'','none','mask',command=lambda x: self.update_panel('mask')).grid(row=0,column=4,padx=5,sticky='W')
             sfm2=ttk.Frame(master=sfm); ttk.Entry(sfm2,textvariable=p.nan,width=5).grid(row=0,column=0,sticky='W'); p.sfm=sfm2
             if option=='mask': p.sfm.grid(row=0,column=5,sticky='W') if p._nan.get()=='mask' else p.sfm.grid_forget()

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
          if p.init==1 and (option in [1,'mask']): #save parameter
             p0=zdata(); p0.dvars=[i.get() for i in p.dvars]; p0.vmin=p.vmin.get(); p0.vmax=p.vmax.get(); p0.scale=p.scale.get()
             p0.transpose=p.transpose.get(); p0.grid=p.grid.get(); p0.bnd=p.bnd.get(); p0.dims=[*p.dims]
          if p.init==0:
             p.vmin=tk.DoubleVar(wd); p.vmax=tk.DoubleVar(wd); p.scale=tk.DoubleVar(wd); p.vmin.set(0); p.vmax.set(0); p.scale.set(1)
             p.transpose=tk.IntVar(wd); p.grid=tk.IntVar(wd); p.bnd=tk.IntVar(wd); p.transpose.set(0); p.grid.set(0); p.bnd.set(0)
             p._nan=tk.StringVar(wd); p.nan=tk.DoubleVar(wd); p._nan.set('none'); p.nan.set(0)
             if self.fmt==3: p.sctr=tk.IntVar(wd); p.srat=tk.DoubleVar(wd); p.ctr=tk.IntVar(wd); p.id=tk.IntVar(wd); p.sctr.set(0); p.srat.set(1); p.ctr.set(0); p.id.set(0)
          if option==0: self.var.set(p.var); self.vars['values']=p.vars

          #update panel
          fid=self.fids[fname]; p.var=self.var.get(); p.cvar=fid[p.var]; p.isource,p.fvar=[1,fid['vsource']] if p.var=='msource' else [0,0]
          p.dims=p.cvar.shape; p.dnames=[*p.cvar.dimensions]
          if p.var in['su2','sv2']: p.dims=[*p.dims,1]; p.dnames=[*p.dnames,'uv']
          p.info='  dim={}'.format(p.dims); self.info.config(text=p.info)
          p.dvars=[tk.StringVar(wd) for i in p.dims];
          p.xs=[*p.dims] if self.fmt==2 else [[*arange(i)] for i in p.dims]

          sfm1=ttk.Frame(master=fm); sfm1.grid(row=0,column=0,sticky='W',pady=5); inode=0
          if not hasattr(self,'hgrid'): self.read_hgrid()
          for n,[dn,ds,dvar] in enumerate(zip(p.dnames,p.dims,p.dvars)):
              if ds in [self.hgrid.np,self.hgrid.ne]: inode=1
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
              if self.fmt==2 and dn.startswith('dim_'):
                 if ds==self.hgrid.np: dn='node'; dvar.set('all')
                 if ds==self.hgrid.ne: dn='elem'; dvar.set('all')
                 if ds!=self.hgrid.np and ds!=self.hgrid.ne: dn='D{}'.format(n)
              sfm11=ttk.Frame(master=sfm1); sfm11.grid(row=0,column=n,sticky='W',pady=5)
              ttk.Label(master=sfm11,text='  '+dn).grid(row=0,column=0,sticky='W'); p.dnames[n]=dn
              mm=ttk.Combobox(sfm11,textvariable=dvar,values=vs,width=dw,); mm.grid(row=0,column=1,sticky='W')
              mm.bind("<<ComboboxSelected>>",lambda x: self.set_anim()); ap.dns.append(mm)

          if self.fmt==2:
             if p.var in ['su2','sv2']:
                ttk.Label(master=sfm1,text='  zoom').grid(row=0,column=len(p.dims)+1,sticky='W')
                ttk.Entry(sfm1,textvariable=p.scale,width=5).grid(row=0,column=len(p.dims)+2,sticky='W',padx=2)
             else:
                if inode==1:
                  ttk.OptionMenu(sfm1,p._nan,'','none','mask',command=lambda x: self.update_panel('mask')).grid(row=0,column=len(p.dims)+2,padx=5,sticky='W')
                  sfm2=ttk.Frame(master=sfm1); ttk.Entry(sfm2,textvariable=p.nan,width=5).grid(row=0,column=0,sticky='W'); p.sfm=sfm2
                  if option=='mask': p.sfm.grid(row=0,column=len(p.dims)+3,sticky='W') if p._nan.get()=='mask' else p.sfm.grid_forget()
                  wd.geometry('420x185')

          #add limit panel
          sfm=ttk.Frame(master=fm); sfm.grid(row=1,column=0,sticky='W')
          ttk.Label(master=sfm,text='  limit').grid(row=0,column=0,sticky='W')
          ttk.Entry(sfm,textvariable=p.vmin,width=6).grid(row=0,column=1,sticky='W',padx=2)
          ttk.Entry(sfm,textvariable=p.vmax,width=6).grid(row=0,column=2,sticky='W')
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
          if p.init==1 and (option in [1,'mask']) and array_equal(array(p0.dims),array(p.dims)):
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
       if self.scaling is not None: self.scaling_window()
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
       if fmt==2: ts=[ds[max([iv-ns,iv0])],]  #back one
       if fmt==6: ts=[ds[min([iv+ns,nrec-1])],] #forward one
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
          gd=self.hgrid; p=self.params[fname]; data=fids[fname]; pfmt=0; fxy=0; mask=None if p._nan.get()=='none' else p.nan.get()
          if p.sflux.get()!='None':
             if not hasattr(gd,'lon'): gd0=read_schism_hgrid(self.run+'hgrid.ll'); gd.lon=gd0.x; gd.lat=gd0.y
             sid=read(self.run+'sflux'+os.sep+p.sflux.get()+'.nc',1); sx=array(sid['lon'][:]); sy=array(sid['lat'][:]); sid.close()
             for i,k in zip(sx,sy): plot(i,k,'-',color='orange',lw=0.5,alpha=1,zorder=1); fxy=1
             for i,k in zip(sx.T,sy.T): plot(i,k,'-',color='orange',lw=0.5,alpha=1,zorder=1)
             if abs(gd.lon-gd.x).mean()>1: pfmt=-1; slimit(gd.lon,gd.lat,data)
          if p.ctr.get()==1:  gd.plot(fmt=1,value=data,clim=[p.vmin.get(),p.vmax.get()],mask=mask,ticks=11,cmap='jet',method=1,xy=fxy,zorder=0)
          if p.grid.get()==1: gd.plot(xy=fxy,zorder=2)
          if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=0.5,alpha=0.5,xy=fxy,zorder=2)
          p.hp=gca(); slimit(gd.x,gd.y,data)
       elif ((self.fmt==2 and (p.var in ['su2','sv2'])) or (self.fmt==1 and fname=='uv3D.th.nc')) and p.data.ndim==2 and p.data.shape[1]==2: #vector (su2,sv2) in hotstart, or uv3d.th.nc
          if not hasattr(self,'hgrid'): self.read_hgrid()
          gd=self.hgrid
          if p.grid.get()==1: gd.plot()
          if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=0.5,alpha=0.5)
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
            gd=self.hgrid; srat=p.srat.get(); isource=p.isource if hasattr(p,'isource') else 0
            sind=p.isc if p.var in ['source_elem','msource','vsource'] else p.isk; xi,yi=gd.xctr[sind],gd.yctr[sind]; eid=arange(len(sind))
            data=ones(p.data.shape) if p.var in ['source_elem','sink_elem'] else p.data.copy(); data2=p.data2 if isource==1 else data.copy()
            fpn=data!=-9999; xi,yi,data,eid,data2=xi[fpn],yi[fpn],data[fpn],eid[fpn],data2[fpn] #remove -9999 values
            if data.max()<=0: data=-data #plot negative values (vsink)
            fpn=data>0; xi,yi,data,eid,data2=xi[fpn],yi[fpn],data[fpn],eid[fpn],data2[fpn] #only keep data>0
            if data.size==0: print('no valid points found!'); close(hf); return

            #plot and label
            if p.ctr.get()==0:
               hg=scatter(xi,yi,s=data*srat,c='r')
            else:
               if isource==1 and srat<0 and p.ctr.get()==1:  
                  hg=scatter(xi,yi,s=-data2*srat,c=data,cmap='jet'); data=data2
               else:
                  hg=scatter(xi,yi,s=srat*10,c=data,cmap='jet')
            p.hp=gca(); slimit(gd.x,gd.y,data); pfmt=2
            if p.grid.get()==1: gd.plot()
            if p.bnd.get()==1:  gd.plot_bnd(c='k',lw=0.3)
            if p.id.get()==1: #plot source id
               if hasattr(p,'fmt')  and hasattr(p,'xm') and p.fmt==pfmt: fp=(xi>=p.xm[0])*(xi<=p.xm[1])*(yi>=p.ym[0])*(yi<=p.ym[1]); xi,yi,eid=xi[fp],yi[fp],eid[fp]
               for xii,yii,eidi in zip(xi,yi,eid): text(xii,yii,'{}'.format(eidi),fontsize=7)

            #legend
            if p.ctr.get()==0 or (isource==1 and srat<0 and p.ctr.get()==1):
               v1,v2=data.min(),data.max();  m1,m2=int(log10(v1)),int(log10(v2))
               m1=max([0,m1]) if m2>=0 else m2; ms=[i for i in arange(m1,m2+1) if (10.0**i>=v1) and (10.0**i<=v2)]
               if len(ms)==0:
                  hl=legend(*hg.legend_elements("sizes", num=[(v1+v2)/2])) #legend
               else:
                  hl=legend(*hg.legend_elements("sizes", num=[abs(srat)*10.0**i for i in ms])) #legend
                  for i,m in enumerate(ms): hl.texts[i].set_text('$10^{'+str(m)+'}$')  #set legend value
            #colorbar
            if p.ctr.get()==1 or (isource==1 and srat<0 and p.ctr.get()==1):
               cm.ScalarMappable.set_clim(hg,vmin=vm[0],vmax=vm[1])
               hc=colorbar(fraction=0.05,aspect=50,spacing='proportional',pad=0.02); hc.set_ticks(linspace(*vm,11)); hc.ax.set_ylim(vm)
       elif self.fmt in [1,2,3,4]: # bnd, nudge, hotstart, source.nc
          vm=[p.vmin.get(),p.vmax.get()]
          if p.data.ndim==1:
              i0=p.ax[0]; xi,xn=p.xs[i0], p.dnames[i0]
              if not hasattr(self,'hgrid'):self.read_hgrid()
              if fname.endswith('_nu.nc') and xn=='node' and p.ctr.get()==1:
                  gd=self.hgrid; mask=None if p._nan.get()=='none' else p.nan.get()
                  if not hasattr(p,'sindn'): p.sindn=array(read(self.run+os.path.sep+fname,1).variables['map_to_global_node'])-1
                  vi=zeros(gd.np); vi[p.sindn]=p.data
                  gd.plot(fmt=1,value=vi,clim=[p.vmin.get(),p.vmax.get()],mask=mask,ticks=11,cmap='jet',method=1); p.hp=gca()
                  if p.grid.get()==1: gd.plot()
                  if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=0.5,alpha=0.5)
                  pfmt=6; slimit(gd.x,gd.y,p.data)
              elif self.fmt==2 and (xn in ['node', 'elem','node/elem','dim_{}'.format(self.hgrid.np),'dim_{}'.format(self.hgrid.ne)]): #schism grid plot
                  gd=self.hgrid; mask=None if p._nan.get()=='none' else p.nan.get()
                  gd.plot(fmt=1,value=p.data,clim=[p.vmin.get(),p.vmax.get()],mask=mask,ticks=11,cmap='jet',method=1); p.hp=gca()
                  if p.grid.get()==1: gd.plot()
                  if p.bnd.get()==1:  gd.plot_bnd(c='rg',lw=0.5,alpha=0.5)
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
                 hg=contourf(yi,xi,p.data,levels=linspace(*vm,50),extend='both',cmap='jet')
                 xlabel(yn); ylabel(xn); pfmt=5; slimit(yi,xi,p.data)
              else:
                 hg=contourf(xi,yi,p.data.T,levels=linspace(*vm,50),extend='both',cmap='jet')
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
             if hasattr(p,'ym'): y1,y2=p.ym; y2=y1+max([y2-y1,1e-10]); setp(p.hp,ylim=[y1,y2])
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
          dind=','.join([':' if (i in ['all','mean','min','max','sum']) else str(i) for i in dns]); isource=p.isource if hasattr(p,'isource') else 0
          if isource==1: dns2=[dns[0],dns[2]]; dind2=','.join([':' if (i in ['all','mean','min','max','sum']) else str(i) for i in dns2])
          if p.var in ['su2','sv2']: dind=dind[:-2]
          if p.var in ['su2','sv2'] and dns[-1]!='0': #deal with vector in hotstart
             fid=fids[fname]; exec('p.data=array(concatenate((fid["su2"][{}][...,None],fid["sv2"][{}][...,None]),axis=-1))'.format(dind,dind))
          else:
             exec('p.data=array(cvar[{}])'.format(dind))
             if isource==1: exec('p.data2=array(p.fvar[{}])'.format(dind2))
          if (p.vmin.get()==0 and p.vmax.get()==0) or (hasattr(p,'itr') and p.itr!=p.dvars[-1].get()) or isnan(p.vmin.get()+p.vmax.get()):
              p.vmin.set(p.data.min()); p.vmax.set(p.data.max())
          if self.fmt==1: p.itr=p.dvars[-1].get()
          p.info='  dim={}, [{}, {}]'.format(p.dims,'{:15f}'.format(p.data.min()).strip(),'{:15f}'.format(p.data.max()).strip())

          #sanity check
          errmsg=' found in "{}, {}[{}]" !!'.format(fname,p.var,dind)
          if sum(isnan(p.data))!=0: print('nan value'+errmsg); p.data=None; return #check nan values
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
          if isource==1: #for source.nc
             isht=0
             for n, dn in enumerate(dns2):
                 if dn in ['mean','min','max','sum']: exec('p.data2=p.data2.{}(axis={})'.format(dn,n-isht))
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
       [fnames.append(i) for i in snames if i=='surface_restore.nc']  #surface_restore.nc
       [fnames.append('source_input') for i in snames if i=='source_sink.in']           #source_input
       [snames.remove(i) for i in snames if i in ['source_sink.in','vsource.th','vsink.th','msource.th']]  #remove source_input
       snames=[i for i in snames if not i.startswith('vsource.')]                                        #remove vsource.th
       [fnames.append(i) for i in snames if i.endswith('.nc') and (i not in ['.source.nc',*fnames]) and (not i.startswith('.th_'))]  #other nc files
       fnames.extend(unique(['*.th' for i in snames if i.endswith('.th')]))          #*.th
       [fnames.append(i) for i in snames if i.endswith('.ll')]                       #hgrid.gr3
       [fnames.append(i) for i in snames if i.endswith('.gr3') and (i not in fnames)]#gr3
       [fnames.append(i) for i in snames if i.endswith('.prop')]                     #prop
       mc=[i for i in snames if i.endswith('.ic') and ('hvar' in i)]                 #ic files
       [fnames.append(i) for i in snames if i.endswith('.ic') and (i not in mc)]; fnames.extend(mc)   #ic
       self.fnames=fnames; self.thfiles=[i for i in snames if i.endswith('.th')]     #*.th
       self.sflux=[i[:-3] for i in os.listdir(self.run+'sflux') if i.endswith('nc')] if fexist(self.run+'sflux') else None
       # self.StartT=0 #update this later

   def cmd_window(self):
        import tkinter as tk
        from tkinter import ttk,filedialog,font
        def cmd_from_file():
            fname=filedialog.askopenfilename(initialdir=self.run, title = "choose file")
            if len(fname)==0 or fname.strip()=='': return
            fid=open(fname,'r'); txt.insert('1.0',''.join([*fid.readlines(),'\n'])); fid.close()
        cw=tk.Toplevel(self.window); cw.geometry("400x300"); cw.title('command input')
        cw.rowconfigure(0,minsize=150, weight=1); cw.columnconfigure(0,minsize=2, weight=1)
        txt=tk.Text(master=cw,width=150,height=14); txt.grid(row=0,column=0,pady=2,padx=2,sticky='nsew')
        sfm=ttk.Frame(master=cw); sfm.grid(row=1,column=0,padx=10)
        rbn=ttk.Button(sfm, text= "run",command=lambda: self.cmd_exec(txt.get('1.0',tk.END)))
        fbn=ttk.Button(sfm, text= "file",command=cmd_from_file); dbn=ttk.Button(sfm, text= "clean",command=lambda: txt.delete("1.0",tk.END))
        fbn.grid(row=0,column=0,padx=10); dbn.grid(row=0,column=1,padx=10); rbn.grid(row=0,column=2,padx=10)
        cw.update(); xm=max([txt.winfo_width(),rbn.winfo_width()]); ym=txt.winfo_height()+rbn.winfo_height()+12
        txt.insert('1.0',self.cmd) if hasattr(self,'cmd') else txt.insert('1.0','#'+'control vars: [fname,p,fid,hf,ax,gd]')
        print('control vars: [fname,p,fid,hf,ax,gd]')
        if self.scaling is not None:
           srat=float(self.scaling); fs=int(10*srat); xm=int(xm*srat); ym=int(ym*srat)
           DF=font.Font(family="Arial",size=fs); ttk.Style().configure("Custom.TButton", font=("Arial", fs))
           txt.config(font=DF); rbn.config(style="Custom.TButton"); fbn.config(style="Custom.TButton"); dbn.config(style="Custom.TButton")
        cw.geometry('{}x{}'.format(xm,ym)); cw.update()

   def cmd_exec(self,cmd):
        self.cmd=cmd
        fname=self.fname; p,fid=None,None
        if fname in self.params: p=self.params[fname]
        if fname in self.fids: fid=self.fids[fname]
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

   def scaling_window(self):
       import tkinter as tk
       from tkinter import ttk,font

       wd,scaling=self.window,self.scaling
       def get_elem(fm):
           xs=[]; [xs.extend(get_elem(i)) if isinstance(i,ttk.Frame) else xs.append(i) for i in fm.winfo_children()]
           return xs
       xs=get_elem(wd); srat=float(scaling); fs=int(srat*10)  #scaling factor
       wd.geometry('{}x{}'.format(int(srat*400),int(srat*175)))
       DF=font.Font(family="Arial",size=fs); ttk.Style().configure("Custom.TButton", font=("Arial", fs))
       wd.option_add("*TCombobox*Listbox.font", ("Arial", fs)) #font size in dropdown list
       [i.config(font=DF) for i in xs if isinstance(i,ttk.Label) or isinstance(i,tk.Button) or isinstance(i,ttk.Entry) or isinstance(i,tk.Checkbutton) or isinstance(i,ttk.Combobox)]
       [i.config(style="Custom.TButton") for i in xs if isinstance(i,ttk.Button) or isinstance(i,ttk.OptionMenu) or isinstance(i,ttk.Menubutton)]
       [i['menu'].config(font=DF) for i in xs if isinstance(i,ttk.OptionMenu)]
       wd.update() #menu is treated seperately

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
       if self.scaling is not None: self.scaling_window()

if __name__=="__main__":
    pass
