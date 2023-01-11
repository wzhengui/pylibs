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
    def VINFO(self):
        return get_VINFO(self)

    def plot(self,ax=None,method=0,fmt=0,value=None,mask=None,ec=None,fc=None,
             lw=0.1,levels=None,ticks=None,xlim=None,ylim=None,clim=None,extend='both',cb=True,**args):
        '''
        plot grid with default color value (grid depth)
        method=0: using tricontourf; method=1: using PolyCollection (old method)
        fmt=0: plot grid only; fmt=1: plot filled contours; fmt=2: plot contour lines
        value: color value size(np,or ne)
        mask: size(ne); only plot elements (mask=True))
        ec: color of grid line;  fc: element color; lw: grid line width
        levels=100: number of colors for depths; levels=array([v1,v2,...]): depths for plot
        ticks=[v1,v2,...]: colorbar ticks; ticks=10: number of ticks
        clim=[vmin,vmax]: value range for plot/colorbar
        cb=False: not add colorbar
        '''

        if ec is None: ec='None'
        if fc is None: fc='None'
        if ax is None: ax=gca()

        if method==0:
           fp3=self.i34==3; fp4=self.i34==4
           # if mask is not None: fp3=fp3*mask; fp4=fp4*mask
           if (fmt==0)|(ec!='None'): #compute lines of grid
              #tri
              tri=self.elnode[fp3,:3]; tri=c_[tri,tri[:,0]]
              x3=self.x[tri]; y3=self.y[tri]
              x3=c_[x3,ones([sum(fp3),1])*nan]; x3=reshape(x3,x3.size)
              y3=c_[y3,ones([sum(fp3),1])*nan]; y3=reshape(y3,y3.size)
              #quad
              quad=self.elnode[fp4,:]; quad=c_[quad,quad[:,0]]
              x4=self.x[quad]; y4=self.y[quad]
              x4=c_[x4,ones([sum(fp4),1])*nan]; x4=reshape(x4,x4.size)
              y4=c_[y4,ones([sum(fp4),1])*nan]; y4=reshape(y4,y4.size)

           if fmt==0:
              if ec=='None': ec='k'
              hg=plot(r_[x3,x4],r_[y3,y4],lw=lw,color=ec);
           elif fmt in [1,2]:
              tri=r_[self.elnode[(fp3|fp4),:3],c_[self.elnode[fp4,0],self.elnode[fp4,2:]]]
              #determine value
              if value is None:
                 value=self.dp
              else:
                 if len(value)==self.ne:
                    value=self.interp_elem_to_node(value=value)
                 elif len(value)!=self.np:
                    sys.exit('value has wrong size: {}'.format(value.shape))

              #detemine clim
              if clim is None:
                 fpn=~isnan(value); vmin,vmax=min(value[fpn]),max(value[fpn])
                 if vmin==vmax: vmax=vmax+1e-6
              else:
                 vmin,vmax=clim

              #detemine levels
              if levels is None: levels=51
              if not hasattr(levels,'__len__'): levels=linspace(vmin,vmax,int(levels))

              #set mask
              if sum(isnan(value))!=0: tri=tri[~isnan(value[tri].sum(axis=1))]

              if (vmax-vmin)/(abs(vmax)+abs(vmin))<1e-10:
                 if fmt==1: hg=tricontourf(self.x,self.y,tri,value,vmin=vmin,vmax=vmax,extend=extend,**args)
                 if fmt==2: hg=tricontour(self.x,self.y,tri,value,vmin=vmin,vmax=vmax,extend=extend,**args)
              else:
                 if fmt==1: hg=tricontourf(self.x,self.y,tri,value,levels=levels,vmin=vmin,vmax=vmax,extend=extend,**args)
                 if fmt==2: hg=tricontour(self.x,self.y,tri,value,levels=levels,vmin=vmin,vmax=vmax,extend=extend,**args)

              #add colobar
              cm.ScalarMappable.set_clim(hg,vmin=vmin,vmax=vmax)
              if cb and fmt==1:
                 #----new method
                 hc=colorbar(hg); self.hc=hc
                 if ticks is not None:
                    if not hasattr(ticks,'__len__'):
                       hc.set_ticks(linspace(vmin,vmax,int(ticks)))
                    else:
                       hc.set_ticks(ticks)
                 #----old method
                 #hc=colorbar(hg); self.hc=hc;
                 #if ticks is not None: hc.set_ticks(ticks)
                 #hc.set_clim([vmin,vmax]);

              #plot grid
              if ec!='None': hg=plot(r_[x3,x4],r_[y3,y4],lw=lw,color=ec);

           self.hg=hg; #show(block=False
           if xlim is not None: setp(ax,xlim=xlim)
           if ylim is not None: setp(ax,ylim=ylim)
           if mpl.get_backend().lower() in ['qt5agg','qtagg']:
              acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
              if 'bp' not in ats: self.bp=schism_bpfile()
              if 'query' not in ats: self.query_pt()
              if 'bnd' not in ats: self.create_bnd()
           return hg
        elif method==1:
           #creat polygon
           xy4=c_[self.x,self.y][self.elnode];
           xy4=array([s[0:-1,:] if (i34==3 and len(s)==4) else s for s,i34 in zip(xy4,self.i34)])

           #elem value
           if value is None:
              if not hasattr(self,'dpe'): self.compute_ctr()
              value=self.dpe
           else:
              if len(value)==self.np:
                 value=self.interp_node_to_elem(value=value)
              elif len(value)!=self.ne:
                 sys.exit('value has wrong size: {}'.format(value.shape))

           # apply mask
           if mask is not None: xy4=xy4[mask]; value=value[mask]
              #ind=nonzero(mask)[0]
              #xy4=xy4[ind];
              #value=value[ind]

           #get clim
           if clim is None: clim=[min(value),max(value)]

           #plot
           if fmt==0:
              hg=mpl.collections.PolyCollection(xy4,lw=lw,edgecolor=ec,facecolor=fc,antialiased=False,**args)
           else:
              hg=mpl.collections.PolyCollection(xy4,lw=lw,edgecolor=ec,array=value,clim=clim,antialiased=False,**args)

              #add colorbar
              if cb:
                 hc=colorbar(hg); self.hc=hc;
                 if ticks is not None: hc.set_ticks(ticks)
                 hc.set_clim(clim);

           #add to figure
           ax.add_collection(hg)
           ax.autoscale_view()
           self.hg=hg; #show(block=False)
           if xlim is not None: setp(ax,xlim=xlim)
           if ylim is not None: setp(ax,ylim=ylim)
           return hg
    def plot_grid(self,**args):
        '''
        alias for plot_grid()
        '''
        return self.plot_grid(**args)

    def plot_bnd(self,c='k',lw=0.5,ax=None,**args):
        '''
          plot schims grid boundary

          gd.plot_bnd(): plot bnd
          gd.plot_bnd(c='rb'): open bnd in red,land bnd in blue

        '''
        if ax!=None: sca(ax)
        if len(c)==1: c=c*2
        if not hasattr(self,'nob'): self.compute_bnd()

        #get indices for bnds
        sindo=[]
        for i in arange(self.nob):
            sindo=r_[sindo,-1,self.iobn[i]]
        sindo=array(sindo).astype('int'); fpn=sindo==-1
        bx1=self.x[sindo]; by1=self.y[sindo]
        bx1[fpn]=nan; by1[fpn]=nan

        sindl=[]
        for i in arange(self.nlb):
            if self.island[i]==0:
               sindl=r_[sindl,-1,self.ilbn[i]]
            else:
               sindl=r_[sindl,-1,self.ilbn[i],self.ilbn[i][0]]
        sindl=array(sindl).astype('int'); fpn=sindl==-1
        bx2=self.x[sindl]; by2=self.y[sindl]
        bx2[fpn]=nan; by2[fpn]=nan

        hb1=plot(bx1,by1,c[0],lw=lw,**args)
        hb2=plot(bx2,by2,c[-1],lw=lw,**args)
        #show(block=False)
        self.hb=[hb1,hb2]
        if mpl.get_backend().lower() in ['qt5agg','qtagg']:
           acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
           if 'bp' not in ats: self.bp=schism_bpfile()
           if 'query' not in ats: self.query_pt()
           if 'bnd' not in ats: self.create_bnd()

    def read_hgrid(self,fname,*args):
        #attribute tracking the file originally read, mainly used for savez and save_pkl
        self.source_file = fname

        fid=open(fname,'r'); lines=fid.readlines(); fid.close()

        #read ne and np; lx,ly and dp
        self.ne,self.np=array(lines[1].split()[0:2]).astype('int')
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
        interpolate node values to element values
            default is self.dp => self.dpe
        '''
        #interpolate
        dp=self.dp if (value is None) else value
        fp3=self.i34==3; fp4=~fp3; dpe=zeros(self.ne)
        dpe[fp3]=dp[self.elnode[fp3,:3]].mean(axis=1)
        dpe[fp4]=dp[self.elnode[fp4]].mean(axis=1)
        return dpe

    def interp_elem_to_node(self,value=None,fmt=0,p=1):
        '''
        interpolate element values to nodes
        if value not given, dpe is used
        fmt=0: simple avarage; fmt=1: inverse distance (power=p)
        fmt=2: maximum of surrounding nodal values
        fmt=3: minimum of surrounding nodal values
        '''
        #element values
        if not hasattr(self,'nne'): self.compute_nne()
        if (value is None) and (not hasattr(self,'dpe')): self.compute_ctr()
        v0=self.dpe if (value is None) else value

        #interpolation
        vs=v0[self.ine]
        if fmt==0:
           w=self.ine!=-1; tw=w.sum(axis=1)
           if sum(isnan(value))!=0:
              vs[~w]=0; v=vs.sum(axis=1)/tw
           else:
              v=(w*vs).sum(axis=1)/tw
        if fmt==2: vs[self.ine==-1]=v0.min()-1; v=vs.max(axis=1)
        if fmt==3: vs[self.ine==-1]=v0.max()+1; v=vs.min(axis=1)
        if fmt==1:
              dist=abs((self.xctr[self.ine]+1j*self.yctr[self.ine])-(self.x+1j*self.y)[:,None])
              w=1/(dist**p); w[self.ine==-1]=0; tw=w.sum(axis=1); v=(w*vs).sum(axis=1)/tw
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
        if (not hasattr(self,'ine'))  or fmt==1: self.compute_nne()
        if (not hasattr(self,'ic3'))  or fmt==1: self.compute_ic3()

    def compute_ctr(self):
        '''
        compute element center information: (xctr,yctr,dpe)
        '''
        if not hasattr(self,'xctr'):
           fp3=self.i34==3; fp4=~fp3; self.xctr,self.yctr,self.dpe=zeros([3,self.ne])
           self.xctr[fp3]=self.x[self.elnode[fp3,:3]].mean(axis=1)
           self.xctr[fp4]=self.x[self.elnode[fp4,:]].mean(axis=1)
           self.yctr[fp3]=self.y[self.elnode[fp3,:3]].mean(axis=1)
           self.yctr[fp4]=self.y[self.elnode[fp4,:]].mean(axis=1)
           self.dpe[fp3]=self.dp[self.elnode[fp3,:3]].mean(axis=1)
           self.dpe[fp4]=self.dp[self.elnode[fp4,:]].mean(axis=1)
        return self.dpe

    def compute_area(self):
        fp=self.elnode[:,-1]<0;
        x1=self.x[self.elnode[:,0]]; y1=self.y[self.elnode[:,0]];
        x2=self.x[self.elnode[:,1]]; y2=self.y[self.elnode[:,1]];
        x3=self.x[self.elnode[:,2]]; y3=self.y[self.elnode[:,2]];
        x4=self.x[self.elnode[:,3]]; y4=self.y[self.elnode[:,3]]; x4[fp]=x1[fp]; y4[fp]=y1[fp]
        self.area=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)+(x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
        return self.area

    def compute_gradient(self,fmt=0,value=None,outfmt=0):
        '''
        Compute gradient on each element first, then convert to node value
          fmt: interploation method (see details in interp_elem_to_node())
              (0: simple avarage; fmt=1: inverse distance; fmt=2: maximum of surrounding nodal values)
          value=None: gd.dp is used;    value: array of [np,] or [ne]
          outfmt=0: return (dpdx,dpdy,dpdxy); outfmt=1: return (dpdx,dpdy,dpdxy,dpedx,dpedy,dpedxy)
        '''
        if not hasattr(self,'area'): self.compute_area()
        if not hasattr(self,'dpe'): self.compute_ctr()

        #get node value
        v0=self.dp if value is None else value
        if len(v0)==self.ne: v0=self.interp_elem_to_node(value=v0)

        #get pts
        fp=self.elnode[:,-1]<0; fpn=~fp;
        x1=self.x[self.elnode[:,0]]; y1=self.y[self.elnode[:,0]]; v1=v0[self.elnode[:,0]]
        x2=self.x[self.elnode[:,1]]; y2=self.y[self.elnode[:,1]]; v2=v0[self.elnode[:,1]]
        x3=self.x[self.elnode[:,2]]; y3=self.y[self.elnode[:,2]]; v3=v0[self.elnode[:,2]]
        x4=self.x[self.elnode[:,3]]; y4=self.y[self.elnode[:,3]]; v4=v0[self.elnode[:,3]]
        x4[fp]=x1[fp]; y4[fp]=y1[fp]; v4[fp]=v1[fp]
        a1=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2
        a2=((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2

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
        fmt=2: compute (ns,isidenode,isdel), and (xcj,ycj,dps,distj)
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
           return self.ns,self.isidenode,self.isdel

    def compute_bnd(self,bxy=None):
        '''
        compute boundary information. If bxy is provided, define open/land boundries

        bxy: endpoint coordinates of open boundaries. Examples:
            1). bxy=[x1,x2,y1,y2]  #only one open boundary
            2). bxy=[[x1,x2,y1,y2],[x1,x2,y1,y2],...]  #multiple open boundaries
            3). bxy="bpfile", with paired build points sequentially defining each boundaries
        '''
        print('computing grid boundaries')
        if not hasattr(self,'isdel') or not hasattr(self,'isidenode'): self.compute_side(fmt=1)

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
        sinds=dict(zip(isn[:,0],arange(nbs))) #dict for sides
        ifb=ones(nbs).astype('int'); nb=0; nbn=[]; ibn=[]
        while(sum(ifb)!=0):
            #start points
            id0=isn[nonzero(ifb==1)[0][0],0]; id=isn[sinds[id0],1]; ibni=[id0,id]
            ifb[sinds[id0]]=0; ifb[sinds[id]]=0;
            while True:
                id=isn[sinds[id],1]; ifb[sinds[id]]=0
                if(id==id0): break
                ibni.append(id)
            nb=nb+1; nbn.append(len(ibni)); ibn.append(array(ibni))

        #sort bnd
        nbn=array(nbn); ibn=array(ibn,dtype='O'); fps=flipud(argsort(nbn)); nbn,ibn=nbn[fps],ibn[fps]
        if ibn.shape[0]==1: ibn=ibn.astype('int')

        #find the outline
        island=ones(nb).astype('int')
        for i in arange(nb):
            px=self.x[ibn[i].astype('int')]; i0=nonzero(px==px.min())[0][0]
            sid=ibn[i][array([(i0-1)%nbn[i],i0,(i0+1)%nbn[i]])].astype('int')
            if signa(self.x[sid],self.y[sid])>0: island[i]=0; bid=i; break

        #put outline bnd ahead
        if bid!=0:
           island=ones(nb).astype('int'); island[0]=0
           nbn=array([nbn[bid],*nbn[:bid],*nbn[(bid+1):]])
           ibn=array([ibn[bid],*ibn[:bid],*ibn[(bid+1):]],dtype='O')

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

        #define open/land/island boundaries
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

    def compute_nne(self,**args):
        '''
        alias for compute_node_ball()
        '''
        return self.compute_node_ball(**args)

    def compute_node_ball(self):
        '''
        compute nodal ball information: nne,mnei,indel,ine
        where:
             nne:   number of elements in nodal ball
             mnei:  maximum number of elements in nodal ball
             indel: indices for each nodal ball
             ine:   indices for each nodal ball, but in maxtrix " shape=[np,max(nne)]"
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

    def save(self, fname=None,**args):
        '''
        save hgrid as (*.npz, *.pkl, *.gr3, *.ll or *.ic)
          examples:
                 1). gd.save('grid.npz') or gd.save('grid')
                 2). gd.save('grid.pkl')
                 3). gd.save('hgrid.gr3') or gd.save('hgrid.ll') or gd.save('temp.ic',value=5)
        '''
        if fname is None: fname = '{}.npz'.format(os.path.splitext(self.source_file)[0])
        if fname.endswith('.gr3') or fname.endswith('.ll') or fname.endswith('.ic'):
           self.write_hgrid(fname,**args)
        else:
           s=zdata(); s.hgrid=self; savez(fname,s,**args)

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
        fstr=('{:d} '+fmt+' \n')*self.ne
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

    def split_quads(self,angle_min=60,angle_max=120,fname=None):
        '''
        1). split the quads that have angle (<angle_min or >angle_max), add append the connectivity in the end
        2). fname: output a new grid "fname" if fname!=None
        '''
        if not hasattr(self,'index_bad_quad'): self.check_quads(angle_min,angle_max)

        #compute (angle_max-angle_min) in splitted triangle
        qind=self.index_bad_quad;
        x=self.x[self.elnode[qind,:]]; y=self.y[self.elnode[qind,:]];

        #compute difference between internal angles
        for i in arange(4):
            id1=mod(i-1+4,4); id2=i; id3=mod(i+1,4)
            x1=x[:,id1]; x2=x[:,id2]; x3=x[:,id3];
            y1=y[:,id1]; y2=y[:,id2]; y3=y[:,id3];

            a1=angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2))
            a2=angle((x2-x3)+1j*(y2-y3))-angle((x1-x3)+1j*(y1-y3))
            a3=angle((x3-x1)+1j*(y3-y1))-angle((x2-x1)+1j*(y2-y1))
            a1=mod(a1*180/pi+360,360);a2=mod(a2*180/pi+360,360);a3=mod(a3*180/pi+360,360);

            #compute amax-amin
            a=c_[a1,a2,a3];
            Ai=a.max(axis=1)-a.min(axis=1)
            if i==0:
                A=Ai
            else:
                A=c_[A,Ai]

        #split quads
        flag=sign(A[:,0]+A[:,2]-A[:,1]-A[:,3])

        ne=self.ne; nea=len(self.index_bad_quad);
        self.elnode=r_[self.elnode,ones([nea,4])-3].astype('int');
        for i in arange(nea):
            ind=self.index_bad_quad[i]
            nds=self.elnode[ind,:].copy();
            if flag[i]>=0:
                self.elnode[ind,:]=r_[nds[[0,1,2]],-2]; self.i34[ind]=3
                self.elnode[ne+i,:]=r_[nds[[2,3,0]],-2]
            else:
                self.elnode[ind,:]=r_[nds[[1,2,3]],-2]; self.i34[ind]=3
                self.elnode[ne+i,:]=r_[nds[[3,0,1]],-2]

        self.ne=ne+nea
        self.i34=r_[self.i34,ones(nea)*3].astype('int');
        self.elnode=self.elnode.astype('int')

        #write new grids
        if fname is not None: self.write_hgrid(fname)

    def check_quads(self,angle_min=60,angle_max=120,fname='bad_quad.bp'):
        '''
        1). check the quality of quads, violation when internal angle < angle_min, or >angle_max
        2). the locations of bad quads are saved in file "fname"
        '''

        qind=nonzero(self.i34==4)[0];
        x=self.x[self.elnode[qind,:]]; y=self.y[self.elnode[qind,:]];

        #compute internal angle
        a=[];
        for i in arange(4):
            id1=mod(i-1+4,4); id2=i; id3=mod(i+1,4)
            x1=x[:,id1]; x2=x[:,id2]; x3=x[:,id3];
            y1=y[:,id1]; y2=y[:,id2]; y3=y[:,id3];

            ai=angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2))
            a.append(ai*180/pi);
        a=array(a).T; a=mod(a+360,360)

        #check violation
        for i in arange(4):
            if i==0:
                fp=(a[:,i]<=angle_min)|(a[:,i]>=angle_max)
            else:
                fp=fp|(a[:,i]<=angle_min)|(a[:,i]>=angle_max)

        self.index_bad_quad=qind[nonzero(fp)[0]];

        #output bad_quad location as bp file
        if not hasattr(self,'xctr'): self.compute_ctr()
        qxi=self.xctr[self.index_bad_quad]; qyi=self.yctr[self.index_bad_quad]
        sbp=schism_bpfile(); sbp.nsta=len(qxi); sbp.x=qxi; sbp.y=qyi; sbp.z=zeros(sbp.nsta); sbp.write_bpfile(fname)

    def plot_bad_quads(self,color='r',ms=12,*args):
        #plot grid with bad quads
        if not hasattr(self,'index_bad_quad'): self.check_quads()
        if not hasattr(self,'xctr'): self.compute_ctr()

        qxi=self.xctr[self.index_bad_quad]; qyi=self.yctr[self.index_bad_quad]
        self.plot_grid()
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

    def check_skew_elems(self,angle_min=5,fname='skew_element.bp',fmt=0):
        '''
        1) check schism grid's skewness with angle<=angle_min
        2) the locations of skew elements are (xskew,yskew), and also save in file "fname"
        Inputs:
            angle_min: skew element if one of element's internal angles is smaller than angle_min
            fname=None: not save skew_element.bp; fname!=None: save skew_element.bp
            fmt=1: return indices of skew elements
        '''

        if not hasattr(self,'dpe'): self.compute_ctr()

        #for triangles
        fp=nonzero(self.i34==3)[0]; x=self.x[self.elnode[fp,:3]]; y=self.y[self.elnode[fp,:3]]; sind=[]
        for i in arange(3):
            x1=x[:,i]; x2=x[:,(i+1)%3]; x3=x[:,(i+2)%3]; y1=y[:,i]; y2=y[:,(i+1)%3]; y3=y[:,(i+2)%3]
            ai=abs(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))*180/pi
            sind.extend(nonzero(ai<=angle_min)[0])
        sind3=fp[unique(sind).astype('int')]

        #for quads
        fp=nonzero(self.i34==4)[0]; x=self.x[self.elnode[fp,:]]; y=self.y[self.elnode[fp,:]]; sind=[]
        for i in arange(4):
            x1=x[:,i]; x2=x[:,(i+1)%4]; x3=x[:,(i+2)%4]; y1=y[:,i]; y2=y[:,(i+1)%4]; y3=y[:,(i+2)%4]
            ai=abs(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))*180/pi
            sind.extend(nonzero(ai<=angle_min)[0])
        sind4=fp[unique(sind).astype('int')]; sindw=r_[sind3,sind4]

        #combine and save
        if fname is not None:
            self.xskew,self.yskew,self.zskew=self.xctr[sindw],self.yctr[sindw],self.dpe[sindw]
            sbp=schism_bpfile(); sbp.nsta=len(self.xskew); sbp.x,sbp.y,sbp.z=self.xskew,self.yskew,self.zskew; sbp.write_bpfile(fname)
        if fmt==1: return sindw

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
        npt=len(pxy); sindp=arange(npt)
        if not hasattr(self,'bndinfo'): self.compute_bnd()
        if not hasattr(self.bndinfo,'nb'): self.compute_bnd()
        for i in arange(self.bndinfo.nb):
            fpb=self.bndinfo.ibn[i].astype('int'); fp=inside_polygon(pxy[sindp],self.x[fpb],self.y[fpb])==1
            sindp=sindp[fp] if self.bndinfo.island[i]==0 else sindp[~fp]
            if len(sindp)==0: break
        sind=zeros(npt).astype('int'); sind[sindp]=1
        return sind

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
               if S.nlb<S.nob:
                  bid=S.bid[-1]; pid=nonzero(S.ibn[bid]==S.pt[-1])[0][0]
                  S.nlb=S.nlb+1; ibni=r_[S.ibn[bid][pid:],S.ibn[bid][0]]; S.ilbn.append(ibni)
                  hlb=plot(self.x[ibni],self.y[ibni],'g-'); S.hlb.append(hlb)

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
               if S.npt%2==1 and pid0>=pid: return  #the 2nd pt is ahead of the 1st pt
            if bid not in S.bid: S.ibn[bid]=r_[S.ibn[bid][pid:],S.ibn[bid][:pid]] #reorder boundary points

            #new bnd pt
            S.pt.append(ip); S.bid.append(bid); S.npt=S.npt+1
            hp=plot(self.x[ip],self.y[ip],'ro'); S.hp.append(hp)

            #new open bnd
            if S.npt%2==0:
               S.nob=S.nob+1; ibni=S.ibn[bid][pid0:(pid+1)]; S.iobn.append(ibni)
               hob=plot(self.x[ibni],self.y[ibni],'r-'); S.hob.append(hob)

            #new land bnd
            if S.npt>2 and S.npt%2==1 and bid0==bid:
               S.nlb=S.nlb+1; ibni=S.ibn[bid][pid0:(pid+1)]; S.ilbn.append(ibni)
               hlb=plot(self.x[ibni],self.y[ibni],'g-'); S.hlb.append(hlb)

            #add a new land bnd to the end of the segment
            if S.npt>=2 and bid0!=bid:
               S.nlb=S.nlb+1; ibni=r_[S.ibn[bid0][pid0:],S.ibn[bid0][0]]; S.ilbn.append(ibni)
               hlb=plot(self.x[r_[ibni,S.ibn[bid0][0]]],self.y[r_[ibni,S.ibn[bid0][0]]],'g-'); S.hlb.append(hlb)
            gcf().canvas.draw()

        def remove_pt(x,y):
            if S.npt==0: return
            bid=S.bid[-1]; pid=nonzero(S.ibn[bid]==S.pt[-1])[0][0]

            #remove bnd pt
            S.hp[-1][0].remove(); S.hp.pop(); S.pt.pop(); S.bid.pop(); S.npt=S.npt-1

            #remove open bnd
            if S.npt%2==1: S.hob[-1][0].remove(); S.hob.pop(); S.nob=S.nob-1; S.iobn.pop()

            #remove land bnd
            if (S.nlb>S.nob) or (S.nlb==S.nob and S.npt%2==0 and S.npt>0):
               S.hlb[-1][0].remove(); S.hlb.pop(); S.ilbn.pop(); S.nlb=S.nlb-1
            gcf().canvas.draw()

        #add bnd icon
        if mpl._pylab_helpers.Gcf.get_active() is None: self.plot_grid()
        acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
        abn=acs[nonzero(ats=='bnd')[0][0]] if 'bnd' in ats else gcf().canvas.toolbar.addAction('bnd')

        #add bndinfo capsule
        if not hasattr(self,'bndinfo'): self.bndinfo=zdata()
        S=self.bndinfo; S.hp=[]; S.hob=[]; S.hlb=[]; S.nob=0; S.iobn=[]; S.nlb=0; S.ilbn=[]; S.npt=0; S.pt=[]; S.bid=[]

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

class schism_bpfile:
    def __init__(self):
        self.nsta=0; self.x=array([]); self.y=array([]); self.z=array([]);
        self.station=[]; self.hp=[]; self.ht=[]
        self.edit()

    @property
    def VINFO(self):
        return get_VINFO(self)

    def read_reg(self,fname):
        self.read_bpfile(fname,fmt=1)

    def read_bpfile(self,fname,fmt=0):
        #read file content
        lines=[i.strip().split() for i in open(fname,'r').readlines()]
        stations=[i.strip().split('!')[-1] for i in open(fname,'r').readlines()[2:] if ('!' in i)]
        if fmt==0:
            self.nsta=int(lines[1][0])
            if self.nsta==0: return
            fc=lambda x: x if len(x)==4 else [*x[:4],x[4][1:]]
            data=array([fc(line) for line in lines[2:(2+self.nsta)]])

            self.x=data[:,1].astype(float)
            self.y=data[:,2].astype(float)
            self.z=data[:,3].astype(float)
        elif fmt==1:
            self.nsta=int(lines[2][0])
            if self.nsta==0: return
            data=squeeze(array([lines[3:]])).astype('float')
            self.x=data[:,0]
            self.y=data[:,1]
            self.z=zeros(self.nsta)
        else:
            sys.exit('unknow format')

        #get unique station data.
        if len(stations)==self.nsta:
           self.station=array(stations)
        else:
           self.station=array(['{}'.format(i) for i in arange(self.nsta)])

    def save(self,fname):
        '''
        If fname.endswith('reg'), save points as ACE/gredit *.reg file. Otherwise, save as *.bp file
        '''
        if fname.endswith('.reg'):
           self.write_bpfile(fname,fmt=1)
        else:
           self.write_bpfile(fname,fmt=0)

    def write_reg(self,fname):
        self.write_bpfile(fname,fmt=1)

    def write_bpfile(self,fname,fmt=0):
        '''
        fmt=0: write ACE/gredit *.bp file
        fmt=1: write ACE/gredit *.reg file
        '''

        fid=open(fname,'w+')
        #write header
        if hasattr(self,'note'): fid.write('ACE/gredit: {}'.format(self.note))
        if fmt==0: fid.write('bpfile in ACE/gredit format\n{}\n'.format(self.nsta))
        if fmt==1: fid.write('Region in ACE/gredit format\n1\n{} 1\n'.format(self.nsta))

        #get station names
        stations=[i+1 for i in arange(self.nsta)]
        if hasattr(self,'station') and len(self.station)==self.nsta: stations=self.station

        #write pts
        for i in arange(self.nsta):
            if fmt==0: fid.write('{:<d} {:<.8f} {:<.8f} {:<.8f} !{}\n'.format(i+1,self.x[i],self.y[i],self.z[i],stations[i]))
            if fmt==1: fid.write('{:<.8f} {:<.8f}\n'.format(self.x[i],self.y[i]))
        fid.close()

    def get_unique_pts(self,fmt=0):
        '''
        compute unique pts
            fmt=0: compute ux,uy,uz,ustation of the point
            fmt=1: replace (x,y,z,station) by (ux,uy,uz,ustation)
        '''
        #get unique locations
        upxy,sind=unique(self.x+1j*self.y,return_index=True); sind=sort(sind)
        self.ux=self.x[sind]; self.uy=self.y[sind]
        self.uz=self.z[sind]; self.ustation=self.station[sind]
        if fmt==1: self.x,self.y,self.z,self.station,self.nsta=self.ux,self.uy,self.uz,self.ustation,len(self.ux)
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
        px,py=proj(prj0=prj0,prj1=prj1,x=self.x,y=self.y,lon0=lon0,lat0=lat0)
        if fmt==0: self.x,self.y=px,py
        return [px,py]

    def write_shapefile(self,fname,prj='epsg:4326'):
        self.shp_bp=zdata()
        self.shp_bp.type='POINT'
        self.shp_bp.xy=c_[self.x,self.y]
        self.shp_bp.prj=get_prj_file(prj)

        if hasattr(self,'station'):
            self.shp_bp.attname=['station']
            self.shp_bp.attvalue=self.station
        write_shapefile_data(fname,self.shp_bp)

    def plot_station(self,ax=None,color='r',marker='.',ls=None,label=True,fmt=0,**args):
        '''
        plot points on current figure
          fmt=0: plot all points
          fmt=1: plot unique points
        '''

        self.edit()
        #pre-processing
        if ls is None: ls='None'
        lc=color if label else 'None'
        if not None: ax=gca()
        if fmt==0: sx,sy,sz,stations=self.x,self.y,self.z,self.station
        if fmt==1: sx,sy,sz,stations=self.get_unique_pts()

        #plot
        self.hp=[]; self.ht=[]
        for i,station in enumerate(stations):
            hpi=plot(sx[i],sy[i],marker=marker,color=color,linestyle=ls,**args); self.hp.append(hpi)
            hti=text(sx[i],sy[i],station,color=lc); self.ht.append(hti)
        #show(block=False)
        return [self.hp,self.ht]

    def compute_acor(self,gd):
        #compute areal coordinates, and gd is the schism grid
        self.ie,self.ip,self.acor=gd.compute_acor(c_[self.x,self.y])
        return self.ie,self.ip,self.acor

    def edit(self):
        def connect_actions():
            self.cidmove=gcf().canvas.mpl_connect('motion_notify_event', onmove)
            self.cidpress=gcf().canvas.mpl_connect('button_press_event', onclick)
            if self.nsta!=0 and len(self.hp)==0: self.plot_station()
            acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]); ac=acs[nonzero(ats=='Pan')[0][0]]
            if not ac.isChecked(): ac.trigger()
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
            if dlk==1 and btn==2:
               gcf().canvas.mpl_disconnect(self.cidpress)
               gcf().canvas.mpl_disconnect(self.cidmove)
               acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs]); ac=acs[nonzero(ats=='Pan')[0][0]]
               if ac.isChecked(): ac.trigger()
               gcf().canvas.draw()

        def add_pt(x,y):
            self.nsta=self.nsta+1; self.station=[*self.station,'{}'.format(self.nsta)]
            self.x=r_[self.x,x]; self.y=r_[self.y,y]; self.z=r_[self.z,0.0]

            #plot point
            if len(self.hp)!=0:
                hp=self.hp[-1][0]; ht=self.ht[-1]
                color=hp.get_color(); mk=hp.get_marker(); ms=hp.get_markersize(); ls=hp.get_linestyle()
                fs=ht.get_fontsize(); fw=ht.get_fontweight(); fc=ht.get_color()
            else:
                color='r'; mk='.'; ms=6.0; ls='None'; fs=10; fw='normal'; fc='r'
            hpi=plot(x,y,marker=mk,markersize=ms,color=color,linestyle=ls); self.hp.append(hpi)
            hti=text(x,y,self.station[-1],color=fc,fontsize=fs,fontweight=fw); self.ht.append(hti)
            gcf().canvas.draw()

        def remove_pt(x,y):
            if self.nsta==0: return
            distp=squeeze(abs((self.x-x)+1j*(self.y-y))); sid=nonzero(distp==distp.min())[0][0]
            color='r'; mk='.'; ms=6.0; ls='None'; fs=10; fw='normal'; fc='r'
            for i in arange(sid,self.nsta):
                if i==self.nsta-1:
                   self.hp[-1][0].remove(); self.ht[-1].remove()
                   del self.hp[-1]; del self.ht[-1]
                else:
                   xi=self.x[i+1]; yi=self.y[i+1]
                   self.x[i]=xi; self.y[i]=yi; self.station[i]='{}'.format(i+1)
                   self.hp[i][0].remove(); self.ht[i].remove()
                   hpi=plot(xi,yi,marker=mk,markersize=ms,color=color,linestyle=ls); self.hp[i]=hpi
                   hti=text(xi,yi,self.station[i],color=fc,fontsize=fs,fontweight=fw); self.ht[i]=hti
            self.x=self.x[:-1]; self.y=self.y[:-1]; self.z=self.z[:-1]; self.station=self.station[:-1]; self.nsta=self.nsta-1
            gcf().canvas.draw()

        def move_pt(xi,yi):
            distp=squeeze(abs((self.x-xi)+1j*(self.y-yi))); sid=nonzero(distp==distp.min())[0][0]
            color='r'; mk='.'; ms=6.0; ls='None'; fs=10; fw='normal'; fc='r'
            self.x[sid]=xi; self.y[sid]=yi
            self.hp[sid][0].remove(); self.ht[sid].remove()
            hpi=plot(xi,yi,marker=mk,markersize=ms,color=color,linestyle=ls); self.hp[sid]=hpi
            hti=text(xi,yi,self.station[sid],color=fc,fontsize=fs,fontweight=fw); self.ht[sid]=hti
            gcf().canvas.draw()

        if mpl._pylab_helpers.Gcf.get_active() is not None:
            acs=gcf().canvas.toolbar.actions(); ats=array([i.iconText() for i in acs])
            abp=acs[nonzero(ats=='bp')[0][0]] if 'bp' in ats else gcf().canvas.toolbar.addAction('bp')
            #if not abp.isCheckable(): abp.setCheckable(True)

            #disconnect and clean previous bpfile
            if hasattr(abp,'bp'):
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
    bp=schism_bpfile(); bp.read_bpfile(fname,fmt=fmt)
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

def read_schism_reg(fname):
    '''
    read schism *.reg file created by ACE/gredit
    '''
    return read_schism_bpfile(fname,fmt=1)

def save_schism_grid(fname='grid',path='.',fmt=0):
    '''
    read and save path/{hgrid.gr3,hgrid.ll,vgrid.in}
       fname: save name
       path:  directory whether grids exist
       fmt=0: not save grid's full geometry; fmt=1: save
    '''
    grd=path+'/hgrid.gr3'; grd0=path+'/hgrid.ll'; vrd='{}/vgrid.in'.format(path); S=zdata()
    if os.path.exists(grd):
       gd=read_schism_hgrid(grd)
       if os.path.exists(grd0): gd0=read_schism_hgrid(grd0); gd.lon,gd.lat=gd0.x,gd0.y
       if fmt==1: gd.compute_all(); gd.compute_bnd()
       S.hgrid=gd
    if os.path.exists(vrd): S.vgrid=read_schism_vgrid(vrd)
    if (not hasattr(S,'hgrid')) and (not hasattr(S,'vgrid')): sys.exit('not found: {}, {}'.format(grd,vrd))
    savez(fname,S)
    return S

class schism_vgrid:
    def __init__(self):
        pass

    @property
    def VINFO(self):
        return get_VINFO(self)

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
        compute schism zcor (ivcor=1)
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
    def write_vgrid(self,fname='vgrid.in',fmt=0):
        '''
        write schism vertical grid
            fmt=0: write vgrid.in in latest format of ivcor=1 (one line per lelvel)
            fmt=1: write vgrid.in in old format of ivcor=1    (one line per node)
        '''
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
    compute schism zcor (ivcor=1)
        sigma: sigma cooridinate (dim=[np,nvrt])
        dp: depth at nodes (dim=[np] or [1])
        eta: surface elevation (dim=[np] or [1])
        fmt: output format of zcor
            fmt=0: bottom depths byeond kbp are extended
            fmt=1: bottom depths byeond kbp are nan
        kbp: index of bottom layer (not necessary, just to speed up if provided for ivcor=1)
        method=1 and ivcor=2: return zcor and kbp
        ifix=1 and ivcor=2: using traditional sigma in shallow if error raise
    '''

    if ivcor==1:
        np=sigma.shape[0]
        if not hasattr(dp,'__len__'):  dp=ones(np)*dp
        if not hasattr(eta,'__len__'): eta=ones(np)*eta

        #get kbp
        if kbp is None:
            kbp=array([nonzero(abs(i+1)<1e-10)[0][-1] for i in sigma])

        #thickness of water column
        hw=dp+eta

        #add elevation
        zcor=hw[:,None]*sigma+eta[:,None]
        fpz=hw<0; zcor[fpz]=-dp[fpz][:,None]

        #change format
        if fmt==1:
            for i in arange(np):
                zcor[i,:kbp[i]]=nan
        return zcor
    elif ivcor==2:
        #get dimension of pts
        if not hasattr(dp,'__len__'):
            np=1; dp=array([dp])
        else:
            np=len(dp)
        if not hasattr(eta,'__len__'): eta=ones(np)*eta
        zcor=ones([vd.nvrt,np])*nan

        cs=(1-vd.theta_b)*sinh(vd.theta_f*vd.sigma)/sinh(vd.theta_f)+ \
            vd.theta_b*(tanh(vd.theta_f*(vd.sigma+0.5))-tanh(vd.theta_f*0.5))/2/tanh(vd.theta_f*0.5)
        #for sigma layer: depth<=h_c
        hmod=dp.copy(); fp=hmod>vd.h_s; hmod[fp]=vd.h_s
        fps=hmod<=vd.h_c
        zcor[(vd.kz-1):,fps]=vd.sigma[:,None]*(hmod[fps][None,:]+eta[fps][None,:])+eta[fps][None,:]

        #depth>h_c
        fpc=eta<=(-vd.h_c-(hmod-vd.h_c)*vd.theta_f/sinh(vd.theta_f))
        if sum(fpc)>0:
            if ifix==0: sys.exit('Pls choose a larger h_c: {}'.format(vd.h_c))
            if ifix==1: zcor[(vd.kz-1):,~fps]=eta[~fps][None,:]+(eta[~fps][None,:]+hmod[~fps][None,:])*vd.sigma[:,None]
        else:
            zcor[(vd.kz-1):,~fps]=eta[~fps][None,:]*(1+vd.sigma[:,None])+vd.h_c*vd.sigma[:,None]+cs[:,None]*(hmod[~fps]-vd.h_c)

        #for z layer
        kbp=-ones(np).astype('int'); kbp[dp<=vd.h_s]=vd.kz-1
        fpz=dp>vd.h_s; sind=nonzero(fpz)[0]
        for i in sind:
            for k in arange(0,vd.kz-1):
                if (-dp[i]>=vd.ztot[k])*(-dp[i]<=vd.ztot[k+1]):
                    kbp[i]=k;
                    break
            #check
            if kbp[i]==-1:
                sys.exit('can not find a bottom level for node')
            elif kbp[i]<0 or kbp[i]>=(vd.kz-1):
                sys.exit('impossible kbp,kz: {}, {}'.format(kbp[i],vd.kz))

            #assign values
            zcor[kbp[i],i]=-dp[i]
            for k in arange(kbp[i]+1,vd.kz-1):
                zcor[k,i]=vd.ztot[k]
        zcor=zcor.T; vd.kbp=kbp

        #change format
        if fmt==0:
            for i in arange(np):
                zcor[i,:kbp[i]]=zcor[i,kbp[i]]
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
    tr=mpl.tri.Triangulation(x,y); gd=schism_grid()
    gd.np,gd.ne=np,len(tr.triangles); gd.x,gd.y,gd.dp=x,y,z
    gd.elnode=c_[tr.triangles,-2*ones([gd.ne,1])].astype('int'); gd.i34=3*ones(gd.ne).astype('int')

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

    #find max/min side or angle values
    angles,sides=[],[];  fp3=gd.i34==3; fp4=gd.i34==4
    id1,id2,id3=ones([3,gd.ne]).astype('int'); sid=arange(gd.ne)
    for i in arange(4):
        id1[fp3]=i%3; id2[fp3]=(i+1)%3; id3[fp3]=(i+2)%3
        id1[fp4]=i%4; id2[fp4]=(i+1)%4; id3[fp4]=(i+2)%4
        x1=gd.x[gd.elnode[sid,id1]]; x2=gd.x[gd.elnode[sid,id2]]; x3=gd.x[gd.elnode[sid,id3]]
        y1=gd.y[gd.elnode[sid,id1]]; y2=gd.y[gd.elnode[sid,id2]]; y3=gd.y[gd.elnode[sid,id3]]
        ai=abs(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))*180/pi; angles.append(ai)
        si=abs((x1+1j*y1)-(x2+1j*y2)); sides.append(si)
    angles=array(angles).T; sides=array(sides)
    mangle=angles.min(axis=1); sidemin=sides.min(axis=0); sidemax=sides.max(axis=0)

    #filter illegal elements
    gd.compute_area(); gd.compute_ctr()
    fangle=nonzero(mangle<angle_min)[0] if (angle_min is not None) else array([])
    farea=nonzero(gd.area>area_max)[0] if (area_max is not None) else array([])
    fside_max=nonzero(sidemax>side_max)[0] if (side_max is not None) else array([])
    fside_min=nonzero(sidemin<side_min)[0] if (side_min is not None) else array([])
    sindp=r_[fangle,farea,fside_max,fside_min].astype('int')

    #filter elements inside region
    if (reg_in is not None) and len(sindp)!=0:
        bp=read_schism_bpfile(reg_in,fmt=1)
        fpr=inside_polygon(c_[gd.xctr[sindp],gd.yctr[sindp]],bp.x,bp.y)==1; sindp=sindp[fpr]

    #filter elements outside region
    if reg_out is not None:
        bp=read_schism_bpfile(reg_out,fmt=1)
        sindo=nonzero(inside_polygon(c_[gd.xctr,gd.yctr],bp.x,bp.y)==0)[0]; sindp=r_[sindp,sindo]

    sind=setdiff1d(arange(gd.ne),sindp)

    #add back elements with dangling pts
    ips=setdiff1d(arange(gd.np),unique(gd.elnode[sind].ravel()))
    if len(ips)!=0:
        gd.compute_nne(); sinde=[]
        for ip in ips:
            ies=gd.indel[ip]
            if method==0: ai=sidemax[ies]; sinde.append(ies[nonzero(ai==min(ai))[0][0]])
            if method==1: ai=mangle[ies]; sinde.append(ies[nonzero(ai==max(ai))[0][0]])
        sind=sort(r_[sind,array(sinde)])

    #delete elements
    gd.ne,gd.i34,gd.elnode=len(sind),gd.i34[sind],gd.elnode[sind]
    gd.area,gd.xctr,gd.yctr,gd.dpe=gd.area[sind],gd.xctr[sind],gd.yctr[sind],gd.dpe[sind]
    return gd

def combine_icm_output(rundir='.',sname='icm.nc',fmt=0,outfmt=0):
    '''
    combine schism-icm station outputs
      rundir: run directory
      sname: save name for combined output
      fmt=1: copy subdomain output file, and then read
      outfmt=0: float32;  outfmt=1: float64
    '''

    outdir=rundir+'/outputs/'
    bp=read_schism_bpfile(rundir+'/istation.in')
    fid=Dataset(outdir+sname,'w',format='NETCDF4'); fvar=fid.variables; fdim=fid.dimensions

    #get fnames information
    fnames=[outdir+i for i in os.listdir(outdir) if (i.startswith('icm_') and i.endswith('.nc'))]; nts=[]
    if fmt==1: [copyfile(i,i+'.copy') for i in fnames]; fnames=[i+'.copy' for i in fnames] #copy output files
    for i in fnames: C=ReadNC(i,1); nts.append(C.variables['time'].size); C.close() #time length

    #combine station output
    for n,fname in enumerate(fnames):
        C=ReadNC(fname,1); cvar=C.variables; cdim=C.dimensions; sind=array(cvar['istation'][:])-1
        if n==0:
           #def dim
           nt=min(nts)
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
                  fid.createVariable(cn,np.float32,cdn,fill_value=False)
               else:
                  fid.createVariable(cn,cvar[cn].dtype,cdn,fill_value=False)
           fid.createVariable('station',str,['nstation'],fill_value=False)
           fvar['time'][:]=array(cvar['time'][:nt])/86400 #set time

        #set variables
        for i,cn in enumerate(cvar):
            cdn=[*cvar[cn].dimensions]; cds=[cdim[k].size for k in cdn]
            if cdn[0]=='nstation': fvar[cn][sind]=array(cvar[cn][:])
            if len(cds)==3 and cds[0]==1: fvar[cn][sind]=array(cvar[cn][0,:,:nt])
            if len(cds)==3 and cds[0]!=1: fvar[cn][sind]=array(cvar[cn][...,:nt]).transpose([1,0,2])
        fvar['station'][:]=bp.station
        C.close()
        if fmt==1: os.remove(fname)
    fid.close()

def combine_schism_hotstart(outdir='.',fmt=0,irec=None):
    '''
    combine schism hotstart
       outdir: schism output directory
       fmt=0: combine the last hotstart; fmt=1: combine all hotstart
       irec: step number of hotstrat (fmt is disabled when irec is set)
    '''

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
        fid=Dataset(outdir+'/'+fname,'w',format='NETCDF4');  #open file
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

def read_schism_output(run,varname,xyz,stacks=None,ifs=0,nspool=1,sname=None,fname=None,grid=None,fmt=0,extend=0):
    '''
    extract time series of SCHISM results @xyz or transects @xy (works for scribe IO and combined oldIO)
       run:     run directory where (grid.npz or hgrid.gr3) and outputs are located
       varname: variables to be extracted; accept shortname(s) or fullname(s) (elev, hvel, horizontalVelX, NO3, ICM_NO3, etc. )
       xyz:     c_[x,y,z], or bpfile, or c_[x,y]
       fmt:     (optional) 0: read time series @xyz;     1: read transect @xy
       stacks:  (optional) output stacks to be extract; all avaiable stacks will be extracted if not specified
       ifs=0:   (optional) extract results @xyz refers to free surface (default); ifs=1: refer to fixed levels
       nspool:  (optional) sub-sampling frequency within each stack (npsool=1 means all records)
       sname:   (optional) variable name for save
       fname:   (optional) save the results as fname.npz
       grid:    (optional) grid=read_schism_hgrid('hgrid.gr3'); grid.compute_all(); used to speed up
       extend:  (optional) 0: extend bottom value beyond;  1: assign nan for value beyond bottom
    '''

    #get schism outputs information
    bdir=run+'/outputs'; modules,outfmt,dstacks,dvars,dvars_2d =get_schism_output_info(bdir,1)

    #variables
    if isinstance(varname,str): varname=[varname]
    if isinstance(sname,str): sname=[sname]

    #read grid
    if grid is not None: gd=grid
    if (grid is None) and os.path.exists(run+'/grid.npz'): gd=loadz(run+'/grid.npz').hgrid
    if (grid is None) and (not os.path.exists(run+'/grid.npz')): gd=read_schism_hgrid(run+'/hgrid.gr3')

    #read station coordinates (xyz)
    if isinstance(xyz,str): bp=read_schism_bpfile(xyz); xyz=c_[bp.x,bp.y,bp.z]
    if array(xyz).ndim==1: xyz=array(xyz)[None,:]
    lx,ly=xyz.T[:2]; npt=len(lx); lz=xyz.T[2] if xyz.shape[1]==3 else zeros(npt)
    pie,pip,pacor=gd.compute_acor(c_[lx,ly]); pip,pacor=pip.T,pacor.T

    #extract time series@xyz
    S=zdata(); mtime=[]; mdata=[[] for i in varname]
    stacks=dstacks if (stacks is None) else [*array(stacks).ravel()] #check outputs stacks
    for istack in stacks:
        print('reading stack: {}'.format(istack))
        if outfmt==0:
           if fmt==0: Z0=ReadNC('{}/zCoordinates_{}.nc'.format(bdir,istack),1); Z=Z0.variables['zCoordinates']
           C0=ReadNC('{}/out2d_{}.nc'.format(bdir,istack),1)
        else:
           C0=ReadNC('{}/schout_{}.nc'.format(bdir,istack),1); Z=C0.variables['zcor']
        mti0=array(C0.variables['time'])/86400; nt=len(mti0); mti=mti0[::nspool]; nrec=len(mti); mtime.extend(mti)
        if istack==stacks[0]: nvrt=C0.dimensions['nSCHISM_vgrid_layers'].size
        zii=None

        for m,varnamei in enumerate(varname):
            vs=[]; svars=get_schism_var_info(varnamei,modules,fmt=outfmt) #get variable information
            for n,[vari,svar] in enumerate(svars):
                if svar in dvars_2d: #2D
                    C=C0.variables[svar]; np=C.shape[1]
                    if np==gd.np: vi=array([(array([C[i,j] for j in pip])*pacor).sum(axis=0) for i in arange(0,nt,nspool)])
                    if np==gd.ne: vi=array([C[i,pie] for i in arange(0,nt,nspool)])
                    vs.append(vi)
                else: #3D
                    #read zcoor,and extend zcoor downward
                    if (zii is None) and fmt==0:
                       if n==0: zii=(array([[array(Z[i,:,k])[pip] for k in arange(nvrt)] for i in arange(0,nt,nspool)])*pacor[None,None,...]).sum(axis=2).transpose([1,0,2])
                       for k in arange(nvrt-1): z1=zii[nvrt-k-2]; z2=zii[nvrt-k-1]; z1[abs(z1)>1e8]=z2[abs(z1)>1e8]
                       if ifs==0: zii=zii-zii[-1][None,...]

                    #read data for the whole vertical
                    C1=ReadNC('{}/{}_{}.nc'.format(bdir,svar,istack),1) if outfmt==0 else C0
                    C=C1.variables[svar]; np=C.shape[1]; nd=C.ndim
                    if np==gd.np and nd==3: vii=(array([[array(C[i,:,k])[pip]  for k in arange(nvrt)] for i in arange(0,nt,nspool)])*pacor[None,None,...]).sum(axis=2)
                    if np==gd.np and nd==4: vii=(array([[array(C[i,:,k])[pip]  for k in arange(nvrt)] for i in arange(0,nt,nspool)])*pacor[None,None,...,None]).sum(axis=2)
                    if np==gd.ne: vii=array([[C[i,pie,k] for k in arange(nvrt)] for i in arange(0,nt,nspool)])
                    vii=vii.transpose([1,0,2]) if nd==3 else vii.transpose([1,0,2,3])
                    if extend==0:
                       for k in arange(nvrt-1): z1=vii[nvrt-k-2]; z2=vii[nvrt-k-1]; z1[abs(z1)>1e8]=z2[abs(z1)>1e8] #extend value at bottom
                    else:
                       fpn=abs(vii)>1e8; vii[fpn]=nan
                    if outfmt==0: C1.close()

                    #interp in the vertical
                    if fmt==0: #time series
                       for ivs in arange(2):
                           if ivs==1 and nd==3: continue
                           viii=vii if nd==3 else vii[:,:,:,ivs]
                           vi=ones([nrec,npt]); zm=-tile(lz,[nrec,1]); dz=1e-10
                           ziii=zii[-1]-dz; fpz=zm>ziii; zm[fpz]=ziii[fpz]
                           ziii=zii[0] +dz; fpz=zm<ziii; zm[fpz]=ziii[fpz]
                           for k in arange(nvrt-1):
                               z1=zii[k]; z2=zii[k+1]; v1=viii[k]; v2=viii[k+1]; fpz=(zm>z1)*(zm<=z2)
                               if sum(fpz)!=0: vi[fpz]=v1[fpz]+(v2[fpz]-v1[fpz])*(zm[fpz]-z1[fpz])/(z2[fpz]-z1[fpz])
                           vs.append(vi)
                    else: #transect
                       if nd==3: vs.append(vii)
                       if nd==4: vs.append(vii[...,0]); vs.append(vii[...,1])
            mdata[m].extend(array(vs).transpose([1,2,0])) if (fmt==0 or (svar in dvars_2d)) else mdata[m].extend(array(vs).transpose([2,3,1,0]))
        C0.close()
        if fmt==0 and outfmt==0: Z0.close()

    #save data
    if sname is None: sname=varname
    for m,k in enumerate(sname):
        if array(mdata[m]).ndim==3: S.__dict__[k]=squeeze(array(mdata[m]).transpose([1,0,2]))
        if array(mdata[m]).ndim==4: S.__dict__[k]=squeeze(array(mdata[m]).transpose([1,0,2,3]))
    S.time=array(mtime)
    if fname is not None: savez(fname,S)
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
          stacks=unique([int(i[:-3].split('_')[-1]) for i in os.listdir(run) if (i.startswith('out2d_') and i.endswith('.nc'))])
          C=ReadNC(run+'/out2d_{}.nc'.format(stacks[0]),1); svars_2d=setdiff1d(array([*C.variables]),dvars); C.close()
          svars=setdiff1d(r_[ovars,svars_2d],dvars)
       else:
          stacks=unique([int(i[:-3].split('_')[-1]) for i in os.listdir(run) if (i.startswith('schout_') and i.endswith('.nc'))])
          fname=[i for i in os.listdir(run) if (i.startswith('schout_') and i.endswith('_{}.nc'.format(stacks[0])))][0]
          C=ReadNC(run+'/{}'.format(fname),1); svars=setdiff1d(array([*C.variables]),dvars)
          svars_2d=array([i for i in svars if C.variables[i].ndim==2]); C.close()

       #get SCHISM modules
       M=get_schism_var_info(fmt=outfmt).__dict__
       modules=array([k for k in M if len(intersect1d([*M[k].values()],svars))!=0])
       return [modules,outfmt,stacks,svars,svars_2d]
    else:
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

    if fmt==1: Hydro=oHydro
    mdict={'Hydro':Hydro,'ICM':ICM,'COS':COS,'SED':SED}; S=zdata(); C=[]
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

if __name__=="__main__":
    pass
