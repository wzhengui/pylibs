#/usr/bin/env python3
from pylib import *

class schism_grid:
    def __init__(self):
        pass

    def plot_grid(self,ax=None,method=0,fmt=0,value=None,mask=None,ec=None,fc=None,
             lw=0.1,levels=None,ticks=None,clim=None,extend='both',cb=True,**args):
        '''
        plot grid with default color value (grid depth)
        method=0: using tricontourf; method=1: using PolyCollection (old method)
        fmt=0: plot grid only; fmt=1: plot color
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
           elif fmt==1:
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
                 vmin,vmax=min(value),max(value)
              else:
                 vmin,vmax=clim

              #detemine levels
              if levels is None: levels=51
              if squeeze(array([levels])).size==1:
                 levels=linspace(vmin,vmax,int(levels))

              #set mask
              if sum(isnan(value))!=0: tri=tri[~isnan(value[tri].sum(axis=1))]

              if vmin==vmax:
                 hg=tricontourf(self.x,self.y,tri,value,vmin=vmin,vmax=vmax,extend=extend,**args)
              else:
                 hg=tricontourf(self.x,self.y,tri,value,levels=levels,vmin=vmin,vmax=vmax,extend=extend,**args)

              #add colobar
              cm.ScalarMappable.set_clim(hg,vmin=vmin,vmax=vmax)
              if cb:
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

           self.hg=hg; #show(block=False)
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
           return hg

    def compute_bnd(self):
        pass

    def create_bnd(self,figsize=[8,9]):

        #compute xm and ym
        xm=[self.x.min(),self.x.max()]; ym=[self.y.min(),self.y.max()]
        dx,dy=0.01*diff(xm),0.01*diff(ym); xm=[xm[0]-dx,xm[1]+dx]; ym=[ym[0]-dy,ym[1]+dy]

        #add all bnd pts
        self.binfo=npz_data()
        bind=[]; [bind.extend(i) for i in self.iobn]; [bind.extend(i) for i in self.ilbn]
        self.binfo.bind=array(bind); self.binfo.x=self.x[self.binfo.bind]; self.binfo.y=self.y[self.binfo.bind]
        self.binfo.obp=array([]).astype('int');  self.binfo.lbp=array([]).astype('int'); self.binfo.hp=[]

        #plot grid
        figure(figsize=[8,9])
        ax0=axes([0.01,0.01,0.98,0.94])
        self.plot_bnd(marker='.',ms=3,color='k')
        setp(gca(),xticklabels=[],yticklabels=[],xlim=xm,ylim=ym)
        move_figure(gcf(),0,0)

        #add active buttion
        def add_pt_open_bnd(*args,gd=self,ax=ax0):

            sca(ax)
            while True:
                xy0=ginput(1)
                if len(xy0)==0:
                    break
                else:
                    xi0,yi0=xy0[0]

                #find the nearest pts
                dist=abs((gd.binfo.x-xi0)+1j*(gd.binfo.y-yi0)); sind=nonzero(dist==min(dist))[0][0]
                sid=gd.binfo.bind[sind]; xi=gd.x[sid]; yi=gd.y[sid]

                gd.binfo.obp=r_[gd.binfo.obp,sid]

                #plot pts
                hp=plot(xi,yi,'r.',ms=6)
                show(); gd.binfo.hp.append(hp[0])

        def remove_pt_open_bnd(*args,gd=self,ax=ax0):
            if len(gd.binfo.obp)!=0:
                gd.binfo.obp=gd.binfo.obp[:-1]
                sca(ax); gd.binfo.hp[-1].remove(); gd.binfo.hp.pop()

        ax1=axes([0.02,0.96,0.2,0.03]); ax2=axes([0.7,0.96,0.25,0.03])
        hb1=Button(ax1,'add open boundary',color='r'); hb1.on_clicked(add_pt_open_bnd)
        hb2=Button(ax2,'remove boundary pts',color='gray'); hb2.on_clicked(remove_pt_open_bnd)

        sca(ax0); self.ha=gca(); self.hb=[hb1,hb2]
        return

    def plot_bnd(self,c='k',lw=1,ax=None,**args):
        '''
          plot schims grid boundary

          gd.plot_bnd(): plot bnd
          gd.plot_bnd(c='rb'): open bnd in red,land bnd in blue

        '''
        if ax!=None: sca(ax)
        if len(c)==1: c=c*3

        #get indices for bnds
        sindo=[]
        for i in arange(self.nob):
            sindo=r_[sindo,-1,self.iobn[i]]
        sindo=array(sindo).astype('int'); fpn=sindo==-1
        bx1=self.x[sindo]; by1=self.y[sindo]
        bx1[fpn]=nan; by1[fpn]=nan

        sindl=[]
        for i in arange(self.nlb):
            sindl=r_[sindl,-1,self.ilbn[i]]
        sindl=array(sindl).astype('int'); fpn=sindl==-1
        bx2=self.x[sindl]; by2=self.y[sindl]
        bx2[fpn]=nan; by2[fpn]=nan

        hb1=plot(bx1,by1,c[0],lw=lw,**args)
        hb2=plot(bx2,by2,c[-1],lw=lw,**args)
        #show(block=False)
        self.hb=[hb1,hb2]

    def read_hgrid(self,fname,*args):
        with open(fname,'r') as fid:
            lines=fid.readlines()

        #read ne and np
        num=array(lines[1].split()[0:2]).astype('int')
        self.ne=num[0]; self.np=num[1]

        #read lx,ly and dp
        num=[]
        for i in arange(self.np):
            num.append(array(lines[2+i].split()[1:4]));
        num=array(num).astype('float')
        self.x=num[:,0]
        self.y=num[:,1]
        self.dp=num[:,2]

        if len(lines)<(2+self.np+self.ne):
            return

        #read elnode and i34
        num=[]
        for i in arange(self.ne):
            num.append(lines[2+self.np+i].split())
        num=array([s if len(s)==6 else [*s,'-1'] for s in num])
        num=num.astype('int')

        self.i34=num[:,1]
        self.elnode=num[:,2:]-1

        #compute ns
        self.compute_ns()

        if len(lines)<(4+self.np+self.ne):
            return

        #read obnd info
        n=2+self.np+self.ne; num=array(lines[n].split()[0]).astype('int'); n=n+2;
        self.nob=num

        self.nobn=[];
        self.iobn=[];
        for i in arange(self.nob):
            num=array(lines[n].split()[0]).astype('int')
            self.nobn.append(num)
            num=[];
            for m in arange(self.nobn[i]):
                num.append(lines[n+m+1].split()[0])
            self.iobn.append(array(num).astype('int')-1)
            n=n+self.nobn[i]+1;
        self.nobn=array(self.nobn);
        self.iobn=array(self.iobn,dtype='O')
        if len(self.iobn)==1: self.iobn=self.iobn.astype('int')


        #read lbnd info
        num=array(lines[n].split()[0]).astype('int'); n=n+2;
        self.nlb=num

        self.nlbn=[];
        self.ilbn=[];
        self.island=[];
        for i in arange(self.nlb):
            #if 'island' in lines[n]:
            #    self.island.append(1)
            #else:
            #    self.island.append(0)
            num=array(lines[n].split()[0]).astype('int')
            self.nlbn.append(num)
            num=[];
            for m in arange(self.nlbn[i]):
                num.append(lines[n+m+1].split()[0])
            self.ilbn.append(array(num).astype('int')-1)

            #update island flag method
            if num[0]==num[-1]:
                self.island.append(1)
            else:
                self.island.append(0)

            n=n+self.nlbn[i]+1;
        self.island=array(self.island);
        self.nlbn=array(self.nlbn);
        self.ilbn=array(self.ilbn,dtype='O');
        if len(self.ilbn)==1: self.ilbn=self.ilbn.astype('int')

    def read_prop(self,fname):
        '''
        read schism prop, and return the values
        '''
        pdata=loadtxt(fname); 
        pvalue=pdata[:,1] if pdata.ndim==2 else pdata[None,:][:,1]

        return pvalue

    def interp_node_to_elem(self,value=None):
        '''
        interpolate node values to element values
            default is self.dp => self.dpe
        '''
        #get node value
        if value is None:
            dp=self.dp
        else:
            dp=value

        #interpolate
        fp1=self.i34==3; fp2=self.i34==4;
        dpe=zeros(self.ne)*nan
        dpe[fp1]=mean(dp[self.elnode[fp1,0:3]],axis=1)
        dpe[fp2]=mean(dp[self.elnode[fp2,:]],axis=1)

        return dpe

    def interp_elem_to_node(self,value=None,method=0,p=1):
        '''
        interpolate element values to nodes
        if value not given, dpe is used
        method=0: simple avarage; method=1: inverse distance (power=p)
        '''
        #-specify element values
        if value is None:
            if not hasattr(self,'dpe'): self.compute_ctr()
            v0=self.dpe
        else:
            v0=value;

        #compute node ball
        if not hasattr(self,'nne'): self.compute_node_ball()

        #interpolation
        v=[];
        for i in arange(self.np):
            ind=self.ine[i];
            if method==0: #aveaging
                vi=sum(v0[ind])/self.nne[i];
            else: #inverse distance
                W=1/((self.xctr[ind]-self.x[i])**2+(self.yctr[ind]-self.y[i])**2)**(p/2); #weight
                vi=sum(W*v0[ind])/sum(W)
            v.append(vi)
        v=array(v)
        return v

    def compute_ctr(self):
        '''
        compute element center XYZ
        '''
        if not hasattr(self,'xctr'):
           fp1=self.i34==3; fp2=self.i34==4;
           self.xctr=zeros(self.ne)*nan
           self.yctr=zeros(self.ne)*nan

           self.xctr[fp1]=mean(self.x[self.elnode[fp1,0:3]],axis=1)
           self.yctr[fp1]=mean(self.y[self.elnode[fp1,0:3]],axis=1)
           self.xctr[fp2]=mean(self.x[self.elnode[fp2,:]],axis=1)
           self.yctr[fp2]=mean(self.y[self.elnode[fp2,:]],axis=1)

        self.dpe=self.interp_node_to_elem()
        return self.dpe

    def compute_area(self):
        fp=self.elnode[:,-1]<0;
        x1=self.x[self.elnode[:,0]]; y1=self.y[self.elnode[:,0]];
        x2=self.x[self.elnode[:,1]]; y2=self.y[self.elnode[:,1]];
        x3=self.x[self.elnode[:,2]]; y3=self.y[self.elnode[:,2]];
        x4=self.x[self.elnode[:,3]]; y4=self.y[self.elnode[:,3]]; x4[fp]=x1[fp]; y4[fp]=y1[fp]
        self.area=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)+(x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
        return self.area

    def compute_gradient(self):
        if not hasattr(self,'area'): self.compute_area()
        if not hasattr(self,'dpe'): self.compute_ctr()
        #get pts
        fp=self.elnode[:,-1]<0; fpn=~fp;
        x1=self.x[self.elnode[:,0]]; y1=self.y[self.elnode[:,0]]; v1=self.dp[self.elnode[:,0]]
        x2=self.x[self.elnode[:,1]]; y2=self.y[self.elnode[:,1]]; v2=self.dp[self.elnode[:,1]]
        x3=self.x[self.elnode[:,2]]; y3=self.y[self.elnode[:,2]]; v3=self.dp[self.elnode[:,2]]
        x4=self.x[self.elnode[:,3]]; y4=self.y[self.elnode[:,3]]; v4=self.dp[self.elnode[:,3]]
        x4[fp]=x1[fp]; y4[fp]=y1[fp]; v4[fp]=v1[fp]
        a1=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2
        a2=((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2

        #compute gradients
        self.dpedx=(v1*(y2-y3)+v2*(y3-y1)+v3*(y1-y2))/(2*a1)
        self.dpedy=((x3-x2)*v1+(x1-x3)*v2+(x2-x1)*v3)/(2*a1)
        self.dpedxy=sqrt(self.dpedx**2+self.dpedy**2)

        #modify quads
        dpedx2=(v1[fpn]*(y3[fpn]-y4[fpn])+v3[fpn]*(y4[fpn]-y1[fpn])+v4[fpn]*(y1[fpn]-y3[fpn]))/(2*a2[fpn])
        dpedy2=((x4[fpn]-x3[fpn])*v1[fpn]+(x1[fpn]-x4[fpn])*v3[fpn]+(x3[fpn]-x1[fpn])*v4[fpn])/(2*a2[fpn])
        dpedxy2=sqrt(dpedx2**2+dpedy2**2)

        self.dpedx[fpn]=(self.dpedx[fpn]+dpedx2)/2
        self.dpedy[fpn]=(self.dpedy[fpn]+dpedy2)/2
        self.dpedxy[fpn]=(self.dpedxy[fpn]+dpedxy2)/2

        #get node value------
        self.dpdx=self.interp_elem_to_node(value=self.dpedx)
        self.dpdy=self.interp_elem_to_node(value=self.dpedy)
        self.dpdxy=self.interp_elem_to_node(value=self.dpedxy)

        return self.dpdx,self.dpdy,self.dpdxy

    def compute_ns(self):
        '''
        compute number of hgrid sides
        '''
        #collect sides
        fp3=self.i34==3; fp4=self.i34==4; sind=array([[],[]]).astype('int').T
        for i in arange(3): sind=r_[sind,c_[self.elnode[fp3,mod(i+3,3)],self.elnode[fp3,mod(i+4,3)]]]
        for i in arange(4): sind=r_[sind,c_[self.elnode[fp4,mod(i+4,4)],self.elnode[fp4,mod(i+5,4)]]]

        #sort side
        sind=sort(sind,axis=1).T; sid=unique(sind[0]+1j*sind[1]); self.ns=len(sid)
        return self.ns

    def compute_node_ball(self):
        '''
        compute nodal ball information
        '''
        nne=zeros(self.np).astype('int');
        ine=[[] for i in arange(self.np)];
        for i in arange(self.ne):
            inds=self.elnode[i,:self.i34[i]];
            nne[inds]=nne[inds]+1
            [ine[indi].append(i) for indi in inds]
        self.nne=nne
        self.ine=array([array(ine[i]) for i in arange(self.np)],dtype='O');
        return self.nne, self.ine

    def compute_acor(self,pxy):
        '''
        compute acor coodinate for points pxy[npt,2]

        usage: ie,ip,acor=compute_acor(c_[xi,yi]), where xi and yi are array of coordinates
        outputs: ie[npt],ip[npt,3],acor[npt,3]
               ie:  the element number
               ip:  the nodal indices of the ie
               acor: the area coordinate
        '''
        #compute parents element
        pie,pip=self.inside_grid(pxy,fmt=1); fpn=pie!=-1

        #cooridate for triangles, and areas
        x1,x2,x3=self.x[pip[fpn]].T; y1,y2,y3=self.y[pip[fpn]].T; x,y=pxy[fpn].T
        A=signa(c_[x1,x2,x3],c_[y1,y2,y3]); A1=signa(c_[x,x2,x3],c_[y,y2,y3])
        A2=signa(c_[x1,x,x3],c_[y1,y,y3]);  #A3=signa(c_[x1,x2,x],c_[y1,y2,y])

        #area coordinate
        pacor=zeros(pip.shape); pacor[fpn]=c_[A1/A,A2/A,1-(A1+A2)/A]
        return [pie,pip,pacor]

    def interp(self,xyi):
        #interpolate to get depth at xyi
        ie,ip,acor=self.compute_acor(xyi)
        dpi=(self.dp[ip]*acor).sum(axis=1)
        return dpi

    def write_hgrid(self,fname,value=None,elnode=1,bndfile=None,Info=None):
        '''
        write *.gr3 file
            fname: file name
            value: depth value to be outputted
                   value=const: uniform value in space
                   value=dp[np]: specify depth value
                   value=None:  grid's default depth self.dp is used
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

        #write *gr3
        with open(fname,'w+') as fid:
            fid.write('!grd info:{}\n'.format(Info))
            fid.write('{} {}\n'.format(self.ne,self.np))
            for i in arange(self.np):
                fid.write('{:<d} {:<.8f} {:<.8f} {:<.8f}\n'.format(i+1,self.x[i],self.y[i],dp[i]))
            if elnode!=0:
                for i in arange(self.ne):
                    if self.i34[i]==3: fid.write('{:<d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))
                    if self.i34[i]==4: fid.write('{:<d} {:d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))

            #write bnd information
            if bndfile is not None: fid.writelines(open(bndfile,'r').readlines())

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

    def split_quads(self,angle_min=60,angle_max=120,fname='new.gr3'):
        '''
        1). split the quads that have angle (<angle_min or >angle_max), add append the connectivity in the end
        2). output a new grid "fname"
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
        self.write_hgrid(fname)


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
        #with open("{}".format(fname),'w+') as fid:
        #    for i in arange(len(qxi)):
        #        fid.write('{} {} 0\n'.format(qxi[i],qyi[i]))

    def plot_bad_quads(self,color='r',ms=12,*args):
        #plot grid with bad quads
        if not hasattr(self,'index_bad_quad'): self.check_quads()
        if not hasattr(self,'xctr'): self.compute_ctr()

        qxi=self.xctr[self.index_bad_quad]; qyi=self.yctr[self.index_bad_quad]
        self.plot_grid()
        plot(qxi,qyi,'.',color=color,ms=ms,*args)
        #show(block=False)
        pass

    def proj(self,prj0,prj1='epsg:4326',x=None,y=None,lon0=None,lat0=None):
        '''
        transform the projection of schism grid's coordinates
        Inputs:
            prj0: projection name of schism grid
            prj1: target projection name; default is 'epsg:4326'
            x,y: values of grid coordiantes; default is (gd.x, gd.y)
            lon0,lat0: lon&lat of cpp projection center; needed only if 'cpp' in [prj0,prj1]
                       if ("ll"=>"cpp") and (lon0 or lat0 is not provided): lon0=mean(x); lat0=mean(y)
        '''
        if (x is None) or (y is None): x=self.x; y=self.y
        x1,y2=proj(prj0=prj0,prj1=prj1,x=x,y=y,lon0=lon0,lat0=lat0)
        return [x1,y2]

    def check_skew_elems(self,angle_min=5,fname='skew_element.bp'):
        '''
        1) check schism grid's skewness with angle<=angle_min
        2) the locations of skew elements are (xskew,yskew), and also save in file "fname"
        '''

        if not hasattr(self,'dpe'): self.compute_ctr()

        #for triangles
        fp=self.i34==3; x=self.x[self.elnode[fp,:3]]; y=self.y[self.elnode[fp,:3]]; xctr=self.xctr[fp]; yctr=self.yctr[fp]; zctr=self.dpe[fp]
        sind=[]
        for i in arange(3):
            id1=i; id2=(i+1)%3; id3=(i+2)%3
            x1=x[:,id1]; x2=x[:,id2]; x3=x[:,id3]
            y1=y[:,id1]; y2=y[:,id2]; y3=y[:,id3]
            ai=abs(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))*180/pi
            sindi=nonzero(ai<=angle_min)[0]
            if len(sindi)!=0: sind.extend(sindi)
        sind=array(sind)
        if len(sind)!=0:
            XS3=xctr[sind]; YS3=yctr[sind]; ZS3=zctr[sind]
        else:
            XS3=array([]); YS3=array([]); ZS3=array([])

        #for quads
        fp=self.i34==4; x=self.x[self.elnode[fp,:]]; y=self.y[self.elnode[fp,:]]; xctr=self.xctr[fp]; yctr=self.yctr[fp]; zctr=self.dpe[fp]
        sind=[]
        for i in arange(4):
            id1=i; id2=(i+1)%4; id3=(i+2)%4
            x1=x[:,id1]; x2=x[:,id2]; x3=x[:,id3]
            y1=y[:,id1]; y2=y[:,id2]; y3=y[:,id3]
            ai=abs(angle((x1-x2)+1j*(y1-y2))-angle((x3-x2)+1j*(y3-y2)))*180/pi
            sindi=nonzero(ai<=angle_min)[0]
            if len(sindi)!=0: sind.extend(sindi)
        sind=array(sind)
        if len(sind)!=0:
            XS4=xctr[sind]; YS4=yctr[sind]; ZS4=zctr[sind]
        else:
            XS4=array([]); YS4=array([]); ZS4=array([])

        #combine and save
        self.xskew=r_[XS3,XS4]; self.yskew=r_[YS3,YS4]; zskew=r_[ZS3,ZS4]
        sbp=schism_bpfile(); sbp.nsta=len(self.xskew); sbp.x=self.xskew; sbp.y=self.yskew; sbp.z=zskew; sbp.write_bpfile(fname)

    def inside_grid(self,pxy,fmt=0):
        '''
          compute element indices that pxy[npt,2] resides. '-1' means outside of the grid domain
          usage:
               sind=inside_grid(pxy)
               sind,ptr=inside_grid(pxy,fmt=1)
               ptr is triange indice (=1:1st triangle; =2: 2nd triangle), used for computing area coordinates
        '''
        #first triangles
        sindp=self.elnode[:,:3]; px=self.x[sindp]; py=self.y[sindp]
        pie=inside_polygon(pxy,px.T,py.T,fmt=1); fp1=pie!=-1; sind2=nonzero(~fp1)[0]
        if fmt==1: pip=-ones([len(pxy),3]).astype('int'); pip[fp1]=sindp[pie[fp1]]

        #2nd triangles
        if len(sind2)!=0:
           fp34=self.i34==4; sindp=self.elnode[fp34,:][:,array([0,2,3])]; px=self.x[sindp]; py=self.y[sindp]
           pie2=inside_polygon(pxy[sind2],px.T,py.T,fmt=1); fp2=pie2!=-1; pie[sind2[fp2]]=nonzero(fp34)[0][pie2[fp2]]
           if fmt==1: pip[sind2[fp2]]=sindp[pie2[fp2]]

        if fmt==0:
           return pie
        else:
           return [pie,pip]

    def write_shapefile_bnd(self,fname,prjname='epsg:4326'):
        self.shp_bnd=npz_data()
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
        self.shp_bnd.prj=get_prj_file(prjname)
        write_shapefile_data(fname,self.shp_bnd)

    def write_shapefile_node(self,fname,prjname='epsg:4326'):
        self.shp_node=npz_data()
        self.shp_node.type='POINT'
        self.shp_node.xy=c_[self.x,self.y]
        self.shp_node.attname=['id_node']
        self.shp_node.attvalue=arange(self.np)+1;
        self.shp_node.prj=get_prj_file(prjname)
        write_shapefile_data(fname,self.shp_node)

    def write_shapefile_element(self,fname,prjname='epsg:4326'):
        self.shp_elem=npz_data()
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
        self.shp_elem.prj=get_prj_file(prjname)
        write_shapefile_data(fname,self.shp_elem)

class schism_grid_ll(schism_grid):
    def __init__(self):
        pass

    def read_hgrid(self,fname,gr3=None):
        with open(fname,'r') as fid:
            lines=fid.readlines()

        #read ne and np
        num=array(lines[1].split()[0:2]).astype('int')
        self.ne=num[0]; self.np=num[1]

        #read lx,ly and dp
        num=[]
        for i in arange(self.np):
            num.append(array(lines[2+i].split()[1:4]));
        num=array(num).astype('float')
        self.x=num[:,0]
        self.y=num[:,1]
        self.dp=num[:,2]

        #read parents' elnode and bndinfo
        if gr3!=None:
           pattrs=['i34','elnode','nob','nobn','iobn','nlb','nlbn','ilbn','island']
           for pattr in pattrs:
               if hasattr(gr3,pattr):
                  exec('self.'+pattr+'=gr3.'+pattr)

class schism_bpfile:
    def __init__(self):
        self.nsta=0; self.x=array([]); self.y=array([]); self.z=array([]);
        self.station=[]; self.hp=[]; self.ht=[]
        try:
            if mpl._pylab_helpers.Gcf.get_active() is not None:
                abp=gcf().canvas.toolbar.addAction('bp'); abp.triggered.connect(self.edit_bp)
        except:
            pass

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

    def write_bpfile(self,fname,fmt=0):
        '''
        fmt=0: write ACE/gredit *.bp file
        fmt=1: write ACE/gredit *.reg file
        '''

        fid=open(fname,'w+')
        #write header
        if hasattr(self,'note'): fid.write('ACE/gredit: {}'.format(self.note))
        if fmt==0: fid.write('\n{}\n'.format(self.nsta))
        if fmt==1: fid.write('\n1\n{} 1\n'.format(self.nsta))

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

    def write_shapefile(self,fname,prjname='epsg:4326'):
        self.shp_bp=npz_data()
        self.shp_bp.type='POINT'
        self.shp_bp.xy=c_[self.x,self.y]
        self.shp_bp.prj=get_prj_file(prjname)

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

    def edit_bp(self):
        self.cid=gcf().canvas.mpl_connect('button_press_event', self.onclick)
        print('edit bpfile (double click): left=add pt, right=remove pt; middle=stop')

    def onclick(self,sp):
        dlk=int(sp.dblclick); btn=int(sp.button); bx=sp.xdata; by=sp.ydata
        #double click
        if dlk==1 and btn==1: self.add_pt(bx,by)
        if dlk==1 and btn==3: self.remove_pt(bx,by)
        if btn==2: gcf().canvas.mpl_disconnect(self.cid)

    def add_pt(self,x,y):
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

    def remove_pt(self,x,y):
        dist=abs((self.x-x)+1j*(self.y-y))
        sind=nonzero(dist==min(dist))[0][0]; fp=setdiff1d(arange(self.nsta),sind)
        self.nsta=self.nsta-1; self.x=self.x[fp]; self.y=self.y[fp]; self.z=self.z[fp]
        self.station=array(self.station)[fp]
        self.hp[sind][0].remove()
        self.ht[sind].remove()
        del self.hp[sind]; del self.ht[sind]

def read_schism_hgrid(fname):
    #read_schism_hgrid(fname):
    gd=schism_grid()
    gd.read_hgrid(fname)
#    gd.plot_grid()
    return gd

def read_schism_hgrid_ll(fname,gr3=None):
    #read hgrid.ll
    #gr3=read_schism_grid('hgrid.gr3')
    gd=schism_grid_ll()
    gd.read_hgrid(fname,gr3)
    return gd

def read_schism_bpfile(fname,fmt=0):
    '''
    read schism *bp (fmt=0) or *.reg (fmt=1) file created by ACE/gredit
    '''
    bp=schism_bpfile();
    bp.read_bpfile(fname,fmt=fmt)
    return bp

def save_schism_grid(fname='grid',path='.'):
    '''
    read and save path/{hgrid.gr3,vgrid.in}
    '''
    gname='{}/hgrid.gr3'.format(path); vname='{}/vgrid.in'.format(path); S=npz_data();
    if os.path.exists(gname): S.hgrid=read_schism_hgrid(gname)
    if os.path.exists(vname): S.vgrid=read_schism_vgrid(vname)
    if (not hasattr(S,'hgrid')) and (not hasattr(S,'vgrid')): sys.exit('not found: {}, {}'.format(gname,vname))
    save_npz(fname,S)
    return S

class schism_vgrid:
    def __init__(self):
        pass

    def read_vgrid(self,fname):
        #read schism vgrid
        fid=open(fname,'r'); lines=fid.readlines(); fid.close()

        self.ivcor=int(lines[0].strip().split()[0]); self.nvrt=int(lines[1].strip().split()[0])
        if self.ivcor==1:
            #read vgrid info
            lines=lines[2:]
            self.kbp=array([int(i.split()[1])-1 for i in lines]); self.np=len(self.kbp)
            self.sigma=-ones([self.np,self.nvrt])
            for i,line in enumerate(lines):
                self.sigma[i,self.kbp[i]:]=array(line.strip().split()[2:]).astype('float')
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

def read_schism_vgrid(fname):
    '''
    read schism vgrid information working only for ivcor=1
    '''
    vd=schism_vgrid()
    vd.read_vgrid(fname)
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

def create_schism_vgrid(fname='vgrid.in',nvrt=10,fmt=0,h_c=10,theta_b=0.5,theta_f=1.0):
    '''
    create a simple schism pure S vgrid for ivcor=2
    '''
    if fmt==0:
        fid=open(fname,'w+')
        fid.write('2 !ivcor\n{} 1 1.e6 !nvrt\nZ levels\n1 -1.e6\nS levels\n'.format(nvrt))
        fid.write('{} {} {} !h_c,theta_b,theta_f\n'.format(h_c,theta_b,theta_f))
        w=[fid.write('{:5}  {:8.5f}\n'.format(i+1,k)) for i,k in zip(arange(nvrt),linspace(-1,0,nvrt))]
        fid.close()
    else:
        sys.exit('fmt!=0 not working yet')

def getglob(dirpath='.',fmt=0):
    '''
    get global information about schism run (ne,ns,np,nvrt,nproc,ntracers,ntrs)
    dirpath: run directory or outputs directory
    fmt=0: default is 0; fmt(!=0) are for eariler schism versions
    '''
   
    rstr,bdir=srank(0,dirpath=dirpath,fmt=1)
    fname='{}/local_to_global_{}'.format(bdir,rstr) #local_to_global_0000 or local_to_global_000000
    
    #get fname
    #if (fname is None) and os.path.exists('local_to_global_0000'): fname='local_to_global_0000'
    #if (fname is None) and os.path.exists('local_to_global_000000'): fname='local_to_global_000000'
    #if (fname is None) and os.path.exists('./outputs/local_to_global_0000'): fname='./outputs/local_to_global_0000'
    #if (fname is None) and os.path.exists('./outputs/local_to_global_000000'): fname='./outputs/local_to_global_000000'
    #if fname is None: sys.exit('fname unknown')

    #get info
    S=npz_data()
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
    S=npz_data()
    ne=int(lines[0].strip()); np=int(lines[ne+1].strip()); ns=int(lines[ne+np+2].strip())
    S.ielg=array([i.strip().split() for i in lines[1:(ne+1)]])[:,1].astype('int')-1
    S.iplg=array([i.strip().split() for i in lines[(ne+2):(ne+np+2)]])[:,1].astype('int')-1
    S.islg=array([i.strip().split() for i in lines[(ne+np+3):(ne+np+ns+3)]])[:,1].astype('int')-1

    #find line for np,ne
    for i in arange(ne+np+ns+3,len(lines)): 
        sline=lines[i].strip().split()
        if len(sline)!=2: continue
        if int(sline[0])==np and int(sline[1])==ne: nd=i; break; 

    slines=array([i.strip().split() if len(i.split())==5 else [*i.strip().split(),'-1'] for i in lines[(nd+np+1):(nd+np+ne+1)]]).astype('int')
    i34=slines[:,0].astype('int'); elnode=slines[:,1:].astype('int')-1

    S.ne,S.np,S.ns,S.i34,S.elnode=ne,np,ns,i34,elnode
    return S

def read_schism_param(fname,*args):
  with open(fname,'r') as fid:
    lines=fid.readlines()

  param={}
  for line in lines:
    line=line.strip()
    if len(line)==0 or line[0]=='!' or line[0]=='&': continue
    ind=line.find('!');
    if(ind!=-1): line=line[0:ind];
    ind=line.find('=');
    keyi=line[:ind].strip();
    vali=line[(ind+1):].strip();
    param[keyi]=vali
    if((len(args)>0) and (args[0]==1)):
       if vali.lstrip('-').replace('.','',1).isdigit(): param[keyi]=float(vali)
  return param;

def write_schism_param(fname,param):
    pkeys=sorted(param.keys())
    with open(fname,'w+') as fid:
        for i in range(len(pkeys)):
           fid.write('{:10}= {:}\n'.format(pkeys[i],param[pkeys[i]]))

def sms2gr3(fname_2dm,fname_gr3='new.gr3'):
    #2dm to gr3 format: sms2gr3(fname_2dm,fname_gr3='new.gr3')
    gd=schism_grid()

    #read 2dm file
    with open(fname_2dm,'r') as fid:
        lines=fid.readlines()

    # parse every line
    import re
    enum=[]; i34=[]; elnode=[];
    pnum=[]; xyz=[];
    for line in lines:
        #---triangle--
        m=re.match('^E3T (\d+) (\d+) (\d+) (\d+)',line)
        if m!=None:
            enum.append(m.groups()[0])
            elnode.append([*m.groups()[1:],'-1'])
            i34.append(3)
            continue

        #---quads--
        m=re.match('^E4Q (\d+) (\d+) (\d+) (\d+) (\d+)',line)
        if m!=None:
            enum.append(m.groups()[0])
            elnode.append(m.groups()[1:])
            i34.append(4)
            continue

        #----node---
        m=re.match('^ND (\d+) (.*)\n',line)
        if m!=None:
            pnum.append(m.groups()[0])
            xyz.append(m.groups()[1].split())
            continue

    #str2num
    enum=array(enum).astype('int')
    i34=array(i34)
    elnode=array(elnode).astype('int')-1
    ind=argsort(enum)
    enum=enum[ind]; i34=i34[ind]; elnode=elnode[ind,:]

    pnum=array(pnum).astype('int')
    xyz=array(xyz).astype('float64')
    ind=argsort(pnum)
    pnum=pnum[ind]; xyz=xyz[ind,:]

    #assign grid attribute and write grid
    gd.ne=len(enum); gd.np=len(pnum)
    gd.i34=i34; gd.elnode=elnode;
    gd.x=xyz[:,0]; gd.y=xyz[:,1];  gd.dp=xyz[:,2];

    #write grid
    gd.write_hgrid(fname_gr3);
    return gd

if __name__=="__main__":
    pass

