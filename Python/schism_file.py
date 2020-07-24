#/usr/bin/env python3
from pylib import *

class schism_grid(object):
    def __init__(self):
        pass

    def plot_grid(self,ax=None,method=0,plotz=0,value=None,mask=None,ec='k',fc='None',
             lw=0.1,levels=None,ticks=None,clim=None,extend='both',cb=True,**args):
        #code for plot grid with default color value (grid depth)
        #method=0: using tricontourf; method=1: using PolyCollection (old method)
        #plotz=0: plot grid only; plotz=1: plot color 
        #value: color value size(np,or ne)
        #mask: size(ne); only plot elements (mask=True))
        #ec: color of grid line;  fc: element color; lw: grid line width
        #levels=100: number of colors for depths; levels=array([v1,v2,...]): depths for plot 
        #ticks=[v1,v2,...]: colorbar ticks; ticks=10: number of ticks
        #clim=[vmin,vmax]: value range for plot/colorbar
        #cb=False: not add colorbar

        if ec==None: ec='None'
        if ax==None: ax=gca();
        if method==0: 
           fp3=self.i34==3; fp4=self.i34==4
           if mask is not None: fp3=fp3*mask; fp4=fp4*mask
           if (plotz==0)|(ec!='None'): #compute lines of grid
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

           if plotz==0:
              hg=plot(r_[x3,x4],r_[y3,y4],lw=lw,color=ec);
           elif plotz==1:
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

              hg=tricontourf(self.x,self.y,tri,value,levels=levels,vmin=vmin,vmax=vmax,extend=extend,**args)

              #add colobar
              if cb:
                 #----new method
                 cm.ScalarMappable.set_clim(hg,vmin=vmin,vmax=vmax)
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

           self.hg=hg
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
           if plotz==0:
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
           self.hg=hg;
           return hg

    def plot_bnd(self,PlotSep=0,c='k',**args):
        if PlotSep==0:
            bi=[]
            for i in arange(self.nob):
                bi=r_[bi,self.iobn[i]]
                bi=r_[bi,self.ilbn[i]]
                bi=bi.astype('int')
                bx=self.x[bi]; by=self.y[bi];
            hb=plt.fill(bx,by,fc='none',ec=c,**args)
            self.hb=hb
        else:
            hob=[];
            for i in arange(self.nob):
                bi=self.iobn[i]
                bi=bi.astype('int')
                bx=self.x[bi]; by=self.y[bi];
                hb=plot(bx,by,color='r',**args)
                hob.append(hb)
            self.hob=hob

            hlb=[];
            for i in arange(self.nob):
                bi=self.ilbn[i]
                bi=bi.astype('int')
                bx=self.x[bi]; by=self.y[bi];
                hb=plot(bx,by,color='g',**args)
                hlb.append(hb)
            self.hlb=hlb
            pass

        hi=[];
        for i in arange(self.nob,self.nlb):
            ibi=self.ilbn[i]
            ibi=ibi.astype('int')
            ibx=self.x[ibi]; iby=self.y[ibi];
            hii=plt.fill(ibx,iby,fc='none',ec=c,**args)
            hi.append(hii)
        self.hi=hi;

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
        self.compute_ns(method=1)

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
        self.iobn=array(self.iobn);

        #read lbnd info
        num=array(lines[n].split()[0]).astype('int'); n=n+2;
        self.nlb=num

        self.nlbn=[];
        self.ilbn=[];
        self.island=[];
        for i in arange(self.nlb):
            if 'island' in lines[n]:
                self.island.append(1)
            else:
                self.island.append(0)
            num=array(lines[n].split()[0]).astype('int')
            self.nlbn.append(num)
            num=[];
            for m in arange(self.nlbn[i]):
                num.append(lines[n+m+1].split()[0])
            self.ilbn.append(array(num).astype('int')-1)

            n=n+self.nlbn[i]+1;
        self.island=array(self.island);
        self.nlbn=array(self.nlbn);
        self.ilbn=array(self.ilbn);

        #old method-----------------
        #        lines0=lines;
        #        lines=remove_tail(lines)
        #        #read obnd info
        #        n=2+self.np+self.ne; num=str2num(lines[n])
        #        self.nob=num.astype('int'); n=n+1
        #        self.nobn=zeros(self.nob,dtype='int')
        #        self.iobn=[]
        #        n=n+1;
        #        for i in arange(self.nob):
        #            num=str2num(lines[n]).astype('int')
        #            self.nobn[i]=num[0]
        #            n=n+1; num=str2num(lines[n:(n+self.nobn[i])]).astype('int')-1
        #            self.iobn.append(squeeze(num))
        #            n=n+self.nobn[i]
        #        self.iobn=array(self.iobn)
        #
        #        #read lbnd info
        #        num=str2num(lines[n])
        #        self.nlb=num.astype('int'); n=n+1
        #        self.nlbn=zeros(self.nlb,dtype='int')
        #        self.island=zeros(self.nlb,dtype='int')
        #        self.ilbn=[]
        #        n=n+1;
        #        for i in arange(self.nlb):
        #            num=str2num(lines[n]).astype('int')
        #            if 'island' in lines0[n]:
        #                self.island[i]=1
        #            self.nlbn[i]=num[0]
        #            n=n+1; num=str2num(lines[n:(n+self.nlbn[i])]).astype('int')-1
        #            self.ilbn.append(squeeze(num))
        #            n=n+self.nlbn[i]
        #        self.ilbn=array(self.ilbn)

    def interp_node_to_elem(self,value=None):
        #interpolate node values to elements

        #--get node value
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
        #interpolate element values to nodes
        #if value not given, dpe is used
        #method=0: simple avarage; method=1: inverse distance (power=p)

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

    def compute_ctr(self,value=None):
        self.xctr=zeros(self.ne)*nan
        self.yctr=zeros(self.ne)*nan
        #self.dpe=zeros(self.ne)*nan

        fp1=self.i34==3; fp2=self.i34==4;
        self.xctr[fp1]=mean(self.x[self.elnode[fp1,0:3]],axis=1)
        self.yctr[fp1]=mean(self.y[self.elnode[fp1,0:3]],axis=1)
        self.xctr[fp2]=mean(self.x[self.elnode[fp2,:]],axis=1)
        self.yctr[fp2]=mean(self.y[self.elnode[fp2,:]],axis=1)

        if value is None:
            self.dpe=self.interp_node_to_elem()
        else:
            self.dpe=self.interp_node_to_elem(value=value)

        #self.dpe[fp1]=mean(self.dp[self.elnode[fp1,0:3]],axis=1)
        #self.dpe[fp2]=mean(self.dp[self.elnode[fp2,:]],axis=1)
        #self.dpe=array([self.dp[e[:-1]].mean() if (i34==3 and len(e)==4) else self.dp[e].mean() \
        #               for e,i34 in zip(self.elnode,self.i34)]);

    def compute_area(self):
        fp=self.elnode[:,-1]<0;
        x1=self.x[self.elnode[:,0]]; y1=self.y[self.elnode[:,0]];
        x2=self.x[self.elnode[:,1]]; y2=self.y[self.elnode[:,1]];
        x3=self.x[self.elnode[:,2]]; y3=self.y[self.elnode[:,2]];
        x4=self.x[self.elnode[:,3]]; y4=self.y[self.elnode[:,3]]; x4[fp]=x1[fp]; y4[fp]=y1[fp]
        self.area=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)+(x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2

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

    def compute_ns(self,method=1):
        #compute number of grid sides

        if method==1:
            #triangles
            ind=nonzero(self.i34==3)[0]
            if len(ind)>0:
                for i in arange(3):
                    i1=mod(i+3,3); i2=mod(i+4,3);
                    if i==0:
                        x=c_[self.elnode[ind,i1],self.elnode[ind,i2]]
                    else:
                        x=r_[x,c_[self.elnode[ind,i1],self.elnode[ind,i2]]];

            #quadlateral
            ind=nonzero(self.i34==4)[0]
            if len(ind)>0:
                for i in arange(4):
                    i1=mod(i+4,4); i2=mod(i+5,4);
                    x=r_[x,c_[self.elnode[ind,i1],self.elnode[ind,i2]]];

            #-sort sides-----------
            y=sort(x,axis=1)
            uy=unique(y,axis=0)

            self.ns=uy.shape[0]
        elif method==2:
            #--method 2--------------------
            #open bondary
            nobs=0
            for i in arange(self.nob):
                nobs=nobs+self.nobn[i]-1

            nlbs=0
            #land boundary
            for i in arange(self.nlb):
                if self.ilbn[i][0]==self.ilbn[i][-1]:
                    nlbs=nlbs+self.nlbn[i]-2+self.island[i]
                else:
                    nlbs=nlbs+self.nlbn[i]-1+self.island[i]

            #compute ns
            fp1=self.i34==3; fp2=self.i34==4;
            ns=int((sum(fp1)*3+sum(fp2)*4+nobs+nlbs)/2);
            self.ns=ns

    def compute_node_ball(self):
        nne=zeros(self.np).astype('int');
        ine=[[] for i in arange(self.np)];
        for i in arange(self.ne):
            inds=self.elnode[i,:self.i34[i]];
            nne[inds]=nne[inds]+1
            [ine[indi].append(i) for indi in inds]
        self.nne=nne
        self.ine=array([array(ine[i]) for i in arange(self.np)]);
   
    def compute_acor(self,pxy,N=100):
        #compute acor coodinate for points pts(xi,yi)
        
        #compute the corresponding residing elem 
        npt=pxy.shape[0]; pind=arange(npt); ie=array([]).astype('int');
        while(len(pind)>0):
             #get the first N pts)
             if len(pind)<=N:
                pindi=pind
             else:
                pindi=pind[:N] 
             pind=setdiff1d(pind,pindi)
             ie=r_[ie,self.inside_grid(pxy[pindi])]
         
        ie=ie.astype('int') 
        #compute area coordinate
        ip0=self.elnode[ie]; i34=self.i34[ie]
        x1=self.x[ip0][:,0]; x2=self.x[ip0][:,1]; x3=self.x[ip0][:,2]; x4=self.x[ip0][:,3]; x=pxy[:,0]
        y1=self.y[ip0][:,0]; y2=self.y[ip0][:,1]; y3=self.y[ip0][:,2]; y4=self.y[ip0][:,3]; y=pxy[:,1]

        A1=((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))/2
        A11=((x2-x)*(y3-y)-(x3-x)*(y2-y))/2
        A12=((x-x1)*(y3-y1)-(x3-x1)*(y-y1))/2
        A13=((x2-x1)*(y-y1)-(x-x1)*(y2-y1))/2
        fA1=(abs(A11)+abs(A12)+abs(A13)-abs(A1))/abs(A1);

        A2=((x3-x1)*(y4-y1)-(x4-x1)*(y3-y1))/2
        A21=((x3-x)*(y4-y)-(x4-x)*(y3-y))/2
        A22=((x-x1)*(y4-y1)-(x4-x1)*(y-y1))/2
        A23=((x3-x1)*(y-y1)-(x-x1)*(y3-y1))/2
        fA2=(abs(A21)+abs(A22)+abs(A23)-abs(A2))/abs(A2);

        ip=[];acor=[];
        for i in arange(npt):
            if abs(fA1[i])<1e-5: #pt in 1st triange
               a1=A11[i]/A1[i]; a2=A12[i]/A1[i]; a3=1-a1-a2;
               if (a1<0)|(a2<0): sys.exit('check pt: {}, {}, {}, {}, {},{}'.format(i,i34[i],ie[i],A1[i],A11[i],A12[i]))
               ip.append(ip0[i,:3]); acor.append(array([a1,a2,a3]))
            else:
               a1=A21[i]/A2[i]; a2=A22[i]/A2[i]; a3=1-a1-a2;
               if (a1<0)|(a2<0)|(i34[i]==3)|(abs(fA2[i])>=1e-5): sys.exit('check pt: {}, {}, {}, {}, {},{}'.format(i,i34[i],ie[i],A1[i],A11[i],A12[i]))
               ip.append(ip0[i,array([0,2,3])]); acor.append(array([a1,a2,a3]))
        ip=array(ip); acor=array(acor)

        #save acor
        return ie,ip,acor             

    def interp(self,xyi,*args):
        return interpolate.griddata(c_[self.x,self.y],self.dp,xyi,*args);

    def write_hgrid(self,fname,elnode=1,bnd=0,Info=None):
        with open(fname,'w+') as fid:
            fid.write('!grd info:{}\n'.format(Info))
            fid.write('{} {}\n'.format(self.ne,self.np))
            for i in arange(self.np):
                #fid.write('{:<6d} {:<16.8f} {:<16.8f} {:<10.3f}\n'.format(i+1,self.x[i],self.y[i],self.dp[i]))
                fid.write('{:<d} {:<.8f} {:<.8f} {:<.8f}\n'.format(i+1,self.x[i],self.y[i],self.dp[i]))
            if elnode!=0:
                for i in arange(self.ne):
                    #if self.i34[i]==3: fid.write('{:<6d} {:2d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))
                    #if self.i34[i]==4: fid.write('{:<6d} {:2d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))
                    if self.i34[i]==3: fid.write('{:<d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))
                    if self.i34[i]==4: fid.write('{:<d} {:d} {:d} {:d} {:d} {:d}\n'.format(i+1,self.i34[i],*self.elnode[i,:]+1))

    def split_quads(self,angle_min=60,angle_max=120,fname='new.gr3'):
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


    def check_quads(self,angle_min=60,angle_max=120,fname='bad_quad.xyz'):
        #check the quality of quads, violation when internal angle < angle_min, or >angle_max

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

        #output bad_quad location as xyz file
        if not hasattr(self,'xctr'): self.compute_ctr()
        qxi=self.xctr[self.index_bad_quad]; qyi=self.yctr[self.index_bad_quad];
        with open("{}".format(fname),'w+') as fid:
            for i in arange(len(qxi)):
                fid.write('{} {} 0\n'.format(qxi[i],qyi[i]))

    def plot_bad_quads(self,color='r',ms=12,*args):
        #plot grid with bad quads
        if not hasattr(self,'index_bad_quad'): self.check_quads()
        if not hasattr(self,'xctr'): self.compute_ctr()

        qxi=self.xctr[self.index_bad_quad]; qyi=self.yctr[self.index_bad_quad]
        self.plot_grid()
        plot(qxi,qyi,'.',color=color,ms=ms,*args)
        pass;

    def inside_grid(self,pxy):
        x0=pxy[:,0][:,None]; y0=pxy[:,1][:,None]

        # for all trignales and first triganges of quads--
        x1=self.x[self.elnode[:,0]][None,:]; x2=self.x[self.elnode[:,1]][None,:]; x3=self.x[self.elnode[:,2]][None,:]
        y1=self.y[self.elnode[:,0]][None,:]; y2=self.y[self.elnode[:,1]][None,:]; y3=self.y[self.elnode[:,2]][None,:]

        a1=(x0-x3)*(y2-y3)-(x2-x3)*(y0-y3)
        a2=(x1-x3)*(y0-y3)-(x0-x3)*(y1-y3)
        a3=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)

        fp=(a1>=0)*(a2>=0)*(a3>=0);
        pn=[];
        for i in arange(pxy.shape[0]):
            ind=nonzero(fp[i,:])[0]
            if len(ind)==0:
                pn.append(None)
            else:
                pn.append(ind[0])
        pn=array(pn)

        # for the 2nd trignale of quads--
        fp=self.i34==4; ind0=nonzero(fp)[0];
        fpn=pn==None; indn=nonzero(fpn)[0];  x0=x0[indn]; y0=y0[indn] #remaining pts
        if len(ind0)!=0 and len(indn)!=0 :
            x1=self.x[self.elnode[fp,0]][None,:]; x2=self.x[self.elnode[fp,2]][None,:]; x3=self.x[self.elnode[fp,3]][None,:]
            y1=self.y[self.elnode[fp,0]][None,:]; y2=self.y[self.elnode[fp,2]][None,:]; y3=self.y[self.elnode[fp,3]][None,:]

            a1=(x0-x3)*(y2-y3)-(x2-x3)*(y0-y3)
            a2=(x1-x3)*(y0-y3)-(x0-x3)*(y1-y3)
            a3=(x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)

            fp=(a1>=0)*(a2>=0)*(a3>=0);
            pn2=[];
            for i in arange(len(indn)):
                ind=nonzero(fp[i,:])[0]
                if len(ind)==0:
                    pn2.append(None)
                else:
                    pn2.append(ind[0])
            pn2=array(pn2)
            pn[fpn]=ind0[pn2]

        return pn.astype('int')

    def write_shapefile_bnd(self,fname,prjname='epsg:4326'):
        self.shp_bnd=npz_data()
        self.shp_bnd.type='POLYLINE'
        for i in arange(self.nob):
            ind=self.iobn[i]
            xyi=c_[self.x[ind],self.y[ind]];
            xyi=insert(xyi,0,nan,axis=0);
            if i==0:
                xy=xyi
            else:
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
        #gr3=read_schism_grid('hgrid.gr3')
        with open(fname,'r') as fid:
            lines0=fid.readlines()
        lines=remove_tail(lines0)
        #read ne and np
        num=str2num(lines[1]).astype('int')
        self.ne=num[0]; self.np=num[1]

        #read lx,ly and dp
        n=2; num=str2num(lines[n:(n+self.np)]);
        self.x=num[:,1]
        self.y=num[:,2]
        self.dp=num[:,3]

        #read parents' elnode and bndinfo
        if gr3!=None:
           pattrs=['i34','elnode','nob','nobn','iobn','nlb','nlbn','ilbn','island']
           for pattr in pattrs:
               if hasattr(gr3,pattr):
                  exec('self.'+pattr+'=gr3.'+pattr)

class schism_bpfile(object):
    def __init__(self):
        pass

    def read_bpfile(self,fname):
        lines=[line.split() for line in open(fname,'r').read().strip().split('\n')]
        self.nsta=int(lines[1][0])

        fc=lambda x: x if len(x)==4 else [*x[:4],x[4][1:]]
        data=array([fc(line) for line in lines[2:(2+self.nsta)]])

        self.x=data[:,1].astype(float)
        self.y=data[:,2].astype(float)
        self.z=data[:,3].astype(float)
       
        upxy,ind=unique(self.x+1j*self.y,return_index=True)
        ux=self.x[ind]; uy=self.y[ind];uz=self.z[ind]
        self.ux=ux; self.uy=uy;self.uz=uz;
        #get unique station data.
        if data.shape[1]==5: 
           self.station=data[:,4]
           self.ustation=self.station[ind]; 
        else:
           self.station=array(['{}'.format(i) for i in arange(self.nsta)])
           self.ustation=self.station[ind]

    def write_bpfile(self,fname):
        with open(fname,'w+') as fid:
            fid.write('!\n{}\n'.format(self.nsta))
            for i in arange(self.nsta):
                if hasattr(self,'station'):
                    fid.write('{:<d} {:<.8f} {:<.8f} {:<.8f} !{}\n'.format(i+1,self.x[i],self.y[i],self.z[i],self.station[i]))
                else:
                    fid.write('{:<d} {:<.8f} {:<.8f} {:<.8f}\n'.format(i+1,self.x[i],self.y[i],self.z[i]))

    def plot_station(self,ax=None,ls='',label=True,**args):
        if not None: ax=gca()
        hp=plot(self.ux,self.uy,linestyle=ls,**args)
        self.hp=hp
        if label:
           ht=[];
           for i in arange(len(self.ustation)):
               hti=text(self.ux[i],self.uy[i],self.ustation[i],color='r')
               ht.append(hti)
           self.ht=array(ht)

    def compute_acor(self,gd): 
        #compute areal coordinates, and gd is the schism grid
        self.ie,self.ip,self.acor=gd.compute_acor(c_[self.x,self.y])

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

def read_schism_bpfile(fname):
    bp=schism_bpfile();
    bp.read_bpfile(fname)
    return bp

def read_schism_vgrid(fname,gd,node=None,eta=0,flag=0):
    #read vgrid
    with open(fname,'r') as fid:
        lines=fid.readlines()
    R=re.findall(r'\d+',lines[0]); ivcor=int(R[0])
    R=re.findall(r'\d+',lines[1]); nvrt=int(R[0])

    sigma=ones([gd.np,nvrt])*nan;
    if ivcor==1:
        for i in arange(gd.np):
            R=loadtxt(StringIO(lines[i+2]));
            kbp=int(R[1]-1);
            sigma[i,kbp:]=R[2:]
        dpi=tile(gd.dp,[nvrt,1]).T
        zcor=sigma*dpi
    else:
        zcor='not work for ivcor/=1 yet'

    # if only subset of nodes is needed
    if node is not None:
        zcor=zcor[node,:]

    #extend values in the bottom
    if flag!=0:
        mzcor=nanmin(zcor,axis=1)
        sz=zcor.shape
        for i in arange(zcor.shape[1]):
            fp=isnan(zcor[:,i]);
            zcor[fp,i]=mzcor[fp]

    return zcor

def getglob(fname,flag=0):
    glob=npz_data()
    with open(fname,'r') as fid:
      line=fid.readline().split()
      if flag==0: #old format (ne,np,nvrt,nproc,ntrs)
        glob.ne=int(line[0])
        glob.np=int(line[1])
        glob.nvrt=int(line[2])
        glob.nproc=int(line[3])
        glob.ntrs=array(line[4:]).astype('int')
        glob.ntracers=sum(glob.ntrs)
      elif flag==1: #new format (ns,ne,np,nvrt,nproc,ntrs)
        glob.ns=int(line[0])
        glob.ne=int(line[1])
        glob.np=int(line[2])
        glob.nvrt=int(line[3])
        glob.nproc=int(line[4])
        glob.ntrs=array(line[5:]).astype('int')
        glob.ntracers=sum(glob.ntrs)
      elif flag==2: #new format (ns,ne,np,nvrt,nproc,ntracers,ntrs)
        glob.ns=int(line[0])
        glob.ne=int(line[1])
        glob.np=int(line[2])
        glob.nvrt=int(line[3])
        glob.nproc=int(line[4])
        glob.ntracers=int(line[5])
        glob.ntrs=array(line[6:]).astype('int')
    return glob

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

##----read vertical grid--------------------------------------------------------
#    hgrid=r'D:\Work\E3SM\ChesBay\Hycom\hgrid.gr3'
#    vgrid=r'D:\Work\E3SM\ChesBay\Hycom\vgrid.in'
#
#    gd=read_schism_hgrid(hgrid)
#    zcor=read_schism_vgrid(vgrid,gd,node=gd.iobn[0],flag=1)

#----check quads and plot------------------------------------------------------
#    fname=r'D:\Work\E3SM\ChesBay\Grid\ChesBay_1.gr3'
#    fname1=r'D:\Work\E3SM\ChesBay\Grid\ChesBay_1a.gr3'
#    gd=read_schism_hgrid(fname);
#    gd.check_quads(70,110,r'D:\Work\E3SM\ChesBay\Grid\AA.xyz');
#    gd.split_quads(70,110,fname1)


#    gd=read_schism_hgrid(fname);
#    gd.check_quads(angle_min=70,angle_max=110);
#    gd.plot_bad_quads();


#---transform 2dm to gr3 format, and do projection of grid---------------------
#    fname_2dm=r'D:\Work\E3SM\ChesBay\Grid\ChesBay_1.2dm'
#    fname=r'D:\Work\E3SM\ChesBay\Grid\A.gr3'
#    fname_ll=r'D:\Work\E3SM\ChesBay\Grid\A.ll'
#
#    sms2gr3(fname_2dm,fname);
#    proj(fname,0,'epsg:26918',fname_ll,0,'epsg:4326');
#------------------------------------------------------------------------------
#    g=read_schism_hgrid(fname)

#    fname=r'C:\Users\Zhengui\Desktop\Python\learn\Station.bp_SED';
#    bp=read_schism_bpfile(fname)
#    bp.ux,bp.uy=transform(Proj(init='epsg:26910'),Proj(init='epsg:4326'),bp.ux,bp.uy);

#------------------read grid and plot------------------------------------------
#    fname=r'C:\Users\Zhengui\Desktop\Python\learn\hgrid.gr3'
#    g=read_schism_hgrid(fname)
#    fname=r'C:\Users\Zhengui\Desktop\Python\learn\T.gr3'
#    g.write_hgrid(fname,1)

#    fname=r'C:\Users\Zhengui\Desktop\Python\learn\hgrid.ll'
#    g2=read_schism_hgrid_ll(fname,g)
#
#    g2.plot_bnd(); bp.plot_station(c='g',marker='.',ms=10);
#    [setp(hti,color='r',fontsize=10,fontweight='bold') for hti in bp.ht]

#    fname=r'C:\Users\Zhengui\Desktop\Python\learn\bk.gr3'
#    g=read_schism_hgrid(fname)
#    dpe=[g.dp[e[:-1]].mean() if (i34==3 and len(e)==4) else g.dp[e].mean() for e,i34 in zip(g.elnode,g.i34)];
#    g.plot_grid(plotz=1,cmap='hsv',ec='None',clim=[0,25])
#    g.plot_bnd(lw=1)


#    fname2=r'C:\Users\Zhengui\Desktop\Python\learn\hgrid.ll'
#    g2=read_schism_hgrid_ll(fname2,g)
