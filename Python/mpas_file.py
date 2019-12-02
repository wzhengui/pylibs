#!/usr/bin/env python3
from pylib import *



class mpas_grid(object):
    def __init__(self):
        pass

    def plot_grid(self,ax=None,ec='k',lw=0.1,fc='None',plotz=0,value=None,**args):
        if ax==None: ax=gca();
        x6=self.lonVertex[self.verticesOnCell]*180/pi; y6=self.latVertex[self.verticesOnCell]*180/pi

        #make sure grid lons close
        pind=nonzero((x6.max(axis=1)-x6.min(axis=1))>180)[0]
        for pindi in pind:
            if abs(mean(y6[pindi]))>88: continue
            fp=x6[pindi,:]>90; x6[pindi,fp]=x6[pindi,fp]-360

        xy6=c_[x6[:,:,None],y6[:,:,None]]
        if plotz==0:
            hg=mpl.collections.PolyCollection(xy6,lw=lw,edgecolor=ec,facecolor=fc,antialiased=False,**args)
        else:
            hg=mpl.collections.PolyCollection(xy6,lw=lw,edgecolor=ec,array=value,antialiased=False,**args)
            hc=colorbar(hg);
            self.hc=hc;

        ax.add_collection(hg)
        ax.autoscale_view()
        self.hg=hg

        return hg

    def read_grid(self,fname,*args):
        #read mpas grid'
        C=ReadNC(fname); cvars=[cvar for cvar in C.variables]
        svars=array(('latCell','lonCell','xCell','yCell','zCell','latEdge','lonEdge',
                     'xEdge','yEdge','zEdge','latVertex','lonVertex','xVertex','yVertex',
                     'zVertex','cellsOnEdge','nEdgesOnCell','nEdgesOnEdge','edgesOnCell',
                     'edgesOnEdge','dvEdge','dcEdge','angleEdge','areaCell','cellsOnCell',
                     'verticesOnCell','verticesOnEdge','edgesOnVertex','cellsOnVertex'))
        fvars=array(('cellsOnEdge','cellsOnCell','cellsOnVertex','edgesOnCell','edgesOnEdge',
                     'edgesOnVertex','verticesOnCell','verticesOnEdge'))
        for svar in svars:
            if svar not in cvars: continue
            if svar in fvars:
                exec("self.{}=array(C.variables['{}'])-1".format(svar,svar))
            else:
                exec("self.{}=array(C.variables['{}'])".format(svar,svar))

        #add new parameters, and process grid
        self.nCell=len(self.xCell); self.nEdge=len(self.xEdge); self.nVertex=len(self.xVertex)
        fp=(self.verticesOnCell==self.nVertex)|(self.verticesOnCell==-1); cind=nonzero(fp)
        for i,j in zip(*cind):
            self.verticesOnCell[i,j]=self.verticesOnCell[i,j-1]

def read_mpas_grid(fname):
    gd=mpas_grid()
    gd.read_grid(fname)
    return gd


if __name__=="__main__":
    pass;
    #gd=read_mpas_grid('init.nc')
    #gd.plot_grid(lw=0.2,fc='r');

