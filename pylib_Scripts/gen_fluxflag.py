#!/usr/bin/env python3
#write fluxflag.prop in regions (*.reg) for computing fluxes 
#note: *.reg is restrcited to 4 sorted pts with the 1st two pts 
#      define upstream, and the 2nd two pts defines downstream
from pylib import *

#--------------------------------------------------------
#inputs
#--------------------------------------------------------
grd='../grid.npz'      #model grid (*.npz, or *.gr3)
regions=['DSJ.reg','SJR.reg','TSL.reg','SCR.reg']    #regions

#--------------------------------------------------------
#read grid
#--------------------------------------------------------
gd=loadz(grd).hgrid if grd.endswith('.npz') else read_schism_hgrid(grd)
#--------------------------------------------------------
#assign different values in regions
#--------------------------------------------------------
pvi=-ones(gd.ne).astype('int'); gd.compute_ctr() 
for m,region in enumerate(regions):
    #read region info
    bp=read_schism_bpfile(region,fmt=1)
    if bp.nsta!=4: sys.exit(f'{region}''s npt!=4')
    x1,x2,x3,x4=bp.x; y1,y2,y3,y4=bp.y 

    #middle pts
    mx1=(x1+x4)/2; mx2=(x2+x3)/2
    my1=(y1+y4)/2; my2=(y2+y3)/2

    #for lower region
    px=array([mx1,mx2,x3,x4]); py=array([my1,my2,y3,y4])
    pvi[inside_polygon(c_[gd.xctr,gd.yctr],px,py)==1]=m
  
    #for upper region
    px=array([mx1,x1,x2,mx2]); py=array([my1,y1,y2,my2])
    pvi[inside_polygon(c_[gd.xctr,gd.yctr],px,py)==1]=m+1
gd.write_prop('fluxflag.prop',value=pvi,fmt='{:3d}')
