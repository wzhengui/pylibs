#!/usr/bin/env python3
from pylib import *
close('all')

#----------------------------------------------------------------
#inputs
#----------------------------------------------------------------
fnames=['a.asc', 'b.asc' ]  #DEM files (*npz, *.asc, *.tif, *.tiff); For tif/tiff file, center position is assumed
sname='DEM_contour' #name of shapefile to be outputted
levels=[-10,-5, 0]  #contour values (note ocean depth is usually negative)

#optional
nproc=None  #number of threads in parallel computing (default is 10)
reg=None    #region where contours will be extracted (*.bp, *.reg, xy)
#----------------------------------------------------------------
#compute and output contours
#----------------------------------------------------------------
compute_dem_contour(fnames,levels,sname,reg=reg,nproc=nproc)
