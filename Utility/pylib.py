#---------------------------------------------------------------------
#import system lib
#---------------------------------------------------------------------
import os,sys

Libs=['pylib','mylib','schism_file']
if not set(Libs).issubset(set(sys.modules.keys())):
   #---------------------------------------------------------------------
   #spyder pylab library
   #C:/Program Files/Python/Python3[57]/Lib/site-packages/matplotlib/pylab.py
   #---------------------------------------------------------------------
   pversion=sys.version.split(' ')[0]
   #print(pversion)

   #---------------------------------------------------------------------
   #libraries of packages
   #---------------------------------------------------------------------
   #matplotlib
   import matplotlib as mpl
   if os.getenv('HOSTNAME') is not None:
      if 'frontera' in os.getenv('HOSTNAME'): mpl.use('tkagg')
   from matplotlib import pyplot as plt
   from matplotlib import cbook, mlab
   from matplotlib.dates import *
   from matplotlib.pyplot import *

   #numpy
   import numpy as np
   from numpy import *
   from numpy.random import *
   from numpy.linalg import *
   import numpy.ma as ma
   #from numpy.fft import *

   #scipy
   import scipy as sp
   from scipy import (optimize,interpolate,signal)
   from scipy.fftpack import fft, ifft
   #from scipy import (optimize,interpolate,io,signal)

   #pandas
   import pandas as pd

   #misc
   import re
   import datetime
   #from io import StringIO
   #import imp
   #import importlib as imp

   #proj
   from pyproj import Transformer
   #from pyproj import Proj, transform

   #netcdf
   from netCDF4 import Dataset

   #excel
   try:
      import xlsxwriter as xw
   except:
      pass

   #mpi4py
   try:
      from mpi4py import MPI
   except:
       pass

   #url download
   try:
      import urllib
   except:
      pass

   #reload
   try:
     from importlib import reload
   except:
     pass

   #---------------------------------------------------------------------
   #libraries of self-defined modules
   #---------------------------------------------------------------------
   import mylib
   from mylib import (get_xtick,close_data_loop,datenum,
        loadz,npz_data,save_npz,find_continuous_sections,
        smooth,daytime_length,move_figure,lpfilt,mdivide,signa,
        inside_polygon,command_outputs,near_pts,proj,
        get_prj_file,mfft,read_shapefile_data,write_shapefile_data,
        ReadNC,WriteNC,harmonic_fit,harmonic_analysis,get_hycom,
        get_stat,get_subplot_position,load_bathymetry,plot_taylor_diagram,
        convert_dem_format)
        #convert_matfile_format,

   import schism_file
   from schism_file import (read_schism_hgrid, read_schism_hgrid_ll,read_schism_bpfile,getglob,
        schism_grid,schism_vgrid,schism_bpfile,sms2gr3,read_schism_vgrid,save_schism_grid,
        compute_zcor,read_schism_param,write_schism_param,read_schism_local_to_global,
        create_schism_vgrid)

   if os.getenv('HOME')!=None:
       sys.path.append(os.getenv('HOME'))

   #sys.modules['loadz'] = mylib #in case oldmodule name used
   #sys.modules['read_schism_file'] = schism_file #in case oldmodule name used

   #import mpas_file
   #from mpas_file import (read_mpas_grid)


