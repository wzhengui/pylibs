#---------------------------------------------------------------------
#import system lib
#---------------------------------------------------------------------
import os,sys
from glob import glob

Libs=['pylib','mylib','schism_file']
if not set(Libs).issubset(set(sys.modules.keys())):
   pversion=sys.version.split(' ')[0] #print(pversion)

   #---------------------------------------------------------------------
   #load pylib libraries of packages
   #---------------------------------------------------------------------
   #matplotlib
   import matplotlib as mpl
   HNAME=str(os.getenv('HOSTNAME')); TNAME=str(os.getenv('TACC_SYSTEM'))
   if ('frontera' in HNAME) or ('stampede2' in HNAME): mpl.use('tkagg')
   if ('frontera' in TNAME) or ('stampede2' in TNAME): mpl.use('tkagg')
   from matplotlib import pyplot as plt
   from matplotlib.dates import date2num, datestr2num,num2date
   if hasattr(mpl.dates,'set_epoch'):
      try:
         mpl.dates.set_epoch('0000-12-31')
      except:
         pass
   from matplotlib.pyplot import *

   import platform
   if platform.system().lower()=='windows':
      try:
         if get_ipython().__class__.__name__!='ZMQInteractiveShell': mpl.use('Qt5Agg')
      except:
         pass

   #numpy
   import numpy as np
   from numpy import *
   from numpy.random import rand,randn
   from numpy.linalg import *

   #scipy
   import scipy as sp
   from scipy import interpolate
   from scipy.fftpack import fft, ifft

   #mpi4py
   try:
      from mpi4py import MPI
   except:
      pass
      #from pyUtility.mylib import parallel_jobs
      #MPI=parallel_jobs()

   #url download
   try:
      import urllib
      from urllib.request import urlretrieve as urlsave
      import ssl
      try:
          _create_unverified_https_context = ssl._create_unverified_context
      except AttributeError:
          # Legacy Python that doesn't verify HTTPS certificates by default
          pass
      else:
          # Handle target environment that doesn't support HTTPS verification
          ssl._create_default_https_context = _create_unverified_https_context
   except:
      pass

   #------------------------------------------------
   #old import
   #------------------------------------------------
   #from numpy.random import *
   #import numpy.ma as ma
   #from matplotlib import cbook, mlab
   #from matplotlib.dates import *
   #netcdf
   #from netCDF4 import Dataset
   #misc
   #import datetime
   #import pandas as pd
   #sympy
   #try:
   #  from sympy import init_session as sym_init
   #except:
   #  pass
   #pickle
   #import pickle
   #import copy
   #from copy import copy as scopy
   #from copy import deepcopy as dcopy

   #---------------------------------------------------------------------
   #libraries of self-defined modules
   #---------------------------------------------------------------------
   import pyUtility.mylib as mylib
   sys.modules['mylib'] = mylib
   from pyUtility.mylib import (get_xtick,close_data_loop,datenum,quickdatenum,
        get_INFO,loadz,zdata,savez,find_cs,npz2mat,convert_matfile,
        smooth,daytime_length,move_figure,bpfilt,lpfilt,mdivide,signa,
        inside_polygon,mdist,command_outputs,near_pts,proj,proj_pts,rewrite,rewrite_input,
        get_prj_file,mfft,interp_vertical,read_shapefile_data,write_shapefile_data,
        ReadNC,WriteNC,harmonic_fit,harmonic_analysis,get_hycom,compute_contour,EOF,
        get_stat,get_subplot_position,get_subplot_position2,load_bathymetry,plot_taylor_diagram,
        convert_dem_format,get_hpc_command,least_square_fit,read_yaml,read_excel,
        write_excel,rtext,mklink,cindex,resize,savefig,pplot,blit_manager,read,add_xtick,
        get_qnode,modify_figure,parallel_jobs,fig_IFNO,ceqstate,subdomain_index)

   import pyUtility.schism_file as schism_file
   sys.modules['schism_file'] = schism_file
   from pyUtility.schism_file import (read_schism_hgrid, read_schism_bpfile,getglob,
        schism_grid,schism_vgrid,schism_bpfile,sms2grd,read_schism_vgrid,save_schism_grid,
        compute_zcor,read_schism_param,write_schism_param,read_schism_local_to_global,
        create_schism_vgrid,srank,grd2sms,scatter_to_schism_grid,delete_schism_grid_element,
        read_schism_prop,read_schism_reg,interp_schism_3d,get_schism_var_info,
        read_schism_output,change_schism_param,get_schism_output_info,get_schism_grid_subdomain,
        get_schism_output_subset,combine_schism_hotstart,combine_icm_output,read_schism_slab,
        convert_schism_source,schism_view,schism_check,zcor_to_schism_grid,compute_schism_volume)

   if os.getenv('HOME')!=None:
       sys.path.append(os.getenv('HOME'))

   #sys.modules['loadz'] = mylib #in case oldmodule name used
   #sys.modules['read_schism_file'] = schism_file #in case oldmodule name used

   #import mpas_file
   #from mpas_file import (read_mpas_grid)

   #---------------------------------------------------------------------
   #alias
   #---------------------------------------------------------------------
   from os.path import exists as fexist
   from pyUtility.mylib import savez as save_npz; mylib.save_npz=savez
   from pyUtility.mylib import zdata as npz_data; mylib.npz_data=zdata
   from pyUtility.mylib import least_square_fit as lsq; mylib.least_square_fit=lsq
   from pyUtility.mylib import move_figure as mvfig
   from pyUtility.mylib import modify_figure as mfig
   from pyUtility.mylib import find_cs as find_continuous_sections; mylib.find_continuous_sections=find_cs
   from pyUtility.mylib import convert_matfile as mat2npz
   from pyUtility.mylib import convert_matfile as loadm
   from pyUtility.mylib import convert_dem_format as read_dem
   from pyUtility.mylib import read_shapefile_data as read_shp
   from pyUtility.mylib import write_shapefile_data as write_shp
   from pyUtility.schism_file import read_schism_hgrid as read_grd
   from pyUtility.schism_file import read_schism_bpfile as read_bp
   from pyUtility.schism_file import change_schism_param as chparam
