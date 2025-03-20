#Copyright 2021, Zhengui WANG
#Apache License, Version 2.0; http://www.apache.org/licenses/LICENSE-2.0
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
   #temp fix: try to be compatible with higher numpy version
   if ('numpy.core' in sys.modules) and ('numpy._core' not in sys.modules):
      sys.modules['numpy._core']=sys.modules['numpy.core']
      nms=[i.split('.')[-1] for i in sys.modules if i.startswith('numpy.core.')]
      for i in nms: sys.modules['numpy._core.'+i]=sys.modules['numpy.core.'+i]

   #scipy
   import scipy as sp
   from scipy import interpolate
   from scipy.fftpack import fft, ifft

   #mpi4py
   #try:
   #   from mpi4py import MPI
   #except:
   #   pass
   #   #from src.mylib import parallel_jobs
   #   #MPI=parallel_jobs()

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
   path_pylib=os.path.dirname(__file__)
   if os.path.exists(os.path.dirname(__file__)+'/pylibs/src'):
      import pylibs.src.mylib as mylib; sys.modules['src']=sys.modules['pylibs.src']
      path_src=path_pylib+'/pylibs/src'; path_scripts=path_pylib+'/pylibs/scripts'
   else:
      import src.mylib as mylib; path_src=path_pylib+'/src'; path_scripts=path_pylib+'/scripts'
   sys.modules['mylib'] = mylib
   from mylib import (xtick,get_xtick,close_data_loop,datenum,quickdatenum,
        add_basemap,get_INFO,loadz,zdata,savez,find_cs,npz2mat,read_mat,sort_all,
        cmean,smooth,doy,daytime_length,move_figure,bpfilt,lpfilt,mdivide,signa,sub_lines,sub_polygons,
        inside,inside_polygon,mdist,command_outputs,near_pts,proj,proj_pts,rewrite,rewrite_input,
        get_prj_file,mfft,interp_vertical,read_shapefile_data,write_shapefile_data,
        ReadNC,WriteNC,harmonic_fit,harmonic_analysis,get_hycom,compute_contour,EOF,
        get_stat,get_subplot_position,get_subplot_position2,load_dem,plot_taylor_diagram,
        read_dem,get_hpc_command,least_square_fit,read_yaml,read_excel, write_excel,rtext,
        mklink,sindex,pindex,nindex,cindex,resize,savefig,pplot,blit_manager,read,add_xtick,
        get_qnode,modify_figure,parallel_jobs,fig_IFNO,ceqstate,subdomain_index,interp,
        nargout,pause,isnumber,ncfile)

   if os.path.exists(os.path.dirname(__file__)+'/pylibs/src'):
      import pylibs.src.schism_file as schism_file
   else:
      import src.schism_file as schism_file
   sys.modules['schism_file'] = schism_file
   from schism_file import (read_schism_hgrid, read_schism_bpfile,getglob,
        schism_grid,schism_vgrid,schism_bpfile,sms2grd,read_schism_vgrid,save_schism_grid,
        compute_zcor,read_schism_param,write_schism_param,read_schism_local_to_global,
        create_schism_vgrid,srank,grd2sms,scatter_to_schism_grid,delete_schism_grid_element,
        read_schism_prop,read_schism_reg,interp_schism_3d,get_schism_var_info,check_schism_ihot,
        read_schism_output,change_schism_param,get_schism_output_info,get_schism_grid_subdomain,
        get_schism_output_subset,combine_schism_hotstart,combine_icm_output,read_schism_slab,
        convert_schism_source,schism_view,schism_check,zcor_to_schism_grid,compute_schism_volume,
        read_schism_grid,schism_transect)

   if os.getenv('HOME')!=None:
       sys.path.append(os.getenv('HOME'))

   #sys.modules['loadz'] = mylib #in case oldmodule name used
   #sys.modules['read_schism_file'] = schism_file #in case oldmodule name used

   #import mpas_file
   #from mpas_file import (read_mpas_grid)

   #old module names
   sys.modules['pyUtility']=sys.modules['src']

   #---------------------------------------------------------------------
   #alias
   #---------------------------------------------------------------------
   from os.path import exists as fexist
   from numpy import array_equal as eq
   from src.mylib import savez as save_npz; mylib.save_npz=savez
   from src.mylib import zdata as npz_data; mylib.npz_data=zdata
   from src.mylib import least_square_fit as lsq; mylib.least_square_fit=lsq
   from src.mylib import move_figure as mvfig
   from src.mylib import modify_figure as mfig
   from src.mylib import find_cs as find_continuous_sections; mylib.find_continuous_sections=find_cs
   from src.mylib import read_mat as mat2npz
   from src.mylib import read_mat as convert_matfile
   from src.mylib import read_dem as convert_dem_format
   from src.mylib import load_dem as load_bathymetry
   from src.mylib import read_shapefile_data as read_shp
   from src.mylib import write_shapefile_data as write_shp
   from src.mylib import harmonic_analysis as HA
   from src.schism_file import read_schism_hgrid as read_hgrid
   from src.schism_file import read_schism_grid as read_grd
   from src.schism_file import read_schism_grid as grd
   from src.schism_file import read_schism_bpfile as read_bp
   from src.schism_file import change_schism_param as chparam
