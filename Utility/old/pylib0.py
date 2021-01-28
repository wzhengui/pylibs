#!/usr/bin/evn python3
from pylab2 import *
import os,sys
import scipy as sp
from scipy import (optimize,interpolate,io,signal)
from scipy.fftpack import fft
from pyproj import Proj, transform
import importlib as imp
from mpi4py import MPI
from netCDF4 import Dataset
from io import StringIO
import skill_metrics as sm
import re

#my library
from str2num import (str2num,remove_tail)
import convert_matfile_format as cmat
from date_proc import datenum
from loadz import loadz,npz_data,save_npz
from misc import (wipe,reload,smooth,clear_globals,DaytimeLength,move_figure,lpfilt,mdivide,signa,
     inside_polygon,command_outputs,near_pts,proj,close_data_loop,get_prj_file,
     mfft)
from shpfile import read_shapefile_data,write_shapefile_data
from read_schism_file import (read_schism_hgrid, read_schism_hgrid_ll,read_schism_bpfile,getglob,
     schism_grid,schism_bpfile,sms2gr3,read_schism_vgrid,read_schism_param,write_schism_param)
from netcdf import ReadNC, WriteNC

if os.getenv('HOME')!=None:
    sys.path.append(os.getenv('HOME'))
