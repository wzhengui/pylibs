#!/usr/bin/evn python3
from pylib import *


def ReadNC(fname,med=1,order=0):
    #ReadNC(fname,med=1)
    #if med=1: return netcdf.Dateset(fname)
    #if med=2: reorgnaized Dataset similar to npz_data
    #order=1: only works for med=2
    #order=0: variable dimension order read not changed for python format
    #order=1: variable dimension order read reversed follwoing in matlab/fortran format
    C=Dataset(fname);

    if med==1:
        return C
    else:
        ncdims=[i for i in C.dimensions]
        ncvars=[i for i in C.variables]
        F=npz_data();
        F.file_format=C.file_format
        F.dimname=ncdims
        F.dims=[C.dimensions[i].size for i in ncdims]
        F.vars=ncvars
        for i in ncvars:
            fi=npz_data();
            dimi=C.variables[i].dimensions;
            fi.dimname=dimi
            fi.dims=[C.dimensions[j].size for j in dimi]
            fi.val=C.variables[i][:]
            fi.attrs=C.variables[i].ncattrs()
            for j in C.variables[i].ncattrs():
                ncattri=C.variables[i].getncattr(j);
                exec('fi.{}=ncattri'.format(j))

            if order==1:
                fi.dimname=list(flipud(fi.dimname))
                fi.dims=list(flipud(fi.dims))
                nm=flipud(arange(ndim(fi.val)));
                fi.val=fi.val.transpose(nm)

            exec('F.{}=fi'.format(i))

        return F

def WriteNC(C,fname,med=1,order=0):
    #WriteNC(C,fname,med=1)
    #C is data source
    #if med=1, C has netcdf.Dataset format
    #if med=2, C has different format
    #order=0: variable dimension order written not changed for python format
    #order=1: variable dimension order written reversed follwoing in matlab/fortran format
    if med==1:
        #----write NC files-------------
        fid=Dataset(fname,'w',format=C.file_format); #C.file_format
        ncdims=[i for i in C.dimensions]
        ncvars=[i for i in C.variables]
        for dimi in ncdims:
            fid.createDimension(dimi,C.dimensions[dimi].size)
        if order==0:
            for vari in ncvars:
                vid=fid.createVariable(vari,C.variables[vari].dtype,C.variables[vari].dimensions)
                for attri in C.variables[vari].ncattrs():
                    vid.setncattr(attri,C.variables[vari].getncattr(attri))
                fid.variables[vari][:]=C.variables[vari][:]
        elif order==1:
            for vari in ncvars:
                vid=fid.createVariable(vari,C.variables[vari].dtype,flipud(C.variables[vari].dimensions))
                for attri in C.variables[vari].ncattrs():
                    vid.setncattr(attri,C.variables[vari].getncattr(attri))
                nm=flipud(arange(ndim(C.variables[vari][:])));
                fid.variables[vari][:]=C.variables[vari][:].transpose(nm)

        fid.close()
    else:
        #----write NC files-------------
        fid=Dataset(fname,'w',format=C.file_format); #C.file_format
        for i in range(len(C.dims)):
            fid.createDimension(C.dimname[i],C.dims[i])

        if order==0:
            for vari in C.vars:
                vi=eval('C.{}'.format(vari));
                vid=fid.createVariable(vari,vi.val.dtype,vi.dimname)
                for j in vi.attrs:
                    attri=eval('vi.{}'.format(j))
                    vid.setncattr(j,attri)
                fid.variables[vari][:]=vi.val
        elif order==1:
            for vari in C.vars:
                vi=eval('C.{}'.format(vari));
                vid=fid.createVariable(vari,vi.val.dtype,flipud(vi.dimname))
                for j in vi.attrs:
                    attri=eval('vi.{}'.format(j))
                    vid.setncattr(j,attri)
                if ndim(vi.val)>=2:
                    nm=flipud(arange(ndim(vi.val)));
                    fid.variables[vari][:]=vi.val.transpose(nm)
                else:
                    fid.variables[vari][:]=vi.val


        fid.close()

if __name__=='__main__':

    F2=ReadNC(r'D:\Work\E3SM\E3SMScript\run4ie\sflux\sflux_air_1.002.nc',2)
    WriteNC(F2,'T2.nc',2)

    F=ReadNC('T2.nc')

#    pass
#    # read NC files
#    C=Dataset('sflux_air_1.002.nc')
#
#    ncdims=[i for i in C.dimensions]
#    ncvars=[i for i in C.variables]
#
#
#    [print("{}".format(i)) for i in ncdims]
#    [print("{}".format(i)) for i in ncvars]
#
#    #----write NC files-------------
#    fid=Dataset('test.nc','w',format='NETCDF3_CLASSIC'); #C.file_format
#
#    for dimi in ncdims:
#        fid.createDimension(dimi,C.dimensions[dimi].size)
#
#    for vari in ncvars:
#        vid=fid.createVariable(vari,C.variables[vari].dtype,C.variables[vari].dimensions)
#        for attri in C.variables[vari].ncattrs():
#           vid.setncattr(attri,C.variables[vari].getncattr(attri))
#        fid.variables[vari][:]=C.variables[vari][:]
#    fid.close()
#
#    ## check results
#    F=Dataset('test.nc');
