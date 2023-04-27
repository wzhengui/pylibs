import multiprocessing as mp
import pathlib

import numpy as np
from netCDF4 import Dataset

from pylib import inside_polygon, schism_grid

def subset_on_stack(stack, nidxs, eidxs):

    for fname in ncfiles:

        ds = Dataset(f'outputs/{fname}_{stack}.nc')

        #get dimensions size
        nLevels = ds.dimensions['nSCHISM_vgrid_layers'].size
        nMax_face_nodes = ds.dimensions['nMaxSCHISM_hgrid_face_nodes'].size
        one = 1
        two = 2

        fout = Dataset(f'./{path}/{fname}_{stack}.nc', 'w', format='NETCDF3_CLASSIC')
        fout.createDimension('nSCHISM_hgrid_node', gd.np)
        fout.createDimension('nSCHISM_hgrid_face', gd.ne)
        fout.createDimension('nSCHISM_hgrid_edge', gd.ns)
        fout.createDimension('nMaxSCHISM_hgrid_face_nodes', nMax_face_nodes)
        fout.createDimension('nSCHISM_vgrid_layers', nLevels)
        fout.createDimension('one', one)
        fout.createDimension('two', two)
        fout.createDimension('time', None)

        #time
        fout.createVariable('time', 'f', ('time',))
        fout['time'][:] = ds.variables['time'][:]
        fout['time'].i23d = 0

        if fname == 'out2d':
            fout.createVariable('SCHISM_hgrid_edge_nodes', 'i', ('nSCHISM_hgrid_edge', 'two'))
            fout['SCHISM_hgrid_edge_nodes'][:] = isidenode + 1

            fout.createVariable('SCHISM_hgrid_face_nodes', 'i', ('nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'))
            fout['SCHISM_hgrid_face_nodes'][:] = gd.elnode + 1

        for var in ds.variables:

            #dims = ds.variables[var].ndim

            if 'nSCHISM_hgrid_node' in ds.variables[var].dimensions:
                print(f'Processing on stack {stack}, var {var}.')
                dims = ds.variables[var].ndim
                if dims == 1:
                    fout.createVariable(var, ds.variables[var].dtype, ds.variables[var].dimensions)
                    fout.variables[var][:] = ds.variables[var][nidxs]
                    #for attr in ds.variables[var].attrs.keys():
                    #    fout.variables[var].setncattr(attr, ds.variables[var].attrs.get(attr))

                elif dims == 2:
                    fout.createVariable(var, ds.variables[var].dtype, ds.variables[var].dimensions)
                    fout.variables[var][:, :] = ds.variables[var][:][:, nidxs]

                elif dims == 3:
                    fout.createVariable(var, ds.variables[var].dtype, ds.variables[var].dimensions)
                    fout.variables[var][:, :, :] = ds.variables[var][:][:, nidxs, :]

            elif 'nSCHISM_hgrid_face' in ds.variables[var].dimensions and 'time' in ds.variables[var].dimensions:
                    fout.createVariable(var, ds.variables[var].dtype, ds.variables[var].dimensions)
                    fout.variables[var][:, :] = ds.variables[var][:][:, eidxs]
            else:
                continue

            for attr in ds.variables[var].ncattrs():
                fout.variables[var].setncattr(attr, ds.variables[var].getncattr(attr))
        ds.close()
        fout.close()

if __name__ == '__main__':

    #input 1: stack
    stack_start = 15
    stack_end = 42

    #input 2: bbox
    lon_min = -92.0  #-82.0
    lon_max = -88.0  #-80.0
    lat_min = 29   #25.0
    lat_max = 31   #27.5

    #input 3: choose which nc files to be subsetted 
    #ncfiles = ['out2d', 'zCoordinates', 'salinity', 'temperature', 'horizontalVelX', 'horizontalVelY']
    ncfiles = ['out2d', 'zCoordinates', 'horizontalVelX', 'horizontalVelY']
    #ncfiles = ['out2d', 'zCoordinates']
    #ncfiles = ['horizontalVelX', 'horizontalVelY']

    #input 4: directory to save subsetting results
    dirname = 'subset_NC'
    path = pathlib.Path(dirname)
    if path.exists():
        print('Directory exists!')
    else:
        print(f'Create a new directory {dirname}')
        path.mkdir(parents=True, exist_ok=True)

    #input 5: save sub-grid
    save_subgrid = True

    #build polygon
    px = np.array([lon_min, lon_max, lon_max, lon_min])
    py = np.array([lat_max, lat_max, lat_min, lat_min])

    #build new hgrid
    ds = Dataset(f'./outputs/out2d_1.nc')

    gd = schism_grid()

    gd.elnode = ds.variables['SCHISM_hgrid_face_nodes'][:]-1
    gd.x = ds.variables['SCHISM_hgrid_node_x'][:]
    gd.y = ds.variables['SCHISM_hgrid_node_y'][:]
    gd.dp = ds.variables['depth'][:]
    gd.ne = len(gd.elnode)
    gd.np = len(gd.x)
    gd.i34 = np.ones(gd.ne).astype('int')
    fp3 = gd.elnode[:, -1] < 0
    gd.i34[fp3] = 3
    gd.i34[~fp3] = 4

    gd.compute_ctr()

    #indexes of elements inside the box
    eidxs = np.nonzero(inside_polygon(np.c_[gd.xctr, gd.yctr], px, py) == 1)[0]
    gd.elnode = gd.elnode[eidxs]

    nidxs = np.unique(gd.elnode)
    fpn = nidxs>=0
    nidxs = nidxs[fpn]
    gd.x = gd.x[nidxs]
    gd.y = gd.y[nidxs]
    gd.dp = gd.dp[nidxs]
    gd.ne = len(eidxs)
    gd.np = len(nidxs)
    gd.i34 = gd.i34[eidxs]

    #construct new element connectivity
    node2node = dict(zip(nidxs, np.arange(gd.np)))
    for i in np.arange(gd.elnode.size):
        if gd.elnode.ravel()[i]<0: continue
        gd.elnode.ravel()[i] = node2node[gd.elnode.ravel()[i]]

    #compute side
    gd.ns, isidenode, isdel = gd.compute_side(fmt=1)
    
    #save grid
    if save_subgrid:
        gd.save('hgrid_sub.gr3')

    ds.close()

    #create stack list
    stacks = [i for i in np.arange(stack_start, stack_end+1)]
    npool = len(stacks) if len(stacks) < mp.cpu_count() else mp.cpu_count()-1
    #npool = len(stacks) if len(stacks) < mp.cpu_count()/2 else int(mp.cpu_count()/2)
    print(f'npool is {npool}')

    pool = mp.Pool(npool)
    pool.starmap(subset_on_stack, [(i, nidxs, eidxs) for i in stacks])
    #pool.join()
    pool.close()
    del pool

