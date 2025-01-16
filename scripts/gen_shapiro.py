#!/usr/bin/env python3
#region files have WGS84 coordinate in this script
#you can find the examples of region files used in script here:/sciclone/data10/kpark07/regfiles

from pylib import *

def write_shapiro(
    grd,
    shapiro_max=0.5, threshold_slope=0.5,
    depths=None, shapiro_vals1=None,
    s_regions=None, s_shapiro_maxs=None, s_threshold_slopes=None,
    regions=None, shapiro_vals2=None, i_set_add_s=None,
    fname='shapiro.gr3'
 ):
    '''
    write shapiro fileter strength value based on depth and specfied regions

    Input:
        grd:     grid name (*.gr3 or *.npz, where *.npz is python format)
        depths:  two hgrid depth(m) to distinguish river, land and the transition zone
        mvalues: lower and upper limits of shapiro values

        regions: list of region file names. e.g. regions=('GoME_1.reg','GoME_2.reg')
        s_regions: list of subregions file names for different parameters of shapiro filter
        s_shapiro_maxs: maximum shapiro values for subregions
        s_threshold_slopes: slope threshold for subregions
        rvalues: list of values for each region:  e.g. (0.1,0.5)
        i_set_add_s: identifier for setting or adding value,  0: set value; 1: add value
    '''
    #read hgrid
    if grd.endswith('.npz'):
        gd=loadz(grd).hgrid
    else:
        gd=read_schism_hgrid(grd)

    # compute bathymetry gradient on each node
    _, _, slope = gd.compute_gradient(fmt=2)

    # compute shapiro coefficients
    shapiro=shapiro_max*tanh(2*slope/threshold_slope)
    if s_regions is not None:
       for s_shapiro_max, s_threshold_slope, s_region in zip(s_shapiro_maxs,s_threshold_slopes,s_regions):
           print(f'managing {s_shapiro_max} shapiro with {s_threshold_slope} in {s_region}')
           bp=read_schism_bpfile(s_region)
           px,py=proj_pts(bp.x,bp.y,'epsg:4326','epsg:26918')
           sind=inside_polygon(c_[gd.x,gd.y], px,py).astype('bool')
           shapiro[sind]=s_shapiro_max*tanh(2*slope[sind]/s_threshold_slope)
    # further tweaks on shallow waters
    if len(depths) != len(shapiro_vals1):
        raise Exception(f'lengths of depths {len(depths)} and shapiro_vals1 {len(shapiro_vals1)} inconsistent')
    fp = gd.dp < depths[-1]
    shapiro[fp] = maximum(shapiro[fp], interp(gd.dp[fp], depths, shapiro_vals1))

    #set or add values in regions
    if regions is not None:
        for i_set_add, rvalue, region in zip(i_set_add_s, shapiro_vals2, regions):
            bp=read_schism_bpfile(region)
            px,py=proj_pts(bp.x,bp.y,'epsg:4326','epsg:26918') 
            sind=inside_polygon(c_[gd.x,gd.y], px,py).astype('bool')
            
            if i_set_add==0:
                print(f'setting {rvalue} shapiro in {region}')
                fp=sind
                shapiro[fp]=rvalue
            else:
                print(f'adding {rvalue} shapiro in {region}')
                sum(sind)
                sind2=(gd.dp>depths[0])  # additional condition: deeper water, dp > -1 m
                fp=(sind & sind2)
                shapiro[fp]=shapiro[fp]+rvalue
    #Edit values that is higher than maximum shaprio value
    sind=(shapiro>shapiro_max)
    shapiro[sind]=shapiro_max
    #save shapiro.gr3
    gd.dp=shapiro
    gd.write_hgrid(fname)

if __name__=="__main__":
    outfilename = './shapiro.gr3'

    if os.path.exists(outfilename):
        os.remove(outfilename)

    write_shapiro(
        grd='./hgrid_utm.gr3',  # grid name (*.gr3 or *.npz, where *.npz is python format)
        shapiro_max=0.5,
        threshold_slope=0.5,
        depths=[-99999, 20, 50],  # tweaks in shallow waters
        shapiro_vals1=[0.2, 0.2, 0.05],  # tweaks in shallow waters
        s_regions=['./LD1.bp','./LD2.bp','./LD3.bp','./SD1.bp','./SD2.bp','./NM1.bp','./NM2.bp','./NM3.bp','./INC1.bp','./INC2.bp','./INC3.bp','./INC4.bp','./INC5.bp','./CA1.bp','./CA2.bp','./CA3.bp','./CA4.bp'],  # subregions for different shapiro parameters 
        s_shapiro_maxs=[0.005,0.005,0.005,\
                        0.3,0.3,\
                        0.5,0.5,0.5,\
                        0.5,0.5,0.5,0.5,0.5,\
                        0.5,0.5,0.5,0.5], # maximum shapiro values for subregions
        s_threshold_slopes=[0.8,0.8,0.8,\
                            0.7,0.7,\
                            0.5,0.5,0.5,\
                            0.5,0.5,0.5,0.5,0.5,\
                            0.5,0.5,0.5,0.5],# slope threshold for subregions
        regions=['./INC1.bp','./INC3.bp','./INC4.bp','./INC5.bp','./CA2.bp','./CA3.bp','./CA4.bp'],    # regions for set or add values
        shapiro_vals2=[0.15,0.2,0.15,0.15,0.1,0.1,0.1],  # tweaks in regions, the order matters
        i_set_add_s=[1,0,1,1,1,1,1,1],  # 0: set; 1: add
        fname=outfilename
    )

