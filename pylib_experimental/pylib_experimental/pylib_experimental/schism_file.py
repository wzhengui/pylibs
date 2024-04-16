#/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import copy
import unittest

from pylib import schism_grid, zdata, WriteNC


# ---------------------experimental functions and classes-------------------------------------

def cread_schism_hgrid(fname):
    '''
    A wrapper function to read SCHISM hgrid.gr3 file using c++ binding,
    then copy the data to pylib's grid object.
    Usefull when the grid is large and the python reading is slow.
    '''
    import hgrid_pybind

    hgrid_obj = hgrid_pybind.HGrid(fname, True)  # read c++ HGrid object from file and optionally get side information

    gd=schism_grid()  # initialize empty pylib's grid object
    gd.source_file = str(fname)

    # copy the member variables from the c++ HGrid object to pylib's grid object
    gd.np = hgrid_obj.np
    gd.ne = hgrid_obj.ne
    gd.x = hgrid_obj.x
    gd.y = hgrid_obj.y
    gd.dp = hgrid_obj.z
    gd.elnode = hgrid_obj.elements
    gd.i34 = hgrid_obj.i34
    gd.ns = hgrid_obj.ns

    # need more testing
    # gd.iobn = np.array(hgrid_obj.openBoundaryNodes, dtype=object)
    # gd.nobn = np.array(hgrid_obj.nobn)
    # gd.nob = np.array(hgrid_obj.nob)
    # gd.ilbn = np.array(hgrid_obj.landBoundaryNodes, dtype=object)
    # gd.nlbn = np.array(hgrid_obj.nlbn)
    # gd.nlb = np.array(hgrid_obj.nlb)
    # gd.island = np.array(hgrid_obj.island)

    return gd

def read_schism_vgrid_cached(vg_filename, overwrite_cache=False):
    vg_cache_fname = os.path.splitext(vg_filename)[0] + '.pkl'
    if not overwrite_cache:
        try:
            with open(vg_cache_fname, 'rb') as handle:
                vg = pickle.load(handle)
        except Exception as e:
            print(f'{e}\nError reading cache file {vg_cache_fname}.\nReading from original file.')
            # read original file
            vg = read_schism_vgrid(vg_filename)
            with open(vg_cache_fname, 'wb') as handle:
                pickle.dump(vg, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        vg = read_schism_vgrid(vg_filename)
        with open(vg_cache_fname, 'wb') as handle:
            pickle.dump(vg, handle, protocol=pickle.HIGHEST_PROTOCOL)
    return vg

def combine_dataframes(A, B, weights=[1.0, 1.0]):
    import numbers

    AB = copy.deepcopy(A)

    # warn if B's time period does not contain A's time
    if B.index[0] > A.index[0] or B.index[-1] < A.index[-1]:
        print('Warning: B\'s time period does not contain A\'s time period')
        print('In the interpolated B, NaN will be filled for the time period that is not covered by B')
        print('and the NaN in the interpolated B will be treated as 0 when summing up')

    # Interpolate B to A's index
    B_interpolated = B.reindex(A.index).interpolate(method='time')

    # Find the columns that are in B but not in A
    new_columns = B_interpolated.columns.difference(A.columns)

    # append these columns to A
    AB = pd.concat([A, B_interpolated[new_columns]], axis=1)

    # Find the common columns in A and B
    common_columns = A.columns.intersection(B_interpolated.columns)

    # Sum the values from both dataframes for the common columns
    if isinstance(weights[0], numbers.Number):
        AB[common_columns] = weights[0] * A[common_columns] + weights[1] * B_interpolated[common_columns]
    elif isinstance(weights[0], TimeHistory):
        vs_A = weights[0].df[common_columns]
        vs_B = weights[1].df[common_columns]
        vs_B = vs_B.reindex(vs_A.index).interpolate(method='time')
        weights_A = vs_A / (vs_A + vs_B)
        weights_B = vs_B / (vs_A + vs_B)
        AB[common_columns] = weights_A * A[common_columns] + weights_B * B_interpolated[common_columns]
    else:
        raise ValueError('weights must be a list of two numbers or two TimeHistory objects')

    return AB

class TimeHistory:
    """Class for handling SCHISM's *.th file format.
    A *.th file is a simple ascii file containing time series data,
    the 1st column is time in seconds or days,
    the 2nd to last column is data, so each column is a time series.
    However, the *.th file itself lacks information about the time units and start time
    and this class is to bind those information in a dataframe.
    The class also sets the time column as the index of the dataframe
    and assign each data column with a meaningful name (e.g., station name)
    to facilitate data queries and processing.
    """
    @property
    def time(self):
        # original time in *th's format
        seconds = (self.df.index - self.df.index[0]).total_seconds().values
        t = seconds / self.sec_per_time_unit
        return t

    @property
    def delta_t(self): # time step
        return self.time[1] - self.time[0]

    @property
    def n_time(self):
        return self.df.shape[0]

    @property
    def n_station(self):
        return self.df.shape[1]

    @property
    def stations(self):
        return self.df.columns.values.tolist()

    @property
    def data(self):
        # original data excluding time
        return self.df.values


    def __init__(self, start_time_str="2000-01-01 00:00:00", data_array=None, columns=None, th_unit='seconds'):
        """
        Initialize a TimeHistory object from a data array,
        assuming the first column is time, and the rest are data.
        Note that you need to provide the start time and the time unit.
        """

        # list of main attributes
        self.df = pd.DataFrame(data_array)
        self.meaningful_columns = False  # some functions only work when column names are meaningful,
        self.th_unit = th_unit
                                         # e.g., station names or element ids, but not datetime + column numbers
        self.sec_per_time_unit = None

        # set time unit
        unit_dict = {'seconds': 1, 'minutes': 60, 'hours': 3600, 'days': 86400, 'weeks': 604800, 'years': 31536000}
        self.sec_per_time_unit = unit_dict[th_unit]  # only related to the time column of a *.th file

        # set column names, which usually are datetime + station ids
        if type(columns) is list:  # from a user-specified list
            if len(columns) == self.df.shape[1]:
                self.df.columns = [str(x) for x in columns]
                self.meaningful_columns = True
            elif len(columns) == self.df.shape[1]-1:
                print('number of column labels does not match the data array, assuming the first column of the data is time')
                self.df.columns = [str(x) for x in ['datetime'] + columns]
                self.meaningful_columns = True
            else:
                raise Exception('number of columns does not match')
        elif columns is None:  # first col is time and the rest are column numbers
            self.df.columns = ['datetime'] + [str(x) for x in range(1, self.df.shape[1])]
            self.meaningful_columns = False
        else:
            raise Exception('unknown columns type')

        # Lay some ground rules
        # force the first column's name to "datetime"
        self.df.rename(columns={0: 'datetime'}, inplace=True)
        second_series = self.df['datetime'].values * self.sec_per_time_unit
        # convert the first column to pandas DatetimeIndex
        time_stamps = pd.to_timedelta(second_series, unit='s') + pd.to_datetime(start_time_str)
        self.df['datetime'] = time_stamps
        # set datetime as index
        self.df.set_index('datetime', inplace=True)
        # force column names to be string
        self.df.columns = self.df.columns.map(str)

    @classmethod
    def from_file(cls, file_name, start_time_str="2000-01-01 00:00:00", th_unit='seconds', columns=None):
        """
        Initialize from a file.
        Note that the *.th file doen't have information about the time units and start time,
        and columns names, so you need to provide them.
        """
        data = np.loadtxt(file_name)
        return cls(data_array=data, th_unit=th_unit, start_time_str=start_time_str, columns=columns)

    def __getitem__(self, selector):
        """Subset the TimeHistory object by column names"""

        # parse row and col selectors from selector
        if type(selector) is str:
            selector = [selector]
        elif isinstance(selector, np.ndarray):
            if len(selector.shape) == 1:  # 1D array of column names
                selector = selector.astype(str).tolist()
            else:
                raise IndexError("Column names must be a 1D array")

        if type(selector) is list and all(isinstance(x, str) for x in selector):  # subset by column names
            column_names = selector
            subset_data = np.array(self.df[column_names])
            return TimeHistory(
                start_time_str=self.df.index[0],
                data_array=np.c_[self.time, subset_data],
                columns=column_names,
                th_unit=self.th_unit
            )
        elif isinstance(selector, tuple):  # subset by row and column indices; column index does not include time
            if len(selector) != 2:
                raise IndexError("Only 2D indexing is supported")
            row_idx, col_idx = selector
            subset_data = self.data[row_idx, col_idx]
            subset_time_str = self.df.index[row_idx]
            subset_time = self.time[row_idx]
        elif isinstance(selector, slice):  # subset by row index only, i.e., by time
            subset_df = self.df.loc[selector]
            subset_data = subset_df.values
            subset_time_str = subset_df.index
            subset_time = np.array(subset_df.index - subset_df.index[0]).astype('timedelta64[s]').astype(float) / self.sec_per_time_unit
            col_idx = range(len(self.df.columns))

            # if isinstance(row_idx, slice) or isinstance(col_idx, slice):
            #     # Handle slices here
            #     rows = self.data[row_idx] if isinstance(row_idx, int) else [self.data[i] for i in range(*row_idx.indices(len(self.data)))]
            #     if isinstance(col_idx, int):
            #         return [row[col_idx] for row in rows]
            #     else:
            #         return [row[col_idx] for row in rows for col_idx in range(*col_idx.indices(len(row)))]
            # else:
            #     # Regular integer indexing
            #     data = self.data[row_idx][col_idx]

            return TimeHistory(
                start_time_str=subset_time_str[0],
                data_array=np.c_[subset_time, subset_data],
                columns=self.df.columns[col_idx].tolist()
            )
        else:
            raise IndexError("Unknown type of index")

    def __add__(self, other, weights=[1.0, 1.0]):
        """
        Add two TimeHistory objects together
        Interpolate other to self's time stamps;
        other attributes also inherit from self.
        When combining columns with the same name,
        the default weights are 1.0 for both self and other,
        meaning that the values are summed up.
        if you want to average the values, set weights=[0.5, 0.5]
        You can also specify two TimeHistory objects as the surrogates of the weights,
        e.g., when combining two sets of msource, you can set weights=[vsource1, vsource2],
        i.e., a weighted average based on the volume of the two sources,
        """

        A = copy.deepcopy(self)
        B = other

        A.df = combine_dataframes(A.df, B.df, weights=weights)

        return A

    def __eq__(self, other) -> bool:
        """Check if two TimeHistory objects are equal"""
        for att in ['sec_per_time_unit']:
            if getattr(self, att) != getattr(other, att):
                print(f'{att} is not equal')
                return False

        # dataframes, strict test
        # return self.df.equals(other.df)

        # test if the data are equal within a tolerance
        if np.any(self.df.columns != other.df.columns):
            print('column labels are not equal')
            return False
        return np.allclose(self.df, other.df, rtol=0.001, atol=0.0001)

    def writer(self, file_name, np_savetxt_args={'fmt':'%.4f', 'delimiter':' ', 'newline':'\n'}):
        # assemble data array in *.th format and write to file
        np.savetxt(file_name, np.c_[self.time, self.data], **np_savetxt_args)


class SourceSinkIn():
    def __init__(self, ele_groups=[[], []]):
        self.ele_groups = ele_groups  # 0: source; 1: sink

    @property
    def n_group(self):
        return len(self.ele_groups)

    @property
    def np_group(self):
        return [len(x) for x in self.ele_groups]

    @property
    def ip_group(self):
        return [np.array(x) for x in self.ele_groups]

    @property
    def n_source(self):
        return self.np_group[0]

    @property
    def n_sink(self):
        return self.np_group[1]

    @classmethod
    def from_file(cls, filename):
        ele_groups = [[], []]
        with open(filename, 'r') as file:
            for k in range(0, 2):  # 0: source; 1: sink
                num_points = int(file.readline().strip().split()[0])
                for _ in range(num_points):
                    ele_groups[k].append(int(file.readline().strip().split()[0]))
                # blank line between groups
                if k == 0:
                    file.readline()
        source_sink_in = cls(ele_groups=ele_groups)
        source_sink_in.print_info()

        return source_sink_in

    def print_info(self):
        print(f"nsource: {self.n_source}")
        if self.n_source > 0:
            print(f"first and last ele: {self.ele_groups[0][0]}, {self.ele_groups[0][-1]}")

        print(f"nsink: {self.n_sink}")
        if self.n_sink > 0:
            print(f"first and last ele: {self.ele_groups[1][0]}, {self.ele_groups[1][-1]}")

    def writer(self, filename=None):
        with open(filename, 'w') as fout:
            for k in range(0, self.n_group):
                print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                fout.write(f"{self.np_group[k]}\n")
                for i in range(0, self.np_group[k]):
                    fout.write(f"{self.ele_groups[k][i]}\n")
                fout.write("\n")  # empty line

    def __eq__(self, other) -> bool:
        for g1, g2 in zip(self.ele_groups, other.ele_groups):
            if g1.tolist() != g2.tolist():
                return False
        return True


class source_sink:
    """
    Class for handling all source/sink inputs:

    source_sink.in,
    vsource.th,
    vsink.th,
    msource.th

    These files are always used together,
    so a class is created to bind and process them,
    for example, adding or comparing two source_sink objects

    In addition, there are no complete time info in the *.th files,
    so the TimeHistory class is also used to bind the time info and
    provide additional functions.
    """
    @property
    def source_eles(self):  # element index starts from 1
        return self.source_sink_in.ele_groups[0]

    @property
    def sink_eles(self):  # element index starts from 1
        return self.source_sink_in.ele_groups[1]

    @property
    def nsource(self):
        return self.source_sink_in.n_source

    @property
    def nsink(self):
        return self.source_sink_in.n_sink

    def __init__(self, vsource:TimeHistory, vsink:TimeHistory, msource:list):
        """initialize from TimeHistory objects,
        vsource: TimeHistory object for volume source
        vsink: TimeHistory object for volume sink
        msource: list of TimeHistory objects for mass source

        Note that msource is different from SCHISM's native msource.th
        because saving all tracers in one array is not convenient for
        subsequent processing; instead, each tracer is saved in a separate
        TimeHistory object in a list
        """

        # list of main attributes
        self.vsource = vsource
        self.vsink = vsink
        self.msource = msource
        self.ntracers = None
        self.source_sink_in = None

        # if vsource, vsink, and msources are properly set,
        # then ntracers and source_sink_in can be decided without additional inputs
        if vsource is not None:
            self.ntracers = len(msource)
            source_eles = vsource.df.columns.astype(int).values  # index starts from 1
        else:
            self.ntracers = 0
            source_eles = []

        if vsink is not None:
            sink_eles = vsink.df.columns.astype(int).values
        else:
            sink_eles = []

        self.source_sink_in = SourceSinkIn(ele_groups=[source_eles, sink_eles])

        self.sanity_check()

    @classmethod
    def dummy(cls, start_time_str='2000-01-01 00:00:00', timestamps=[0.0, 86400.0*365*100],
                 source_eles=[], sink_eles=[], ntracers=2):
        """create a dummy source_sink object for testing purpose"""

        # initialize a set of source/sink files from scratch
        nt = len(timestamps)
        nsources = len(source_eles)
        nsinks = len(sink_eles)
        vsource = None
        vsink = None
        msource = [None] * ntracers

        if nsources > 0:
            vsource = TimeHistory(start_time_str=start_time_str,
                                  data_array=np.c_[np.array(timestamps), np.zeros([nt, nsources])],
                                  columns=source_eles)
            # dummy temperature, set to -9999, i.e., ambient temperature
            msource[0] = TimeHistory(start_time_str=start_time_str,
                                     data_array=np.c_[np.array(timestamps), -9999*np.ones([nt, nsources])],
                                     columns=source_eles)
            # dummy salinity, set to 0
            msource[1] = TimeHistory(start_time_str=start_time_str,
                                     data_array=np.c_[np.array(timestamps), np.zeros([nt, nsources])],
                                     columns=source_eles)

        if nsinks > 0:
            vsink = TimeHistory(start_time_str=start_time_str,
                                data_array=np.c_[np.array(timestamps), np.zeros([nt, nsinks])],
                                columns=sink_eles)

        return cls(vsource, vsink, msource)

    @classmethod
    def from_files(cls, source_dir, start_time_str='2000-01-01 00:00:00'):
        '''
        Initialize from existing source/sink files under the source_dir.
        Note that these files don't have start time information,
        and you need to provide it.
        '''
        vsource = None
        msource = None
        vsink = None

        # ele_groups are defined in source_sink.in
        source_sink_in = SourceSinkIn.from_file(filename=f"{source_dir}/source_sink.in")

        if source_sink_in.n_source > 0:
            print('reading vsource\n')
            vsource = TimeHistory.from_file(
                f"{source_dir}/vsource.th",
                start_time_str=start_time_str,
                columns=source_sink_in.ele_groups[0],  # source_eles, index starts from 1
            )

            print('reading msource\n')
            msource_total = np.loadtxt(f"{source_dir}/msource.th")
            msource_t = msource_total[:, 0]
            ntracers = (msource_total.shape[1] - 1) / vsource.n_station
            if (int(ntracers) != ntracers):
                raise ValueError("Number of tracers must be an integer, vsource and msource don't match")
            else:
                ntracers = int(ntracers)

            # Split msource_total into separate TimeHistory objects,
            # one for each tracer. This facilitates subsequent processing.
            msource = [None] * ntracers
            for i in range(ntracers):
                idx1 = i * vsource.n_station + 1
                idx2 = (i + 1) * vsource.n_station + 1
                msource[i] = TimeHistory(
                    data_array=np.c_[msource_t, msource_total[:, idx1:idx2]],
                    start_time_str=start_time_str,
                    columns=source_sink_in.ele_groups[0]
                )

        if source_sink_in.n_sink > 0:
            print('reading vsink\n')
            vsink = TimeHistory.from_file(
                f"{source_dir}/vsink.th",
                start_time_str=start_time_str,
                columns=source_sink_in.ele_groups[1]
            )

        source_sink = cls(vsource, vsink, msource)
        source_sink.sanity_check()
        return source_sink
    
    @classmethod
    def from_ncfile(cls, ncfile, start_time_str='2000-01-01 00:00:00'):
        '''
        Initialize from an existing source.nc file.
        Note that this file doesn't have start time information,
        and you need to provide it.
        '''
        import netCDF4 as nc

        # read source.nc
        with nc.Dataset(ncfile, 'r') as f:
            # source
            source_elem = f.variables['source_elem'][:]
            # vsource
            vs_data = f.variables['vsource'][:]
            dt_vs = f.variables['time_step_vsource'][:][0]
            time_vs = np.arange(0, dt_vs*vs_data.shape[0], dt_vs)
            vsource = TimeHistory(data_array=np.c_[time_vs, vs_data], start_time_str=start_time_str, columns=source_elem.tolist())
            # msource
            ms_data = f.variables[f'msource'][:]
            ntracers = ms_data.shape[1]
            dt_ms = f.variables[f'time_step_msource'][:][0]
            time_ms = np.arange(0, dt_ms*ms_data.shape[0], dt_ms)
            msources = [None] * ntracers
            for i in range(ntracers):
                msources[i] = TimeHistory(data_array=np.c_[time_ms, ms_data[:, i, :]], start_time_str=start_time_str, columns=source_elem.tolist())
            # vsink
            sink_elem = f.variables['sink_elem'][:]
            vsink_data = f.variables['vsink'][:]
            dt_vsink = f.variables['time_step_vsink'][:][0]
            time_vsink = np.arange(0, dt_vsink*vsink_data.shape[0], dt_vsink)
            vsink = TimeHistory(data_array=np.c_[time_vsink, vsink_data], start_time_str=start_time_str, columns=sink_elem.tolist())

        source_sink = cls(vsource, vsink, msources)
        source_sink.sanity_check()
        return source_sink

    def subset_by_time(self, start_time_str, end_time_str):
        '''
        Subset source/sink files by time.
        '''
        time_slice = slice(start_time_str, end_time_str)
        if self.vsource is not None:
            subset_vsource = self.vsource[time_slice]
            subset_msource = [x[time_slice] for x in self.msource]
        else:
            subset_vsource = None
            subset_msource = None

        if self.vsink is not None:
            subset_vsink = self.vsink[time_slice]
        else:
            subset_vsink = None

        return source_sink(subset_vsource, subset_vsink, subset_msource)

    def subset_by_idx(self, source_idx, sink_idx):
        '''
        Subset source/sink files by index.
        '''
        if self.vsource is not None:
            subset_vsource = self.vsource[:, source_idx]
            subset_msource = [x[:, source_idx] for x in self.msource]
        else:
            subset_vsource = None
            subset_msource = None

        if self.vsink is not None:
            subset_vsink = self.vsink[:, sink_idx]
        else:
            subset_vsink = None

        return source_sink(subset_vsource, subset_vsink, subset_msource)

    def subset_by_ele(self, source_eles=[], sink_eles=[]):
        '''subset source/sink files by element ids (index starts from 1)'''
        if self.vsource is not None:
            if source_eles == []:  # no subsetting
                subset_vsource = self.vsource
                subset_msource = self.msource
            else:
                subset_vsource = self.vsource[source_eles]
                subset_msource = [x[source_eles] for x in self.msource]
        else:
            subset_vsource = None
            subset_msource = None

        if self.vsink is not None:
            if sink_eles == []:
                subset_vsink = self.vsink  # no subsetting
            else:
                subset_vsink = self.vsink[sink_eles]
        else:
            subset_vsink = None

        return source_sink(subset_vsource, subset_vsink, subset_msource)


    def clip_by_polygons(self, hgrid, polygons_xy=[]):
        '''
        Select source/sink elements by polygons.
        An hgrid of schism_grid type is required to get element coordinates.
        polygons: a list of 2D np arrays of (x, y) coordinates of each polygon
        '''
        hgrid.compute_ctr()

        # select source and sink elements
        inside_source = np.zeros(self.nsource, dtype=bool)
        inside_sink = np.zeros(self.nsink, dtype=bool)
        for polygon_xy in polygons_xy:
            ele_xy = np.c_[hgrid.xctr[self.source_eles-1], hgrid.yctr[self.source_eles-1]]
            inside_source += inside_polygon(ele_xy, polygon_xy[:, 0], polygon_xy[:, 1]).astype(bool)

            ele_xy = np.c_[hgrid.xctr[self.sink_eles-1], hgrid.yctr[self.sink_eles-1]]
            inside_sink += inside_polygon(ele_xy, polygon_xy[:, 0], polygon_xy[:, 1]).astype(bool)

        inside_ss = self.subset_by_idx(inside_source, inside_sink)
        outside_ss = self.subset_by_idx(~inside_source, ~inside_sink)

        return inside_ss, outside_ss

    def writer(self, output_dir):
        '''
        Write source/sink files to the output_dir.
        '''
        os.makedirs(output_dir, exist_ok=True)

        self.source_sink_in.writer(f"{output_dir}/source_sink.in")
        if self.vsource is not None:
            self.vsource.writer(f"{output_dir}/vsource.th")

            msource_total = self.msource[0].time
            for i in range(self.ntracers):
                msource_total = np.c_[msource_total, self.msource[i].df.values]
            np.savetxt(f"{output_dir}/msource.th", msource_total)

        if self.vsink is not None:
            self.vsink.writer(f"{output_dir}/vsink.th")

        # additional outputs in *.nc format
        self.nc_writer(output_dir=output_dir)

    def diag_writer(self, hgrid, output_dir):
        '''writer for diagnostic files'''
        hgrid.compute_ctr()

        if self.vsource:
            np.savetxt(
                f'{output_dir}/sources.xyz',
                np.c_[hgrid.xctr[self.source_eles-1], hgrid.yctr[self.source_eles-1], self.vsource.df.mean().values]
            )
        if self.vsink:
            np.savetxt(
                f'{output_dir}/sinks.xyz',
                np.c_[hgrid.xctr[self.sink_eles-1], hgrid.yctr[self.sink_eles-1], self.vsink.df.mean().values]
            )

    def nc_writer(self, output_dir=None):
        if output_dir is None:
            raise Exception("output_dir is required.")
        os.makedirs(output_dir, exist_ok=True)

        # create netcdf data
        C=zdata(); C.vars=[]; C.file_format='NETCDF4'

        # create dummy source/sink if they are empty
        if self.nsource == 0:
            dummy_ss = source_sink.dummy(start_time_str=self.vsink.df.index[0], timestamps=self.vsink.time, source_eles=['1'], sink_eles=['1'], ntracers=2)
            vsource = dummy_ss.vsource
            msource = dummy_ss.msource
            vsink = self.vsink
            nsource = 1
            nsink = self.nsink
            ntracers = 2
        else:
            dummy_ss = source_sink.dummy(start_time_str=self.vsource.df.index[0], timestamps=self.vsource.time, source_eles=['1'], sink_eles=['1'], ntracers=self.ntracers)
            vsource = self.vsource
            msource = self.msource
            vsink = dummy_ss.vsink
            nsource = self.nsource
            nsink = 1
            ntracers = self.ntracers

        C.dimname=['nsources', 'nsinks', 'ntracers', 'time_msource','time_vsource','time_vsink','one']
        C.dims=[nsource, nsink, ntracers, msource[0].n_time, vsource.n_time, vsink.n_time, 1]

        C.vars.extend(['source_elem','vsource','msource'])
        vi=zdata(); vi.dimname=('nsources',); vi.val=vsource.df.columns.values.astype(int); C.source_elem=vi
        vi=zdata(); vi.dimname=('time_vsource','nsources'); vi.val=vsource.data; C.vsource=vi
        msource_data = np.stack([x.data for x in msource], axis=1)  # cast into a 3D array of shape (nt, ntracers, nsources)
        vi=zdata(); vi.dimname=('time_msource','ntracers','nsources'); vi.val=msource_data; C.msource=vi

        C.vars.extend(['sink_elem','vsink'])
        vi=zdata(); vi.dimname=('nsinks',); vi.val=vsink.df.columns.values.astype(int); C.sink_elem=vi
        vi=zdata(); vi.dimname=('time_vsink','nsinks',); vi.val=vsink.data; C.vsink=vi

        C.vars.extend(['time_step_vsource','time_step_msource','time_step_vsink'])
        vi=zdata(); vi.dimname=('one',); vi.val=vsource.delta_t; C.time_step_vsource=vi
        vi=zdata(); vi.dimname=('one',); vi.val=msource[0].delta_t; C.time_step_msource=vi
        vi=zdata(); vi.dimname=('one',); vi.val=vsink.delta_t; C.time_step_vsink=vi

        WriteNC(f'{output_dir}/source.nc', C)

    def sanity_check(self):
        # check consistency of source_sink_in and vsource/vsink/msource
        if self.vsource is None and self.vsink is None:
            raise ValueError("vsource and vsink cannot be both None")

        if self.vsource is not None:
            if len(self.msource) != self.ntracers:
                raise Exception('inconsistent number of tracers')
            if self.nsource != self.msource[0].n_station:
                raise Exception('inconsistent number of msource stations')
            if self.nsource != self.vsource.n_station:
                raise Exception('inconsistent number of vsource stations')
            if np.min(self.vsource.df.values, axis=None) < 0:
                raise Exception('vsource must be non-negative')

        if self.vsink is not None:
            if self.nsink != self.vsink.n_station:
                raise Exception('inconsistent number of sink stations')
            if np.max(self.vsink.df.values, axis=None) > 0:
                raise Exception('vsink must be non-positive')

    def __add__(self, other):
        '''
        Add source/sink other to source/sink self,
        retaining self's time stamps.
        '''
        A = self
        B = other

        # sanity check
        if A.nsource == 0 and B.nsource == 0 and A.nsink == 0 and B.nsink == 0:
            raise Exception('both source and sink are empty')

        # most cases are trivial unless both A and B have source
        if A.nsource == 0 and B.nsource == 0:  # neither has source
            vsource = None
            msource = None
        elif A.nsource == 0:  # only B has source
            vsource = B.vsource
            msource = B.msource
        elif B.nsource == 0:  # only A has source
            vsource = A.vsource
            msource = A.msource
        else:  # both have source
            vsource = A.vsource + B.vsource  # using TimeHistory.__add__
            msource = [None] * A.ntracers
            for i in range(A.ntracers):  # also using TimeHistory.__add__, but with weights
                msource[i] = A.msource[i].__add__(B.msource[i], weights=[A.vsource, B.vsource])

        # most cases are trivial unless both A and B have sink
        if A.nsink == 0 and B.nsink == 0:  # neither has sink
            vsink = None
        elif A.nsink == 0:  # only B has sink
            vsink = B.vsink
        elif B.nsink == 0:  # only A has sink
            vsink = A.vsink
        else:  # both have sink
            vsink = A.vsink + B.vsink  # using TimeHistory.__add__

        return type(self)(vsource=vsource, vsink=vsink, msource=msource)

    def __eq__(self, other):
        for att in ['source_sink_in', 'vsource', 'vsink']:
            if getattr(self, att) != getattr(other, att):
                print(f'{att} not equal')
                return False
        for i, [ms_A, ms_B] in enumerate(zip(self.msource, other.msource)):
            if ms_A != ms_B:
                print('msource {i} not equal')
                return False
        return True

# ---------------------------- unit test ----------------------------
class test_cread_schism_hgrid(unittest.TestCase):
    def test_read_schism_hgrid(self):
        print('\n\n*************** test_read_schism_hgrid ****************')

        gd0 = schism_grid('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I16/hg.gr3')
        gd = cread_schism_hgrid('/sciclone/schism10/feiye/STOFS3D-v7/Inputs/I16/hg.gr3')

        gd0.compute_bnd()
        gd.compute_bnd()

        self.assertEqual(gd0.ne, gd.ne)
        self.assertEqual(gd0.np, gd.np)
        self.assertEqual(gd0.ns, gd.ns)
        self.assertEqual(gd0.nob, gd.nob)
        
        print(gd)
        print(gd.x)
        print(gd.y)
        print(f'number of elements: {gd.ne}')
        print(f'number of nodes: {gd.np}')
        print(f'number of sides: {gd.ns}')
        print(f'number of open boundaries: {gd.nob}')
class test_add_source_sink(unittest.TestCase):
    def test_add_source_sink(self):
        print('\n\n*************** test_add_source_sink ****************')
        # read from prepared sample files
        A = source_sink.from_files('/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Source_sink/Test_data/source_sink_sample1')
        B = source_sink.from_files('/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Source_sink/Test_data/source_sink_sample2')
        AB = source_sink.from_files('/sciclone/data10/feiye/SCHISM_REPOSITORY/schism/src/Utility/Pre-Processing/STOFS-3D-Atl-shadow-VIMS/Pre_processing/Source_sink/Test_data/source_sink_sample12')

        A_add_B = A + B

        self.assertEqual(A_add_B, AB)

class test_combine_dataframes(unittest.TestCase):
    def setUp(self):
        # create some random time series for test
        index_A = pd.date_range(start='2023-01-01', end='2023-12-31', freq='D')
        index_B = pd.date_range(start='2023-01-01', end='2023-10-31', freq='3D')  # lower frequency
        self.df_A = pd.DataFrame(np.random.rand(len(index_A), 3), index=index_A, columns=['A', 'B', 'C'])
        self.df_B = pd.DataFrame(np.random.rand(len(index_B), 3), index=index_B, columns=['B', 'C', 'D'])

    def test_combine_dataframes(self):
        print('\n\n*************** test combine dataframes ****************')
        df_combined = combine_dataframes(self.df_A, self.df_B)

        print("First few rows of DataFrame A:")
        print(self.df_A.head())
        print("First few rows of DataFrame B:")
        print(self.df_B.head())
        print("First few rows of Combined DataFrame:")
        print(df_combined.head())

        # Check the shape of the resulting dataframe
        self.assertEqual(df_combined.shape[0], self.df_A.shape[0])  # check row number
        self.assertEqual(df_combined.shape[1], 4)  # check column number

        # Check the column names of the resulting dataframe
        self.assertTrue(all(np.isin(['A', 'B', 'C', 'D'], df_combined.columns)))  # check column names

        # Check that the values of the resulting dataframe are correct
        common_columns = self.df_A.columns.intersection(self.df_B.columns)
        for col in common_columns:
            self.assertTrue(all(np.isclose(df_combined[col].loc[self.df_B.index], self.df_A[col].loc[self.df_B.index] + self.df_B[col], atol=1e-5)))

if __name__ == "__main__":
    # run unit tests
    unittest.main()
