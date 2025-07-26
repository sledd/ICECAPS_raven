from netCDF4 import Dataset
import glob
import datetime
import numpy as np

def load_netcdf(filepath, in_vars):
    """ open netcdf file, load variables (in_vars should be a list) and
        output dictionary of variables"""

    out_vars = {}

    with Dataset(filepath, mode = 'r') as open_netcdf:
        if len(in_vars)>0:
            for var in in_vars:
                out_vars[var] = open_netcdf.variables[var][:]
        else:
            for var in open_netcdf.variables.keys():
                out_vars[var] = open_netcdf.variables[var][:]
    
        try:
            out_vars['time_unit'] = open_netcdf.variables['time'].units
        except:
            None

    return out_vars


def datetimes_to_seconds(dates_in):
    """ convert datetime objects to seconds since 1970/1/1 (arbitrary)
    dates_in should be a list or array of datetime object
    """
    seconds_out = np.asarray([(t-datetime.datetime(1970,1,1)).total_seconds() for t in dates_in])
    return seconds_out


def load_raven_sleigh_data(varnames_list):

    asfs_filepath = '/psd3data/arctic/raven_process/asfs/2_level_product/'
    asfs_lev2_files = glob.glob(asfs_filepath+'seb.level2.0.melt-sleighsfs.10min.*.nc')
    
    asfs_lev2_files.sort()


    # lev1_slow_vars = ['up_short_hemisp_qc','up_long_hemisp_qc','down_short_hemisp_qc','down_long_hemisp_qc','subsurface_heat_flux_A_qc',
    #               'subsurface_heat_flux_B_qc','skin_temp_surface_qc','temp_qc','snow_depth_qc','zenith_true_qc',
    #               'down_short_diffuse','down_short_direct','up_short_hemisp','up_long_hemisp','down_short_hemisp','down_long_hemisp',
    #               'snow_depth','temp','brightness_temp_surface','skin_temp_surface','subsurface_heat_flux_A','subsurface_heat_flux_B',
    #               'subsurface_heat_flux_C','zenith_true','snow_gpr_dist','base_time','time']

    asfs_data_lev2 = {}

    for fname in asfs_lev2_files[:]:
        # print(fname)
    
        fdic = load_netcdf(fname, varnames_list)
        fstart_time = datetime.datetime.strptime(fdic['time_unit'], 'seconds since %Y-%m-%dT%H:%M:%S.000000')
        fdic['dates'] = np.asarray([fstart_time+datetime.timedelta(seconds=int(m)) for m in fdic['time']])

        # varnames_list.remove('base_time')
        # varnames_list.remove('time')
        
        for var in varnames_list+['dates']:
            if var in ['base_time','time']:
                continue
            if var not in asfs_data_lev2:
                asfs_data_lev2[var] = fdic[var]
            else:
                asfs_data_lev2[var] = np.ma.concatenate( (asfs_data_lev2[var], fdic[var]), axis=0 )

    if 'down_short_hemisp' in varnames_list and 'up_short_hemisp' in varnames_list:
        asfs_data_lev2['net_short_hemisp'] = asfs_data_lev2['down_short_hemisp'] - asfs_data_lev2['up_short_hemisp']
        asfs_data_lev2['albedo'] =  asfs_data_lev2['up_short_hemisp']/asfs_data_lev2['down_short_hemisp']
        
    if 'down_short_diffuse' in varnames_list and 'down_short_direct' in varnames_list:
        asfs_data_lev2['diffuse_frac'] = asfs_data_lev2['down_short_diffuse']/asfs_data_lev2['down_short_direct']

    with np.load('/home/asledd/ICECAPS/Raven_SW-correct-tskin_emis985_all-times_20250113.npz') as npz:
        asfs_data_lev2['skin_temp'] = np.ma.MaskedArray(**npz)

    return asfs_data_lev2