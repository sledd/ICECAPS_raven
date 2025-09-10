from netCDF4 import Dataset
import glob
import datetime
import numpy as np
from scipy import interpolate
import pandas as pd

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


def open_corrected_simba_temps(file_version='v4.3_20250815_sfcdif+4.0_ksi16.0'):
    acceptable_file_versions = ['v3_20250724', 'v2_20250723', 'v3_20250724_sfcdif+2', 'v2_20250723_sfcdif+2',
                               'v4_20250726_sfcdif+4.0', 'v4_20250726_sfcdif+0.0',
                               'v4.1_20250729_sfcdif+4.0', 'v4.1_20250729_sfcdif+0.0',
                               'v4.2_20250730_sfcdif+4.0_ksi20.0', 'v4.2_20250730_sfcdif+4.0_ksi15.0',
                               'v4.2_20250730_sfcdif+6.0_ksi15.0',
                               'v4.3_20250815_sfcdif+4.0_ksi15.0', 'v4.3_20250815_sfcdif+6.0_ksi15.0',
                                'v4.3_20250815_sfcdif+4.0_ksi16.0']

    if file_version not in acceptable_file_versions:
        print("Please choose of these versions: ", acceptable_file_versions)
        return 
        
    varnames = ['temp','initial_height','time']
    
    if file_version=='v2_20250723':
        fn = '/home/asledd/ICECAPS/firnprofile_solar-correction_rolling-variance-max-sfc_created_20250722.nc'
    elif 'v4.1_' in file_version:
        fn = f'/home/asledd/ICECAPS/firnprofile_solar-correction_rolling-variance-max-sfc-{file_version}_created_20250729.nc'
        varnames += ['detected_surface']
    elif 'v4.2_' in file_version or 'v4.3_' in file_version:
        fn = f'/home/asledd/ICECAPS/firnprofile_solar-correction_rolling-variance-max-sfc-{file_version}_created_20250730.nc'
        varnames += ['detected_surface']
    elif 'v4_' in file_version:
        fn = f'/home/asledd/ICECAPS/firnprofile_solar-correction_rolling-variance-max-sfc-{file_version}_created_20250726.nc'
    else:
        fn = f'/home/asledd/ICECAPS/firnprofile_solar-correction_rolling-variance-max-sfc-{file_version}_created_20250725.nc'
    
    season_data_solar = {}

    fdic = load_netcdf(fn, varnames)
    
    fdic['dates'] = np.asarray([datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(s)) for s in fdic['time']])
    season_data_solar['seconds'] = fdic['time']

    new_varnames = {'temp':'temperature', 'dates':'dates','initial_height':'height', 'detected_surface':'snow_surface'}
    for varname in fdic.keys():
        if varname in new_varnames:
            season_data_solar[new_varnames[varname]] = fdic[varname]
    

    return season_data_solar


def interpolate_and_mask(var_data, var_secs, goal_secs):
    """
    interpolate var_data to have same time steps as goal_secs; mask data outside of bounds
    """
    var_data = np.ma.masked_invalid(var_data)
    f_var = interpolate.interp1d(var_secs[~var_data.mask], var_data[~var_data.mask], fill_value=-999, bounds_error=False)
    interped_var = f_var(goal_secs)
    interped_var = np.ma.masked_equal(interped_var, -999)
    return interped_var

def make_temps_rel_to_sfc(heights_rel_to_simba, temps_rel_to_simba, sfc_rel_to_simba):

    depths_rel_to_sfc = np.repeat(heights_rel_to_simba.copy()[np.newaxis,:], temps_rel_to_simba.shape[0], axis=0)
    sfc_2d = np.ma.repeat(sfc_rel_to_simba[:,np.newaxis], temps_rel_to_simba.shape[1], axis=1)
    depths_rel_to_sfc = depths_rel_to_sfc-sfc_2d


    temporary_temp = np.ma.masked_where(depths_rel_to_sfc>0, temps_rel_to_simba) # temporary temperature that's only at or below the surface

    m = 100. # I don't remember why this is 100!
    temps_rel_to_sfc = []
    
    for i in range(temporary_temp.shape[0]):
        temp_i = temporary_temp[i,:-1]
        sub_t = temp_i[~temp_i.mask]
    
        sub_t_filled = np.ma.concatenate( (sub_t, np.full(int(m-sub_t.count()), -999 )) )
        temps_rel_to_sfc.append(sub_t_filled)
        
    temps_rel_to_sfc = np.ma.masked_equal(temps_rel_to_sfc, -999)

    ## new height variable that decreases downward
    sub_depths = np.arange(0,temps_rel_to_sfc.shape[1])
    sub_depths *= -2

    return sub_depths, np.ma.masked_invalid(temps_rel_to_sfc)


def get_simba_time_of_day(dates_in):
    simba_hours = np.array([d.hour for d in dates_in])
    simba_mins = np.array([d.minute for d in dates_in])
    
    ## time of day
    simba_tod = simba_hours+simba_mins/60.
    
    return simba_tod

def get_diurnal_cycle(dates_in, data_in, stat='mean'):
    
    tod = get_simba_time_of_day(dates_in)
    if stat=='mean':
        daily_data_out = np.array([data_in[np.ma.where(tod==t)[0],:].mean(axis=0) for t in np.unique(tod)])
    elif stat=='std':
        daily_data_out = np.array([data_in[np.ma.where(tod==t)[0],:].std(axis=0) for t in np.unique(tod)])

    return daily_data_out



def load_uncorrected_simba():
    simba_filepath = '/psd3data/arctic/sledd/Raven_simba/'
    fn = 'firnprofile-combined-corrected-at38.sled.level2.beta.15min.20240517-20240821_created_June26_2025.nc'
    
    season_data_unc = {}
    if 'firnprofile-combined-corrected' in fn:
        varnames = ['temperature','initial_height','time']
        fdic = load_netcdf(simba_filepath+fn, varnames)
    else:
        varnames = ['temp','initial_height','time']
        fdic = load_netcdf(fn, varnames)
    
    if 'firnprofile-combined-corrected' in fn:
        fstart_time = datetime.datetime.strptime(fdic['time_unit'], 'minutes since %Y-%m-%d')
        fdic['dates'] = np.asarray([fstart_time+datetime.timedelta(minutes=int(m)) for m in fdic['time']])
    else:    
        fdic['dates'] = np.asarray([datetime.datetime(1970,1,1)+datetime.timedelta(seconds=int(s)) for s in fdic['time']])
    season_data_unc['seconds'] = fdic['time']
    
    for var in ['temperature','dates','initial_height']:
        if var=='temperature' and fn=='firnprofile_maxheight_solar-correction_created_20250417.nc':
            season_data_unc['temperature'] = fdic['temp']
        else:
            season_data_unc[var] = fdic[var]
    
    if 'firnprofile-combined-corrected' in fn:
        season_data_unc['height'] = season_data_unc['initial_height'][::-1]
        season_data_unc['temperature'] = season_data_unc['temperature'][:,::-1]
    else:
        season_data_unc['height'] = season_data_unc['initial_height']
        
    return season_data_unc

def load_ogre_surfaces():
    ogre_surface_out = {}

    for ogre in ['rav1_south','rav2_origin','rav3_west']:
        df = pd.read_csv('/psd3data/arctic/sledd/Raven_OGRE/'+ogre+'.csv')
        
        ogre_date = np.array([dt.strptime(d, '%m/%d/%y %H:%M') for d in df['Datetime'].values])       
        ogre_sfc = df[' Average Ant. Height [m]'].values
        ogre_surface_out[ogre] = {'dates':ogre_date, 'sfc':ogre_sfc}
    return ogre_surface_out

def detrend_temps_rel_to_sfc(temps_rel_to_sfc):
    x = np.arange(temps_rel_to_sfc.shape[0])

    slopes_with_depth = []
    
    for i in range(temps_rel_to_sfc.shape[1]):
        y = temps_rel_to_sfc[:,i]
    
        try:
            coeff = np.ma.polyfit(x, y, 2)
        except:
            detrend_temps = np.ma.concatenate( (detrend_temps, np.full(x.shape, np.nan)[:, np.newaxis]), axis=1 )
            slopes_with_depth.append(np.nan)
            continue
        # no need to use the original x values here just for visualizing the polynomial
        x_poly = np.linspace(x.min(), x.max())
        y_poly = np.polyval(coeff, x_poly)
        
        # we need the original x values here, so we can remove the trend from all points
        trend = np.polyval(coeff, x)
        t_detrend = y - trend
        # note that simply subtracting the trend might not be enough for other data sets
        slopes_with_depth.append(coeff[0])
            
        if i==0:
            detrend_temps = t_detrend[:, np.newaxis]
        else:
            detrend_temps = np.ma.concatenate( (detrend_temps, t_detrend[:, np.newaxis]), axis=1 )
    
    
    return np.ma.masked_invalid(detrend_temps)

def interpolate_snowpit_densities(init_density, final_density, out_seconds):

    ## datetimes are from install and demob
    st_seconds = datetimes_to_seconds([datetime.datetime(2024,5,15)])
    en_seconds = datetimes_to_seconds([datetime.datetime(2024,8,19)])

    f_density = interpolate.interp1d([st_seconds[0], en_seconds[0]],[init_density, final_density])
    density_timeseries = f_density(out_seconds)

    return density_timeseries

def calculate_storage(temps_rel_to_sfc, seconds_in, sub_depths, d_layer=20):

    c_i = 2.1*10**-3 # MJ kg-1 C-1
    d = 0.02 # m
    
    # from snow pits at install and demob
    # now averaged in upper 20/30 cm
    # should just the average be used instead of range?
    if d_layer==20:
        rho_lo = interpolate_snowpit_densities(238, 229, simba_seconds)
        rho_hi = interpolate_snowpit_densities(328, 379, simba_seconds)
    elif d_lay==30:
        rho_lo = interpolate_snowpit_densities(242, 226, simba_seconds)
        rho_hi = interpolate_snowpit_densities(325, 387, simba_seconds)

    else:
        print('go calculate average density over new d_layer or choose 20 or 30 cm')
        return
    avg_temps = np.ma.copy(temps_rel_to_sfc)[:,np.where(sub_depths>=-d_layer)[0]].mean(axis=1)
    avg_temp_diff = avg_temps[2:] - avg_temps[:-2]
    time_diff = 2.*(seconds_in[2:] - seconds_in[:-2])
    
    # 10^6 to convert from MJ to J
    storage_avg_lo = 10.**6*c_i*rho_lo[1:-1]*0.2*avg_temp_diff/time_diff
    storage_avg_hi = 10.**6*c_i*rho_hi[1:-1]*0.2*avg_temp_diff/time_diff

    return {'high_density':storage_avg_hi, 'low_density':storage_avg_lo}


def get_gpr_melt():
    df = pd.read_csv('/psd3data/arctic/sledd/Melt_Index_despiked_GPR_Raven2024.csv')
    
    gpr_dates = np.array([datetime.datetime.strptime(d, '%Y-%m-%d %H:%M:%S') for d in df['date'].values])       
    gpr_melt = df['Melt_Index_Normal'].values
    return {'dates':gpr_dates, 'melt':gpr_melt}

def get_gpr_melt_dates(gpr_melt_data, simba_dates, melt_threshold=0.8):

    melty_dates = np.ma.masked_where(gpr_melt_data<melt_threshold, simba_dates)
    melty_day_mn = np.unique([' '.join((str(d.day), str(d.month))) for d in melty_dates[~melty_dates.mask]])
    all_day_mn = np.asarray([' '.join((str(d.day), str(d.month))) for d in simba_dates])
    
    ## inefficient but whatever
    simba_melty_indices = np.zeros(simba_dates.shape)
    
    for md in melty_day_mn:    
        simba_melty_indices[np.where(all_day_mn==md)] = 1

    return simba_melty_indices

def get_subset_dTdz_cycle(melt_indices, data_in, melt_or_dry, dates_in):
    
    melt_masking = {'melt':1, 'dry':0}

    if melt_or_dry not in melt_masking.keys():
        print('Give melt or dry to determine subset')
        return

    melt_val = melt_masking[melt_or_dry]
    
    data_in_subset = data_in[np.where(melt_indices==melt_val)[0],:]
    simba_dates_subset = dates_in[np.where(melt_indices==melt_val)[0]]
    
    daily_avgs_subset = get_diurnal_cycle(simba_dates_subset, data_in_subset, stat='mean')
    return daily_avgs_subset
