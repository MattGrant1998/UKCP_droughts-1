import iris
import numpy as np
# import iris.experimental.equalise_cubes
import pre_processing
import os
import numpy.ma as ma
import matplotlib.pyplot as plt
import iris.coord_categorisation as cat
import iris.quickplot as qplt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import colors
import iris.plot as iplt
import cartopy.crs as ccrs
import csv
from iris.coord_categorisation import add_categorised_coord
import time

# move index_numbers and period names to pre_processing - will require some refactoring
INDEX_NUMBERS = [
    3,
    6,
    12
]

PERIOD_NAMES = {
    3: 'trimester',
    6: 'semester',
    12: '12-month',
}

INDEX_NAME = ['DSI-3', 'DSI-6', 'DSI-12']


def equalise_coordinate_bounds(cubelist, coord):
    cube1 = cubelist[0]
    for cube in cubelist[1:]:
        cube.coord(coord).bounds = cube1.coord(coord).bounds

    return cubelist

def calculate_average_climatology(climatology_cube):
    """
    Computes the average climatology over the 1980-2000 period, for each month.

    Args:
         climatology_cube (iris.cube.Cube):
            cube of precipitation from 1980-2000

    Returns:
        average_climatology (iris.cube.Cube):
            cube with each entry the average precipitation for the
            respective month and grid point over the time period 1980-2000
    """
    average_climatology_list = iris.cube.CubeList()

    for i in range(12):
        climatology_month_i = climatology_cube[i::12, ...]
        average_climatology_i = climatology_month_i.collapsed('time', iris.analysis.MEAN)
        average_climatology_list.append(average_climatology_i)

    average_climatology = average_climatology_list.merge_cube()

    return average_climatology


def add_n_month_co(period_length, period_name, cube):
    # This function adds n-month-period coordinate to the cube
    # Cube is aggregated by means over rolling window of length period_length (e.g. semester means)
    # e.g. if period length is 6, add coordinate semester_number to cube
    # where semester N correspond to the window ending in month N
    # e.g. 1 = period ending in Jan = August to Jan

    year_period = int(cube.coord('year').points[-1] - cube.coord('year').points[0])
    pts = np.arange(period_length, 13).tolist() + np.arange(1, 13).tolist() * (year_period)

    newcoord = iris.coords.AuxCoord(pts, standard_name=None, long_name=None, var_name=None,
                                    units='1', bounds=None, attributes=None, coord_system=None)

    cube.add_aux_coord(newcoord, 0)
    cube.coord('unknown').rename(period_name + '_number')


def add_year_co(period_length, cube):
    # This function adds integer year coordinate to substitute rational ones.
    # Cube is aggregated by means over rolling window of length period_length (e.g. semester means)
    # Period window is associated to a year if it is ending in such year (e.g. (August '67, Jan '68) has semester coordinate 1 and year coordinate '68 )

    pts = cube.coord('year').points
    integer_years = np.modf(pts)[0] == 0
    pts[~integer_years] = np.ceil(pts[~integer_years])

    newcoord = iris.coords.AuxCoord(pts, standard_name=None, long_name=None, var_name=None,
                                    units='1', bounds=None, attributes=None, coord_system=None)

    cube.add_aux_coord(newcoord, 0)
    cube.coord('unknown').rename('year_int')


def get_anomaly2(cube, clim, period_length, period_name):
    # Input:,
    # cube: 3d cube with aux coord year and month_number,
    # clim: numpy array of length 12 with the climatological averages for every n-period ending in a month
    # period_length: integer larger or equal to 1
    # period_name: string

    # Output: n_monthly mean anomaly cube with respect to the given climatology

    lastyear = int(cube.coord('year').points[-1])
    firstyear = int(cube.coord('year').points[0])
    no_of_years = lastyear - firstyear

    if period_length > 1:  # n-monthly mean anomaly

        # time dimension this cube= original number months - window length +1
        roll_mean = cube.rolling_window('month_number', iris.analysis.MEAN, period_length)

        # add semester as aux_coord to the cubes of interest:
        add_n_month_co(period_length, period_name, roll_mean)

        # add integer year as aux_coord to the cubes of interest:
        add_year_co(period_length, roll_mean)

        first_year_clim = clim[period_length - 1:]
        other_years_clim = np.tile(clim, (no_of_years - 1, 1, 1))
        all_years_clim = np.concatenate((first_year_clim, other_years_clim), axis=0)
        anomaly = roll_mean - all_years_clim
        # # define the anomaly for the first year (they are only 12-n+1 months)
        # anomaly0 = [roll_mean[0:12 - period_length + 1] - clim[period_length - 1:]]
        #
        # for y in range(firstyear + 1, lastyear + 1):
        #     anomaly0 += [roll_mean.extract(iris.Constraint(year_int=lambda cell: cell == y)) - clim]

    elif period_length == 1:
        # compute monthly precipitation anomaly wrt reference period
        all_years_clim = np.tile(clim, (no_of_years, 1, 1))
        anomaly = cube - all_years_clim
        # for y in range(firstyear, lastyear + 1):
        #     year_constraint = iris.Constraint(year=y)
        #     cube_year = cube.extract(year_constraint)
        #
        #     # clim.data runs from dec to nov, so need to roll it back so it goes jan to dec
        #     clim_data = np.roll(clim.data, -1)
        #
        #     if y == firstyear:
        #         clim_year = clim_data[11]
        #     elif y == lastyear:
        #         clim_year = clim_data[:11]
        #     else:
        #         clim_year = clim_data
        #
        #     # extract the year and take difference with reference mean
        #     anomaly_cube_y = cube_year - clim_year
        #
        #     if y == firstyear:
        #         scalar_time_coord = anomaly_cube_y.coord('time')
        #         dim_time_coord = iris.coords.DimCoord(scalar_time_coord.points,
        #                                               standard_name=scalar_time_coord.standard_name,
        #                                               units=scalar_time_coord.units)
        #         anomaly_cube_y.remove_coord(scalar_time_coord)
        #         anomaly_cube_y.add_aux_coord(dim_time_coord, 0)
        #         iris.util.promote_aux_coord_to_dim_coord(anomaly_cube_y, 'time')
        #         print(anomaly_cube_y)
        #
        #     anomaly0.append(anomaly_cube_y)
        # print(cube_year - clim_year)

    else:
        print('Period length should be a positive integer')
        return

    # cube_list = iris.cube.CubeList(anomaly0)
    # anomaly = anomaly0.concatenate_cube()
    return anomaly


def get_dsi(n_anomaly_cube, anomaly_cube, n):
    # input:
    # n_anomaly: cube with n-monthly anomalies
    # anomaly: cube with monthly anomalies
    # annual_mean_cube: cube with average yearly precipitation over locations
    # n: window length

    # Output: drought severity index cube
    print('Computing DSI...')
    dsi0 = n_anomaly_cube.copy()

    dsi = dsi0.data
    anomaly = anomaly_cube.data[n - 1:]  # avoid first n-1 months for which we cannot compute n-monthly mean.

    #     print(len(anomaly))
    n_anomaly = n_anomaly_cube.data
    #     print(len(dsi))
    #     print(len(n_anomaly))

    dsi[:, :, :] = 0
    for mo in range(n_anomaly_cube.shape[0]):
        dsi_thismonth = dsi[mo, :, :]
        dsi_lastmonth = dsi[mo - 1, :, :]
        anomaly_mo = anomaly[mo, :, :]
        anomaly_mo_1 = anomaly[mo - 1, :, :]
        n_anomaly_mo = n_anomaly[mo, :, :]

        id = np.where((dsi_lastmonth == 0) & (n_anomaly_mo < 0) & (anomaly_mo < 0))
        dsi_thismonth[id] = - anomaly_mo[id]
        id2 = np.where((dsi_lastmonth > 0) & (n_anomaly_mo < 0) & (dsi_lastmonth > anomaly_mo))
        dsi_thismonth[id2] = dsi_lastmonth[id2] - anomaly_mo[id2]
        dsi[mo, :, :] = dsi_thismonth

    dsi0.data = dsi
    return dsi0


def calculate_standardised_dsi(precip_cube, hist_precip_cube, dsi_number):
    period_name = PERIOD_NAMES[dsi_number]

    climatology_years = hist_precip_cube.coord('year').points
    firstyear_past = int(climatology_years[1])
    lastyear_past = int(climatology_years[-1])

    # calculate the average precipitation over historic period
    average_historic_months = hist_precip_cube.aggregated_by('month_number', iris.analysis.MEAN)
    average_historic_months.rename(f' monthly climatological precipitation from {firstyear_past} to {lastyear_past})')

    # find the monthly anomalies from the average climatology
    monthly_anomaly = get_anomaly2(precip_cube, average_historic_months.data, 1, 'monthly')
    monthly_anomaly.rename(f' monthly rainfall anomaly (wrt {firstyear_past}-{lastyear_past})')

    # Calculate the n-month historical average precipitation, where n is dsi_number
    historic_n_months = hist_precip_cube.rolling_window('month_number', iris.analysis.MEAN, dsi_number)
    add_n_month_co(dsi_number, period_name, historic_n_months)  # add n month coordinate to aggregate over later
    add_year_co(dsi_number, historic_n_months)  # add integer year coordinate
    average_historic_n_months = historic_n_months.aggregated_by(f'{period_name}_number', iris.analysis.MEAN)

    # find the anomaly over the n month period
    monthly_n_anomaly = get_anomaly2(precip_cube, average_historic_n_months.data, dsi_number, period_name)

    dsi_cube = get_dsi(monthly_n_anomaly, monthly_anomaly, dsi_number)

    # Standardise the dsi by dividing by the mean annual rainfall at each grid point
    annual_mean = precip_cube.aggregated_by('year', iris.analysis.SUM).collapsed('time', iris.analysis.MEAN)
    standardised_dsi = dsi_cube / annual_mean * 100

    # Redefine metadata for the dsi cube
    standardised_dsi.rename(f' DSI-{dsi_number} (wrt {firstyear_past}-{lastyear_past})')
    standardised_dsi.units = '%'

    # Remove month_number coordinate as it is does not make sense after calculating dsi
    standardised_dsi.remove_coord('month_number')

    return standardised_dsi

# global_precip = iris.load_cube('/scratch/mgrant/UKCP/global/full_time_slice/TS1/monthly/27/pr/27_198012-200011_pr.nc')
# calculate_standardised_dsi(global_precip, global_precip, 3)

def reapply_mask(cube, mask_cube):
    mask = mask_cube.data.mask
    masked_cube = iris.util.mask_cube(cube, mask)

    return masked_cube


def calculate_models_dsi(model, domain, PPE=False):
    if PPE:
        model = str(model).zfill(2)

    hist_precip_filepath = f'{pre_processing.data_files}/uk/{domain}/{model}/TS1/' \
                           f'monthly/pr/{model}_198012-200011_pr.nc'
    hist_precip = iris.load_cube(hist_precip_filepath)
    if 'yyyymm' in [coord.name() for coord in hist_precip.aux_coords]:
        hist_precip.remove_coord('yyyymm')

    for timeslice in pre_processing.SLICES:
        years = pre_processing.SLICES[timeslice]
        firstyear = years[0] - 1
        lastyear = years[-1]
        precip_filepath = f'{pre_processing.data_files}/uk/{domain}/{model}/{timeslice}' \
                          f'/monthly/pr/{model}_{firstyear}12-{lastyear}11_pr.nc'
        precip = iris.load_cube(precip_filepath)
        if 'yyyymm' in [coord.name() for coord in precip.aux_coords]:
            precip.remove_coord('yyyymm')

        for dsi_number in INDEX_NUMBERS:
            print('Calculating DSI for...')
            print('Model: ', model)
            print('Timeslice: ', timeslice)
            print('DSI number: ', dsi_number)

            res = pre_processing.RESOLUTION[domain]
            firstmonth = dsi_number - 1
            dsi = calculate_standardised_dsi(precip, hist_precip, dsi_number)
            # print(dsi)
            # mask = iris.load_cube(f'/scratch/mgrant/UKCP/data/obs/{res}/rainfall/1980-2000/'
            #                       f'obs_{res}_198012-200011_rainfall.nc')
            # mask = mask[firstmonth:]
            # dsi_masked = reapply_mask(dsi, mask)
            # print(dsi_masked)
            dsi_filepath = f'{pre_processing.data_files}/uk/{domain}/{model}' \
                           f'/{timeslice}/monthly/dsi-{dsi_number}/'
            if not os.path.exists(dsi_filepath):
                os.makedirs(dsi_filepath)

            firstyear = int(dsi.coord('year_int').points[0])
            lastyear = int(dsi.coord('year_int').points[-1])
            firstmonth_string = str(firstmonth).zfill(2)
            dsi_filename = f'{model}_{firstyear}{firstmonth_string}-{lastyear}11_dsi-{dsi_number}.nc'
            iris.save(dsi, dsi_filepath + dsi_filename)

# calculate_models_dsi(1, 'local', PPE=True)

def calculate_dsi_full_global_timeseries():
    for model in pre_processing.CMIP_MODELS:
        print(model)
        member = pre_processing.MEMBERS[model]
        full_timeslice_filepath = f'/project/ukcp18/data/land-gcm/uk/60km/rcp85/{member}/pr/mon/' \
                                  f'latest/pr_rcp85_land-gcm_uk_60km_{member}_mon_189912-209911.nc'

        hist_filepath = f'{pre_processing.data_files}/uk/global/{model}/' \
                        f'TS1/monthly/pr/{model}_198012-200011_pr.nc'

        full_timeslice_cube = iris.load_cube(full_timeslice_filepath)
        full_timeslice_cube = full_timeslice_cube[0, ...] # remove ensemble member coord
        UK_mask = iris.load_cube('/project/ukcp/extra/lsm_land-gcm_uk_60km.nc')
        full_timeslice_masked = pre_processing.apply_landsea_mask(full_timeslice_cube, UK_mask, 'global')

        if 'yyyymm' in [coord.name() for coord in full_timeslice_masked.aux_coords]:
            full_timeslice_masked.remove_coord('yyyymm')

        hist_cube = iris.load_cube(hist_filepath)

        for dsi_number in INDEX_NUMBERS:
            print(dsi_number)
            firstmonth = dsi_number - 1
            dsi = calculate_standardised_dsi(full_timeslice_masked, hist_cube, dsi_number)
            dsi_filepath = f'{pre_processing.data_files}/uk/global/{model}' \
                           f'/full_timeslice/monthly/dsi-{dsi_number}/'

            if not os.path.exists(dsi_filepath):
                os.makedirs(dsi_filepath)

            firstmonth_string = str(firstmonth).zfill(2)
            dsi_filename = f'{model}_1900{firstmonth_string}-209911_dsi-{dsi_number}.nc'

            iris.save(dsi, dsi_filepath + dsi_filename)

def load_data(domain, model, index_number):
    firstmonth = str(index_number - 1).zfill(2)
    model_path = f'/data/users/mgrant/ukcp/droughts/data/{domain}/{model}/TS1/monthly/dsi-{index_number}/{model}_1981{firstmonth}-200011_dsi-{index_number}.nc'
    obs_path = f'/data/users/mgrant/ukcp/droughts/data/obs/1km/dsi-{index_number}/1980-2000/obs_1km_198102-200011_dsi-{index_number}.nc'

    return model_path, obs_path

def calculate_obs_dsi():
    for slice in pre_processing.ALT_TIMESLICES:
        firstyear = slice[0]
        lastyear = slice[1]
        for res in pre_processing.OBS_RES:
            precip_filepath = f'{pre_processing.data_files}/uk/obs/{res}/rainfall/' \
                              f'{firstyear}-{lastyear}/obs_{res}_{firstyear}12-{lastyear}11_rainfall.nc'
            precip = iris.load_cube(precip_filepath)
            print(precip)
            for dsi_number in INDEX_NUMBERS:
                print('Calculating DSI for...')
                print(f'Model: obs {res}')
                print('DSI number: ', dsi_number)

                firstmonth = dsi_number - 1
                dsi = calculate_standardised_dsi(precip, precip, dsi_number)
                # mask = iris.load_cube(f'/scratch/mgrant/UKCP/data/obs/{res}/rainfall/'
                #                       f'1980-2000/obs_{res}_198012-200011_rainfall.nc')
                # mask = mask[firstmonth:]
                # dsi_masked = reapply_mask(dsi, mask)
                print(dsi)
                dsi_filepath = f'{pre_processing.data_files}/uk/obs/{res}/dsi-{dsi_number}/{firstyear}-{lastyear}/'

                if not os.path.exists(dsi_filepath):
                    os.makedirs(dsi_filepath)

                # firstyear = int(dsi.coord('year_int').points[0])
                # lastyear = int(dsi.coord('year_int').points[-1])
                firstmonth_string = str(firstmonth).zfill(2)
                dsi_filename = f'obs_{res}_{firstyear + 1}{firstmonth_string}-{lastyear}11_dsi-{dsi_number}.nc'
                iris.save(dsi, dsi_filepath + dsi_filename)


def equalise_ensemble_cube_and_obs_coord_bounds(ensemble_cube, domain, dsi_number):
    firstmonth = (str(dsi_number - 1)).zfill(2)
    # ensemble_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{ensemble}_multimodel_{stat}/'
    #                                     f'{timeslice}/monthly/dsi-{dsi_number}/{ensemble}_multimodel_{stat}_'
    #                                     f'1981{firstmonth}-200011_dsi-{dsi_number}.nc')

    # ensemble_stddev_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{ensemble}_multimodel_{stat}/'
    #                                       f'{timeslice}/monthly/dsi-{dsi_number}/{ensemble}_multimodel_{stat}_'
    #                                       f'1981{firstmonth}-200011_dsi-{dsi_number}.nc')

    obs_res = pre_processing.RESOLUTION[domain]
    obs_cube = iris.load_cube(f'{pre_processing.data_files}/uk/obs/{obs_res}/dsi-{dsi_number}/1980-2000/'
                              f'obs_{obs_res}_1981{firstmonth}-200011_dsi-{dsi_number}.nc')

    cubelist = iris.cube.CubeList((obs_cube, ensemble_cube))

    if 'grid_latitude' in [coord.name() for coord in obs_cube.dim_coords]:
        equalised_bounds_cubelist = equalise_coordinate_bounds(
            equalise_coordinate_bounds(cubelist, 'grid_latitude'), 'grid_longitude')

    elif 'projection_x_coordinate' in [coord.name() for coord in obs_cube.dim_coords]:
        equalised_bounds_cubelist = equalise_coordinate_bounds(
            equalise_coordinate_bounds(cubelist, 'projection_x_coordinate'), 'projection_y_coordinate')

    equalised_bounds_ensemble_cube = equalised_bounds_cubelist[1]
    return equalised_bounds_ensemble_cube


def create_ensemble_cubes_and_ensemble_stat_cubes(country, domain, timeslice, dsi_number, ensemble):
    """
    :param country:
    :param domain:
    :param timeslice:
    :param dsi_number:
    :param ensemble: Takes vales of 'full', 'PPE', or 'CMIP'
    :return:
    """
    model_list = iris.cube.CubeList()
    years = pre_processing.SLICES[timeslice]
    firstyear = years[0]
    lastyear = years[-1]
    firstmonth = str(dsi_number - 1).zfill(2)

    ensemble_list = pre_processing.ENSEMBLE[ensemble]
    # if ensemble == 'PPE':
    #     ensemble_list = pre_processing.PPE_MEMBERS
    # elif ensemble == 'CMIP':
    #     ensemble_list = pre_processing.CMIP_MODELS
    # elif ensemble == 'full':
    #     ensemble_list = pre_processing.PPE_MEMBERS + pre_processing.CMIP_MODELS

    model_one = ensemble_list[0]

    if isinstance(model_one, int):
        model_one = str(model_one).zfill(2)
    print(country)
    filepath_one = f'{pre_processing.data_files}/{country}/{domain}/{model_one}/{timeslice}/monthly/dsi-{dsi_number}/' \
                   f'{model_one}_{firstyear}{firstmonth}-{lastyear}11_dsi-{dsi_number}.nc'
    cube_one = iris.load_cube(filepath_one)

    # Add ensemble member coord
    if 'ensemble_member' not in [coord.name() for coord in cube_one.aux_coords]:
        ensemble_member_coord = iris.coords.AuxCoord(0, long_name='ensemble_member')
        cube_one.add_aux_coord(ensemble_member_coord)
    else:
        cube_one.coord('ensemble_member').points = 0

    model_list.append(cube_one)

    for n, model in enumerate(ensemble_list[1:]):
        if isinstance(model, int):
            model = str(model).zfill(2)

        filepath = f'{pre_processing.data_files}/{country}/{domain}/{model}/{timeslice}/monthly/dsi-{dsi_number}/' \
                   f'{model}_{firstyear}{firstmonth}-{lastyear}11_dsi-{dsi_number}.nc'
        print(filepath)
        cube = iris.load_cube(filepath)
        print('loaded: ', model)
        # if ensemble == 'full' and domain == 'local' and not isinstance(model, int):
        #     cube = cube[:, :-1, :]

        cube_copy = cube_one.copy()
        cube_copy.data = cube.data

        cube_copy.coord('ensemble_member').points = n+1
        # if domain == 'global':
        #     cube_copy.coord('ensemble_member').points = cube.coord('ensemble_member').points
        # elif domain == 'local':
        #     print(cube_copy)
        #     ensemble_member_coord = iris.coords.AuxCoord(n+1, long_name='ensemble_member')
        #     cube_copy.add_aux_coord(ensemble_member_coord)
        #     cube.remove_coord('forecast_reference_time')
        # if domain == 'local':
        #     cube.remove_coord('forecast_reference_time')

        # print(cube.coord('time'))
        # if n > 0:
        #     cube.coord('time').points = model_list[0].coord('time').points
        # print(cube_copy)
        model_list.append(cube_copy)
        # print(cube.coord('time'))
        # print(cube_copy)
    iris.util.equalise_attributes(model_list)
    # iris.util.equalise_cubes(model_list, ['time'])
    print(model_list)
    models_merged = model_list.merge_cube()

    # Save a cube of all the ensmeble members
    ensemble_cube_path = f'{pre_processing.data_files}/{country}/{domain}/ensemble_cubes/' \
                         f'{ensemble}/{timeslice}/monthly/dsi-{dsi_number}/'
    ensemble_cube_name = f'{ensemble}_ensemble_{firstyear}{firstmonth}-{lastyear}11_dsi-{dsi_number}.nc'
    if not os.path.exists(ensemble_cube_path):
        os.makedirs(ensemble_cube_path)

    iris.save(models_merged, ensemble_cube_path + ensemble_cube_name)


def find_ensemble_stats(country, domain, timeslice, dsi_number, ensemble):
    years = pre_processing.SLICES[timeslice]
    firstyear = years[0]
    lastyear = years[-1]
    firstmonth = str(dsi_number - 1).zfill(2)

    filepath_in = f'{pre_processing.data_files}/{country}/{domain}/ensemble_cubes/' \
                  f'{ensemble}/{timeslice}/monthly/dsi-{dsi_number}/'
    filename_in = f'{ensemble}_ensemble_{firstyear}{firstmonth}-{lastyear}11_dsi-{dsi_number}.nc'

    ensemble_cube = iris.load_cube(filepath_in + filename_in)
    # Find means and std-devs of the ensemble cubes
    stats = [
        'mean',
        # 'std-dev',
        # 'median'
    ]

    # master_filepath = f'{pre_processing.data_files}/uk/{domain}/multimodel/{timeslice}/monthly/dsi-{dsi_number}/'

    for stat in stats:
        filepath = f'{pre_processing.data_files}/{country}/{domain}/{ensemble}_multimodel_{stat}/' \
                   f'{timeslice}/monthly/dsi-{dsi_number}/'
        if not os.path.exists(filepath):
            os.makedirs(filepath)

        filename = f'/{ensemble}_multimodel_{stat}_{firstyear}{firstmonth}-{lastyear}11_dsi-{dsi_number}.nc'
        if stat == 'mean':
            cube = ensemble_cube.collapsed('ensemble_member', iris.analysis.MEAN)
        elif stat == 'std-dev':
            cube = ensemble_cube.collapsed('ensemble_member', iris.analysis.STD_DEV)
        elif stat == 'median':
            cube = ensemble_cube.collapsed('ensemble_member', iris.analysis.MEDIAN)

        equalised_bounds_cube = equalise_ensemble_cube_and_obs_coord_bounds(cube, domain, dsi_number)
        iris.save(equalised_bounds_cube, filepath + filename)



# multimodel_stats('global', 'TS1', 3, ensemble='CMIP')

# multimodel_stats('global', 'TS1', 3, ensemble='PPE')

def main():
    # for domain in pre_processing.DOMAIN:
    #     for model in pre_processing.MODELS:
    #         calculate_models_dsi(model, domain)
    # for domain in pre_processing.DOMAIN:
    #     for member in pre_processing.PPE_MEMBERS:
    #         calculate_models_dsi(member, domain, PPE=True)
    # calculate_obs_dsi()
    # calculate_dsi_full_global_timeseries()
    # for number in INDEX_NUMBERS:
    #     for slice in pre_processing.ALT_TIMESLICES:
    #         var = 'dsi-' + str(number)
    #         firstyear = slice[0]
    #         lastyear = slice[1]
    #         pre_processing.create_2p2km_obs(var, firstyear, lastyear)
    # pre_processing.create_2p2km_obs('dsi-3', 1980, 2000)
    for country in pre_processing.COUNTRIES:
        for dom in pre_processing.DOMAIN:
            print(dom)
            for slice in pre_processing.SLICES:
                print(slice)
                for number in INDEX_NUMBERS:
                    print(number)
                    for ensemble in pre_processing.ENSEMBLE:
                        print(ensemble)
                        # create_ensemble_cubes_and_ensemble_stat_cubes(country, dom, slice, number, ensemble)
                        find_ensemble_stats(country, dom, slice, number, ensemble)

if __name__ == '__main__':
    main()

# test_cube = iris.load_cube('/net/spice/scratch/mgrant/UKCP/local/full_time_slice/TS1/monthly/mi-bd773/pr/mi-bd773_198012-200011_pr.nc')
# test_cube = test_cube[:, 320:330, 260:270]
# # qplt.pcolormesh(test_cube[0, ...])
# # plt.show()
# dsi = calculate_standardised_dsi(test_cube, test_cube, 12)
#
# print(dsi)
# print(int(dsi.coord('year_int').points[0]))
# print(dsi.coord('month_number').points[0])
