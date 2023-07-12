import numpy

import DSI
import pre_processing
import iris
import iris.plot as iplt
import matplotlib
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy.ma as ma
import pre_processing
import iris.quickplot as qplt
import os
import numpy as np
from matplotlib import cm
import math
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import warnings
import itertools
import tkinter as tk
import csv


def create_custom_diverging_colormap(color1, color2):
    # define top and bottom colormaps
    color1_r = color1 + '_r'
    top = cm.get_cmap(color1_r, 128)  # r means reversed version
    bottom = cm.get_cmap(color2, 128)  # combine it all
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))
    colorbar = ListedColormap(newcolors, name=color1 + color2)
    return colorbar


def load_model_cube(country, dom, model, timeslice, drought_index, index_number):
    if isinstance(model, int):
        model = str(model).zfill(2)

    firstmonth = (str(index_number - 1)).zfill(2)
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]

    index_cube = iris.load_cube(f'{pre_processing.data_files}/{country}/{dom}/{model}/{timeslice}/'
                                f'monthly/{drought_index}-{index_number}/{model}_'
                                f'{firstyear}{firstmonth}-{lastyear}11_{drought_index}-{index_number}.nc')

    # Mask local dsi cube - ideally this happens in pre-processing but this will do for now
    # if dom == 'local':
    #     cube_with_mask = iris.load_cube(f'{pre_processing.data_files}/uk/obs/2.2km/rainfall/'
    #                                     f'1980-2000/obs_2.2km_198012-200011_rainfall.nc')
    #     lats = index_cube.coord('grid_latitude').points
    #     longs = index_cube.coord('grid_longitude').points
    #     lsm_lat_constraint = iris.Constraint(grid_latitude=lambda cell: lats[0] <= cell <= lats[-1])
    #     lsm_long_constraint = iris.Constraint(grid_longitude=lambda cell: longs[0] <= cell <= longs[-1])
    #     cube_with_mask = cube_with_mask.extract(lsm_lat_constraint & lsm_long_constraint)
    #     firstmonth = index_number - 1
    #     mask = cube_with_mask.data.mask[firstmonth:, ...]
    #     index_cube = iris.util.mask_cube(index_cube, mask)

    return index_cube


# cube = load_model_cube('local', 'MRI-CGCM3', 'TS1', 'dsi', 3)
# qplt.pcolormesh(cube[0,...])
# plt.show()

def find_temporal_mean(dom, model, drought_index):
    index_cubelist = iris.cube.CubeList()
    for index_number in DSI.INDEX_NUMBERS:
        firstmonth = (str(index_number - 1)).zfill(2)
        obs = iris.load_cube(f'{pre_processing.data_files}/uk/obs/1km/{drought_index}-{index_number}/'
                             f'obs_1km_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')
        obs_temporal_mean = obs.collapsed('time', iris.analysis.MEAN)
        index_cubelist.append(obs_temporal_mean)

        for timeslice in pre_processing.SLICES:
            index_cube = load_model_cube(dom, model, timeslice, drought_index, index_number)
            index_temporal_mean = index_cube.collapsed('time', iris.analysis.MEAN)
            index_cubelist.append(index_temporal_mean)

    return index_cubelist


def regrid_obs_and_comparison_to_domain(obs_cube, comparison_cube, domain_cube):
    regrid_scheme = iris.analysis.Nearest(extrapolation_mode='extrapolate')

    regridded_obs = obs_cube.regrid(domain_cube, regrid_scheme)
    regridded_comparison = comparison_cube.regrid(domain_cube, regrid_scheme)
    regridded_cubes = iris.cube.CubeList((regridded_obs, regridded_comparison, domain_cube))

    return regridded_cubes


def calculate_time_averaged_RMSE(obs_cube, model_cube):
    temporally_averaged_obs = obs_cube.collapsed('time', iris.analysis.MEAN)
    temporally_averaged_model = model_cube.collapsed('time', iris.analysis.MEAN)

    RMSE = ((temporally_averaged_obs - temporally_averaged_model) ** 2) ** 0.5

    return RMSE


# def constrain_local_to_global_domain(local_cube, global_cube):
#     global_lats = global_cube.coord('latitude')
#     global_longs = global_cube.coord('longitude')
#
#     min_lat = np.min(global_lats.points)
#     min_long = np.min(global_longs.points)
#     max_lat = np.max(global_lats.points)
#     max_long = np.max(global_longs.points)
#     print(max_long, min_long, max_lat, min_lat)
#     lat_constraint = iris.Constraint(grid_latitude=lambda cell: min_lat <= cell <= max_lat)
#     long_constraint = iris.Constraint(grid_longitude=lambda cell: min_long <= cell <= max_long)
#
#     coord_constraint = lat_constraint & long_constraint
#
#     constrained_local_cube = local_cube.extract(coord_constraint)
#
#     return constrained_local_cube
def equalise_lat_longs(cube1, cube2):
    # if 'latitude' in [coord.name() for coord in (cube1.dim_coords and cube2.dim_coords)] or \
    #         [coord.name() for coord in (cube1.dim_coords and cube2.dim_coords)]:
    #     x_coord = 'longitude'
    #     y_coord = 'latitude'
    # # elif 'projection_y_coordinate' in [coord.name() for coord in (cube1.dim_coords and cube2.dim_coords)]:
    # #     x_coord = 'projection_x_coordinate'
    # #     y_coord = 'projection_y_coordinate'
    # else:
    #     raise(ValueError(f'{cube1.name()} or {cube2.name()} have incompatible x-y coords'))
    # print(cube1.coord('latitude'))
    # print(cube2.coord('latitude'))
    if 'grid_latitude' in [coord.name() for coord in cube1.dim_coords and cube2.dim_coords]:
        latitude = 'grid_latitude'
        longitude = 'grid_longitude'
    elif 'latitude' in [coord.name() for coord in cube1.aux_coords and cube2.aux_coords]:
        latitude = 'latitude'
        longitude = 'longitude'

    lat_diff = cube1.coord(latitude).points - cube2.coord(latitude).points
    long_diff = cube1.coord(longitude).points - cube2.coord(longitude).points
    max_lat_diff = np.max(abs(lat_diff))
    max_long_diff = np.max(abs(long_diff))

    if max_lat_diff > 1 * 10 ** (-5):
        raise ValueError(f'The latitudes between the observation and model cubes are different. '
                         f'Max difference is {max_lat_diff}.')
    else:
        cube2.coord(latitude).points = cube1.coord(latitude).points
        cube2.coord(latitude).bounds = cube1.coord(latitude).bounds

    if max_long_diff > 1 * 10 ** (-5):
        raise ValueError(f'The longitudes between the observation and model cubes are different. '
                         f'Max difference is {max_long_diff}.')
    else:
        cube2.coord(longitude).points = cube1.coord(longitude).points
        cube2.coord(longitude).bounds = cube1.coord(longitude).bounds

    return cube1, cube2


def equalise_cube_boundaries(cube1, cube2):
    if 'projection_x_coordinate' in [coord.name() for coord in cube1.dim_coords and cube2.dim_coords]:
        cube1_x = cube1.coord('projection_x_coordinate').points
        cube1_xmin = cube1_x[0]
        cube1_xmax = cube1_x[-1]

        cube1_y = cube1.coord('projection_y_coordinate').points
        cube1_ymin = cube1_y[0]
        cube1_ymax = cube1_y[-1]

        cube2_x = cube2.coord('projection_x_coordinate').points
        cube2_xmin = cube2_x[0]
        cube2_xmax = cube2_x[-1]

        cube2_y = cube2.coord('projection_y_coordinate').points
        cube2_ymin = cube2_y[0]
        cube2_ymax = cube2_y[-1]

        xmin = min(cube1_xmin, cube2_xmin)
        xmax = max(cube1_xmax, cube2_xmax)

        ymin = min(cube1_ymin, cube2_ymin)
        ymax = max(cube1_ymax, cube2_ymax)

        x_constraint = iris.Constraint(projection_x_coordinate=lambda cell: xmin <= cell <= xmax)
        y_constraint = iris.Constraint(projection_y_coordinate=lambda cell: ymin <= cell <= ymax)

    elif 'grid_latitude' in [coord.name() for coord in cube1.dim_coords and cube2.dim_coords]:
        cube1_x = cube1.coord('grid_longitude').points
        cube1_xmin = cube1_x[0]
        cube1_xmax = cube1_x[-1]

        cube1_y = cube1.coord('grid_latitude').points
        cube1_ymin = cube1_y[0]
        cube1_ymax = cube1_y[-1]

        cube2_x = cube2.coord('grid_longitude').points
        cube2_xmin = cube2_x[0]
        cube2_xmax = cube2_x[-1]

        cube2_y = cube2.coord('grid_latitude').points
        cube2_ymin = cube2_y[0]
        cube2_ymax = cube2_y[-1]

        xmin = min(cube1_xmin, cube2_xmin)
        xmax = max(cube1_xmax, cube2_xmax)

        ymin = min(cube1_ymin, cube2_ymin)
        ymax = max(cube1_ymax, cube2_ymax)

        x_constraint = iris.Constraint(projection_x_coordinate=lambda cell: xmin <= cell <= xmax)
        y_constraint = iris.Constraint(projection_y_coordinate=lambda cell: ymin <= cell <= ymax)

    else:
        raise (AttributeError(f'{cube1.name()} and {cube2.name()} cubes do not have matching coordinate systmes or '
                              f'have incompatible coordinate systems'))

    cube1_constrained = cube1.extract(x_constraint & y_constraint)
    cube2_constrained = cube2.extract(x_constraint & y_constraint)

    # Equalise the bounds of the cubes
    # if 'projection_x_coordinate' in [coord.name() for coord in cube1.dim_coords and cube2.dim_coords]:
    #     cube2_constrained.coord('projection_x_coordinate').bounds = \
    #         cube1_constrained.coord('projection_x_coordinate').bounds
    #     cube2_constrained.coord('projection_y_coordinate').bounds = \
    #         cube1_constrained.coord('projection_y_coordinate').bounds
    #
    # elif 'grid_latitude' in [coord.name() for coord in cube1.dim_coords and cube2.dim_coords]:
    #     cube2_constrained.coord('grid_latitude').bounds = \
    #         cube1_constrained.coord('grid_latitude').bounds
    #     cube2_constrained.coord('grid_longitude').bounds = \
    #         cube1_constrained.coord('grid_longitude').bounds

    return cube1_constrained, cube2_constrained


# obs_cube = iris.load_cube('/data/users/mgrant/ukcp/droughts/data/obs/12km/dsi-3/1980-2000/obs_12km_198102-200011_dsi-3.nc')
# reg_cube = iris.load_cube('/data/users/mgrant/ukcp/droughts/data/regional/MRI-CGCM3/TS1/monthly/dsi-3/MRI-CGCM3_198102-200011_dsi-3.nc')
# print(obs_cube)
# print(reg_cube)
# const_obs, const_reg = equalise_cube_boundaries(obs_cube, reg_cube)
# print(const_obs)
# print(const_reg)
def create_difference_cube(model_cube, obs_cube, stat, percent=False, regrid='None'):
    warnings.filterwarnings("ignore", category=UserWarning)
    # firstmonth = (str(index_number - 1)).zfill(2)
    #
    # if regrid == 'None':
    #     res = pre_processing.RESOLUTION[domain]
    # else:
    #     res = pre_processing.RESOLUTION[regrid]
    #
    # model_cube = load_model_cube(domain, model, 'TS1', drought_index, index_number)
    # obs_cube = iris.load_cube(f'{pre_processing.data_files}/uk/obs/{res}/{drought_index}-{index_number}/1980-2000/'
    #                           f'obs_{res}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

    if stat == 'mean':
        av_model_cube = model_cube.collapsed('time', iris.analysis.MEAN)
        av_obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)
    elif stat == 'std-dev':
        av_model_cube = model_cube.collapsed('time', iris.analysis.STD_DEV)
        av_obs_cube = obs_cube.collapsed('time', iris.analysis.STD_DEV)
    else:
        raise TypeError(f'Stat input must be either mean or std-dev, not {stat}')

    # print(av_model_cube)
    # print(av_obs_cube)
    if not regrid == 'None':
        regrid_scheme = iris.analysis.Nearest(extrapolation_mode='extrapolate')
        regrid_cube = load_model_cube(regrid, 'MRI-CGCM3', 'TS1', 'dsi', 3)
        av_model_cube = av_model_cube.regrid(regrid_cube, regrid_scheme)

    # av_model_cube.remove_coord('month_number') # need to run dsi code again with removing coord

    if 'month_number' in [coord.name() for coord in av_obs_cube.aux_coords]:
        av_obs_cube.remove_coord('month_number')
    # plt.figure()
    # plt.subplot(121)
    # qplt.pcolormesh(av_obs_cube)
    #
    # plt.subplot(122)
    # qplt.pcolormesh(av_model_cube)
    # plt.show()

    # print(av_obs_cube.coord('latitude'))
    # print(av_model_cube.coord('latitude'))
    if ('latitude' or 'grid_latitude') in \
            [coord.name() for coord in av_model_cube.aux_coords or av_model_cube.dim_coords]:
        # if regrid == 'None' and domain == 'global':
        av_model_cube, av_obs_cube = equalise_lat_longs(av_model_cube, av_obs_cube)

    if percent:
        diff_cube = ((av_model_cube - av_obs_cube) / av_obs_cube) * 100
    else:
        diff_cube = av_model_cube - av_obs_cube

    return diff_cube


def save_difference_cube(
        country, model, drought_index, index_number, domain, stat, percent=False, regrid='None'
):
    firstmonth = str(index_number - 1).zfill(2)

    if isinstance(model, int):
        model = str(model).zfill(2)

    if regrid == 'None':
        res = pre_processing.RESOLUTION[domain]
    else:
        res = pre_processing.RESOLUTION[regrid]

    if percent:
        comparison = 'anomaly_from_obs'
    else:
        comparison = 'error'

    model_cube = iris.load_cube(f'{pre_processing.data_files}/{country}/{domain}/{model}/TS1/'
                                f'monthly/{drought_index}-{index_number}/{model}_'
                                f'1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

    obs_cube = iris.load_cube(f'{pre_processing.data_files}/{country}/obs/{res}/{drought_index}-{index_number}/'
                              f'1980-2000/obs_{res}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

    diff_cube = create_difference_cube(model_cube, obs_cube, stat, percent, regrid=regrid)

    difference_filepath = f'{pre_processing.data_files}/{country}/verification/{comparison}/{domain}/{model}' \
                          f'/TS1/monthly/{drought_index}-{index_number}/{stat}/'

    if not regrid == 'None':
        difference_filepath += f'{regrid}_grid/'

    if not os.path.exists(difference_filepath):
        os.makedirs(difference_filepath)

    difference_filename = f'{country}_{model}_1981{firstmonth}-200011_' \
                          f'{drought_index}-{index_number}_{comparison}.nc'
    iris.save(diff_cube, difference_filepath + difference_filename)

    return diff_cube


def create_all_models_difference_plots(
        country, domain, drought_index, index_number, ensemble, stat, regrid='None', with_means=False
):
    firstmonth = str(index_number).zfill(2)
    difference_cubelist = iris.cube.CubeList(())

    if with_means:
        model_list = pre_processing.ENSEMBLE[ensemble] + pre_processing.MEANS
    else:
        model_list = pre_processing.ENSEMBLE[ensemble]

    for model in model_list:
        difference_cube = iris.load_cube(f'{pre_processing.data_files}/{country}/anomaly_from_obs/{model}'
                                         f'/TS1/monthly/{drought_index}_{index_number}/{stat}/{country}_{model}_1981'
                                         f'{firstmonth}-200011_{drought_index}-{index_number}_difference_from_obs.nc')

        difference_cubelist.append(difference_cube)

    num_plots = len(difference_cubelist)
    num_cols = min(num_plots, 4)
    num_rows = math.ceil(num_plots / num_cols)
    fig_height = 4 * num_rows

    if fig_height == 4:
        fig_height += 1
        fig_width = 16
    else:
        fig_width = 12

    cbar_max = 75
    # 12 for PPE
    # 10/12 seems to work with full
    plt.figure(figsize=(fig_width, fig_height), dpi=400)
    # for i in range(i_len):
    #     for j in range(4):
    if stat == 'std-dev':
        stat_name = 'Std Dev'
    else:
        stat_name = stat[0].upper() + stat[1:]

    for i in range(num_rows):
        for j in range(num_cols):
            index = (i * num_cols) + j
            if index < num_plots:
                plt.subplot2grid((num_rows, num_cols), (i, j))

                # plt.subplot2grid((i_len, 4), (i, j))
                if regrid != 'None' and regrid != domain:
                    plt.suptitle(f'Difference in {stat_name} Between {domain[0].upper() + domain[1:]} '
                                 f'Models and Observations for '
                                 f'{drought_index.upper()}-{index_number} on {regrid[0].upper() + regrid[1:]} Grid',
                                 fontsize=20)
                else:
                    plt.suptitle(
                        f'Difference {stat_name} Between {domain[0].upper() + domain[1:]} Models and Observations for '
                        f'{drought_index.upper()}-{index_number}', fontsize=20)
                # try:
                # iplt.pcolormesh(difference_cubelist[(i * 4) + j], vmin=-cbar_max, vmax=cbar_max, cmap='BrBG_r')
                # model = model_list[(i * 4) + j]

                iplt.pcolormesh(difference_cubelist[index], vmin=-cbar_max, vmax=cbar_max, cmap='BrBG_r')
                model = model_list[index]

                plt.subplots_adjust(left=0, hspace=0.3, wspace=0)

                plt.gca().coastlines(resolution='50m')

                plt.title(model, fontsize=15)

                # if j == 3:
                plt.colorbar(label='Difference in DSI (%)', fraction=0.046, pad=0.04, extend='both')
            # except:
            #     e = IndexError
            #     print("e",e)

    model_file = {
        'full': 'all_models',
        'PPE': 'all_PPE_models',
        'CMIP': 'all_CMIP_models'
    }

    if not regrid == 'None' and not regrid == domain:
        filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/{model_file[ensemble]}]/' \
                   f'{domain}_vs_obs_on_{regrid}_grid/{stat}/'
        filename = f'difference_in_{drought_index}-{index_number}_between_' \
                   f'{domain}_and_obs_for_{model_file[ensemble]}_on_{regrid}_grid'
    else:
        filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/{model_file[ensemble]}/{domain}_vs_obs/{stat}/'
        filename = f'difference_in_{drought_index}-{index_number}_between_{domain}_and_obs_for_{model_file[ensemble]}'

    if with_means:
        filename += '_and_multimodle_means'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename)
    plt.close()


def CMIP_ensemble_difference_plots(drought_index, index_number):
    matplotlib.use('Agg')
    cubelist = iris.cube.CubeList(())

    for domain in pre_processing.DOMAIN:
        cube = create_difference_cube('CMIP_ensemble', drought_index, index_number, domain)

        cubelist.append(cube)

    plt.figure(figsize=(9, 12))

    plt.suptitle(f'Percentage Increase From Observations \n to CMIP Ensemble Means for '
                 f'{drought_index.upper()}-{index_number}', fontsize=20)
    cbar_max = 75
    cbar_min = -75
    cmap = 'BrBG_r'
    extend = 'both'
    filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/CMIP_ensemble/{drought_index}-{index_number}/'
    filename = f'difference_in_{drought_index}-{index_number}_between_CMIP_ensemble_and_obs'

    num = len(pre_processing.DOMAIN)
    ensemble_names = ['CMIP', 'PPE', 'Full']
    for i in range(num):
        domain = pre_processing.DOMAIN[i]
        plt.subplot2grid((1, num), (0, i))
        iplt.pcolormesh(cubelist[i], vmax=cbar_max, vmin=cbar_min, cmap=cmap)
        plt.subplots_adjust(left=0, hspace=0.3, wspace=0)
        plt.gca().coastlines(resolution='50m')
        plt.title(f'{domain[0].upper() + domain[1:]}', fontsize=12)
        plt.colorbar(label='Difference in DSI (%)', fraction=0.046, pad=0.04, extend=extend)

    plt.show()

# CMIP_ensemble_difference_plots('dsi', 3)
    # if not os.path.exists(filepath):
    #     os.makedirs(filepath)
    #
    # plt.savefig(filepath + filename)
    # plt.close()

# create_all_models_difference_plots('local', 'dsi', 12, 'CMIP', regrid='None', with_means=True)
def difference_plots_of_ensemble_means_difference_from_obs(drought_index, index_number, highlight_low_values=False):
    matplotlib.use('Agg')
    cubelist = iris.cube.CubeList(())

    for ensemble in pre_processing.MEANS:
        for domain in pre_processing.DOMAIN:
            cube = create_difference_cube(ensemble, drought_index, index_number, domain)

            if highlight_low_values:
                threshold = 10
                cube.data = np.ma.masked_greater(abs(cube.data), threshold)

            cubelist.append(cube)

    plt.figure(figsize=(9, 12))

    if highlight_low_values:
        cbar_max = threshold
        cbar_min = 0
        cmap = 'pink'
        extend = 'max'
        plt.suptitle(f'Percentage Difference Between Observations \n and Ensemble Means for '
                     f'{drought_index.upper()}-{index_number} for Differences < {threshold}%', fontsize=20)
        filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/ensembles/{drought_index}-{index_number}/'
        filename = f'difference_in_{drought_index}-{index_number}_between_ensembles_and_obs_when_difference<{threshold}'
    else:
        plt.suptitle(f'Percentage Increase From Observations \n to Ensemble Means for '
                     f'{drought_index.upper()}-{index_number}', fontsize=20)
        cbar_max = 75
        cbar_min = -75
        cmap = 'BrBG_r'
        extend = 'both'
        filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/ensembles/{drought_index}-{index_number}/'
        filename = f'difference_in_{drought_index}-{index_number}_between_ensembles_and_obs'

    num_rows = len(pre_processing.DOMAIN)
    num_cols = len(pre_processing.MEANS)
    ensemble_names = ['CMIP', 'PPE', 'Full']
    for i in range(num_rows):
        for j in range(num_cols):
            index = (i * num_cols) + j
            ensemble = ensemble_names[i]
            domain = pre_processing.DOMAIN[j]
            plt.subplot2grid((num_rows, num_cols), (i, j))
            iplt.pcolormesh(cubelist[index], vmax=cbar_max, vmin=cbar_min, cmap=cmap)
            plt.subplots_adjust(left=0, hspace=0.3, wspace=0)
            plt.gca().coastlines(resolution='50m')
            plt.title(f'{ensemble} {domain[0].upper() + domain[1:]} Ensemble', fontsize=12)
            plt.colorbar(label='Difference in DSI (%)', fraction=0.046, pad=0.04, extend=extend)

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename)
    plt.close()


def create_difference_plots_of_mean_and_stddev_with_obs_comparison(
        model, drought_index, index_number, domain, stat, regrid='None'
):
    matplotlib.use('Agg')
    warnings.filterwarnings("ignore", category=UserWarning)
    firstmonth = (str(index_number - 1)).zfill(2)

    if isinstance(model, int):
        model = str(model).zfill(2)

    if regrid == 'None':
        res = pre_processing.RESOLUTION[domain]
    else:
        res = pre_processing.RESOLUTION[regrid]

    model_cube = load_model_cube('uk', domain, model, 'TS1', drought_index, index_number)
    obs_cube = iris.load_cube(
        f'{pre_processing.data_files}/uk/obs/{res}/{drought_index}-{index_number}/1980-2000/'
        f'obs_{res}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

    if stat == 'mean':
        collapsed_model_cube = model_cube.collapsed('time', iris.analysis.MEAN)
        collapsed_obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)
    elif stat == 'std-dev':
        collapsed_model_cube = model_cube.collapsed('time', iris.analysis.STD_DEV)
        collapsed_obs_cube = obs_cube.collapsed('time', iris.analysis.STD_DEV)
    else:
        raise TypeError(f'Stat input must be either mean or std-dev, not {stat}')

    diff_cube = iris.load_cube(f'{pre_processing.data_files}/uk/anomaly_from_obs/{model}'
                               f'/TS1/monthly/{drought_index}_{index_number}/{stat}/uk_{model}_1981'
                               f'{firstmonth}-200011_{drought_index}-{index_number}_difference_from_obs.nc')

    mean_difference = np.nanmean(abs(diff_cube.data))

    max_model = np.max(collapsed_model_cube.data)
    max_obs = np.max(collapsed_obs_cube.data)

    min_model = np.min(collapsed_model_cube.data)
    min_obs = np.min(collapsed_obs_cube.data)

    max_dsi = max(max_model, max_obs)
    min_dsi = min(min_model, min_obs)

    # max_difference = np.max(abs(percentage_diff_cube.data))
    # if 'projection_x_coordinate' in [coord.name() for coord in collapsed_obs_cube.dim_coords]:
    #     obs_coords = pre_processing.COORDS['x_y']
    # elif 'grid_latitude' in [coord.name() for coord in collapsed_obs_cube.dim_coords]:
    #     obs_coords = pre_processing.COORDS['lat_long']
    #
    # if 'projection_x_coordinate' in [coord.name() for coord in collapsed_model_cube.dim_coords]:
    #     model_coords = pre_processing.COORDS['x_y']
    # elif 'grid_latitude' in [coord.name() for coord in collapsed_model_cube.dim_coords]:
    #     model_coords = pre_processing.COORDS['lat_long']
    # else:
    #     print('wrong!')

    # mean_model_dsi = collapsed_model_cube.collapsed(model_coords, iris.analysis.MEAN)
    # mean_obs_dsi = collapsed_obs_cube.collapsed(obs_coords, iris.analysis.MEAN)

    plt.figure(figsize=(18, 10))

    if not regrid == 'None' and not regrid == domain:
        plt.suptitle(
            f'{domain[0].upper() + domain[1:]} {model} vs Observations for {drought_index.upper()}-{index_number}'
            f'from 1980 to 2000 on {regrid[0].upper() + regrid[1:]} Grid',
            fontsize=25)
    else:
        plt.suptitle(
            f'{domain[0].upper() + domain[1:]} {model} vs Observations for '
            f'{drought_index.upper()}-{index_number} from 1980 to 2000',
            fontsize=25)

    plt.subplot(131)
    iplt.pcolormesh(collapsed_obs_cube, vmax=max_dsi, vmin=min_dsi, cmap='copper_r')
    plt.title(f'Observations', fontsize=20)
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(fraction=0.046, pad=0.04, label='DSI (%)')
    # plt.figtext(0.22, 0.15, f"average {drought_index}-{index_number} = {mean_obs_dsi.data:.2f}",
    #             ha='center', fontsize=15)

    plt.subplot(132)
    iplt.pcolormesh(collapsed_model_cube, vmax=max_dsi, vmin=min_dsi, cmap='copper_r')
    plt.title(f'{domain[0].upper() + domain[1:]} {model}', fontsize=20)
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(fraction=0.046, pad=0.04, label='DSI (%)')
    # plt.figtext(0.5, 0.15, f"average {drought_index}-{index_number} = {mean_model_dsi.data:.2f}",
    #             ha='center', fontsize=15)

    plt.subplot(133)
    iplt.pcolormesh(diff_cube, vmax=100, vmin=-100, cmap='BrBG_r')
    plt.title(f'Percentage Increase Between \n {model} and Observations', fontsize=20)
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(fraction=0.046, pad=0.04, label='Difference (%)', extend='both')
    plt.figtext(0.78, 0.15, f"absolute mean percentage difference = {mean_difference:.2f} %", ha='center', fontsize=15)

    # plt.show()
    if not regrid == 'None' and not regrid == domain:
        filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/{model}/' \
                   f'{domain}_vs_obs_on_{regrid}_grid/'
        filename = f'difference_in_{drought_index}-{index_number}_between_{domain}_and_obs_for_{model}_on_{regrid}_grid'
    else:
        filepath = f'{pre_processing.plots_files}/spatial_maps/difference_plots/{model}/{domain}_vs_obs/'
        filename = f'difference_in_{drought_index}-{index_number}_between_{domain}_and_obs_for_{model}'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename)
    plt.close()


def table_of_stat_difference(country, stat, comparison):
    summary_dict = {}
    for number in DSI.INDEX_NUMBERS:
        firstmonth = str(number - 1).zfill(2)
        domain_dict = {}
        for domain in pre_processing.DOMAIN:
            model_dict = {}
            for model in pre_processing.CMIP_MODELS + pre_processing.MEANS + pre_processing.MEDIANS:
                diff_cube = iris.load_cube(f'{pre_processing.data_files}/{country}/verification/{comparison}/'
                                           f'{domain}/{model}/TS1/monthly/dsi-{number}/{stat}/{country}_{model}_1981'
                                           f'{firstmonth}-200011_dsi-{number}_{comparison}.nc')

                mean = np.nanmean(abs(diff_cube.data))

                model_dict[model] = mean
            domain_dict[domain] = model_dict
        summary_dict[f'dsi-{number}'] = domain_dict

    if comparison == 'error':
        comparison_table = 'MAE'
    else:
        comparison_table = comparison

    filepath = f'{pre_processing.data_files}/{country}/verification/{comparison_table}/summary_tables/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'{country}_{comparison_table}_summary_table.csv'

    if stat == 'std-dev':
        stat_name = 'Std Dev'
    else:
        stat_name = stat[0].upper() + stat[1:]

    with open(filepath + filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Drought Index', 'Domain', 'Model', comparison_table])
        writer.writeheader()
        index = 1
        for index_key, index_value in summary_dict.items():
            for dom_key, dom_value in index_value.items():
                for model_key, value in dom_value.items():
                    writer.writerow({'Drought Index': index_key,
                                     'Domain': dom_key,
                                     'Model': model_key,
                                     comparison_table: round(value, 2)})
                    index += 1




def create_obs_comparison_plots(drought_index, index_number, regrid='None'):
    warnings.filterwarnings("ignore", category=UserWarning)
    firstmonth = (str(index_number - 1)).zfill(2)

    cubelist = iris.cube.CubeList(())
    max_dsi = 0
    # min_dsi = 1000

    for res in pre_processing.OBS_RES:
        obs_cube = iris.load_cube(f'{pre_processing.data_files}/uk/obs/{res}/{drought_index}-{index_number}/1980-2000/'
                                  f'obs_{res}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

        if res != regrid and regrid != 'None':
            target_grid = iris.load_cube(f'{pre_processing.data_files}/uk/obs/{regrid}/'
                                         f'{drought_index}-{index_number}/1980-2000/'
                                         f'obs_{regrid}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

            regrid_scheme = iris.analysis.Nearest()
            obs_cube = obs_cube.regrid(target_grid, regrid_scheme)

        mean_obs_cube = obs_cube.collapsed('time', iris.analysis.MEAN)

        dim_coords = mean_obs_cube.coords(dim_coords=True)
        coord_names = []

        cubelist.append(mean_obs_cube)

        # for coord in dim_coords:
        #     coord_names.append(coord.name())
        #
        # cube_max = mean_obs_cube.collapsed(coord_names, iris.analysis.MAX).data
        # # cube_min = mean_obs_cube.collapsed(coord_names, iris.analysis.MIN).data
        #
        # if cube_max > max_dsi:
        #     max_dsi = cube_max

        # if cube_min < min_dsi:
        #     min_dsi = cube_min
    max_dsi = {
        3: 11.5,
        6: 16.5,
        12: 22.5
    }

    plt.figure(figsize=(12, 8))
    plt.suptitle(f'{drought_index.upper()}-{index_number} of Observations at all Resolutions', fontsize=24)
    # ticks = [1, 3, 5, 7, 9, 11]

    for i, cube in enumerate(cubelist):
        res = pre_processing.OBS_RES[i]
        plt.subplot2grid((1, 3), (0, i))
        im = iplt.pcolormesh(cube, vmax=max_dsi[index_number], vmin=0, cmap='cubehelix')
        plt.title(f'{res} Observations', fontsize=20)
        plt.gca().coastlines(resolution='50m')
        # cb = plt.colorbar(im, cmap='tab20c_r', ticks=ticks, fraction=0.046, pad=0.04, label='DSI (%)')
        # cb.ax.tick_params()
        plt.colorbar(fraction=0.066, pad=0.04, label='DSI (%)')

    if not regrid == 'None':
        filepath = f'{pre_processing.plots_files}/spatial_maps/obs_comparison/' \
                   f'{regrid}_grid/{drought_index}_{index_number}/'
        filename = f'spatial_maps_of_observed_{drought_index}-{index_number}_regridded_to_{regrid}_grid'
    else:
        filepath = f'{pre_processing.plots_files}/spatial_maps/obs_comparison/' \
                   f'native_grids/{drought_index}_{index_number}/'
        filename = f'spatial_maps_of_observed_{drought_index}-{index_number}'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename)


# create_obs_comparison_plots('dsi', 3, regrid='60km')
# create_difference_plots_with_obs_comparison('full_multimodel_mean', 'dsi', 12, 'local')
# create_difference_plots('MRI-CGCM3', 'dsi', 3, 'regional', regrid='global')

def create_RMSE_plots(model, drought_index, index_number, domain):
    firstmonth = (str(index_number - 1)).zfill(2)

    if domain == 'global':
        comparison = 'local'
    else:
        comparison = 'global'

    domain_cube = load_model_cube(domain, model, 'TS1', drought_index, index_number)
    comparison_cube = load_model_cube(comparison, model, 'TS1', drought_index, index_number)

    # if domain == 'global':
    #     obs_res = ''
    obs_cube = iris.load_cube(f'{pre_processing.data_files}/uk/obs/1km/{drought_index}-{index_number}/'
                              f'obs_1km_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')

    regridded_cubes = regrid_obs_and_comparison_to_domain(obs_cube, comparison_cube, domain_cube)
    obs_vs_comparison_RMSE = calculate_time_averaged_RMSE(regridded_cubes[0], regridded_cubes[1])
    obs_vs_domain_RMSE = calculate_time_averaged_RMSE(regridded_cubes[0], regridded_cubes[2])
    max_comparison_RMSE = np.max(obs_vs_comparison_RMSE.data)
    max_domain_RSME = np.max(obs_vs_domain_RMSE.data)
    max_value = max(max_comparison_RMSE, max_domain_RSME)

    average_comparison_RMSE = np.round(np.mean(obs_vs_comparison_RMSE.data), 3)
    average_domain_RMSE = np.round(np.mean(obs_vs_domain_RMSE.data), 3)

    plt.figure()
    plt.suptitle(f'{drought_index.upper()}-{index_number} RMSE of {model} against Observations')
    plt.subplot(121)
    iplt.pcolormesh(obs_vs_comparison_RMSE, vmax=max_value)
    plt.title(f'{comparison[0].upper() + comparison[1:]} vs Observations \n (mean = {average_comparison_RMSE})')
    plt.colorbar(fraction=0.046, pad=0.04)

    plt.subplot(122)
    iplt.pcolormesh(obs_vs_domain_RMSE, vmax=max_value)
    plt.title(f'{domain[0].upper() + domain[1:]} vs Observations \n (mean = {average_domain_RMSE})')
    plt.colorbar(fraction=0.046, pad=0.04)

    filepath = f'{pre_processing.plots_files}/spatial_maps/RMSE/{model}/{domain}_grid/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'{drought_index}-{index_number}_RMSE_of_local_and_global_against_obs_for_{model}_on_{domain}_grid'

    plt.savefig(filepath + filename)
    plt.close()


def create_spatial_plots(dom, model, drought_index):
    cubelist = find_temporal_mean(dom, model, drought_index)
    plt.figure(figsize=(18, 24))
    plt.suptitle(f'Average Drought Severity Index over 20-year Timeslice for {model} on {dom} domain', fontsize=25)
    index_number = [3, 6, 12]
    years = ['(1980-2000)', '(1980-2000)', '(2020-2040)', '(2060-2080)']

    # cube_max = np.zeros(12)
    # for i in range(len(cubelist)):
    #     cube = cubelist[i]
    #     cube_max[i] = np.max(cube.data)
    #
    # colorbar_max = np.zeros(3)
    # for i in range(3):
    #     colorbar_max[i] = np.max(cube_max[i*4:(i+1)*4])

    for j in range(4):
        for i in range(3):
            if j == 0:
                model_title = 'Observed'
            else:
                model_title = model

            plt.subplot2grid((3, 4), (i, j))
            # plt.axes(projection=ccrs.PlateCarree()).coastlines()
            # iplt.pcolormesh(cubelist[(i*4)+j], cmap='copper_r', vmax=colorbar_max[i])
            iplt.pcolormesh(cubelist[(i * 4) + j], cmap='copper_r')
            plt.gca().coastlines(resolution='50m')
            plt.title(f'{model_title} {drought_index.upper()}-{index_number[i]} \n {years[j]}', fontsize=20)

            # if j == 3:
            plt.colorbar(label='DSI (%)', fraction=0.046, pad=0.04)

    # plt.subplots_adjust(left=0.1,
    #                     bottom=0.1,
    #                     right=0.9,
    #                     top=0.9,
    #                     wspace=0,
    #                     hspace=0.2)

    # plt.show()
    filepath = f'{pre_processing.plots_files}/spatial_maps/{dom}_vs_obs_{drought_index}/{model}/'
    if not os.path.exists(filepath):
        os.makedirs(filepath)
    figname = f'temporally_averaged_{drought_index}_across_{model}_timeslices.png'
    plt.savefig(filepath + figname, dpi=600)
    plt.close()


def find_spatial_mean(dom, model, drought_index, index_number):
    timeseries_list = iris.cube.CubeList()
    for timeslice in pre_processing.SLICES:
        index_cube = load_model_cube(dom, model, timeslice, drought_index, index_number)
        # if dom == 'global':
        #     index_cube.remove_coord('latitude')
        #     index_cube.remove_coord('longitude')
        #     print(index_cube)
        # elif dom == 'local':
        #     grid_areas = iris.analysis.cartography.area_weights(index_cube)

        coords = pre_processing.COORDS[dom]
        index_cube_timeseries = index_cube.collapsed(coords, iris.analysis.MEAN)

        timeseries_list.append(index_cube_timeseries)

    return timeseries_list


def create_cubelist_for_mean_timeseries_comparison(model, time_slice, drought_index, index_number):
    timeseries_list = iris.cube.CubeList()
    # years = pre_processing.SLICES[time_slice]
    # firstyear = years[0]
    # lastyear = years[-1]
    firstmonth = str(index_number - 1).zfill(2)

    for dom in pre_processing.DOMAIN:
        model_cube = load_model_cube(dom, model, time_slice, drought_index, index_number)
        coords = pre_processing.COORDS[dom]
        # grid_areas = iris.analysis.cartography.area_weights(model_cube)
        # print(grid_areas)
        mean_model_cube = model_cube.collapsed(coords, iris.analysis.MEAN)
        timeseries_list.append(mean_model_cube)

    if time_slice == 'TS1':
        obs_cube = iris.load_cube(f'{pre_processing.data_files}/uk/obs/1km/{drought_index}-{index_number}/1980-2000/'
                                  f'obs_1km_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')
        coords = ('latitude', 'longitude')
        # grid_areas = iris.analysis.cartography.area_weights(obs_cube)
        mean_obs_cube = obs_cube.collapsed(coords, iris.analysis.MEAN)
        timeseries_list.append(mean_obs_cube)

    # Equalise the time coord as differences occur due to different calendars
    # for cube in timeseries_list[1:]:
    #     cube.coord('time').points = timeseries_list[0].coord('time').points

    for cube in timeseries_list:
        print(cube)
    time_coord1 = timeseries_list[0].coord('time').points
    time_coord2 = timeseries_list[1].coord('time').points
    print(time_coord1 - time_coord2)

    return timeseries_list


def mean_timeseries_comparison(model, time_slice, drought_index, index_number):
    cubes = create_cubelist_for_mean_timeseries_comparison(model, time_slice, drought_index, index_number)
    # model_cube = load_model_cube('global', model, time_slice, drought_index, index_number)
    # coords = pre_processing.COORDS['global']
    # mean_cube = model_cube.collapsed(coords, iris.analysis.MEAN)
    labels = pre_processing.DOMAIN
    if time_slice == 'TS1':
        labels.append('observations')

    fig, ax = plt.subplots()
    for n, cube in enumerate(cubes):
        iplt.plot(cube, label=labels[n])

    ax.set_title(f'{drought_index.upper()}-{index_number} UK averaged timeseries')
    plt.legend()

    # plt.show()


def create_all_timeslices_timeseries_plots(dom, model, drought_index, index_number):
    timeseries_list = find_spatial_mean(dom, model, drought_index, index_number)
    years = {0: '(1980-2000)', 1: '(2020-2040)', 2: '(2060-2080)'}

    # Plot timeseries
    plt.figure(figsize=(24, 6))
    plt.suptitle(f'{drought_index.upper()} average timeseries across UK for {model}')
    ax1 = plt.subplot2grid((1, 3), (0, 0))
    iplt.plot(timeseries_list[0], 'salmon')

    plt.title(f'{drought_index.upper()}-{index_number} {years[0]}')
    plt.subplot2grid((1, 3), (0, 1), sharey=ax1)
    iplt.plot(timeseries_list[1], 'salmon')
    plt.title(f'{drought_index.upper()}-{index_number} {years[1]}')

    plt.subplot2grid((1, 3), (0, 2), sharey=ax1)
    iplt.plot(timeseries_list[2], 'salmon')
    plt.title(f'{drought_index.upper()}-{index_number} {years[2]}')

    # Save plots
    filepath = f'{pre_processing.plots_files}/timeseries/{drought_index}-{index_number}/{dom}/{model}/'
    if not os.path.exists(filepath):
        os.makedirs(filepath)
    figname = f'spatially_averaged_timeseries_{drought_index}-{index_number}_across_{model}_timeslices.png'
    plt.savefig(filepath + figname, dpi=600)
    plt.close()


def create_all_spatial_plots():
    for drought_index in pre_processing.DROUGHT_INDICES:
        # for dom in pre_processing.DOMAIN:
        #     for model in pre_processing.CMIP_MODELS:
        #         create_spatial_plots(dom, model, drought_index)

        for number in DSI.INDEX_NUMBERS:
            create_obs_comparison_plots(drought_index, number, regrid='None')
            create_obs_comparison_plots(drought_index, number, regrid='60km')


def create_all_timeseries_plots():
    for drought_index in pre_processing.DROUGHT_INDICES:
        for index_number in DSI.INDEX_NUMBERS:
            for dom in pre_processing.DOMAIN:
                for model in pre_processing.CMIP_MODELS:
                    create_all_timeslices_timeseries_plots(dom, model, drought_index, index_number)


def create_all_RMSE_plots():
    for model in pre_processing.CMIP_MODELS:
        print(model)
        for drought_index in pre_processing.DROUGHT_INDICES:
            print(drought_index)
            for index_number in DSI.INDEX_NUMBERS:
                print(index_number)
                for dom in pre_processing.DOMAIN:
                    print(dom)
                    # create_RMSE_plots(model, drought_index, index_number, dom)


import time


def save_all_difference_cubes():
    combinations = itertools.product(
        pre_processing.COUNTRIES,
        pre_processing.MEDIANS,
        pre_processing.DROUGHT_INDICES,
        DSI.INDEX_NUMBERS,
        pre_processing.DOMAIN,
        ['mean']
    )

    for country, model, drought_index, index_number, domain, stat in combinations:
        save_difference_cube(country, model, drought_index, index_number, domain, stat)


def create_all_tables():
    # combinations = itertools.product(
    #     pre_processing.COUNTRIES,
    # )

    for country in pre_processing.COUNTRIES:
        table_of_stat_difference(country, 'mean', 'error')



def create_all_difference_plots():
    for drought_index in pre_processing.DROUGHT_INDICES:
        print(drought_index)
        for index_number in DSI.INDEX_NUMBERS:
            print(index_number)
            difference_plots_of_ensemble_means_difference_from_obs(drought_index, index_number)
            difference_plots_of_ensemble_means_difference_from_obs(drought_index, index_number,
                                                                   highlight_low_values=True)
            # create_difference_plots(model, drought_index, index_number, 'regional')
            for dom in pre_processing.DOMAIN:
                # print(dom)
                for model in pre_processing.CMIP_MODELS + pre_processing.PPE_MEMBERS + pre_processing.MEANS:
                    print(model)
                    print('Non-regrid')
                    # create_difference_plots_with_obs_comparison(model, drought_index, index_number, dom)

                    print('Global regrid')
                    # create_difference_plots_with_obs_comparison(model, drought_index, index_number, dom, regrid='global')

                    # time.sleep(0.1)
                    # # Create the main Tkinter window
                    # root = tk.Tk()
                    # # Schedule the plot function to run in the main thread
                    # root.after(0, create_difference_plots_with_spatial_comparison(model, drought_index, index_number, dom))
                    # # Start the Tkinter event loop
                    # root.mainloop()
                    # root.destroy()

                    # root.after(0, create_difference_plots_with_spatial_comparison(model, drought_index, index_number, dom, regrid='global'))
                    # # Start the Tkinter event loop
                    # root.mainloop()
                    # root.destroy()

                # for ensemble in pre_processing.ENSEMBLE:
                #     # create_all_models_difference_plots(dom, drought_index, index_number, ensemble)
                #     # create_all_models_difference_plots(dom, drought_index, index_number, ensemble, regrid='global')
                #     create_all_models_difference_plots(dom, drought_index, index_number, ensemble, with_means=True)


def main():
    # create_all_spatial_plots()
    # create_all_timeseries_plots()
    # create_all_RMSE_plots()
    save_all_difference_cubes()
    create_all_tables()
    # create_all_difference_plots()
    # create_difference_plots('MRI-CGCM3', 'dsi', 12, 'global')


if __name__ == '__main__':
    main()

# mean_diff_Linear = []
# mean_diff_Nearest = []
# for number in DSI.INDEX_NUMBERS:
#     for dom in pre_processing.DOMAIN:
#         for model in pre_processing.MODELS:
#             print(model)
#             mean_diff_Linear.append(
#                 create_difference_plots(model, 'dsi', number, dom, regrid='global', regrid_scheme='Linear'))
#             mean_diff_Nearest.append(
#                 create_difference_plots(model, 'dsi', number, dom, regrid='global', regrid_scheme='Nearest'))
#         for member in pre_processing.PPE_MEMBERS:
#             print(member)
#             mean_diff_Linear.append(
#                 create_difference_plots(member, 'dsi', number, dom, regrid='global', regrid_scheme='Linear'))
#             mean_diff_Nearest.append(
#                 create_difference_plots(member, 'dsi', number, dom, regrid='global', regrid_scheme='Nearest'))
# mean_diff_Nearest = np.array(mean_diff_Nearest)
# mean_diff_Linear = np.array(mean_diff_Linear)
# diff = mean_diff_Nearest - mean_diff_Linear
# print('Linear = ', mean_diff_Linear)
# print('Nearest = ', mean_diff_Nearest)
# print('Difference = ', diff)
# print('Nearest better = ', np.size(np.where(diff<0)), 'Linear better = ', np.size(np.where(diff>0)))
# print('nearest mean dev = ', np.mean(mean_diff_Nearest), 'linear mean dev = ', np.mean(mean_diff_Linear))
