import iris
import numpy as np
import DSI
import os
import pre_processing
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import math
import itertools

COUNTRIES = [
    'Scotland',
    'England',
    'Wales',
    'Northern_Ireland'
]

MASK_GRID = {
    'global': 'n216_normgrid',
    'regional': 'rcm_12kmrotpol',
    'local': 'cpm_2p2kmrotpol'
}
#
# dsi_gcm = iris.load_cube('/data/users/mgrant/ukcp/droughts/data/regional/MRI-CGCM3/TS1/'
#                          'monthly/dsi-3/MRI-CGCM3_198102-200011_dsi-3.nc')
# Scotland_native_mask = iris.load_cube('/project/ukcp18/shapefilebased_weightcubes/rcm_12kmrotpol_countries_Scotland.nc')
# Scotland_mask_const = pre_processing.spatially_constrain_cube(Scotland_native_mask)
# qplt.pcolormesh(Scotland_native_mask)
# plt.show()
# print(Scotland_native_mask.coord('grid_longitude').points)

# x_constraint = iris.Constraint(projection_x_coordinate=lambda cell: x_min <= cell <= x_max)
# y_constraint = iris.Constraint(projection_y_coordinate=lambda cell: -3.5 <= cell <= y_max)
# regrid_scheme = iris.analysis.Linear()
# Scotland_gcm_mask = Scotland_native_mask.regrid(dsi_gcm, regrid_scheme)
# scotland_gcm_mask = np.where(Scotland_gcm_mask.data == 0, True, False)
# scotland_mask = np.tile(scotland_gcm_mask, (238, 1, 1))
#
# dsi_gcm_scotland = dsi_gcm.copy()
# dsi_mask = dsi_gcm.data.mask
#
# combined_mask = dsi_mask | scotland_mask
#
# dsi_gcm_scotland.data.mask = combined_mask


def make_model_country_cube(country, domain, model, timeslice, drought_index, index_number):
    firstmonth = str(index_number - 1).zfill(2)
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]

    if isinstance(model, int):
        model = str(model).zfill(2)

    if model == 'obs':
        res = pre_processing.DOM_OBS_RES[domain]
        uk_cube_filepath = f'{pre_processing.data_files}/uk/obs/{res}/{drought_index}-{index_number}/1980-2000/' \
                           f'obs_{res}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc'
    elif model == 'obs_2.2km':
        uk_cube_filepath = f'{pre_processing.data_files}/uk/obs/2.2km/{drought_index}-{index_number}/1980-2000/' \
                           f'obs_2.2km_1981{firstmonth}-200011_{drought_index}-{index_number}.nc'
    else:
        uk_cube_filepath = f'{pre_processing.data_files}/uk/{domain}/{model}/' \
                           f'{timeslice}/monthly/{drought_index}-{index_number}/' \
                           f'{model}_{firstyear}{firstmonth}-{lastyear}11_{drought_index}-{index_number}.nc'

    uk_cube = iris.load_cube(uk_cube_filepath)
    country_cube = uk_cube.copy()

    native_country_mask_filepath = f'/project/ukcp18/shapefilebased_weightcubes/{MASK_GRID[domain]}_countries_{country}.nc'
    native_country_mask_cube = iris.load_cube(native_country_mask_filepath)

    # if domain == 'global':
    regrid_scheme = iris.analysis.Linear()
    country_mask_cube = native_country_mask_cube.regrid(uk_cube, regrid_scheme)
    # else:
    #     country_mask_cube = pre_processing.spatially_constrain_cube(country_mask_cube)

    uk_mask = uk_cube.data.mask
    country_mask_at_one_time = np.where(country_mask_cube.data <= 0.5, True, False)

    time_length = len(uk_cube.coord('time').points)
    country_mask = np.tile(country_mask_at_one_time, (time_length, 1, 1))

    combined_mask = uk_mask | country_mask

    country_cube.data.mask = combined_mask
    if model == 'obs':
        filepath_out = f'{pre_processing.data_files}/{country}/obs/{res}/{drought_index}-{index_number}/1980-2000/'
        filename_out = f'1obs_{res}_1981{firstmonth}-200011_{drought_index}-{index_number}.nc'
    elif model == 'obs_2.2km':
        filepath_out = f'{pre_processing.data_files}/{country}/obs/2.2km/{drought_index}-{index_number}/1980-2000/'
        filename_out = f'obs_2.2km_1981{firstmonth}-200011_{drought_index}-{index_number}.nc'
    else:
        filepath_out = f'{pre_processing.data_files}/{country}/{domain}/{model}/' \
                       f'{timeslice}/monthly/{drought_index}-{index_number}/'
        filename_out = f'{model}_{firstyear}{firstmonth}-{lastyear}11_{drought_index}-{index_number}.nc'

    if not os.path.exists(filepath_out):
        os.makedirs(filepath_out)

    iris.save(country_cube, filepath_out + filename_out)


# make_model_country_cube('Northern_Ireland', 'global', 'ACCESS1-3', 'TS1', 'dsi', 6)
# def make_all_country_cubes():
#     for country in COUNTRIES:
#         for domain in pre_processing.DOMAIN:
#             for model in pre_processing.CMIP_MODELS + pre_processing.PPE_MEMBERS + pre_processing.MEANS:
#                 for timeslice in pre_processing.SLICES:
#                     for index in pre_processing.DROUGHT_INDICES:
#                         for number in DSI.INDEX_NUMBERS:
#                             make_model_country_cube(country, domain, model, timeslice, index, number)

def make_all_country_cubes():
    combinations = itertools.product(
         COUNTRIES,
         pre_processing.DOMAIN,
         pre_processing.MODELS + pre_processing.MEANS + pre_processing.MEDIANS,
         pre_processing.SLICES,
         pre_processing.DROUGHT_INDICES,
         DSI.INDEX_NUMBERS
    )

    for country, domain, model, timeslice, index, number in combinations:
        make_model_country_cube(country, domain, model, timeslice, index, number)

def main():
    make_all_country_cubes()


if __name__ == '__main__':
    main()

