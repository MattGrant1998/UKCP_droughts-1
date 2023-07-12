import pymannkendall as mk
import iris
import numpy as np
import DSI
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import pre_processing
import os
import statsmodels.api as sm
from statsmodels.stats.stattools import durbin_watson as DW


def perform_DW_test(cube):
    dL = 1.75
    dU = 1.79
    data = cube.data
    DW_array = np.zeros_like(data[0, ...])
    for i in range(data.shape[1]):
        for j in range(data.shape[2]):
            if DW_array.mask[i, j] == True:
                continue
            DSI = data[:, i, j]
            time = cube.coord('time').points
            time = sm.add_constant(time)
            reg = sm.OLS(DSI, time).fit()
            dw_result = DW(resids=np.array(reg.resid))

            if dw_result < dL:
                # Reject null hypothesis - serial correlation could be an issue
                DW_array[i, j] = -1
            elif (dw_result >= dL) and (dw_result <= dU):
                # Results are inconclusive
                DW_array[i, j] = 0
            elif dw_result > dU:
                # Do not reject hypothesis - serial correlation not a problem
                DW_array[i, j] = 1
            else:
                print('mk test returned unexpected result')

    DW_cube = cube[0, ...].copy()
    DW_cube.data = DW_array

    return DW_cube

def compute_MK_trendtest_cube(drought_index_cube, HR_mod=True):
    drought_index_data = drought_index_cube.data
    trend_array = np.zeros_like(drought_index_data[0, ...])
    for i in range(drought_index_data.shape[1]):
        for j in range(drought_index_data.shape[2]):
            if trend_array.mask[i, j] == True:
                continue

            if HR_mod:
                mk_result = mk.hamed_rao_modification_test(drought_index_data[:, i, j])
            else:
                mk_result = mk.original_test(drought_index_data[:, i, j])

            if mk_result[0] == 'increasing':
                trend_array[i, j] = 1
            elif mk_result[0] == 'decreasing':
                trend_array[i, j] = -1
            elif mk_result[0] == 'no trend':
                trend_array[i, j] = 0
            else:
                print('mk test returned unexpected result')

    trend_cube = drought_index_cube[0,...].copy()
    trend_cube.data = trend_array

    return trend_cube

drought_index_filepath = f'{pre_processing.data_files}/global/MPI-ESM-LR/' \
                             f'TS1/monthly/dsi-3/'

drought_index_filename = f'MPI-ESM-LR_198102-' \
                             f'200011_dsi-3.nc'
cube = iris.load_cube(drought_index_filepath + drought_index_filename)

# dw_cube = perform_DW_test(cube)
# qplt.pcolormesh(dw_cube)
# plt.show()

trendtest_array = np.zeros_like(cube.data[0,...])
for i in range(cube.shape[1]):
    for j in range(cube.shape[2]):
        mk_result = mk.original_test(cube.data[:, i, j])
        print(mk_result)

def save_DW_test_cube(index, index_number, dom, model, timeslice):
    firstmonth = str(index_number - 1).zfill(2)
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]

    drought_index_filepath = f'{pre_processing.data_files}/{dom}/{model}/' \
                             f'{timeslice}/monthly/{index}-{index_number}/'

    drought_index_filename = f'{model}_{firstyear}{firstmonth}-' \
                             f'{lastyear}11_{index}-{index_number}.nc'

    DWtest_filepath = f'{pre_processing.stat_tests_files}/DW_test/{dom}/{model}/' \
                         f'{timeslice}/monthly/{index}-{index_number}/'

    DWtest_filename = f'{model}_{firstyear}{firstmonth}-' \
                         f'{lastyear}11_{index}-{index_number}_DW_test.nc'

    if not os.path.exists(DWtest_filepath):
        os.makedirs(DWtest_filepath)

    drought_index_cube = iris.load_cube(drought_index_filepath + drought_index_filename)
    DWtest_cube = perform_DW_test(drought_index_cube)
    iris.save(DWtest_cube, DWtest_filepath + DWtest_filename)


def save_obs_DW_test_cube(index, index_number, res):
    firstmonth = str(index_number - 1).zfill(2)
    drought_index_filepath = f'{pre_processing.data_files}/obs/{res}/{index}-{index_number}/'
    drought_index_filename = f'obs_{res}_1981{firstmonth}-200011_{index}-{index_number}.nc'
    DWtest_filepath = f'{pre_processing.stat_tests_files}/DW_test/obs/{res}/{index}-{index_number}/'
    DWtest_filename = f'obs_1981{firstmonth}-200011_{index}-{index_number}_DW_test.nc'

    if not os.path.exists(DWtest_filepath):
        os.makedirs(DWtest_filepath)

    drought_index_cube = iris.load_cube(drought_index_filepath + drought_index_filename)
    DWtest_cube = perform_DW_test(drought_index_cube)
    iris.save(DWtest_cube, DWtest_filepath + DWtest_filename)

    # else:
    #     print('Skipping: ')
    #     print('Observations')
    #     print(index, '-', index_number)

def save_obs_MK_trendtest_cube(index, index_number, res, HR_mod=True):
    firstmonth = str(index_number - 1).zfill(2)
    drought_index_filepath = f'{pre_processing.data_files}/obs/1km/{index}-{index_number}/'
    drought_index_filename = f'obs_{res}_1981{firstmonth}-200011_{index}-{index_number}.nc'

    if HR_mod:
        trendtest_filepath = f'{pre_processing.stat_tests_files}/MK_trendtest_HRmod/' \
                             f'obs/{res}/{index}-{index_number}/1980-2000/'
        trendtest_filename = f'obs_1981{firstmonth}-200011_{index}-{index_number}_MK_trendtest.nc'
    else:
        trendtest_filepath = f'{pre_processing.stat_tests_files}/MK_trendtest/obs/{res}/{index}-{index_number}/1980-2000/'
        trendtest_filename = f'obs_1981{firstmonth}-200011_{index}-{index_number}_MK_trendtest.nc'

    if not os.path.exists(trendtest_filepath):
        os.makedirs(trendtest_filepath)

    drought_index_cube = iris.load_cube(drought_index_filepath + drought_index_filename)
    trendtest_cube = compute_MK_trendtest_cube(drought_index_cube, HR_mod)
    iris.save(trendtest_cube, trendtest_filepath + trendtest_filename)

    # else:
    #     print('Skipping: ')
    #     print('Observations')
    #     print(index, '-', index_number)


def save_models_MK_trendtest_cube(index, index_number, dom, model, timeslice, HR_mod=True):
    firstmonth = str(index_number - 1).zfill(2)
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]

    drought_index_filepath = f'{pre_processing.data_files}/{dom}/{model}/' \
                             f'{timeslice}/monthly/{index}-{index_number}/'

    drought_index_filename = f'{model}_{firstyear}{firstmonth}-' \
                             f'{lastyear}11_{index}-{index_number}.nc'

    if HR_mod:
        trendtest_filepath = f'{pre_processing.stat_tests_files}/MK_trendtest_HRmod/{dom}/{model}/' \
                             f'{timeslice}/monthly/{index}-{index_number}/'

        trendtest_filename = f'{model}_{firstyear}{firstmonth}-' \
                             f'{lastyear}11_{index}-{index_number}_MK_trendtest.nc'
    else:
        trendtest_filepath = f'{pre_processing.stat_tests_files}/MK_trendtest/{dom}/{model}/' \
                             f'{timeslice}/monthly/{index}-{index_number}/'

        trendtest_filename = f'{model}_{firstyear}{firstmonth}-' \
                             f'{lastyear}11_{index}-{index_number}_MK_trendtest.nc'

    if not os.path.exists(trendtest_filepath):
        os.makedirs(trendtest_filepath)

    drought_index_cube = iris.load_cube(drought_index_filepath + drought_index_filename)
    trendtest_cube = compute_MK_trendtest_cube(drought_index_cube, HR_mod)
    iris.save(trendtest_cube, trendtest_filepath + trendtest_filename)

    # else:
    #     print('Skipping: ')
    #     print(index, '-', index_number)
    #     print(dom)
    #     print(model)
    #     print(timeslice)


def MKtrendtest_selected_length_global(index, index_number, model, firstyear, lastyear):
    firstmonth = str(index_number - 1).zfill(2)
    drought_index_filepath = f'{pre_processing.data_files}/global/{model}/full_timeslice/monthly/' \
                             f'{index}-{index_number}/{model}_1900{firstmonth}-209911_{index}-{index_number}.nc'
    drought_index_cube = iris.load_cube(drought_index_filepath)

    # Slice the drought index cube to be for the years of relevance only
    year_constraint = iris.Constraint(year=lambda cell: firstyear - 1 <= cell <= lastyear)
    drought_index_yearsliced = drought_index_cube.extract(year_constraint)
    drought_index_timesliced = drought_index_yearsliced[11:-1, ...]

    trendtest_filepath = f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}/' \
                         f'other_timeslices/{firstyear}-{lastyear}/monthly/{index}-{index_number}/'

    trendtest_filename = f'{model}_{firstyear}{firstmonth}-{lastyear}11_{index}-{index_number}_MKtrendtest.nc'

    if not os.path.exists(trendtest_filepath):
        os.makedirs(trendtest_filepath)


    trendtest_cube = compute_MK_trendtest_cube(drought_index_timesliced)
    trendtest_cube.rename(f'DSI MK trend test from {firstyear}-{lastyear}')
    iris.save(trendtest_cube, trendtest_filepath + trendtest_filename)


def main():
    print('Finding MK trend test for...')
    # save_obs_MK_trendtest_cube('dsi', 12)
    # save_models_MK_trendtest_cube('dsi', 6, 'local', 'ACCESS1-3', 'TS1')
    for index in pre_processing.DROUGHT_INDICES:
        for index_number in DSI.INDEX_NUMBERS:
            print('Drought index: ', index, '-', index_number)
            print('Domain: Obs')
            # save_obs_MK_trendtest_cube(index, index_number)
            for model in pre_processing.CMIP_MODELS:
                print('Model: ', model)
                # MKtrendtest_selected_length_global(index, index_number, model, 1981, 2099)
                # MKtrendtest_selected_length_global(index, index_number, model, 1901, 2000)
                # MKtrendtest_selected_length_global(index, index_number, model, 1981, 2010)
                # MKtrendtest_selected_length_global(index, index_number, model, 1971, 2000)
                # MKtrendtest_selected_length_global(index, index_number, model, 1961, 2010)
                # MKtrendtest_selected_length_global(index, index_number, model, 1900, 2015)
                for dom in pre_processing.DOMAIN:
                    print('Domain: ', dom)
                    save_models_MK_trendtest_cube(index, index_number, dom, model, 'TS1')
                    # save_DW_test_cube(index, index_number, dom, model, 'TS1')
                for timeslice in pre_processing.SLICES:
                    print('Time Slice: ', timeslice)
                    # save_models_MK_trendtest_cube(index, index_number, dom, model, timeslice)
                    # save_DW_test_cube(index, index_number, 'global', model, timeslice)
                    # save_DW_test_cube(index, index_number, 'local', model, timeslice)

            # for res in pre_processing.OBS_RES:
            #     save_obs_DW_test_cube(index, index_number, res)
            #     save_obs_MK_trendtest_cube(index, index_number, res)


            # save_models_MK_trendtest_cube('dsi', 3, 'global', 'MRI-CGCM3', 'TS1')
# save_obs_MK_trendtest_cube('dsi', 6)
# if __name__ == '__main__':
#     main()
