import pymannkendall as mk
import iris
import numpy as np
import DSI
import iris.plot as iplt
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import pre_processing
import os
import iris.quickplot as qplt

def mask_local_cube(local_cube, index_number):
    cube_with_mask = iris.load_cube('/scratch/mgrant/UKCP/data/obs/2.2km/'
                                    'rainfall/obs_2.2km_198012-200011_rainfall.nc')
    lats = local_cube.coord('grid_latitude').points
    longs = local_cube.coord('grid_longitude').points
    lsm_lat_constraint = iris.Constraint(grid_latitude=lambda cell: lats[0] <= cell <= lats[-1])
    lsm_long_constraint = iris.Constraint(grid_longitude=lambda cell: longs[0] <= cell <= longs[-1])
    cube_with_mask = cube_with_mask.extract(lsm_lat_constraint & lsm_long_constraint)
    firstmonth = index_number - 1
    mask = cube_with_mask.data.mask[0, ...]
    masked_local_cube = iris.util.mask_cube(local_cube, mask)

    return masked_local_cube

def plot_DW_test(domain, model, drought_index, index_number, timeslice):
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]

    firstmonth = str(index_number - 1).zfill(2)

    cube = iris.load_cube(f'{pre_processing.stat_tests_files}/DW_test/{domain}/{model}'
                          f'/{timeslice}/monthly/{drought_index}-{index_number}/{model}_1981{firstmonth}'
                          f'-200011_{drought_index}-{index_number}_DW_test.nc')

    # for i, cube in enumerate(cubelist):
    #     constrained_cube = pre_processing.spatially_constrain_cube(cube)
    #     cubelist[i] = constrained_cube

    plt.figure()
    colors = (["#7BC8F6", "#cacaca", "brown"])
    cmap = ListedColormap(colors)
    bounds = [-1, 0, 1, 2]
    norm = BoundaryNorm(bounds, cmap.N)
    # ticks = [-0.5, 0.5, 1.5]
    # labels = ['Reject NH', 'Inconclusive', 'Do Not Reject NH']

    im = iplt.pcolormesh(cube, cmap=cmap, norm=norm)
    plt.title(f'Durbin-Watson Test for {drought_index}-{index_number} from \n '
              f'{firstyear}-{lastyear} for {domain} {model}')
    plt.gca().coastlines(resolution='50m')
    # cb = plt.colorbar(im, cmap=cmap, norm=norm, ticks=ticks, fraction=0.046, pad=0.04)
    plt.colorbar()
    # cb.set_ticklabels(labels)
    # cb.ax.tick_params()

    filepath = f'{pre_processing.plots_files}/DW_test/obs_vs_all_domains/'
    filename = f'{drought_index}-{index_number}_DW_trendtest_plots_for_{domain}_{model}.png'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename, dpi=300)
    plt.close()

# plot_obs_and_all_domains_MK_trends('global', 'MRI-CGCM3', 'dsi', 3, 'TS1')

def plot_obs_and_all_domains_MK_trends(model, drought_index, index_number):
    firstmonth = str(index_number - 1).zfill(2)
    obs_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/obs/1km/{drought_index}-{index_number}/'
                              f'obs_1981{firstmonth}-200011_{drought_index}-{index_number}_MK_trendtest.nc')

    local_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/local/{model}'
                                 f'/TS1/monthly/{drought_index}-{index_number}/{model}_1981{firstmonth}'
                                 f'-200011_{drought_index}-{index_number}_MK_trendtest.nc')

    regional_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/regional/{model}'
                                   f'/TS1/monthly/{drought_index}-{index_number}/{model}_1981{firstmonth}'
                                   f'-200011_{drought_index}-{index_number}_MK_trendtest.nc')

    global_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}'
                                 f'/TS1/monthly/{drought_index}-{index_number}/{model}_1981{firstmonth}'
                                 f'-200011_{drought_index}-{index_number}_MK_trendtest.nc')

    cubelist = iris.cube.CubeList((obs_cube, local_cube, regional_cube, global_cube))

    for i, cube in enumerate(cubelist):
        constrained_cube = pre_processing.spatially_constrain_cube(cube)
        cubelist[i] = constrained_cube

    plt.figure(figsize=(18, 18))
    plt.subplots_adjust(left=0.01, right=0.8, bottom=0.1, top=0.9, wspace=0.02, hspace=0.1)
    plt.suptitle(f'Mann-Kendall Trend Test for {drought_index.upper()}-{index_number}'
                 f' in Observations and {model} \n Over the Historic Period (1980-2000)',
                 fontsize=25)
    colors = (["#7BC8F6", "#cacaca", "brown"])
    cmap = ListedColormap(colors)
    bounds = [-1, 0, 1, 2]
    norm = BoundaryNorm(bounds, cmap.N)
    ticks = [-0.5, 0.5, 1.5]
    labels = ['decreasing', 'no trend', 'increasing']
    plot_titles = ['Observations', 'Local', 'Regional', 'Global']
    for i in range(2):
        for j in range(2):
            ax = plt.subplot2grid((2, 2), (i, j))
            pos = (i * 2) + j
            im = iplt.pcolormesh(cubelist[pos], cmap=cmap, norm=norm)
            plt.title(plot_titles[pos], fontsize=20)
            plt.gca().coastlines(resolution='50m')
            cb = plt.colorbar(im, cmap=cmap, norm=norm, ticks=ticks, fraction=0.046, pad=0.04)
            cb.set_ticklabels(labels)
            cb.ax.tick_params(labelsize=20)

    filepath = f'{pre_processing.plots_files}/MK_test/obs_vs_all_domains/{model}/{drought_index}-{index_number}/'
    filename = f'{drought_index}-{index_number}_MK_trendtest_plots_for_{model}_and_obs.png'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename, dpi=300)
    plt.close()


def plot_global_MKtrends_for_various_hist_timeslices(model, drought_index, index_number):
    firstmonth = str(index_number - 1).zfill(2)
    TS1_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/'
                              f'{model}/TS1/monthly/{drought_index}-{index_number}/'
                              f'{model}_1981{firstmonth}-200011_{drought_index}-{index_number}_MK_trendtest.nc')

    thirty_year_cube1 = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}/'
                                       f'other_timeslices/1981-2010/monthly/{drought_index}-{index_number}/'
                                       f'{model}_1981{firstmonth}-201011_{drought_index}-{index_number}_MKtrendtest.nc')

    thirty_year_cube2 = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}/'
                                       f'other_timeslices/1971-2000/monthly/{drought_index}-{index_number}/'
                                       f'{model}_1971{firstmonth}-200011_{drought_index}-{index_number}_MKtrendtest.nc')

    forty_year_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}/'
                                     f'other_timeslices/1971-2010/monthly/{drought_index}-{index_number}/'
                                     f'{model}_1971{firstmonth}-201011_{drought_index}-{index_number}_MKtrendtest.nc')

    fifty_year_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}/'
                                     f'other_timeslices/1961-2010/monthly/{drought_index}-{index_number}/'
                                     f'{model}_1961{firstmonth}-201011_{drought_index}-{index_number}_MKtrendtest.nc')

    full_hist_cube = iris.load_cube(f'{pre_processing.stat_tests_files}/MK_trendtest/global/{model}/'
                                    f'other_timeslices/1900-2015/monthly/{drought_index}-{index_number}/'
                                    f'{model}_1900{firstmonth}-201511_{drought_index}-{index_number}_MKtrendtest.nc')

    cubelist = iris.cube.CubeList((TS1_cube, thirty_year_cube1, thirty_year_cube2,
                                   forty_year_cube, fifty_year_cube, full_hist_cube))
    plt.figure(figsize=(36, 18))
    plt.suptitle(f'Mann-Kendall Trend Test for {drought_index.upper()}-{index_number}'
                 f' in Global Model of {model} Over the Various Time Slices',
                 fontsize=25)
    colors = (["#7BC8F6", "#cacaca", "brown"])
    cmap = ListedColormap(colors)
    bounds = [-1, 0, 1, 2]
    norm = BoundaryNorm(bounds, cmap.N)
    ticks = [-0.5, 0.5, 1.5]
    labels = ['decreasing', 'no trend', 'increasing']
    plot_titles = ['1981-2000', '1981-2010', '1971-2000', '1971-2010', '1961-2010', '1900-2015']
    for i in range(2):
        for j in range(3):
            n = (i * 3) + j
            ax = plt.subplot2grid((2, 3), (i, j))
            im = iplt.pcolormesh(cubelist[n], cmap=cmap, norm=norm)
            plt.title(plot_titles[n], fontsize=20)
            plt.gca().coastlines(resolution='50m')
            cb = plt.colorbar(im, cmap=cmap, norm=norm, ticks=ticks, fraction=0.046, pad=0.04)
            cb.set_ticklabels(labels)
            cb.ax.tick_params(labelsize=20)

    filepath = f'{pre_processing.plots_files}/MK_test/global_multiple_timeslices/{model}/{drought_index}-{index_number}/'
    filename = f'{drought_index}-{index_number}_MK_trendtest_plots_for_global_{model}_various_timeslices.png'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename, dpi=600)
    plt.close()
# plot_obs_and_all_domains_MK_trends('ACCESS1-3', 'dsi', 12)

# plot_obs_local_global_MK_trends('MRI-CGCM3', 'dsi', 3)
def main():
    for model in pre_processing.CMIP_MODELS:
        for index in pre_processing.DROUGHT_INDICES:
            # plot_obs_local_global_MK_trends(model, index, 12)
            for index_number in DSI.INDEX_NUMBERS:
                # plot_obs_and_all_domains_MK_trends(model, index, index_number)
                # plot_global_MKtrends_for_various_hist_timeslices(model, index, index_number)
                for dom in pre_processing.DOMAIN:
                    plot_DW_test(dom, model, index, index_number, 'TS1')

if __name__ == '__main__':
    main()
