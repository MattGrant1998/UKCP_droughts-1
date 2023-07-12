import iris
import numpy as np
import DSI
import os
import pre_processing
import iris.plot as iplt
import matplotlib.pyplot as plt
import math
import itertools

def create_anomaly_cube(domain, model, drought_index, index_number, timeslice1, timeslice2):
    """
    Creates a cube of the percentage increase in DSI-n from the timeslice1 mean to the timeslice2 mean.
    :param domain:
    :param model:
    :param drought_index:
    :param index_number:
    :param timeslice1:
    :param timeslice2:
    :return:
    """
    if isinstance(model, int):
        model = str(model).zfill(2)

    firstmonth = (str(index_number - 1)).zfill(2)
    timeslice1_firstyear = pre_processing.SLICES[timeslice1][0]
    timeslice1_lastyear = pre_processing.SLICES[timeslice1][-1]

    timeslice2_firstyear = pre_processing.SLICES[timeslice2][0]
    timeslice2_lastyear = pre_processing.SLICES[timeslice2][-1]

    cube1 = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/{timeslice1}/monthly/'
                           f'{drought_index}-{index_number}/{model}_{timeslice1_firstyear}{firstmonth}-'
                           f'{timeslice1_lastyear}11_{drought_index}-{index_number}.nc')

    cube2 = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/{timeslice2}/monthly/'
                           f'{drought_index}-{index_number}/{model}_{timeslice2_firstyear}{firstmonth}-'
                           f'{timeslice2_lastyear}11_{drought_index}-{index_number}.nc')

    cube1_mean = cube1.collapsed('time', iris.analysis.MEAN)
    cube2_mean = cube2.collapsed('time', iris.analysis.MEAN)

    anomaly_cube = ((cube2_mean - cube1_mean) / cube1_mean) * 100

    filepath = f'{pre_processing.data_files}/uk/{domain}/{model}/' \
               f'{timeslice2}_minus_{timeslice1}/monthly/{drought_index}-{index_number}/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'{model}_{timeslice2}_minus_{timeslice1}_{drought_index}-{index_number}.nc'

    iris.save(anomaly_cube, filepath + filename)

def find_difference_between_CMIP_and_PPE_means(domain, drought_index, index_number):
    CMIP_anomaly_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/CMIP_multimodel_mean/TS3_minus_TS1'
                                       f'/monthly/{drought_index}-{index_number}/'
                                       f'CMIP_multimodel_mean_TS3_minus_TS1_{drought_index}-{index_number}.nc')

    PPE_anomaly_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/PPE_multimodel_mean/TS3_minus_TS1'
                                      f'/monthly/{drought_index}-{index_number}/'
                                      f'PPE_multimodel_mean_TS3_minus_TS1_{drought_index}-{index_number}.nc')

    difference_cube = CMIP_anomaly_cube - PPE_anomaly_cube

    difference_filepath = f'{pre_processing.data_files}/uk/{domain}/CMIP_PPE_anomaly_difference/TS3_minus_TS1' \
                          f'/monthly/{drought_index}-{index_number}/'

    if not os.path.exists(difference_filepath):
        os.makedirs(difference_filepath)

    difference_filename = f'CMIP_PPE_anomaly_difference_TS3_minus_TS1_{drought_index}-{index_number}.nc'

    iris.save(difference_cube, difference_filepath + difference_filename)


def plot_individual_TS3_TS1_anomaly_cubes(domain, model, drought_index, index_number):
    if isinstance(model, int):
        model = str(model).zfill(2)
    anomaly_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/'
                                  f'TS3_minus_TS1/monthly/{drought_index}-{index_number}/'
                                  f'{model}_TS3_minus_TS1_{drought_index}-{index_number}.nc')

    # max_anomaly = np.max(anomaly_cube.data)
    # min_anomaly = abs(np.min(anomaly_cube.data))
    # cbar_max = max(min_anomaly, max_anomaly)
    # if domain == 'global':
    #     cbar_max = 3
    # else:
    #     cbar_max = 4
    plt.figure(figsize=(6, 6), dpi=200)
    iplt.pcolormesh(anomaly_cube, vmin=-100, vmax=100, cmap='BrBG_r')
    plt.title(f'{domain[0].upper() + domain[1:]} {model} Increase in \n {drought_index.upper()}-{index_number} '
              f'from the 1980-2000 Mean to the 2060-2080 Mean')
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(label='Increase in DSI (%)', extend='both')
    plt.gca().set_aspect('equal')
    plt.subplots_adjust(left=0, right=0.9, bottom=0.1, top=0.9)

    filepath = f'{pre_processing.plots_files}/spatial_maps/anomaly_plots/individual_plots/{domain}/' \
               f'{model}/TS3_minus_TS1/{drought_index}-{index_number}/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'percentage_anomaly_in_{drought_index}-{index_number}_' \
               f'from_1980-2000_to_2060-2080_for_{domain}_{model}.png'

    plt.savefig(filepath + filename)
    plt.close()


def plot_domain_comparison_anomalies(model, index_number):
    if isinstance(model, int):
        model = str(model).zfill(2)

    cubelist = []
    for domain in pre_processing.DOMAIN:
        anomaly_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/'
                                      f'TS3_minus_TS1/monthly/dsi-{index_number}/'
                                      f'{model}_TS3_minus_TS1_dsi-{index_number}.nc')
        cubelist.append(anomaly_cube)
    plt.figure(figsize=(14, 8))
    domain_names = ['GCM', 'RCM', 'CPM']

    if model == 'CMIP_multimodel_mean':
        model_name = 'CMIP Ensemble'
    elif model == 'PPE_multimodel_mean':
        model_name = 'PPE Ensemble'
    elif model == 'full_multimodel_mean':
        model_name = 'CMIP + PPE Ensemble'
    else:
        model_name = model

    if model == 'CMIP_PPE_anomaly_difference':
        plt.suptitle(f'Difference Between the \n CMIP and PPE Anomalies for \n DSI-{index_number}', fontsize=38)
    else:
        plt.suptitle(f'Anomaly Between 1980-2000 and 2060-2080 \n for the {model_name} \n DSI-{index_number}', fontsize=38)
    for i in range(3):
        plt.subplot2grid((1, 3), (0, i))

        cbar_max = 100
        iplt.pcolormesh(cubelist[i], vmin=-cbar_max, vmax=cbar_max, cmap='BrBG_r')

        plt.gca().coastlines(resolution='50m')

        plt.title(domain_names[i], fontsize=30)
        plt.tight_layout()
        plt.subplots_adjust(wspace=0.5, right=0.9)
        cb = plt.colorbar(label='Difference in DSI (%)', fraction=0.066, pad=0.04, extend='both')
        cb.set_label(label='Difference in DSI (%)', fontsize=23)
        cb.ax.tick_params(labelsize=18)
        # set_ticklabels(cb.get_ticklabels(), fontsize=15)

    filepath = f'{pre_processing.plots_files}/spatial_maps/anomaly_plots/model/all_domains/' \
               f'{model}/TS3_minus_TS1/dsi-{index_number}/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'percentage_anomaly_in_dsi-{index_number}_' \
               f'from_1980-2000_to_2060-2080_for_all_domains_{model}.png'

    plt.savefig(filepath + filename)
    plt.close()

# plot_domain_comparison_anomalies('CMIP_multimodel_mean', 3)
# plot_individual_TS3_TS1_anomaly_cubes('global', 'PPE_multimodel_mean', 'dsi', 3)
def plot_all_models_anomaly_cubes(domain, drought_index, index_number, ensemble, with_means=False):
    anomaly_cubelist = iris.cube.CubeList(())

    if with_means:
        model_list = pre_processing.ENSEMBLE[ensemble] + pre_processing.MEANS
    else:
        model_list = pre_processing.ENSEMBLE[ensemble]

    for model in model_list:
        model_filepath = f'{pre_processing.data_files}/uk/{domain}/{model}/' \
                   f'TS3_minus_TS1/monthly/{drought_index}-{index_number}/'
        model_filename = f'{model}_TS3_minus_TS1_{drought_index}-{index_number}.nc'
        anomaly_cube = iris.load_cube(model_filepath + model_filename)
        # anomaly_cube = create_anomaly_cube(domain, model, drought_index, index_number, 'TS1', 'TS3')
        anomaly_cubelist.append(anomaly_cube)

    num_plots = len(anomaly_cubelist)
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
    for i in range(num_rows):
        for j in range(num_cols):
            index = (i * num_cols) + j
            if index < num_plots:
                plt.subplot2grid((num_rows, num_cols), (i, j))
                plt.suptitle(f'Difference Between {domain[0].upper() + domain[1:]} Models '
                             f'{drought_index.upper()}-{index_number} Between 1980-2000 Mean and 2060-2080 Mean',
                             fontsize=20)

                iplt.pcolormesh(anomaly_cubelist[index], vmin=-cbar_max, vmax=cbar_max, cmap='BrBG_r')
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

    filepath = f'{pre_processing.plots_files}/spatial_maps/anomaly_plots/{model_file[ensemble]}/{domain}/' \
               f'{model}/TS3_minus_TS1/{drought_index}-{index_number}/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'percentage_anomaly_in_{drought_index}-{index_number}_' \
               f'from_1980-2000_to_2060-2080_for_{domain}_{model_file[ensemble]}.png'
    if with_means:
        filename += '_and_multimodle_means'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename)
    plt.close()


def create_TS3_TS1_anomaly_plots_with_comparison(domain, model, drought_index, index_number, regrid='None'):
    firstmonth = (str(index_number - 1)).zfill(2)

    if isinstance(model, int):
        model = str(model).zfill(2)

    TS1_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/TS1/monthly/'
                              f'{drought_index}-{index_number}/{model}_1981{firstmonth}-'
                              f'200011_{drought_index}-{index_number}.nc')

    TS3_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/TS3/monthly/'
                              f'{drought_index}-{index_number}/{model}_2061{firstmonth}-'
                              f'208011_{drought_index}-{index_number}.nc')

    anomaly_cube = iris.load_cube(f'{pre_processing.data_files}/uk/{domain}/{model}/'
                                  f'TS3_minus_TS1/monthly/{drought_index}-{index_number}/'
                                  f'{model}_TS3_minus_TS1_{drought_index}-{index_number}.nc')

    av_TS1_cube = TS1_cube.collapsed('time', iris.analysis.MEAN)
    av_TS3_cube = TS3_cube.collapsed('time', iris.analysis.MEAN)

    mean_difference = np.nanmean(anomaly_cube.data)

    max_dsi = max(np.nanmax(av_TS1_cube.data), np.nanmax(av_TS3_cube.data))
    min_dsi = min(np.nanmin(av_TS1_cube.data), np.nanmin(av_TS3_cube.data))

    plt.figure(figsize=(18, 10))

    if not regrid == 'None' and not regrid == domain:
        plt.suptitle(
            f'Anomaly comparison for the 1980-2000 and 2060-2080 Average {drought_index.upper()}-{index_number} '
            f'for {domain[0].upper() + domain[1:]} {model} on {regrid[0].upper() + regrid[1:]} Grid',
            fontsize=25)
    else:
        plt.suptitle(
            f'Anomaly comparison for the 1980-2000 and 2060-2080 Average {drought_index.upper()}-{index_number} '
            f'for {domain[0].upper() + domain[1:]} {model}',
            fontsize=25)

    plt.subplot(131)
    iplt.pcolormesh(av_TS1_cube, vmax=max_dsi, vmin=min_dsi, cmap='copper_r')
    plt.title(f'1980-2000', fontsize=20)
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(fraction=0.046, pad=0.04, label='DSI (%)')
    # plt.figtext(0.22, 0.15, f"average {drought_index}-{index_number} = {mean_obs_dsi.data:.2f}",
    #             ha='center', fontsize=15)

    plt.subplot(132)
    iplt.pcolormesh(av_TS3_cube, vmax=max_dsi, vmin=min_dsi, cmap='copper_r')
    plt.title(f'2060-2080', fontsize=20)
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(fraction=0.046, pad=0.04, label='DSI (%)')
    # plt.figtext(0.5, 0.15, f"average {drought_index}-{index_number} = {mean_model_dsi.data:.2f}",
    #             ha='center', fontsize=15)

    plt.subplot(133)
    iplt.pcolormesh(anomaly_cube, vmax=100, vmin=-100, cmap='BrBG_r')
    plt.title(f'Anomaly', fontsize=20)
    plt.gca().coastlines(resolution='50m')
    plt.colorbar(fraction=0.046, pad=0.04, label='Difference (%)', extend='both')
    plt.figtext(0.78, 0.15, f"average anomaly = {mean_difference:.2f} %", ha='center', fontsize=15)

    # plt.show()
    if not regrid == 'None' and not regrid == domain:
        filepath = f'{pre_processing.plots_files}/spatial_maps/anomaly_plots/comparison/{domain}/{model}/' \
                   f'{drought_index}-{index_number}/TS3_minus_TS1/{regrid}_grid/'
        filename = f'TS3_minus_TS1_anomalies_for_{drought_index}-{index_number}_{domain}_{model}_on_{regrid}_grid'
    else:
        filepath = f'{pre_processing.plots_files}/spatial_maps/anomaly_plots/comparison/{domain}/{model}/' \
                   f'{drought_index}-{index_number}/TS3_minus_TS1/native_grid/'
        filename = f'TS3_minus_TS1_anomalies_for_{drought_index}-{index_number}_{domain}_{model}'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    plt.savefig(filepath + filename)
    plt.close()

# create_TS3_TS1_anomaly_plots_with_comparison('MRI-CGCM3', 'dsi', 3, 'global', regrid='None')

def create_TS3_minus_TS1_anomalies():
    print('------------------- Calculating TS3 - TS1 anomalies -----------------------')
    for domain in pre_processing.DOMAIN:
        print('Domain: ', domain)
        for model in pre_processing.MEANS:
            print('Model: ', model)
            for index_number in DSI.INDEX_NUMBERS:
                print('Index: ', 'dsi-', index_number)
                find_difference_between_CMIP_and_PPE_means(domain, 'dsi', index_number)
                # create_anomaly_cube(domain, model, 'dsi', index_number, 'TS1', 'TS3')


def plot_all_anomaly_plots():
    for domain in pre_processing.DOMAIN:
        for model in pre_processing.CMIP_MODELS + pre_processing.PPE_MEMBERS + pre_processing.MEANS:
            for number in DSI.INDEX_NUMBERS:
                # plot_individual_TS3_TS1_anomaly_cubes(domain, model, 'dsi', number)
                create_TS3_TS1_anomaly_plots_with_comparison(domain, model, 'dsi', number)
                # for ensemble in pre_processing.ENSEMBLE:
                    # plot_all_models_anomaly_cubes(domain, 'dsi', number, ensemble, with_means=False)
                    # plot_all_models_anomaly_cubes(domain, 'dsi', number, ensemble, with_means=True)

def plot_all_domain_comparisons_anomalies():
    for number in DSI.INDEX_NUMBERS:
        plot_domain_comparison_anomalies('CMIP_PPE_anomaly_difference', number)
    # combinations = itertools.product(
    #     pre_processing.CMIP_MODELS + pre_processing.PPE_MEMBERS + pre_processing.MEANS,
    #     DSI.INDEX_NUMBERS
    # )
    #
    # for model, index_number in combinations:
    #     plot_domain_comparison_anomalies(model, index_number)

def main():
    # create_TS3_minus_TS1_anomalies()
    # plot_all_anomaly_plots()
    plot_all_domain_comparisons_anomalies()

if __name__ == '__main__':
    main()
