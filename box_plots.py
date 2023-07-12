import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import numpy as np
import os
import DSI
import drought_index_plots
import make_country_cubes
import pre_processing
import pickle
import seaborn as sns
import matplotlib
import pandas as pd
from matplotlib.lines import Line2D
import itertools


def find_spatial_coord_names(cube):
    if ('projection_x_coordinate' and 'projection_y_coordinate') in [coord.name() for coord in cube.dim_coords]:
        spatial_coords = ['projection_y_coordinate', 'projection_x_coordinate']
    elif ('grid_latitude' and 'grid_longitude') in [coord.name() for coord in cube.dim_coords]:
        spatial_coords = ['grid_latitude', 'grid_longitude']
    else:
        raise ValueError(
            f'Dimension coordinates of the cube do not include either gird_longitude/latitude or projection '
            f'x/y. Likely incorrect cube is being passed through the function. '
            f'Dimension coordinates of the cube are: {cube.dim_coords}')

    return spatial_coords


def find_cube_mean(cube):
    # Load in Ensemble data
    coords = cube.dim_coords
    mean_cube = cube.collapsed(coords, iris.analysis.MEAN)

    return mean_cube


def add_ensemble_member_coord(cube, model):
    if 'ensemble_member' not in [coord.name() for coord in cube.dim_coords or cube.aux_coords]:
        if len(model) > 2:
            ensemble_member = pre_processing.MEMBERS[model]
        else:
            ensemble_member = int(model)

        ensemble_member_coord = iris.coords.AuxCoord(ensemble_member, standard_name=None, units='1',
                                                     long_name='ensemble_member', var_name='ensemble_member')
        cube.add_aux_coord(ensemble_member_coord)


def remove_unnecessary_aux_coords(cube):
    coord_names = []
    for coord in cube.aux_coords or cube.dim_coords:
        if not coord.name() == 'ensemble_member':
            coord_names.append(coord.name())

    for coord in coord_names:
        cube.remove_coord(coord)
    # if 'ensemble_member_id' in [coord.name() for coord in cube.aux_coords]:
    #     cube.remove_coord('ensemble_member_id')
    # if 'forecast_reference_time' in [coord.name() for coord in cube.aux_coords]:
    #     cube.remove_coord('forecast_reference_time')
    # if 'time' in [coord.name() for coord in cube.aux_coords]:
    #     cube.remove_coord('time')


def create_dictionary_of_means_for_each_ensemble_set(country, domain, timeslice, drought_index, index_number):
    firstmonth = str(index_number - 1).zfill(2)
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]

    means_dict = {}
    for ensemble in pre_processing.ENSEMBLE:
        models = pre_processing.ENSEMBLE[ensemble]
        means_dict[ensemble] = iris.cube.CubeList(())
        for model in models:
            if isinstance(model, int):
                model = str(model).zfill(2)

            filepath = f'{pre_processing.data_files}/{country}/{domain}/{model}/{timeslice}/monthly/' \
                       f'{drought_index}-{index_number}/'
            filename = f'{model}_{firstyear}{firstmonth}-{lastyear}11_{drought_index}-{index_number}.nc'
            cube = iris.load_cube(filepath + filename)
            cube_mean = find_cube_mean(cube)
            add_ensemble_member_coord(cube_mean, model)
            remove_unnecessary_aux_coords(cube_mean)
            means_dict[ensemble].append(cube_mean)

        iris.util.equalise_attributes(means_dict[ensemble])
        means_dict[ensemble] = means_dict[ensemble].merge_cube()

    dict_filepath = f'/scratch/mgrant/UKCP/droughts/data/{country}/boxplots_data/ensemble_mean_dictionary/'
    if not os.path.exists(dict_filepath):
        os.makedirs(dict_filepath)
    dict_filename = f'ensemble_mean_dictionary_{domain}_{timeslice}_{drought_index}-{index_number}.pickle'

    with open(dict_filepath + dict_filename, 'wb') as file:
        pickle.dump(means_dict, file)


def save_all_ensemble_dictionaries():
    for country in pre_processing.COUNTRIES:
        for domain in pre_processing.DOMAIN:
            for timeslice in pre_processing.SLICES:
                for drought_index in pre_processing.DROUGHT_INDICES:
                    for index_number in DSI.INDEX_NUMBERS:
                        create_dictionary_of_means_for_each_ensemble_set(
                            country, domain, timeslice, drought_index, index_number
                        )


# def save_all_ensemble_dictionaries():
#     combinations = itertools.product(pre_processing.COUNTRIES,
#                                      pre_processing.DOMAIN,
#                                      pre_processing.SLICES,
#                                      pre_processing.DROUGHT_INDICES,
#                                      DSI.INDEX_NUMBERS)
#
#     for country, domain, timeslice, drought_index, index_number in combinations:
#         create_dictionary_of_means_for_each_ensemble_set(country, domain, timeslice, drought_index, index_number)


# save_all_ensemble_dictionaries()

# cube = iris.load_cube(f'{pre_processing.data_files}/uk/regional/IPSL-CM5A-MR/TS1/monthly/dsi-3/IPSL-CM5A-MR_198102-200011_dsi-3.nc')
# cube_mean = find_cube_mean(cube)
# print(cube_mean)


def change_domain_labels_in_df(dataframe):
    replacement_map = {
        'global': 'GCM',
        'regional': 'RCM',
        'local': 'CPM'
    }

    dataframe['Domain'] = dataframe['Domain'].replace(replacement_map)

    return dataframe


def make_box_plots(country, timeslice, drought_index, index_number, for_paper=False):

    matplotlib.use('Agg')
    firstyear = pre_processing.SLICES[timeslice][0]
    lastyear = pre_processing.SLICES[timeslice][-1]
    arrays = []
    labels = []
    # raw_data = []
    if timeslice == 'TS1':
        firstmonth = str(index_number - 1).zfill(2)
        obs_cube = iris.load_cube(f'{pre_processing.data_files}/{country}/obs/1km/{drought_index}-{index_number}/'
                                  f'1980-2000/obs_1km_1981{firstmonth}-200011_{drought_index}-{index_number}.nc')
        obs_mean = find_cube_mean(obs_cube)
        obs_mean_data = obs_mean.data

    for domain in pre_processing.DOMAIN:
        print(domain)
        dict_filepath = f'/scratch/mgrant/UKCP/droughts/data/{country}/boxplots_data/ensemble_mean_dictionary/'
        dict_filename = f'ensemble_mean_dictionary_{domain}_{timeslice}_{drought_index}-{index_number}.pickle'
        with open(dict_filepath + dict_filename, 'rb') as file:
            domain_ensemble_dict = pickle.load(file)
        for ensemble in reversed(pre_processing.ENSEMBLE.keys()):
            data = domain_ensemble_dict[ensemble].data
            print(timeslice, drought_index, '-', index_number, ensemble)
            print(data)
            arrays.append(data)

            # raw_data = np.concatenate((raw_data, data))
            # raw_data.extend(data)
            if ensemble == 'full':
                ensemble_name = 'CMIP + PPE'
            else:
                ensemble_name = ensemble

            if for_paper:
                domain_dict = {
                    'global': 'GCM',
                    'regional': 'RCM',
                    'local': 'CPM'
                }
                labels.append(ensemble_name + ' ' + domain_dict[domain])
            else:
                labels.append(domain[0].upper() + domain[1:] + ' ' + ensemble_name)


    colors = {
        'global': 'goldenrod',
        'regional': 'red',
        'local': 'seagreen'
    }

    # for array in arrays:
    #     print(np.shape(array))
    colors_edge = [colors['global']] * 3 + [colors['regional']] * 3 + [colors['local']] * 3
    colors_whiskers = [colors['global']] * 6 + [colors['regional']] * 6 + [colors['local']] * 6
    colors_median = ['darkgoldenrod'] * 3 + ['darkred'] * 3 + ['green'] * 3

    # # Define the number of box plots per group
    # plots_per_group = 3
    #
    # # Calculate the total number of groups
    # total_groups = len(data) // plots_per_group
    #
    # # Generate positions for the box plots
    # positions = np.arange(total_groups) * plots_per_group + 1
    positions = [2, 4, 6, 10, 12, 14, 18, 20, 22]

    ppe_markers = ["o"]
    cmip_markers = ["P", "^", "s", "X"]
    marker_groups = {
        4: cmip_markers,
        12: ppe_markers * 12,
        16: ppe_markers * 12 + cmip_markers
    }

    markers = {
        4: cmip_markers,
        12: ppe_markers,
        16: ppe_markers + cmip_markers
    }

    cmip_group = ['cmip'] * 4
    ppe_group = ['ppe'] * 12
    size_groups = {
        4: cmip_group,
        12: ppe_group,
        16: ppe_group + cmip_group
    }

    cmip_size = 160
    ppe_size = 70
    size = {
        4: (cmip_size, ppe_size),
        12: (ppe_size, cmip_size),
        16: (cmip_size, ppe_size)
    }

    fig, ax = plt.subplots(figsize=(16, 10))
    print('Started making boxplots')
    # ax = fig.add_axes([0, 0, 1, 1])
    if country == 'uk':
        country_name = f'the {country.upper()}'
    else:
        country_name = country
    plt.suptitle(f'Boxplots of the Mean {drought_index.upper()} values across \n'
                 f'{country_name} from {firstyear}-{lastyear}', fontsize=35)


    bp = ax.boxplot(arrays, positions=positions, labels=labels, patch_artist=True, showfliers=False, whis=(0, 100))
    plt.xlim([1, 23])

    for i, data in enumerate(arrays):
        arrays[i] = data[~np.isnan(data)]

    for patch, color in zip(bp['boxes'], colors_edge):
        patch.set_edgecolor(color)
        patch.set_facecolor('white')

    for whisker, color in zip(bp['whiskers'], colors_whiskers):
        whisker.set(color=color)

    for cap, color in zip(bp['caps'], colors_whiskers):
        cap.set(color=color)

    for median, color in zip(bp['medians'], colors_median):
        median.set(color=color)

    # sns.stripplot(data=arrays, x=[x_position] * len(data), hue=labels, doge=True, jitter=False, color='black', ax=ax)
    for i in range(len(arrays)):
        facecolor = ['white'] * len(arrays[i])

        if i in range(0, 3):
            edgecolor = [colors['global']] * len(arrays[i])
            if i != 0:
                facecolor[0] = colors['global']
        elif i in range(3, 6):
            edgecolor = [colors['regional']] * len(arrays[i])
            if i != 3:
                facecolor[0] = colors['regional']
        elif i in range(6, 9):
            edgecolor = [colors['local']] * len(arrays[i])
            if i != 6:
                facecolor[0] = colors['local']
        else:
            raise ValueError(f'{i} outside of expected range of 9')

        data = pd.DataFrame(
            dict(y=arrays[i],
                 x=np.repeat(positions[i] - 0.5, len(arrays[i])),
                 style=marker_groups[len(arrays[i])],
                 # m=markers[len(arrays[i])],
                 e=edgecolor,
                 f=facecolor,
                 s_groups=size_groups[len(arrays[i])],)
                 # s=size[len(arrays[i])])
        )

        sns.scatterplot(
            data=data, x='x', y='y', style=marker_groups[len(arrays[i])], markers=markers[len(arrays[i])], legend=False,
            size=size_groups[len(arrays[i])], sizes=size[len(arrays[i])], c=data['f'], edgecolor=data['e'], linewidth=1,
        )

        # y_position = arrays[i]
        # x_position = np.repeat(positions[i] - 0.5, len(y_position))
        # for x, y, marker in zip(x_position, y_position, markers[len(y_position)]):
        #     if marker == 'o':
        #         size = 10
        #     else:
        #         size = 30
        #
        #     plt.scatter(
        #         x_position, y_position, c='None', s=size,
        #         edgecolor=color, alpha=0.2, marker=marker
        #     )

    # x_values = [i + 1 for i in range(len(raw_data))]  # Generate x-coordinates for data points
    # ax.scatter(obs_mean_data, color='black', alpha=0.5)

    ax.set_xlabel(None)
    ax.set_ylabel('DSI (%)', fontsize=25, labelpad=20)

    legend_markers = cmip_markers + ['o'] + ['o']
    legend_labels = pre_processing.CMIP_MODELS + ['PPE Member', 'Driving PPE Member']
    legend_edgecolors = ['k'] * 6
    legend_facecolors = ['white'] * 5 + ['k']
    handles = [Line2D([], [], marker=marker, linestyle='None', label=label,
                      markerfacecolor=facecolor, markersize=15, markeredgecolor=edgecolor)
               for marker, label, facecolor, edgecolor in
               zip(legend_markers, legend_labels, legend_facecolors, legend_edgecolors)]
    plt.legend(fontsize=20, handles=handles, ncol=1, loc='upper left', bbox_to_anchor=(1, 1))
    # handles, labels = plt.gca().get_legend_handles_labels()
    #
    # for handle in handles:
    #     handle.set_markers(20)



    if timeslice == 'TS1':
        # Add a horizontal line where the obs is
        line_value = obs_mean_data  # Specify the y-value for the line
        ax.axhline(line_value, color='darkblue', linestyle='--')
        ax.text(0.01, line_value, 'Obs', fontsize=20, color='darkblue', ha='left', va='center')

    plt.xticks(rotation=45, fontsize=25)
    plt.yticks(fontsize=20)
    plt.margins(0.2)
    plt.setp(ax.get_xticklabels(), ha="right")
    plt.subplots_adjust(bottom=0.2)
    # plt.show()

    print('Plots made')
    fig_filepath = f'{pre_processing.plots_files}/boxplots/{country}/{timeslice}/{drought_index}-{index_number}/'
    plt.subplots_adjust(bottom=0.25, top=0.85, right=0.74)

    if not os.path.exists(fig_filepath):
        os.makedirs(fig_filepath)

    fig_filename = f'boxplot_{country}_{timeslice}_{drought_index}-{index_number}.png'

    if for_paper:
        fig_filename = fig_filename[:-4] + '_paper_version.png'

    plt.savefig(fig_filepath + fig_filename)
    print('figure_saved')
    plt.close()

# make_box_plots('uk', 'TS1', 'dsi', 6, for_paper=True)
def plot_all_boxplots():
    combinations = itertools.product(
        pre_processing.COUNTRIES,
        pre_processing.SLICES,
        DSI.INDEX_NUMBERS,
    )

    for country, timeslice, number in combinations:
        make_box_plots(country, timeslice, 'dsi', number, for_paper=True)
    # for country in pre_processing.COUNTRIES:
    #     for timeslice in pre_processing.SLICES:
    #         for drought_index in pre_processing.DROUGHT_INDICES:
    #             for number in DSI.INDEX_NUMBERS:
    #                 make_box_plots(country, timeslice, drought_index, number, for_paper=True)


# def plot_all_boxplots():
#     combinations = itertools.product(
#         pre_processing.COUNTRIES,
#         pre_processing.SLICES,
#         pre_processing.DROUGHT_INDICES,
#         DSI.INDEX_NUMBERS
#     )
#
#     for country, timeslice, drought_index, index_number in combinations:
#         make_box_plots(country, timeslice, drought_index, index_number)


def main():
    # save_all_ensemble_dictionaries()
    plot_all_boxplots()


if __name__ == '__main__':
    main()

# make_box_plots('TS1', 'dsi', 3)