import numpy as np
import iris
import itertools
import pre_processing
import matplotlib.pyplot as plt
import DSI
import os
import csv

def find_coords(cube):
    if 'grid_latitude' in [coord.name() for coord in cube.dim_coords]:
        coords = ('grid_latitude', 'grid_longitude')
    elif 'projection_x_coordinate' in [coord.name() for coord in cube.dim_coords]:
        coords = ('projection_x_coordinate', 'projection_y_coordinate')

    return coords


def flatten_cube_data(cube):
    data = cube.data
    flattened_data = np.ndarray.flatten(data)
    return flattened_data


def find_pattern_correlation(temporal_mean_model_cube, temporal_mean_obs_cube):
    flattened_model_data = flatten_cube_data(temporal_mean_model_cube)
    flattened_obs_data = flatten_cube_data(temporal_mean_obs_cube)

    pattern_correlation_matrix = np.ma.corrcoef(flattened_model_data, flattened_obs_data)
    pattern_correlation = pattern_correlation_matrix[0, 1]

    return pattern_correlation


def find_normalised_spatial_standard_deviation(temporal_mean_model_cube, temporal_mean_obs_cube):
    model_coords = find_coords(temporal_mean_model_cube)
    obs_coords = find_coords(temporal_mean_obs_cube)

    model_spatial_SD = temporal_mean_model_cube.collapsed(model_coords, iris.analysis.STD_DEV)
    obs_spatial_SD = temporal_mean_obs_cube.collapsed(obs_coords, iris.analysis.STD_DEV)

    normalised_SD_cube = model_spatial_SD / obs_spatial_SD
    normalised_SD = normalised_SD_cube.data

    return normalised_SD


def calculate_skill_score(model_cube, obs_cube):
    max_pattern_correlation = 1

    temporal_mean_model = model_cube.collapsed('time', iris.analysis.MEAN)
    temporal_mean_obs = obs_cube.collapsed('time', iris.analysis.MEAN)

    normalised_SD = find_normalised_spatial_standard_deviation(temporal_mean_model, temporal_mean_obs)
    pattern_correlation = find_pattern_correlation(temporal_mean_model, temporal_mean_obs)

    skill_numerator = 4 * (1 + pattern_correlation)

    skill_denominator = (normalised_SD + 1 / normalised_SD) ** 2 * (1 + max_pattern_correlation)

    skill_score = skill_numerator / skill_denominator

    return skill_score


def create_skill_score_dict(country):
    skill_score_dict = {}
    for number in DSI.INDEX_NUMBERS:
        firstmonth = str(number - 1).zfill(2)
        domain_dict = {}
        for domain in pre_processing.DOMAIN:
            model_dict = {}
            res = pre_processing.RESOLUTION[domain]
            obs_cube = iris.load_cube(f'/net/data/users/mgrant/ukcp/droughts/data/{country}/obs/{res}/'
                                      f'dsi-{number}/1980-2000/obs_{res}_1981{firstmonth}-200011_dsi-{number}.nc')
            for model in pre_processing.CMIP_MODELS + pre_processing.MEANS + pre_processing.MEDIANS:
                model_cube = iris.load_cube(f'/net/data/users/mgrant/ukcp/droughts/data/{country}/{domain}/{model}/TS1/'
                                            f'monthly/dsi-{number}/{model}_1981{firstmonth}-200011_dsi-{number}.nc')

                skill_score = calculate_skill_score(model_cube, obs_cube)
                model_dict[model] = skill_score
            domain_dict[domain] = model_dict
        skill_score_dict[f'dsi-{number}'] = domain_dict

    return skill_score_dict


def create_and_save_skill_score_csv(country):
    skill_score_dict = create_skill_score_dict(country)

    filepath = f'{pre_processing.data_files}/{country}/verification/SDPC_skill_score/summary_tables/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'{country}_SDPC_skill_score_summary_table.csv'

    with open(filepath + filename, 'w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=['Drought Index', 'Domain', 'Model', 'SDPC Skill Score'])
        writer.writeheader()
        index = 1
        for index_key, index_value in skill_score_dict.items():
            for dom_key, dom_value in index_value.items():
                for model_key, value in dom_value.items():
                    writer.writerow({'Drought Index': index_key,
                                     'Domain': dom_key,
                                     'Model': model_key,
                                     'SDPC Skill Score': round(value, 2)})
                    index += 1


def main():
    for country in pre_processing.COUNTRIES:
        create_and_save_skill_score_csv(country)

if __name__ == '__main__':
    main()
