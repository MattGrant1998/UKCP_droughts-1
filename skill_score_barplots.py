import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

import DSI
import pre_processing
import itertools


skill_scores = [
    'anomaly_from_obs',
    'MAE'
    'SDPC_skill_score'
]

comparisons = [
    'global_vs_downscaled',
    'CMIP_vs_PPE'
]

model_list = {
    'global_vs_downscaled': pre_processing.CMIP_MODELS + ['CMIP Mean'],
    'CMIP_vs_PPE': ['CMIP Mean', 'PPE Mean', 'CMIP + PPE Mean', 'CMIP Median', 'PPE Median', 'CMIP + PPE Median']
}


def change_ensemble_average_labels_in_df(dataframe):
    replacement_map = {
        'CMIP_multimodel_mean': 'CMIP Mean',
        'PPE_multimodel_mean': 'PPE Mean',
        'full_multimodel_mean': 'CMIP + PPE Mean',
        'CMIP_multimodel_median': 'CMIP Median',
        'PPE_multimodel_median': 'PPE Median',
        'full_multimodel_median': 'CMIP + PPE Median'
    }

    dataframe['Model'] = dataframe['Model'].replace(replacement_map)

    return dataframe


def change_domain_labels_in_df(dataframe):
    replacement_map = {
        'global': 'GCM',
        'regional': 'RCM',
        'local': 'CPM'
    }

    dataframe['Domain'] = dataframe['Domain'].replace(replacement_map)

    return dataframe

def take_reciprocal_of_skill_score(dataframe, score_title):
    # Take the reciprocal of each entry in the score_title column
    dataframe[score_title] = dataframe[score_title].rdiv(1)

    dataframe.rename(columns={score_title:f'1 / ({score_title})'}, inplace=True)

    return dataframe


def create_plot_of_skill_scores(score_type, country, number, comparison, average_type):
    skill_score_df = pd.read_csv(f'/data/users/mgrant/ukcp/droughts/data/{country}/verification/'
                                 f'{score_type}/summary_tables/{country}_{score_type}_summary_table.csv')
    skill_score_df = change_ensemble_average_labels_in_df(skill_score_df)
    DSInumber_df = skill_score_df[skill_score_df['Drought Index'] == f'dsi-{number}']

    score_title = list(DSInumber_df.columns.values)[3]

    # if score_title == 'Anomaly from Obs':
    #     # Take reciprocal of anomaly so that the lowest anomaly shows highest on the bar chart
    #     DSInumber_df = take_reciprocal_of_skill_score(DSInumber_df, score_title)
    #     score_title = list(DSInumber_df.columns.values)[3]

    plt.figure()
    if country == 'uk':
        Country = 'the UK'
    else:
        Country = country

    # if score_type == 'SDPC_skill_score':
    #     score_title = 'SDPC Skill Score'
    # elif score_type == 'anomaly_from_obs':
    #     score_title = 'Anomaly from Obs'
    # else:
    #     raise ValueError(f'{score_type} not accepted as argument in the function')


    plt.suptitle(f'{score_title} for DSI-{number} Across {Country}', fontsize=20)

    if comparison == 'global_vs_downscaled':
        x_list = ['MRI-CGCM3', 'MPI-ESM-LR', 'ACCESS1-3', 'IPSL-CM5A-MR', 'CMIP Mean']
        x_df = DSInumber_df.loc[DSInumber_df["Model"].isin(x_list)]['Model']
        hue_df = DSInumber_df['Domain']
        colors = ['bisque', 'salmon', 'paleturquoise']
    elif comparison == 'CMIP_vs_PPE':
        hue_list = [f'CMIP {average_type}', f'PPE {average_type}', f'CMIP + PPE {average_type}']

        x_df = DSInumber_df.loc[DSInumber_df["Model"].isin(hue_list)]['Domain']
        hue_df = DSInumber_df.loc[DSInumber_df["Model"].isin(hue_list)]['Model']
        colors = ['powderblue', 'deepskyblue', 'royalblue']
    # print(comparison)
    # print(x_df)
    # print(hue_df)
    ax = sns.barplot(
        data=DSInumber_df, x=x_df, y=DSInumber_df[score_title],
        hue=hue_df, palette=colors
    )

    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3)
    plt.xticks(rotation=30)

    if score_type == 'anomaly_from_obs':
        ax.set_ylabel(str(list(DSInumber_df.columns.values)[3]) + ' (%)', fontsize=15)
    else:
        ax.set_ylabel(str(list(DSInumber_df.columns.values)[3]), fontsize=15)

    if score_type == 'SDPC_skill_score':
        ax.set_ylim(0, 1)

    ax.set_xlabel('')
    plt.setp(ax.get_xticklabels(), ha="right", fontsize=15)
    plt.subplots_adjust(bottom=0.25, top=0.85)

    if comparison == 'global_vs_downscaled':
        filepath = f'{pre_processing.plots_files}/verification_barcharts/{comparison}/' \
                   f'{score_type}/{country}/dsi-{number}/'
    elif comparison == 'CMIP_vs_PPE':
        filepath = f'{pre_processing.plots_files}/verification_barcharts/{comparison}/' \
                   f'{score_type}/{average_type}/{country}/dsi-{number}/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'{country}_{score_type}_dsi-{number}_barchart.png'
    plt.savefig(filepath + filename)
    plt.close()


# create_plot_of_skill_scores('SDPC_skill_score', 'Scotland', 3, 'global_vs_downscaled')
# create_plot_of_skill_scores('SDPC_skill_score', 'uk', 3, 'CMIP_vs_PPE')


def create_plot_of_skill_scores_all_dsi(score_type, country, comparison, average_type='Mean', for_paper=False):
    plt.figure(figsize=(24, 6))
    skill_score_df = pd.read_csv(f'/data/users/mgrant/ukcp/droughts/data/{country}/verification/'
                                 f'{score_type}/summary_tables/{country}_{score_type}_summary_table.csv')
    skill_score_df = change_ensemble_average_labels_in_df(skill_score_df)
    score_title = list(skill_score_df.columns.values)[3]

    models = model_list[comparison]
    df_for_finding_max = skill_score_df.loc[skill_score_df["Model"].isin(models)]
    max_error = df_for_finding_max[score_title].max()

    print('max error: ', max_error)
    if for_paper:
        skill_score_df = change_domain_labels_in_df(skill_score_df)
        legend_fontsize = 25
        title_fontsize = 35
        ylabel_fontsize = 30
        xtick_fontsize = 25
    else:
        legend_fontsize = 15
        title_fontsize = 25
        ylabel_fontsize = 20
        xtick_fontsize = 15

    for i, number in enumerate(DSI.INDEX_NUMBERS):
        DSInumber_df = skill_score_df[skill_score_df['Drought Index'] == f'dsi-{number}']

        # if score_title == 'Anomaly from Obs':
        #     # Take reciprocal of anomaly so that the lowest anomaly shows highest on the bar chart
        #     DSInumber_df = take_reciprocal_of_skill_score(DSInumber_df, score_title)
        #     score_title = list(DSInumber_df.columns.values)[3]

        if comparison == 'global_vs_downscaled':
            x_list = ['MRI-CGCM3', 'MPI-ESM-LR', 'ACCESS1-3', 'IPSL-CM5A-MR', 'CMIP Mean']
            x_df = DSInumber_df.loc[DSInumber_df["Model"].isin(x_list)]['Model']
            hue_df = DSInumber_df['Domain']
            colors = ['bisque', 'salmon', 'paleturquoise']
        elif comparison == 'CMIP_vs_PPE':
            hue_list = ['CMIP ' + average_type, 'PPE ' + average_type, 'CMIP + PPE ' + average_type]
            x_df = DSInumber_df.loc[DSInumber_df["Model"].isin(hue_list)]['Domain']
            hue_df = DSInumber_df.loc[DSInumber_df["Model"].isin(hue_list)]['Model']
            print(x_df, hue_df)
            colors = ['powderblue', 'deepskyblue', 'royalblue']


        # print(comparison)
        # print(x_df)
        # print(hue_df)
        plt.subplot2grid((1, 3), (0, i))
        ax = sns.barplot(
            data=DSInumber_df, x=x_df, y=DSInumber_df[score_title],
            hue=hue_df, palette=colors
        )

        if number == 12:
            if comparison == 'global_vs_downscaled':
                h_alignment = 1.5
            else:
                h_alignment = 1.95
            ax.legend(loc='upper right', bbox_to_anchor=(h_alignment, 1.1), fontsize=legend_fontsize)
        else:
            plt.legend([], [], frameon=False)
        # ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3, fontsize=15)
        plt.xticks(rotation=30)

        if country == 'uk':
            Country = 'the UK'
        else:
            Country = country

        if score_type == 'MAE' and comparison != 'CMIP_vs_PPE':
            plt.title(f'{score_title} for DSI-{number} Across {Country}', fontsize=title_fontsize, y=1.1)
        elif score_type == 'SDPC_skill_score':
            plt.title(f'{score_title} for \n DSI-{number} Across {Country}', fontsize=title_fontsize, y=1.1)
        else:
            plt.title(f'{score_title} for DSI-{number} \n Across {Country}', fontsize=title_fontsize, y=1.1)

        if score_type == 'SDPC_skill_score':
            y_max = 1
            plt.locator_params(axis='y', nbins=2)
        else:
            y_max = max_error + 0.2
            print('y_max: ', y_max)
            plt.locator_params(axis='y', nbins=6)


        ax.set_ylim(0, y_max)

        if score_type == 'anomaly_from_obs':
            ax.set_ylabel(score_title + ' (%)', fontsize=ylabel_fontsize)
        else:
            ax.set_ylabel(score_title, fontsize=ylabel_fontsize)

        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), ha="right", fontsize=xtick_fontsize)
        plt.setp(ax.get_yticklabels(), fontsize=25)

        plt.subplots_adjust(bottom=0.25, top=0.85, right=0.95, left=0.1)

        # plt.subplots_adjust(hspace=1.5)
        plt.tight_layout()
    # plt.suptitle(f'DSI {score_title} Across {Country}', fontsize=20)

    filepath = f'{pre_processing.plots_files}/verification_barcharts/{comparison}/' \
               f'{score_type}/{country}/all_dsi/'
    if comparison == 'global_vs_downscaled':
        filepath = f'{pre_processing.plots_files}/verification_barcharts/{comparison}/' \
                   f'{score_type}/{country}/all_dsi/'
    elif comparison == 'CMIP_vs_PPE':
        filepath = f'{pre_processing.plots_files}/verification_barcharts/{comparison}/' \
                   f'{score_type}/{average_type}/{country}/all_dsi/'

    if not os.path.exists(filepath):
        os.makedirs(filepath)

    filename = f'{country}_{score_type}_all_dsi_barchart.png'
    plt.savefig(filepath + filename)
    plt.close()

create_plot_of_skill_scores_all_dsi('MAE', 'Scotland', 'CMIP_vs_PPE', average_type='Median', for_paper=True)
def create_all_skill_score_plots():
    combinations = itertools.product(
        ['MAE', 'SDPC_skill_score'],
        pre_processing.COUNTRIES,
        DSI.INDEX_NUMBERS,
        comparisons
    )

    for skill_score, country, number, comparison in combinations:
        create_plot_of_skill_scores(skill_score, country, number, comparison, average_type='Median')
        create_plot_of_skill_scores_all_dsi(skill_score, country, comparison, average_type='Median', for_paper=True)


def main():
    create_all_skill_score_plots()


# if __name__ == '__main__':
#     main()

