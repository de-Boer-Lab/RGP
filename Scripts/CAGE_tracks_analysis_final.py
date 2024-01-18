#June 15, 2023
#Fig 4D and Supplemental Figure 7
#This script will make density plots for cage track analysis
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy import stats
import logging
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--path_to_predictions",
                    help="Relative path of enformer predictions for all cage tracks(.npy)", type=str, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()
os.mkdir(hparams.output_path)
### Read in files
genomic_all = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
dinuc_local_all = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]

dinuc_local_all_regions = dinuc_local_all.transpose(2,0,1).reshape(-1,dinuc_local_all.shape[1])
genomic_all_regions = genomic_all.transpose(2,0,1).reshape(-1,genomic_all.shape[1])

print(dinuc_local_all_regions.shape)
genomic_mean = genomic_all_regions.mean(axis = 0)
genomic_std = genomic_all_regions.std(axis = 0)
genomic_z_score = (genomic_all_regions - genomic_mean)/genomic_std
dinuc_local_z_score = (dinuc_local_all_regions - genomic_mean)/genomic_std


# density plots of CAGE correlations no filtering based on z score
corr_genomic = pd.DataFrame(genomic_z_score).corr()
upper_genomic = pd.DataFrame(corr_genomic.where(np.triu(np.ones(corr_genomic.shape), k=1).astype(np.bool)))
all_corr_genomic = upper_genomic.values.tolist()
flat_list = [item for sublist in all_corr_genomic for item in sublist]
no_nan_genomic = [item for item in flat_list if not(pd.isnull(item)) == True]

corr_dinuc = pd.DataFrame(dinuc_local_z_score).corr()
upper_dinuc = pd.DataFrame(corr_dinuc.where(np.triu(np.ones(corr_dinuc.shape), k=1).astype(np.bool)))
all_corr_dinuc = upper_dinuc.values.tolist()
flat_list_dinuc = [item for sublist in all_corr_dinuc for item in sublist]
no_nan_dinuc = [item for item in flat_list_dinuc if not(pd.isnull(item)) == True]

all_data = pd.DataFrame(no_nan_genomic, columns = ['genomic'])
all_data['di-local'] = no_nan_dinuc
all_data.to_csv(hparams.output_path + '/all_CAGE_data.csv')

all_data_melt = all_data.melt()
all_data_melt.columns = ['data_type', 'correlations']

green = (0/255, 102/255, 0/255)
purple = (128/255, 0/255, 128/255)
blue = (0/255, 0/255, 238/255)

palette = {'genomic':purple,
           'random':blue,
           'di-local': green}

sns.set(font_scale = 1.5)
sns.set_theme(style='white')
#plot density plot
fig, ax = plt.subplots(figsize=(4,4))
sns.displot(data=all_data_melt, x="correlations", hue="data_type", kind="kde", palette=palette)
plt.xlabel('Track Correlations')
plt.ylabel('Density')
ax.legend(title='')
plt.savefig(hparams.output_path + '/distplot_correlation_all_tracks_genomic_naive_zscore.pdf', dpi=400, bbox_inches='tight')


def z_score_filt_plot(genomic, naive, z_score, index, axes):

    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)

    palette = {'genomic':purple,
               'random':blue,
               'di-local': green}

    #### Filtering by Z-score

    dinuc_filtered = pd.DataFrame(naive[(np.any(naive > z_score, axis = 1))])
    genomic_filtered = pd.DataFrame(genomic[(np.any(genomic > z_score, axis = 1))])
    print("z-value")
    print(z_score)
    print("genomic_filtered")
    print(genomic_filtered.shape)
    print("genomic_filtered")
    print(dinuc_filtered.shape)

    corr_genomic = pd.DataFrame(genomic_filtered).corr()
    upper_genomic = pd.DataFrame(corr_genomic.where(np.triu(np.ones(corr_genomic.shape), k=1).astype(np.bool)))
    all_corr_genomic = upper_genomic.values.tolist()
    flat_list = [item for sublist in all_corr_genomic for item in sublist]
    no_nan_genomic = [item for item in flat_list if not(pd.isnull(item)) == True]

    corr_dinuc = pd.DataFrame(dinuc_filtered).corr()
    upper_dinuc = pd.DataFrame(corr_dinuc.where(np.triu(np.ones(corr_dinuc.shape), k=1).astype(np.bool)))
    all_corr_dinuc = upper_dinuc.values.tolist()
    flat_list_dinuc = [item for sublist in all_corr_dinuc for item in sublist]
    no_nan_dinuc = [item for item in flat_list_dinuc if not(pd.isnull(item)) == True]

    all_data = pd.DataFrame(no_nan_genomic, columns = ['genomic'])
    all_data['di-local'] = no_nan_dinuc
    all_data.to_csv(hparams.output_path + '/z_filt_data_' + str(z_score) + '.csv')
    print("Rank sum test")
    stat, pvalue = stats.ranksums(all_data['genomic'], all_data['di-local'])
    print(pvalue)
    print(stat)
    label_p = '$\it{stat}$ = ' + str(stat)
    all_data_melt = all_data.melt()
    all_data_melt.columns = ['data_type', 'correlations']

    sns.kdeplot(ax=axes[index], data=all_data_melt, x="correlations", hue="data_type", palette=palette)
    axes[index].set_xlabel('Track Correlations')
    axes[index].set_ylabel('Density')
    axes[index].legend(labels = ['evolved (n = ' + str(genomic_filtered.shape[0]) + ')' , 'naive (n = ' + str(dinuc_filtered.shape[0]) + ')'], loc='upper left', fontsize=15)
    axes[index].set_title('Z-score > '+ str(z_score))
    axes[index].title.set_size(15)
    axes[index].tick_params(labelsize=15)
    #axes[index].annotate(label_p, xy=(.1, 0.95), xycoords=axes[index].transAxes, fontsize=15)

#create 5 panel cage track-track correlation plot with different z score cutoffs
sns.set_theme(style='white')
fig, axes = plt.subplots(1, 5, figsize=(20, 5))
z_score_filt_plot(genomic_z_score, dinuc_local_z_score, 1, 0, axes)
z_score_filt_plot(genomic_z_score, dinuc_local_z_score, 3, 1, axes)
z_score_filt_plot(genomic_z_score, dinuc_local_z_score, 10, 2, axes)
z_score_filt_plot(genomic_z_score, dinuc_local_z_score, 50, 3, axes)
z_score_filt_plot(genomic_z_score, dinuc_local_z_score, 100, 4, axes)
plt.savefig(hparams.output_path + '/distplot_correlation_all_tracks_z_scores_genomic_naive.pdf', dpi=300, bbox_inches="tight")
