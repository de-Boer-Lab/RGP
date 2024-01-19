##June 15, 2023
#Fig 3A - simplfied version (only plots track for genomic, random and local dinucleotide shuffled)
#This script will make plots for the enformer prediction track of interest for all 1000 regions
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns
from numpy import genfromtxt
import os
import logging
import sys
import argparse
from scipy.stats.stats import pearsonr

bases = ['A', 'T', 'C', 'G']
#enformer prediction window
SEQUENCE_LENGTH = 114_688

parser = argparse.ArgumentParser()

parser.add_argument("--path_to_predictions",
                    help="Relative path of enformer predictions (.npy)", type=str, required=True)
parser.add_argument("--track_to_plot",
                    help="Cell type and assay", type=str, required=True)
parser.add_argument("--track_index",
                    help="track to use for the cell type and assay you want to plot", type=int, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()
# os.mkdir(hparams.output_path)
path_to_predictions = '/scratch/st-cdeboer-1/iluthra/randomDNA/all_iPSC_related_H1/'
### Read in files
random_all = np.load(path_to_predictions +'/random_sequences.npy')[:,:,0:1000]
genomic_all = np.load(path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
di_local_all = np.load(path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]


def plot_track(ax, track, title, max_current):
    sns.color_palette("pastel")
    ax.fill_between(np.linspace(0, SEQUENCE_LENGTH, num=len(track)), track,zorder=0)
    sns.despine(top=True, right=True, bottom=True)
    ax.set_title(title)
    ax.set_xlabel('')
    ax.set_ylim(0, max_current)

#plot all 1000 predicted sequences
for i in (np.arange(0, 1000, 1).tolist()):
    track_genomic = genomic_all[:,40,i]
    track_random = random_all[:,40,i]
    track_dilocal = di_local_all[:,40,i]

    corr_coeff = pearsonr(track_genomic[:,], track_dilocal[:,])
    if corr_coeff[0] > 0.29:
        if corr_coeff[0] < 0.30:
            print(corr_coeff)
            print(i)
    track_genomic_max = track_genomic.max()
    track_random_max =track_random.max()
    track_dilocal_max = track_dilocal.max()
    #make the y-axis the same for all 3 tracks
    max_current = max([track_genomic_max, track_random_max, track_dilocal_max])

    # Plotting
    fig = plt.figure(figsize = (12, 4))
    gs = fig.add_gridspec(ncols=1, nrows=3, hspace=0.5, wspace=0.3)
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[1,0],sharex=ax0)
    ax2 = fig.add_subplot(gs[2,0],sharex=ax0)

    for ax in [ax0, ax1, ax2]:
        ax.set_xticks(np.arange(0, SEQUENCE_LENGTH, 20_000))
        ax.tick_params(labelbottom=False, bottom=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    ax2.set_xlim(0,SEQUENCE_LENGTH)
    ax2.set_xticks([0, 25_000, 50_000, 75_000,100_000])
    ax2.tick_params(labelbottom=True, bottom=True)
    ax2.set_xticklabels(np.arange(0, 0 + SEQUENCE_LENGTH, 25_000))

    plot_track(ax0, track_genomic, 'Genomic:' + hparams.track_to_plot, max_current)
    plot_track(ax1, track_random , 'Random:' + hparams.track_to_plot, max_current)
    plot_track(ax2, track_dilocal, 'Di-local:' + hparams.track_to_plot, max_current)

    fig.savefig(hparams.output_path + '/region_' + str(i) + '.png',  bbox_inches='tight')
