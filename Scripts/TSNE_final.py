#June 15, 2023
#Fig 3C and Supplementary Figure 8
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy import stats
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

import logging
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--path_to_predictions",
                    help="Relative path of enformer predictions (.npy)", type=str, required=True)
parser.add_argument("--cell_type",
                    help="Cell type for experiment", type=str, required=True)
parser.add_argument("--track_indices",
                    help="tracks to use for the cell type", nargs = '+', type=int, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()
os.mkdir(hparams.output_path)

def tsne_projection_data_annotation(tsne_data, evolved_name, naive_name, palette, output_path):
    #to annotate by label

    annot_evolved = [evolved_name for i in range(896000)]
    annot_naive = [naive_name for i in range(896000)]
    annot = annot_evolved + annot_naive
    tsne_data['annot'] = annot
    #to shuffle order of rows to ensure to proper plotting
    tsne_data = tsne_data.sample(frac = 1)

    fig, axes = plt.subplots(figsize=(5, 5))
    sns.scatterplot(x="comp-1", y="comp-2",hue='annot', palette=palette,  s=3, alpha=0.4, linewidth=0,
                    data=tsne_data).set(title=str(evolved_name) + " and " + str(naive_name) + " T-SNE projection (n = 1000)", rasterized=True)
    plt.legend(loc='upper right', title='')
    plt.xlabel('tSNE1')
    plt.ylabel('tSNE2')
    plt.savefig(output_path + '/TSNE_n100_normalized' + str(evolved_name) + '_' + str(naive_name) + '_' + 'shuffled_order.png', dpi=400, bbox_inches="tight")
    plt.close()


def tsne_projection_data_annotation_all(tsne_data, evolved_name, random_name, dilocal_name, palette, output_path):
    #to annotate by label
    annot_evolved = [evolved_name for i in range(896000)]
    annot_random = [random_name for i in range(896000)]
    annot_dilocal = [dilocal_name for i in range(896000)]

    annot = annot_evolved + annot_random + annot_dilocal
    tsne_data['annot'] = annot
    #to shuffle order of rows
    tsne_data = tsne_data.sample(frac = 1)
    fig, axes = plt.subplots(figsize=(5, 5))
    sns.scatterplot(x="comp-1", y="comp-2",hue='annot', palette=palette,  s=3, alpha=1, linewidth=0,
                    data=tsne_data).set(title="evolved, random and Dilocal T-SNE projection (n = 1000 loci)")
    plt.legend(loc='upper right', title='')
    plt.xlabel('tSNE1', fontsize=15)
    plt.ylabel('tSNE2', fontsize=15)
    plt.savefig(output_path + '/TSNE_n1000_normalized_all_data_log_alpha1.png', dpi=400, bbox_inches="tight")
    plt.close()

def tsne_projection_track_annotation(tsne_data, evolved_name, naive_name, evolved_tracks, naive_tracks, track_index, track_name, cmap, output_path):
    fig, axes = plt.subplots(figsize=(5, 5))
    #to annotate by track valuess
    annot_evolved = evolved_tracks[:,track_index].tolist()
    annot_naive = naive_tracks[:,track_index].tolist()
    annot = annot_evolved + annot_naive
    tsne_data['annot'] = annot

    sns.scatterplot(x="comp-1", y="comp-2",c = tsne_data.annot , cmap = cmap, s=3, alpha=1,linewidth=0,
                    data=tsne_data).set(title=str(evolved_name) + " and " + str(naive_name) + " T-SNE projection (n = 1000)", rasterized=True)
    norm = plt.Normalize(tsne_data.annot.min(), tsne_data.annot.max())
    sm = plt.cm.ScalarMappable(cmap = cmap, norm = norm)
    sm.set_array([])
    # Remove the legend and add a colorbar
    axes.figure.colorbar(sm, label = track_name)
    plt.xlabel('tSNE1')
    plt.ylabel('tSNE2')
    plt.savefig(output_path + '/TSNE_n1000_normalized' + str(evolved_name) + '_' + str(naive_name) + '_' + track_name + '.png', dpi=400, bbox_inches="tight")
    plt.close()

def tsne_projection_split(tsne_data, evolved_name, naive_name, evolved_tracks, naive_tracks, track_index, track_name, cmap, output_path):
    evolved = tsne_data.iloc[0:896000,:]

    evolved['annot'] = np.log2(evolved_tracks[:,track_index].tolist())
    naive = tsne_data.iloc[896000:1792000,:]

    naive['annot'] = np.log2(naive_tracks[:,track_index].tolist())
    #check code with removed
    naive['annot'] = naive['annot'].mask(naive['annot'] < 0, 0)
    evolved['annot'] = evolved['annot'].mask(evolved['annot'] < 0, 0)
    vmin = min(min(evolved.annot), min(naive.annot))
    vmax = max(max(evolved.annot), max(naive.annot))

    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(10, 5))
    sns.scatterplot(ax = ax1, x="comp-1", y="comp-2",c = evolved.annot , cmap = cmap, vmin=vmin, vmax=vmax, s=3, alpha=0.4, linewidth=0,
                    data=evolved, rasterized=True)
    ax1.set(xlabel='tSNE1', ylabel='tSNE2')
    ax1.set_title(evolved_name, fontsize=25)

    sns.scatterplot(ax=ax2, x="comp-1", y="comp-2",c = naive.annot , cmap = cmap,vmin=vmin, vmax=vmax, s=3, alpha=0.4, linewidth=0,
                    data=naive, rasterized=True)
    ax2.set(xlabel='tSNE1', ylabel='tSNE2')
    ax2.set_yticks([])
    ax2.set_ylabel('')
    ax2.set_title(naive_name, fontsize=25)
    fig.subplots_adjust(bottom=0.1, top=0.9, left=0.1, right=0.8, wspace=0.02, hspace=0.02)
    # Remove the legend and add a colorbar
    norm = plt.Normalize(vmin, vmax)
    sm = plt.cm.ScalarMappable(cmap = cmap, norm = norm)
    sm.set_array([])
    # Remove the legend and add a colorbar
    cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
    cb = fig.colorbar(sm, cb_ax, label = track_name)
    tick_list = cb.ax.get_yticks()
    power = [pow(2, item) for item in tick_list]
    power = [ '%.1f' % elem for elem in power ]
    cb.ax.set_yticklabels(power, fontsize=20)
    ax2.xaxis.label.set_size(25)
    ax2.yaxis.label.set_size(25)
    ax1.xaxis.label.set_size(25)
    ax1.yaxis.label.set_size(25)
    ax2.tick_params(labelsize=20)
    ax1.tick_params(labelsize=20)
    plt.savefig(output_path + '/TSNE_n1000_normalized' + evolved_name + '_' + naive_name + '_' + track_name + '_split' + '.png', dpi=400, bbox_inches="tight")
    plt.close()


### Read in files
random_all = np.load(hparams.path_to_predictions +'/random_sequences.npy')[:,:,0:1000]
evolved_all = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
di_local_all = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]

cmap = sns.color_palette("flare", as_cmap=True)
green = (0/255, 102/255, 0/255)
purple = (128/255, 0/255, 128/255)
blue = (0/255, 0/255, 238/255)

palette = {'evolved':purple,
           'random':blue,
           'naive': green}


cell_type_index = hparams.track_indices
# #fill in indcies of interest for current cell types
evolved_cell = evolved_all[:,cell_type_index,:]
evolved_cell_test = np.array(evolved_cell[:,:,list(range(0, 1000, 1))])
evolved_test_regions_cell_melt = evolved_cell_test.transpose(2,0,1).reshape(-1,evolved_cell_test.shape[1])

random_cell = random_all[:,cell_type_index,:]
random_cell_test = np.array(random_cell[:,:,list(range(0, 1000, 1))])
random_test_regions_cell_melt = random_cell_test.transpose(2,0,1).reshape(-1,random_cell_test.shape[1])

dilocal_cell = di_local_all[:,cell_type_index,:]
dilocal_cell_test = np.array(dilocal_cell[:,:,list(range(0, 1000, 1))])
dilocal_test_regions_cell_melt = dilocal_cell_test.transpose(2,0,1).reshape(-1,dilocal_cell_test.shape[1])

#to Calculate TSNE projections

#TSNE plotting evolved vs dilocal#########################################################################
#combine the evolved and dilocal one_matrix
evolved_dilocal_stack = np.vstack((evolved_test_regions_cell_melt,dilocal_test_regions_cell_melt ))
tsne = TSNE(n_components=2, verbose=1, random_state=123, learning_rate = "auto")
evolved_dilocal_stack = StandardScaler().fit_transform(evolved_dilocal_stack)
z = tsne.fit_transform(evolved_dilocal_stack)
evolved_dilocal = pd.DataFrame()
evolved_dilocal["comp-1"] = z[:,0]
evolved_dilocal["comp-2"] = z[:,1]
#only need to calculate this once can read it in next time to save time
evolved_dilocal.to_csv(hparams.output_path + "/genomic_dilocal_normalized_n1000.csv")

evolved_dilocal = pd.read_csv(hparams.output_path + "/genomic_dilocal_normalized_n1000.csv")

tsne_projection_data_annotation(evolved_dilocal, "evolved", "naive", palette, hparams.output_path)

tsne_projection_track_annotation(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 0, "DNase", cmap, hparams.output_path)
tsne_projection_split(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 0, "DNase", cmap, hparams.output_path)

tsne_projection_track_annotation(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 1, "H3K27ac", cmap, hparams.output_path)
tsne_projection_split(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 1, "H3K27ac", cmap, hparams.output_path)

tsne_projection_track_annotation(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 2, "H3K4me1", cmap, hparams.output_path)
tsne_projection_split(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 2, "H3K4me1", cmap, hparams.output_path)

tsne_projection_track_annotation(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, random_test_regions_cell_melt, 3, "H3K4me3", cmap, hparams.output_path)
tsne_projection_split(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 3, "H3K4me3", cmap, hparams.output_path)

tsne_projection_track_annotation(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 4, "H3K27me3", cmap, hparams.output_path)
tsne_projection_split(evolved_dilocal, "evolved", "naive", evolved_test_regions_cell_melt, dilocal_test_regions_cell_melt, 4, "H3K27me3", cmap, hparams.output_path)

###to tsne for all data at same time
all_data_stack = np.vstack((evolved_test_regions_cell_melt,random_test_regions_cell_melt,dilocal_test_regions_cell_melt ))
tsne = TSNE(n_components=2, verbose=1, random_state=123, learning_rate = "auto")
all_data_stack = StandardScaler().fit_transform(all_data_stack)
z = tsne.fit_transform(all_data_stack)
all_data = pd.DataFrame()
all_data["comp-1"] = z[:,0]
all_data["comp-2"] = z[:,1]
all_data.to_csv(hparams.output_path + "/all_data_normalized_n1000.csv")
all_data = pd.read_csv(hparams.output_path + "/all_data_normalized_n1000.csv")
tsne_projection_data_annotation_all(all_data, "evolved", "random", "naive", palette, hparams.output_path)
