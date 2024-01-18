#June 15, 2023
#Fig 4E, Extended Data Figure 9,10

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy import stats
from PIL import Image as im

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
                    help="tracks to use for the cell type should be in this order ['Dnase', 'H3K27ac','H3K4me1','H3K4me3', 'H3K27me3']", nargs = '+', type=int, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()
os.mkdir(hparams.output_path)
### Read in files
evolved_all = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
di_local_all = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]

def scatter_matrix_marks(evolved, naive, evolved_name, naive_name, index_list, cell_type, output_path):
    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)
    print(index_list[0])
    all_marks = pd.DataFrame(
        {'DNase evolved': evolved[:,index_list[0],:].flatten('F'),
         'H3K27ac evolved': evolved[:, index_list[1], :].flatten('F'),
         'H3K4me1 evolved': evolved[:, index_list[2], :].flatten('F'),
         'H3K4me3 evolved': evolved[:, index_list[3], :].flatten('F'),
         'H3K27me3 evolved': evolved[:, index_list[4], :].flatten('F'),
         'DNase naive': naive[:,index_list[0],:].flatten('F'),
         'H3K27ac naive': naive[:, index_list[1], :].flatten('F'),
         'H3K4me1 naive': naive[:, index_list[2], :].flatten('F'),
         'H3K4me3 naive': naive[:, index_list[3], :].flatten('F'),
         'H3K27me3 naive': naive[:, index_list[4], :].flatten('F'),
        })

    print(all_marks.shape)
    fig, axes = plt.subplots(5, 5, figsize=(50, 50))

    subplot('DNase naive', 'DNase evolved',0, 0, axes, all_marks)
    subplot('DNase naive', 'H3K27ac naive', 0, 1, axes, all_marks)
    subplot('DNase naive', 'H3K4me1 naive',0, 2, axes, all_marks)
    subplot( 'DNase naive', 'H3K4me3 naive', 0, 3, axes, all_marks)
    subplot('DNase naive', 'H3K27me3 naive', 0, 4, axes, all_marks)
    #COLUMN 2
    subplot('H3K27ac evolved', 'DNase evolved', 1, 0, axes, all_marks)
    subplot( 'H3K27ac naive', 'H3K27ac evolved', 1, 1, axes, all_marks)
    subplot('H3K27ac naive', 'H3K4me1 naive', 1, 2, axes, all_marks)
    subplot('H3K27ac naive','H3K4me3 naive', 1, 3, axes, all_marks)
    subplot('H3K27ac naive', 'H3K27me3 naive',1, 4, axes, all_marks)
    #COLUMN 3
    subplot('DNase evolved','H3K4me1 evolved', 2, 0, axes, all_marks)
    subplot('H3K27ac evolved', 'H3K4me1 evolved', 2, 1, axes, all_marks)
    subplot('H3K4me1 naive', 'H3K4me1 evolved', 2, 2, axes, all_marks)
    subplot( 'H3K4me1 naive', 'H3K4me3 naive', 2, 3, axes, all_marks)
    subplot( 'H3K4me1 naive', 'H3K27me3 naive',2, 4, axes, all_marks)
    #COLUMN 4
    subplot('DNase evolved', 'H3K4me3 evolved', 3, 0, axes, all_marks)
    subplot('H3K27ac evolved', 'H3K4me3 evolved', 3, 1, axes, all_marks)
    subplot('H3K4me1 evolved', 'H3K4me3 evolved', 3, 2, axes, all_marks)
    subplot( 'H3K4me3 naive', 'H3K4me3 evolved', 3, 3, axes, all_marks)
    subplot('H3K4me3 naive', 'H3K27me3 naive',3, 4, axes, all_marks)
    #COLUMN 5
    subplot('DNase evolved', 'H3K27me3 evolved', 4, 0, axes, all_marks)
    subplot('H3K27ac evolved', 'H3K27me3 evolved', 4, 1, axes, all_marks)
    subplot('H3K4me1 evolved', 'H3K27me3 evolved', 4, 2, axes, all_marks)
    subplot('H3K4me3 evolved', 'H3K27me3 evolved', 4, 3, axes, all_marks)
    subplot('H3K27me3 naive', 'H3K27me3 evolved',4, 4, axes, all_marks)


    plt.savefig(output_path + '/scatter_' + evolved_name + '_' + naive_name + '_all_marks' + cell_type + '.pdf', dpi=300, bbox_inches="tight")

def subplot(x_label, y_label, index_x, index_y, axes, data_frame):
    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)

    sns.scatterplot(ax=axes[index_x, index_y], data=data_frame, x=x_label, y=y_label, alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[index_x,index_y].xaxis.label.set_size(35)
    axes[index_x,index_y].yaxis.label.set_size(35)

    if "evolved" in x_label:
        axes[index_x,index_y].xaxis.label.set_color(purple)
    else:
        axes[index_x,index_y].xaxis.label.set_color(green)
    if "evolved" in y_label:
        axes[index_x,index_y].yaxis.label.set_color(purple)
    else:
        axes[index_x,index_y].yaxis.label.set_color(green)

    axes[index_x,index_y].tick_params(labelsize=30)
    (r, p) = stats.pearsonr(data_frame[x_label], data_frame[y_label])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[index_x,index_y].annotate(label, xy=(.6, 0.9), xycoords=axes[index_x,index_y].transAxes, fontsize=35)

def calc_corr(x_label, y_label, data_frame):
    (r, p) = stats.pearsonr(data_frame[x_label], data_frame[y_label])
    return(r)

def correlation_matrix(evolved, naive, evolved_name, naive_name, index_list, cell_type, output_path):
    all_marks = pd.DataFrame(
        {'DNase evolved': evolved[:,index_list[0],:].flatten(),
         'H3K27ac evolved': evolved[:, index_list[1], :].flatten(),
         'H3K4me1 evolved': evolved[:, index_list[2], :].flatten(),
         'H3K4me3 evolved': evolved[:, index_list[3], :].flatten(),
         'H3K27me3 evolved': evolved[:, index_list[4], :].flatten(),
         'DNase naive': naive[:,index_list[0],:].flatten(),
         'H3K27ac naive': naive[:, index_list[1], :].flatten(),
         'H3K4me1 naive': naive[:, index_list[2], :].flatten(),
         'H3K4me3 naive': naive[:, index_list[3], :].flatten(),
         'H3K27me3 naive': naive[:, index_list[4], :].flatten(),
        })

    print(all_marks.shape)
    fig, axes = plt.subplots(5, 5, figsize=(50, 50))
    data = np.empty([5,5])
    #COLUMN 1
    data[0,0] = calc_corr('DNase evolved', 'DNase naive', all_marks)
    data[0,1] = calc_corr('DNase naive', 'H3K27ac naive', all_marks)
    data[0,2] = calc_corr('DNase naive', 'H3K4me1 naive', all_marks)
    data[0,3] = calc_corr('DNase naive', 'H3K4me3 naive', all_marks)
    data[0,4] = calc_corr('DNase naive', 'H3K27me3 naive', all_marks)
    #COLUMN 2
    data[1,0] = calc_corr('H3K27ac evolved', 'DNase evolved', all_marks)
    data[1,1] = calc_corr('H3K27ac evolved', 'H3K27ac naive', all_marks)
    data[1,2] = calc_corr('H3K27ac naive', 'H3K4me1 naive', all_marks)
    data[1,3] = calc_corr('H3K27ac naive', 'H3K4me3 naive', all_marks)
    data[1,4] = calc_corr('H3K27ac naive', 'H3K27me3 naive', all_marks)
    #COLUMN 3
    data[2,0] = calc_corr('H3K4me1 evolved', 'DNase evolved', all_marks)
    data[2,1] = calc_corr('H3K4me1 evolved', 'H3K27ac evolved', all_marks)
    data[2,2] = calc_corr('H3K4me1 evolved', 'H3K4me1 naive', all_marks)
    data[2,3] = calc_corr('H3K4me1 naive', 'H3K4me3 naive', all_marks)
    data[2,4] = calc_corr('H3K4me1 naive', 'H3K27me3 naive', all_marks)
    #COLUMN 4
    data[3,0] = calc_corr('H3K4me3 evolved', 'DNase evolved', all_marks)
    data[3,1] = calc_corr('H3K4me3 evolved', 'H3K27ac evolved', all_marks)
    data[3,2] = calc_corr('H3K4me3 evolved', 'H3K4me1 evolved', all_marks)
    data[3,3] = calc_corr('H3K4me3 evolved', 'H3K4me3 naive', all_marks)
    data[3,4] = calc_corr('H3K4me3 naive', 'H3K27me3 naive',  all_marks)
    #COLUMN 5
    data[4,0] = calc_corr('H3K27me3 evolved','DNase evolved', all_marks)
    data[4,1] = calc_corr('H3K27me3 evolved', 'H3K27ac evolved', all_marks)
    data[4,2] = calc_corr('H3K27me3 evolved', 'H3K4me1 evolved', all_marks)
    data[4,3] = calc_corr('H3K27me3 evolved', 'H3K4me3 evolved', all_marks)
    data[4,4] = calc_corr('H3K27me3 evolved', 'H3K27me3 naive', all_marks)


    np.fill_diagonal(data, np.NaN)

    final_corrs = pd.DataFrame(data)
    final_corrs.columns = ['DNase', 'H3K27ac','H3K4me1','H3K4me3','H3K27me3']
    final_corrs.index = ['DNase', 'H3K27ac','H3K4me1','H3K4me3','H3K27me3']
    palette = sns.color_palette("magma", as_cmap=True)

    fig, axes = plt.subplots(figsize=(5, 5))
    axes.set_facecolor("lightgrey")
    sns.heatmap(final_corrs, cmap = palette, mask=final_corrs.isnull())
    plt.yticks(rotation=45)
    plt.xticks(rotation=45)
    plt.savefig(output_path + '/heatmap_' + evolved_name + '_' + naive_name + '_all_marks_correlations' + cell_type + '.pdf', dpi=300, bbox_inches="tight")

def barplot(evolved, naive, evolved_name, naive_name, index_list, cell_type, output_path):
    all_marks = pd.DataFrame(
        {'DNase evolved': evolved[:,index_list[0],:].flatten('F'),
         'H3K27ac evolved': evolved[:, index_list[1], :].flatten('F'),
         'H3K4me1 evolved': evolved[:, index_list[2], :].flatten('F'),
         'H3K4me3 evolved': evolved[:, index_list[3], :].flatten('F'),
         'H3K27me3 evolved': evolved[:, index_list[4], :].flatten('F'),
         'DNase naive': naive[:,index_list[0],:].flatten('F'),
         'H3K27ac naive': naive[:, index_list[1], :].flatten('F'),
         'H3K4me1 naive': naive[:, index_list[2], :].flatten('F'),
         'H3K4me3 naive': naive[:, index_list[3], :].flatten('F'),
         'H3K27me3 naive': naive[:, index_list[4], :].flatten('F'),
        })


    data = pd.DataFrame(columns = ['Col1', 'Col2', 'pearson', 'type_between'])

    #Dnase - H3K27ac
    correlation = calc_corr('DNase naive', 'H3K27ac naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K27ac'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K27ac evolved', 'DNase evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K27ac'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    #Dnase - H3K4me1

    correlation = calc_corr('DNase naive', 'H3K4me1 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K4me1'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K4me1 evolved', 'DNase evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K4me1'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    #DNase - H3K4me3

    correlation = calc_corr('DNase naive', 'H3K4me3 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K4me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K4me3 evolved', 'DNase evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K4me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    ## DNase - H3K27me3

    correlation = calc_corr('DNase naive', 'H3K27me3 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K27me3 evolved','DNase evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['DNase'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    #H3K27ac - H3K4me1

    correlation = calc_corr('H3K27ac naive', 'H3K4me1 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K27ac'],'Col2': ['H3K4me1'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K4me1 evolved', 'H3K27ac evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K27ac'],'Col2': ['H3K4me1'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    #H3K27ac - H3K4me3

    correlation = calc_corr('H3K27ac naive', 'H3K4me3 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K27ac'],'Col2': ['H3K4me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K4me3 evolved', 'H3K27ac evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K27ac'],'Col2': ['H3K4me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    #H3K27ac - H3K4me3

    correlation = calc_corr('H3K27ac naive', 'H3K27me3 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K27ac'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K27me3 evolved', 'H3K27ac evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K27ac'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    #all H3K4me1

    correlation = calc_corr('H3K4me1 naive', 'H3K4me3 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K4me1'],'Col2': ['H3K4me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K4me3 evolved', 'H3K4me1 evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K4me1'],'Col2': ['H3K4me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    # #All H3K4me3

    correlation = calc_corr('H3K4me1 naive', 'H3K27me3 naive', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K4me1'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K27me3 evolved', 'H3K4me1 evolved',  all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K4me1'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    # #All H3K27me3
    correlation = calc_corr('H3K4me3 naive', 'H3K27me3 naive',all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K4me3'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['naive']})
    data = pd.concat([data, data_2_add])

    correlation = calc_corr('H3K27me3 evolved', 'H3K4me3 evolved', all_marks)
    data_2_add = pd.DataFrame({'Col1': ['H3K4me3'],'Col2': ['H3K27me3'], 'pearson': [correlation], 'type_between':['evolved']})
    data = pd.concat([data, data_2_add])

    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)

    palette = {'evolved':purple,
               'Random':blue,
               'naive': green}
    data['Col1-Col2'] = data.Col1 + data.Col2
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={'height_ratios':[6,1]})
    #fig, axes = plt.subplots(figsize=(20, 5))
    sns.barplot(ax=axes[0], data=data,x="Col1-Col2", y="pearson", hue="type_between", palette = palette)
    #axes[0].xticks(rotation = 45, ha='right', rotation_mode='anchor')
    axes[0].set_xticks([])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('Correlation (Pearsons r)')
    axes[0].set_title(cell_type, fontsize=25)
    axes[0].yaxis.label.set_size(25)
    axes[0].tick_params(labelsize=25)
    #plt.xticks(rotation = 45, ha='right', rotation_mode='anchor')
    #plt.savefig('/scratch/st-cdeboer-1/iluthra/randomDNA//all_iPSC_related//box_plot_' + evolved_name + '_' + naive_name + '_all_marks_correlations' + cell_type + '.pdf', dpi=300, bbox_inches="tight")

    bubble_data = pd.DataFrame(columns = ['1', '2', '3', '4', '5', '6', '7','8','9','10'], index = ['DNase', 'H3K27ac','H3K4me1','H3K4me3','H3K27me3'])
    print(bubble_data)
    bubble_data['1'] = [1,1,0,0,0]
    bubble_data['2'] = [1,0,1,0,0]
    bubble_data['3'] = [1,0,0,1,0]
    bubble_data['4'] = [1,0,0,0,1]
    bubble_data['5'] = [0,1,1,0,0]
    bubble_data['6'] = [0,1,0,1,0]
    bubble_data['7'] = [0,1,0,0,1]
    bubble_data['8'] = [0,0,1,1,0]
    bubble_data['9'] = [0,0,1,0,1]
    bubble_data['10'] = [0,0,0,1,1]


    melted = bubble_data.melt(ignore_index=False)
    melted['track'] = melted.index
    sns.scatterplot(ax= axes[1], data=melted,x="variable", y="track", hue="value", s = 150, palette = ['lightgrey', 'black'], legend = False, clip_on=False)
    sns.despine(ax= axes[1], bottom = True, left = True)
    axes[1].set_xticks([])
    axes[1].set_xlabel('')
    axes[1].yaxis.label.set_size(25)
    axes[1].tick_params(labelsize=20)

    plt.savefig(output_path + '/bubble_data_' + evolved_name + '_' + naive_name + '_all_marks_correlations' + cell_type + '2.pdf', dpi=300, bbox_inches="tight")

barplot(evolved_all, di_local_all, "evolved", "naive", hparams.track_indices, hparams.cell_type, hparams.output_path)

scatter_matrix_marks(evolved_all, di_local_all, "evolved", "naive", hparams.track_indices, hparams.cell_type, hparams.output_path)

correlation_matrix(evolved_all, di_local_all, "evolved", "Di-Local", hparams.track_indices, hparams.cell_type, hparams.output_path)


barplot(evolved_all, di_local_all, "evolved", "Di-Local", [0,1,2,3,4], "iPS_DF_male_newborn")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [5,6,7,8,9], "foreskin_melanocyte_male_newborn")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [10,11,12,13,14], "neural_progenitor_cell_originated_from_H9")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [15,16,17,18,19], "mesenchymal_stem_cell_originated_from_H1")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [20,21,22,23,24], "iPS_DF_male_newborn_2")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [25,26,27,28,29], "hepatocyte_originated_from_H9")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [30,31,32,33,34], "foreskin_keratinocyte_male_newborn")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [35,36,37,38,39], "mesendoderm_originated_from_H1")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [45,46,47,48,49], "IMR-90")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [50,51,52,53,54], "foreskin_fibroblast_male_newborn")
barplot(evolved_all, di_local_all, "evolved", "Di-Local", [55,56,57,58,59], "neural stem progenitor cell originated from H1-hESC")
