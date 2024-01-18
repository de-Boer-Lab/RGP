#June 15, 2023
#Fig 4C
#This script will make a scatter plot for 3 germ layers of your choice
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
parser.add_argument("--track_indices",
                    help="tracks to use should be ordered as [mesoderm endoderm ectoderm]", nargs = '+', type=int, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()
os.mkdir(hparams.output_path)
### Read in files
random_all = np.load(hparams.path_to_predictions +'/random_sequences.npy')[:,:,0:1000]
evolved_all = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
di_local_all = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]

def corrfunc(x, y, **kws):
    (r, p) = stats.pearsonr(x, y)
    label = 'r = ' + str(round(r, 3))
    ax = plt.gca()
    ax.annotate(label,
                xy=(.4, 1.05), xycoords=ax.transAxes)

#single_scatter(evolved_all[:,25,:], di_local_all[:, 25, :], 'Endoderm - H9', 'Endoderm - H9 di_local')

def single_scatter(evolved, naive, evolved_name, naive_name):
    current_track = pd.DataFrame(
        {evolved_name: evolved.flatten(),
         naive_name: naive.flatten()
        })
    fig, axes = plt.subplots(figsize=(10, 10))
    sns.scatterplot(data=current_track, x=evolved_name, y=naive_name, alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized = True)
    axes.xaxis.label.set_size(16)
    axes.yaxis.label.set_size(16)
    axes.tick_params(labelsize=18)
    (r, p) = stats.pearsonr(current_track[evolved_name], current_track[naive_name])
    label = 'r = ' + str(round(r, 2))
    axes.annotate(label, xy=(.9, 0.95), xycoords=axes.transAxes, fontsize=18)
    plt.savefig(output_path + '/scatter_' + naive_name + '_' + evolved_name + '.pdf', dpi=400, bbox_inches="tight")


def scatter_matrix(evolved, naive, evolved_name, naive_name, track_list, output_path):

    meso_endo_ecto = pd.DataFrame(
        {'Ectoderm evolved': evolved[:,track_list[0],:].flatten('F'),
         'Mesoderm evolved': evolved[:, track_list[1], :].flatten('F'),
         'Endoderm evolved': evolved[:, track_list[2], :].flatten('F'),
         'Ectoderm naive': naive[:, track_list[0], :].flatten('F'),
         'Mesoderm naive': naive[:, track_list[1], :].flatten('F'),
         'Endoderm naive': naive[:, track_list[2], :].flatten('F')
        })
    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)

    print(meso_endo_ecto)
    meso_endo_ecto.to_csv(output_path + '/meso_endo_ecto.csv')
    #sns.set(font_scale = 1.5)
    fig, axes = plt.subplots(3, 3, figsize=(30, 30))
    sns.scatterplot(ax=axes[0, 0], data=meso_endo_ecto, x='Mesoderm evolved', y='Mesoderm naive', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[0,0].xaxis.label.set_size(35)
    axes[0,0].yaxis.label.set_size(35)
    axes[0,0].tick_params(labelsize=30)
    axes[0,0].xaxis.label.set_color(purple)
    axes[0,0].yaxis.label.set_color(green)
    (r, p) = stats.pearsonr(meso_endo_ecto['Mesoderm evolved'], meso_endo_ecto['Mesoderm naive'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[0,0].annotate(label, xy=(.6, 0.9), xycoords=axes[0,0].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[0,0].annotate(label_p, xy=(.6, 0.8), xycoords=axes[0,0].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[1, 0], data=meso_endo_ecto, x='Mesoderm evolved', y='Endoderm evolved', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[1,0].xaxis.label.set_size(35)
    axes[1,0].yaxis.label.set_size(35)
    axes[1,0].tick_params(labelsize=30)
    axes[1,0].xaxis.label.set_color(purple)
    axes[1,0].yaxis.label.set_color(purple)
    (r, p) = stats.pearsonr(meso_endo_ecto['Mesoderm evolved'], meso_endo_ecto['Endoderm evolved'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[1,0].annotate(label, xy=(.6, 0.9), xycoords=axes[1,0].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[1,0].annotate(label_p, xy=(.6, 0.8), xycoords=axes[1,0].transAxes, fontsize=35)


    sns.scatterplot(ax=axes[2, 0], data=meso_endo_ecto, x='Mesoderm evolved', y='Ectoderm evolved', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[2,0].xaxis.label.set_size(35)
    axes[2,0].yaxis.label.set_size(35)
    axes[2,0].tick_params(labelsize=30)
    axes[2,0].xaxis.label.set_color(purple)
    axes[2,0].yaxis.label.set_color(purple)
    (r, p) = stats.pearsonr(meso_endo_ecto['Mesoderm evolved'], meso_endo_ecto['Ectoderm evolved'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[2,0].annotate(label, xy=(.6, 0.9), xycoords=axes[2,0].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[2,0].annotate(label_p, xy=(.6, 0.8), xycoords=axes[2,0].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[0, 1], data=meso_endo_ecto, x='Endoderm naive', y='Mesoderm naive', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[0,1].xaxis.label.set_size(35)
    axes[0,1].yaxis.label.set_size(35)
    axes[0,1].tick_params(labelsize=30)
    axes[0,1].xaxis.label.set_color(green)
    axes[0,1].yaxis.label.set_color(green)
    (r, p) = stats.pearsonr(meso_endo_ecto['Endoderm naive'], meso_endo_ecto['Mesoderm naive'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[0,1].annotate(label, xy=(.6, 0.9), xycoords=axes[0,1].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[0,1].annotate(label_p, xy=(.6, 0.8), xycoords=axes[0,1].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[1, 1], data=meso_endo_ecto, x='Endoderm evolved', y='Endoderm naive', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[1,1].xaxis.label.set_size(35)
    axes[1,1].yaxis.label.set_size(35)
    axes[1,1].tick_params(labelsize=30)
    axes[1,1].xaxis.label.set_color(purple)
    axes[1,1].yaxis.label.set_color(green)
    (r, p) = stats.pearsonr(meso_endo_ecto['Endoderm evolved'], meso_endo_ecto['Endoderm naive'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[1,1].annotate(label, xy=(.6, 0.9), xycoords=axes[1,1].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[1,1].annotate(label_p, xy=(.6, 0.8), xycoords=axes[1,1].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[2, 1], data=meso_endo_ecto, x='Endoderm evolved', y='Ectoderm evolved', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[2,1].xaxis.label.set_size(35)
    axes[2,1].yaxis.label.set_size(35)
    axes[2,1].tick_params(labelsize=30)
    axes[2,1].xaxis.label.set_color(purple)
    axes[2,1].yaxis.label.set_color(purple)
    (r, p) = stats.pearsonr(meso_endo_ecto['Endoderm evolved'], meso_endo_ecto['Ectoderm evolved'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[2,1].annotate(label, xy=(.6, 0.9), xycoords=axes[2,1].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[2,1].annotate(label_p, xy=(.6, 0.8), xycoords=axes[2,1].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[0, 2], data=meso_endo_ecto, x='Ectoderm naive', y='Mesoderm naive', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[0,2].xaxis.label.set_size(35)
    axes[0,2].yaxis.label.set_size(35)
    axes[0,2].tick_params(labelsize=30)
    axes[0,2].xaxis.label.set_color(green)
    axes[0,2].yaxis.label.set_color(green)
    (r, p) = stats.pearsonr(meso_endo_ecto['Ectoderm naive'], meso_endo_ecto['Mesoderm naive'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[0,2].annotate(label, xy=(.6, 0.9), xycoords=axes[0,2].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[0,2].annotate(label_p, xy=(.6, 0.8), xycoords=axes[0,2].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[1, 2], data=meso_endo_ecto, x='Ectoderm naive', y='Endoderm naive', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[1,2].xaxis.label.set_size(35)
    axes[1,2].yaxis.label.set_size(35)
    axes[1,2].tick_params(labelsize=30)
    axes[1,2].xaxis.label.set_color(green)
    axes[1,2].yaxis.label.set_color(green)
    (r, p) = stats.pearsonr(meso_endo_ecto['Ectoderm naive'], meso_endo_ecto['Endoderm naive'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[1,2].annotate(label, xy=(.6, 0.9), xycoords=axes[1,2].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print(label_p)
    axes[1,2].annotate(label_p, xy=(.6, 0.8), xycoords=axes[1,2].transAxes, fontsize=35)

    sns.scatterplot(ax=axes[2, 2], data=meso_endo_ecto, x='Ectoderm evolved', y='Ectoderm naive', alpha=.1, s = 100, color = '#000000', linewidth=0, rasterized=True)
    axes[2,2].xaxis.label.set_size(35)
    axes[2,2].yaxis.label.set_size(35)
    axes[2,2].tick_params(labelsize=30)
    axes[2,2].xaxis.label.set_color(purple)
    axes[2,2].yaxis.label.set_color(green)
    (r, p) = stats.pearsonr(meso_endo_ecto['Ectoderm evolved'], meso_endo_ecto['Ectoderm naive'])
    label = '$\it{r}$ = ' + str(round(r, 3))
    axes[2,2].annotate(label, xy=(.6, 0.9), xycoords=axes[2,2].transAxes, fontsize=35)
    label_p = '$\it{p}$ = ' + str(p)
    print("p-value")
    print(stats.pearsonr(meso_endo_ecto['Ectoderm evolved'], meso_endo_ecto['Ectoderm naive']))
    axes[2,2].annotate(label_p, xy=(.6, 0.8), xycoords=axes[2,2].transAxes, fontsize=35)

    plt.savefig(output_path + '/scatter_' + evolved_name + '_' + naive_name + '_meso_endo_ecto' + '.pdf', dpi=300, bbox_inches="tight")


scatter_matrix(evolved_all, di_local_all, "evolved", "naive", hparams.track_indices, hparams.output_path)

#For creating single scatter plots
single_scatter(evolved_all[:,40,:], evolved_all[:,41,:], "evolved-H1-DNase", "evolved-H1-H3K27ac")
single_scatter(di_local_all[:,40,:], di_local_all[:,41,:], "DiLocal-H1-DNase", "DiLocal-H1-H3K27ac")
single_scatter(random_all[:,40,:], random_all[:,41,:], "Random-H1-DNase", "Random-H1-H3K27ac")

single_scatter(evolved_all[:,40,:], evolved_all[:,42,:], "evolved-H1-DNase", "evolved-H1-H3K4me1")
single_scatter(di_local_all[:,40,:], di_local_all[:,42,:], "DiLocal-H1-DNase", "DiLocal-H1-H3K4me1")
single_scatter(random_all[:,40,:], random_all[:,42,:], "Random-H1-DNase", "Random-H1-H3K4me1")

single_scatter(evolved_all[:,40,:], evolved_all[:,43,:], "evolved-H1-DNase", "evolved-H1-H3K4me3")
single_scatter(di_local_all[:,40,:], di_local_all[:,43,:], "DiLocal-H1-DNase", "DiLocal-H1-H3K4me3")
single_scatter(random_all[:,40,:], random_all[:,43,:], "Random-H1-DNase", "Random-H1-H3K4me3")

single_scatter(evolved_all[:,40,:], evolved_all[:,44,:], "evolved-H1-DNase", "evolved-H1-H3K27me3")
single_scatter(di_local_all[:,40,:], di_local_all[:,44,:], "DiLocal-H1-DNase", "DiLocal-H1-H3K27me3")
single_scatter(random_all[:,40,:], random_all[:,44,:], "Random-H1-DNase", "Random-H1-H3K27me3")
