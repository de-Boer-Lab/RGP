#June 15, 2023
#Fig 3B
#This script will make ecdfs for all nucleotide shuffled enformer predictions
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
                    help="Relative path of enformer predictions (.npy)", type=str, required=True)
parser.add_argument("--cell_type",
                    help="Cell type for experiment", type=str, required=True)
parser.add_argument("--cell_type_index",
                    help="Cell type index", type=int, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()

### Read in files
#only use first 1000 predictions to speed up analysis
random_all = np.load(hparams.path_to_predictions +'/random_sequences.npy')[:,:,0:1000]
genomic_all = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
mono_all = np.load(hparams.path_to_predictions+'/mononuc.npy')[:,:,0:1000]
tri_all = np.load(hparams.path_to_predictions+'/trinuc.npy')[:,:,0:1000]
di_all = np.load(hparams.path_to_predictions+'/dinuc.npy')[:,:,0:1000]
mono_local_all = np.load(hparams.path_to_predictions+'/mononuc_local.npy')[:,:,0:1000]
di_local_all = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]
tri_local_all = np.load(hparams.path_to_predictions+'/trinuc_local.npy')[:,:,0:1000]

#default indices used in enformer_multi.pbs match to these 12 cell types (all iPSC related)
def all_cell_types(output_path):
    cell_types_all = ['iPS_DF_6.9_male_newborn',
    'foreskin_melanocyte_male_newborn',
    'neural_progenitor_cell_originated_from_H9',
    'mesenchymal_stem_cell_originated_from_H1-hESC',
    'iPS_DF_19.11_male_newborn',
    'hepatocyte_originated_from_H9',
    'foreskin_melanocyte_male_newborn',
    'mesendoderm_originated_from_H1-hESC',
    'H1-hESC',
    'IMR-90',
    'foreskin_fibroblast_male_newborn',
    'neural_stem_progenitor_cell_originated_from_H1-hESC']
    index_counter = 0
    for i in range(0, len(cell_types_all)):

        cell_type = cell_types_all[i]
        print("index_counter")
        print(index_counter)
        random = random_all[:,index_counter:index_counter+5,:]
        genomic = genomic_all[:,index_counter:index_counter+5,:]
        mono = mono_all[:,index_counter:index_counter+5,:]
        tri = tri_all[:,index_counter:index_counter+5,:]
        di = di_all[:,index_counter:index_counter+5,:]
        mono_local = mono_local_all[:,index_counter:index_counter+5,:]
        di_local = di_local_all[:,index_counter:index_counter+5,:]
        tri_local = tri_local_all[:,index_counter:index_counter+5,:]

        for i in range(0,5):
            if i == 0:
                track = "DNase"
                v_min = -0.5
                v_max = 1
                c_min = 0.960
                c_max = 1.002
                cx_min = -0.01
                cx_max = 2
            if i == 1:
                track = "H3K27ac"
                v_min = -0.5
                v_max = 8
                c_min = 0.960
                c_max = 1.002
                cx_min = 2.0
                cx_max = 5
            if i == 2:
                track = "H3K4me1"
                v_min = -0.5
                v_max = 6
                c_min = 0.4
                c_max = 1.003
                cx_min = -0.01
                cx_max = 6
            if i == 3:
                track = "H3K4me3"
                v_min = -0.5
                v_max = 4
                c_min = 0.7
                c_max = 1.003
                cx_min = 0.8
                cx_max = 5
            if i == 4:
                track = "H3K27me3"
                v_min = -0.5
                v_max = 4
                c_min = 0.7
                c_max = 1.003
                cx_min = 1
                cx_max = 5

            track_list = pd.DataFrame(
                {'random_' + track: random[:, i, :].flatten('F'),
                 'genomic_' + track: genomic[:, i, :].flatten('F'),
                 'mono_'+ track: mono[:, i, :].flatten('F'),
                 'di_'+ track: di[:, i, :].flatten('F'),
                 'tri_'+ track: tri[:, i, :].flatten('F'),
                 'mono_local_'+ track: mono_local[:, i, :].flatten('F'),
                 'di_local_'+ track: di_local[:, i, :].flatten('F'),
                 'tri_local_'+ track: tri_local[:, i, :].flatten('F')
                })

            sns.set(font_scale = 1)
            sns.set_style(style='white')
            fig, ax = plt.subplots()
            #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
            green = (0/255, 102/255, 0/255)
            purple = (128/255, 0/255, 128/255)
            blue = (0/255, 0/255, 238/255)
            #overlay ecdfs for each sequence set prediction
            sns.ecdfplot(track_list['genomic_' + track], color = purple )
            sns.ecdfplot(track_list['random_'+ track], color = blue)
            sns.ecdfplot(track_list['mono_'+ track], color = (255/255,127/255,80/255))
            sns.ecdfplot(track_list['di_'+ track], color = (255/255,20/255,147/255))
            sns.ecdfplot(track_list['tri_'+ track], color = (0/255,238/255,238/255))
            sns.ecdfplot(track_list['mono_local_'+ track], color = (238/255,59/255,59/255))
            sns.ecdfplot(track_list['di_local_'+ track], color = green)
            sns.ecdfplot(track_list['tri_local_'+ track], color = (171/255,130/255,255/255))
            plt.legend(labels = ['genomic', 'random', 'mono', 'di', 'tri', 'mono_local','di_local', 'tri_local'], loc='lower right', fontsize=10)
            plt.title(cell_type)
            plt.ylim(0, 1.1)
            #ax.yaxis.grid(True)
            plt.xlabel('Enformer ' + track +  ' Prediction')
            ax.tick_params(bottom=True, left=True)
            plt.savefig(output_path + '/ecdfs_shuffle/eCDF_' + cell_type + '_' + track + '.pdf', dpi=400, bbox_inches='tight')

            #make y-axis zoomed in version of cdf
            sns.set(font_scale = 1)
            sns.set_style(style='white')
            fig, ax = plt.subplots()

            sns.ecdfplot(track_list['genomic_' + track], color = purple )
            sns.ecdfplot(track_list['random_'+ track], color = blue)
            sns.ecdfplot(track_list['mono_'+ track], color = (255/255,127/255,80/255))
            sns.ecdfplot(track_list['di_'+ track], color = (255/255,20/255,147/255))
            sns.ecdfplot(track_list['tri_'+ track], color = (0/255,238/255,238/255))
            sns.ecdfplot(track_list['mono_local_'+ track], color = (238/255,59/255,59/255))
            sns.ecdfplot(track_list['di_local_'+ track], color = green)
            sns.ecdfplot(track_list['tri_local_'+ track], color = (171/255,130/255,255/255))
            plt.legend(labels = ['genomic', 'random', 'mono', 'di', 'tri', 'mono_local', 'di_local', 'tri_local'], loc='lower right', fontsize=10)
            plt.ylim(c_min, c_max)
            plt.xlim(cx_min, cx_max)
            plt.title(cell_type)
            plt.xlabel('Enformer ' + track + ' Prediction')
            #ax.yaxis.grid(True)
            ax.tick_params(bottom=True, left=True)
            plt.savefig(output_path + '/eCDF_cropped_' + cell_type + '_' + track + '_.pdf', dpi=400, bbox_inches='tight')

        print("done all tracks for one cell type")
        index_counter = index_counter + 5


def one_cell_type(cell_type, index_start, output_path):
    random = random_all[:,index_start:index_start+5,:]
    genomic = genomic_all[:,index_start:index_start+5,:]
    mono = mono_all[:,index_start:index_start+5,:]
    tri = tri_all[:,index_start:index_start+5,:]
    di = di_all[:,index_start:index_start+5,:]
    mono_local = mono_local_all[:,index_start:index_start+5,:]
    di_local = di_local_all[:,index_start:index_start+5,:]
    tri_local = tri_local_all[:,index_start:index_start+5,:]

    for i in range(0,5):
        if i == 0:
            track = "DNase"
            v_min = -0.5
            v_max = 1
            c_min = 0.960
            c_max = 1.002
            cx_min = -0.01
            cx_max = 2
        if i == 1:
            track = "H3K27ac"
            v_min = -0.5
            v_max = 8
            c_min = 0.960
            c_max = 1.002
            cx_min = 2.0
            cx_max = 5
        if i == 2:
            track = "H3K4me1"
            v_min = -0.5
            v_max = 6
            c_min = 0.4
            c_max = 1.003
            cx_min = -0.01
            cx_max = 6
        if i == 3:
            track = "H3K4me3"
            v_min = -0.5
            v_max = 4
            c_min = 0.7
            c_max = 1.003
            cx_min = 0.8
            cx_max = 5
        if i == 4:
            track = "H3K27me3"
            v_min = -0.5
            v_max = 4
            c_min = 0.7
            c_max = 1.003
            cx_min = 1
            cx_max = 5


        track_list = pd.DataFrame(
            {'random_' + track: random[:, i, :].flatten('F'),
             'genomic_' + track: genomic[:, i, :].flatten('F'),
             'mono_'+ track: mono[:, i, :].flatten('F'),
             'di_'+ track: di[:, i, :].flatten('F'),
             'tri_'+ track: tri[:, i, :].flatten('F'),
             'mono_local_'+ track: mono_local[:, i, :].flatten('F'),
             'di_local_'+ track: di_local[:, i, :].flatten('F'),
             'tri_local_'+ track: tri_local[:, i, :].flatten('F')
            })

        sns.set(font_scale = 1)
        sns.set_style(style='white')
        fig, ax = plt.subplots()
        #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
        green = (0/255, 102/255, 0/255)
        purple = (128/255, 0/255, 128/255)
        blue = (0/255, 0/255, 238/255)

        sns.ecdfplot(track_list['genomic_' + track], color = purple )
        sns.ecdfplot(track_list['random_'+ track], color = blue)
        sns.ecdfplot(track_list['mono_'+ track], color = '#bb486a')
        sns.ecdfplot(track_list['di_'+ track], color = green)
        sns.ecdfplot(track_list['tri_'+ track], color = '#ff9c79')
        sns.ecdfplot(track_list['mono_local_'+ track], color = '#bb486a', linestyle='--')
        sns.ecdfplot(track_list['di_local_'+ track], color = green, linestyle='--')
        sns.ecdfplot(track_list['tri_local_'+ track], color = '#ff9c79', linestyle='--')
        plt.legend(labels = ['genomic', 'random sampled', 'mononucleotide', 'dinucleotide', 'trinucleotide', 'mononucleotide', 'dinucleotide', 'trinucleotide'], loc='lower right', fontsize=10)
        plt.title(cell_type)
        plt.ylim(0, 1.1)
        #ax.yaxis.grid(True)
        plt.xlabel('Enformer ' + track +  ' Prediction')
        ax.tick_params(bottom=True, left=True)
        plt.savefig(output_path + '/eCDF_' + cell_type + '_' + track + '.pdf', dpi=400, bbox_inches='tight')

        #make cropped version of cdf
        sns.set(font_scale = 2)
        sns.set_style(style='white')
        fig, ax = plt.subplots()
        #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))

        sns.ecdfplot(track_list['genomic_' + track], color = purple )
        sns.ecdfplot(track_list['random_'+ track], color = blue)
        sns.ecdfplot(track_list['mono_'+ track], color = '#bb486a')
        sns.ecdfplot(track_list['di_'+ track], color = green)
        sns.ecdfplot(track_list['tri_'+ track], color = '#ff9c79')
        sns.ecdfplot(track_list['mono_local_'+ track], color = '#bb486a', linestyle='--')
        sns.ecdfplot(track_list['di_local_'+ track], color = green, linestyle='--')
        sns.ecdfplot(track_list['tri_local_'+ track], color = '#ff9c79', linestyle='--')
        plt.legend(labels = ['genomic sequence', 'random sampled', 'mononucleotide', 'dinucleotide', 'trinucleotide', 'mononucleotide', 'dinucleotide', 'trinucleotide'], loc='lower right', fontsize=10)
        plt.ylim(c_min, c_max)
        plt.xlim(cx_min, cx_max)
        #plt.title(cell_type)
        plt.xlabel('Enformer ' + cell_type + ' ' + track + ' Prediction')
        #ax.yaxis.grid(True)
        ax.tick_params(bottom=True, left=True)
        plt.savefig(output_path + '/eCDF_cropped_' + cell_type + '_' + track + '_1.pdf', dpi=400, bbox_inches='tight')

#to create ecdfs for 1 cell type
one_cell_type(hparams.cell_type, hparams.cell_type_index, hparams.output_path)
#to create ecdfs for all cell types
all_cell_types(output_path)
