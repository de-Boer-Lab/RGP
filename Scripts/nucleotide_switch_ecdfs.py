#June 15, 2023
#Fig 3D
#This script will make ecdfs for all nucleotide switch enformer predictions
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
os.mkdir(hparams.output_path)
### Read in files

random_all = np.load(hparams.path_to_predictions+'/random_sequences.npy')[:,:,0:1000]
genomic_all = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
AG2GA_all = np.load(hparams.path_to_predictions+'/AG2GA.npy')[:,:,0:1000]
AC2CA_all = np.load(hparams.path_to_predictions+'/AC2CA.npy')[:,:,0:1000]
AT2TA_all = np.load(hparams.path_to_predictions+'/AT2TA.npy')[:,:,0:1000]
CA2AC_all = np.load(hparams.path_to_predictions+'/CA2AC.npy')[:,:,0:1000]
CG2GC_all = np.load(hparams.path_to_predictions+'/CG2GC.npy')[:,:,0:1000]
GA2AG_all = np.load(hparams.path_to_predictions+'/GA2AG.npy')[:,:,0:1000]
GC2CG_all = np.load(hparams.path_to_predictions+'/GC2CG.npy')[:,:,0:1000]
TA2AT_all = np.load(hparams.path_to_predictions+'/TA2AT.npy')[:,:,0:1000]


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

        random = random_all[:,index_counter:index_counter+5,:]
        genomic = genomic_all[:,index_counter:index_counter+5,:]
        AG2GA = AG2GA_all[:,index_counter:index_counter+5,:]
        AC2CA = AC2CA_all[:,index_counter:index_counter+5,:]
        AT2TA = AT2TA_all[:,index_counter:index_counter+5,:]
        CA2AC = CA2AC_all[:,index_counter:index_counter+5,:]
        CG2GC = CG2GC_all[:,index_counter:index_counter+5,:]
        GA2AG = GA2AG_all[:,index_counter:index_counter+5,:]
        GC2CG = GC2CG_all[:,index_counter:index_counter+5,:]
        TA2AT = TA2AT_all[:,index_counter:index_counter+5,:]

        for i in range(0,5):
            print(i)
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
                cx_min = 1
                cx_max = 5
            if i == 2:
                track = "H3K4me1"
                v_min = -0.5
                v_max = 6
                c_min = 0.7
                c_max = 1.003
                cx_min = -0.01
                cx_max = 6
            if i == 3:
                track = "H3K4me3"
                v_min = -0.5
                v_max = 4
                c_min = 0.7
                c_max = 1.003
                cx_min = 0.5
                cx_max = 5
            if i == 4:
                track = "H3K27me3"
                v_min = -0.5
                v_max = 4
                c_min = 0.7
                c_max = 1.003
                cx_min = 0.3
                cx_max = 5


            track_list = pd.DataFrame(
                {'AG2GA_' + track: AG2GA[:, i, :].flatten(),
                 'AC2CA_' + track: AC2CA[:, i, :].flatten(),
                 'AT2TA_' + track: AT2TA[:, i, :].flatten(),
                 'CA2AC_'+ track: CA2AC[:, i, :].flatten(),
                 'CG2GC_'+ track: CG2GC[:, i, :].flatten(),
                 'GA2AG_'+ track: GA2AG[:, i, :].flatten(),
                 'GC2CG_'+ track: GC2CG[:, i, :].flatten(),
                 'TA2AT_'+ track: TA2AT[:, i, :].flatten(),
                 'random_' + track: random[:, i, :].flatten(),
                 'genomic_' + track: genomic[:, i, :].flatten()
                })

            green = (0/255, 102/255, 0/255)
            purple = (128/255, 0/255, 128/255)
            blue = (0/255, 0/255, 238/255)

            sns.set(font_scale = 1.0)
            sns.set_style(style='white')
            fig, ax = plt.subplots()
            #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
            sns.ecdfplot(track_list['genomic_' + track], color = purple)
            sns.ecdfplot(track_list['random_' + track], color = blue)
            sns.ecdfplot(track_list['AG2GA_'+ track], color = (255/255,127/255,80/255))
            sns.ecdfplot(track_list['AC2CA_'+ track], color = (255/255,20/255,147/255))
            sns.ecdfplot(track_list['AT2TA_'+ track], color = (0/255,238/255,238/255))
            sns.ecdfplot(track_list['CA2AC_'+ track], color = (238/255,59/255,59/255))
            sns.ecdfplot(track_list['CG2GC_'+ track], color = (171/255,130/255,255/255))
            sns.ecdfplot(track_list['GA2AG_'+ track], color = (0/255,205/255,0/255))
            sns.ecdfplot(track_list['GC2CG_'+ track], color = (131/255,139/255,131/255))
            sns.ecdfplot(track_list['TA2AT_'+ track], color = (255/255,215/255,0/255))
            plt.legend(labels = ['genomic', 'random', 'AG' + r'$\rightarrow$' + 'GA', 'AC' + r'$\rightarrow$' + 'CA','AT' + r'$\rightarrow$' + 'TA', 'CAr' + r'$\rightarrow$' + 'AC', 'CG' + r'$\rightarrow$' + 'GC','GA' + r'$\rightarrow$' + 'AG', 'GC' + r'$\rightarrow$' + 'CG', 'TA' + r'$\rightarrow$' + 'AT'], loc='lower right', fontsize=10)
            plt.title(cell_type)
            plt.xlabel('Enformer ' + track +  ' Prediction')
            plt.ylim(0, 1.1)
            #ax.yaxis.grid(True)
            ax.tick_params(bottom=True, left=True)
            plt.savefig(output_path + '/eCDF_' + cell_type + '_' + track + '.pdf', dpi=400, bbox_inches='tight')

            #make cropped version of cdf
            sns.set(font_scale = 1.0)
            sns.set_style(style='white')
            fig, ax = plt.subplots()
            #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
            sns.ecdfplot(track_list['genomic_' + track], color = purple)
            sns.ecdfplot(track_list['random_' + track], color = blue)
            sns.ecdfplot(track_list['AG2GA_'+ track], color = (255/255,127/255,80/255))
            sns.ecdfplot(track_list['AC2CA_'+ track], color = (255/255,20/255,147/255))
            sns.ecdfplot(track_list['AT2TA_'+ track], color = (0/255,238/255,238/255))
            sns.ecdfplot(track_list['CA2AC_'+ track], color = (238/255,59/255,59/255))
            sns.ecdfplot(track_list['CG2GC_'+ track], color = (171/255,130/255,255/255))
            sns.ecdfplot(track_list['GA2AG_'+ track], color = (0/255,205/255,0/255))
            sns.ecdfplot(track_list['GC2CG_'+ track], color = (131/255,139/255,131/255))
            sns.ecdfplot(track_list['TA2AT_'+ track], color = (255/255,215/255,0/255))
            plt.legend(labels = ['genomic', 'random', 'AG' + r'$\rightarrow$' + 'GA', 'AC' + r'$\rightarrow$' + 'CA','AT' + r'$\rightarrow$' + 'TA', 'CAr' + r'$\rightarrow$' + 'AC', 'CG' + r'$\rightarrow$' + 'GC','GA' + r'$\rightarrow$' + 'AG', 'GC' + r'$\rightarrow$' + 'CG', 'TA' + r'$\rightarrow$' + 'AT'], loc='lower right', fontsize=10)
            plt.ylim(c_min, c_max)
            plt.xlim(cx_min, cx_max)
            plt.xlabel('Enformer ' + track +  ' Prediction')
            plt.title(cell_type)
            #ax.yaxis.grid(True)
            ax.tick_params(bottom=True, left=True)
            plt.savefig(output_path + '/eCDF' + cell_type + '_' + track + '_cropped.pdf', dpi=400, bbox_inches='tight')

            print("done")
        print("done all tracks for one cell type")
        index_counter = index_counter + 5

def one_cell_type(cell_type, index_counter, output_path):
    random = random_all[:,index_counter:index_counter+5,:]
    genomic = genomic_all[:,index_counter:index_counter+5,:]
    AG2GA = AG2GA_all[:,index_counter:index_counter+5,:]
    AC2CA = AC2CA_all[:,index_counter:index_counter+5,:]
    AT2TA = AT2TA_all[:,index_counter:index_counter+5,:]
    CA2AC = CA2AC_all[:,index_counter:index_counter+5,:]
    CG2GC = CG2GC_all[:,index_counter:index_counter+5,:]
    GA2AG = GA2AG_all[:,index_counter:index_counter+5,:]
    GC2CG = GC2CG_all[:,index_counter:index_counter+5,:]
    TA2AT = TA2AT_all[:,index_counter:index_counter+5,:]

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
            {'AG2GA_' + track: AG2GA[:, i, :].flatten(),
             'AC2CA_' + track: AC2CA[:, i, :].flatten(),
             'AT2TA_' + track: AT2TA[:, i, :].flatten(),
             'CA2AC_'+ track: CA2AC[:, i, :].flatten(),
             'CG2GC_'+ track: CG2GC[:, i, :].flatten(),
             'GA2AG_'+ track: GA2AG[:, i, :].flatten(),
             'GC2CG_'+ track: GC2CG[:, i, :].flatten(),
             'TA2AT_'+ track: TA2AT[:, i, :].flatten(),
             'random_' + track: random[:, i, :].flatten(),
             'genomic_' + track: genomic[:, i, :].flatten()
            })

        sns.set(font_scale = 1)
        sns.set_style(style='white')
        fig, ax = plt.subplots()
        #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
        green = (0/255, 102/255, 0/255)
        purple = (128/255, 0/255, 128/255)
        blue = (0/255, 0/255, 238/255)

        sns.set(font_scale = 1.0)
        sns.set_style(style='white')
        fig, ax = plt.subplots()
        #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
        sns.ecdfplot(track_list['genomic_' + track], color = purple)
        sns.ecdfplot(track_list['random_' + track], color = blue)
        sns.ecdfplot(track_list['AG2GA_'+ track], color = '#ffb46d')
        sns.ecdfplot(track_list['AC2CA_'+ track], color = '#c15083')
        sns.ecdfplot(track_list['AT2TA_'+ track], color = '#87ba5b')
        sns.ecdfplot(track_list['CG2GC_'+ track], color = '#8785bb')
        sns.ecdfplot(track_list['GA2AG_'+ track], color = '#ffb46d', linestyle='--')
        sns.ecdfplot(track_list['CA2AC_'+ track], color = '#c15083',linestyle='--')
        sns.ecdfplot(track_list['TA2AT_'+ track], color = '#87ba5b', linestyle='--')
        sns.ecdfplot(track_list['GC2CG_'+ track], color = '#8785bb', linestyle='--')
        plt.legend(labels = ['genomic sequence', 'random sampled', 'AG' + r'$\rightarrow$' + 'GA', 'AC' + r'$\rightarrow$' + 'CA','AT' + r'$\rightarrow$' + 'TA', 'CG' + r'$\rightarrow$' + 'GC', 'GA' + r'$\rightarrow$' + 'AG', 'CA' + r'$\rightarrow$' + 'AC', 'TA' + r'$\rightarrow$' + 'AT', 'GC' + r'$\rightarrow$' + 'CG'], loc='lower right', fontsize=10)
        plt.ylim(c_min, c_max)
        plt.xlim(cx_min, cx_max)
        plt.title(cell_type)
        plt.xlabel('Enformer ' + track +  ' Prediction')
        plt.ylim(0, 1.1)
        #ax.yaxis.grid(True)
        ax.tick_params(bottom=True, left=True)
        plt.savefig(output_path + '/eCDF_' + cell_type + '_' + track + '.pdf', dpi=400, bbox_inches='tight')

        #make cropped version of cdf
        sns.set(font_scale = 2)
        sns.set_style(style='white')
        fig, ax = plt.subplots()
        #sns.ecdfplot(x="variable", hue="value", data=pd.melt((track_list)))
        sns.ecdfplot(track_list['genomic_' + track], color = purple)
        sns.ecdfplot(track_list['random_' + track], color = blue)
        sns.ecdfplot(track_list['AG2GA_'+ track], color = '#ffb46d')
        sns.ecdfplot(track_list['AC2CA_'+ track], color = '#c15083')
        sns.ecdfplot(track_list['AT2TA_'+ track], color = '#87ba5b')
        sns.ecdfplot(track_list['CG2GC_'+ track], color = '#8785bb')
        sns.ecdfplot(track_list['GA2AG_'+ track], color = '#ffb46d', linestyle='--')
        sns.ecdfplot(track_list['CA2AC_'+ track], color = '#c15083',linestyle='--')
        sns.ecdfplot(track_list['TA2AT_'+ track], color = '#87ba5b', linestyle='--')
        sns.ecdfplot(track_list['GC2CG_'+ track], color = '#8785bb', linestyle='--')
        plt.legend(labels = ['genomic sequence', 'random sampled', 'AG' + r'$\rightarrow$' + 'GA', 'AC' + r'$\rightarrow$' + 'CA','AT' + r'$\rightarrow$' + 'TA', 'CG' + r'$\rightarrow$' + 'GC', 'GA' + r'$\rightarrow$' + 'AG', 'CA' + r'$\rightarrow$' + 'AC', 'TA' + r'$\rightarrow$' + 'AT', 'GC' + r'$\rightarrow$' + 'CG'], loc='lower right', fontsize=10)
        plt.ylim(c_min, c_max)
        plt.xlim(cx_min, cx_max)
        plt.xlabel('Enformer ' + cell_type + ' ' + track + ' Prediction')
        #plt.title(cell_type)
        #ax.yaxis.grid(True)
        ax.tick_params(bottom=True, left=True)
        plt.savefig(output_path + '/eCDF' + cell_type + '_' + track + '_cropped.pdf', dpi=400, bbox_inches='tight')

#to create ecdfs for 1 cell type
one_cell_type(hparams.cell_type, hparams.cell_type_index, hparams.output_path)
#to create ecdfs for all cell types
all_cell_types(hparams.output_path)
