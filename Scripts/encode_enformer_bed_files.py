#June 15, 2023
#This script makes the bed files used to overlap the enformer predicted and ENCODE predicted peaks
#Need to run this before Evolution_plots_final.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
from scipy import stats
from PIL import Image as im
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.mlab as mlab

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
parser.add_argument("--encode_file_paths",
                    help="path to where ENCODE bed files", type=str, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)

hparams, _ = parser.parse_known_args()
os.mkdir(hparams.output_path)

### Read in files
evolved = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
naive = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]

tracks = ['Dnase', 'H3K27ac','H3K4me1','H3K4me3', 'H3K27me3']

def make_bed_files(evolved_all, input_bed_encode, track, cell_type, track_list, output_path):
    #Clip the ENCODE file to just chr start and end
    all_bed_split_128 = pd.read_csv("/project/st-cdeboer-1/iluthra//enformer_random_DNA/enformer_data/ENCODE/genomic_regions_coordinates_n1000_114kb_128bins.bed", sep = '\t')

    if track == "Dnase":
        index = track_list[0]
    if track == "H3K27ac":
        index = track_list[1]
    if track == "H3K4me1":
        index = track_list[2]
    if track == "H3K4me3":
        index = track_list[3]
    if track == "H3K27me3":
        index = track_list[4]

    enformer_predictions = evolved[:,index,:].flatten('F')
    all_bed_split_128[track] = enformer_predictions

    all_bed_split_128.to_csv(output_path + "/genomic_regions_coordinates_n1000_114kb_128bins_" + cell_type + '_' + track + ".bed", index = False, sep = '\t', header = None)

    #make the ENCODE peak file 3 columns
    bed_file_H1 = pd.read_csv("/project/st-cdeboer-1/iluthra//enformer_random_DNA/enformer_data/ENCODE/" + input_bed_encode, sep = '\t', header = None)
    bed_file_coor = bed_file_H1.iloc[:,0:3]
    bed_file_coor.to_csv(output_path + cell_type + '_' + track + "encode.bed", index = False, header = None, sep = '\t')

#need to download corresponding bed files from ENCODE manually
make_bed_files(evolved, "ENCFF983UCL.bed", "Dnase", hparams.cell_type, hparams.track_indices, hparams.output_path)
make_bed_files(evolved, "ENCFF520LUH.bed", "H3K27ac",hparams.cell_type, hparams.track_indices, hparams.output_path)
make_bed_files(evolved, "ENCFF711NNX.bed", "H3K4me1", hparams.cell_type,hparams.track_indices, hparams.output_path)
make_bed_files(evolved, "ENCFF760ZNE.bed", "H3K4me3", hparams.cell_type,hparams.track_indices, hparams.output_path)
make_bed_files(evolved, "ENCFF760MIC.bed", "H3K27me3",hparams.cell_type, hparams.track_indices, hparams.output_path)

# run bedtools with this file and enformer test region file
#bedtools intersect -wao -a genomic_regions_coordinates_n1000_114kb_128bins_h1_dnase.bed -b H1-hESC_track.bed > track_enformer_overlap.bed
