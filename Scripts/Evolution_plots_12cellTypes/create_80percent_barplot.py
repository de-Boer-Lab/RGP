#Fig 4B and and Supplementary Figure 6
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
from scipy import stats
all_cell_80 = pd.read_csv("/scratch/st-cdeboer-1/iluthra/randomDNA/Analysis_outputs_06152023/H1_DNase_rep1_07052023/Evolution_plots_12cellTypes/all_80_percent_output.txt", sep = '\t', header = None)
encode_qc = pd.read_csv("/project/st-cdeboer-1/iluthra/enformer/test_enformer/de-former/Paper_Review_Analysis/Final_clean_scripts_figures/Evolution_plots_12cellTypes/RGP_cell_types_annotation_file_CJ_IL.csv", header = None, delimiter=r"\s+")


all_cell_80.columns = ['cell_type', 'track', '80']
print(encode_qc)
encode_qc.columns = ['track','cell_type', 'track_index', 'index','exp','encode','qc']
merged = pd.merge(all_cell_80, encode_qc, on=['cell_type','track'])
dnase = merged[merged.track == "DNASE"]
print(dnase['80'].median())

H3K27ac = merged[merged.track == "H3K27ac"]
print(H3K27ac['80'].median())

H3K27me3 = merged[merged.track == "H3K27me3"]
print(H3K27me3['80'].median())

H3K4me1 = merged[merged.track == "H3K4me1"]
print(H3K4me1['80'].median())

H3K4me3 = merged[merged.track == "H3K4me3"]
print(H3K4me3['80'].median())
fig, ax = plt.subplots(figsize=(10,10))
sns.barplot(data=merged, x="track", y="80", hue="cell_type")
plt.xlabel('Enformer Chromatin Track')
plt.ylabel('Relative enrichment of naive sequence activity')
plt.savefig('/project/st-cdeboer-1/iluthra/enformer/test_enformer/de-former/Paper_Review_Analysis/Final_clean_scripts_figures/Evolution_plots_12cellTypes/80_percent_barplot.pdf', dpi=400, bbox_inches='tight')


fig, ax = plt.subplots(figsize=(10,10))
ax = sns.boxplot(data=merged, x="track", y="80", boxprops={'alpha': 0.4}, fliersize = 0)
sns.stripplot(data=merged, x="track", y="80", dodge=True, ax=ax)
plt.xlabel('Enformer Chromatin Track')
plt.ylabel('Relative enrichment of naive sequence activity')
plt.savefig('/project/st-cdeboer-1/iluthra/enformer/test_enformer/de-former/Paper_Review_Analysis/Final_clean_scripts_figures/Evolution_plots_12cellTypes/80_percent_boxplot.pdf', dpi=400, bbox_inches='tight')
