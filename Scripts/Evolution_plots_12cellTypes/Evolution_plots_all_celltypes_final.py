#June 15, 2023
#Fig 4A, Supplementary Figure 11,12,13

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

parser = argparse.ArgumentParser()

parser.add_argument("--path_to_predictions",
                    help="Relative path of enformer predictions (.npy)", type=str, required=True)
parser.add_argument("--cell_type",
                    help="Cell type for experiment", type=str, required=True)
parser.add_argument("--track_index",
                    help="track index", nargs = '+', type=int, required=True)
parser.add_argument("--track_name",
                    help="track type", type=str, required=True)
parser.add_argument("--output_path",
                    help="path to output plots", type=str, required=True)
parser.add_argument("--bed_file",
                    help="path to bedfile", type=str, required=True)
hparams, _ = parser.parse_known_args()

def intersect_distribution_panel(intersect_bed, cell_type,track_name, index_x, index_y, output_path):
    encode_enformer_merged = pd.read_csv(output_path + intersect_bed, sep = '\t', header = None)

    encode_enformer_merged.columns = ['binchr', 'binstart','binend','enformer_track_value', 'peakchr','peakstart','peakend','bp_overlap']
    #sort by overlap by biggest to smallest
    encode_enformer_merged = encode_enformer_merged.sort_values(by='bp_overlap', ascending=False)
    encode_enformer_merged = encode_enformer_merged.drop_duplicates(['binchr', 'binstart','binend'], keep='first')

    in_peak = pd.DataFrame(encode_enformer_merged[encode_enformer_merged.iloc[:,7] >= 10])
    not_in = pd.DataFrame(encode_enformer_merged[encode_enformer_merged.iloc[:,7] < 10])

    peak_yes = ['in_peak' for i in range(in_peak.shape[0])]
    peak_no = ['not_in_peak' for i in range(not_in.shape[0])]

    peak_dnase = pd.DataFrame(in_peak.iloc[:,3])
    peak_dnase['annot'] = peak_yes
    peak_dnase.columns = [track_name, 'annot']
    notpeak_dnase = pd.DataFrame(not_in.iloc[:,3])
    notpeak_dnase['annot'] = peak_no
    notpeak_dnase.columns = [track_name, 'annot']
    print("mean")
    print(peak_dnase[track_name].mean())
    print(notpeak_dnase[track_name].mean())

    plot_peak_data = pd.concat([peak_dnase, notpeak_dnase])
    print(plot_peak_data)

    plot_peak_data.columns = [track_name, 'value']
    single_violin_plots(index_x, index_y, axes, plot_peak_data, track_name)


def single_violin_plots(index_x, index_y, axes, data_frame, track_name):
    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)

    my_pal = {"in_peak": "#FEFE62", "not_in_peak": "#88CCEE"}
    #sns.set(font_scale = 2)
    n_in_peak = data_frame[data_frame.value == 'in_peak'].shape[0]
    n_not_in_peak = data_frame[data_frame.value == 'not_in_peak'].shape[0]

    in_peak = data_frame[data_frame.value == 'in_peak']
    not_in_peak = data_frame[data_frame.value == 'not_in_peak']
    #print("Rank sum test")
    stat, pvalue = stats.ranksums(in_peak[track_name], not_in_peak[track_name])

    sns.violinplot(ax=axes[index_x, index_y], x="value", y=track_name, data=data_frame, palette = my_pal)
    axes[index_x,index_y].set(xlabel="", ylabel=track_name)
    axes[index_x,index_y].yaxis.label.set_size(35)
    axes[index_x,index_y].xaxis.label.set_size(35)
    axes[index_x,index_y].tick_params(labelsize=25)
    axes[index_x,index_y].set_xticklabels(['in peak\n' + 'n = ' + str(n_in_peak),'not in peak\n'+ 'n = ' + str(n_not_in_peak)])

    if pvalue == 0.0:
        label_p = '$\it{p}$ < 2.2e-16'
        axes[index_x,index_y].annotate(label_p, xy=(.7, 0.95), xycoords=axes[index_x,index_y].transAxes, fontsize=25)
    else:
        label_p = '$\it{p}$ = ' + str(pvalue)
        axes[index_x,index_y].annotate(label_p, xy=(.7, 0.95), xycoords=axes[index_x,index_y].transAxes, fontsize=25)


def evolution_panels(intersect_bed, track, cell_type, evolved, naive, index, y_max_right, y_max_left, y_min_right, index_x, index_y, axes, which_plot, output_path = '/'):

    evolved = pd.DataFrame({track: evolved[:, index, :].flatten('F')})
    #print("YO")
    #print(evolved.max())
    naive = pd.DataFrame({track: naive[:, index, :].flatten('F')})
    #print(naive.max())

    #calculate ECDF
    ecdf_naive = ECDF(naive[track])
    ecdf_evolved = ECDF(evolved[track])

    encode_enformer_merged = pd.read_csv(output_path + intersect_bed, sep = '\t', header = None)
    #print("here")

    encode_enformer_merged.columns = ['binchr', 'binstart','binend','enformer_track_value', 'peakchr','peakstart','peakend','bp_overlap']
    #sort by overlap by biggest to smallest
    encode_enformer_merged = encode_enformer_merged.sort_values(by='bp_overlap', ascending=False)
    encode_enformer_merged = encode_enformer_merged.drop_duplicates(['binchr', 'binstart','binend'], keep='first')
    #duplicate = encode_enformer_merged[encode_enformer_merged.duplicated(['binstart'], keep=False)]

    in_peak = pd.DataFrame(encode_enformer_merged[encode_enformer_merged.bp_overlap >= 10])
    not_in = pd.DataFrame(encode_enformer_merged[encode_enformer_merged.bp_overlap < 10])
    #print(in_peak)
    #rint(not_in)
    peak_yes = ['in_peak' for i in range(in_peak.shape[0])]
    peak_no = ['not_in_peak' for i in range(not_in.shape[0])]

    peak_current_track = pd.DataFrame(in_peak.enformer_track_value)
    peak_current_track['annot'] = peak_yes
    peak_current_track.columns = [track, 'annot']
    notpeak_current_track = pd.DataFrame(not_in.enformer_track_value)
    notpeak_current_track['annot'] = peak_no
    notpeak_current_track.columns = [track, 'annot']

    plot_peak_data = pd.concat([peak_current_track, notpeak_current_track])

    plot_peak_data.columns = [track, 'annot']
    plot_peak_data = plot_peak_data.sort_values(track)
    plot_peak_data.annot = plot_peak_data.annot.map( {"not_in_peak" : 0, "in_peak" : 1} )
    plot_peak_data['cumulative'] = plot_peak_data.loc[::-1, 'annot'].cumsum()[::-1]
    plot_peak_data = plot_peak_data.reset_index()
    plot_peak_data.columns = ['original_index', track, 'annot', 'cumulative']

    plot_peak_data['bins'] = range(1, len(plot_peak_data) + 1)
    plot_peak_data['bins'] = plot_peak_data['bins'].values[::-1]
    plot_peak_data['ratio'] = (plot_peak_data.cumulative/plot_peak_data.bins)
    #print("80% calculation")

    #print(plot_peak_data)
    #print(plot_peak_data[plot_peak_data.ratio >= 0.8])
    plot_peak_data_80_greater = plot_peak_data[plot_peak_data.ratio >= 0.8].iloc[:1][track].iloc[0]
    #print(plot_peak_data_80_greater)

    #print("Max values")
    #print(max(ecdf_naive.x))
    #print(max(ecdf_evolved.x))
    range_limit = max(max(ecdf_evolved.x),max(ecdf_naive.x))
    #print(range_limit)
    index_list = []
    naive_list = []
    evolved_list = []
    for k in np.arange(0,range_limit,0.01):
        index_list.append(k)
        naive_list.append(1 - ecdf_naive(k))
        evolved_list.append(1 - ecdf_evolved(k))


    cdf_ratio = np.log2(np.array(evolved_list)/np.array(naive_list))

    def closest_value(input_list, input_value):
        arr = np.asarray(input_list)
        i = (np.abs(arr - input_value)).argmin()
        return arr[i]

    tmp = closest_value(index_list, plot_peak_data_80_greater)
    closest_value_index = index_list.index(tmp)

    #print('Value at 80% for ' + cell_type + ' ' + track + ' is '+  str(pow(2, cdf_ratio[closest_value_index])))
    print(cell_type + '\t' + track + '\t' + str(pow(2, cdf_ratio[closest_value_index])))
    #supplemental ecdf curve
    if which_plot == 'ecdf':
        single_ecdf_plots(index_x, index_y, axes, ecdf_naive, ecdf_evolved, track, cell_type)
    if which_plot == 'evolution_plot_all':
        single_evolution_plot(index_x, index_y, axes, plot_peak_data, track, y_max_left, y_max_right, y_min_right, index_list, cdf_ratio, cell_type)
    if which_plot == 'evolution_plot_single':

        single_evolution_plot(index_x, index_y, axes, plot_peak_data, track, y_max_left, y_max_right, y_min_right, index_list, cdf_ratio, cell_type)
        #plt.savefig(output_path + '/evolution_overlay' + '_' + track + 'ratio_axis.pdf', dpi=400, bbox_inches="tight")


def single_ecdf_plots(index_x, index_y, axes, ecdf_naive, ecdf_evolved, track_name, cell_type):
    green = (0/255, 102/255, 0/255)
    purple = (128/255, 0/255, 128/255)
    blue = (0/255, 0/255, 238/255)
    axes[index_x,index_y].plot(ecdf_naive.x, 1 - ecdf_naive.y, color = green)
    axes[index_x,index_y].plot(ecdf_evolved.x, 1 - ecdf_evolved.y, color = purple)
    axes[index_x,index_y].set_ylim(0, 0.06)
    #plt.xlim(-0.01, 19)
    axes[index_x,index_y].xaxis.label.set_size(25)
    axes[index_x,index_y].yaxis.label.set_size(25)
    axes[index_x,index_y].set(xlabel="Enformer " + cell_type + ' ' + track_name + " predictions", ylabel = "Fraction of bins > x")
    axes[index_x,index_y].legend(labels = ['naive', 'evolved'], fontsize=25)
    axes[index_x,index_y].title.set_text(track_name)
    axes[index_x,index_y].title.set_size(25)
    axes[index_x,index_y].tick_params(labelsize=25)
    axes[index_x,index_y].grid()



def single_evolution_plot(index_x, index_y, axes, plot_peak_data, track, y_max_left, y_max_right, y_min_right, index_list, cdf_ratio, cell_type):
    #overlay Evolution subplot

    ax2 = axes[index_x,index_y].twinx()
    axes[index_x,index_y].plot(plot_peak_data[track], plot_peak_data.ratio, color = '#40B0A6')
    axes[index_x,index_y].set_xlabel("Enformer Evolved " + cell_type + ' ' + track + " predictions", fontsize=30)
    axes[index_x,index_y].set_ylabel("Fraction of bins > x in peaks", fontsize=35, color = '#40B0A6')
    #ax1.set_xlim(-0.01, 19)
    axes[index_x,index_y].set_ylim(0, y_max_left)
    axes[index_x,index_y].tick_params(labelsize=25)

    ax2.plot(index_list, cdf_ratio, color = '#E1BE6A')
    #ax2.set_ylabel("log2(evolved/naive)", fontsize=15, color = '#E1BE6A')

    ax2.set_ylabel("Evolved:Naive", color = '#E1BE6A')
    ax2.set_ylim(y_min_right, y_max_right)
    tick_list = ax2.get_yticks()
    power = [(pow(2, item)) for item in tick_list]
    power = [ '%.2f' % elem for elem in power ]
    power = [str(s) + ":1" for s in power]
    #ax2.set_xlim(-0.01, 19)
    ax2.set_yticklabels(power)
    ax2.xaxis.label.set_size(35)
    ax2.yaxis.label.set_size(35)
    ax2.tick_params(labelsize=30)
    axes[index_x,index_y].grid()
    axes[index_x,index_y].title.set_text(track)
    axes[index_x,index_y].title.set_size(25)


### Read in files
evolved = np.load(hparams.path_to_predictions+'/genomic_new.npy')[:,:,0:1000]
naive = np.load(hparams.path_to_predictions+'/dinuc_local.npy')[:,:,0:1000]

#Fig 4 Panel A

fig, axes = plt.subplots(2, 2, figsize=(30, 30))
evolution_panels(hparams.bed_file, hparams.track_name, hparams.cell_type, evolved, naive, hparams.track_index,10, 1.2, -2, 0, 0, axes, 'evolution_plot_single', hparams.output_path)
plt.savefig(hparams.output_path + '/encode_enformer_intersect_' + hparams.track_name + hparams.cell_type + '_track.pdf', dpi=400, bbox_inches="tight")
