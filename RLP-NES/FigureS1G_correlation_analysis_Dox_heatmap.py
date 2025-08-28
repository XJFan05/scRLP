#! usr/bin/env python

"""
Author: Xiaojuan Fan
date: 2022-11-02
E-mail: fanx3@nih.gov
Description: Correlation analysis of FPKM from APEX-seq
Note: 
"""


import argparse
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from string import ascii_letters
import numpy as np
import pandas as pd
import math

def createHelp():
    """
    Create the command line interface of the program.
    """

    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-p', '--input-path', dest='p', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/rsem/merge.NES.TPM.anno.txt', help='input file path')
    parser.add_argument('-o', '--out-path', dest='fnOut', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/rsem/NES.TPM.correlation.png', help='output file')
    op=parser.parse_args()
    return op


if __name__ == '__main__':
    op = createHelp()
    
    data_1 = []
    data_2 = []
    out_list = []
    
    rsem_merge_file = open(op.p).readlines()
    sample_name = rsem_merge_file[0].strip().split('\t')[3:]
   # print (sample_name)
    sample_number = list(range(12))
    pair_order_list = itertools.combinations(sample_number,2)
    # print (list(pair_order_list))
    # for i in list(pair_order_list):
    for j in range(1,len(rsem_merge_file)):
        words = rsem_merge_file[j].strip().split('\t')[3:]
        #print (rsem_merge_file[j])
        TPM_list = [float(x) for x in words]
        if min(TPM_list) > 0:
            TPM_list = [math.log10(x) for x in TPM_list]
            out_list.append(TPM_list)
        #         data_1.append(TPM_list[i[0]])
        #         data_2.append(TPM_list[i[1]])

        # res = stats.spearmanr(data_1,data_2)
        # out_list.append('\t'.join([sample_name[i[0]]+'_'+sample_name[i[1]], str(res.correlation)]))
        # out_file = open(op.fnOut, 'w')
        # out_file.write('\n'.join(out_list) + '\n')

        ####################### two samples comparison
        # fig = plt.figure(figsize = (3,3))
        # ax=fig.add_subplot(111)
        # fig.subplots_adjust(bottom=0.15, right=0.9, top=0.9, left=0.15)
        # sns.set(style="ticks")
        # sns.scatterplot(x=data_1, y=data_2, alpha=0.4, palette="blue", \
        #     legend=False, s=8)
        # # plt.xlim(0.5,10000000)
        # # plt.ylim(0.5,10000000)
        # plt.xscale('log')
        # plt.yscale('log')
        # plt.savefig(op.fnOut+sample_name[i[0]]+'_'+sample_name[i[1]]+'_TPM_correlation.png',dpi=600)
        # plt.show()
        #######################

    ####################### diagonal correlation matrix
    df = pd.DataFrame(data=out_list, columns=sample_name)
    print (df)

    corr = df.corr(method='pearson')

    # Set up the matplotlib figure

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    # sns.heatmap(corr, mask=mask, cmap=cmap, 
    #             square=True, linewidths=.5, cbar_kws={"shrink": .5})

    g=sns.clustermap(corr, figsize=(7,7),
                  cmap='Blues', linewidths=0, square=True,
                  vmin = 0.8,
                  cbar_pos=(0, 0.2, 0.03, 0.4), dendrogram_ratio=(.1, .2),
                  cbar_kws={"orientation": "horizontal", "shrink": 0.2})
    # mask = np.tril(np.ones_like(corr))
    # values = g.ax_heatmap.collections[0].get_array().reshape(corr.shape)
    # new_values = np.ma.array(values, mask=mask)
    # g.ax_heatmap.collections[0].set_array(new_values)
    g.ax_cbar.set_position([0.8, 0.9, 0.15, 0.05])

    #plt.savefig(op.fnOut, dpi=600)
    plt.show()