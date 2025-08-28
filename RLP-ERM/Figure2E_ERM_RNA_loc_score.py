#! usr/bin/env python

"""
Author: Xiaojuan Fan
date: 2025-02-06
E-mail: fanx3@nih.gov
Description: 
Note: 
"""


import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def createHelp():
    """
    Create the command line interface of the program.
    """

    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-p', '--input-path', dest='p', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/', help='input file path')
    parser.add_argument('-c', '--htseq-path', dest='c', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/htseq/', help='input file path')
    parser.add_argument('-o', '--out-path', dest='fnOut', default='', help='output file')
    op=parser.parse_args()
    return op

def confidence_score(df):
    """
    Compute confidence score based on the ratio of editing events to read count in ERM and cytoplasm.
    """
    # Calculate the ratio of editing events to read count for ERM and cytoplasm
    df['editing_efficiency_ERM'] = df['edit_count_ERM'] / df['read_count_ERM']
    df['editing_efficiency_cytoplasm'] = df['edit_count_cytoplasm'] / df['read_count_cytoplasm']
    #print (df['read_count_ERM'])

    # Compute confidence score
    df['confidence_score'] = df['editing_efficiency_ERM'] / df['editing_efficiency_cytoplasm']

    # Handle cases where editing_ratio_cytoplasm is zero
    df.loc[df['editing_efficiency_cytoplasm'] == 0, 'confidence_score'] = 1000

    # Take log2 of the confidence score
    df['confidence_score'] = np.log2(df['confidence_score'])

    return df

def exp_count_extract(htseq_file):
    """
    Extract the expression value from htseq file
    """
    htseq_dic = {}
    for i in range(1, len(htseq_file)):
        words = htseq_file[i].strip().split('\t')
        geneID = words[0]
        read_count = int(words[1])
        htseq_dic[geneID] = read_count
    return htseq_dic

def edit_event_calculate(JACUSA_cluster, htseq_dic, edit_index):
    """
    Calculate the RNA editing events in each edited RNA
    """
    edited_gene_dic = {}
    for i in range(1, len(JACUSA_cluster)):
        words = JACUSA_cluster[i].strip().split('\t')
        gene_info = words[-1]
        chromosome = words[0]
        editing_position = words[1]
        strand = words[5]
        if 'ENSG' in gene_info:
            geneNA = gene_info.split(',')[-1]
            geneID = gene_info.split(',')[0]
            geneType = gene_info.split(',')[1]
            edit_reads = int(words[edit_index].split(',')[-2])
            if edit_reads >= 3:
                ############################### editing events
                try:
                    edited_gene_dic[geneID]['edit_count'] += edit_reads
                except KeyError:
                    edited_gene_dic[geneID] = {
                        'chromosome': chromosome,
                        'editing_position': editing_position,
                        'geneID': geneID,
                        'geneNA': geneNA,
                        'geneType': geneType,
                        'strand': strand,
                        'edit_count': edit_reads,
                        'read_count': htseq_dic.get(geneID, 0)
                    }
    return edited_gene_dic

def read_gene_list(file_path):
    """
    Read gene names from a file, skipping comment lines starting with '#'.
    """
    genes = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):  # Skip empty lines and comments
                genes.append(line)
    return genes

if __name__ == '__main__':
    op = createHelp()
    
    Time_flag = '72h'  
    edit_index = -5   # 6h:-9, 12h:-8, 24h:-7, 48h:-6, 72h:-5
    JACUSA_fileNA_ERM = op.p + 'TadA_ERM_3PABP_time_points_result_cluster_anno_' + Time_flag + '.out'
    JACUSA_fileNA_CYT = op.p + 'TadA_NES_3PABP_0_vs_100_and_2000_result_cluster_anno.out'
    htseq_fileNA_ERM = op.c + 'TadA_ERM_3PABP_long_' + Time_flag + '_htseq.txt'
    htseq_fileNA_CYT = op.c + 'TadA_NES_3PABP_2000_htseq.txt'
    three_sets_genes = op.p + 'Common_RNAs_three_sets.txt'

    # Read data
    JACUSA_ERM_file = open(JACUSA_fileNA_ERM).readlines()
    JACUSA_CYT_file = open(JACUSA_fileNA_CYT).readlines()
    htseq_ERM_file = open(htseq_fileNA_ERM).readlines()
    htseq_CYT_file = open(htseq_fileNA_CYT).readlines()
    htseq_ERM_dic = exp_count_extract(htseq_ERM_file)
    htseq_CYT_dic = exp_count_extract(htseq_CYT_file)


    # Extract editing information
    edited_gene_ERM_dic = edit_event_calculate(JACUSA_ERM_file, htseq_ERM_dic, edit_index)
    edited_gene_CYT_dic = edit_event_calculate(JACUSA_CYT_file, htseq_CYT_dic, -5)

    # Convert edited_gene_dic to DataFrame
    df_ERM = pd.DataFrame.from_dict(edited_gene_ERM_dic, orient='index')
    df_CYT = pd.DataFrame.from_dict(edited_gene_CYT_dic, orient='index')

    # Debugging: Print column names to check for consistency
    print("Columns in df_ERM:", df_ERM.columns)
    print("Columns in df_CYT:", df_CYT.columns)

    # Merge ERM and CYT data
    df = pd.merge(df_ERM, df_CYT, on='geneNA', suffixes=('_ERM', '_cytoplasm'))

    # Debugging: Print column names after merge
    print("Columns in merged df:", df.columns)

    # Calculate confidence score
    df = confidence_score(df)

    # Filter genes with read count > 0 in both ERM and cytoplasm
    print(f"Total genes before filtering: {len(df)}")
    df = df[(df['read_count_ERM'] > 0) & (df['read_count_cytoplasm'] > 0)]
    print(f"Total genes after filtering (read count > 0 in both conditions): {len(df)}")

    # Sort the DataFrame based on confidence score
    df = df.sort_values('confidence_score', ascending=False)

    # Select and reorder the columns for output
    output_columns = [
        'geneID_ERM', 'geneNA', 'geneType_ERM', 'strand_ERM',
        'edit_count_ERM', 'read_count_ERM', 'edit_count_cytoplasm', 'read_count_cytoplasm',
        'editing_efficiency_ERM', 'editing_efficiency_cytoplasm', 'confidence_score'
    ]
    df = df[output_columns]

    # Output the results
    output_file = f"{op.p}/TadA_ERM_PABP_{Time_flag}_confidence_scores.csv"
    df.to_csv(output_file, sep='\t', index=False)

    # Read subset of genes
    subset_genes = read_gene_list(three_sets_genes)
    subset_df = df[df['geneNA'].isin(subset_genes)]

    # Calculate and print mean and percentile confidence scores for subset
    mean_confidence_subset = subset_df['confidence_score'].mean()
    percentiles = subset_df['confidence_score'].quantile([0.75, 0.50, 0.25])
    print(f"Statistics for Common_RNAs_three_sets ({len(subset_df)} genes):")
    print(f"Mean Confidence Score: {mean_confidence_subset:.4f}")
    print(f"75th Percentile (Q3): {percentiles[0.75]:.4f}")
    print(f"50th Percentile (Median): {percentiles[0.50]:.4f}")
    print(f"25th Percentile (Q1): {percentiles[0.25]:.4f}")

    # Plot distributions
    plt.figure(figsize=(2.5, 1.8))  # Slightly larger for two distributions
    # Total confidence score distribution
    sns.histplot(df['confidence_score'], kde=True, bins=40, color='brown', alpha=0.1, 
                 label='All Genes', stat='density',linewidth=0.1)
    # Subset confidence score distribution
    sns.histplot(subset_df['confidence_score'], kde=True, bins=40, color='teal', alpha=0.1, 
                 label=f'Subset ({len(subset_df)} Genes)', stat='density',linewidth=0.1)
    
    plt.xlabel("Confidence Score (log2)", fontsize=8, fontfamily='Arial')
    plt.ylabel("Density", fontsize=8, fontfamily='Arial')

    plt.xlim(-12, 2)
    plt.legend(fontsize=4, frameon=False)
    plt.tight_layout(pad=0.5)

    # Save the plot
    plot_file = f"{op.p}/TadA_ERM_PABP_{Time_flag}_confidence_score_distribution.png"
    plt.savefig(plot_file, dpi=600, bbox_inches='tight')
    plt.show()