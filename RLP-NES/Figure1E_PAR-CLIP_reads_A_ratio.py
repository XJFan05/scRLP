#! usr/bin/env python

"""
Author: Xiaojuan Fan
date: 2024-10-23
E-mail: fanx3@nih.gov
Description: Calculate the read length distribution of PAR-CLIP
Note: 
"""

from collections import defaultdict
from Bio import SeqIO
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import stats
from scipy.stats import mannwhitneyu

def createHelp():
    """
    Create the command line interface of the program.
    """

    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-p', '--path', dest='p', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/MH598_PAR/cutadapt/', help='input file path')
    parser.add_argument('-o', '--out-path', dest='fnOut', default='', help='output file')
    op=parser.parse_args()
    return op

def calculate_read_length_distribution(fasta_file):
    """
    Calculate the read length distribution from a FASTA file.
    """
    read_length_distribution = defaultdict(int)
    total_reads = 0

    # Read each record in the FASTA file
    for record in SeqIO.parse(op.p+fasta_file, "fasta"):
        read_length = len(record.seq)
        read_length_distribution[read_length] += 1
        total_reads += 1

    return dict(read_length_distribution), total_reads

def calculate_A_ratio(fasta_file):
    """
    Calculate the ratio of 'A' nucleotides in each read of a FASTA file.
    """
    A_ratios = []

    # Read each record in the FASTA file
    for record in SeqIO.parse(op.p+fasta_file, "fasta"):
        sequence = str(record.seq).upper()
        total_length = len(sequence)
        A_count = sequence.count('A')  # Count the number of 'A' nucleotides in the sequence
        
        # Calculate the ratio of 'A' in the read
        A_ratio = A_count / total_length if total_length > 0 else 0
        A_ratios.append(A_ratio * 100)  # Convert to percentage
    return A_ratios

def calculate_significance_group_comparison(control_files, target_files, control_labels, target_labels):
    """
    Calculate statistical significance between control group (5 samples) and target group (3 samples).
    Performs overall group comparison and individual sample comparisons.
    """
    print("\n" + "="*60)
    print("STATISTICAL SIGNIFICANCE ANALYSIS")
    print("5 Control Samples vs 3 Target Samples")
    print("="*60)
    
    # Collect all A ratios for each group
    all_control_ratios = []
    all_target_ratios = []
    
    print("\nCONTROL GROUP (n=5):")
    for i, (control_file, label) in enumerate(zip(control_files, control_labels)):
        control_ratios = calculate_A_ratio(control_file)
        all_control_ratios.extend(control_ratios)
        print(f"  {label}: {len(control_ratios)} reads, mean A% = {np.mean(control_ratios):.2f}")
    
    print(f"\nTARGET GROUP (n=3):")
    for i, (target_file, label) in enumerate(zip(target_files, target_labels)):
        target_ratios = calculate_A_ratio(target_file)
        all_target_ratios.extend(target_ratios)
        print(f"  {label}: {len(target_ratios)} reads, mean A% = {np.mean(target_ratios):.2f}")
    
    
    # Perform Mann-Whitney U test for overall comparison
    try:
        mw_statistic, mw_pvalue = mannwhitneyu(all_control_ratios, all_target_ratios, 
                                              alternative='two-sided')
        print(f"Mann-Whitney U test p-value: {mw_pvalue:.2e}")
    except Exception as e:
        print(f"Mann-Whitney U test failed: {e}")
        mw_pvalue = np.nan
        significance = "error"
    
    # Create summary
    summary_results = {
        'Analysis': 'Group_Comparison_5v3',
        'Mann_Whitney_p': mw_pvalue,
    }
    
    # Save results
    summary_df = pd.DataFrame([summary_results])
    
    print(f"\n{'='*60}")
    print("SUMMARY")
    print("="*60)
    print(summary_df.to_string(index=False, float_format='%.4f'))
     
    return summary_df

def plot_polyA_ratio_distribution(fasta_files, labels):
    """
    Plot the cumulative distribution of polyA ratios for multiple FASTA files using seaborn.
    """
    # Create a DataFrame to hold the cumulative distribution data
    plot_data = []

    for fasta_file, label in zip(fasta_files, labels):
        polyA_ratios = calculate_A_ratio(fasta_file)
        if len(polyA_ratios) == 0:
            continue

        # Sort the polyA ratios
        sorted_ratios = sorted(polyA_ratios)
        cumulative_percentages = []
        
        for i, ratio in enumerate(sorted_ratios):
            cumulative_percentage = (i + 1) / len(sorted_ratios) * 100
            cumulative_percentages.append({
                'PolyA_Ratio': ratio, 
                'Cumulative_Percentage': cumulative_percentage, 
                'Sample': label
            })
        
        plot_data.extend(cumulative_percentages)

    # Convert to DataFrame for plotting
    df = pd.DataFrame(plot_data)
    print(df)

    # Plot using seaborn
    plt.figure(figsize=(3, 2))
    sns.lineplot(data=df, x='PolyA_Ratio', y='Cumulative_Percentage', hue='Sample', \
                 palette=colors, linewidth=2.0, legend=False)
    
    # Customize the plot
    plt.xlabel('A Ratio (%)')
    plt.ylabel('Cumulative Percentage of Sequences')

    # Show the plot
    plt.savefig(op.p+file_prefix+'polyA_ratio_cumulative.png', dpi=600)
    plt.show()

if __name__ == '__main__':
    op = createHelp()
    
    file_prefix = 'TadA_NES_'
    
    # REORGANIZED FILE LISTS FOR EASIER PAIRING
    control_files = [
        file_prefix+'ctl_cutadapt.fa',
        file_prefix+'noUV_cutadapt.fa',
        file_prefix+'low_cutadapt.fa',
        file_prefix+'high_cutadapt.fa'
    ]
    
    target_files = [
        file_prefix+'3PABP_ctl_cutadapt.fa',
        file_prefix+'3PABP_noUV_cutadapt.fa',
        file_prefix+'3PABP_low_cutadapt.fa',
        file_prefix+'3PABP_high_cutadapt.fa'
    ]
    
    control_labels = ['ctl', 'noUV', 'low', 'high', '3PABP_ctl']
    target_labels = ['3PABP_noUV', '3PABP_low', '3PABP_high']
    
    # All files for plotting
    fasta_files = control_files + target_files
    labels = ['ctl', 'noUV', 'low', 'high', 'RLP_ctl', 'RLP_noUV', 'RLP_low', 'RLP_high']  # Labels for the legend
    colors = ['#E6E5E5','#D3ECF3', '#FACAE2', '#CCC1D9', '#A0A0A0', '#34C3EF', '#DF4291', '#8C4ED4']

    # PERFORM SIGNIFICANCE TESTING
    print("Starting significance analysis...")
    significance_results = calculate_significance_group_comparison(control_files, target_files, control_labels, target_labels)

    # Plot the cumulative PolyA ratio distributions
    plot_polyA_ratio_distribution(fasta_files, labels)