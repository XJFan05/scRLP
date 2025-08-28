#!/usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2025-02-13
E-mail: fanx3@nih.gov
Description: 
    1) Normalize the edit counts by read count from the htseq file.
    2) Keep only genes that appear at more than one time point.
    3) Standardize each gene's normalized editing events (row-wise z-score normalization).
    4) Perform hierarchical clustering using the correlation metric and assign flat clusters.
    5) Reorder clusters based on the onset of editing (first time point > 0).
    6) Output a TSV file merging heatmap data and gene-cluster mapping in heatmap order.
    7) Print the number of genes per cluster and filtering statistics.
    8) Output a time course matrix and design file with sample information.
"""

import argparse
import os
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list

def createHelp():
    """Create the command line interface."""
    epilog_string = "Report bugs to fanxiaojuan@picb.ac.cn"
    description_string = (
        "Normalizes and standardizes editing events, filters genes present in >1 time point, "
        "clusters genes, orders clusters by editing onset, plots a heatmap, outputs a TSV merging "
        "heatmap data with gene-cluster mapping (in heatmap order), and provides a time course matrix "
        "and design file."
    )
    parser = argparse.ArgumentParser(description=description_string, epilog=epilog_string)
    parser.add_argument('-p', '--editing-path', dest='p', 
                        default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/', 
                        help='Input file path')
    parser.add_argument('-o', '--out-path', dest='fnOut', 
                        default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/Figures/', 
                        help='Output directory path')
    parser.add_argument('-n', '--n-clusters', dest='n_clusters', type=int, default=5,
                        help='Number of clusters for hierarchical clustering')
    return parser.parse_args()

def normalize_editing_events(editing_path):
    """Normalize editing events by read count for each time point."""
    time_points = ['6h', '12h', '24h', '48h', '72h']
    normalized_dfs = []
    all_genes = set()

    for time_point in time_points:
        editing_file = os.path.join(editing_path, f'TadA_ERM_PABP_{time_point}_confidence_scores.csv')
        editing_df = pd.read_csv(editing_file, sep='\t', header=0)
        editing_df = editing_df[editing_df['read_count_ERM'] != 0]
        editing_df['gene'] = editing_df['geneNA']
        editing_df['normalized_editing'] = editing_df['edit_count_ERM'] / editing_df['read_count_ERM']
        editing_df['time_point'] = time_point
        normalized_dfs.append(editing_df)
        all_genes.update(editing_df['gene'])

    normalized_df = pd.concat(normalized_dfs, ignore_index=True)
    missing_entries = [
        {'gene': gene, 'time_point': tp, 'normalized_editing': 0}
        for gene in all_genes for tp in time_points
        if not ((normalized_df['gene'] == gene) & (normalized_df['time_point'] == tp)).any()
    ]
    if missing_entries:
        normalized_df = pd.concat([normalized_df, pd.DataFrame(missing_entries)], ignore_index=True)
    
    return normalized_df

def filter_genes(normalized_df):
    """Filter genes edited in more than one time point."""
    pivot = normalized_df.pivot(index='gene', columns='time_point', values='normalized_editing')
    nonzero_counts = (pivot != 0).sum(axis=1)
    genes_to_keep = nonzero_counts[nonzero_counts > 1].index   #genes edited in more than one time point
    return normalized_df[normalized_df['gene'].isin(genes_to_keep)]

def plot_heatmap(normalized_df, output_path, n_clusters):
    """Plot a heatmap of standardized editing events with ordered clusters."""
    time_points = ['6h', '12h', '24h', '48h', '72h']
    pivot_df = normalized_df.pivot(index='gene', columns='time_point', values='normalized_editing')[time_points]
    pivot_df.index.name = 'gene'  # Ensure index is named correctly
    pivot_df_std = pivot_df.apply(
        lambda row: (row - row.mean()) / row.std() if row.std() != 0 else row - row.mean(), axis=1
    )
    pivot_df_std.index.name = 'gene'  # Reinforce index name

    # Hierarchical clustering
    Z = linkage(pivot_df_std.values, method='average', metric='correlation')
    cluster_labels = fcluster(Z, t=n_clusters, criterion='maxclust')

    # Cluster-to-color mapping
    cluster_color_map = {1: "#E27508", 2: "#F72FA0", 3: "#1996F6", 4: "#963CFD", 5: "#FCB715"}
    row_colors = pd.Series(cluster_labels, index=pivot_df_std.index).map(cluster_color_map)

    # Create clustermap
    g = sns.clustermap(
        pivot_df_std,
        cmap='Blues',
        yticklabels=False,
        row_cluster=True,
        col_cluster=False,
        row_linkage=Z,
        row_colors=row_colors,
        cbar_kws={"shrink": 0.5},
        cbar_pos=(0.785, 0.35, 0.03, 0.25),
        figsize=(4, 5)
    )
    g.ax_heatmap.set_xlabel('Time Point')

    # Adjust layout
    # Get the current position of the heatmap, 4‑tuple (x0, y0, width, height) giving the axes’ bounding box
    pos = g.ax_heatmap.get_position()

    # Compute a new, narrower box. 
    new_pos = [pos.x0 + 0.02,      # shift right slightly
               pos.y0,             # same vertical placement
               pos.width * 0.75,   # 75% of the original width
               pos.height]         # same height
    g.ax_heatmap.set_position(new_pos)

    # Ajust row‐color strip
    row_colors_pos = g.ax_row_colors.get_position()
    g.ax_row_colors.set_position([new_pos[0] - 0.04, row_colors_pos.y0, 0.03, row_colors_pos.height])

    plt.tight_layout(pad=1.5)
    #g.savefig(output_path, dpi=600, bbox_inches='tight')
    plt.show()

    # Extract ordered genes and data
    ordered_genes = pivot_df_std.index[g.dendrogram_row.reordered_ind]
    ordered_data = pivot_df_std.loc[ordered_genes]
    merged_data = ordered_data.copy()
    merged_data['cluster'] = pd.Series(cluster_labels, index=pivot_df_std.index).loc[ordered_genes]
    cluster_counts = pd.Series(cluster_labels).value_counts().sort_index()

    return ordered_data, cluster_counts, merged_data

def main():
    op = createHelp()
    normalized_df = normalize_editing_events(op.p)
    filtered_df = filter_genes(normalized_df)

    # Report filtering statistics
    total_genes = normalized_df['gene'].nunique()
    kept_genes = filtered_df['gene'].nunique()
    print("Total genes before filtering:", total_genes)
    print("Genes kept (edited in >1 time point):", kept_genes)
    print("Genes filtered out:", total_genes - kept_genes)

    # Plot heatmap and get data
    heatmap_out = os.path.join(op.fnOut, 'edit_count_heatmap.png')
    heatmap_data, cluster_counts, merged_data = plot_heatmap(filtered_df, heatmap_out, op.n_clusters)

    # Gene type information
    time_points = ['6h', '12h', '24h', '48h', '72h']
    gene_type_list = [
        pd.read_csv(os.path.join(op.p, f"TadA_ERM_PABP_{tp}_confidence_scores.csv"), sep='\t', header=0)
        [['geneNA', 'geneType_ERM']].assign(time_point=tp)
        for tp in time_points
    ]
    gene_type_all = pd.concat(gene_type_list, ignore_index=True).drop_duplicates(subset=['geneNA'])
    gene_type_all.rename(columns={'geneNA': 'gene', 'geneType_ERM': 'gene_type'}, inplace=True)
    print("gene_type_all:", gene_type_all)  # Debug

    # Merge with gene type
    merged_data_reset = merged_data.reset_index()  # Should create 'gene' column
    print("merged_data_reset before rename:", merged_data_reset)  # Debug
    # Ensure the gene column is correctly named
    if 'gene' not in merged_data_reset.columns and 'time_point' in merged_data_reset.columns:
        merged_data_reset = merged_data_reset.rename(columns={'time_point': 'gene'})
    print("merged_data_reset after rename:", merged_data_reset)  # Debug
    merged_data_with_type = pd.merge(merged_data_reset, gene_type_all[['gene', 'gene_type']], on='gene', how='left')
    merged_data_with_type.set_index('gene', inplace=True)

    # Save merged data
    merged_out = os.path.join(op.fnOut, 'ERM_gene_heatmap_cluster.tsv') # Merge heatmap data with cluster mapping (genes edted in more than 1 sample)
    merged_data_with_type.to_csv(merged_out, sep='\t')
    print(f"Merged gene heatmap, cluster mapping, and gene type data saved to {merged_out}")

    # Print cluster counts
    print("Number of genes in each cluster:")
    print(cluster_counts)

if __name__ == '__main__':
    main()