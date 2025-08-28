import pandas as pd
import numpy as np
from collections import defaultdict

def parse_jacusa_file(filepath, column_index=-5):
    """
    Parse JACUSA2 output file and extract relevant information for both A→G and C→T editing
    """
    data = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip header lines starting with ##
            if line.startswith('##') or line.startswith('chr\t'):
                continue
                
            # Skip empty lines
            if not line:
                continue
                
            parts = line.split('\t')
                
            # Extract relevant columns
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            score = float(parts[4])
            strand = parts[5]
            
            # Find the reference and alternative bases
            if 'TadA' in filepath:  # Check for TadA induced A→G editing
                ref_base = parts[-2]
                alt_base = 'G'  # A→G editing
            elif 'APO' in filepath:
                ref_base = parts[-2]
                alt_base = 'T'  # C→T editing
            else:
                raise Exception("No mutation type was specified")

            
            # The last count column (sample of interest)
            gene_info = parts[-1]
            
            # Find count columns (format: A,C,G,T)
            count_str = parts[column_index] # This is the column with counts in sample of interest. -5 is the index for Dox2000, -6 is the index for Dox100
            counts = [int(x) for x in count_str.split(',')]
            count_column = counts

            target_counts = count_column  # Last column with counts
            a_count, c_count, g_count, t_count = target_counts
            
            # Determine editing type and edited count based on alternative base
            edited_count = 0
            editing_type = None
            
            if alt_base == 'G':  # A→G editing (TadA induced)
                edited_count = g_count
                editing_type = 'A_to_G'
            elif alt_base == 'T':  # C→T editing (APOBEC1 induced)
                edited_count = t_count
                editing_type = 'C_to_T'
            else:
                raise Exception("No mutation type was specified")
            
            total_count = sum(target_counts)
            # Extract gene information
            gene_id = None
            gene_type = None
            gene_name = None
            
            if 'None' not in gene_info:
                gene_parts = gene_info.split(',')
                if len(gene_parts) >= 3:
                    gene_id = gene_parts[0]
                    gene_type = gene_parts[1]
                    gene_name = gene_parts[2]
            
                data.append({
                    'chromosome': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'score': score,
                    'edited_count': edited_count,
                    'editing_type': editing_type,
                    'gene_id': gene_id,
                    'gene_type': gene_type,
                    'gene_name': gene_name
                })

    return pd.DataFrame(data)

def calculate_ecpm_by_editing_type(df):
    """
    Calculate ECPM for both A→G and C→T editing types separately
    """
    results = []
    
    # Process each editing type separately
    for editing_type in ['A_to_G', 'C_to_T']:
        # Filter data for this editing type
        type_df = df[df['editing_type'] == editing_type].copy()
        #print (type_df[5:])
        
        if type_df.empty:
            continue
            
        # Calculate total edited counts for this editing type
        total_edited_counts = type_df['edited_count'].sum()
        #print (total_edited_counts)
        
        if total_edited_counts == 0:
            continue
        
        # Group by gene and calculate statistics
        gene_stats = type_df.groupby(['gene_id', 'gene_type', 'gene_name']).agg({
            'edited_count': 'sum',
            'start': 'count'  # Count of editing sites per gene (Counts how many rows (editing sites) each gene has)
        }).reset_index()
        
        # Rename the count column
        gene_stats.rename(columns={'start': 'edit_events'}, inplace=True)
        #print (gene_stats[5:])
        
        # Calculate ECPM (similar to CPM normalization)
        gene_stats['ECPM'] = (gene_stats['edited_count'] / total_edited_counts) * 1_000_000
        
        # Add editing type information
        gene_stats['editing_type'] = editing_type
        #print (gene_stats[5:])
        
        results.append(gene_stats)
    
    # Combine results from both editing types
    if results:
        combined_results = pd.concat(results, ignore_index=True)
        return combined_results
    else:
        return pd.DataFrame()

def format_output(gene_ecpm_data, output_file):
    """
    Format output in the specified format: geneID, geneType, geneNA, edit_events, ECPM
    """
    # Select and reorder columns
    output_df = gene_ecpm_data[['gene_id', 'gene_type', 'gene_name', 'edit_events', 'edited_count', 'ECPM']].copy()
    #print (output_df[5:])
    
    output_df['ECPM'] = output_df['ECPM'].round(6)

    n_genes = output_df.shape[0]
    print(f"Total genes output (edit_events ≥3 & edited_count ≥9): {n_genes}")
    
    # Sort by ECPM in descending order
    output_df = output_df.sort_values('ECPM', ascending=False)
    
    if output_file:
        output_df.to_csv(output_file, index=False, sep='\t',
                        header=['geneID', 'geneType', 'geneNA', 'edit_events', 'edited_count', 'ECPM'])
        print(f"\nResults saved to: {output_file}")
    
    return output_df

def main(filepath, output_file, analysis_type='combined'):
    """
    Main function to process JACUSA file and calculate ECPM
    """
    print("Parsing JACUSA2 file...")
    df = parse_jacusa_file(filepath, column_index=-5)
    
    if df.empty:
        print("No data found in the file")
        return None
    
    print(f"Found {len(df)} editing sites")
    
    # Count editing types
    editing_counts = df['editing_type'].value_counts()
    print("Editing type distribution:")
    for edit_type, count in editing_counts.items():
        print(f"  {edit_type}: {count} sites")
    
    print(f"Total edited counts: {df['edited_count'].sum()}")
    
    # Calculate ECPM for each editing type separately
    print("\nCalculating ECPM for each editing type separately...")
    gene_ecpm = calculate_ecpm_by_editing_type(df)
    
    if not gene_ecpm.empty:
        print("\n=== ECPM Results by Editing Type ===")
        for editing_type in gene_ecpm['editing_type'].unique():
            print(f"\n--- {editing_type} ---")
            type_data = gene_ecpm[gene_ecpm['editing_type'] == editing_type]
            format_output(type_data, output_file)

            # Validation: check ECPM sum
            ecpm_sum = type_data['ECPM'].sum()
            print(f"\n{editing_type} - Sum of ECPM values: {ecpm_sum:.2f} (should be ~1,000,000)")
    
    return df, gene_ecpm

# Example usage
if __name__ == "__main__":
    filepath = "/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/TadA_NES_3PABP_0_vs_100_and_2000_result_cluster_anno.out"
    
    print("\n" + "="*50)
    
    # For separate analysis (A→G and C→T separately)
    print("=== SEPARATE ANALYSIS ===")
    output_file = "/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/TadA_NES_3PABP_2000_anno_ECPM.tsv"
    sites_data_sep, gene_ecpm_data_sep = main(filepath, analysis_type='separate', output_file=output_file)
