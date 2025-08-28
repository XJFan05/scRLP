#!/usr/bin/env python

"""
Author: Modified from Xiaojuan Fan's editing analysis script
Date: 2025-7-3
Description: Performance-optimized calculation of read depth distribution profile 
at the upstream of transcription termination site for merged single-cell RNA-seq data
Note: Processes bedgraph files and creates read depth profiles upstream of polyA tails
"""

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from collections import defaultdict

def createHelp():
    """
    Create the command line interface of the program.
    """
    epilog_string = "Performance-optimized script for calculating read depth distribution profile"
    description_string = 'Calculate read depth distribution profile upstream of 3\'-polyA tail from bedgraph file'
    parser = argparse.ArgumentParser(description=description_string, epilog=epilog_string)
    parser.add_argument('-i_bedgraph', '--input-bedgraph', dest='bedgraph_file', 
                       default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/MH611_iCell8/STAR/merged_highqual_merged.bedgraph', 
                       type=str, help='input bedgraph file with read depth information')
    parser.add_argument('-i_refgene', '--input-refgene', dest='fnIn_ref', 
                       default='/Users/fanx3/Desktop/database/Human/NCBI/refGene_hg38_112123.txt', 
                       help='input refgene annotation file')
    parser.add_argument('-l', '--PAtail-upstream', dest='l', default=1000, type=int,
                       help='upstream length from transcription termination site to analyze')
    parser.add_argument('-o', '--output-file', dest='fnOut', 
                       default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/MH611_iCell8/STAR/pseudobulk_PAtail_expression_profile.txt', 
                       type=str, help='output file for read depth profile')
    op = parser.parse_args()
    return op

def parse_bedgraph_to_dict(bedgraph_file):
    """
    Parse bedgraph file and convert to per-nucleotide depth dictionary
    """
    print(f"Parsing bedgraph file: {bedgraph_file}")
    start_time = time.time()
    
    depth_dict = defaultdict(dict)
    total_intervals = 0
    
    try:
        with open(bedgraph_file, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chr_name = parts[0]
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                depth = float(parts[3])
                
                # Convert interval to per-nucleotide depth
                for pos in range(start_pos, end_pos):
                    depth_dict[chr_name][pos] = depth
                
                total_intervals += 1
                
                if total_intervals % 100000 == 0:
                    print(f"Processed {total_intervals} intervals...")
    
    except IOError:
        print(f"Error: Cannot read bedgraph file {bedgraph_file}")
        return {}
    
    # Calculate statistics
    total_positions = sum(len(chr_data) for chr_data in depth_dict.values())
    total_depth = sum(sum(chr_data.values()) for chr_data in depth_dict.values())
    avg_depth = total_depth / total_positions if total_positions > 0 else 0
    
    print(f"Bedgraph parsing completed in {time.time() - start_time:.2f} seconds")
    print(f"Processed {total_intervals} intervals covering {total_positions} positions")
    print(f"Total chromosomes: {len(depth_dict)}")
    print(f"Average read depth: {avg_depth:.2f}")
    
    return dict(depth_dict)

def refgene_dic_all_genes(refgene_file):
    """
    Build refgene dictionary for all genes (no pre-filtering needed)
    """
    print(f"Building refgene dictionary from {refgene_file}...")
    ref_dict = {}
    processed_count = 0
    
    try:
        with open(refgene_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                words = line.strip().split('\t')
                if len(words) < 15:
                    continue
                    
                gene_name = words[12]
                chr_name = words[2]
                strand = words[3]
                tss = int(words[4])
                tes = int(words[5])
                
                # Pre-parse and convert to integers
                try:
                    exon_starts = [int(x) for x in words[9].rstrip(',').split(',') if x]
                    exon_ends = [int(x) for x in words[10].rstrip(',').split(',') if x]
                    
                    if len(exon_starts) == len(exon_ends):  # Valid exon structure
                        ref_dict.setdefault(gene_name, []).append([tss, tes, exon_starts, exon_ends, strand, chr_name])
                        processed_count += 1
                except ValueError:
                    continue
        
        print(f'Refgene dictionary built: {len(ref_dict)} genes, {processed_count} transcripts')
        return ref_dict
        
    except IOError:
        print(f"Error: Cannot read refgene file {refgene_file}")
        return {}

def calculate_depth_profiles(depth_dict, ref_dict, upstream_length):
    """
    Calculate read depth distribution profiles upstream of polyA tails
    """
    print("Calculating read depth profiles...")
    result_list_max_depth = []
    
    processed_genes = 0
    genes_with_data = 0
    
    for gene_name, transcripts in ref_dict.items():
        result_list = []
        depth_totals = []
        
        for tss, tes, exon_starts, exon_ends, strand, chr_name in transcripts:
            # Initialize depth array
            depth_array = np.zeros(upstream_length, dtype=np.float32)
            remaining_length = upstream_length
            
            # Get chromosome data once
            chr_data = depth_dict.get(chr_name, {})
            if not chr_data:
                continue
            
            if strand == '+':
                # Reverse iterate through exons for + strand (from 3' to 5')
                for start, end in reversed(list(zip(exon_starts, exon_ends))):
                    if end <= tes and remaining_length > 0:
                        exon_length = end - start + 1
                        
                        if remaining_length <= exon_length:
                            # Process partial exon
                            for i in range(remaining_length):
                                pos = end - i
                                if pos in chr_data:
                                    depth_array[upstream_length - remaining_length + i] = chr_data[pos]
                            remaining_length = 0
                            break
                        else:
                            # Process full exon
                            for i in range(exon_length):
                                pos = end - i
                                if pos in chr_data and remaining_length > 0:
                                    depth_array[upstream_length - remaining_length] = chr_data[pos]
                                    remaining_length -= 1
                                    
            elif strand == '-':
                # Forward iterate through exons for - strand (from 3' to 5')
                for start, end in zip(exon_starts, exon_ends):
                    if start >= tss and remaining_length > 0:
                        exon_length = end - start + 1
                        
                        if remaining_length <= exon_length:
                            # Process partial exon
                            for i in range(remaining_length):
                                pos = start + i
                                if pos in chr_data:
                                    depth_array[upstream_length - remaining_length + i] = chr_data[pos]
                            remaining_length = 0
                            break
                        else:
                            # Process full exon
                            for i in range(exon_length):
                                pos = start + i
                                if pos in chr_data and remaining_length > 0:
                                    depth_array[upstream_length - remaining_length] = chr_data[pos]
                                    remaining_length -= 1
            
            # Reverse array to have 3'-end at position 0
            depth_array = depth_array[::-1]
            
            # Check if we have any read depth
            total_depth = np.sum(depth_array)
            if total_depth > 0:
                depth_totals.append(total_depth)
                result_list.append(depth_array)
        
        # Select transcript with maximum total depth
        if result_list:
            best_idx = np.argmax(depth_totals)
            best_result = result_list[best_idx]
            result_list_max_depth.append(best_result)
            genes_with_data += 1
            
        processed_genes += 1
        
        if processed_genes % 1000 == 0:
            print(f"Processed {processed_genes} genes, {genes_with_data} with data...")
    
    print(f'Total genes processed: {processed_genes}')
    print(f'Genes with read depth data: {genes_with_data}')
    print(f'Transcripts with depth profiles: {len(result_list_max_depth)}')
    return result_list_max_depth

def calculate_profile_statistics(profile_list):
    """
    Calculate average and standard deviation of depth profiles
    """
    if not profile_list:
        return np.array([]), np.array([])
        
    # Convert to numpy array
    profile_matrix = np.array(profile_list, dtype=np.float64)
    
    # Calculate statistics
    profile_ave = np.mean(profile_matrix, axis=0)
    profile_std = np.std(profile_matrix, axis=0)
    
    return profile_ave, profile_std

def main():
    start_time = time.time()
    op = createHelp()
    
    # Validate inputs
    if not os.path.exists(op.bedgraph_file):
        print(f"Error: Bedgraph file not found: {op.bedgraph_file}")
        return 1
    
    if not os.path.exists(op.fnIn_ref):
        print(f"Error: RefGene file not found: {op.fnIn_ref}")
        return 1
    
    # Create output directory
    output_dir = os.path.dirname(op.fnOut)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Parse bedgraph file
    print("Step 1: Parsing bedgraph file...")
    parse_start = time.time()
    depth_dict = parse_bedgraph_to_dict(op.bedgraph_file)
    print(f"Bedgraph parsing completed in {time.time() - parse_start:.2f} seconds")
    
    if not depth_dict:
        print("No read depth data found in bedgraph file!")
        return 1
    
    # Step 2: Build refgene dictionary
    print("Step 2: Building refgene dictionary...")
    refgene_start = time.time()
    ref_dict = refgene_dic_all_genes(op.fnIn_ref)
    print(f'RefGene processing completed in {time.time() - refgene_start:.2f} seconds')
    
    if not ref_dict:
        print("No genes found in refgene annotation!")
        return 1
    
    # Step 3: Calculate read depth profiles
    print("Step 3: Calculating read depth profiles...")
    profile_start = time.time()
    results = calculate_depth_profiles(depth_dict, ref_dict, op.l)
    print(f"Profile calculation completed in {time.time() - profile_start:.2f} seconds")
    
    if not results:
        print("No read depth profiles generated!")
        return 1
    
    # Step 4: Compute average profile
    print("Step 4: Computing average read depth profile...")
    avg_start = time.time()
    profile_ave, profile_std = calculate_profile_statistics(results)
    print(f"Statistics calculation completed in {time.time() - avg_start:.2f} seconds")
    
    # Step 5: Save results
    print("Step 5: Saving results...")
    save_start = time.time()
    
    # Convert to strings efficiently
    profile_ave_str = [f"{x:.6f}" for x in profile_ave]
    profile_std_str = [f"{x:.6f}" for x in profile_std]
    
    with open(op.fnOut, 'w') as fnOut:
        fnOut.write('\t'.join(['profile_pseudobulk'] + profile_ave_str) + '\n')
        fnOut.write('\t'.join(['std_pseudobulk'] + profile_std_str) + '\n')
    
    print(f"Results saved to: {op.fnOut}")
    print(f"Saving completed in {time.time() - save_start:.2f} seconds")
    
    # Print summary statistics
    print(f"\nSummary:")
    print(f"Total transcripts analyzed: {len(results)}")
    print(f"Mean read depth: {np.mean(profile_ave):.6f}")
    print(f"Max read depth: {np.max(profile_ave):.6f}")
    print(f"Read depth at 3'-end (position 0): {profile_ave[0]:.6f}")
    print(f"Read depth at 5'-end (position {op.l-1}): {profile_ave[op.l-1]:.6f}")
    
    total_time = time.time() - start_time
    print(f"\nTotal processing time: {total_time:.2f} seconds")
    print(f"Performance: {len(results)/total_time:.2f} transcripts/second")
    
    return 0

if __name__ == '__main__':
    exit(main())