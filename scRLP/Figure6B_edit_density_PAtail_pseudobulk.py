#!/usr/bin/env python

"""
Author: Modified from Xiaojuan Fan's script
Date: 2025-7-3
Description: Performance-optimized calculation of A-to-G RNA editing frequency profile 
at the upstream of transcription termination site for pseudobulk SMART-seq data
Note: Processes strand-corrected mutation types and creates editing frequency profiles
"""

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
from collections import defaultdict
import multiprocessing as mp
from functools import partial

def createHelp():
    """
    Create the command line interface of the program.
    """
    epilog_string = "Performance-optimized script for calculating A-to-G RNA editing frequency profile"
    description_string = 'Merge single-cell JACUSA2 results and calculate A-to-G editing frequency profile upstream of 3\'-end'
    parser = argparse.ArgumentParser(description=description_string, epilog=epilog_string)
    parser.add_argument('-i_jacusa', '--input-jacusa-dir', dest='jacusa_dir', 
                       default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/MH611_iCell8/JACUSA/', 
                       type=str, help='input JACUSA directory containing single-cell results')
    parser.add_argument('-ID', '--cell-id', dest='ID', 
                       default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/MH611_iCell8/JACUSA/singe_cell_samples_edit_filter_highQual.txt', 
                       type=str, help='input cell ID list for filter')
    parser.add_argument('-i_refgene', '--input-refgene', dest='fnIn_ref', 
                       default='/Users/fanx3/Desktop/database/Human/NCBI/refGene_hg38_112123.txt', 
                       help='input refgene annotation file')
    parser.add_argument('-l', '--PAtail-upstream', dest='l', default=1000, type=int,
                       help='upstream length from transcription termination site to analyze')
    parser.add_argument('-t', '--treatment-column', dest='t', default=7, type=int,
                       help='column number of treatment condition in JACUSA output (0-indexed)')
    parser.add_argument('-min_cells', '--minimum-cells', dest='min_cells', default=0, type=int,
                       help='minimum number of cells contributing to a site')
    parser.add_argument('-n_cores', '--num-cores', dest='n_cores', default=4, type=int,
                       help='number of CPU cores to use for parallel processing')
    parser.add_argument('-batch_size', '--batch-size', dest='batch_size', default=100, type=int,
                       help='number of files to process in each batch')
    parser.add_argument('-o', '--output-file', dest='fnOut', 
                       default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/MH611_iCell8/JACUSA/pseudobulk_PAtail_editing_profile.txt', 
                       type=str, help='output file for editing frequency profile')
    op = parser.parse_args()
    return op

def refgene_dic_optimized(refgene_file, edited_genes):
    """
    Optimized refgene dictionary building with set lookup
    """
    print(f"Building refgene dictionary from {len(edited_genes)} edited genes...")
    ref_dict = {}
    edited_genes_set = set(edited_genes)  # Convert to set for O(1) lookup
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
                if gene_name in edited_genes_set:  # O(1) lookup instead of O(n)
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

def process_jacusa_file_batch(file_batch, treatment_col):
    """
    Process a batch of JACUSA files and extract A-to-G editing sites
    """
    batch_merged_sites = defaultdict(lambda: {
        'edited_cov': 0, 
        'cells': set(), 'gene_info': None
    })
    batch_gene_list = []
    
    for jacusa_file in file_batch:
        cell_id = os.path.basename(jacusa_file).replace('_result_sites_anno.out', '')
        
        try:
            with open(jacusa_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or line.startswith('contig') or line.startswith('chr\t'):
                        continue
                        
                    words = line.strip().split('\t')
                    if len(words) < treatment_col + 1:
                        continue
                        
                    # Get mutation type from last column (already adjusted for strand)
                    mutation_type = words[-1]
                    
                    # Only process A>G mutations
                    if mutation_type != 'A>G':
                        continue
                    
                    # Parse basic information
                    chr_name = words[0]
                    pos = int(words[1])
                    strand = words[5]
                    
                    # Parse treatment counts
                    try:
                        counts = [int(x) for x in words[treatment_col].split(',')]
                        if len(counts) != 4:
                            continue
                            
                        # counts = [A, C, G, T]
                        total_cov = sum(counts)
                        if total_cov == 0:
                            continue
                        
                        # For A>G mutation, edited coverage is G count
                        if strand == '+':
                            edited_cov = counts[2]  # G count (A->G)
                        elif strand == '-':
                            edited_cov = counts[1]
                        else:
                            continue
                        
                        # Extract gene information from second-to-last column
                        gene_info = words[-2]
                        if ',' in gene_info:
                            gene_name = gene_info.split(',')[2]  # Take first part as gene name
                        else:
                            continue
                        
                        if gene_name:
                            batch_gene_list.append(gene_name)
                            
                            site_key = f"{chr_name}:{pos}"
                            
                            # Add to merged data
                            batch_merged_sites[site_key]['edited_cov'] += edited_cov
                            batch_merged_sites[site_key]['cells'].add(cell_id)
                            batch_merged_sites[site_key]['gene_info'] = gene_info
                            
                    except (ValueError, IndexError):
                        continue
                        
        except IOError:
            print(f"Warning: Could not read file {jacusa_file}")
            continue
    
    return dict(batch_merged_sites), batch_gene_list

def merge_jacusa_results_parallel(jacusa_dir, treatment_col, min_cells, n_cores, batch_size, sample_list):
    """
    Parallel processing of JACUSA files with batching
    """
    print("Merging JACUSA2 results across single cells (parallel processing)...")

    # Find all JACUSA result files
    jacusa_files_original = glob.glob(os.path.join(jacusa_dir, "sample_*_result_sites_anno.out"))
    jacusa_files = []
    for sample in jacusa_files_original:
        sample_ID = '_'.join(sample.split('/')[-1].split('_')[0:2])
        #print(sample_ID)
        #if sample_ID in sample_list:
        jacusa_files.append(sample)
    
    print(f"Found {len(jacusa_files)} JACUSA files to merge")
    
    if not jacusa_files:
        print("No JACUSA files found!")
        return {}, []
    
    # Create batches of files
    file_batches = [jacusa_files[i:i + batch_size] for i in range(0, len(jacusa_files), batch_size)]
    print(f"Processing {len(file_batches)} batches with {n_cores} cores")
    
    # Process batches in parallel
    process_func = partial(process_jacusa_file_batch, treatment_col=treatment_col)
    
    with mp.Pool(processes=n_cores) as pool:
        batch_results = pool.map(process_func, file_batches)
    
    # Merge results from all batches
    print("Merging batch results...")
    merged_sites = defaultdict(lambda: {
        'edited_cov': 0, 
        'cells': set(), 'gene_info': None
    })
    all_genes = []
    
    for batch_sites, batch_genes in batch_results:
        all_genes.extend(batch_genes)
        
        for site_key, site_data in batch_sites.items():
            merged_sites[site_key]['edited_cov'] += site_data['edited_cov']
            merged_sites[site_key]['cells'].update(site_data['cells'])
            if merged_sites[site_key]['gene_info'] is None:
                merged_sites[site_key]['gene_info'] = site_data['gene_info']
    
    # Get unique gene list
    edited_genes = list(set(all_genes))
    print(f"Found {len(edited_genes)} unique genes with A-to-G editing")
    
    # Create editing frequency dictionary
    print("Creating editing frequency dictionary...")
    editing_dict = defaultdict(dict)
    valid_sites = 0
    
    for site_key, site_data in merged_sites.items():
        if len(site_data['cells']) >= min_cells:
            chr_name, pos_str = site_key.split(':')
            pos = int(pos_str)
            
            # Calculate editing frequency
            edited_cov = site_data['edited_cov']
            
            editing_freq = edited_cov
            editing_dict[chr_name][pos] = editing_freq
            valid_sites += 1
    
    print(f"Created editing frequency dictionary with {valid_sites} sites from >= {min_cells} cells")
    print(f"Average editing frequency: {np.mean([freq for chr_data in editing_dict.values() for freq in chr_data.values()]):.4f}")
    
    return dict(editing_dict), edited_genes

def coverage_polyA_regions_optimized(editing_dic, ref_dic, upstream_length):
    """
    Optimized editing frequency calculation with numpy arrays
    """
    print("Calculating editing frequency profiles (optimized)...")
    result_list_max_editing = []
    
    processed_genes = 0
    
    for gene_name, transcripts in ref_dic.items():
        result_list = []
        editing_totals = []
        
        for tss, tes, exon_starts, exon_ends, strand, chr_name in transcripts:
            # Initialize editing frequency array
            editing_array = np.zeros(upstream_length, dtype=np.float32)
            remaining_length = upstream_length
            
            # Get chromosome data once
            chr_data = editing_dic.get(chr_name, {})
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
                                    editing_array[upstream_length - remaining_length + i] = chr_data[pos]
                            remaining_length = 0
                            break
                        else:
                            # Process full exon
                            for i in range(exon_length):
                                pos = end - i
                                if pos in chr_data and remaining_length > 0:
                                    editing_array[upstream_length - remaining_length] = chr_data[pos]
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
                                    editing_array[upstream_length - remaining_length + i] = chr_data[pos]
                            remaining_length = 0
                            break
                        else:
                            # Process full exon
                            for i in range(exon_length):
                                pos = start + i
                                if pos in chr_data and remaining_length > 0:
                                    editing_array[upstream_length - remaining_length] = chr_data[pos]
                                    remaining_length -= 1
            
            # Reverse array to have 3'-end at position 0
            editing_array = editing_array[::-1]
            #print (editing_array)
            
            # Check if we have any editing
            total_editing = np.sum(editing_array)
            if total_editing > 0:
                editing_totals.append(total_editing)
                result_list.append(editing_array)
        
        # Select transcript with maximum editing
        if result_list:
            best_idx = np.argmax(editing_totals)
            best_result = result_list[best_idx]
            result_list_max_editing.append(best_result)
            processed_genes += 1
    
    print(f'Total genes with editing analyzed: {processed_genes}')
    print(f'Total transcripts with editing: {len(result_list_max_editing)}')
    return result_list_max_editing

def matrix_average_optimized(profile_list):
    """
    Optimized matrix operations using numpy
    """
    if not profile_list:
        return np.array([]), np.array([])
        
    # Convert to numpy array once
    profile_matrix = np.array(profile_list, dtype=np.float64)
    
    # Vectorized operations
    profile_ave = np.mean(profile_matrix, axis=0)
    profile_std = np.std(profile_matrix, axis=0)
    
    return profile_ave, profile_std

def main():
    start_time = time.time()
    op = createHelp()
    
    # Validate inputs
    if not os.path.exists(op.jacusa_dir):
        print(f"Error: JACUSA directory not found: {op.jacusa_dir}")
        return 1
    
    if not os.path.exists(op.fnIn_ref):
        print(f"Error: RefGene file not found: {op.fnIn_ref}")
        return 1
    
    sampleID_file = open(op.ID).readlines()
    sample_list = [x.strip() for x in sampleID_file]
    #print (sample_list)
    
    # Create output directory
    output_dir = os.path.dirname(op.fnOut)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Parallel merging of JACUSA results
    print("Step 1: Merging JACUSA results and filtering A-to-G sites...")
    merge_start = time.time()
    editing_dic, edited_genes = merge_jacusa_results_parallel(
        op.jacusa_dir, op.t, op.min_cells, op.n_cores, op.batch_size, sample_list)
    print(f"Merging completed in {time.time() - merge_start:.2f} seconds")
    
    if not editing_dic or not edited_genes:
        print("No A-to-G editing data found after merging!")
        return 1
    
    # Step 2: Optimized refgene dictionary building
    print("Step 2: Building refgene dictionary...")
    refgene_start = time.time()
    ref_dic = refgene_dic_optimized(op.fnIn_ref, edited_genes)
    print(f'RefGene processing completed in {time.time() - refgene_start:.2f} seconds')
    print(f'Total genes used: {len(ref_dic)}')
    
    if not ref_dic:
        print("No genes found in refgene!")
        return 1
    
    # Step 3: Calculate editing frequency profiles
    print("Step 3: Calculating editing frequency profiles...")
    coverage_start = time.time()
    results = coverage_polyA_regions_optimized(editing_dic, ref_dic, op.l)
    print(f"Profile calculation completed in {time.time() - coverage_start:.2f} seconds")
    
    if not results:
        print("No editing frequency profiles generated!")
        return 1
    
    # Step 4: Compute average profile
    print("Step 4: Computing average editing frequency profile...")
    avg_start = time.time()
    profile_ave, profile_std = matrix_average_optimized(results)
    print(f"Averaging completed in {time.time() - avg_start:.2f} seconds")
    
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
    print(f"Mean editing frequency: {np.mean(profile_ave):.6f}")
    print(f"Max editing frequency: {np.max(profile_ave):.6f}")
    print(f"Editing frequency at 3'-end (position 0): {profile_ave[0]:.6f}")
    
    total_time = time.time() - start_time
    print(f"\nTotal processing time: {total_time:.2f} seconds")
    print(f"Performance: {len(results)/total_time:.2f} transcripts/second")
    
    return 0

if __name__ == '__main__':
    exit(main())