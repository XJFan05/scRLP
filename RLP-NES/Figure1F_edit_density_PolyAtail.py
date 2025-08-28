#! usr/bin/env python

"""
Author: Xiaojuan Fan
Date: 2025-1-10
E-mail: fanxiaojuan@picb.ac.cn
Description: Calculate the editing events distribution at the upstream of transcription termination site (absolute length, without intron)
"""

"""
Note: remove intron or not; strand or strandless
"""

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt

def createHelp():
    """
    Create the command line interface of the program.
    """
    
    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to calculate the editing events distribution at the upstream of TTS'
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i_vcf', '--input file', dest='fnIn', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/TadA_NES_3PABP_0_vs_100_and_2000_result_sites_anno.out', type=str,help='input vcf file')
    parser.add_argument('-i_refgene', '--input-refgene', dest='fnIn_ref', default='/Users/fanx3/Desktop/database/Human/NCBI/refGene_hg38_112123.txt', help='input refgene')
    parser.add_argument('-l', '--PAtail-upstream', dest='l', default=1000, type=int,help='input upstream length of transcription termination site')
    parser.add_argument('-s', '--strand', dest='s', default='', type=str, help='strand (+ or -)')
    parser.add_argument('-m', '--mutation', dest='m', default='AG', help='CT or AG')
    parser.add_argument('-o', '--output file', dest='fnOut', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/TadA_NES_3PABP_0_vs_100_and_2000_PAtail_profile.txt', type=str,help='output file')
    op=parser.parse_args()
    return op

def refgene_dic(refgene):
    """
    Build refgene dictionary
    """
    ref_dict={}
    for line in refgene:
        if '#' != line[0]:
            words=line.strip().split('\t')
            strand=words[3]
            cds_info = words[13]
            cde_info = words[14]
            geneNA=words[12]
            if geneNA in edited_genes:
                chrNA=words[2]
                TSS=int(words[4])
                TES=int(words[5])
                exon_start_list_include_utr5=[int(a) for a in words[9][0:-1].split(',')]
                exon_end_list_include_utr3=[int(a) for a in words[10][0:-1].split(',')]
                ref_dict.setdefault(geneNA,[]).append([TSS,TES,exon_start_list_include_utr5,\
                                                    exon_end_list_include_utr3,strand,chrNA])
    print ('query dictionary construction finished, time used: '+str(time.time()-start_time))
    return ref_dict

def vcf_dic_build(JACUSA_file):
    """
    Convert JACUSA output file to vcf format
    """
    vcf_dict0={}
    vcf_dict100={}
    vcf_dict2000={}
    gene_list = []
    for i in range(1,len(JACUSA_file)):
        words = JACUSA_file[i].strip().split('\t')
        #print (words)
        if len(words) == len(JACUSA_file[1].split('\t')):  # change the value in different sample
            #print (words)
            chrNA=words[0]
            chrPOS=words[1]
            geneNA = words[-1].split(',')[-1]
            gene_type = words[-1].split(',')[1]
            #print (geneNA)
            gene_list.append(geneNA)
            ctl_editing_list = [int(x) for x in words[6].split(',')]
            editing_100_list = [int(x) for x in words[-6].split(',')]
            editing_2000_list = [int(x) for x in words[-5].split(',')]
            if op.m == 'CT':
                coverage_0 = ctl_editing_list[3]
                coverage_100 = editing_100_list[3]
                coverage_2000 = editing_2000_list[3]
            elif op.m == 'AG':
                coverage_0 = ctl_editing_list[2] 
                coverage_100 = editing_100_list[2]
                coverage_2000 = editing_2000_list[2]
            vcf_dict0.setdefault(chrNA, {})[chrPOS]=coverage_0
            vcf_dict100.setdefault(chrNA, {})[chrPOS]=coverage_100
            vcf_dict2000.setdefault(chrNA, {})[chrPOS]=coverage_2000
    gene_list = list(set(gene_list))
    # print (len(gene_list))
    print ('vcf dictionary construction finished, time used: '+str(time.time()-start_time))
    return vcf_dict0,vcf_dict100,vcf_dict2000,gene_list

def coverage_polyA_regions(vcf_dic):
    """
    Count average coverage among start site and stop site
    """
    #print pos_list_set
    result_list_max_coverage = []
    editing_sites_distance = []
    for key in ref_dic:
        PAtail_coverage_list = [0] * op.l
        PAtail_coverage_total = []
        result_list = []
        #print gene_id
        #print pos_list_set
        ########calculate coverage at upstream of transcription termination site
        remaining_length = op.l
        #if key == 'ACTB':
        for TSS, TES, exon_start_list, exon_end_list, strand, chrNA in ref_dic[key]:
            if strand == '+':
                # Include upstream region from polyA site for the plus strand
                for start, end in reversed(list(zip(exon_start_list, exon_end_list))):
                    #print (start,end,TES)
                    if end <= TES:  # Exons upstream of the polyA site
                        exon_length = end - start + 1
                        if remaining_length <= exon_length:
                            for i in range(end, end - remaining_length, -1):
                                PAtail_coverage_list[op.l - remaining_length] = vcf_dic.get(chrNA, {}).get(str(i), 0)
                                remaining_length -= 1
                            break
                        else:
                            for i in range(end, start - 1, -1):
                                #print (end,start)
                                PAtail_coverage_list[op.l - remaining_length] = vcf_dic.get(chrNA, {}).get(str(i), 0)
                                remaining_length -= 1
            if strand == '-':
                # Include upstream region from polyA site for the plus strand
                for start, end in list(zip(exon_start_list, exon_end_list)):
                    if start >= TSS:  # Exons upstream of the polyA site
                        exon_length = end - start + 1
                        #print (start,end)
                        if remaining_length <= exon_length:
                            for i in range(start, start + remaining_length):
                                PAtail_coverage_list[op.l - remaining_length] = vcf_dic.get(chrNA, {}).get(str(i), 0)
                                remaining_length -= 1
                            break
                        else:
                            for i in range(start, end + 1):
                                PAtail_coverage_list[op.l - remaining_length] = vcf_dic.get(chrNA, {}).get(str(i), 0)
                                remaining_length -= 1
            PAtail_coverage_list = PAtail_coverage_list[::-1]
            #print (PAtail_coverage_list)
            if np.nansum(PAtail_coverage_list) != 0:
                PA_effect_coverage = len(PAtail_coverage_list) - PAtail_coverage_list.count(0)
                PAtail_coverage_total.append(PA_effect_coverage)
                result_list.append(PAtail_coverage_list)
        if result_list != []:
            out_result = result_list[PAtail_coverage_total.index(max(PAtail_coverage_total))]
            result_list_max_coverage.append(out_result)
            
            ########################## Calculate the distances between editing sites in the RNA sequence.
            # Find positions of non-zero editing sites
            editing_positions = [i for i, value in enumerate(out_result) if value > 0]
            # Calculate distances between consecutive editing sites
            distances = [editing_positions[i] - editing_positions[i - 1] for i in range(1, len(editing_positions))]
            editing_sites_distance.extend(distances)
            ################################
    print ('Total out transcripts: ' + str(len(result_list_max_coverage)))
    return result_list_max_coverage,editing_sites_distance

def matrix_average(profile_list):
    """
    Calculate lines' average
    """
    #profile_list = []
    profile_matrix = np.array(profile_list)
    #print profile_matrix.shape
    profile_ave = np.mean(profile_matrix,axis=0)
    profile_std = np.std(profile_matrix,axis=0)
    return (profile_ave,profile_std)


def plot_distance_distribution(distances, output_file):
    """
    Plot the distribution of editing site distances
    """
    #print (distances[:10])
    plt.figure(figsize=(2, 1.5))
    plt.subplots_adjust(left=0.12, bottom=0.15, right=0.9, top=0.9)
    plt.hist(distances, bins=200, edgecolor='white')
    # plt.xlabel('Distance Between Editing Sites (nt)')
    # plt.ylabel('Frequency (Log Scale)')
    plt.xlim(0,140)
    plt.xticks([0,50,100])
    # plt.yticks([0,5000,10000])
    #plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    # plt.savefig(op.fnOut+output_file, dpi=600)
    plt.show()

if __name__ == '__main__':
    start_time=time.time()
    op = createHelp()
    
    editing_file = open(op.fnIn).readlines()
    vcf_dic0,vcf_dic100,vcf_dic2000,edited_genes = vcf_dic_build(editing_file)
    #print (edited_genes)
    #print (vcf_dic2000['chr4']['102888477'])

    print ('Start to build refgene dic...')
    refgene = open(op.fnIn_ref).readlines()
    ref_dic = refgene_dic(refgene)
    # print (ref_dic['HES4'])
    print ('Total genes used: ' + str(len(ref_dic.keys())))

    results_0,distance_0 = coverage_polyA_regions(vcf_dic0)
    results_100,distance_100 = coverage_polyA_regions(vcf_dic100)
    results_2000,distance_2000 = coverage_polyA_regions(vcf_dic2000)
    print (len(results_2000))

    # plot_distance_distribution(distance_2000,'distance_2000.png')
    # plot_distance_distribution(distance_100,'distance_100.png')

    pro_0,std_0 = matrix_average(results_0)
    pro_100,std_100 = matrix_average(results_100)
    pro_2000,std_2000 = matrix_average(results_2000)
    pro_100 = [str(x) for x in pro_100]
    std_100 = [str(x) for x in std_100]
    pro_2000 = [str(x) for x in pro_2000]
    std_2000 = [str(x) for x in std_2000]

    fnOut=open(op.fnOut,'w')
    fnOut.write('\t'.join(['profile_Dox100']+pro_100) + '\n' + '\t'.join(['std_Dox100']+std_100) + '\n' \
            +'\t'.join(['profile_Dox2000']+pro_2000) + '\n' + '\t'.join(['std_Dox2000']+std_2000) + '\n')
    
    # ###################################  plot for editing events distribution before polyA tail
    fig = plt.figure(figsize = (3,1.8))
    ax=plt.subplot(111)
    fig.subplots_adjust(bottom=0.15, right=0.9, top=0.9, left=0.15)
    plt.plot(pro_0,label = 'Dox_0',color='grey',linewidth=0.5) #data_list[0] is transcript average value
    plt.plot(pro_100,label = 'Dox_100',color='skyblue',linewidth=0.5)
    plt.plot(pro_2000,label = 'Dox_2000',color='blue',linewidth=0.5)
    #plt.yscale('log')
    ax.legend(labelspacing=0.5, frameon = False, loc=0, fontsize=8) #bbox_to_anchor=(0.5, -0.1)
    ax=plt.subplot(111)
    #ax.set_xlim(145,160)
    # ax.set_ylim(-0.0002,0.0055)
    # ax.set_yticks([0,0.002,0.004])
    ax.set_xticks([0,200,400,600,800,1000])
    # ax.tick_params(labelbottom='off')
    # plt.axvline(x=15,linewidth=0.2,color='black',linestyle='-')
    # plt.axvline(x=45,linewidth=0.2,color='black',linestyle='-')
    # ###################################

    # ######################################
    # plt.savefig(op.fnOut,dpi=600)
    plt.show()
            
    # print ('Done!')
    # ##################################