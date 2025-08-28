#! usr/bin/env python

"""
Author: Xiaojuan Fan
date: 2023-11-21
E-mail: fanx3@nih.gov
Description: Plot the editing site distribution among mRNA
Note: 
"""


import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def createHelp():
    """
    Create the command line interface of the program.
    """

    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-i', '--input-file', dest='i', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/GJA1_stop_0_vs_100_and_2000_result_cluster_anno.out', help='input file path')
    parser.add_argument('-i_refgene', '--input-refgene', dest='fnIn_ref', default='/Users/fanx3/Desktop/database/Human/NCBI/refGene_hg38_112123.txt', help='input refgene')
    parser.add_argument('-w_UTR5', '--window-UTR5', dest='w_utr5', default=30, type=int,help='input window number in UTR5')
    parser.add_argument('-w_CDS', '--window-CDS', dest='w_cds', default=60, type=int,help='input window number in CDS')
    parser.add_argument('-w_UTR3', '--window-UTR3', dest='w_utr3', default=30, type=int,help='input window number in UTR3')
    parser.add_argument('-over', '--overlap-percent', dest='over', default=0.5, type=float, help='overlap percentage of two adjacent windows')
    parser.add_argument('-m', '--mutation', dest='m', default='AG', help='CT or AG')
    parser.add_argument('-o', '--out-path', dest='fnOut', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/GJA1_stop_0_vs_100_and_2000_result_cluster_metaGenePlot.png', help='output file')
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
                CDS=int(words[6])
                CDE=int(words[7])
                exon_start_list_include_utr5=[int(a) for a in words[9][0:-1].split(',')]
                exon_end_list_include_utr3=[int(a) for a in words[10][0:-1].split(',')]
                ref_dict.setdefault(geneNA,[]).append([TSS,TES,CDS,CDE,exon_start_list_include_utr5,\
                                                    exon_end_list_include_utr3,strand,geneNA,chrNA])
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
        # print (words)
        if len(words) == 16:  # change the value in different sample
            #print (words)
            chrNA=words[0]
            chrPOS=words[1]
            geneNA = words[-1].split(';')[0].split(',')[-1]
            gene_type = words[-1].split(';')[0].split(',')[1]
            if gene_type == 'protein_coding':
                # print (words)
                gene_list.append(geneNA)
                ctl_editing_list = [int(x) for x in words[6].split(',')]
                editing_100_list = [int(x) for x in words[-6].split(',')]
                editing_2000_list = [int(x) for x in words[-5].split(',')]
                if op.m == 'CT':
                    coverage_0 = ctl_editing_list[3] / (ctl_editing_list[1] + ctl_editing_list[3])
                    coverage_100 = editing_100_list[3] / (editing_100_list[1] + editing_100_list[3])
                    coverage_2000 = editing_2000_list[3] / (editing_2000_list[1] + editing_2000_list[3])
                elif op.m == 'AG':
                    if editing_2000_list[0] == 0:
                        editing_2000_list[0] = 1
                    if editing_100_list[0] == 0:
                        editing_100_list[0] = 1
                    coverage_0 = ctl_editing_list[2] / (ctl_editing_list[0] + ctl_editing_list[2])
                    coverage_100 = editing_100_list[2] / (editing_100_list[0] + editing_100_list[2])
                    coverage_2000 = editing_2000_list[2] / (editing_2000_list[0] + editing_2000_list[2])
                vcf_dict0.setdefault(chrNA, {})[chrPOS]=coverage_0
                vcf_dict100.setdefault(chrNA, {})[chrPOS]=coverage_100
                vcf_dict2000.setdefault(chrNA, {})[chrPOS]=coverage_2000
    gene_list = list(set(gene_list))
    # print (len(gene_list))
    print ('vcf dictionary construction finished, time used: '+str(time.time()-start_time))
    return vcf_dict0,vcf_dict100,vcf_dict2000,gene_list

def coverage_count(vcf_dic):
    """
    Count average coverage across UTR5, UTR3, CDS, and full transcript
    """
    result_list_max_coverage = []
    for key in ref_dic:
        # if key == 'SDF4':
        #     print (ref_dic[key])
        result_list = []
        gene_coverage = []
        for pos_list in ref_dic[key]:
            #print pos_list
            query_pos_list=range(pos_list[0],pos_list[1])
            #print query_pos_list
            query_pos_list=[str(a) for a in query_pos_list]
        ########calculate coverage for each nucleotide
            coverage_list=[0]*len(query_pos_list)
            for i in range(len(query_pos_list)):
                try:
                    coverage_list[i]=vcf_dic[pos_list[-1]][query_pos_list[i]]
                    # print (coverage_list[i])
                except:
                    coverage_list[i]=0
            # print (coverage_list)
    ########calculate UTR5 average coverage
            coverage_list_without_intron_UTR5=[]
            ave_win_UTR5_coverage_list = []
            #print pos_list[6]
            #if ave_transcript_coverage != 0:
            # if key == 'SDF4':
            #     print (pos_list)
            if pos_list[6] == '-':
                if (pos_list[1] == pos_list[3]): #no UTR5
                    break
                else:
                    coverage_list_without_intron_UTR5=[]
                    UTR5_exon_start_list = [a for a in pos_list[4] if a > pos_list[3]]
                    UTR5_exon_end_list = [a for a in pos_list[5] if a > pos_list[3]]
                    if len(UTR5_exon_end_list) > len(UTR5_exon_start_list):
                        UTR5_exon_start_list.insert(0,pos_list[3])
                    if len(UTR5_exon_end_list) != len(UTR5_exon_start_list):
                        raise Exception('UTR5 exon start and end not in pair')
                    for ind in range(len(UTR5_exon_start_list)):
                        coverage_list_without_intron_UTR5.extend(coverage_list[UTR5_exon_start_list[ind]-pos_list[0]:UTR5_exon_end_list[ind]-pos_list[0]])
                    coverage_list_without_intron_UTR5 = coverage_list_without_intron_UTR5[::-1]
                    UTR5_len = len(coverage_list_without_intron_UTR5)
                    if UTR5_len < (op.w_utr5):
                        break
                    else:
                        win_size_UTR5 = int(UTR5_len/(op.w_utr5*(1-op.over)+op.over))
                        step_size_UTR5 = int(win_size_UTR5*(1-op.over))
                        ave_win_UTR5_coverage_list = [0]*op.w_utr5
                        for i in range(op.w_utr5):
                            win_start_pos=i*step_size_UTR5
                            win_coverage=sum(coverage_list_without_intron_UTR5[win_start_pos:win_start_pos+win_size_UTR5])/win_size_UTR5
                            ave_win_UTR5_coverage_list[i]=win_coverage
            elif pos_list[6] == '+':
                if pos_list[0] == pos_list[2]: #no UTR5
                    break
                else:
                    UTR5_exon_start_list = [a for a in pos_list[4] if a < pos_list[2]]
                    UTR5_exon_end_list = [a for a in pos_list[5] if a < pos_list[2]]
                    if len(UTR5_exon_end_list) < len(UTR5_exon_start_list):
                        UTR5_exon_end_list.append(pos_list[2])
                    if len(UTR5_exon_end_list) != len(UTR5_exon_start_list):
                        raise Exception('UTR5 exon start and end not in pair')
                    for ind in range(len(UTR5_exon_start_list)):
                        coverage_list_without_intron_UTR5.extend(coverage_list[UTR5_exon_start_list[ind]-pos_list[0]:UTR5_exon_end_list[ind]-pos_list[0]])
                    UTR5_len = len(coverage_list_without_intron_UTR5)
                    if UTR5_len < (op.w_utr5): # UTR5 length is too short to split
                        break
                    else: # With overlap of each window
                        win_size_UTR5 = int(UTR5_len/(op.w_utr5*(1-op.over)+op.over))
                        step_size_UTR5 = int(win_size_UTR5*(1-op.over))
                        ave_win_UTR5_coverage_list = [0]*op.w_utr5
                        for i in range(op.w_utr5):
                            win_start_pos=i*step_size_UTR5
                            win_coverage=sum(coverage_list_without_intron_UTR5[win_start_pos:win_start_pos+win_size_UTR5])/win_size_UTR5
                            ave_win_UTR5_coverage_list[i]=win_coverage
    #------------------------------------------------------------------------------ 
    ########calculate CDS average coverage
            coverage_list_without_intron_CDS=[]
            CDS_exon_start_list = [a for a in pos_list[4] if a >= pos_list[2] and a <= pos_list[3]]
            CDS_exon_end_list = [a for a in pos_list[5] if a >= pos_list[2] and a <= pos_list[3]]
            #print CDS_exon_start_list,CDS_exon_end_list
            if len(CDS_exon_start_list) == 0 or len(CDS_exon_end_list) == 0:
                CDS_exon_start_list = [pos_list[2]]
                CDS_exon_end_list = [pos_list[3]]
            #print CDS_exon_start_list,CDS_exon_end_list
            if CDS_exon_start_list[0] != pos_list[2]:
                CDS_exon_start_list.insert(0,pos_list[2])
            if CDS_exon_end_list[-1] != pos_list[3]:
                CDS_exon_end_list.append(pos_list[3])
            if len(CDS_exon_end_list) != len(CDS_exon_start_list):
                print (CDS_exon_start_list,CDS_exon_end_list)
                raise Exception('CDS exon start and end not in pair')
            for ind in range(len(CDS_exon_start_list)):
                coverage_list_without_intron_CDS.extend(coverage_list[CDS_exon_start_list[ind]-pos_list[0]:CDS_exon_end_list[ind]-pos_list[0]])
            CDS_len = len(coverage_list_without_intron_CDS)
            if CDS_len < (op.w_cds): # CDS length is too short to split
                break
            else: # No overlap among each window
                win_size_CDS = int(CDS_len/(op.w_cds*(1-op.over)+op.over))
                step_size_CDS = int(win_size_CDS*(1-op.over))
                ave_win_CDS_coverage_list = [0]*op.w_cds
                for i in range(op.w_cds):
                    win_start_pos=i*step_size_CDS
                    if pos_list[6] == '-':
                        coverage_list_without_intron_CDS = coverage_list_without_intron_CDS[::-1]
                        win_coverage=sum(coverage_list_without_intron_CDS[win_start_pos:win_start_pos+win_size_CDS])/win_size_CDS
                        ave_win_CDS_coverage_list[i]=win_coverage
                    elif pos_list[6] == '+':
                        coverage_list_without_intron_CDS = coverage_list_without_intron_CDS
                        win_coverage=sum(coverage_list_without_intron_CDS[win_start_pos:win_start_pos+win_size_CDS])/win_size_CDS
                        ave_win_CDS_coverage_list[i]=win_coverage
    #------------------------------------------------------------------------------ 
    ########calculate UTR3 average coverage
            coverage_list_without_intron_UTR3=[]
            ave_win_UTR3_coverage_list = []
            if pos_list[6] == '-':
                if pos_list[0] == pos_list[2]: #no UTR3
                    break
                else:
                    UTR3_exon_start_list = [a for a in pos_list[4] if a < pos_list[2]]
                    UTR3_exon_end_list = [a for a in pos_list[5] if a < pos_list[2]]
                    if len(UTR3_exon_end_list) < len(UTR3_exon_start_list):
                        UTR3_exon_end_list.append(pos_list[2])
                    if len(UTR3_exon_end_list) != len(UTR3_exon_start_list):
                        raise Exception('UTR3 exon start and end not in pair')
                    for ind in range(len(UTR3_exon_start_list)):
                        coverage_list_without_intron_UTR3.extend(coverage_list[UTR3_exon_start_list[ind]-pos_list[0]:UTR3_exon_end_list[ind]-pos_list[0]])
                    coverage_list_without_intron_UTR3 = coverage_list_without_intron_UTR3[::-1]
                    #print coverage_list_without_intron_UTR3
                    UTR3_len = len(coverage_list_without_intron_UTR3)
                    if UTR3_len < (op.w_utr3): # UTR3 length is too short to split
                        break
                    else:
                        win_size_UTR3 = int(UTR3_len/(op.w_utr3*(1-op.over)+op.over))
                        step_size_UTR3 = int(win_size_UTR3*(1-op.over))
                        ave_win_UTR3_coverage_list = [0]*op.w_utr3
                        for i in range(op.w_utr3):
                            win_start_pos=i*step_size_UTR3
                            win_coverage=sum(coverage_list_without_intron_UTR3[win_start_pos:win_start_pos+win_size_UTR3])/win_size_UTR3
                            ave_win_UTR3_coverage_list[i]=win_coverage
            elif pos_list[6] == '+':
                if pos_list[1] == pos_list[3]: #no UTR3
                    break
                else:
                    UTR3_exon_start_list = [a for a in pos_list[4] if a > pos_list[3]]
                    UTR3_exon_end_list = [a for a in pos_list[5] if a > pos_list[3]]
                    if len(UTR3_exon_end_list) > len(UTR3_exon_start_list):
                        UTR3_exon_start_list.insert(0,pos_list[3])
                    if len(UTR3_exon_end_list) != len(UTR3_exon_start_list):
                        raise Exception('UTR3 exon start and end not in pair')
                    for ind in range(len(UTR3_exon_start_list)):
                        coverage_list_without_intron_UTR3.extend(coverage_list[UTR3_exon_start_list[ind]-pos_list[0]:UTR3_exon_end_list[ind]-pos_list[0]])
                    UTR3_len = len(coverage_list_without_intron_UTR3)
                    if UTR3_len < (op.w_utr3): # UTR3 length is too short to split
                        break
                    else:
                        win_size_UTR3 = int(UTR3_len/(op.w_utr3*(1-op.over)+op.over))
                        step_size_UTR3 = int(win_size_UTR3*(1-op.over))
                        ave_win_UTR3_coverage_list = [0]*op.w_utr3
                        for i in range(op.w_utr3):
                            win_start_pos=i*step_size_UTR3
                            win_coverage=sum(coverage_list_without_intron_UTR3[win_start_pos:win_start_pos+win_size_UTR3])/win_size_UTR3
                            ave_win_UTR3_coverage_list[i]=win_coverage
    #------------------------------------------------------------------------------ 
            transcript_list = ave_win_UTR5_coverage_list+ave_win_CDS_coverage_list+ave_win_UTR3_coverage_list
            transcript_length_effect = len(transcript_list) - transcript_list.count(0)
            if transcript_length_effect != 0:
                gene_coverage.append(transcript_length_effect)
                result_list.append(transcript_list)
        if result_list != []:
            result_list_max_coverage.append(result_list[gene_coverage.index(max(gene_coverage))])
    print ('Total used transcript: ' + str(len(result_list_max_coverage)) + ', time used: '+str(time.time()-start_time))
    return result_list_max_coverage

def matrix_average(profile_list):
    """
    Calculate lines' average
    """
    #profile_list = []
    profile_matrix = np.array(profile_list)
    #print profile_matrix.shape
    profile_ave = np.mean(profile_matrix,axis=0)
    profile_std = np.std(profile_matrix,axis=0)
    return (profile_ave)

if __name__ == '__main__':
    start_time=time.time()
    op = createHelp()

    editing_file = open(op.i).readlines()
    vcf_dic0,vcf_dic100,vcf_dic2000,edited_genes = vcf_dic_build(editing_file)

    print ('Start to build refgene dic...')
    refgene = open(op.fnIn_ref).readlines()
    ref_dic = refgene_dic(refgene)
    # print (ref_dic['HES4'])
    print ('Total genes used: ' + str(len(ref_dic.keys())))

    space_UTR5='\t'.join(['']*(op.w_utr5))
    space_UTR3='\t'.join(['']*(op.w_utr3))
    space_CDS='\t'.join(['']*(op.w_cds))
    
    results_0 = coverage_count(vcf_dic0)
    results_100 = coverage_count(vcf_dic100)
    results_2000 = coverage_count(vcf_dic2000)

    pro_0 = matrix_average(results_0)
    pro_100 = matrix_average(results_100)
    pro_2000 = matrix_average(results_2000)

    pro_0 = [0] * (op.w_utr5 + op.w_cds + op.w_utr3)

    data = {
        'Position': list(range(120)),  # Example x-axis positions (0 to 60)
        'Dox_0': pro_0,              # Replace with actual pro_0 data
        'Dox_100': pro_100,          # Replace with actual pro_100 data
        'Dox_2000': pro_2000         # Replace with actual pro_2000 data
    }
    #print (data)

    # Convert to long format for Seaborn
    df = pd.DataFrame(data)
    long_df = pd.melt(df, id_vars=['Position'], var_name='Condition', value_name='Value')

    # Plot using Seaborn
    sns.set(style="ticks", font_scale=0.8)
    plt.figure(figsize=(2.2, 1.5))

    # Create the line plot
    sns.lineplot(
        data=long_df, x='Position', y='Value', hue='Condition',
        palette={'Dox_0': 'grey', 'Dox_100': 'skyblue', 'Dox_2000': 'blue'},
        linewidth=1, legend=False  # Remove the legend
    )

    # Customize the plot
    plt.xticks([0, 30, 90, 120])
    plt.axvline(x=30, linewidth=0.2, color='grey', linestyle='-')
    plt.axvline(x=90, linewidth=0.2, color='grey', linestyle='-')
    #plt.legend(loc='best', frameon=False, fontsize=8, labelspacing=0.5)
    plt.subplots_adjust(bottom=0.15, right=0.9, top=0.9, left=0.2)

    ######################################
    plt.savefig(op.fnOut,dpi=600)
    plt.show()
            
    print ('Done!')
    ##################################