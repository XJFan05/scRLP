#! usr/bin/env python

"""
Author: Xiaojuan Fan
date: 2023-06-21
E-mail: fanx3@nih.gov
Description: Filter and annotation of JACUSA2 result
Note: 
"""

"""
file format:
---------------------------------------------------
#contig	start	end	name	stat	strand	bases11	bases21	bases22	info	filter_info	refBase
chr1	14676	14677	variant	0.4520960848494724	-	0,8,0,6	0,12,0,3	0,11,0,6	*	Y	G
chr1	14679	14680	variant	0.43335239140085946	-	0,15,0,1	0,17,0,0	0,21,0,0	*	*	G
chr1	14762	14763	variant	0.35196580357032303	-	0,26,0,1	0,27,0,0	0,29,0,0	*	*	G
"""


import argparse
import portion as P
import time as time
import os

def createHelp():
    """
    Create the command line interface of the program.
    """

    epilog_string="Any bug is welcome reported to fanxiaojuan@picb.ac.cn"
    description_string='The program is going to '
    parser = argparse.ArgumentParser(description=description_string,epilog=epilog_string)
    parser.add_argument('-p', '--input-path', dest='p', default='/Users/fanx3/Desktop/RNA_localizatioin/single-cell/RNA-seq/JACUSA_TadA/', help='input file path')
    parser.add_argument('-nc', '--ctl-number', dest='c', default=4, help='input ctl sample')
    parser.add_argument('-nt', '--treat-number', dest='t', default=2, help='input treat sample')
    parser.add_argument('-a', '--annotate', dest='a', default='/Users/fanx3/Desktop/database/Human/GENCODE/gencode.v42.annotation.gtf', help='input file path')
    op=parser.parse_args()
    return op

def editing_site_filter(editing_file):
    """
    Filter editing sites
    """
    editing_site_list = [0] * (op.c + op.t)
    editing_cluster = []
    for i in range(0,len(editing_file)):
        if '#' not in editing_file[i]:
            words = editing_file[i].strip().split('\t')
            snp_info = []
            editing_info_list = []
            for k in range(0,op.c):
                mut_info = [int(x) for x in words[6+k].split(',')]
                if mut_type == 'CT':
                    if mut_info[-1] > 0:
                        snp_info.append(1)
                    else:
                        snp_info.append(0)
                elif mut_type =='AG':
                    if mut_info[-2] > 0:
                        snp_info.append(1)
                    else:
                        snp_info.append(0)
                else:
                    raise Exception('Please input correct mutation type CT or AG!')
            if sum(snp_info) < 1:  
                for j in range(0,(op.c+op.t)):
                    editing_info = [int(x) for x in words[6+j].split(',')]
                    #print (editing_info)
                    editing_info_list.extend(editing_info)
                    if mut_type == 'CT':
                        A = (editing_info[1] > 4) and (editing_info[3] > 2)
                        B = (editing_info[0] == 0) and (editing_info[2] == 0)
                        # print (words)
                        # print (editing_info)
                        if A and B:
                            if editing_info[3] / (editing_info[1]+editing_info[3]) > 0.01:
                                editing_site_list[j] += 1
                    elif mut_type == 'AG':
                        A = (editing_info[0] > 4) and (editing_info[2] > 2)
                        B = (editing_info[1] == 0) and (editing_info[3] == 0)
                        if A and B:
                            if editing_info[2] / (editing_info[0]+editing_info[2]) > 0.01:
                                editing_site_list[j] += 1
                    else:
                        raise Exception('Please enter the correct mutation type CT or AG!')
                # z_score = float(words[4])
                #print (editing_info_list)
                if mut_type == 'CT':
                    A = (editing_info_list[-3] > 4) and (editing_info_list[-1] > 2)
                    B = (editing_info_list[-7] > 4) and (editing_info_list[-5] > 2)
                    C = (editing_info_list[-2]+editing_info_list[-4]+editing_info_list[-6]+editing_info_list[-8]==0)
                    if (A and C) or (B and C):
                        if editing_info[3] / (editing_info[1]+editing_info[3]) > 0.01:
                            editing_cluster.append(words)
                if mut_type == 'AG':
                    A = (editing_info_list[-4] > 4) and (editing_info_list[-2] > 2)
                    B = (editing_info_list[-8] > 4) and (editing_info_list[-6] > 2)
                    C = (editing_info_list[-1]+editing_info_list[-3]+editing_info_list[-5]+editing_info_list[-7]==0)
                    if (A and C) or (B and C):
                        if editing_info[2] / (editing_info[0]+editing_info[2]) > 0.01:
                            editing_cluster.append(words)
    return editing_site_list, editing_cluster

def gencode_dic(gencode_file):
    """
    Build gencode database dictionary from GTF file
    """
    gen_dic = {}
    current_gene = None
    exon_starts = []
    exon_ends = []

    for line in gencode_file:
        if line.startswith('#'):
            continue
        words_gen = line.strip().split('\t')
        feature_type = words_gen[2]
        if feature_type == 'gene':
            if current_gene:
                gen_dic.setdefault((chr_no, strand), []).append([TSS_start, TSS_end, gene_id, gene_type, gene_name, exon_starts, exon_ends])
            chr_no = words_gen[0]
            TSS_start = int(words_gen[3])
            TSS_end = int(words_gen[4])
            strand = words_gen[6]
            gene_info = words_gen[8]
            gene_id = gene_info.split('gene_id "')[1].split('"')[0]
            gene_name = gene_info.split('gene_name "')[1].split('"')[0]
            gene_type = gene_info.split('gene_type "')[1].split('"')[0]
            exon_starts = []
            exon_ends = []
            current_gene = gene_id
            if 'readthrough_gene' in gene_info:
                current_gene = None
        elif feature_type == 'exon' and current_gene:
            exon_starts.append(int(words_gen[3]))
            exon_ends.append(int(words_gen[4]))

    if current_gene:
        gen_dic.setdefault((chr_no, strand), []).append([TSS_start, TSS_end, gene_id, gene_type, gene_name, exon_starts, exon_ends])

    return gen_dic

def choose_gene(gene_list, editing_site, strand):
    chosen_gene = None
    closest_distance = float('inf')
    for gene in gene_list:
        TSS_start, TSS_end, gene_id, gene_type, gene_name, exon_starts, exon_ends = gene
        # Check if editing site is within any exon
        for exon_start, exon_end in zip(exon_starts, exon_ends):
            if exon_start <= editing_site <= exon_end:
                chosen_gene = gene
                break
        if chosen_gene:
            break
        # Check proximity to 3'-UTR region
        if gene_type == 'protein_coding' and strand == '+':
            distance_to_3utr = abs(editing_site - TSS_end)
        elif gene_type == 'protein_coding' and strand == '-':
            distance_to_3utr = abs(editing_site - TSS_start)
        else:
            distance_to_3utr = float('inf')
        
        if distance_to_3utr < closest_distance:
            closest_distance = distance_to_3utr
            chosen_gene = gene

    if not chosen_gene:
        for gene in gene_list:
            if gene[3] == 'protein_coding':
                chosen_gene = gene
                break
    if not chosen_gene:
        for gene in gene_list:
            if gene[3] == 'lncRNA':
                chosen_gene = gene
                break
    return chosen_gene[2:5]

def SNV_anno(flat_dic,snv_list):
    """
    Annotate editing site
    """
    editing_anno_list = []
    for i in range(0,len(snv_list)):
        info_snv = snv_list[i]
        chr_snv = info_snv[0]
        strand_snv = info_snv[5]
        snv_site = int(info_snv[2])
        key_snv = (chr_snv,strand_snv)
        gene_list = []
        snv_gene = []
        if key_snv in flat_dic:
            for j in range(0,len(flat_dic[key_snv])):
                trans_info = flat_dic[key_snv][j]
                trans_start = trans_info[0]
                trans_end = trans_info[1]
                #print (trans_info)
                if snv_site > trans_start and snv_site < trans_end:
                    gene_list.append(trans_info)
        if len(gene_list) > 1:
            snv_gene = choose_gene(gene_list,snv_site,strand_snv)
        elif len(gene_list) == 1:
            #print (gene_list)
            snv_gene = gene_list[0][2:5]
        else:
            snv_gene = ['None','None','None']
        #print (snv_gene)
        out_gene = ','.join(snv_gene)
        editing_anno_list.append(info_snv + [out_gene])
    return editing_anno_list

def filter_clust(cluster_list):
    """
    Filter cluster based on editing sites distance
    """
    out_list = []
    cluster_dic = {}
    for i in range(0,len(cluster_list)):
        #print (cluster_list[i])
        cluster_info = cluster_list[i]
        geneNA = cluster_info[-1]
        cluster_dic.setdefault(geneNA,[]).append(cluster_info)
    
    #print (cluster_dic['ENSG00000188157.15,protein_coding,AGRN'])
    cluster_no = 0
    cluster_gene = []
    for key in cluster_dic:
        if len(cluster_dic[key]) > 2:
            cluster_no += 1
            cluster_gene.append(key)
            for i in range(0, len(cluster_dic[key])):
                out_list.append('\t'.join(cluster_dic[key][i]))
    gene_list_file = open(op.p+file_tag+'_edited_genes.txt','w')
    gene_list_file.write('\n'.join(cluster_gene)+'\n')
    print (file_tag + ' total cluster: ' + str(cluster_no))
    return out_list

if __name__ == '__main__':
    op = createHelp()
    time_start = time.time()

    for fileNA in sorted(os.listdir(op.p)):
        if 'GJA1_stop_JACUSA2_result_MH590.out' in fileNA:
            file_tag = fileNA.strip().split('_JACUSA2')[0]
            print (file_tag + ' start ...')

            if 'APO' in file_tag:
                mut_type = 'CT'
            elif ('TadA' in file_tag) or ('GJA1' in file_tag):
                mut_type = 'AG'
            else:
                raise Exception('Editing type should be CT or AG')
            
            JACUSA_out = open(op.p+fileNA).readlines()
            editing_site_list, editing_cluster = editing_site_filter(JACUSA_out)

            print ('Total editing sites identified: ')
            print (editing_site_list)
            print ('Effective editing sites: ' + str(len(editing_cluster)))

            gencode_file = open(op.a).readlines()
            gtf_dic = gencode_dic(gencode_file)
            print ('GTF dic done!')
            print ('Time used: ' + str(time.time() - time_start) + ' s.')

            editing_sites_anno = SNV_anno(gtf_dic, editing_cluster)
            print ('Total annotated editing sites: ' + str(len(editing_sites_anno)))

            out_file_sites = open(op.p+file_tag+'_0_vs_100_and_2000_result_sites_anno.out','w')
            out_file_sites.write(JACUSA_out[0].strip() + '\t' + 'geneNA' + '\n')
            for sites in editing_sites_anno:
                sites_line = '\t'.join(sites)
                out_file_sites.write(sites_line + '\n')

            # cluster_out = cluster_anno(editing_sites_anno)
            print ('Editing sites in clusters before filter: ' + str(len(editing_sites_anno)))

            cluster_filter = filter_clust(editing_sites_anno)

            print ('Output editing sites in clusters: ' + str(len(cluster_filter)))

            out_file = open(op.p+file_tag+'_0_vs_100_and_2000_result_cluster_anno.out','w')
            out_file.write(JACUSA_out[0].strip() + '\t' + 'geneNA' + '\n')
            out_file.write('\n'.join(cluster_filter) + '\n')

            print ('Done!')
            print ('Time used: ' + str(time.time() - time_start) + ' s.')