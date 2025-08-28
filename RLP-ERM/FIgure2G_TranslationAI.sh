#!/bin/bash
#$ -N translationAI
#$ -cwd
#$ -pe threaded 32
#$ -l mem_total=32G
#$ -e err.translationAI.txt
#$ -o out.translationAI.txt
###########################################################################
dataDir='/home/fanx3/translationAI'
python3 CDS_CDE_predictor2.py 2000 l1,l2,l3,l4,l5 10,10 ${dataDir}/refSeq/refSeq_hg19_compl_1.0_1.fa CDS,CDE 0