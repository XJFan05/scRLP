#!/bin/bash
#SBATCH --partition=defq
#SBATCH --job-name=ribo_seq
#SBATCH --time=10:00:00
#SBATCH --ntasks=32
#SBATCH -e err.ribo_seq.txt
#SBATCH -o out.ribo_seq.txt
#############################################
module load star/2.7.9a
module load samtools/1.20
module load sratoolkit/3.0.2
module load bowtie/1.2.2
module load R/4.4.1
resDir='/home/fanx3/RNA_localization/single_cell_localization/ribosome_footprint/GSM4425991'
StarIndexDir='/home/fanx3/database/Homo_sapiens/GENCODE/StarIndex/'
BowtieIndexDir='/home/fanx3/database/Homo_sapiens/rRNA/BowtieIndex/'
RiboTaperDir='/home/fanx3/program/src/RiboTaper'
BedtoolsDir='/home/fanx3/program/src/conda/bin/'
human_ref='/home/fanx3/database/Homo_sapiens/GENCODE/gencode.v42.annotation.txt'

for file in ${resDir}/raw_data/*.fastq
do
base=`basename ${file} .fastq`
#echo ${base}
fastq-dump ${resDir}/raw_data/${base} -O ${resDir}/fastq/
fastqc -t 30 -o ${resDir}/fastqc/ --noextract ${resDir}/raw_data/${base}.fastq

cutadapt -e 0.1 -m 20 -j 35 -u 1 -a 'CTGTAGGCACCATCAAT' --discard-untrimmed -o ${resDir}/cutadapt/${base}_cutadapt.fastq ${resDir}/raw_data/${base}.fastq >${resDir}/cutadapt/${base}_log.txt
fastqc -t 30 -o ${resDir}/fastqc/ --noextract ${resDir}/cutadapt/${base}_cutadapt.fastq

bowtie --threads 30 --seedlen=23 --un ${resDir}/rRNA/${base}.rRNA.unmapped.fq ${BowtieIndexDir}/human_rRNA ${resDir}/cutadapt/${base}_cutadapt.fastq >${resDir}/rRNA/${base}.rRNA.sam
fastqc -t 30 -o ${resDir}/fastqc/ --noextract ${resDir}/rRNA/${base}.rRNA.unmapped.fq

STAR --runMode alignReads --runThreadN 35 --genomeDir ${StarIndexDir} \
     --sjdbGTFfile ${StarIndexDir}/gencode.v42.annotation.gtf \
     --alignEndsType Extend5pOfRead1 \
     --outSAMtype BAM Unsorted \
     --outSAMattributes All \
     --outSAMattrRGline ID:${base} SM:${base} \
     --readFilesIn ${resDir}/rRNA/${base}.rRNA.unmapped.fq \
     --outReadsUnmapped Fastx --outSJfilterReads Unique  \
     --outFileNamePrefix ${resDir}/STAR/${base}. \
     --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverReadLmax 0.04
fastqc -t 30 -o ${resDir}/fastqc/ --noextract ${resDir}/star/${base}.Unmapped.out.mate1
samtools sort -@ 35 -o ${resDir}/STAR/${base}.Aligned.sortedByCoord.out.bam ${resDir}/STAR/${base}.Aligned.out.bam
samtools index -@ 35 -b ${resDir}/STAR/${base}.Aligned.sortedByCoord.out.bam
##################################
mkdir ${resDir}/metaplots/${base}
bash create_metaplots.bash ${resDir}/STAR/${base}.Aligned.sortedByCoord.out.bam ${RiboTaperDir}/Human_annotation_hg38/start_stops_FAR.bed ${resDir}/metaplots/${base}/${base} ${BedtoolsDir} ${RiboTaperDir}/scripts/
python /home/fanx3/code/IRES_element/src/ribosome_footprint/reads_filter.py -bam ${resDir}/star/${base}.Aligned.sortedByCoord.out.bam -l 28 -o ${resDir}/star/${base}.Aligned.sortedByCoord.out.filter.bam
done
###############################


