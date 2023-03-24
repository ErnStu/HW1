#!/bin/bash

threads=6
g_ref="../refs/mm10"

#check if index files exist
if [ -f ${g_ref}.1.ht2 ]; then
    echo "genome index is already done"
else 
    hisat2-build -p ${threads} ../refs/mm.fa ${g_ref}
fi

#FASTQC analysis and FASTQ trimming
for i in ../inputs/*_1.fastq.gz
do
  R1=${i}
  R2="../inputs/"$(basename ${i} _1.fastq.gz)"_2.fastq.gz"
  fastqc -t ${threads} ${R1} ${R2} -o ../outputs/fastqc
  trim_galore -j ${threads} --paired --quality 20 --length 20 --illumina -o ../outputs/ ${R1} ${R2}
done

#MULTIQC on raw files
multiqc -o ../outputs/fastqc ../outputs/fastqc

#FASTQC analysis on trimmed files
for i in ../outputs/*_1_val_1.fq.gz
do
  R1=${i}
  R2="../outputs/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
  fastqc -t ${threads} ${R1} ${R2} -o ../outputs/fastqc_clean
done

#MULTIQC on trimmed files
multiqc -o ../outputs/fastqc_clean ../outputs/fastqc_clean

#Mapping, QC and quantification
for i in ../outputs/*1_val_1.fq.gz
do
  base=$(basename $i _1_val_1.fq.gz)
  hisat2 -p ${threads} --dta -x ../refs/mm10 -1 ../outputs/${base}_1_val_1.fq.gz -2 ../outputs/${base}_2_val_2.fq.gz -S ../outputs/${base}.sam
done

for i in ../outputs/*.sam
do
  base=$(basename $i .sam)
  samtools view -F 1796 -bS -@ ${threads} ../outputs/${base}.sam | samtools sort -@ ${threads} -o ../outputs/tmp.bam
  samtools sort -@ ${threads} -n ../outputs/tmp.bam -o ../outputs/tmp2.bam
  samtools fixmate -@ ${threads} -m -r ../outputs/tmp2.bam ../outputs/tmp3.bam
  samtools sort -@ ${threads} ../outputs/tmp3.bam -o ../outputs/tmp2.bam
  samtools markdup -@ ${threads} -r -s ../outputs/tmp2.bam ../outputs/tmp3.bam
  samtools sort -@ ${threads} ../outputs/tmp3.bam -o ../outputs/${base}.bam
done

rm ../outputs/*.sam
rm ../outputs/tmp.bam ../outputs/tmp2.bam ../outputs/tmp3.bam

gtf="../refs/gtf_mm10.gtf"

for i in ../outputs/*.bam;
do
  ID=$(basename $i .bam)
  stringtie -p ${threads} -G ${gtf} -eB -o ../outputs/stringtie/${ID}/${ID}.gtf -A ../outputs/stringtie/${ID}/${ID}.tab ${i}
done

#Correlation diagrams
multiBamSummary bins --outFileName ../results/mapped.npz --binSize 1000 -p ${threads} --outRawCounts ../results/raw_counts.tsv -b outputs/*.bam
plotCorrelation --removeOutliers -in ../results/mapped.npz -c pearson -p heatmap -o ../results/mapped_data_heatmap.pdf
plotCorrelation --removeOutliers -in ../results/mapped.npz -c pearson -p scatterplot -o ../results/mapped_data_scatter.pdf
plotPCA -in ../results/mapped.npz -o ../results/mapped_data_PCA.pdf