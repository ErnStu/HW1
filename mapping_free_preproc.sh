#!/bin/bash

threads=6

#check if index files exist
if [ -f refs/transcriptome/seq.bin ]; then
    echo "genome index is already done"
else 
    salmon index -p ${threads} -t refs/transcriptome.fa -i refs/transcriptome
fi

#FASTQC analysis and FASTQ trimming
for i in inputs/*_1.fastq.gz
do
  R1=${i}
  R2="inputs/"$(basename ${i} _1.fastq.gz)"_2.fastq.gz"
  fastqc -t ${threads} ${R1} ${R2} -o ../outputs/fastqc
  trim_galore -j ${threads} --paired --quality 20 --length 20 --illumina -o outputs/ ${R1} ${R2}
done

#MULTIQC on raw files
multiqc -o outputs/fastqc outputs/fastqc

#FASTQC analysis on trimmed files
for i in outputs/*_1_val_1.fq.gz
do
  R1=${i}
  R2="outputs/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
  fastqc -t ${threads} ${R1} ${R2} -o outputs/fastqc_clean
done

#MULTIQC on trimmed files
multiqc -o outputs/fastqc_clean outputs/fastqc_clean

#Quantification
for i in outputs/*_1_val_1.fq.gz
do
  R1=${i}
  R2="outputs/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
  OUT="outputs/salmon/"$(basename ${i} _1_val_1.fq.gz)"_transcriptome_quant"
  salmon quant -p ${threads} -i refs/transcriptome -l IU -1 ${R1} -2 ${R2} --validateMappings -o ${OUT}
done