#!/bin/sh

#Calculate (and print to the screen) the number of sequences in mm10 reference genome.

grep -c "^>" ../refs/mm10.fa.gz

#Calculate (and print to the screen) the number of reads in each sample.

echo $(zcat ../inputs/SRR8985047.fastq.gz|wc -l)/4|bc

echo $(zcat ../inputs/SRR8985048.fastq.gz|wc -l)/4|bc

echo $(zcat ../inputs/SRR8985051.fastq.gz|wc -l)/4|bc

echo $(zcat ../inputs/SRR8985052.fastq.gz|wc -l)/4|bc

#Calculate the number of protein-coding genes in your genome.

gunzip -c ../refs/GRCm38.gtf.gz | grep -w "gene" | grep -w "protein_coding" | wc -l