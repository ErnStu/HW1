#!/bin/sh

#Download the reference genome in FASTA format and put it into references directory.

wget -O ../refs/mm10.fa.gz "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz"

#Download the reference transcriptome in FASTA format and put it into references directory.

wget -O ../refs/mm10_transcriptome.fa.gz "https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"

#Download GTF/GFF3 file for mm10 reference genome.

wget -O ../refs/GRCm38.gtf.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz

#Download raw FASTQ files.

prefetch SRR8985047 SRR8985048 SRR8985051 SRR8985052

fastq-dump SRR8985047 SRR8985048 SRR8985051 SRR8985052

gzip SRR8985047 SRR8985048 SRR8985051 SRR8985052