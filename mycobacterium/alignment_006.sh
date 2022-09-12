#!/usr/bin/bash

cd ../data/

bwa mem reference/GCF_000013925.1_ASM1392v2_genomic.fna fastq/trimmed_fastq/P7741_R1_trimmed.fastq \
	fastq/trimmed_fastq/P7741_R2_trimmed.fastq > ../results/sam/P7741.aligned.sam
