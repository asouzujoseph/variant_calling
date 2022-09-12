#!/usr/bin/bash

cd ../results/sam
## convert to bam file
samtools view -S -b P7741.aligned.sam > P7741.aligned.bam

# sort BAM file by coordinates
samtools sort -o P7741.aligned.sorted.bam P7741.aligned.bam

## get statistics about the sorted bam file
samtools flagstat P7741.aligned.sorted.bam


