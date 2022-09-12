#!/usr/bin/bash
cd ../results/sam/

## first index the sorted bam file before visulaization wth any tool / software

samtools index P7741.aligned.sorted.bam

## Visualize with tview( a  bcftools component)
samtools tview P7741.aligned.sorted.bam ../../data/reference/GCF_000013925.1_ASM1392v2_genomic.fna
