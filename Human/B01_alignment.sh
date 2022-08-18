#!/usr/bin/env bash

cd ..
cd ..
mkdir results/alignments

bwa mem  \
data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
data/fastq/mother_R1.fastq.gz \
data/fastq/mother_R2.fastq.gz \
> results/alignments/mother.sam

