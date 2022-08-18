#!/usr/bin/bash env
cd ..
cd ..
cd results/alignments
samtools view -bh mother.sorted.sam > mother.bam
