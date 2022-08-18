#!/usr/bin/env bash
cd ..
cd ..
cd results/alignments

for SAMPLE in mother father son
do
    samtools index "$SAMPLE".rg.md.bam
done  

