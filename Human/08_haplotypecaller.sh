#!/usr/bin/env bash

cd ../../
mkdir -p results/variants
for SAMPLE in mother father son
do
    gatk HaplotypeCaller \
    --reference A_refs/data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
    --input results/bqsr/"$SAMPLE".recal.bam \
    --output results/variants/"$SAMPLE".HC.g.vcf \
    --bam-output results/variants/"$SAMPLE".phased.bam \
    --intervals chr20 \
    --ERC GVCF 
done

