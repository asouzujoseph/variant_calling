#!/usr/bin/env bash

cd ../../
mkdir -p results/bqsr
for SAMPLE in mother father son
do
gatk BaseRecalibrator \
--reference A_refs/data/reference/Homo_sapiens.GRCh38.dna.chromosome.20.fa \
--input C_all/alignments/"$SAMPLE".rg.md.bam \
--known-sites A_refs/data/variants/GCF.38.filtered.renamed.vcf \
--known-sites A_refs/data/variants/1000g_gold_standard.indels.filtered.vcf \
--output results/bqsr/"$SAMPLE".recal.table

gatk ApplyBQSR \
--input C_all/alignments/"$SAMPLE".rg.md.bam \
--bqsr-recal-file results/bqsr/"$SAMPLE".recal.table \
--output results/bqsr/"$SAMPLE".recal.bam
done

