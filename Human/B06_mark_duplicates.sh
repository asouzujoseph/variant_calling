#!/usr/bin/env bash
cd ..
cd ..

cd results/alignments

gatk MarkDuplicates \
--INPUT mother.rg.bam \
--OUTPUT mother.rg.md.bam \
--METRICS_FILE marked_dup_metrics_mother.txt 

