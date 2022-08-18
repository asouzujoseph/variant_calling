#!/usr/bin/env bash
cd ..
cd ..
cd results/alignments

for SAMPLE in mother father son
do
    gatk MarkDuplicates \
    --INPUT "$SAMPLE".rg.bam \
    --OUTPUT "$SAMPLE".rg.md.bam \
    --METRICS_FILE marked_dup_metrics_"$SAMPLE".txt 
done

