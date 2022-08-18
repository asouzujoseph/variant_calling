#!/usr/bin/bash env
cd ..
cd ..
cd results/alignments

gatk AddOrReplaceReadGroups \
--INPUT mother.bam \
--OUTPUT mother.rg.bam \
--RGLB lib1 \
--RGPU H0164.2.ALXX140820 \
--RGPL ILLUMINA \
--RGSM mother \
--RGID H0164.2
