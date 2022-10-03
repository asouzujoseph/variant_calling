#!/usr/bin/bash

cd ../../

for Sample in mother father son
do
	gatk VariantsToTable \
	--variant results/variants/"$Sample".HC.g.vcf \
	--fields CHROM -F POS -F TYPE -GF GT \
	--output results/variants/"$Sample".HC.table
done
