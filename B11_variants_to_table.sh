#!/usr/bin/env/bash

cd ..
cd ..
cd results

gatk VariantsToTable \
--variant variants/mother.HC.vcf \
--fields CHROM -F POS -F TYPE -GF GT \
--output variants/mother.HC.table

grep -c "SNP" variants/mother.HC.table
grep -c "INDEL" variants/mother.HC.table
