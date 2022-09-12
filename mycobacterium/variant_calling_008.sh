#!/usr/bin/bash

# -f flags the path to the reference genome
cd ../results/sam/

## Calculate the read coverage of positions in the genome. 
## the flag -O b tells bcftools to generate a bcf format output file
## -o specifies where to write the output file
## -f specifies the path to the reference genome

bcftools mpileup -O b -o P7741_raw.bcf \
-f ../../data/reference/GCF_000013925.1_ASM1392v2_genomic.fna P7741.aligned.sorted.bam  

## Identify single nucleotide variants 
### The flag --ploidy specifies the number of chromosomes
## -m allows for multiallelic and rare-variant calling
## -v tells the program to output variant sites only
## - o specifies the output file
bcftools call --ploidy 1 -m -v -o P7741_variants.vcf P7741_raw.bcf

#### Filter and report the SNV in variant calling format (VCF)
vcfutils.pl varFilter P7741_variants.vcf > P7741_final_variants.vcf
