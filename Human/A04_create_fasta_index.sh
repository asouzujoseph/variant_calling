#!/usr/bin/env bash
cd ..
cd ..
cd data/reference 

samtools faidx Homo_sapiens.GRCh38.dna.chromosome.20.fa
gatk CreateSequenceDictionary --REFERENCE Homo_sapiens.GRCh38.dna.chromosome.20.fa

