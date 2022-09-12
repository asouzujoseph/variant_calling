#!/usr/bin/bash
cd ../data/fastq/

trimmomatic PE -threads 8 P7741_R1.fastq P7741_R2.fastq \
		P7741_R1_trimmed.fastq P7741_R1unp_trimmed.fastq \
		P7741_R2_trimmed.fastq P7741_R2unp_trimmed.fastq \
		SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15

## MINLEN = discard any reads that do not have at least 25 bases remaining after the trimming step
## Sliding window of size 4 that will remove bases if their phred score is below 20
## the three additional numbers (2:40:15) will tell Trimmomatic how to handle sequence matches to the adapter sequence.
