#!/usr/bin/bash
cd ../data/fastq
fastqc *.fastq
mkdir -p ../../results/fastqc_untrimmed_reads
mv *.zip ../../results/fastqc_untrimmed_reads
mv *.html ../../results/fastqc_untrimmed_reads
cd  ../../results/fastqc_untrimmed_reads/*.zip:

for file in *.zip:
do
unzip $file
done
