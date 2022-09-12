#!/usr/bin/bash

set -e  
cd ../data/fastq/
echo "Running FastQC ..."
fastqc *.fastq

mkdir -p ../../results/fastqc_untrimmed_reads
echo "Saving FastQC results..."
mv *.zip ../../results/fastqc_untrimmed_reads/
mv *.html ../../results/fastqc_untrimmed_reads/

cd ../../results/fastqc_untrimmed_reads/
echo "Unzipping ..."
for filename in *.zip
do
	unzip $filename
done

echo "saving summary..."
cat */summary.txt > ../fastqc_summaries.txt
