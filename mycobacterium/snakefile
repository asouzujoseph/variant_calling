## snakemake workflow for variant calling on mycobacterium
#SAMPLES=["P7741","GCA_000013925.2_ASM1392v2"]

rule all:
	input:
		"vcf/P7741_final_variants.vcf"

rule download_reads:
	output:
		"raw_reads/{sample}.fastq.gz"

	shell:
		"""
		wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR333/ERR3335404/{wildcards.sample}.fastq.gz -O {output}
		"""

rule download_genome:
	output:
		"ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna.gz"
	shell:
		"wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Mycobacterium_ulcerans/latest_assembly_versions/GCA_000013925.2_ASM1392v2/GCA_000013925.2_ASM1392v2_genomic.fna.gz -O {output}"

rule unzip_reads:
	input:
		r1="raw_reads/{sample}_R1.fastq.gz",
		r2="raw_reads/{sample}_R2.fastq.gz"
	output:
		r1="raw_reads/{sample}_R1.fastq",
		r2="raw_reads/{sample}_R2.fastq"
	shell:
		"""
		gunzip {input.r1}
		gunzip {input.r2}
		"""

rule unzip_genome:
	input:
		"ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna.gz"
	output:
		"ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna"
	shell:
		"""
		gunzip {input} 
		"""
				
rule fastQC:
	input:
		"raw_reads/{sample}.fastq"
	output:
		"results/fastqc_untrimmed/{sample}_fastqc.html",
		"results/fastqc_untrimmed/{sample}_fastqc.zip"
	shell:
		"""
		fastqc {input}
		mv raw_reads/{wildcards.sample}*fastqc.html {output[0]}
		mv raw_reads/{wildcards.sample}*.zip {output[1]}
		"""

rule trimming:
	input:
		r1="raw_reads/{sample}_R1.fastq",
		r2="raw_reads/{sample}_R2.fastq" 
	output:
		out1="results/trimmed/{sample}_R1_trimmed.fastq",
		out1UN="results/trimmed/{sample}_R1_unpaired_trimmed.fastq",
		out2="results/trimmed/{sample}_R2_trimmed.fastq",
		out2UN="results/trimmed/{sample}_R2_unpaired_trimmed.fastq"
	threads:
		4
	shell:
		"""
		trimmomatic PE -threads {threads} {input.r1} {input.r2} \
		{output.out1} {output.out1UN} \
		{output.out2} {output.out2UN} \
		SLIDINGWINDOW:4:20 MINLEN:25 \
	        ILLUMINACLIP:TruSeq3-PE.fa:2:40:15
		"""

rule index_ref_genome:
	input: 
		"ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna"
	output:
		"ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna.gz.bwt"
	shell:
		"bwa index {input}"

rule mapping:
	input:
		ref="ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna",
		r1="results/trimmed/{sample}_R1_trimmed.fastq",
		r2="results/trimmed/{sample}_R2_trimmed.fastq"
	output:
		"results/sam/{sample}.aligned.sam"
	shell:
		"""
		bwa mem {input.ref} {input.r1} {input.r2} > {output}
		"""

rule bam_Sort:
	input:
		"results/sam/{sample}.aligned.sam"
	output:
		bam="results/sortedBam/{sample}.aligned.sorted.bam",
		stats="results/sortedBam/{sample}.aligned.sorted.bam_stats.txt"
	shell:
		"""		
		samtools view -S -b {input} | samtools sort -o {output.bam}
		samtools flagstat {output.bam} > {output.stats}
		"""

rule readCoverage:
	input: 
		ref="ref_genome/GCA_000013925.2_ASM1392v2_genomic.fna",
		bam="results/sortedBam/{sample}.aligned.sorted.bam"
	output:
		"vcf/{sample}_raw.bcf"
	shell:
		"""
		bcftools mpileup -O b -o {output} -f {input.ref} {input.bam}
		"""

rule callFilterVariants:
	input:
		"vcf/{sample}_raw.bcf"
	output:
		var1="vcf/{sample}_variants.vcf",
		var2="vcf/{sample}_final_variants.vcf"
	shell:
		"""
		bcftools call --ploidy 1 -m -v -o {output.var1} {input}
		vcfutils.pl varFilter {output.var1} > {output.var2}
		"""

