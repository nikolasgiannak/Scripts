#!/bin/bash

# Set the input and output directories
indir="input_data"
outdir="output_data"

# Create the output directory if it does not already exist
if [ ! -d "$outdir" ]; then
  mkdir "$outdir"
fi

# Loop through each sample in the input directory
for sample in $indir/*
do
  # Set the sample name and read files
  sample_name=$(basename "$sample")
  read1="$sample/${sample_name}_R1.fastq.gz"
  read2="$sample/${sample_name}_R2.fastq.gz"

  # Trim adapters from the reads
  trimmomatic PE "$read1" "$read2" "${sample_name}_R1_trimmed.fastq.gz" "${sample_name}_R1_unpaired.fastq.gz" "${sample_name}_R2_trimmed.fastq.gz" "${sample_name}_R2_unpaired.fastq.gz" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10

  # Perform QC analysis on the trimmed reads
  fastqc "${sample_name}_R1_trimmed.fastq.gz" -o "$outdir"
  fastqc "${sample_name}_R2_trimmed.fastq.gz" -o "$outdir"

  # Align the trimmed reads to the reference genome using hisat2
  hisat2 -x /path/to/reference/genome -1 "${sample_name}_R1_trimmed.fastq.gz" -2 "${sample_name}_R2_trimmed.fastq.gz" -S "${sample_name}.sam"

  # Convert the SAM file to BAM and sort it
  samtools view -bS "${sample_name}.sam" | samtools sort - "${sample_name}.sorted"

  # Index the sorted BAM file
  samtools index "${sample_name}.sorted.bam"

  # Count the number of reads aligned to each gene in the reference genome
  htseq-count -f bam -r pos -s no -t gene -i ID "${sample_name}.sorted.bam" /path/to/annotation.gtf > "${sample_name}_counts.txt"
done

# Concatenate the count tables for all samples into a single file
cat *_counts.txt > counts.txt