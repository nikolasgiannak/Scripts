##ATAC-seq analysis


## Collecting and filtering Egr2 ChIP-seq data.

## Mapping quality filtering.
 for i in *bam; do sample=$(basename "$i" .bam); samtools view --threads 20 -q 10 -b $i > ${sample}_mapq.bam; done

## Duplicated read removal.
 for i in *_mapq.bam  ; do sample=$(basename "$i" _mapq.bam); samtools rmdup -s $i ${sample}_rmdup.bam; done

## Remove reads fall into blacklist regions.
 for i in *_rmdup.bam; do sample=$(basename "$i" _rmdup.bam); bedtools intersect -v -abam $i -b /data10/genomes/blacklists/mm10-blacklist.bed > ${sample}_bl.bam; done

## Sort bams.
 for i in *bl.bam; do sample=$(basename "$i" .bam); samtools sort $i -o ${sample}_sorted.bam --threads 20; done

## Cleared all temp files and olny kept "_bl.bam"-s. Rename "_bl.bam"s to normal bam.
rm *bl.bam *mapq.bam *rmdup.bam
rm mm_AM_Egr2KO_rep1.bam mm_AM_Egr2KO_rep2.bam mm_AM_WT_rep1.bam mm_AM_WT_rep2.bam
for i in *_bl_sorted.bam; do sample=$(basename "$i" _bl_sorted.bam); mv $i ${sample}.bam; done

## Index bams.
for i in *.bam; do samtools index $i; done

## RPKM normalization.
for i in *bam; do sample=$(basename "$i" .bam); bamCoverage --bam $i --outFileFormat bigwig -o ${sample}.bw --extendReads 120 --binSize 10 --smoothLength 40 --normalizeUsing RPKM --numberOfProcessors 10; done

## Peak calling.
for i in  *.bam; do sample=$(basename "$i" .bam); macs2 callpeak -t $i -f BAM -g 2652783500 --outdir beds/MACS2 -n ${sample} --extsize 120 -q 0.05; done
for i in  *.bam; do sample=$(basename "$i" .bam); macs2 callpeak -t $i -f BAM -g 2652783500 --outdir beds/MACS2_FDR20 -n ${sample} --extsize 120 -q 0.2; done
