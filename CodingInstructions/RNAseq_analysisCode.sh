Code for RNA-seq analysis WT vs KO  (3 replicates each)

#login to server, changes upon different user, here I include my domain
ssh ngiannakis@ngsdeb.med.unideb.hu -p 11240

#enter password

#print path of the working directory starting from the root
pwd

#create directory where you will make the analysis
mkdir RNA_seq_analysis

#create txt file, where you create links of the fastq files you will analyse and name them, include the whole name of the path 
vi getSymlink.txt

#inside the txt file, and in case of comparison between WT to KO samples (3 replicates for each condition), write the following piece of code, creates links  of the fast files and names them
ln -s /path_tof_file/…/file1.fastq.gz mm_WT_rep1.fastq.gz
ln -s  /path_tof_file/…/file2.fastq.gz mm_WT_rep2.fastq.gz
ln -s  /path_tof_file/…/file3.fastq.gz mm_WT_rep3.fastq.gz
ln -s  /path_tof_file/…/file4.fastq.gz mm_KO_rep1.fastq.gz
ln -s  /path_tof_file/…/file5.fastq.gz mm_KO_rep2.fastq.gz
ln -s  /path_tof_file/…/file6.fastq.gz mm_KO_rep3.fastq.gz

#save and  quit

# alignment with HISAT2, output is a bam file
for i in *fastq.gz; do sample=$(basename "$i" fastq.gz); hisat2 -x /data10/genomes/index/mm10 --dta --known-splicesite-infile /data10/genomes/Mus_musculus/UCSC/mm10/Sequence/Hista2Index/genes.ss --met-file ${sample}hisat.stat -p 40 -U $i | samtools sort -@ 40 -o ${sample}hisat.bam -; done

## Estimate transcript abundance by sample.
#create directory for string tie
#gtf files have been previously downloaded
mkdir stringtie
ls *.bam | xargs basename -a | cut -f1 -d "." | parallel mkdir stringtie/{}
ls -d stringtie/* | parallel echo {/} | parallel --progress 'stringtie -G /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -e -B -o ./stringtie/{}/transcripts.gtf -A ./stringtie/{}/gene_abundances.tsv -C ./stringtie/{}/cov_ref.gtf ./{}.hisat.bam &> ./stringtie/{}/stringtie.log'

cd stringtie
find | grep transc | xargs stringtie --merge -G /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o merge_transcripts.gtf

#change directory to stringtie
cd stringtie

#quantify reads
featureCounts -a merge_transcripts.gtf -o counts_for_DEGs.tsv -T 8 -O ../mm_WT_rep1.hisat.bam ../mm_WT_rep2.hisat.bam ../mm_WT_rep3.hisat.bam ../mm_KO_rep1.hisat.bam ../mm_KO_rep2.hisat.bam ../mm_KO_rep3.hisat.bam

#retrieve specific columns from a file and write the result to a new file
cut -f 1,7,8,9,10,11,12 counts_for_DEGs.tsv > gene_count_matrix.tsv

#create sample annotation file needed for the Deseq2 analysis
vi sampleannotation.csv

#write the following
file,condition,type
CTRL_1,CTRL,single-end
CTRL_2,CTRL,single-end
CTRL_3,CTRL,single-end
KO_1,KO,single-end
KO_2,KO,single-end
KO_3,KO,single-end

#save and quite , :wq!

#open R for downstream analysis, DEG analysis
R

#get directory, set directory
getwd()
setwd()

#load libraries needed for downstream analysis, in case you haven’t installed packages before do it now and then install libraries
#install.packages(“ “) is the command
library("DESeq2"); library("BiocParallel"); library("IHW"); library("ggplot2"); library("dplyr"); library("sqldf"); library("ggrepel"); library("tidyverse"); library("stringr")

## register CPUs to do parallel calculations.
register(MulticoreParam(4))

##set count data and sample annotation file names.
count1 = "gene_count_matrix.tsv"
annot1 = "sampleannotation.csv"

## Read counts and column data.
cts_all = as.matrix(read.csv(count1, sep = "\t", row.names = "Geneid"))

vec = c("CTRL", "CTRL", "CTRL",
        "KO", "KO", "KO")

annotation_gtf = gtf_annot = read_tsv("../merge_transcripts.gtf", skip = 2, col_names = F)
annotation_gtf = annotation_gtf %>% filter(X3 == "transcript") %>% select(X9)
annotation_gtf = str_split_fixed(annotation_gtf$X9, ";", 4)[, c(1, 3)]
annotation_gtf = data.frame(annotation_gtf)
annotation_gtf = annotation_gtf[, c(1, 2)]
annotation_gtf$X1 = gsub("gene_id ", "", annotation_gtf$X1)
annotation_gtf$X1 = gsub("\"", "", annotation_gtf$X1)
annotation_gtf$X2 = gsub("gene_name ", "", annotation_gtf$X2)
annotation_gtf$X2 = gsub("\"", "", annotation_gtf$X2)
annotation_gtf = annotation_gtf[!duplicated(annotation_gtf), ]

coldata_all = read.csv(annot1, row.names = 1)

    a1 = vec[1]
    b1 = vec[3 + 1]
    coldata = coldata_all %>% filter(condition %in% c(a1, b1))
    labels = paste(rep(c(a1, b1), each = 3), seq(1, 3, by = 1), sep = "_")
    rownames(coldata) = labels

    cts = cts_all
#create a DESeqDataSet object with the count and cloumn data, with design to compare WT with ndx1-4.
dds = DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)

#low count filter.
dds = dds[rowSums(counts(dds)) > 1, ]

#run DESeq analysis.
dds = DESeq(dds)

#calcuate differentially expressed genes with FDR < 0.1 and log foldchange > |0.5|
resIHW01 = results(dds, alpha = 0.1, lfcThreshold = 0.5)
resIHW01 = lfcShrink(dds = dds, coef = 2, res = resIHW01)

#plot DESeq2 MA plot.
pdf(paste("DEGs-MAplot-IHW_", a1, "_vs_", b1, ".pdf", sep = ""), width = 6, height = 4, paper = 'special')
plotMA(resIHW01, ylim = c(-6, 6))
dev.off()

#subset significant genes based on adjusted p-values.
resSig = subset(resIHW01, padj < 0.1, abs(log2FoldChange) > 0.5)

#normalized counts.
norm_counts = counts(dds, normalized=TRUE) %>% as.data.frame() %>% rownames_to_column(., "rownames")


#reform  output and annotate genes.
df1 = data.frame(Gene = rownames(resIHW01), baseMean = resIHW01$baseMean, log2FoldChange = resIHW01$log2FoldChange, stat = resIHW01$stat, pvalue = resIHW01$pvalue, padj = resIHW01$padj)
df1 = mutate(df1, sig = ifelse((df1$padj < 0.1 & abs(df1$log2FoldChange) > 0.5), "Sig", "Not Sig"), direction = ifelse(sig == "Sig" & log2FoldChange >= 0.5, "UP", ifelse(sig == "Sig" & log2FoldChange <= -0.5, "DOWN", "Unchanged")))
cts2 = data.frame(cts, id = rownames(cts))
out1 = sqldf('select * from df1 left outer join cts2 on df1.Gene = cts2.id')[, 1:14]
out2 = sqldf('select * from out1 left outer join annotation_gtf on out1.Gene = annotation_gtf.X1')[, c(1:14, 16)]
out2 = out2 %>% left_join(norm_counts, by = c("Gene" = "rownames"))


#write table with DEGs
write.table(out2, file = paste("DEGs-list-IHW_", a1, "_vs_", b1, “.csv”, sep = ""), col.names = T, row.names = F, quote = F, sep = "\t")

#plot individual genes.
#nikolas_gene.
d <- plotCounts(dds, gene = "MSTRG.number”, intgroup = "condition", returnData = TRUE)
d$condition = factor(d$condition, levels = c("WT", “KO”))
theme_set(theme_bw())
gg1 = ggplot(d, aes(x = condition, y = count)) +
geom_point(position = position_jitter(w = 0.1, h = 0)) + 
ggtitle(“NIKOLAS-GENE gene expression") + 
theme(legend.title = element_blank(), legend.text = element_text(size = 14, face = "bold"), axis.text.x = element_text(size = 14, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title = element_text(size = 14, face = "bold")) 

# ggsave(plot = gg1, file = "DEGs-nikolas_woohoo_cambridge.pdf", width = 4, height = 4)










