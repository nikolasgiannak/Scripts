## Gene set enrichment analysis 

Let's use a real example https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197576

How to download the files from ftp https://www.ncbi.nlm.nih.gov/geo/info/download.html

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197576/suppl/

Alternative use GEOquery https://bioconductor.org/packages/release/bioc/html/GEOquery.html

```{bash eval=FALSE}
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197576/suppl/GSE197576_raw_gene_counts_matrix.tsv.gz

csvtk pretty -t  GSE197576_raw_gene_counts_matrix.tsv.gz | less -S

csvtk headers -t   GSE197576_raw_gene_counts_matrix.tsv.gz

csvtk cut -t -f1,2,3,12,13 GSE197576_raw_gene_counts_matrix.tsv.gz| head

csvtk cut -t -f1,2,3,12,13 GSE197576_raw_gene_counts_matrix.tsv.gz > raw_counts.tsv
```

Get csvtk at https://github.com/shenwei356/csvtk

### read the data into R and make a DESeq2 object 

follow the tutorial http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

```{r}
library(dplyr)
library(readr)
library(here)
library(DESeq2)

raw_counts<- read_tsv(here("data/raw_counts.tsv"))

raw_counts_mat<- raw_counts[, -1] %>% as.matrix

head(raw_counts_mat)

rownames(raw_counts_mat)<- raw_counts$gene
head(raw_counts_mat)

```

Make a sample sheet 

```{r}
coldata<- data.frame(condition = c("normoxia", "normoxia", "hypoxia", "hypoxia"))

rownames(coldata)<- colnames(raw_counts_mat)

coldata
```

Make a DEseq2 object

```{r}
all(rownames(coldata) == colnames(raw_counts_mat))

dds <- DESeqDataSetFromMatrix(countData = raw_counts_mat,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "hypoxia", "normoxia"))

res %>%
  as.data.frame() %>%
  arrange((padj), desc(log2FoldChange)) %>%
  head(n=30)


significant_genes<- res %>%
  as.data.frame() %>%
  filter(padj <=0.01, abs(log2FoldChange) >= 2) %>% 
  rownames()


significant_genes
```

## pathway analysis

https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html

### over-representation test

```{r}
library(clusterProfiler)

#convert gene symbol to Entrez ID for 

significant_genes_map<- clusterProfiler::bitr(geneID = significant_genes,
                      fromType="SYMBOL", toType="ENTREZID",
                      OrgDb="org.Hs.eg.db")

head(significant_genes_map)

## background genes are genes that are detected in the RNAseq experiment 
background_genes<- res %>% 
  as.data.frame() %>% 
  filter(baseMean != 0) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pull(gene)


res_df<- res %>% 
  as.data.frame() %>% 
  filter(baseMean != 0) %>%
  tibble::rownames_to_column(var = "gene")

background_genes_map<- bitr(geneID = background_genes, 
                            fromType="SYMBOL", 
                            toType="ENTREZID",
                      OrgDb="org.Hs.eg.db")
```

GO term enrichment 

Gene Ontology(GO) defines concepts/classes used to describe gene function, and relationships between these concepts. It classifies functions along three aspects:

MF: Molecular Function
molecular activities of gene products

CC: Cellular Component
where gene products are active

BP: Biological Process
pathways and larger processes made up of the activities of multiple gene products

GO terms are organized in a directed acyclic graph, where edges between terms represent parent-child relationship.

```{r}
ego <- enrichGO(gene          = significant_genes_map$ENTREZID,
                universe      = background_genes_map$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

library(enrichplot)
barplot(ego, showCategory=20) 
dotplot(ego)
```


H: hallmark gene sets
C1: positional gene sets
C2: curated gene sets
C3: motif gene sets
C4: computational gene sets
C5: GO gene sets
C6: oncogenic signatures
C7: immunologic signatures


```{r}
# install.packages("msigdbr")
library(msigdbr)

m_df <- msigdbr(species = "Homo sapiens")
head(m_df)

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)


table(m_t2g$gs_name)
head(m_t2g)

em <- enricher(significant_genes_map$ENTREZID, TERM2GENE=m_t2g, 
               universe = background_genes_map$ENTREZID )
head(em)
```

### Gene set enrichment analysis

```{r}
## you need all the genes and pre-rank them by p-value
## rank all the genes by signed fold change * -log10pvalue.

res_df<- res_df %>% 
  mutate(signed_rank_stats = sign(log2FoldChange) * -log10(pvalue)) %>%
  left_join(background_genes_map, by= c("gene" = "SYMBOL")) %>%
  arrange(desc(signed_rank_stats))

gene_list<- res_df$signed_rank_stats
names(gene_list)<- res_df$ENTREZID

em2 <- GSEA(gene_list, TERM2GENE=m_t2g)

## change the inf to big numbers
res_df<- res_df %>%
  mutate(negative_log10pvalue = -log10(pvalue)) %>%
  mutate(negative_log10pvalue = ifelse(is.infinite(negative_log10pvalue), 1000, negative_log10pvalue)) %>%
  mutate(signed_rank_stats = sign(log2FoldChange) * negative_log10pvalue)

gene_list<- res_df$signed_rank_stats
names(gene_list)<- res_df$ENTREZID


em2 <- GSEA(gene_list, TERM2GENE=m_t2g)
head(em2)

em2@result %>% View()
```

### visualization 


```{r}
p1<- gseaplot(em2, geneSetID = "HALLMARK_G2M_CHECKPOINT", 
              by = "runningScore", title = "HALLMARK_G2M_CHECKPOINT")

p2 <- gseaplot(em2, geneSetID = "HALLMARK_HYPOXIA", 
               by = "runningScore", title = "HALLMARK_HYPOXIA")

p1/p2
```


important thread on background gene selection https://twitter.com/mdziemann/status/1626407797939384320 by Mark Ziemann

Further reading https://twitter.com/tangming2005/status/1671873310257295360