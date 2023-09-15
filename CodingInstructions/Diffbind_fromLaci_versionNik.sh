library(DiffBind); library(tidyverse); library(hexbin)

## Set working directory.
setwd("/Users/nlab/Desktop/hlszlaszlo/Projects/0002-Petros/20180206_ATAC-seq/analysis_hlsz/MACS2/Consensus_set_Andreas_WT_only_H_vs_L_v2/")

## Samples.
samples = read_csv("csv/samples3.csv")

## All consensus.
ATAC.dba = dba(sampleSheet = samples, peakFormat = "narrow")
ATAC.dba.count = dba.count(ATAC.dba, bLog = T, score = DBA_SCORE_RPKM, filter = 10)

## Consensus set.
consensus1 = data.frame(dba.peakset(ATAC.dba.count, bRetrieve = T))
write.table(consensus1, file = 'MDM_ATAC_consensus_WT.tsv', sep = '\t', quote = F, row.names = FALSE, col.names = T)
write.table(consensus1[, 1:3], file = 'MDM_ATAC_consensus_WT.bed', sep = '\t', quote = F, row.names = FALSE, col.names = F)

## D1H D4L
samples_sub = samples[c(1, 2, 5, 7, 8), ]
ATAC.dba = dba(sampleSheet = samples_sub, peakFormat = "narrow")
ATAC.dba.count = dba.count(ATAC.dba, bLog = T, score = DBA_SCORE_TMM_READS_FULL, filter = 10)

## Consensus set.
consensus1 = data.frame(dba.peakset(ATAC.dba.count, bRetrieve = T))
write.table(consensus1, file = 'MDM_ATAC_consensus_D1H_D4L.tsv', sep = '\t', quote = F, row.names = FALSE, col.names = T)
write.table(consensus1[, 1:3], file = 'MDM_ATAC_consensus_D1H_D4L.bed', sep = '\t', quote = F, row.names = FALSE, col.names = F)

## Set contrast.
ATAC.dba.count.contrast = dba.contrast(ATAC.dba.count,
                                       group1 = ATAC.dba.count$masks$MDM_WT_D1_high,
                                       group2 = ATAC.dba.count$masks$MDM_WT_D4_low,
                                       name1 = "D1_High",
                                       name2 = "D4_Low")

ATAC.dba.count.analyze = dba.analyze(ATAC.dba.count.contrast, bCorPlot = FALSE, method = DBA_DESEQ2, bParallel = TRUE, bFullLibrarySize = FALSE)
ATAC.dba.count.analyze.DB = as.tibble(dba.report(ATAC.dba.count.analyze, th = 1))

maframe1 = ATAC.dba.count.analyze.DB %>%
    mutate(sig = ifelse(FDR < 0.05 & abs(Fold) > 1.5 , "Sig", "NotSig"),
           direction = ifelse(sig == "Sig" & Fold > 1.5, "D1_H_spec", ifelse(sig == "Sig" & Fold < -1.5, "D4_L_spec", "Common")),
           facet1 = "D1_H", facet2 = "D4_L")

colnames(maframe1) = c("seqnames",
                       "start",
                       "end",
                       "width",
                       "strand",
                       "Conc",
                       "Conc_1",
                       "Conc_2",
                       "Fold",
                       "p.value",
                       "FDR",
                       "sig",
                       "direction",
                       "facet1",
                       "facet2")

write.table(maframe1, file = paste("DB_", "D1_H", "__vs__", "D4_L", "_all.tsv", sep = ""), sep = '\t', quote = F, row.names = FALSE, col.names = TRUE)
write.table(filter(maframe1, direction == "Common"), file = paste("DB_", "D1_H", "__vs__", "D4_L", "_common.bed", sep = ""), sep = '\t', quote = F, row.names = FALSE, col.names = FALSE)
write.table(filter(maframe1, direction == "D1_H_spec"), file = paste("DB_", "D1_H", "__vs__", "D1_H", "_spec.bed", sep = ""), sep = '\t', quote = F, row.names = FALSE, col.names = FALSE)
write.table(filter(maframe1, direction == "D4_L_spec"), file = paste("DB_", "D1_H", "__vs__", "D4_L", "_spec.bed", sep = ""), sep = '\t', quote = F, row.names = FALSE, col.names = FALSE)

mergedf1 = data.frame()
mergedf1 = rbind(mergedf1, maframe1)
labels1 = data.frame()
labels1[1, "facet1"] = "D4_L"
labels1[1, "WT_count"] = as.numeric(table(maframe1$direction)["D1_H_spec"])
labels1[1, "KO_count"] = as.numeric(table(maframe1$direction)["D4_L_spec"])

mergedf1$facet1 = factor(mergedf1$facet1, levels = c("D1_H"))

mergedf1$direction = factor(mergedf1$direction, levels = c("D1_H_spec", "D4_L_spec"))

theme_set(theme_bw(base_size = 12))
gg1 = ggplot(filter(mergedf1, sig == "NotSig"), aes(x = Conc, y = Fold)) +
    geom_hex(alpha = 0.8) +
    geom_point(data = filter(mergedf1, sig == "Sig"), aes(x = Conc, y = Fold, color = direction), size = 0.5) +
    scale_x_continuous("Normalized Read Count", limits = c(-1, 11), breaks = seq(0, 10, by = 2)) +
    scale_y_continuous("log(Fold Change)", limits = c(-6, 6), breaks = seq(-6, 6, by = 2)) +
    scale_fill_gradient2(low = "#f0f0f0", mid = "#f0f0f0", high = "#2171b5", midpoint = 0.05) +
    scale_color_manual(values = c(“#7eb837", "#377eb8"),
                       name="",
                       breaks = c("WT_spec", "KO_spec"),
                       labels = c(paste("Wild type specific (n = ", as.numeric(table(maframe1$direction)["WT_spec"]), ")", sep = ""),
                                  paste("BACH1-KO specific (n = ", as.numeric(table(maframe1$direction)["KO_spec"]), ")", sep = ""))) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 1.5, linetype = "dotted") +
    geom_hline(yintercept = -1.5, linetype = "dotted") +
    theme(legend.position = "none") +
    facet_wrap(~facet1, ncol = 1) +
    annotate(geom = "text", x = -1, y = 5, label = labels1$WT_count, color = "#7eb837", hjust = 0) +
    annotate(geom = "text", x = -1, y = -5, label = labels1$KO_count, color = “#377eb8”, hjust = 0)

ggsave(gg1, filename = “nikolas_DifferentialBinding_WT_D1H_vs_D4L_hex.pdf”, units = "in", width = 3, height = 3)


