

install.packages("ggplot2")
install.packages("ggrepel")
install.packages("dplyr")
BiocManager::install("EnhancedVolcano")
install.packages("DEGreport")
library(dplyr) 
library(ggplot2) 
library(ggrepel)
library(affy)
library(EnhancedVolcano)
library(magrittr)
library(DEGreport)
theme_set(theme_bw(base_size = 15))#font size
data<-read.csv("/Users/nikolasgiannakis/Desktop/Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/OVcA_Volcano.csv", header=TRUE, sep=",") 
EH_volcano_data  = data.frame(gene = data$Gene, log2_fold_change = data$log2FoldChange, adj_pval = data$neg_log_pval)


EnhancedVolcano(EH_volcano_data, lab = as.character(EH_volcano_data$gene),
                x = "log2_fold_change", y = "adj_pval",
                title = "Healthy vs OvCa")






EnhancedVolcano(data3, lab = rownames(data3),
                x = "log2FoldChange", y = "neg_log_pval",
                selectLab = c("Arg2", "Nr4a1", "Nr4a2", "Gpr35","Akap3","Vegfa"), boxedLabels = TRUE,drawConnectors = TRUE,
                xlim = c(-2, 2),
                xlab  = bquote(~Log[2]~ 'fold change'), 
                ylim = c(0,50),
                ylab=  bquote(~-Log[10]~ adjusted~italic(P)),
                title = 'RvD2 treated vs WT BMDMs',
                pCutoff = 0.05, FCcutoff = 0.58,
                pointSize = 0.9, 
                labSize = 3.0,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8, 
                hline = c(0.05), hlineCol = c('grey0'), hlineType = 'dotted', hlineWidth = 0.5,
                vline = c(-0.58,0.58), vlineCol = c("grey0", "grey0"), vlineType = "dotted", vlineWidth = 0.5,
                legend= c("NS", "Log2 FC", "Adjusted p-value", "Adjusted p-value & Log2 FC"), legendPosition = "bottom", legendIconSize = 3.0, legendLabSize = 10,
                gridlines.major = FALSE, gridlines.minor = FALSE)


## Volcano plot
p= ggplot(data3) +  geom_point(size=0.3,aes(x=log2FoldChange, y=neg_log_pval)) 
+
  ggtitle("Test data") +
  xlab("log2 fold change") + 
  ylab("-log10   adjustedp-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  +
  scale_color_brewer(palette = "Set1")

res_tableOE_ordered <- data2[order(data2$neglogPvalue), ] 
res_tableOE_ordered <- data2[order(data2$p), ]
top20_sigOE_genes <- rownames(res_tableOE_ordered[1:20, ])
top20_sigOE_genes

res_tableOE_ordered$genelabels <- ""
res_tableOE_ordered$genelabels[1:10] <- rownames(res_tableOE_ordered)[1:10]

View(res_tableOE_ordered)
data2<-read.csv("/Users/nikolasgiannakis/Desktop/log2fc_comparison_wt_pparg_nash.csv", header=TRUE, sep=";") 

theplot = ggplot(data2) +  geom_point(size=0.3, aes(x = log2FoldChange_wt_ppargko, y = log2FoldChange_wt_nash)) + scale_color_brewer(palette = "Set1")
theplot

ggplot(data2) +
  geom_point(size=0.3, aes(x = log2FoldChange_wt_ppargko, y = log2FoldChange_wt_nash), colour = TEST5))+
  scale_color_brewer(palette = "Set1") +
  geom_text_repel(
    data2 = dplyr::filter(data2, abs(neglogPvalue) >= 9.64 & abs(LogFC.Ly6Chigh.D1.vs.Ly6Clow.D4.) >= 8.4), #label for the text and filtering based on FC or p-value
    aes(x = LogFC.Ly6Chigh.D1.vs.Ly6Clow.D4., y = neglogPvalue, label = Gene_Symbol),
    size = 6,
    max.iter = 5000) +
  ggtitle("Mov10 overexpression") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) +
  scale_color_brewer(palette = "Set1")

rownames(res_tableOE_ordered)
pp=ggplot(data3, aes(x = log2FoldChange , y = neg_log_pval)) +
  3p=pp + geom_point(size = 0.3)  + #select your x and y axis data, add color (use string characters)
  geom_vline(xintercept = 5) + 
  geom_hline(yintercept = 20) + 
  scale_color_brewer(palette = "Set1") +
  ggtitle("") + #add a title to your plot
  scale_x_continuous("Log2 FC(Healthy vs OvCa Neu)", limits = c(-10, 10)) + #label x axis
  scale_y_continuous(("-Log (P-Value)"), limits = c(0, 20))#data points size
geom_text_repel(
  data = dplyr::filter(data, abs(neglogPvalue) >= 9.64 | abs(D4low.vs.blood.mono) >= 8.4), #label for the text and filtering based on FC or p-value
  aes(label = "Gene_Symbol"),
  size = 6,
  max.iter = 5000) +
  scale_color_brewer(palette = "Set1") +
  ggtitle("") + #add a title to your plot
  scale_x_continuous("Log2 FC(Day 4 Ly6Clow vs Blood monocytes)", limits = c(-10, 10)) + #label x axis
  scale_y_continuous(bquote("-Log (P-Value)"), limits = c(0, 20)) + #label y axis
  geom_hline(yintercept = 0, alpha = 0.5) + 
  geom_vline(xintercept = 0, alpha = 0.5) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none???) )