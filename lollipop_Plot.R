library(ggplot2)
library(dplyr)

data<-read.csv("/Users/nikolasgiannakis/Desktop/Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/GO_analysis/go_common_3_top_hits.csv", sep=",", header=T)
data2<-read.csv("/Users/nikolasgiannakis/Desktop/Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/GO_analysis/go_common_3_top_hits.csv", sep=",", header=T)


data %>%
  arrange(neg_log_FDR) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  # mutate(name=factor(Gene_ontology, levels=Gene_ontology)) %>%   # This trick update the factor levels
  ggplot(aes(x=Gene_ontology, y=neg_log_FDR)) +
  geom_segment( aes(x=Gene_ontology, xend=Gene_ontology, y=0, yend=neg_log_FDR), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )
# Step 1: Call the pdf command to start the plot

pdf(file = "/Users/nikolasgiannakis/Desktop/Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/GO_analysis/common_3.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
# Horizontal version
ggplot(data, aes(x=Gene_ontology, y=neg_log_FDR)) +
  geom_segment( aes(x=Gene_ontology, xend=Gene_ontology, y=0, yend=neg_log_FDR), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + xlab("Gene ontology") +ylab("Negative log10(FDR)")

dev.off()

pdf(file = "/Users/nikolasgiannakis/Desktop/Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/GO_analysis/common_OvCa_BrCa.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4) # The height of the plot in inches
# Horizontal version
ggplot(data2, aes(x=Gene_ontology, y=neg_log_FDR)) +
  geom_segment( aes(x=Gene_ontology, xend=Gene_ontology, y=0, yend=neg_log_FDR), color="skyblue") +
  geom_point( color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  ) + xlab("Gene ontology") +ylab("Negative log10(p-value)")

dev.off()