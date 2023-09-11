
BiocManager::install("CEMiTool")
library("CEMiTool")


expr = read.csv("/Users/nikolasgiannakis/Desktop/Cabin_GV_Het_KO.csv", sep=",", header=T)
rownames(expr) <- expr[, 1] # place the "Row.names" column in the df rownames
expr[, 1] <- NULL # delete the "Row.names" column
data(expr)


expr = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/CemiTool/H_BrCa_OvCa_COPD_cemitool.csv", sep=",", header=T)
rownames(expr) <- expr[, 1] # place the "Row.names" column in the df rownames
expr[, 1] <- NULL # delete the "Row.names" column
data(expr)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# For comparison H vs BrCa only
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

expr = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/CemiTool/H_BrCa/H_BrCa_cemitool.csv", sep=",", header=T)
rownames(expr) <- expr[, 1] # place the "Row.names" column in the df rownames
expr[, 1] <- NULL # delete the "Row.names" column

sample_annot = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/CemiTool/H_BrCa/sample_annot.csv", sep=",", header=T)

cem <- cemitool(expr, sample_annot, set_beta = 20, sample_name_column="SampleName", class_column="Class")

#SET BETA VALUE TO 7. This is the smallest beta value that the cemitool was not giving me an error,
# while calculating  the cem

# options for cem generation 
# a. filter = TRUE
# b. apply_vst = TRUE
# c. n_genes = 2000
# d. set_beta = 7
# e. min_ngen = 20
# f. merge_similar = FALSE
cem <- cemitool(expr, sample_annot, filter = TRUE, apply_vst = TRUE, n_genes = 2000, set_beta = 7,min_ngen = 20,   merge_similar = FALSE, sample_name_column="SampleName", class_column="Class")

# Here in our example it finds 7 modules
nmodules(cem)

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")


# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]
plots[2]

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)


# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

#For HUMAN samples
int_df <- read.delim("/Users/nikolasgiannakis/Desktop/ligand_receptors_human_for_CeMiTool.txt")
head(int_df)
# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]
plots[2]
plots[3]

# create report as html document
generate_report(cem, directory="./Neu_Report_H_BrCa_v20230703_v3")

# write analysis results into files
write_files(cem, directory="./Neu_Tables_H_BrCa")

# save all plots
save_plots(cem, "all", directory="./Neu_Plots_H_BrCa", force = TRUE)

#####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# For comparison H vs OvCa only
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

expr = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/CemiTool/H_OvCa/H_OvCa_cemitool.csv", sep=",", header=T)
rownames(expr) <- expr[, 1] # place the "Row.names" column in the df rownames
expr[, 1] <- NULL # delete the "Row.names" column

sample_annot = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/CemiTool/H_OvCa/sample_annot.csv", sep=",", header=T)

cem <- cemitool(expr, sample_annot, set_beta = 20, sample_name_column="SampleName", class_column="Class")

# Here in our example it finds 7 modules
nmodules(cem)

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")


# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]
plots[2]
plots[3]
plots[4]
plots[5]

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)


# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

#For HUMAN samples
int_df <- read.delim("/Users/nikolasgiannakis/Desktop/ligand_receptors_human_for_CeMiTool.txt")
head(int_df)
# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]
plots[2]
plots[3]
plots[4]

# create report as html document
generate_report(cem, directory="./Neu_Report_H_OvCa")

# write analysis results into files
write_files(cem, directory="./Neu_Tables_H_OvCa")

# save all plots
save_plots(cem, "all", directory="./Neu_Plots_H_OvCa", force = TRUE)

#####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################









sample_annot = read.csv("/Users/nikolasgiannakis/Desktop/sample_annot.csv", sep=",", header=T)
sample_annot = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/CemiTool/sample_annot.csv", sep=",", header=T)

cem <- cemitool(expr, apply_vst = TRUE, sample_annot, plot_diagnostics = FALSE, filter_pval = 0.05, cor_method = "spearman")



cem <- cemitool(expr, sample_annot, plot_diagnostics = FALSE, filter_pval = 0.05, cor_method = "spearman")

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")


# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]
plots[2]

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)


# perform over representation analysis
cem <- mod_ora(cem, gmt_in)

library(tidyverse)
df <- readRDS("/Users/nikolasgiannakis/Desktop/mouse_lr_pair.rds")
write.csv(df, "ligand_receptors2_mouse.csv", row.names=FALSE)

df <- readRDS("/Users/nikolasgiannakis/Desktop/human_lr_pair.rds")
write.table(df, "ligand_receptors_human.csv", row.names=FALSE)


# read interactions
#int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")


#For HUMAN samples
int_df <- read.delim("/Users/nikolasgiannakis/Desktop/ligand_receptors_human_for_CeMiTool.txt")
head(int_df)
# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]
plots[2]
plots[3]

# create report as html document
generate_report(cem, directory="./Neu_Report")

# write analysis results into files
write_files(cem, directory="./Neu_Tables")

# save all plots
save_plots(cem, "all", directory="./Neu_Plots", force = TRUE)




#I will check onlyt the Het modules
expr = read.csv("/Users/nikolasgiannakis/Desktop/cabin_het_ko.csv", sep=",", header=T)
rownames(expr) <- expr[, 1] # place the "Row.names" column in the df rownames
expr[, 1] <- NULL # delete the "Row.names" column

cem <- cemitool(expr[,2:5], apply_vst = TRUE)

nmodules(cem)

#Start analysis from scratch
expr_new = read.csv("/Users/nikolasgiannakis/Desktop/Cabin_GV_Het_KO.csv", sep=",", header=T)
rownames(expr_new) <- expr_new[, 1] # place the "Row.names" column in the df rownames
expr_new[, 1] <- NULL # delete the "Row.names" column


sample_annot_new = read.csv("/Users/nikolasgiannakis/Desktop/sample_annot.csv", sep=",", header=T)

cem_new <- cemitool(expr_new, sample_annot_new, plot_diagnostics = FALSE, filter_pval = 0.05, cor_method = "spearman")

# Here in our example it finds 7 modules
nmodules(cem_new)

find_modules(cem_new)

# generate heatmap of gene set enrichment analysis
cem_new <- mod_gsea(cem_new)
cem <- plot_gsea(cem_new)
show_plot(cem_new, "gsea")
mod_gsea(cem_new)


cem_new <- plot_profile(cem_new)
plots <- show_plot(cem_new, "profile")
plots[1]
plots[2]
plots[3]
plots[4]

library(ggplot2)
interactions_data(cem_new) <- int_df # add interactions
cem_new <- plot_interactions(cem_new) # generate plot
plots_inter <- show_plot(cem_new, "interaction") # view the plot for the first module
plots_inter[1]


# create report as html document
generate_report(cem_new, directory="./Report_new_analysis")

# write analysis results into files
write_files(cem_new, directory="./Tables_new_analysis")

# save all plots
save_plots(cem_new, "all", directory="./Plots_new_analysis")


# Repeat the whole analysis but instead of spearman corrleation, I use pearson, and p_val_cut_off = 0.05

cem_new <- cemitool(expr_new, sample_annot_new, plot_diagnostics = FALSE, min_ngen = 20, filter_pval = 0.05, cor_method = "spearman")

# Here in our example it finds 7 modules
nmodules(cem_new)

find_modules(cem_new)

# generate heatmap of gene set enrichment analysis
cem_new <- mod_gsea(cem_new)
warnings()
cem_new <- plot_gsea(cem_new)
show_plot(cem_new, "gsea")
mod_gsea(cem_new)


cem_new <- plot_profile(cem_new)
plots <- show_plot(cem_new, "profile")
plots[1]
plots[2]
plots[3]
plots[4]
plots[5]

library(ggplot2)
interactions_data(cem_new) <- int_df # add interactions
cem_new <- plot_interactions(cem_new) # generate plot
plots_inter <- show_plot(cem_new, "interaction") # view the plot for the first module
plots_inter[1]


# create report as html document
generate_report(cem_new, directory="./Report_new_analysis2")

# write analysis results into files
write_files(cem_new, directory="./Tables_new_analysis2")

# save all plots
save_plots(cem_new, "all", directory="./Plots_new_analysis2")


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# For Lance Hira GV Het vs KO
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################

expr = read.csv("/Users/nikolasgiannakis/Desktop/For_Lance/Hira_GV/Hira_GV.csv", sep=",", header=T)
rownames(expr) <- expr[, 1] # place the "Row.names" column in the df rownames
expr[, 1] <- NULL # delete the "Row.names" column

sample_annot = read.csv("/Users/nikolasgiannakis/Desktop/For_Lance/Hira_GV/sample_annot.csv", sep=",", header=T)

# At the beginning I run this command, but it couldn't find beta parameter.
# For this reason I run a diagnostics report to check the Beta X R2 plot

cem <- cemitool(expr, sample_annot, plot_diagnostics = TRUE, filter_pval = 0.05)

diagnostic_report(cem_new, force = TRUE)

# As at beta =20, it approximates 0.8 of R2, I will use this as beta
# beta = 20
# pval = 0.05
#correlation method
# It uses here pearson correlation, not spearman, unless set otherwise
cem <- cemitool(expr, sample_annot, set_beta = 20, filter_pval = 0.05)

# Here in our example it finds 7 modules
nmodules(cem)

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")


# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]
plots[2]
plots[3]


# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)


# perform over representation analysis
cem <- mod_ora(cem, gmt_in)


# FOR MOUSE SAMPLES
int_df <- read.delim("/Users/nikolasgiannakis/Desktop/ligand_receptors_mouse_for_CeMiTool.txt")
head(int_df)
# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]
plots[2]


# create report as html document
generate_report(cem, directory="./Hira_GV_pval_0_05_report")

# write analysis results into files
write_files(cem, directory="./Tables_Hira_GV_pval_0_05")

# save all plots
save_plots(cem, "all", directory="./Plots_Hira_GV_pval_0_05", force = TRUE)



#I will repeat the analysis with pval = 0.1, which is the proposed
# I tried with out setting beta, but gave error.
cem <- cemitool(expr, sample_annot, plot_diagnostics = TRUE)

# I repeat with beta = 20
cem <- cemitool(expr, sample_annot, set_beta = 20)
# Here in our example it finds 7 modules
nmodules(cem)

# generate heatmap of gene set enrichment analysis
cem <- mod_gsea(cem)
cem <- plot_gsea(cem)
show_plot(cem, "gsea")


# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]
plots[2]
plots[3]


# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)


# perform over representation analysis
cem <- mod_ora(cem, gmt_in)


# FOR MOUSE SAMPLES
int_df <- read.delim("/Users/nikolasgiannakis/Desktop/ligand_receptors_mouse_for_CeMiTool.txt")
head(int_df)
# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]
plots[2]


# create report as html document
generate_report(cem, directory="./Hira_GV_pval_0_1_report")

# write analysis results into files
write_files(cem, directory="./Tables_Hira_GV_pval_0_1")

# save all plots
save_plots(cem, "all", directory="./Plots_Hira_GV_pval_0_1", force = TRUE)



#####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################



