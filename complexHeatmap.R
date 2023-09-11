#Download package and load library ComplexHeatmap
#It can be done in two ways, either the classic install package and load library or
#use dev tools and install github repository

#Version 1 of installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#Version 2 of installation
install.packages("devtools")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

#Load data in Rstudio

data<-read.csv("/Users/nikolasgiannakis/Desktop/Neu_DEGs_expression.csv", sep=",", header=T)

data2<-read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/Big_Heatmap_all_Neu_DEGs/Neu_norm_exprs_values.csv", header=TRUE, sep=",") 

data3 = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/Shaul_2016_N1_N2/Shaul_N1_N2_genes_H_BrCa_OvCa_clean.csv", header=TRUE, sep=",") 
#Load the metadata for annotation 
metadata<-read.csv("/Users/nikolasgiannakis/Desktop/Neu_DEGs_metadata.csv", sep=",", header=T)

metadata2<-read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/RNAseq_data_heatlhy_ov_br_copd/My_analysis/Big_Heatmap_all_Neu_DEGs/neu_norm_exprs_metadata.csv", sep=",", header=T)

metadata3 = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/Shaul_2016_N1_N2/Shaul_N1_N1_metadata_H_BrCa_OvCa.csv",  header=TRUE, sep=",") 







##########################################
##########################################
ann = data.frame(metadata$Tissue)
colnames(ann) = c("Tissue")
colours = list( "Tissue" = c("Healthy" = "deepskyblue2", "OvCa" = "firebrick3",  "BrCa" = "coral", "COPD" = "darkgreen"))
#  ,   "dodgerblue3", "firebrick3"
#   "Day" = c("1" = "gold", "2" = "darkgreen", "4" = "darkorange"))
colAnn = HeatmapAnnotation(df = ann,
                           which = 'col',
                           col = colours,
                           annotation_width = unit(c(1, 4), 'cm'),
                           gap = unit(2, 'mm'))

##########################################
##########################################
ann = data.frame(metadata2$Condition)
colnames(ann) = c("Condition")
colours = list( "Condition" = c("Healthy" = "deepskyblue3", "OvCa" = "brown4",  "BrCa" = "coral3", "COPD" = "darkgreen"))
#  ,   "dodgerblue3", "firebrick3"
#   "Day" = c("1" = "gold", "2" = "darkgreen", "4" = "darkorange"))
colAnn = HeatmapAnnotation(df = ann,
                           which = 'col',
                           col = colours,
                           annotation_width = unit(c(1, 4), 'cm'),
                           gap = unit(2, 'mm'))






##########################################
##########################################
#Set the numeric matrix that will be used for the heatmap visualization

mdata<-data.matrix(data[,2:ncol(data)])
head(mdata)
dim(mdata)
#Set the row names argument to be the column that has the gene names.
#Here in this test it is the 1st column of the data 

rownames(mdata) = data[,1]


scaled_data = t(scale(t(mdata[, 1:32])))
rownames(scaled_data) = data[,1]
###############################################

mdata2<-data.matrix(data2[,2:ncol(data2)])
head(mdata2)
dim(mdata2)
#Set the row names argument to be the column that has the gene names.
#Here in this test it is the 1st column of the data 

rownames(mdata2) = data2[,1]


scaled_data2 = t(scale(t(mdata2[, 1:32])))
rownames(scaled_data2) = data2[,1]


##########################################
##########################################

#Plot a basic Heatmap

Heatmap(mdata)

#Create a color key  and plot a basic heatmap

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue3", "white", "darkgreen"))
col_fun(seq(-3, 3))

col_fun = colorRamp2(c(-3, 0, 3), c("darkblue", "white", "darkgreen"))
col_fun(seq(-3, 3))


Heatmap(scaled_data, name = "mat", col = col_fun)

#set heatmap borders

Heatmap(scaled_data, name = "mat", border_gp = gpar(col = "black", lty = 2))

#or alternatively set cells borders

Heatmap(mdata, name = "mat", rect_gp = gpar(col = "white", lwd = 2))

# hide column dendrogram

Heatmap(scaled_data, name = "mat", show_column_dend = FALSE) 

#manually set the dimensions of the dendrogram related to the heamap's clustering

Heatmap(mdata, name = "mat",col=col_fun, column_dend_height = unit(2, "cm"), 
        row_dend_width = unit(2, "cm"))


#################################################
#################################################
#################################################
##########CLUSTERING#############################

#NO rows' clustering
Heatmap(scaled_data, name = "mat", col=col_fun, cluster_rows = FALSE) # turn off row clustering

#Hierarchical clustering based on pre-defined distance metrics
#   "pearson" , "spearman" , "kendall"

Heatmap(scaled_data, name = "mat", col=col_fun, clustering_distance_rows = "pearson",
        column_title = "pre-defined distance method (1 - pearson)")

#change fontsize for gene names column
Heatmap(scaled_data, name = "z-score", col=col_fun, row_names_gp = gpar(fontsize = 20))

#################################################
#################################################
#################################################
#################################################



###########SPLIT#################################
##################<<<<K-MEANS>>>>###############
#####################CLUSTERING#################

#CLUSTER ROWS
Heatmap(scaled_data, name = "z-score",col=col_fun, row_km = 2)

#CLUSTER COLUMNS
Heatmap(mdata, name = "z-score",col=col_fun, column_km = 3)

#SPLIT BY COLUMN AND ROW
Heatmap(scaled_data, name = "z-score",col=col_fun, row_km = 3, column_km = 3, row_title = "Cluster %s")

#IF we want to name the ordered clusters manually and in this occasion
# for 3 clusters we do the following:

Heatmap(mdata, name = "z-score",col=col_fun, row_km = 3, column_km = 3,
        row_title= c("Fabulous", "You bitch", "Miss you"),
        row_title_rot = 0)

####################################################
##########END OF SPLIT#############################
####################################################

####################################################
#########<<<<<<<<<<<ADD>>>>>>>>>>>##################
#############<<<<<ANNOTATION>>>>>>##################

Heatmap(scaled_data, name = "z-score",col=col_fun,
        top_annotation = colAnn,
        column_km = 4,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:4),
                                                         labels = c("Cluster 1 ", "Cluster 2 ", "Cluster 3","NikoCluster"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10))),
        row_km = 4)

Heatmap(scaled_data, name = "z-score",col=col_fun,
        top_annotation = colAnn,
        column_km = 4,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 1:4),
                                                         labels = c("Cluster 1 ", "Cluster 2 ", "Cluster 3","Cluster 4","Cluster 5 ", "Cluster 6 ", "Cluster 7"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10))),
        row_km = 7)


#Remove dendrograms from heatmap and column names
set.seed(3)
ht = Heatmap(scaled_data2, name = "z-score",col=col_fun, show_column_dend = TRUE, show_row_dend = FALSE, column_title = "", use_raster = FALSE,
             show_row_names=FALSE,
             show_column_names = FALSE,
             top_annotation = colAnn,
             column_km = 5,
             left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "white"),labels_rot = 0,
                                                              labels = c("Cluster 1 ", "Cluster 2 ", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"), 
                                                              labels_gp = gpar(col = "black", fontsize = 6))),
             row_km = 6,
             width = ncol(scaled_data2)*unit(3, "mm"), 
             height = nrow(scaled_data2)*unit(0.01, "mm"))
ht = draw(ht)
col_fun = colorRamp2(c(-2, 0, 2), c("darkblue", "white", "brown4"))
# col_fun(seq(-1, 1))

r.dend <- row_dend(ht)  #Extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

> # loop to extract genes for each cluster.
  for (i in 1:length(row_order(ht))){
    if (i == 1) {
      clu <- t(t(row.names(scaled_data2[row_order(ht)[[i]],])))
      out <- cbind(clu, paste("cluster", i, sep=""))
      colnames(out) <- c("GeneID", "Cluster")
    } 
    else {
      clu <- t(t(row.names(scaled_data2[row_order(ht)[[i]],])))
      clu <- cbind(clu, paste("cluster", i, sep=""))
      out <- rbind(out, clu)
    }
  }

out

#export
write.table(out, file= "gene_clusters.txt", sep="\t", quote=F, row.names=FALSE)





#No clustering the columns
Heatmap(scaled_data2, name = "z-score",col=col_fun, show_row_names=FALSE,
        show_column_names = FALSE,
        top_annotation = colAnn,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "white"),
                                                         labels = c("Cluster 1 ", "Cluster 2", "Cluster 3", "cluster 4"), 
                                                         labels_gp = gpar(col = "black", fontsize = 10))),
        row_km = 4,
        cluster_columns = FALSE)

Heatmap(scaled_data, name = "z-score",col=col_fun,
        top_annotation = colAnn,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                         labels = c("Cluster 1 ", "Cluster 2", "Cluster 3", "cluster 4"), 
                                                         labels_gp = gpar(col = "white", fontsize = 10))),
        row_km = 4,
        cluster_columns = FALSE)


Heatmap(scaled_data, name = "z-score",col=col_fun,
        top_annotation = colAnn,
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:4),
                                                         labels = c("Cluster 1 ", "Cluster 2", "Cluster 3", "cluster 4", "Cluster 5", "Cluster 6"), 
                                                         labels_gp = gpar(col = "white", fontsize = 7))),
        row_km = 6,
        cluster_columns = FALSE
)



#############################################################################
#############################################################################
#################      S H A U L N1 and N2 genes     ########################
#############################################################################


data3 = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/Shaul_2016_N1_N2/Shaul_N1_N2_genes_H_BrCa_OvCa_clean.csv", header=TRUE, sep=",") 
metadata3 = read.csv("/Users/nikolasgiannakis/Desktop/AAA_Ovarian_cancer_p/Shaul_2016_N1_N2/Shaul_N1_N1_metadata_H_BrCa_OvCa.csv",  header=TRUE, sep=",") 

rownames(data3) = data3$Gene


ann = data.frame(metadata3$Condition)
colnames(ann) = c("Condition")
colours = list( "Condition" = c("Healthy" = "deepskyblue2", "OvCa" = "firebrick3",  "BrCa" = "coral"))

#     "Healthy" = "deepskyblue3", "OvCa" = "brown4",  "BrCa" = "coral3"
#  ,   "dodgerblue3", "firebrick3"
#   "Day" = c("1" = "gold", "2" = "darkgreen", "4" = "darkorange"))
colAnn = HeatmapAnnotation(df = ann,
                           which = 'col',
                           col = colours,
                           annotation_width = unit(c(1, 4), 'cm'),
                           gap = unit(2, 'mm'))


mdata3 <- data.matrix(data3[,2:ncol(data3)])
head(mdata3)
dim(mdata3)
rownames(mdata3) = data3[,1]
scaled_data3 = t(scale(t(mdata3[, 1:28])))
rownames(scaled_data3) = data3[,1]
rownames(scaled_data3) = data3$Gene



Heatmap(scaled_data3, name = "mat", border_gp = gpar(col = "black", lty = 2))

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("blue3", "white", "darkgreen"))
col_fun(seq(-3, 3))

col_fun = colorRamp2(c(-3, 0, 3), c("slategray4", "white", "tomato4"))
col_fun(seq(-3, 3))
#No clustering the columns
set.seed(3)

column_dend_height = unit(2, "cm"), 
row_dend_width = unit(2, "cm")

ht = Heatmap(scaled_data3, name = "z-score",col=col_fun, show_row_names=FALSE,
             show_column_names = FALSE, show_row_dend = FALSE,
             top_annotation = colAnn,
             left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = "white"), 
                                                              labels_rot =  0,
                                                              labels = c("Cluster 1\n n = 14/75 ", "Cluster 2\n n = 36/75", "Cluster 3\n n = 25/75"), 
                                                              labels_gp = gpar(col = "black", fontsize = 8))),
             row_km = 3, column_split = metadata3$Condition,
             cluster_columns = FALSE,
             column_dend_height = unit(0.5, "cm"), 
             row_dend_width = unit(0.5, "cm"))

#row_labels = scaled_data3$)

ht = draw(ht)

# col_fun(seq(-1, 1))

r.dend <- row_dend(ht)  #Extract row dendrogram
rcl.list <- row_order(ht)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters

# loop to extract genes for each cluster.
for (i in 1:length(row_order(ht))){
  if (i == 1) {
    clu <- t(t(row.names(scaled_data3[row_order(ht)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } 
  else {
    clu <- t(t(row.names(scaled_data3[row_order(ht)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}

out

#export
write.table(out, file= "N1_N2_Shaul_gene_clusters.txt", sep="\t", quote=F, row.names=FALSE)


