if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Version 3.16 of BiocManager is not used with the latest version of R. Latest version of R needs BiocManager version 3.17
# But here I have not installed the latest version as I don't want to install all the packages again
BiocManager::install(version = "3.16")

# Firstly, I installl MungeSumStats as it is a dependency for echodata
BiocManager::install("MungeSumstats")

#  I also need rlang version 1.1.0 to load MungeSumstats ( I had version 1.0.6)
install.packages("rlang", version="1.1.0")
library(rlang)
library(MungeSumstats)
# Then, I install echodata as it is a dependency for echoconda
remotes::install_github("RajLabMSSM/echodata")
library(echodata)

# After, I install echoconda because it's a dependency for scKirby
remotes::install_github("RajLabMSSM/echoconda")
library(echoconda)

# scKirby also need the installation od sceasy
devtools::install_github("cellgeni/sceasy")
library(sceasy)
# Now hopefully I will install scKirby
remotes::install_github("neurogenomics/scKirby")
library(scKirby)

install.packages('reticulate')
library(reticulate)

#install loomR AND hdf5r packages that we may need for loading the loom file
devtools::install_github(repo = "hhoeflin/hdf5r")
library(hdf5r)
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)

getwd()
setwd("/Users/nikolasgiannakis/Desktop/")
getwd()
# Connect to loom file
loom_file <- connect(filename = "dev_all.loom", mode = "r+", skip.validate = TRUE)

seurat_object <- CreateSeuratObject(counts = as.matrix(loom_file$exprs))



library(Seurat)
library(loomR)
library(SeuratWrappers)
library(reticulate)
library(SeuratDisk)
loom_file <- "/Users/nikolasgiannakis/Desktop/dev_all.loom"
seurat_looom = LoadLoom(
  loom_file,
  assay = NULL,
  cells = "CellID",
  features = "Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE
)




loom_data= read.loom(loom_file)


# List all Conda environments
conda_envs <- conda_list()

# Display the Conda environments
print(conda_envs)


library(devtools)

#use_python("/Users/nikolasgiannakis/Library/r-miniconda/envs/r-reticulate/bin/python")

use_condaenv("/Users/nikolasgiannakis/Library/r-miniconda/envs/r-reticulate")
reticulate::py_install("loompy")

# Import the necessary Python modules
# Uninstall the current version of NumPy
# Uninstall the current version of NumPy (optional)
py_run_string("import pip; pip.main(['uninstall', '-y', 'numpy'])")

# Install the desired version of NumPy (e.g., 1.23.3)
py_run_string("import pip; pip.main(['install', 'numpy==1.23.3'])")

#py_install("numpy==1.24", pip = TRUE)
#py_install("numpy<1.24", pip = TRUE)

numpy <- import("numpy")
numpy$version$version
# Get the path of the NumPy module
numpy_path <- py_module_path(numpy)

# Print the path
print(numpy_path)
h5py <- import("h5py")

#loompy = reticulate::py_install("loompy")
loompy <- import("loompy")

# Read the loom file using loompy in Python
# Read the loom file using loompy in Python
loom_file <- "/Users/nikolasgiannakis/Desktop/dev_all.loom"
dataset <- connect(filename = "dev_all.loom", mode = "r+", skip.validate = TRUE)

# Access the count matrix
#counts <- dataset$[:, :]$astype(numpy$float32)
#counts <- dataset[[""]]$astype(numpy$float32)
#counts <- dataset$matrix$astype(numpy$float32)
counts <- dataset$matrix$astype(numpy$float32)

counts <- dataset$[:, :]$astype(numpy$float32)
# Convert loom data to Seurat object
seurat_object <- CreateSeuratObject(counts = counts)

# Optionally, you can assign other information from the loom file to the Seurat object
seurat_object$cellnames <- dataset$col_attrs$Name
seurat_object$genenames <- dataset$row_attrs$Name


# Perform additional preprocessing and analysis on the Seurat object as desired
# ...

# Access the data in the Seurat object
# For example, to access the expression matrix:
expression_matrix <- GetAssayData(seurat_object)

# To access the metadata or cell attributes:
metadata <- seurat_object$meta.data
























#  loom <- loomfile
library(SeuratObject)
library(tidyr)
obj = example_obj("loom")
obj2 <- loom_to_seurat(obj) # it gives error with the Csparse validateIndices

Assays(obj2)
Graphs(obj2)
Neighbors(obj2)
Reductions(obj2)

Stdev(obj2)
head(obj2)

obj2@graphs$RNA_snn
obj2@images
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
obj2[["percent.mt"]] <- PercentageFeatureSet(obj2, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(obj2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(obj2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# How many cells are in each cluster
table(Idents(obj2))
# How many cells are in each replicate?
table(obj2$replicate)

DefaultAssay(obj2) <- 'RNA'
obj2 <- NormalizeData(obj2) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')


head(Loadings(obj2, reduction = "pca")[, 1:5])
DimPlot(obj2, label = TRUE)

DimHeatmap(obj2$RNA_snn_res.0.8, dims = 1, cells = 500, balanced = TRUE)




#to solve the issue I do the following
install.packages('installr')
library(installr)
install.Rtools(keep_install_file=TRUE)

remove.packages("Matrix")
install.packages("Matrix", version ="1.5.3")
library(Matrix)


# Now I repeat the previous step. Here to say that loom_to_seurat is a function of the scKirby package but because i couldn't
# use the suggested function of the package that convertt s from loom to seurat ( ingest_data ), as it was giving me and error
# with the filetype_dict., I used the specific function for conversino from loom to seurat and it worked.



obj2 <- loom_to_seurat(obj)
library(Seurat)

# perform visualization and clustering steps
obj2 <- NormalizeData(obj2)
obj2 <- FindVariableFeatures(obj2)
obj2 <- ScaleData(obj2)
obj2 <- RunPCA(obj2, verbose = FALSE, approx=FALSE)
obj2 <- FindNeighbors(obj2, dims = 1:30)
obj2 <- FindClusters(obj2, resolution = 3, verbose = FALSE)
obj2 <- RunUMAP(obj2, dims = 1:30)
DimPlot(obj2, label = TRUE)

DimHeatmap(obj2, dims = 1, cells = 500, balanced = TRUE)









saveRDS(obj2, file = "devDataLinarssonLab.rds")
getwd()
setwd("/Users/nikolasgiannakis/Desktop/LinarssonLab/")
seuratObdevDataLinarssonLab = readRDS(file = "devDataLinarssonLab.rds")









# now i want to save this seurat file and i will use t he seuratdisk package
library(Seurat)
library(SeuratDisk)
devtools::install_github('satijalab/seurat-data')
library(SeuratData)
SaveH5Seurat(obj2, overwrite = TRUE)

#    print(loom)
#                print(loom$filename)
#Convert loom object to seurat
seurat <- ingest_data(obj=obj, output_type = "Seurat") 

loom_to_seurat <- function(obj,
                           verbose=TRUE){
  messager("+ loom ==> Seurat",v=verbose)
  obj2 <- Seurat::as.Seurat(obj)
  return(obj2)
}

dev_all.loom <- connect(filename = "dev_all.loom", mode = "r+", skip.validate = TRUE)



sceasy::convertFormat('dev_all.loom', from="loom", to="sce",
                      outFile='dev_all.rds')



Rpackages_20230718 = tibble::tibble(
  Package = names(installed.packages()[,3]),
  Version = unname(installed.packages()[,3])
)

write.table(Rpackages_20230718, "Rpackages_20230718.csv")
