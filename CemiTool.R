
library(Seurat)
library(SeuratData)
library(parallel)
library(MASS)
library(patchwork)
library(SeuratWrappers)
#I havent been able to install monocle3.     library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(slingshot)
library(scales)
library(viridis)
library(viridisLite)
library(RColorBrewer)
BiocManager::install("dittoSeq")
library(dittoSeq)
install.packages("ggplot.multistats")
library(ggplot.multistats)
BiocManager::install('multtest')
library(multtest)
install.packages('metap')
library(metap)

load(file = "/Users/nikolasgiannakis/Desktop/20211108_ImmuneSOX2negCells_named.RData")


#Load packages/libraries for cell-cell interaction analysis
install.packages('NMF')
library(NMF)
install.packages("circlize")
library(circlize)
library(remotes)
library(ComplexHeatmap)

options(download.file.method = "wget")
devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

#Extract the input files for the CellChat from the Seurat object
data.input = GetAssayData(Immune_SOX2negCells2_named, assay = "RNA", slot = "data")
labels = Idents(Immune_SOX2negCells_named)
names(labels)
meta = data.frame(group = labels, row.names = names(labels))

#Create a CellChat object using data matrix as input
cellchat = createCellChat (object = data.input, meta = meta, group.by = "group")


#cellchat <- addMeta(cellchat, meta = meta)
#ellchat <- setIdent(cellchat, ident.use = "labels")

levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 

#CellChatDB <- CellChatDB.human 
CellChatDB = CellChatDB.mouse# use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

write.csv(CellChatDB$interaction, file = "CellChatDB_mouse_Interaction_v2.csv")
cellchat@DB <- CellChatDB
#Preprocessing the expression data for cell-cell communication analysis
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)

cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

#Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#Here at first instance I had 10, now I changed it to 20, This option is not at the tutorial for one Seurat object
cellchat <- filterCommunication(cellchat, min.cells = 20)
#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#visualize the aggregated cell-cell communication network.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
cellchat@net$weight
#Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group.
mat <- cellchat@net$weight
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("CXCL")
pathways.show <- c("TGFb") 
pathways.show <- c("VEGF") 
pathways.show <- c("CCL") 
pathways.show <- c("IFN-I") 
pathways.show <- c("ANGPTL") 
pathways.show <- c("COLLAGEN") 
pathways.show <- c("MHC-I") 
pathways.show = c("NOTCH")
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

#Compute and visualize the network centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "net") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
#The original command for centrality's calculation is the one down, not above this line
#cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 10, font.size = 10)


#Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest

#Maybe I will need to run this line to make the heatmap script to work.
df$labels = factor(df$labels,levels = levels(object@idents))

gg2 <- netAnalysis_signalingRole_(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")

pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat, signaling = pathways.show)

pair.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
CXCLpair.show <- pair.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = CXCLpair.show, vertex.receiver = vertex.receiver)
