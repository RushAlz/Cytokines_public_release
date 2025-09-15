# tittle: Intercellular communication of the TGFb, COMPLEMENT, and SPP1 pathways in microglia 12 and 13
# author: Loren dos Santos

library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(kableExtra)
library(readxl) 
library(CellChat)
library(patchwork)
library(dplyr)
library(NMF) 
library(ggalluvial) 

# HDF5 file connection and indexing
hfile <- Connect(h5file)
hfile$index()
hfile[["assays"]]
hfile[["assays/SCT"]]
hfile[["reductions"]]
hfile$close_all()

seurat_obj <- LoadH5Seurat(h5file, assays = "SCT", reductions = "umap", graphs = NA) 
rownames_seurat_obj = rownames(seurat_obj@meta.data)

table(seurat_obj@meta.data$cogdx_3grp, useNA = "always")
seurat_obj@meta.data$cogdx_3grp[seurat_obj@meta.data$cogdx_3grp == 0] <- "NCI" # Non cognitive impairment - control
seurat_obj@meta.data$cogdx_3grp[seurat_obj@meta.data$cogdx_3grp == 1] <- "MCI" # Mild cognitive impairment 
seurat_obj@meta.data$cogdx_3grp[seurat_obj@meta.data$cogdx_3grp == 2] <- "AD" # 
seurat_obj@meta.data$cogdx_3grp <- as.factor(seurat_obj@meta.data$cogdx_3grp)
class(seurat_obj@meta.data$cogdx_3grp)

# All cells
Idents(seurat_obj) <- "cogdx_3grp" # set the identity to subset the data
table(Idents(seurat_obj)) 
meta_all <- seurat_obj@meta.data 
table(meta_all$cogdx_3grp)
Idents(seurat_obj) <- "state"
table(Idents(seurat_obj))

##############################
# Create a CellChat object
##############################
cellchat <- createCellChat(object = seurat_obj, 
                           group.by = "state",
                           assay = "SCT") 

# Preparing database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# Use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat) # Sub Define data based on genes expressed in cells

# identify overexpressed ligand-receptor pairs
# future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project gene expression data onto PPI network (optional)
cellchat <- smoothData(cellchat, adj = PPI.human)

# Significant - inference of cell-cell communication network

# Compute the communication probability and infer cellular communication network 
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1) 

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame 
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
all_interactions_cytokines <- subsetCommunication(cellchat, slot.name = "net")

# Visualization
# Visualize dominant senders (sources) and receivers (targets) in a 2D space
netAnalysis_signalingRole_scatter(cellchat)
netVisual_heatmap(cellchat, measure = "weight")

# Identify the major outgoing signaling events
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 6, slot.name = "netP", height = 15) 
ht1

# Identify the major incoming signaling events
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 6, height = 15) 
ht2

# Analysis of global communication patterns

# Infer the number of outgoing communication patterns

selectK(cellchat, pattern ="outgoing")
nPatterns = 4 
cellchat <- identifyCommunicationPatterns(cellchat, pattern ="outgoing", k = nPatterns, height = 17)

# Visualize the direct associations of cell groups and signaling pathways:
netAnalysis_dot(cellchat, slot.name = "netP", pattern ="outgoing", cutoff = NULL, color.use = NULL, dot.size =  c(1, 3), font.size = 7)

# Identify and visualize incoming communication patterns of target cells
selectK(cellchat, pattern = "incoming")

nPatterns = 3 
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, height = 17)

netAnalysis_river(cellchat, pattern = "incoming")

netAnalysis_dot(cellchat, pattern = "incoming", font.size = 7)

# Pathways Mic.12 and Mic.13

# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = c(5,6), targets.use = c(5,6), slot.name = "netP", legend.pos.x = 10)

# LR Mic.12 and Mic.13
netVisual_chord_gene(cellchat, sources.use = c(5,6), targets.use = c(5,6), lab.cex = 0.5,legend.pos.y = 30)

# LR - pathways: TGFb; SPP1; COMPLEMENT - Mic 12;13
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(5,6), targets.use = c(5,6), signaling = c("TGFb", "COMPLEMENT", "SPP1"),legend.pos.x = 8)

pathways.show <- c("TGFb") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# LR contribuition 
netAnalysis_contribution(cellchat, signaling = pathways.show) #TGFb

# Centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# LR - pathway
netVisual_chord_gene(cellchat, signaling = "TGFb", title.name = "TGFb signaling pathway network")


# Hierarchy plot
vertex.receiver <- c(5,6)  
netVisual_aggregate(cellchat, signaling = "TGFb", layout = "hierarchy", 
                    vertex.receiver = vertex.receiver)

# Circle plot
pathways.show <- c("COMPLEMENT") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")


# LR Contribuition
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# LR - pathway
netVisual_chord_gene(cellchat, signaling = "COMPLEMENT", title.name = "COMPLEMENT signaling pathway network")

# Hierarchy plot
# Definir os grupos de células que serão os receptores (alvos)
vertex.receiver <- c(5,6) 
netVisual_aggregate(cellchat, signaling = "COMPLEMENT", layout = "hierarchy", 
                    vertex.receiver = vertex.receiver)

# Circle plot
pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# LR Contribuition
netAnalysis_contribution(cellchat, signaling = pathways.show)

# Centrality scores
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# LR - pathway
netVisual_chord_gene(cellchat, signaling = "SPP1", title.name = "SPP1 signaling pathway network")

# Hierarchy plot
vertex.receiver <- c(5,6)  

netVisual_aggregate(cellchat, signaling = "SPP1", layout = "hierarchy", 
                    vertex.receiver = vertex.receiver)

# Significant LR interactions 
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = c(5,6), targets.use = c(5,6), signaling = c("TGFb", "COMPLEMENT", "SPP1"), remove.isolate = FALSE) 






