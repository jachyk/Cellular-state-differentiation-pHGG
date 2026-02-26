# Note that this workflow is based on the Seurat vignette
# "Analysis, visualization, and integration of spatial datasets with Seurat"
# (vignettes/spatial_vignette.Rmd, Seurat team) and was adapted for this study.
# Citation: Hao et al. Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature Biotechnology (2023) (Seurat V5)


#----- Load packages -----#

library(Seurat)      
library(ggplot2)
library(hdf5r)
library(dplyr)
library(matrixStats)


#----- Define directories -----#

Exp_dir <- "results/spatial_mapping"
dir.create(Exp_dir, recursive = TRUE, showWarnings = FALSE)

plot_dir_pca <- file.path(Exp_dir, "pcaproject")
plot_dir_cca <- file.path(Exp_dir, "cca")
dir.create(plot_dir_pca, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir_cca, recursive = TRUE, showWarnings = FALSE)

ST_source <-"samplewise_ST_data/"

#----- Load reference single-cell dataset -----#
# this is created in >> Script 04 << 
SC_directory <- "/Pseudotime_binned_dataset.rds"
integrated <- readRDS(SC_directory)

#control plot 
DimPlot(
  integrated,
  reduction = "umap.snn",
  group.by = "Pseudobin",
  shuffle = TRUE,
  pt.size = 1.2)


#----- Load individual spatial dataset -----#
# this workflow uses anchor transfer for one 
# sample at a time
# for higher throughput see loop below

ST_directory <- "XX.RData"
load(ST_directory)   
#control plot 
SpatialFeaturePlot(st, features = "nCount_Spatial", pt.size.factor = 2)

nrow(st@meta.data)
length(rownames(st))


#----- Normalize spatial and reference datasets -----#

st <- SCTransform(
  st,
  assay = "Spatial",
  verbose = FALSE,
  ncells = nrow(st@meta.data),
  n_genes = length(rownames(st)),
  new.assay.name = "SCT"
) %>%
  RunPCA(verbose = FALSE)
#normalization step was done for every unique sample 
#and thus left out in the loop below

#save changes in a new object 
DST <- st
DefaultAssay(DST) <- "SCT"


integrated <- SCTransform(
  integrated,
  ncells = 65246,
  verbose = TRUE,
  new.assay.name = "SCT",
  n_genes = length(rownames(integrated))
) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction.key = "SCT_UMAP")

DefaultAssay(integrated) <- "SCT"


#----- Transfer anchors and prediction (for one sample) -----#

anchors <- FindTransferAnchors(
  reference = integrated,
  query = DST,
  normalization.method = "SCT",
  reduction = "pcaproject",
  k.anchor = 60,
  k.score = 50,
  query.assay = "SCT",
  reference.assay = "SCT",
  n.trees = 100,
  dims = 1:50)

predictions.assay <- TransferData(
  anchorset = anchors,
  refdata = integrated$Pseudobin,
  prediction.assay = TRUE,
  weight.reduction = DST[["pca"]],
  dims = 1:30)

DST[["predictions"]] <- predictions.assay


#------Visualization -----------------------------#
#----- Cluster definitions for visualization -----#

Clusters <- c(
  "AC-Tr 1","AC-Tr 2","AC-Tr 3","AC-Tr 4","AC-Tr 5",
  "OPC-Tr 1","OPC-Tr 2","OPC-Tr 3","OPC-Tr 4",
  "NPC-Tr 1","NPC-Tr 2","NPC-Tr 3","NPC-Tr 4",
  "Reactive-Tr 1","Reactive-Tr 2","Reactive-Tr 3","Reactive-Tr 4","Reactive-Tr 5",
  "Cycling-Tr 1","Cycling-Tr 2","Cycling-Tr 3",
  "HOG","INP-like","MOL","Undiff","Perivascular",
  "Lymphoid","MT","Myeloid","Endothelial","Cycling Myeloid"
)

DefaultAssay(DST) <- "predictions"

sample <- "WT_1"

SpatialFeaturePlot(
  DST,
  features = "AC-Tr 1",
  pt.size.factor = 1.9,
  crop = TRUE
)


#----- Plot prediction maps (for one sample) -----#

for (i in Clusters) {
  try({
    message(i)
    png(
      filename = file.path(
        Exp_dir,
        paste0(i, "_", sample, ".png")
      ),
      width = 600,
      height = 500)
    print(
      SpatialFeaturePlot(
        DST,
        features = i,
        pt.size.factor = 1.8,
        crop = TRUE))
    dev.off()
  }, silent = TRUE)
}



#----looped version for faster workflow ---#
#----- Sample definitions -----------------#
# for looped version make sure your data directory is properly structured
# and sample data is named accordingly 

Infant <- c("IHG_1","IHG_2 A","IHG_2 B","IHG_3 A","IHG_3 B","IHG_4","IHG_5")
WT     <- c("WT_1","WT_2","WT_3 A","WT_3 B","WT_4","WT_5","WT_6","WT_7")
DMG    <- c("DMG_1","DMG_2","DMG_3","DMG_4")
DHG    <- c("DHG_1","DHG_2","DHG_3","DHG_4")

samples <- c(WT, Infant, DMG, DHG)

affiliation <- c(
  "WT","WT","WT","WT","WT",
  "Infant","Infant","Infant","Infant","Infant",
  "DMG","DMG","DMG","DMG",
  "DHG","DHG","DHG","DHG"
)


#----- Spatial mapping loop: PCA projection -----#

j <- 0
for (s in samples) {
  j <- j + 1
  a <- affiliation[j]
  load(paste0(ST_source, s, "_ST.RData"))               
  DST <- st
  DefaultAssay(DST) <- "SCT"
  anchors <- FindTransferAnchors(
    reference = integrated,
    query = DST,
    normalization.method = "SCT",
    reduction = "pcaproject",
    k.anchor = 60,
    k.score = 50,
    query.assay = "SCT",
    reference.assay = "SCT",
    n.trees = 100,
    dims = 1:50
  )
  
  predictions.assay <- TransferData(
    anchorset = anchors,
    refdata = integrated$Pseudobin,
    prediction.assay = TRUE,
    weight.reduction = DST[["pca"]],
    dims = 1:30
  )
  
  DST[["predictions"]] <- predictions.assay
  DefaultAssay(DST) <- "predictions"
  
  sample <- paste0(s, a)
  
  for (i in Clusters) {
    
    try({
      
      png(
        filename = file.path(
          plot_dir_pca,
          paste0(i, "_PCAproject_", sample, ".png")
        ),
        width = 600,
        height = 500
      )
      
      print(
        SpatialFeaturePlot(
          DST,
          features = i,
          pt.size.factor = 1.8,
          crop = TRUE
        )
      )
      
      dev.off()
      
    }, silent = TRUE)
  }
}


#----- Spatial mapping loop: CCA projection -----#

H <- 0
for (s in samples) {
  H <- H + 1
  a <- affiliation[H]
  load(paste0(ST_source, s, "_ST.RData"))
  DST <- st
  DefaultAssay(DST) <- "SCT"
  
  anchors <- FindTransferAnchors(
    reference = integrated,
    query = DST,
    normalization.method = "SCT",
    reduction = "cca",
    k.anchor = 60,
    k.score = 50,
    query.assay = "SCT",
    reference.assay = "SCT",
    n.trees = 100,
    dims = 1:50
  )
  
  predictions.assay <- TransferData(
    anchorset = anchors,
    refdata = integrated$Pseudobin,
    prediction.assay = TRUE,
    weight.reduction = DST[["pca"]],
    dims = 1:30
  )
  
  DST[["predictions"]] <- predictions.assay
  DefaultAssay(DST) <- "predictions"
  
  sample <- paste0(s, a)
  
  for (i in Clusters) {
    
    try({
      
      png(
        filename = file.path(
          plot_dir_cca,
          paste0(i, "_CCA_", sample, ".png")
        ),
        width = 600,
        height = 500
      )
      
      print(
        SpatialFeaturePlot(
          DST,
          features = i,
          pt.size.factor = 1.8,
          crop = TRUE
        )
      )
      
      dev.off()
      
    }, silent = TRUE)
  }
}
