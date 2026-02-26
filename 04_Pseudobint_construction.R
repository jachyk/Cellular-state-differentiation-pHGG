# Discretization of Monocle3 pseudotime trajectories into pseudobins
# and transfer of trajectory-resolved states to Seurat objects
# for downstream comparative and spatial analyses.

#----- Load packages -----#

library(monocle3)
library(SingleCellExperiment)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(openxlsx)


#----- Load single-cell and Monocle3 objects -----#
# import neural subset (clusterresolution 0.7)
pHGG <- readRDS("neural_cells.rds")


cds <- load_monocle_objects(
  directory_path = "XXX/"
)


#----- Step 01: Transfer pseudotime values to Seurat metadata -----#

cds$all_pseudotime <- pseudotime(cds)
pseudotime_values <- cds$all_pseudotime

# check for coherent cell codes
cell_names_monocle <- colnames(cds)
cell_names_seurat  <- colnames(pHGG)

# create dataframe for subsequent insertion into 
# seurat object metadata 
pseudotime_df <- data.frame(
  Cell = cell_names_monocle,
  Pseudotime = pseudotime_values)
rownames(pseudotime_df) <- cell_names_monocle

# filter dataframe to include only common cells
common_cells <- intersect(
  cell_names_seurat,
  cell_names_monocle)
pseudotime_df_filtered <-
  pseudotime_df[common_cells, , drop = FALSE]

# Add pseudotime to Seurat object metadata 
pHGG <- AddMetaData(
  pHGG,
  metadata = pseudotime_df_filtered)

# Filter out infinite pseudotime values 
finite_pseudotime_indices <-
  which(is.finite(pHGG@meta.data$Pseudotime))
pHGG <- subset(pHGG,cells = finite_pseudotime_indices)

# control
head(pHGG@meta.data)


#----- Step 02: Bin pseudotime values for each trajectory -----# 

#----- Define trajectory labels in metadata -------------------#

Idents(pHGG) <- "SNNGRAPH_res.0.7"
pHGG@meta.data$Trajectories <- pHGG@meta.data[["SNNGRAPH_res.0.7"]]
unique(pHGG@meta.data$Trajectories)

# set trajectories by renaming duplicated cluster labels 
# in looped level-logic 
for (n in 1:length(levels(pHGG@meta.data$Trajectories))) {
  if (levels(pHGG@meta.data$Trajectories)[n] == "0")
    levels(pHGG@meta.data$Trajectories)[n] <- "OPC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "1")
    levels(pHGG@meta.data$Trajectories)[n] <- "Undiff"
  if (levels(pHGG@meta.data$Trajectories)[n] == "2")
    levels(pHGG@meta.data$Trajectories)[n] <- "Reactive"
  if (levels(pHGG@meta.data$Trajectories)[n] == "3")
    levels(pHGG@meta.data$Trajectories)[n] <- "AC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "4")
    levels(pHGG@meta.data$Trajectories)[n] <- "NPC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "5")
    levels(pHGG@meta.data$Trajectories)[n] <- "Reactive"
  if (levels(pHGG@meta.data$Trajectories)[n] == "6")
    levels(pHGG@meta.data$Trajectories)[n] <- "AC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "7")
    levels(pHGG@meta.data$Trajectories)[n] <- "Reactive"
  if (levels(pHGG@meta.data$Trajectories)[n] == "8")
    levels(pHGG@meta.data$Trajectories)[n] <- "Cycling"
  if (levels(pHGG@meta.data$Trajectories)[n] == "9")
    levels(pHGG@meta.data$Trajectories)[n] <- "AC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "10")
    levels(pHGG@meta.data$Trajectories)[n] <- "NPC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "11")
    levels(pHGG@meta.data$Trajectories)[n] <- "Cycling"
  if (levels(pHGG@meta.data$Trajectories)[n] == "12")
    levels(pHGG@meta.data$Trajectories)[n] <- "Cycling"
  if (levels(pHGG@meta.data$Trajectories)[n] == "13")
    levels(pHGG@meta.data$Trajectories)[n] <- "INP-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "14")
    levels(pHGG@meta.data$Trajectories)[n] <- "Reactive"
  if (levels(pHGG@meta.data$Trajectories)[n] == "15")
    levels(pHGG@meta.data$Trajectories)[n] <- "Int"
  if (levels(pHGG@meta.data$Trajectories)[n] == "16")
    levels(pHGG@meta.data$Trajectories)[n] <- "MOL"
  if (levels(pHGG@meta.data$Trajectories)[n] == "17")
    levels(pHGG@meta.data$Trajectories)[n] <- "OPC-like"
  if (levels(pHGG@meta.data$Trajectories)[n] == "18")
    levels(pHGG@meta.data$Trajectories)[n] <- "OPC-like"
}


# control plot 
DimPlot(
  pHGG,
  reduction = "umap.snn",
  group.by = "Trajectories",
  shuffle = TRUE)

# control if Trajectories are properly available
Idents(pHGG)<- "Trajectories"
unique(Idents(pHGG))



#----- AC_trajectory -----#### 
AC_pHGG <- subset(pHGG, idents=c("AC-like")) 

# control plot (input)
DimPlot(AC_pHGG, reduction = "umap.snn", group.by = "Trajectories",shuffle = T,ncol = 1, label = T, label.size = 6) 

# create Pseudotime_bin_clusters
bin_width <- 5

# fetch highest Pseudotime values
max_pseudotime <- max(AC_pHGG@meta.data$Pseudotime, na.rm = TRUE)

# ensure the upper limit is included in the bins as well
bin_breaks <- seq(0, ceiling(max_pseudotime/bin_width) * bin_width, by = bin_width)

AC_pHGG@meta.data$Pseudotime_bin <- cut(AC_pHGG@meta.data$Pseudotime, 
                                        breaks = bin_breaks,
                                        include.lowest = TRUE,
                                        labels = FALSE)
# include last bin (marked as NA so far) 
AC_pHGG@meta.data$Pseudotime_bin[is.na(AC_pHGG@meta.data$Pseudotime_bin) & AC_pHGG@meta.data$Pseudotime_bin == max_pseudotime] <- length(bin_breaks) - 1

# set the metadata  
AC_pHGG@meta.data$Pseudotime_bin <- as.character(AC_pHGG@meta.data$Pseudotime_bin)
pseudotime_levels <- sort(unique(AC_pHGG@meta.data$Pseudotime_bin))
AC_pHGG@meta.data$Pseudotime_bin <- factor(x = AC_pHGG@meta.data$Pseudotime_bin, levels = pseudotime_levels)

# control plot (output)
DimPlot(AC_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = T, label.size = 6) #,cols = color_sample)

# rename Pseudotime_bins 
unique(AC_pHGG@meta.data$Pseudotime_bin)
for (n in 1:length(levels(AC_pHGG@meta.data$Pseudotime_bin))) {
  if (levels(AC_pHGG@meta.data$Pseudotime_bin) [n] == "1")
    levels(AC_pHGG@meta.data$Pseudotime_bin) [n]<-"AC_Tr 1"
  if (levels(AC_pHGG@meta.data$Pseudotime_bin) [n] == "2")
    levels(AC_pHGG@meta.data$Pseudotime_bin) [n]<-"AC_Tr 2"
  if (levels(AC_pHGG@meta.data$Pseudotime_bin) [n] == "3")
    levels(AC_pHGG@meta.data$Pseudotime_bin )[n]<-"AC_Tr 3"
  if (levels(AC_pHGG@meta.data$Pseudotime_bin) [n] == "4")
    levels(AC_pHGG@meta.data$Pseudotime_bin )[n]<-"AC_Tr 4"
  if (levels(AC_pHGG@meta.data$Pseudotime_bin) [n] == "5")
    levels(AC_pHGG@meta.data$Pseudotime_bin) [n]<-"AC_Tr 5"
}

# control plot (output2)
color_bin <- c("AC_Tr 1"="#060459",
               "AC_Tr 2"= "#4C30A7",
               "AC_Tr 3"="#8E3BBC",
               "AC_Tr 4"="#D75437",
               "AC_Tr 5" = "#F2C72C")

DimPlot(AC_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = F, label.size = 6,cols = color_bin)


#------ OPC_Trajectory ------####
OPC_pHGG <- subset(pHGG, idents=c("OPC-like")) 

# control plot (input)
DimPlot(OPC_pHGG, reduction = "umap.snn", group.by = "Trajectories",shuffle = T,ncol = 1, label = T, label.size = 6)

#  Pseudotime_bin_clusters
bin_width <- 4

# fetch highest Pseudotime values
max_pseudotime <- max(OPC_pHGG@meta.data$Pseudotime, na.rm = TRUE)

# ensure the upper limit is included in the bins as well
bin_breaks <- seq(0, ceiling(max_pseudotime/bin_width) * bin_width, by = bin_width)

OPC_pHGG@meta.data$Pseudotime_bin <- cut(OPC_pHGG@meta.data$Pseudotime, 
                                         breaks = bin_breaks,
                                         include.lowest = TRUE,
                                         labels = FALSE)

# include last bin (marked as NA so far) 
OPC_pHGG@meta.data$Pseudotime_bin[is.na(OPC_pHGG@meta.data$Pseudotime_bin) & OPC_pHGG@meta.data$Pseudotime_bin == max_pseudotime] <- length(bin_breaks) - 1

# set the metadata 
OPC_pHGG@meta.data$Pseudotime_bin <- as.character(OPC_pHGG@meta.data$Pseudotime_bin)
pseudotime_levels <- sort(unique(OPC_pHGG@meta.data$Pseudotime_bin))
OPC_pHGG@meta.data$Pseudotime_bin <- factor(x = OPC_pHGG@meta.data$Pseudotime_bin, levels = pseudotime_levels)

# control plot (output)
DimPlot(OPC_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = T, label.size = 6) #,cols = color_sample)

#rename Pseudotime_bins 
unique(OPC_pHGG@meta.data$Pseudotime_bin)
for (n in 1:length(levels(OPC_pHGG@meta.data$Pseudotime_bin))) {
  if (levels(OPC_pHGG@meta.data$Pseudotime_bin) [n] == "1")
    levels(OPC_pHGG@meta.data$Pseudotime_bin) [n]<-"OPC_Tr 1"
  if (levels(OPC_pHGG@meta.data$Pseudotime_bin) [n] == "2")
    levels(OPC_pHGG@meta.data$Pseudotime_bin) [n]<-"OPC_Tr 2"
  if (levels(OPC_pHGG@meta.data$Pseudotime_bin) [n] == "3")
    levels(OPC_pHGG@meta.data$Pseudotime_bin )[n]<-"OPC_Tr 3"
  if (levels(OPC_pHGG@meta.data$Pseudotime_bin) [n] == "4")
    levels(OPC_pHGG@meta.data$Pseudotime_bin )[n]<-"OPC_Tr 4"
}

#control plot (output2)
color_bin <- c("OPC_Tr 1"= "#4C30A7",
               "OPC_Tr 2"="#8E3BBC",
               "OPC_Tr 3"="#D75437",
               "OPC_Tr 4" = "#F2C72C")

DimPlot(OPC_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = F, label.size = 6,cols = color_bin)



#----- NPC trajectory -----####
NPC_pHGG <- subset(pHGG, idents=c("NPC-like")) #for NPC_Trajectory

# control plot (input)
DimPlot(NPC_pHGG, reduction = "umap.snn", group.by = "Trajectories",shuffle = T,ncol = 1, label = T, label.size = 6) #,cols = color_sample)

# create Pseudotime_bin_clusters
bin_width <- 4

# fetch highest Pseudotime values
max_pseudotime <- max(NPC_pHGG@meta.data$Pseudotime, na.rm = TRUE)

# ensure the upper limit is included in the bins as well
bin_breaks <- seq(0, ceiling(max_pseudotime/bin_width) * bin_width, by = bin_width)

NPC_pHGG@meta.data$Pseudotime_bin <- cut(NPC_pHGG@meta.data$Pseudotime, 
                                         breaks = bin_breaks,
                                         include.lowest = TRUE,
                                         labels = FALSE)

#include last bin 
NPC_pHGG@meta.data$Pseudotime_bin[is.na(NPC_pHGG@meta.data$Pseudotime_bin) & NPC_pHGG@meta.data$Pseudotime_bin == max_pseudotime] <- length(bin_breaks) - 1

#set the metadata 
NPC_pHGG@meta.data$Pseudotime_bin <- as.character(NPC_pHGG@meta.data$Pseudotime_bin)
pseudotime_levels <- sort(unique(NPC_pHGG@meta.data$Pseudotime_bin))
NPC_pHGG@meta.data$Pseudotime_bin <- factor(x = NPC_pHGG@meta.data$Pseudotime_bin, levels = pseudotime_levels)

# control plot (output)
DimPlot(NPC_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = T, label.size = 6)

#rename Pseudotime_bins 
unique(NPC_pHGG@meta.data$Pseudotime_bin)
for (n in 1:length(levels(NPC_pHGG@meta.data$Pseudotime_bin))) {
  if (levels(NPC_pHGG@meta.data$Pseudotime_bin) [n] == "1")
    levels(NPC_pHGG@meta.data$Pseudotime_bin) [n]<-"NPC_Tr 1"
  if (levels(NPC_pHGG@meta.data$Pseudotime_bin) [n] == "2")
    levels(NPC_pHGG@meta.data$Pseudotime_bin) [n]<-"NPC_Tr 2"
  if (levels(NPC_pHGG@meta.data$Pseudotime_bin) [n] == "3")
    levels(NPC_pHGG@meta.data$Pseudotime_bin )[n]<-"NPC_Tr 3"
  if (levels(NPC_pHGG@meta.data$Pseudotime_bin) [n] == "4")
    levels(NPC_pHGG@meta.data$Pseudotime_bin )[n]<-"NPC_Tr 4"
}


#control plot (output2)

color_bin <- c("NPC_Tr 1"= "#4C30A7",
               "NPC_Tr 2"="#8E3BBC",
               "NPC_Tr 3"="#D75437",
               "NPC_Tr 4" = "#F2C72C")

DimPlot(NPC_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = F, label.size = 6,cols = color_bin)




#------ Reactive Trajectory ------####
Reactive_pHGG <- subset(pHGG, idents=c("Reactive")) #for Reactive_Trajectory

# control plot (input)
DimPlot(Reactive_pHGG, reduction = "umap.snn", group.by = "Trajectories",shuffle = T,ncol = 1, label = T, label.size = 6)

# create Pseudotime_bin_clusters
bin_width <- 5

# fetch highest Pseudotime values
max_pseudotime <- max(Reactive_pHGG@meta.data$Pseudotime, na.rm = TRUE)

# ensure the upper limit is included in the bins as well
bin_breaks <- seq(0, ceiling(max_pseudotime/bin_width) * bin_width, by = bin_width)

Reactive_pHGG@meta.data$Pseudotime_bin <- cut(Reactive_pHGG@meta.data$Pseudotime, 
                                              breaks = bin_breaks,
                                              include.lowest = TRUE,
                                              labels = FALSE)

#include last bin  
Reactive_pHGG@meta.data$Pseudotime_bin[is.na(Reactive_pHGG@meta.data$Pseudotime_bin) & Reactive_pHGG@meta.data$Pseudotime_bin == max_pseudotime] <- length(bin_breaks) - 1

#set the metadata 
Reactive_pHGG@meta.data$Pseudotime_bin <- as.character(Reactive_pHGG@meta.data$Pseudotime_bin)
pseudotime_levels <- sort(unique(Reactive_pHGG@meta.data$Pseudotime_bin))
Reactive_pHGG@meta.data$Pseudotime_bin <- factor(x = Reactive_pHGG@meta.data$Pseudotime_bin, levels = pseudotime_levels)

# control plot (output)
DimPlot(Reactive_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = T, label.size = 6) #,cols = color_sample)

#rename Pseudotime_bins 
unique(Reactive_pHGG@meta.data$Pseudotime_bin)
for (n in 1:length(levels(Reactive_pHGG@meta.data$Pseudotime_bin))) {
  if (levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n] == "1")
    levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n]<-"Reactive_Tr 1"
  if (levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n] == "2")
    levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n]<-"Reactive_Tr 2"
  if (levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n] == "3")
    levels(Reactive_pHGG@meta.data$Pseudotime_bin )[n]<-"Reactive_Tr 3"
  if (levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n] == "4")
    levels(Reactive_pHGG@meta.data$Pseudotime_bin )[n]<-"Reactive_Tr 4"
  if (levels(Reactive_pHGG@meta.data$Pseudotime_bin) [n] == "5")
    levels(Reactive_pHGG@meta.data$Pseudotime_bin )[n]<-"Reactive_Tr 5"
}

#control plot (output2)

color_bin <- c("Reactive_Tr 1"= "#4C30A7",
               "Reactive_Tr 2"="#8E3BBC",
               "Reactive_Tr 3"="#D75437",
               "Reactive_Tr 4" = "#F2C72C",
               "Reactive_Tr 5" = "#FDEC58")

DimPlot(Reactive_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = F, label.size = 6,cols = color_bin)


#----- Cycling cells -----####
Cycling_pHGG <- subset(pHGG, idents=c("Cycling")) #for Cycling_Trajectory

# control plot (input)
DimPlot(Cycling_pHGG, reduction = "umap.snn", group.by = "Trajectories",shuffle = T,ncol = 1, label = T, label.size = 6) #,cols = color_sample)

#create Pseudotime_bin_clusters
bin_width <- 8
#fetch highest Pseudotime values
max_pseudotime <- max(Cycling_pHGG@meta.data$Pseudotime, na.rm = TRUE)
#ensure the upper limit is included in the bins as well
bin_breaks <- seq(0, ceiling(max_pseudotime/bin_width) * bin_width, by = bin_width)


Cycling_pHGG@meta.data$Pseudotime_bin <- cut(Cycling_pHGG@meta.data$Pseudotime, 
                                             breaks = bin_breaks,
                                             include.lowest = TRUE,
                                             labels = FALSE)

#include last bin
Cycling_pHGG@meta.data$Pseudotime_bin[is.na(Cycling_pHGG@meta.data$Pseudotime_bin) & Cycling_pHGG@meta.data$Pseudotime_bin == max_pseudotime] <- length(bin_breaks) - 1
#set the metadata
Cycling_pHGG@meta.data$Pseudotime_bin <- as.character(Cycling_pHGG@meta.data$Pseudotime_bin)

pseudotime_levels <- sort(unique(Cycling_pHGG@meta.data$Pseudotime_bin))

Cycling_pHGG@meta.data$Pseudotime_bin <- factor(x = Cycling_pHGG@meta.data$Pseudotime_bin, levels = pseudotime_levels)

#control plot (output1) 
DimPlot(Cycling_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = T, label.size = 6) #,cols = color_sample)

#rename Pseudotime_bins 
unique(Cycling_pHGG@meta.data$Pseudotime_bin)
for (n in 1:length(levels(Cycling_pHGG@meta.data$Pseudotime_bin))) {
  if (levels(Cycling_pHGG@meta.data$Pseudotime_bin) [n] == "1")
    levels(Cycling_pHGG@meta.data$Pseudotime_bin) [n]<-"Cycling_Tr 1"
  if (levels(Cycling_pHGG@meta.data$Pseudotime_bin) [n] == "2")
    levels(Cycling_pHGG@meta.data$Pseudotime_bin) [n]<-"Cycling_Tr 2"
  if (levels(Cycling_pHGG@meta.data$Pseudotime_bin) [n] == "3")
    levels(Cycling_pHGG@meta.data$Pseudotime_bin )[n]<-"Cycling_Tr 3"
}

#control plot (output2) 

color_bin <- c("Cycling_Tr 1"= "#4C30A7",
               "Cycling_Tr 2"="#8F48F1",
               "Cycling_Tr 3" = "#F08A37")

DimPlot(Cycling_pHGG, reduction = "umap.snn", group.by = "Pseudotime_bin",shuffle = T,ncol = 1, label = F, label.size = 6,cols = color_bin)



#------ Step 03: Add new cluster information to metadata ------#####

unique(pHGG@meta.data$Trajectories)
Idents(pHGG)<- "Trajectories"
pHGG@meta.data$PB_clustering<-pHGG@meta.data[["Trajectories"]]

#--- add OPC trajectory
cells_in_cluster <- WhichCells(pHGG, idents = "OPC-like")

#check for same length befor importing
length(pHGG@meta.data$PB_clustering[pHGG@meta.data$Cell %in% cells_in_cluster])
length(OPC_pHGG@meta.data$Pseudotime_bin)

Indices<-which(pHGG@meta.data$Cell %in% OPC_pHGG@meta.data$Cell)
pHGG@meta.data$PB_clustering<- as.character(pHGG@meta.data$PB_clustering)
pHGG@meta.data$PB_clustering[Indices] <- as.character(OPC_pHGG@meta.data$Pseudotime_bin)


# control plot (OPC PT transfer)
DimPlot(pHGG, reduction = "umap.snn", group.by = "PB_clustering",shuffle = T,ncol = 1, label = F, label.size = 6)


#--- add NPC trajectory 
cells_in_cluster <- WhichCells(pHGG, idents = "NPC-like")

#check for same length befor importing
length(pHGG@meta.data$PB_clustering[pHGG@meta.data$Cell %in% cells_in_cluster])
length(NPC_pHGG@meta.data$Pseudotime_bin)
Indices<-which(pHGG@meta.data$Cell %in% NPC_pHGG@meta.data$Cell)
pHGG@meta.data$PB_clustering<- as.character(pHGG@meta.data$PB_clustering)
pHGG@meta.data$PB_clustering[Indices] <- as.character(NPC_pHGG@meta.data$Pseudotime_bin)

# control plot (NPC PT transfer)
DimPlot(pHGG, reduction = "umap.snn", group.by = "PB_clustering",shuffle = T,ncol = 1, label = F, label.size = 6)


#--- add AC trajectory
cells_in_cluster <- WhichCells(pHGG, idents = "AC-like")

#check for same length befor importing
length(pHGG@meta.data$PB_clustering[pHGG@meta.data$Cell %in% cells_in_cluster])
length(AC_pHGG@meta.data$Pseudotime_bin)
Indices<-which(pHGG@meta.data$Cell %in% AC_pHGG@meta.data$Cell)
pHGG@meta.data$PB_clustering<- as.character(pHGG@meta.data$PB_clustering)
pHGG@meta.data$PB_clustering[Indices] <- as.character(AC_pHGG@meta.data$Pseudotime_bin)


#control plot (AC PT transfer)
DimPlot(pHGG, reduction = "umap.snn", group.by = "PB_clustering",shuffle = T,ncol = 1, label = F, label.size = 6)


#--- add Reactive trajectory 
cells_in_cluster <- WhichCells(pHGG, idents = "Reactive")

#check for same length befor importing
length(pHGG@meta.data$PB_clustering[pHGG@meta.data$Cell %in% cells_in_cluster])
length(Reactive_pHGG@meta.data$Pseudotime_bin)
Indices<-which(pHGG@meta.data$Cell %in% Reactive_pHGG@meta.data$Cell)
pHGG@meta.data$PB_clustering<- as.character(pHGG@meta.data$PB_clustering)
pHGG@meta.data$PB_clustering[Indices] <- as.character(Reactive_pHGG@meta.data$Pseudotime_bin)

# control plot (reactive PT transfer)
DimPlot(pHGG, reduction = "umap.snn", group.by = "PB_clustering",shuffle = T,ncol = 1, label = F, label.size = 6)

#--- add Cycling trajectory
cells_in_cluster <- WhichCells(pHGG, idents = "Cycling")
#check for same length befor importing
length(pHGG@meta.data$PB_clustering[pHGG@meta.data$Cell %in% cells_in_cluster])
length(Cycling_pHGG@meta.data$Pseudotime_bin)
Indices<-which(pHGG@meta.data$Cell %in% Cycling_pHGG@meta.data$Cell)
pHGG@meta.data$PB_clustering<- as.character(pHGG@meta.data$PB_clustering)
pHGG@meta.data$PB_clustering[Indices] <- as.character(Cycling_pHGG@meta.data$Pseudotime_bin)

# control plot (Cycling PT transfer)
DimPlot(pHGG, reduction = "umap.snn", group.by = "PB_clustering",shuffle = T,ncol = 1, label = F, label.size = 6)



# set colors for plotting 
CC_binned <- c("AC_Tr 1"="#51233D",
               "AC_Tr 2"= "#7250C2",
               "AC_Tr 3"="#984343",
               "AC_Tr 4"="#E07A61",
               "AC_Tr 5" = "#F1BB32",
               "OPC_Tr 1"= "#1A5851",
               "OPC_Tr 2"="#2A867C",
               "OPC_Tr 3"="#3F7A5E",
               "OPC_Tr 4" = "#3DBDB2",
               "NPC_Tr 1"= "#144775",
               "NPC_Tr 2"="#3746B7",
               "NPC_Tr 3"="#6471D2",
               "NPC_Tr 4" = "#70B6ED",
               "Reactive_Tr 1"= "#4E3570",
               "Reactive_Tr 2"="#862828",
               "Reactive_Tr 3"="#B83F6C",
               "Reactive_Tr 4" = "#F69F49",
               "Reactive_Tr 5" = "#D49F85",
               "Cycling_Tr 1"= "#3B2869",
               "Cycling_Tr 2"="#A926AC",
               "Cycling_Tr 3" = "#EE8EF0",
               "HOG" = "#765DB5",
               "INP-like"= "#625CA0",
               "MOL"="#22CA97",
               "Undiff"="#252A2E"
)

DimPlot(pHGG, reduction = "umap.snn", group.by = "PB_clustering",shuffle = T,ncol = 1, label = F, label.size = 6,pt.size = 1.2, cols = CC_binned)


#----- Step 04: Merge neural cells with immune and stromal cells ----#
#import integrated object (clusterresolution 0.5) for merging

integrated <- readRDS("/integrated_pHGG_res.05.rds")

Idents(integrated)<-"integrated_snn_res.0.5" #set to right clustering the neural object is derived from
integrated@meta.data[["integrated_snn_res.0.5"]]
integrated@meta.data$Pseudobin<-integrated@meta.data[["integrated_snn_res.0.5"]]
unique(integrated@meta.data$Pseudobin)
levels(integrated@meta.data$Pseudobin)


# Pseuditime-binned (PB) cells of the integrated object display neural cells 
# which have been subsetted for analysis
# The goal is to overwrite these exact cells 
for (n in 1:length(levels(integrated@meta.data$Pseudobin))) {
  if (levels(integrated@meta.data$Pseudobin) [n] == "0")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "1")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "2")
    levels(integrated@meta.data$Pseudobin )[n]<-"Myeloid"
  if (levels(integrated@meta.data$Pseudobin) [n] == "3")
    levels(integrated@meta.data$Pseudobin )[n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "4")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "5")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "6")
    levels(integrated@meta.data$Pseudobin) [n]<-"MT"
  if (levels(integrated@meta.data$Pseudobin) [n] == "7")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "8")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "9")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "10")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin )[n] == "11")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "12")
    levels(integrated@meta.data$Pseudobin) [n]<-"Myeloid"
  if (levels(integrated@meta.data$Pseudobin) [n] == "13")
    levels(integrated@meta.data$Pseudobin) [n]<-"Perivascular"
  if (levels(integrated@meta.data$Pseudobin) [n] == "14")
    levels(integrated@meta.data$Pseudobin) [n]<-"Endothelial"
  if (levels(integrated@meta.data$Pseudobin) [n] == "15")
    levels(integrated@meta.data$Pseudobin) [n]<-"MOL"
  if (levels(integrated@meta.data$Pseudobin) [n] == "16")
    levels(integrated@meta.data$Pseudobin) [n]<-"Lymphoid"
  if (levels(integrated@meta.data$Pseudobin) [n] == "17")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "18")
    levels(integrated@meta.data$Pseudobin) [n]<-"PB"
  if (levels(integrated@meta.data$Pseudobin) [n] == "19")
    levels(integrated@meta.data$Pseudobin) [n]<-"Cycling Myeloid"
}

# control plot 
DimPlot(integrated, reduction = "umap.snn", group.by = "Pseudobin",shuffle = T,ncol = 1, label = F, label.size = 6, pt.size = 1.2,cols = CC_binned) 
integrated@meta.data$Cell <- rownames(integrated@meta.data)

# extract PB cells 
Idents(integrated)<-"Pseudobin" 
cells_in_cluster <- WhichCells(integrated, idents = "PB")

# check for same length before importing cell identities 
length(integrated@meta.data$Pseudobin[integrated@meta.data$Cell %in% cells_in_cluster])
length(pHGG@meta.data$PB_clustering)

# match cells between the objects for integration
Indices<-which(integrated@meta.data$Cell %in% pHGG@meta.data$Cell)
integrated@meta.data$Pseudobin<- as.character(integrated@meta.data$Pseudobin)
integrated@meta.data$Pseudobin[Indices] <- as.character(pHGG@meta.data$PB_clustering)

#check if integration was successfull 
DimPlot(integrated, reduction = "umap.snn", group.by = "Pseudobin",shuffle = T,ncol = 1, label = F, label.size = 6, pt.size = 1.2,cols = CC_binned) 

saveRDS(object = integrated,file = "integrated_Pseudobin.rds")

