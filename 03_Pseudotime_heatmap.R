#----- Load packages -----#

library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)


#----- Load data and add pseudotime metadata -----#
# import neural subcluster
pHGG <- readRDS("/neural_cells.rds")

# import pseudotime values per cell, generated in >> Script 01 <<
DF_pseudotime <- read.xlsx("/DF_pseudo_all_cells.xlsx")

pseudotime_df <- data.frame(
  ID = DF_pseudotime$cells,
  Pseudotime = DF_pseudotime$all3_pseudotime)

pseudotime_df <- pseudotime_df[match(Cells(pHGG), pseudotime_df$ID),]

pHGG <- AddMetaData(
  object = pHGG,
  metadata = pseudotime_df$Pseudotime,
  col.name = "Pseudotime")


#----- Marker gene options representing 4 lineages-----#

OPC <- c("PTPRZ1","C1orf61","MOXD1","KLHDC8A","IGFBP2","SPRY4",
         "EGR1","CSPG4","COL20A1","OLIG2","OLIG1","TNR","FYN",
         "GPR17","PLP1","MBP")

AC <- c("CLU","GFAP","MOXD1","NTRK2","GLI3","SOX9","ID4",
        "TTYH1","F3","SPARCL1","AQP4","SLC1A2","GJA1","APLNR")

OPC <- c("PDGFRA","CSPG4","OLIG2","OLIG1","NTRK3","SOX10",
         "TNR","SLC44A1","NTRK2","PLP1","GPR37L1",
         "GPR17","MBP","CDK18","MOG")

NPC <- c("SOX2","SOX4","DCX","NNAT","PTPRN2","STXBP1","NDRG4",
         "TMEM59L","UCHL1","DLX1","DLX2","DLX5","SCG2",
         "OLFM1","SNAP25","NSF","SYP","SYT1")

marker <- na.omit(marker)


#----- Lineage options -----#

lineage1 <- c("1","4","10","13")   # NPC
lineage2 <- c("1","0","17","18")   # OPC
lineage3 <- c("1","6","3","9","5") # AC
lineage4 <- c("1","2","5","7")     # Reactive


#----- Subset according to lineage and subtype -----#

pHGG_sub <- subset(pHGG,subset = SNNGRAPH_res.0.7 %in% lineage3)

#subtype options (strings) are stored in the metadata
pHGG_sub <- subset(pHGG_sub,subset = subtype %in% "WT")

# check if correct cells have been subsetted
DimPlot(
  pHGG_sub,
  reduction = "umap.snn",
  group.by = "SNNGRAPH_res.0.7",
  shuffle = TRUE)


#----- Prepare expression matrix ordered by pseudotime for plotting-----#

expression_data <- FetchData(pHGG_sub, vars = marker)

pseudotime <- pHGG_sub$Pseudotime
expression_data$Pseudotime <- pseudotime
# check and compare data length 
length(pseudotime)
length(expression_data[,1])
# order cells by pseudotime
order_cells <- order(pseudotime)
expression_data <- expression_data[order_cells, ]
pseudotime <- pseudotime[order_cells]
# check and compare data order
head(pseudotime)
expression_data[1:5,1:5]

# transpose the expression data to match the desired format for plotting
expression_data <- t(expression_data)
Pseudotime_lin <- expression_data[nrow(expression_data), ]
expression_data_filtered <- expression_data[-nrow(expression_data), ]

#----- Heatmap: plot marker dynamics along pseudotime -----#
# see also below: looped version for faster workflow #

heatmap_dir <- file.path(output_dir, "marker_heatmaps")
ptscale_dir <- file.path(output_dir, "pseudotime_scales")

dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ptscale_dir, recursive = TRUE, showWarnings = FALSE)


col_fun2 <- colorRamp2(c(0,2,6), hcl_palette = "Viridis")

#name plot after respective lineage
tiff(
  filename = file.path(
    heatmap_dir,
    "AC1_Pseudoheatmap_all.tiff"
  ),
  width = 600,
  height = 800
)

Heatmap(
  expression_data_filtered,
  name = "Expression",
  column_title = "Pseudotime WT",
  row_title = "Marker Genes",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_order = order(pseudotime),
  row_names_gp = gpar(fontsize = 17),
  col = col_fun2,
  use_raster = TRUE
)

dev.off()


#----- Heatmap: pseudotime scale for comparability  -----#

col_fun <- colorRamp2(c(0,7,23.33), hcl_palette = "Inferno")

tiff(
  filename = file.path(
    ptscale_dir,
    "Pseudotime_OPC_all.tiff"
  ),
  width = 600,
  height = 200
)

Heatmap(
  t(Pseudotime_lin),
  name = "Pseudotime",
  column_title = "Pseudotime Infant",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  col = col_fun
)

dev.off()

#----- looped version for faster workflow ---- #
#----- Subtype-wise visualization ------------#
# make sure to change the plotting directories
# before to avoid accidental overwriting

# set lineage for which the loop should be started 
pHGG_sub <- subset(
  pHGG,
  subset = SNNGRAPH_res.0.7 %in% lineage3
)

DimPlot(
  pHGG_sub,
  reduction = "umap.snn",
  group.by = "SNNGRAPH_res.0.7",
  shuffle = TRUE
)

subtypes <- c("WT","Infant","DHG","DMG")

for (i in subtypes) {
  
  message(i)
  
  pHGG_sub1 <- subset(pHGG_sub, subset = subtype %in% i)
  expression_data <- FetchData(pHGG_sub1, vars = marker)
  pseudotime <- pHGG_sub1$Pseudotime
  expression_data$Pseudotime <- pseudotime
  order_cells <- order(pseudotime)
  expression_data <- expression_data[order_cells, ]
  pseudotime <- pseudotime[order_cells]
  expression_data <- t(expression_data)
  Pseudotime_lin <- expression_data[nrow(expression_data), ]
  expression_data_filtered <- expression_data[-nrow(expression_data), ]
  
  col_fun2 <- colorRamp2(c(0,2,6), hcl_palette = "Viridis")
  
  tiff(
    filename = file.path(
      heatmap_dir,
      paste0("Pseudoheatmap_lineage3_", i, ".tiff")
    ),
    width = 500,
    height = 600
  )
  
  print(
    Heatmap(
      expression_data_filtered,
      name = "Expression",
      column_title = paste0("Pseudotime ", i),
      row_title = "Marker Genes",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = TRUE,
      column_order = order(pseudotime),
      row_names_gp = gpar(fontsize = 16),
      col = col_fun2
    )
  )
  
  dev.off()
  
  
  # pseudotime scale for comparability 
  
  col_fun <- colorRamp2(c(0,7,23.33), hcl_palette = "Inferno")
  marks <- c(0,5,10,15,20,25)
  PT_lin <- t(Pseudotime_lin)
  colnames(PT_lin) <-
    ifelse(
      round(PT_lin[1,]) %in% marks,
      as.character(round(PT_lin[1,])),
      "")
  
  tiff(
    filename = file.path(
      ptscale_dir,
      paste0("Pseudotime_scale_lineage3_", i, ".tiff")
    ),
    width = 600,
    height = 200
  )
  
  print(
    Heatmap(
      PT_lin,
      name = "Pseudotime",
      column_title = paste0("Pseudotime ", i),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_column_names = TRUE,
      show_row_names = FALSE,
      column_names_side = "bottom",
      col = col_fun
    )
  )
  
  dev.off()
}
