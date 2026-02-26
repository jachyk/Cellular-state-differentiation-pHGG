#----- Load packages -----#

library(monocle3)
library(ggplot2)
library(dplyr)
library(reshape2)

#----- color palette for plotting ----'
color_palette2 <- c(
  "DMG" = "#264653",
  "DHG" = "#e9c46a",
  "WT"  = "#E76F51",
  "IHG" = "#2a9d8f"
)



#----- Load Monocle3 object -----#

cds <- readRDS("objects/final_monocle3_cds.rds")

# control plot 
plot_cells(
  cds,
  color_cells_by = "SNNGRAPH_res.0.7",
  show_trajectory_graph = TRUE,
  label_branch_points = TRUE,
  label_roots = TRUE
)


#----- Subtype-specific visualization of pseudotime trajectories -----#
# insert specific string for each subtype which are stored in the metadata
# the new lineage_cds subset contains cells of only one subtype
lineage_cds <- cds[, colData(cds)$subtype %in% "DHG"]

plot_cells(
  lineage_cds,
  color_cells_by = "subtype",
  show_trajectory_graph = TRUE
)

plot_cells(
  lineage_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)



#-----  Define lineage clusters for downstream analysis and subsetting -----#

lineage1<-c("1","4","10","13")    #NPC 
lineage2<-c("1","0","17","18")    #OPC
lineage3<-c("1","6","3","9","5")  #AC
lineage4<-c("1","2","5","7")      #Reactive

#subset object by desired lineage

lineage_cds <- cds[, colData(cds)$SNNGRAPH_res.0.7 %in% lineage1] 

plot_cells(
  lineage_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = FALSE
)



#----- Plot gene expression along pseudotime -----#


# Define genes of interest
# e.g. following astrocyte related genes
genes1 <- c(
  "GLI3","SOX9","AQP4",
  "SLC1A2","F3","SPARCL1","FGFR3"
)

# convert data into dataframe format for plotting
cell_metadata <- as.data.frame(colData(lineage_cds))
pseudotime_condition_data <-
  cell_metadata[, c("all3_pseudotime","subtype")]


dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

for (gene in genes1) {
  
  expression_data <-
    as.data.frame(as.matrix(counts(lineage_cds)[gene, ]))
  colnames(expression_data) <- gene
  combined_data <-
    cbind(pseudotime_condition_data, expression_data)
  long_data <- melt(
    combined_data,
    id.vars = c("all3_pseudotime","subtype"),
    variable.name = "gene",
    value.name = "expression"
  )
  
  p <- ggplot(
    long_data,
    aes(
      x = all3_pseudotime,
      y = expression,
      color = subtype)) +
    geom_smooth(method = "loess", linewidth = 1.5, se = FALSE) +
    scale_colour_manual(values = color_palette2) +
    theme_minimal()
  
  ggsave(
    filename = paste0("results/figures/", gene, "_pseudotime.pdf"),
    plot = p,
    width = 6,
    height = 4)
}


#----- Visualize cellular density along pseudotime -----#

ggplot(pseudotime_condition_data, aes(x=all_pseudotime, group = subtype, colour = subtype))+
  geom_vline(xintercept = 10, linetype = "longdash", colour= "#495057", linewidth = 1.3) +
  geom_density(adjust = 1.4, alpha = .5, linewidth = 2.1) +
  theme_minimal()+ scale_colour_manual(values = c("#e9c46a","#264653","#2a9d8f","#d62828"))+ 
  labs(x = "", y = "", title = "")+theme(axis.text  = element_text(size=18))+ theme(legend.position = "none")
