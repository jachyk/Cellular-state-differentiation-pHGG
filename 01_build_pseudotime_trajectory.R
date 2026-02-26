# This workflow follows the general trajectory inference strategy
# described in the Monocle3 vignette "Constructing single-cell trajectories"
# with adaptations for integrated pHGG single-cell data.

#Citation: Cole Trapnell*, Davide Cacchiarelli*, Jonna Grimsby, Prapti Pokharel, Shuqiang Li, Michael Morse, Niall J. Lennon, Kenneth J. Livak, Tarjei S. Mikkelsen, and John L. Rinn.
#The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells.
#Nature Biotechnology, 2014.

#----- Load packages -----#
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(openxlsx)


#----- Step 01: Load neural Seurat object -----#

integrated <- readRDS("data/tumor_subcluster3.rds")


#----- Step 02: Convert Seurat object to Monocle3 CDS -----#

# import UMAP cell embeddings for consistency 
umap_coords <- Embeddings(integrated, "umap.snn")

cds <- as.cell_data_set(integrated, assay = "RNA")
cds <- estimate_size_factors(cds)

reducedDim(cds, "UMAP") <- umap_coords

# import gene information
rowData(cds)$gene_short_name <- rownames(integrated[["RNA"]])


#----- Step 03: Cluster cells -----#

cds <- cluster_cells(
  cds,
  reduction_method = "UMAP",
  resolution = 0.0004,
  num_iter = 10
)


#----- Step 04: Learn trajectory graph -----#

cds <- learn_graph(
  cds,
  learn_graph_control = list(
    minimal_branch_len = 8,
    prune_graph = TRUE,
    orthogonal_proj_tip = FALSE
  )
)

#----- Step 05: Define principal root node -----#
# helper function derived from the Monocle3 vignette 
# calculates pseudotime origin of
# cells annotated as "undifferentiated"  
get_earliest_principal_node <- function(cds, cluster = c("1")) {
  
  cell_ids <- which(colData(cds)$SNNGRAPH_res.0.7 %in% cluster)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$
    pr_graph_cell_proj_closest_vertex
  
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[
      as.numeric(names(which.max(
        table(closest_vertex[cell_ids,])
      )))
    ]
  
  root_pr_nodes
}


#----- Step 06: Calculate pseudotime -----#

cds <- order_cells(
  cds,
  root_pr_nodes = get_earliest_principal_node(cds)
)

# store pseudotime into cds metadata
cds$all3_pseudotime <- pseudotime(cds)


#----- Step 07: Save for downstream analysis -----#

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

saveRDS(
  cds,
  "results/objects/final_monocle3_cds.rds"
)


#----- Step 08: Export pseudotime metadata and prepare plotting tables -----#

DF_pseudo <- as.data.frame(colData(cds))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.xlsx(
  DF_pseudo,
  "results/tables/DF_pseudotime_all_cells.xlsx"
)
