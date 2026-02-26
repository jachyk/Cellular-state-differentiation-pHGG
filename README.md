# Cellular-state-differentiation-pHGG
This repository contains all code developed for Choyka et al. "Developmental trajectories follow a spatially organized hierarchy and contribute to subtype-specific cell identity programs in pediatric high-grade glioma" (unpublished)


**Contents**
- 01_Build_pseudotime_trajectory.R: Setup and create pseudotime object using Monocle3.
- 02_Pseudotime_expression_dynamics.R: Analyze gene expression and cellular density across pseudotime for each subtype. 
- 03_Pseudotime_heatmap_visualization.R: Plot subtypewise pseudotime-dependent gene expression in heatmaps.
- 04_Pseudobin_construction.R: Discretize pseudotime values in Seurat.
- 05_Pseudotime_anchor_transfer.R: Transfer Pseudotime from single-cell-workflows to spatial transcriptomics.
- 06_SPATA2_spatial_trajectory_analysis.R: Validate pseudotime trajectories with spatial trajectory analysis. 
- 07_Bulk_expression_analysis.R: Validate correlation with statistics in external cohort. 
- 08_CNA_inference.R: inference of copy-number-alterations and malignancy scoring 
- 09_Spatial_NMF_clustering: Calculation of spatial NMF clustering (upload in progress) 
