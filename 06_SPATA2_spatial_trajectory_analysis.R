# Spatial trajectory screening is based on 
# the SPATA2 vignette with project-specific adaptations
# Citation: Kueckelhaus, J., Frerich, S., Kada-Benotmane, J. et al. 
# Inferring histology-associated gene expression gradients in spatial transcriptomic studies. Nat Commun 15, 7280 (2024)

#----- Load required libraries -----#

library(SPATA2)
library(SPATAData)
library(ggplot2)
library(dplyr)
library(shiny)



#----- Step 01: Define project directories -----#


project_dir <- "path/to/project"

data_dir <- file.path(project_dir, "data")
figure_dir <- file.path(project_dir, "figures")

sample_id   <- "sample_name"
sample_file <- "SpataObject.rds"

save_dir <- file.path(data_dir, sample_id)



#----- Step 02: Load SPATA object -----#
# one object per sample
spataObj <- loadSpataObject(
  file.path(save_dir, sample_file)
)


# control plot 
plotSurfaceComparison(
  object = spataObj,
  color_by = c("MKI67", "AQP4"),
  display_image = TRUE)



#----- Step 03: Create spatial trajectories (interactive) -----#
# Opens SPATA2 Shiny interface
# Set name and safe trajectories properly within the shiny interface
# Trajectories have been build as shown in Fig. 6
spataObj <- createSpatialTrajectories(spataObj)



#----- Step 04: Save updated SPATA object -----#

spataObj <- setSpataDir(
  object = spataObj,
  dir = file.path(save_dir, sample_file))

saveSpataObject(spataObj)



#----- Step 05: Plot overview of all trajectories -----#

png(file.path(save_dir, "all_spatial_trajectories.png"),height = 800,width  = 800)

plotSpatialTrajectories(
  object = spataObj,
  color_by = "seurat_clusters",
  pt_alpha = 0.4,
  pt_alpha2 = 1)
dev.off()


#----- Step 06: Retrieve spatial trajectory IDs -----#
traj_ids <- getSpatialTrajectoryIds(spataObj)

print(traj_ids)
length(traj_ids)
# Select trajectories to analyse, e.g. 1st
traj_id <- traj_ids[1]  



#----- Step 08: Ridgeplot screening – neuronal program -----#

genes_neuronal <- c(
  "DCX","SOX4","DLX1","DLX2","DLX5","DLX6",
  "NRXN3","SNAP25","SYP"
)

png(
  file.path(
    figure_dir,
    paste0("trajectory_", traj_id, "_ridgeplot.png")),height = 1200,width  = 1000)

plot(
  plotTrajectoryRidgeplot(
    object = spataObj,
    id = traj_id,
    variables = genes_neuronal,
    alpha = 0.6,
    line_size = 1.1,
    clrp = "plasma",
    smooth_span = 1.1,
    smooth_method = "loess",
    method_gs = "zscore"
  ) +
    legendNone() +
    theme(
      text = element_text(size = 17),
      axis.text.y.right = element_text(angle = 45)
    ))
dev.off()



#-----  Validation ridgeplot as shown in FIG.6 -----#

# for DMG 2
plot(
  plotTrajectoryRidgeplot(
    object = spataObj,
    id = "OPCP",
    variables = c(
      "CSPG4","SOX10","NTRK2","OLIG2",
      "TNR","GPR17","PLP1","MBP"
    ),
    alpha = 0.6,
    line_size = 1.1,
    clrp = "plasma",
    smooth_span = 0.9,
    smooth_method = "loess",
    method_gs = "zscore"
  ) +
    legendNone() +
    theme(
      text = element_text(size = 17),
      axis.text.y.right = element_text(angle = 45)
    )
)



# for WT2 
plot(
  plotTrajectoryRidgeplot(
    object = spataObj,
    id = "AC_3",
    variables = c(
      "TTYH1","GLI3","SOX9","SOX2",
      "SLC1A2","GAP43","APLNR","AQP4"
    ),
    alpha = 0.6,
    line_size = 1.1,
    clrp = "plasma",
    smooth_span = 0.9,
    smooth_method = "loess",
    method_gs = "zscore"
  ) +
    legendNone() +
    theme(
      text = element_text(size = 17),
      axis.text.y.right = element_text(angle = 45)
    )
)

