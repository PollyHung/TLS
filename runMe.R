library(magrittr)
library(annotateTLS)
library(dplyr)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(harmony)
library(scCustomize)

# ## Provide Basic Paths
# path <- "~/Desktop/TLS/data/raw-spatial/E-MTAB-13084/spatial"
# samples <- list.files(path)
#
# ## Run Preprocessing
# seurat <- preprocess(samples = samples,
#                      path = path)
#
# ## Plot :D
# setwd(path)
# FeaturePlot_scCustom(seurat, features = c("TLS", "B.cell", "T.cell", "fDC"))
# ggsave(filename = "FeaturePlot.png", width = 7, height = 6, dpi = 600, units = "in")
# DimPlot_scCustom(seurat, group.by = "orig.ident", colors_use = "Set2")
# ggsave(filename = "DimPlot.png", width = 5, height = 3, dpi = 600, units = "in")
# SpatialFeaturePlot(seurat, features = "TLS", pt.size.factor = 2, ncol = 4, crop = FALSE)
# ggsave(filename = "SpatialFeaturePlot-TLS.png", width = 16, height = 8, dpi = 600, units = "in")
#
# ## KNN!! (finally reaching this step)
#
#








