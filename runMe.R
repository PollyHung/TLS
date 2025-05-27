library(magrittr)
library(annotateTLS)
library(dplyr)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(harmony)
library(scCustomize)

# devtools::install_github("PollyHung/annotateTLS")

# ## Provide Basic Paths
# path <- "~/Desktop/TLS/data/raw-spatial/E-MTAB-13530/spatial"
# samples <- list.dirs(path, recursive = F, full.names = F)
#
# ## Run Preprocessing（¯﹃¯）
# seurat <- preprocess(samples = samples, path = path)
#
# ## Plot d(d＇∀＇)
# setwd(path)
# VlnPlot(seurat, features = c("nFeature_SCT", "nCount_SCT"), group.by = "orig.ident", pt.size = 0)
# ggsave(filename = file.path(path, "InterSampleDifference.png"), width = 8, height = 4, dpi = 600, units = "in")
# FeaturePlot_scCustom(seurat, features = c("TLS", "B.cell", "T.cell", "fDC"))
# ggsave(filename = file.path(path, "FeaturePlot.png"), width = 7, height = 6, dpi = 600, units = "in")
# DimPlot_scCustom(seurat, group.by = "orig.ident", colors_use = "Set2")
# ggsave(filename = file.path(path, "DimPlot.png"), width = 5, height = 3, dpi = 600, units = "in")
# SpatialFeaturePlot(seurat, features = "TLS", pt.size.factor = 2, ncol = 4, crop = FALSE, min.cutoff = 0)
# ggsave(filename = file.path(path, "SpatialFeaturePlot-TLS.png"), width = 16, height = 8, dpi = 600, units = "in")
#
# ## KNN!! (finally reaching this step)
# seurat2 <- runKNN(seurat = seurat, exp_threshold = 0.96, min_spots = 5)
# saveRDS(seurat2, "seurat.rds")
#
# ## Plot 2 d(d＇∀＇)
# SpatialDimPlot(seurat2, group.by = "TLS_identity", pt.size.factor = 2, ncol = 4, crop = FALSE)
# ggsave(filename = file.path(path, "SpatialFeaturePlot-TLS-identity.png"), width = 16, height = 8, dpi = 600, units = "in")
# SpatialDimPlot(seurat2, group.by = "TLS_identity", pt.size.factor = 2, ncol = 4, crop = FALSE, image.alpha = 0)
# ggsave(filename = file.path(path, "SpatialFeaturePlot-TLS-identity-clean.png"), width = 16, height = 8, dpi = 600, units = "in")
# SpatialDimPlot(seurat2, group.by = "TLS_identity", pt.size.factor = 0, ncol = 4, crop = FALSE)
# ggsave(filename = file.path(path, "SpatialFeaturePlot-H&E.png"), width = 16, height = 8, dpi = 600, units = "in")














