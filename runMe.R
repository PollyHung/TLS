library(magrittr)
library(annotateTLS)
library(dplyr)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(harmony)
library(scCustomize)

# devtools::install_github("PollyHung/annotateTLS")

## Provide Basic Paths
path <- "~/Desktop/TLS/data/raw-spatial/E-MTAB-13084/spatial"
samples <- list.dirs(path, recursive = F, full.names = F)

## Run Preprocessing（¯﹃¯）
seurat <- preprocess(samples = samples, path = path)

## Plot d(d＇∀＇)
setwd(path)
VlnPlot(seurat, features = c("nFeature_SCT", "nCount_SCT"), group.by = "orig.ident", pt.size = 0)
ggsave(filename = "InterSampleDifference.png", width = 8, height = 4, dpi = 600, units = "in")
FeaturePlot_scCustom(seurat, features = c("TLS", "B.cell", "T.cell", "fDC"))
ggsave(filename = "FeaturePlot.png", width = 7, height = 6, dpi = 600, units = "in")
DimPlot_scCustom(seurat, group.by = "orig.ident", colors_use = "Set2")
ggsave(filename = "DimPlot.png", width = 5, height = 3, dpi = 600, units = "in")
SpatialFeaturePlot(seurat, features = "TLS", pt.size.factor = 2, ncol = 4, crop = FALSE)
ggsave(filename = "SpatialFeaturePlot-TLS.png", width = 16, height = 8, dpi = 600, units = "in")

## KNN!! (finally reaching this step)










