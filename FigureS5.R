## import package
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(sctransform)
library(harmony)
library(dplyr)
library(plyr)

## load all the objects
Hypo_integrated <- readRDS("~/Yi_Huang/objects/Hypo_integrated.rds")
AgRP_subclusters <- readRDS("~/Yi_Huang/objects/AgRP_subclusters.rds")
Hypo_neuron <- readRDS("~/Yi_Huang/objects/Hypo_neuron.rds")
AvpRorb_subclusters <- readRDS("~/Yi_Huang/objects/AvpRorb_subclusters.rds")
Astrocyte <- readRDS("~/Yi_Huang/objects/astrocyte_subclusters.rds")
ODC <- readRDS("~/Yi_Huang/objects/ODC_subclusters.rds")
POMC <- readRDS("~/Yi_Huang/objects/POMC_subclusters.rds")
Tanycyte <- readRDS("~/Yi_Huang/objects/Tanycyte_subclusters.rds")

## figureS5
# figureS5A
DefaultAssay(Hypo_integrated) <- "SCT"
FeaturePlot(Hypo_integrated, features = "Ntsr2", cols = c("lightgrey", "darkgreen")) + theme_void()
FeaturePlot(Hypo_integrated, features = "Slc7a10", cols = c("lightgrey", "darkgreen")) + theme_void()
FeaturePlot(Hypo_integrated, features = "Gfap", cols = c("lightgrey", "darkgreen")) + theme_void()

# figureS5B
DimPlot(Tanycyte, cols = c("pink", "hotpink3", "purple")) + theme_void()
FeaturePlot(Tanycyte, features = c("Dnah12"), cols = c("lightgrey", "purple")) + theme_void()
FeaturePlot(Tanycyte, features = c("Ccdc153"), cols = c("lightgrey", "purple")) + theme_void()
FeaturePlot(Tanycyte, features = c("Frzb"), cols = c("lightgrey", "purple")) + theme_void()
FeaturePlot(Tanycyte, features = c("Vcan"), cols = c("lightgrey", "purple")) + theme_void()
FeaturePlot(Tanycyte, features = c("Col25a1"), cols = c("lightgrey", "purple")) + theme_void()
FeaturePlot(Tanycyte, features = c("Scn7a"), cols = c("lightgrey", "purple")) + theme_void()

# figureS5C
marker <- FindAllMarkers(Tanycyte, logfc.threshold = 0.25, assay = "SCT")

marker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(Tanycyte) <- "SCT"
DoHeatmap(Tanycyte, features = top10$gene, group.colors = c("orange", "black"), label = FALSE) + 
  NoLegend() +
  scale_fill_gradientn(colors=c("#5CACEE","white","red"))