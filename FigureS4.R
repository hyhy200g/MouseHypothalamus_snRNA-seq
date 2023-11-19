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

## FigureS4
# figureS4a
marker <- FindAllMarkers(AgRP_subclusters, logfc.threshold = 0.25, assay = "SCT")

marker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(AgRP_subclusters) <- "SCT"
DoHeatmap(AgRP_subclusters, features = top10$gene, group.colors = c("orange", "black"), label = FALSE) + 
  NoLegend() +
  scale_fill_gradientn(colors=c("#5CACEE","white","red"))

# figureS4B
AgRP_subclusters_numi <- AgRP_subclusters[["RNA"]]@counts

normalize <- function(x){
  return(log(x/sum(x)*10000 + 1))
}
AgRP_subcluster_norm <- apply(AgRP_subclusters_numi, 2, normalize)
AgRP_subcluster_norm[1:5, 1:5]

AgRP_subcluster_expr <- data.frame(cell = colnames(AgRP_subcluster_norm),
                                   Agrp = AgRP_subcluster_norm["Agrp",],
                                   Npy = AgRP_subcluster_norm["Npy",],
                                   subcluster = AgRP_subclusters@meta.data$seurat_clusters)

ggplot(AgRP_subcluster_expr, aes(x=subcluster, y=Agrp, fill=subcluster)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="black", size=3) +
  scale_fill_manual(values = c("pink", "grey", "#EEC591")) + coord_flip()

ggplot(AgRP_subcluster_expr, aes(x=subcluster, y=Npy, fill=subcluster)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="black", size=3) +
  scale_fill_manual(values = c("pink", "grey", "#EEC591")) + coord_flip()

# figureS4C
p <- FeaturePlot(AgRP_subclusters, features = c("Agrp", "Lepr"), cols = c("lightgrey", "red", "blue"),
                 blend = T, blend.threshold = 0)

# figureS4D
FeaturePlot(Hypo_neuron, features = c("Avp", "Rorb"), cols = c("lightgrey", "blue", "red"),
            blend = T, blend.threshold = 0)

# figureS4E
FeaturePlot(Hypo_neuron, features = "Pomc", cols = c("lightgrey", "hotpink4")) + theme_void()

# figureS4F
DimPlot(POMC, cols = c("purple", "darkgreen", "blue")) + theme_void()

POMC_group <- data.frame(subcluster = POMC@meta.data$seurat_clusters,
                         group = POMC@meta.data$group)
ggplot(POMC_group, aes(group, fill=subcluster)) +
  geom_bar(position = "fill", width = 0.9) +
  theme_classic() +
  labs(y = "Cell proportion") +
  scale_fill_manual(values = c("purple",
                               "darkgreen",
                               "blue"))

# figureS4G
marker <- FindAllMarkers(POMC, logfc.threshold = 0.25, assay = "SCT")

marker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(POMC) <- "SCT"
DoHeatmap(POMC, features = top10$gene, group.colors = c("orange", "black"), label = FALSE) + 
  NoLegend() +
  scale_fill_gradientn(colors=c("#5CACEE","white","red"))

# figureS4H
POMC_subclusters_umi <- POMC[["RNA"]]@counts

normalize <- function(x){
  return(log(x/sum(x)*10000 + 1))
}
POMC_subcluster_norm <- apply(POMC_subclusters_umi, 2, normalize)
POMC_subcluster_norm[1:5, 1:5]

POMC_subcluster_expr <- data.frame(cell = colnames(POMC_subcluster_norm),
                                   Pomc = POMC_subcluster_norm["Pomc",],
                                   Xist = POMC_subcluster_norm["Xist",],
                                   Brd1 = POMC_subcluster_norm["Brd1",],
                                   group = POMC@meta.data$group)

ggplot(POMC_subcluster_expr, aes(x=group, y=Pomc, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="black", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))

ggplot(POMC_subcluster_expr, aes(x=group, y=Xist, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="black", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))

ggplot(POMC_subcluster_expr, aes(x=group, y=Brd1, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="black", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))