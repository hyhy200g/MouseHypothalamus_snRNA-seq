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

## Figure4
# figure4A
DefaultAssay(Hypo_neuron) <- "SCT"
FeaturePlot(Hypo_neuron, features = "Agrp", cols = c("lightgrey", "pink")) + theme_void()
FeaturePlot(Hypo_neuron, features = "Npy", cols = c("lightgrey", "pink")) + theme_void()

# figure4B
FeaturePlot(AgRP_subclusters, features = "Agrp", cols = c("lightgrey", "pink")) + theme_void()
FeaturePlot(AgRP_subclusters, features = "Npy", cols = c("lightgrey", "pink")) + theme_void()
DimPlot(AgRP_subclusters, cols = c("#D39394",
                                   "#71838F",
                                   "#E5CC96")) + theme_void()

# figure4C
AgRP_sub_group <- data.frame(group = AgRP_subclusters@meta.data$group)
AgRP_sub_group$subcluster <- AgRP_subclusters@meta.data$seurat_clusters
table(AgRP_sub_group$group, AgRP_sub_group$subcluster)

ggplot(AgRP_sub_group, aes(group, fill=subcluster)) +
  geom_bar(position = "fill", width = 0.9) +
  theme_classic() +
  labs(y = "Cell proportion") +
  scale_fill_manual(values = c("#D39394",
                               "#71838F",
                               "#E5CC96")) +
  RotatedAxis()

# figure4D
AgRP_subclusters_numi <- AgRP_subclusters[["RNA"]]@counts
range(AgRP_subclusters_numi["Xist",]) 


normalize <- function(x){
  return(log(x/sum(x)*10000 + 1))
}
AgRP_subcluster_norm <- apply(AgRP_subclusters_numi, 2, normalize)
AgRP_subcluster_norm[1:5, 1:5]

AgRP_subcluster_expr <- data.frame(cell = colnames(AgRP_subcluster_norm),
                                   Agrp = AgRP_subcluster_norm["Agrp",],
                                   Xist = AgRP_subcluster_norm["Xist",],
                                   Brd1 = AgRP_subcluster_norm["Brd1",],
                                   Hdac3 = AgRP_subcluster_norm["Hdac3",],
                                   group = AgRP_subclusters@meta.data$group)
AgRP_subcluster_expr$group <- factor(AgRP_subcluster_expr$group)
levels(AgRP_subcluster_expr$group) <- c("female_HFD", "female_LFD", "male_HFD", "male_LFD")

ggplot(AgRP_subcluster_expr, aes(x=group, y=Agrp, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))

ggplot(AgRP_subcluster_expr, aes(x=group, y=Xist, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange")) 

ggplot(AgRP_subcluster_expr, aes(x=group, y=Brd1, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange")) +
  geom_jitter(color="lightgrey") 

ggplot(AgRP_subcluster_expr, aes(x=group, y=Hdac3, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange")) +
  geom_jitter(color="lightgrey") 

# figure4E
DimPlot(AvpRorb_subclusters,
        reduction = "umap",
        cols = c("#F2AA4CFF",
                 "#00203FFF"),
        pt.size = 1) +
  theme_classic() +
  theme_void()

# figure4F
marker <- FindAllMarkers(AvpRorb_subclusters, logfc.threshold = 0.25, assay = "SCT")

marker %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(AvpRorb_subclusters) <- "SCT"
DoHeatmap(AvpRorb_subclusters, features = top10$gene, group.colors = c("orange", "black"), label = FALSE) + 
  NoLegend() +
  scale_fill_gradientn(colors=c("#5CACEE","white","red"))

# figure4G
##featureplot
DefaultAssay(AvpRorb_subclusters) <- "SCT"

#### Avp
p0 <- Seurat::FeaturePlot(AvpRorb_subclusters,
                          features = "Avp",
                          cols = c("gray","purple"),
                          label = F,
                          pt.size = 0.5)
p0


#### Rorb
p1 <- Seurat::FeaturePlot(AvpRorb_subclusters,
                          features = "Rorb",
                          cols = c("gray","purple"),
                          label = F,
                          pt.size = 0.5)
p1


#### Igfbp5
p2 <- Seurat::FeaturePlot(AvpRorb_subclusters,
                          features = "Igfbp5",
                          cols = c("gray","purple"),
                          label = F,
                          pt.size = 0.5)
p2


#### Gm42418
p3 <- Seurat::FeaturePlot(AvpRorb_subclusters,
                          features = "Gm42418",
                          cols = c("gray","purple"),
                          label = F,
                          pt.size = 0.5)
p3

##violin plot
AvpRorb_subclusters_numi <- AvpRorb_subclusters[["RNA"]]@counts

normalize <- function(x){
  return(log(x/sum(x)*10000 + 1))
}

AvpRorb_subcluster_norm <- apply(AvpRorb_subclusters_numi, 2, normalize)

AvpRorb_subcluster_expr <- data.frame(cell = colnames(AvpRorb_subcluster_norm),
                                      Avp = AvpRorb_subcluster_norm["Avp",],
                                      Rorb = AvpRorb_subcluster_norm["Rorb",],
                                      Igfbp5 = AvpRorb_subcluster_norm["Igfbp5",],
                                      Gm42418 = AvpRorb_subcluster_norm["Gm42418",],
                                      group = AvpRorb_subclusters@meta.data$group)
AvpRorb_subcluster_expr$group <- factor(AvpRorb_subcluster_expr$group)
levels(AvpRorb_subcluster_expr$group) <- c("female_HFD", "female_LFD", "male_HFD", "male_LFD")

ggplot(AvpRorb_subcluster_expr, aes(x=group, y=Avp, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))

ggplot(AvpRorb_subcluster_expr, aes(x=group, y=Rorb, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))

ggplot(AvpRorb_subcluster_expr, aes(x=group, y=Igfbp5, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))

ggplot(AvpRorb_subcluster_expr, aes(x=group, y=Gm42418, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange"))


