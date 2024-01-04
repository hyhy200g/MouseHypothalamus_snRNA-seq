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

## Figure2
# figure2B
Hypo_integrated_umap <- as.data.frame(Hypo_integrated@reductions$umap@cell.embeddings)
Hypo_integrated_umap$Celltype <- Hypo_integrated@meta.data$Celltype
Hypo_integrated_Celltype <- as.data.frame(Hypo_integrated@meta.data$Celltype)
ggplot(data = Hypo_integrated_umap, mapping = aes(x = UMAP_1, y = UMAP_2, color = Celltype)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  theme_bw() +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c('Astrocytes'="#526E2DFF",
                                'Endothelial cells'="#f6bb86", 
                                'Interneurons'="#82491EFF",
                                'B cells'="#69C8ECFF",
                                'Microglial/Macrophages'="#24325FFF",
                                'Neurons'="#E89242FF", 
                                'Oligodendrocytes'="#917C5DFF", 
                                'Stromal cells'="#9698dc",
                                'Tanycytes'= "#4f8c9d"))

# figure2E
Hypo_mHFD <- subset(Hypo_integrated, sex=="male" & diet=="HFD")
Hypo_fHFD <- subset(Hypo_integrated, sex=="female" & diet=="HFD")
Hypo_mLFD <- subset(Hypo_integrated, sex=="male" & diet=="LFD")
Hypo_fLFD <- subset(Hypo_integrated, sex=="female" & diet=="LFD")

Hypo_fHFD_freq <- as.data.frame(table(Hypo_fHFD@meta.data$Celltype))
ggdonutchart(Hypo_fHFD_freq, "Freq", fill="Var1",
             ggtheme = theme_void(),
             palette = c("#526E2DFF",
                         "#f6bb86",
                         "#82491EFF",
                         "#69C8ECFF",
                         "#24325FFF",
                         "#E89242FF",
                         "#917C5DFF",
                         "#9698dc",
                         "#4f8c9d")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) 

Hypo_fLFD_freq <- as.data.frame(table(Hypo_fLFD@meta.data$Celltype))
ggdonutchart(Hypo_fLFD_freq, "Freq", fill="Var1",
             ggtheme = theme_void(),
             palette = c("#526E2DFF",
                         "#f6bb86",
                         "#82491EFF",
                         "#69C8ECFF",
                         "#24325FFF",
                         "#E89242FF",
                         "#917C5DFF",
                         "#9698dc",
                         "#4f8c9d")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) 

Hypo_mHFD_freq <- as.data.frame(table(Hypo_mHFD@meta.data$Celltype))
ggdonutchart(Hypo_mHFD_freq, "Freq", fill="Var1",
             ggtheme = theme_void(),
             palette = c("#526E2DFF",
                         "#f6bb86",
                         "#82491EFF",
                         "#69C8ECFF",
                         "#24325FFF",
                         "#E89242FF",
                         "#917C5DFF",
                         "#9698dc",
                         "#4f8c9d")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) 

Hypo_mLFD_freq <- as.data.frame(table(Hypo_mLFD@meta.data$Celltype))
ggdonutchart(Hypo_mLFD_freq, "Freq", fill="Var1",
             ggtheme = theme_void(),
             palette = c("#526E2DFF",
                         "#f6bb86",
                         "#82491EFF",
                         "#69C8ECFF",
                         "#24325FFF",
                         "#E89242FF",
                         "#917C5DFF",
                         "#9698dc",
                         "#4f8c9d")) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank()) 

# figure2C
DefaultAssay(Hypo_integrated) <- "SCT"

FeaturePlot(Hypo_integrated, features="Snhg11", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Rbfox3", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Ntsr2", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Col3a1", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Siglech", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Mrc1", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Pdgfra", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Plp1", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Col23a1", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Flt1", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Cga", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()
FeaturePlot(Hypo_integrated, features="Cd79a", cols = c("lightgrey", "#E89242FF"), pt.size = 0.3) + theme_void()

# figure2D
DotPlot(Hypo_integrated, features = c("Meg3", "Snhg11", "Slc1a3", "Ntsr2", "C1qb",
                                      "Mrc1", "Pecam1", "Flt1", "Col23a1", "Dnah12",
                                      "Col1a1", "Col3a1", "Cga", "Chga", "Plp1",
                                      "Sirt2", "Ebf1", "Cd79a"), cols = c("lightgrey", "orange"),
        dot.scale = 6, cluster.idents = F) + RotatedAxis() + ggplot2::coord_flip()
