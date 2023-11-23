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

## Figure3
#figure3A
FeaturePlot(Hypo_neuron, features = c("Slc32a1"), cols = c("lightgrey", "#FF4500")) + theme_void()
FeaturePlot(Hypo_neuron, features = c("Gad1"), cols = c("lightgrey", "#FF4500")) + theme_void()
FeaturePlot(Hypo_neuron, features = c("Gad2"), cols = c("lightgrey", "#FF4500")) + theme_void()
FeaturePlot(Hypo_neuron, features = c("Slc17a6"), cols = c("lightgrey", "#FFA500")) + theme_void()
FeaturePlot(Hypo_neuron, features = c("Hdc"), cols = c("lightgrey", "#8B8B7A")) + theme_void()

#figure3B
Hypo_neuron_umap <- as.data.frame(Hypo_neuron@reductions$umap@cell.embeddings)
Hypo_neuron_umap$subcluster  <- Hypo_neuron@meta.data$Celltype
ggplot(data = Hypo_neuron_umap, mapping = aes(x = UMAP_1, y = UMAP_2, color = subcluster)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  theme_bw() +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("Lhx1/Coro6"="#526E2DFF",
                                "Lhx1"="#B58900FF",
                                "Coro6/Trh"="#9CC590FF",
                                "Coro6"="#24325FFF",
                                "Trh/Tac1/Bdnf" ="#859900FF",
                                "Tac2-1"="#C6B0A2FF",
                                "Tac2-2"="#CB4B16FF",
                                "Tac1/Foxp2"="#E7999DFF",
                                "Tcf4/Nr4a2"="#E3DDB1FF",
                                "Ebf3/Synpr"="#E89242FF",
                                "Synpr"="#A6EEE6FF",
                                "Col11a1"="#A68872FF",
                                "Gad2/Slc32a1"="#CC9EA4FF",
                                "Gal/Gad2/Slc32a1"="#917C5DFF",
                                "Gal"= "#82491EFF",
                                "Gal/Ghrh"="#83CBB0FF",
                                "Vip/Vipr2" ="#54974AFF",
                                "Avp/Rorb"="#D7C1A9FF",
                                "Avp"="#f6bb86",
                                "Oxt"="#6C71C4FF",
                                "Sst"="#B69C8AFF",
                                "Col25a1/Six3"="#B496C2FF",
                                "Nr5a1"="2AA198FF",
                                "Hcrt"="#876041FF",
                                "POMC"="#9698dc",
                                "AgRP"="#D97176FF",
                                "Npw"="#17692CFF",
                                "Hdc"="#FACCBDFF",
                                "Pmch/Cartpt"="#69C8ECFF",
                                "unassigned"="#4f8c9d"))

#figure3C
Hypo_neuron@meta.data$Celltype <- factor(Hypo_neuron@meta.data$Celltype,
                                         levels = c("unassigned", "Pmch/Cartpt", "Hdc", "Npw", "AgRP",
                                                    "POMC", "Hcrt", "Nr5a1", "Col25a1/Six3", "Sst",
                                                    "Oxt", "Avp", "Avp/Rorb", "Vip/Vipr2", "Gal/Ghrh",
                                                    "Gal", "Gal/Gad2/Slc32a1", "Gad2/Slc32a1", "Col11a1",
                                                    "Synpr", "Ebf3/Synpr", "Tcf4/Nr4a2", "Tac1/Foxp2", 
                                                    "Tac2-2", "Tac2-1", "Trh/Tac1/Bdnf", "Coro6", 
                                                    "Coro6/Trh", "Lhx1", "Lhx1/Coro6"))

DotPlot(Hypo_neuron, features = c("Lhx1", "Coro6", "Trh", "Tac1", "Foxp2",
                                  "Bdnf", "Tac2", "Tcf4", "Nr4a2", "Ebf3",
                                  "Synpr", "Col11a1", "Gad2", "Slc32a1", "Gal",
                                  "Ghrh", "Vip", "Vipr2", "Avp", "Rorb",
                                  "Oxt", "Sst", "Col25a1", "Six3", "Nr5a1",
                                  "Hcrt", "Pomc", "Agrp", "Npw", "Hdc",
                                  "Pmch", "Cartpt"), cols = c("lightgrey", "red"),
        dot.scale = 5, cluster.idents = F, group.by = "Celltype") + RotatedAxis()

#figure3D
