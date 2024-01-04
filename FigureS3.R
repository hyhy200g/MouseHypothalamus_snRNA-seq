## import package
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(sctransform)
library(harmony)
library(dplyr)
library(plyr)
library(DAseq)


## load all the objects
Hypo_integrated <- readRDS("~/Yi_Huang/objects/Hypo_integrated.rds")
AgRP_subclusters <- readRDS("~/Yi_Huang/objects/AgRP_subclusters.rds")
Hypo_neuron <- readRDS("~/Yi_Huang/objects/Hypo_neuron.rds")
AvpRorb_subclusters <- readRDS("~/Yi_Huang/objects/AvpRorb_subclusters.rds")
Astrocyte <- readRDS("~/Yi_Huang/objects/astrocyte_subclusters.rds")
ODC <- readRDS("~/Yi_Huang/objects/ODC_subclusters.rds")
POMC <- readRDS("~/Yi_Huang/objects/POMC_subclusters.rds")
Tanycyte <- readRDS("~/Yi_Huang/objects/Tanycyte_subclusters.rds")

## FigureS3
# figureS3A
table(Hypo_neuron@meta.data$Celltype)
meta.data <- Hypo_neuron@meta.data
meta.data$region <- "Hypo"
meta.data$region[meta.data$Celltype == "Lhx1"] <- "SCN"
meta.data$region[meta.data$Celltype == "Trh/Tac1/Bdnf"] <- "PVH"
meta.data$region[meta.data$Celltype == "Tac1/Foxp2"] <- "PVH"
meta.data$region[meta.data$Celltype == "Col11a1"] <- "SPN"
meta.data$region[meta.data$Celltype == "Gad2/Slc32a1"] <- "CeA"
meta.data$region[meta.data$Celltype == "Gal/Gad2/Slc32a1"] <- "CeA"
meta.data$region[meta.data$Celltype == "Gal"] <- "ARC(CeA)" 
meta.data$region[meta.data$Celltype == "Gal/Ghrh"] <- "ARC"
meta.data$region[meta.data$Celltype == "Vip/Vipr2"] <- "SCN"
meta.data$region[meta.data$Celltype == "Avp/Rorb"] <- "SCN"
meta.data$region[meta.data$Celltype == "Avp"] <- "PVN(SCN)" 
meta.data$region[meta.data$Celltype == "Oxt"] <- "PVN"
meta.data$region[meta.data$Celltype == "Sst"] <- "CeA(LHA)" 
meta.data$region[meta.data$Celltype == "Col25a1/Six3"] <- "SCN"
meta.data$region[meta.data$Celltype == "Nr5a1"] <- "VMH"
meta.data$region[meta.data$Celltype == "Hcrt"] <- "LH(LHA)"
meta.data$region[meta.data$Celltype == "POMC"] <- "ARC"
meta.data$region[meta.data$Celltype == "AgRP"] <- "ARC"
meta.data$region[meta.data$Celltype == "Npw"] <- "DMH"
meta.data$region[meta.data$Celltype == "Hdc"] <- "TMN"
meta.data$region[meta.data$Celltype == "Pmch/Cartpt"] <- "LHA"
table(meta.data$region)

Hypo_neuron@meta.data <- meta.data
getPalette = colorRampPalette(brewer.pal(9, "Set1"))(14)
DimPlot(Hypo_neuron, reduction = "umap", pt.size = 0.1, group.by = "region", label = TRUE) +
  scale_color_manual(values = c('ARC'="#526E2DFF",
                                'ARC(CeA)'="#f6bb86", 
                                'CeA'="#0000FF",
                                'CeA(LHA)'="#69C8ECFF",
                                'SCN'="#EEE685",
                                'LH(LHA)'="#E89242FF", 
                                'LHA'="#917C5DFF", 
                                'PVH'="#9698dc",
                                'PVN'= "#4f8c9d",
                                'Hypo' = "#D1D1D1",
                                'PVN(SCN)' = "#7CFC00",
                                'SPN' = "#CD69C9",
                                'TMN' = "#FF8C69",
                                'VMH' = "#CD0000",
                                'DMH' = "#FF69B4")) +
  theme_void() +
  theme(title = element_blank())

ggsave(device = "tiff",
       file = "~/R/project/data_output/mouse.neuron.region.label.tiff",
       width = 10, height = 8,
       units = "in", dpi = 300)

DimPlot(Hypo_neuron, reduction = "umap", pt.size = 0.1, group.by = "region", label = FALSE) +
  scale_color_manual(values = c('ARC'="#526E2DFF",
                                'ARC(CeA)'="#f6bb86", 
                                'CeA'="#0000FF",
                                'CeA(LHA)'="#69C8ECFF",
                                'SCN'="#EEE685",
                                'LH(LHA)'="#E89242FF", 
                                'LHA'="#917C5DFF", 
                                'PVH'="#9698dc",
                                'PVN'= "#4f8c9d",
                                'Hypo' = "#D1D1D1",
                                'PVN(SCN)' = "#7CFC00",
                                'SPN' = "#CD69C9",
                                'TMN' = "#FF8C69",
                                'VMH' = "#CD0000",
                                'DMH' = "#FF69B4")) +
  theme_void() +
  theme(title = element_blank())

ggsave(device = "tiff",
       file = "~/R/project/data_output/mouse.neuron.region.nolabel.tiff",
       width = 10, height = 8,
       units = "in", dpi = 300)
save(Hypo_neuron,
     file = "~/R/project/data/mouse_neural/Hypo_neuron_region.rds")

# figureS3B
Neuron_sample <- data.frame(group = Hypo_neuron@meta.data$group,
                            subcluster = Hypo_neuron@meta.data$Celltype,
                            Cell = colnames(Hypo_neuron))
table(Neuron_sample$sample, Neuron_sample$subcluster)
Neuron_sample$subcluster <- factor(Neuron_sample$subcluster)
levels(Neuron_sample$subcluster)
Neuron_sample$subcluster <- factor(Neuron_sample$subcluster, 
                                   levels = c("Lhx1/Coro6", 
                                              "Lhx1",
                                              "Coro6/Trh",
                                              "Coro6",
                                              "Trh/Tac1/Bdnf",
                                              "Tac2-1",
                                              "Tac2-2",
                                              "Tac1/Foxp2",
                                              "Tcf4/Nr4a2",
                                              "Ebf3/Synpr",
                                              "Synpr",
                                              "Col11a1",
                                              "Gad2/Slc32a1",
                                              "Gal/Gad2/Slc32a1",
                                              "Gal",
                                              "Gal/Ghrh",
                                              "Vip/Vipr2",
                                              "Avp/Rorb",
                                              "Avp",
                                              "Oxt",
                                              "Sst",
                                              "Col25a1/Six3",
                                              "Nr5a1",
                                              "Hcrt",
                                              "POMC",
                                              "AgRP",
                                              "Npw",
                                              "Hdc",
                                              "Pmch/Cartpt",
                                              "unassigned"))

ggplot(Neuron_sample, aes(group, fill=subcluster)) +
  geom_bar(position = "fill", width = 0.9) +
  theme_classic() +
  labs(y = "Cell proportion") +
  scale_fill_manual(values = c("Lhx1/Coro6"="#526E2DFF",
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



# figureS3C
DefaultAssay(Hypo_neuron) <- "SCT"

Hypo_neuron@meta.data$Celltype <- factor(Hypo_neuron@meta.data$Celltype,
                                         levels = c("unassigned", "Pmch/Cartpt", "Hdc", "Npw", "AgRP",
                                                    "POMC", "Hcrt", "Nr5a1", "Col25a1/Six3", "Sst",
                                                    "Oxt", "Avp", "Avp/Rorb", "Vip/Vipr2", "Gal/Ghrh",
                                                    "Gal", "Gal/Gad2/Slc32a1", "Gad2/Slc32a1", "Col11a1",
                                                    "Synpr", "Ebf3/Synpr", "Tcf4/Nr4a2", "Tac1/Foxp2", 
                                                    "Tac2-2", "Tac2-1", "Trh/Tac1/Bdnf", "Coro6", 
                                                    "Coro6/Trh", "Lhx1", "Lhx1/Coro6"))

DotPlot(Hypo_neuron, features = c("Adcyap1", "Agrp", "Avp", "Bdnf", "Cartpt",
                                  "Cck", "Crh", "Gal", "Ghrh", "Gnrh1",
                                  "Grp", "Hcrt", "Nmu", "Nms", "Nmb",
                                  "Npvf", "Npw", "Npy", "Nts", "Oxt",
                                  "Pdyn", "Penk", "Pmch", "Pnoc", "Pomc",
                                  "Prok2", "Rln1", "Sst", "Tac1", "Tac2",
                                  "Trh", "Vip"), cols = c("lightgrey", "red"),
        dot.scale = 5, cluster.idents = F, group.by = "Celltype") + RotatedAxis()

# figureS3D
#### fLFD vs mLFD ####
run_DAseq_comparison(Hypo_neuron,
                     groupA = "female_LFD", 
                     groupB = "male_LFD",
                     labels.1 = c("m39f","m47f"), 
                     labels.2 = c("m36m","m44m"), 
                     path = "./",
                     outputname = "fLFD vs mLFD neuron.pdf")

#### mHFD vs mLFD ####
run_DAseq_comparison(Hypo_neuron,
                     groupA = "male_HFD", 
                     groupB = "male_LFD",
                     labels.1 = c("m50m","m71m"), 
                     labels.2 = c("m36m","m44m"), 
                     path = "./",
                     outputname = "mHFD vs mLFD neuron.pdf")

#### fHFD vs fLFD ####
run_DAseq_comparison(Hypo_neuron,
                     groupA = "female_HFD", 
                     groupB = "female_LFD",
                     labels.1 = c("m48f","m6f"), 
                     labels.2 = c("m39f","m47f"), 
                     path = "./",
                     outputname = "fHFD vs fLFD neuron.pdf")
