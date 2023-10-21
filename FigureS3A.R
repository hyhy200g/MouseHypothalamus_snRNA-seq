library(Seurat)
library(RColorBrewer)
Hypo_neuron <- load("~/R/project/data/mouse_neural/Hypo_neuron.rds")
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




