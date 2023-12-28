#### code for FigureS6
library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)


### load the mouse object ###
setwd("~/Downloads/Yi Huang et al objects/")

fHFDcci <- readRDS(file = "fHFDcci.rds")
fLFDcci <- readRDS(file = "fLFDcci.rds")
mHFDcci <- readRDS(file = "mHFDcci.rds")
mLFDcci <- readRDS(file = "mLFDcci.rds")

celltype.colors = c( "AgRP"="#D97176FF",
                     'Astrocytes'="#808000", # repeated colors, changed in the full dataset
                     "Avp"="#f6bb86",
                     "Avp/Rorb"="#D7C1A9FF",
                     'B cells'="#93FFE8", # repeated colors, changed in the full dataset
                     "Col11a1"="#A68872FF",
                     "Col25a1/Six3"="#B496C2FF",
                     "Coro6"="#24325FFF",
                     "Coro6/Trh"="#9CC590FF",
                     "Ebf3/Synpr"="#E89242FF",
                     'Endothelial cells'="#FFDAB9", # repeated colors, changed in the full dataset
                     "Gad2/Slc32a1"="#CC9EA4FF",
                     "Gal"= "#82491EFF",
                     "Gal/Gad2/Slc32a1"="#917C5DFF",
                     "Gal/Ghrh"="#83CBB0FF",
                     "Hcrt"="#876041FF",
                     "Hdc"="#FACCBDFF",
                     'Interneurons'="#A97142", # repeated colors, changed in the full dataset
                     "Lhx1"="#B58900FF",
                     "Lhx1/Coro6"="#526E2DFF",
                     'Microglial/Macrophages'="#191970", ##repeated colors, changed in the full dataset
                     "Npw"="#17692CFF",
                     "Nr5a1"="#2AA198FF",
                     'Oligodendrocytes'="#AB784E",  # repeated colors, changed in the full dataset
                     "Oxt"="#6C71C4FF",
                     "Pmch/Cartpt"="#69C8ECFF",
                     "POMC"="#9698dc",
                     "Sst"="#B69C8AFF", 
                     'Stromal cells'="#967BB6", # repeated colors, changed in the full dataset
                     "Synpr"="#A6EEE6FF",
                     "Tac1/Foxp2"="#E7999DFF",
                     "Tac2-1"="#C6B0A2FF",
                     "Tac2-2"="#CB4B16FF",
                     'Tanycytes'= "#4863A0",# repeated colors, changed in the full dataset
                     "Tcf4/Nr4a2"="#E3DDB1FF",
                     "Trh/Tac1/Bdnf" ="#859900FF",
                     "unassigned"="#4f8c9d",
                     "Vip/Vipr2" ="#54974AFF"
)


####  FigureS6a #### 
object.list <- list(femaleHFD = fHFDcci, 
                    maleHFD = mHFDcci)


fHFDvsmHFD <- mergeCellChat(object.list, 
                          add.names = names(object.list))

fHFDvsmHFD
# An object of class CellChat created from a merged object with multiple datasets 
# 760 signaling genes.
# 19363 cells.

netVisual_heatmap(fHFDvsmHFD,
                  measure = "weight",
                  comparison = c("femaleHFD","maleHFD"),
                  color.use = celltype.colors,
                  title.name = "mHFD vs fHFD")



#### FigureS6b ####
object.list <- list(maleHFD = mHFDcci, 
                    maleLFD = mLFDcci)


mHFDvsmLFD <- mergeCellChat(object.list, 
                          add.names = names(object.list))
mHFDvsmLFD
# An object of class CellChat created from a merged object with multiple datasets 
# 760 signaling genes.
# 19435 cells.

netVisual_heatmap(mHFDvsmLFD,
                  measure = "weight",
                  title.name = "mHFD vs mLFD",
                  color.use = celltype.colors,
                  font.size = 15,
                  font.size.title = 15,
                  row.show = c("AgRP","Astrocytes","POMC","Avp/Rorb","Hdc"),
                  col.show = c("AgRP","Astrocytes","POMC","Avp/Rorb","Hdc"),
                  comparison = c("maleLFD","maleHFD") )


#### FigureS6c ####

object.list <- list(femaleHFD = fHFDcci, 
                    femaleLFD = fLFDcci)


fHFDvsfLFD <- mergeCellChat(object.list, 
                          add.names = names(object.list))

fHFDvsfLFD
# An object of class CellChat created from a merged object with multiple datasets 
# 760 signaling genes.
# 19159 cells.

netVisual_heatmap(fHFDvsfLFD,
                  measure = "weight",
                  title.name = "fHFD vs fLFD",
                  color.use = celltype.colors,
                  font.size = 15,
                  font.size.title = 15,
                  row.show = c("AgRP","Astrocytes","POMC","Avp/Rorb","Hdc"),
                  col.show = c("AgRP","Astrocytes","POMC","Avp/Rorb","Hdc"),
                  comparison = c("femaleLFD","femaleHFD") )


#### FigureS6d ####
p <- rankNet(mHFDvsmLFD, 
               mode = "comparison", 
               stacked = T, 
               color.use = c("#917C5DFF",
                             "#E89242FF"),
               sources.use = c(1,2),
               targets.use = c(1,2),
               font.size = 20,
               do.stat = TRUE)
p


#### FigureS6e ####
netVisual_bubble(mHFDvsmLFD, 
                 sources.use = c(1), 
                 color.text.use = T,
                 color.text = c("#917C5DFF",
                                "#E89242FF"),
                 n.colors = 10,
                 targets.use = c(1:2),  
                 comparison = c(1,2), 
                 angle.x = 45)

#### FigureS6f ####
netVisual_bubble(mHFDvsmLFD, 
                 sources.use = c(2), 
                 color.text.use = T,
                 color.text = c("#917C5DFF",
                                "#E89242FF"),
                 n.colors = 10,
                 targets.use = c(1:2),  
                 comparison = c(1,2), 
                 angle.x = 45)

#### FigureS6g ####
gg1 <- rankNet(fHFDvsfLFD, 
               mode = "comparison", 
               stacked = T, 
               color.use = c('#82491EFF',
                             '#526E2DFF'),
               sources.use = c(1,2),
               targets.use = c(1,2),
               font.size = 20,
               do.stat = TRUE)
gg1


#### FigureS6h ####
netVisual_bubble(fHFDvsfLFD, 
                 sources.use = c(1), 
                 color.text.use = T,
                 color.text = c('#82491EFF',
                                '#526E2DFF'),
                 n.colors = 10,
                 targets.use = c(1:2),  
                 comparison = c(1,2), 
                 angle.x = 45)

#### FigureS6i####
netVisual_bubble(fHFDvsfLFD, 
                 sources.use = c(2), 
                 color.text.use = T,
                 color.text = c('#82491EFF',
                                '#526E2DFF'),
                 n.colors = 10,
                 targets.use = c(1:2),  
                 comparison = c(1,2), 
                 angle.x = 45)

