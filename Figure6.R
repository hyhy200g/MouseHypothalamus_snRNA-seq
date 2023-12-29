#### code for Figure 6

library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)


#### load the objects ####
setwd("~/Downloads/Yi Huang et al objects/")

fHFDcci <- readRDS(file = "fHFDcci.rds")
fLFDcci <- readRDS(file = "fLFDcci.rds")
mHFDcci <- readRDS(file = "mHFDcci.rds")
mLFDcci <- readRDS(file = "mLFDcci.rds")

object.list <- list(femaleHFD = fHFDcci, 
                    femaleLFD = fLFDcci,
                    maleHFD = mHFDcci,
                    maleLFD = mLFDcci)


cellchat <- mergeCellChat(object.list, 
                          add.names = names(object.list))
cellchat
# An object of class CellChat created from a merged object with multiple datasets 
# 760 signaling genes.
# 38594 cells.

group.color <- c( 'femaleHFD'='#82491EFF',
                  'femaleLFD'='#526E2DFF',
                  'maleHFD'="#917C5DFF",
                  'maleLFD'="#E89242FF"
)

#### Figure 6a ####
gg1 <- compareInteractions(cellchat, 
                           measure = "count",
                           show.legend = F,
                           group = c(1,2,3,4),
                           color.use = c('#82491EFF',
                                         '#526E2DFF',
                                         "#917C5DFF",
                                         "#E89242FF"),
                           angle.x = 90,
                           size.text = 20)
gg1



#### Figure 6b ####

###fHFD vs mHFD comparison ###
object.list <- list(femaleHFD = fHFDcci, 
                    maleHFD = mHFDcci)
fHFDvsmHFD <- mergeCellChat(object.list, 
                          add.names = names(object.list))

fHFDvsmHFD
# An object of class CellChat created from a merged object with multiple datasets 
# 760 signaling genes.
# 19363 cells.


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

str(fHFDvsmHFD)


netVisual_heatmap(fHFDvsmHFD,
                  measure = "weight",
                  title.name = "mHFD vs fHFD",
                  color.use = celltype.colors,
                  font.size = 15,
                  font.size.title = 15,
                  row.show = c("AgRP","Astrocytes","POMC","Avp/Rorb","Hdc"),
                  col.show = c("AgRP","Astrocytes","POMC","Avp/Rorb","Hdc"),
                  comparison = c("femaleHFD","maleHFD") )



#### Figure 6e ####
gg1 <- rankNet(cellchat, 
               mode = "comparison", 
               stacked = T, 
               color.use = c('#82491EFF','#917C5DFF'),
               sources.use = c(1,2),
               targets.use = c(1,2),
               font.size = 20,
               do.stat = TRUE)
gg1


#### Figure 6f ####
# LR AgRP to AgRP and Astrocytes
netVisual_bubble(cellchat, 
                 sources.use = c(1), 
                 color.text.use = T,
                 color.text = c('#82491EFF','#917C5DFF'),
                 n.colors = 10,
                 targets.use = c(1:2),  
                 comparison = c(1,2), 
                 angle.x = 45)


#### Figure 6g ####
# LR Astrocytes to AgRP and Astrocytes
netVisual_bubble(cellchat, 
                 sources.use = c(2), 
                 color.text.use = T,
                 color.text = c('#82491EFF','#917C5DFF'),
                 n.colors = 10,
                 targets.use = c(1:2),  
                 comparison = c(1,2), 
                 angle.x = 45)

