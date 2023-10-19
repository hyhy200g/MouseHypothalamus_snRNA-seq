setwd("~/R/project/entropy/CytoTRACE")
library(CytoTRACE)
library(SingCellaR)
library(monocle3)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Seurat)
library(reshape2)
library(ggbeeswarm)
library(tidyverse)

######## ODC_merged dataset ########
####################################
ODC_merged <- load(file = "./data/ODC_subclustering_sct_v2_clean.rdata")
## CytoTRACE result
mat_ODC <- as.matrix(ODC_merged@assays[["RNA"]]@counts) 
results_ODC <- CytoTRACE(mat = mat_ODC)
cytotrace_score <- as.data.frame(results_ODC["CytoTRACE"])
head(cytotrace_score)
cytotrace_score$Cell <- rownames(cytotrace_score)
plotCytoGenes(results_ODC, numOfGenes = 10, outputDir = "./data_output")
plotCytoTRACE(results_ODC)
save(results_ODC,
     file = "./data_output/ODC_subclustering_sct_v2_clean_cytotrace.rdata")

## umap result
ODC_umap.result <- as.data.frame(ODC_merged@reductions[["umap"]]@cell.embeddings)
ODC_umap.result$Cell <- rownames(ODC_umap.result)
ODC_umap.result <- merge(ODC_umap.result, cytotrace_score, by = "Cell")
ggplot(data = ODC_umap.result, mapping = aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE)) +
  geom_point(size = 0.3) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  scale_color_viridis_c("CytoTRACE",
                        option = "plasma") +
  theme_bw() +
  theme_classic()

ggsave(filename = "ODC_merged_umap.tiff",
       height = 3,
       width = 4)

save(ODC_umap.result,
     file = "./data_output/ODC_cytotrace_umap.rdata")


# prepare the cell annotation file

ODC_cluster <- as.data.frame(ODC_merged@meta.data[["seurat_clusters"]])
ODC_Cell <- colnames(mat_ODC)
ODC_cluster$Cell <- ODC_Cell
colnames(ODC_cluster)[1] <- "Cluster"
write.csv(ODC_cluster,
          file = "./data/ODC_cellAnnotation.csv",
          quote = F,
          row.names = F,
          col.names = T
)

# prepare the umap file
ODC_umap.result <- as.data.frame(ODC_merged@reductions[["umap"]]@cell.embeddings)
ODC_umap.result$Cell <- rownames(ODC_umap.result)
ODC_umap <- data.frame(Cell = ODC_umap.result$Cell,
                         UMAP1 = ODC_umap.result$UMAP_1,
                         UMAP2 = ODC_umap.result$UMAP_2)
write.csv(ODC_umap,
          file = "./data/Fibro_umap.csv",
          quote = F,
          row.names = F,
          col.names = T
)










