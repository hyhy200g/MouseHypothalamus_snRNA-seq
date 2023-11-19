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

# figure5A
FeaturePlot(Astrocyte, features = "Gfap", cols = c("lightgrey", "darkgreen")) + theme_void()
FeaturePlot(Astrocyte, features = "Slc7a10", cols = c("lightgrey", "darkgreen")) + theme_void()
DimPlot(Astrocyte,
        reduction = "umap",
        cols = c("#fdcf9e",
                 "#4e659b"),
        pt.size = 1) +
  theme_classic() +
  theme_void()

# figure5B
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
devtools::install_github('cole-trapnell-lab/monocle3')

# figure5C
FeaturePlot(ODC, features = "Pdgfra", cols = c("lightgrey", "#8B7D6B")) + theme_void()
FeaturePlot(ODC, features = "Plp1", cols = c("lightgrey", "#8B7D6B")) + theme_void()

# figure5D
ODC@meta.data$Celltype <- factor(ODC@meta.data$Celltype, levels = c("ODC-1", "ODC-2", "ODC-3", "ODC-4", "ODC-5"))
levels(ODC@meta.data$Celltype)
DimPlot(ODC, cols = c("#E89242FF",
                      "#fdcf9e",
                      "lightgreen",
                      "#4f8c9d",
                      "#4e659b"),
        group.by = "Celltype") +
  theme_void() +
  theme(title = element_blank())
ODC_mat <- as.matrix(ODC@assays[["RNA"]]@counts)
ODC_result <- CytoTRACE(mat = ODC_mat)
ODC_score <- as.data.frame(ODC_result["CytoTRACE"]) %>%
  mutate(Cell = rownames(ODC_score))
ODC_umap <- as.data.frame(ODC@reductions$umap@cell.embeddings) %>%
  mutate(Cell = rownames(ODC_umap))
ODC_df <- merge(ODC_umap, ODC_score, by = "Cell")
ggplot(data = ODC_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = CytoTRACE)) +
  geom_point(size = 0.3) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  scale_color_viridis_c("CytoTRACE",
                        option = "plasma") +
  theme_bw() +
  theme_void()

# figure5E
ODC_numi <- ODC[["RNA"]]@counts

normalize <- function(x){
  return(log(x/sum(x)*10000 + 1))
}

ODC_norm <- apply(ODC_numi, 2, normalize)

ODC_expr <- data.frame(cell = colnames(ODC_norm),
                       Pdgfra = ODC_norm["Pdgfra",],
                       Plp1 = ODC_norm["Plp1",],
                       group = ODC@meta.data$Celltype)
ODC_expr$group <- factor(ODC_expr$group)
levels(ODC_expr$group) <- c("ODC-1", "ODC-2", "ODC-3", "ODC-4", "ODC-5")

ggplot(ODC_expr, aes(x=group, y=Plp1, fill=group)) +
  geom_violin(scale = "width") +
  theme_classic() +
  theme(legend.position = "none") +
  labs(y="Normalized UMI") +
  stat_summary(fun = "mean", geom = "point", color="yellow", size=3) +
  scale_fill_manual(values = c("brown", "darkgreen", "grey", "orange", "blue"))

# Figure5D_Cytotrace
library(CytoTRACE)

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