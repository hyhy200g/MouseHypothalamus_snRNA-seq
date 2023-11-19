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

## FigureS2
# figureS2A
UMIs <- data.frame(UMIs = Hypo_integrated@meta.data$nCount_RNA)
mean(UMIs$UMIs)
ggplot(UMIs, aes(UMIs)) +
  geom_histogram(bins = 225, color = "grey", alpha = 0.6) +
  geom_vline(xintercept = mean(UMIs$UMIs), linetype = "dashed", color = "red") +
  theme_classic()

Genes <- data.frame(Genes = Hypo_integrated@meta.data$nFeature_RNA)
ggplot(Genes, aes(x= Genes, fill = count)) +
  geom_histogram(bins = 75, color = "black", alpha = 1, position = "identity", fill = "#4876FF") +
  geom_vline(xintercept = mean(Genes$Genes), linetype = "dashed", color = "red") +
  theme_classic()

mito <- data.frame(mito = Hypo_integrated@meta.data$percent.mt)
ggplot(mito, aes(x= mito, fill = count)) +
  geom_histogram(bins = 60, color = "black", alpha = 1, position = "identity", fill = "darkgreen") +
  geom_vline(xintercept = mean(mito$mito), linetype = "dashed", color = "red") +
  theme_classic() 

# figureS2B
DimPlot(Hypo_integrated, group.by = "orig.ident")
DimPlot(Hypo_integrated, group.by = "batch")
DimPlot(Hypo_integrated, group.by = "sex")
DimPlot(Hypo_integrated, group.by = "group")

# figureS2C
Hypo_sample <- data.frame(sample = Hypo_integrated@meta.data$orig.ident,
                          cluster = Hypo_integrated@meta.data$Celltype)
Hypo_sample$sample <- factor(Hypo_sample$sample, levels = c("m48f", "m6f", "m39f", "m47f",
                                                            "m50m", "m71m", "m36m", "m44m"))
ggplot(Hypo_sample, aes(sample, fill=cluster)) +
  geom_bar(position = "fill", width = 0.9) +
  theme_classic() +
  labs(y = "Cell proportion") +
  scale_fill_manual(values = c('Astrocytes'="#526E2DFF",
                               'Endothelial cells'="#f6bb86", 
                               'Interneurons'="#82491EFF",
                               'B cells'="#69C8ECFF",
                               'Microglial/Macrophages'="#24325FFF",
                               'Neurons'="#E89242FF", 
                               'Oligodendrocytes'="#917C5DFF", 
                               'Stromal cells'="#9698dc",
                               'Tanycytes'= "#4f8c9d"))


# figureS2E
#### import packages ####
library(symphony)
library(harmony)
library(singlecellmethods)
library(irlba)
library(tidyverse)
library(data.table)
library(matrixStats)
library(Matrix)
library(plyr)
library(dplyr)
library(rstan)
library(Seurat)
library(ggplot2)
library(ggthemes)
library(ggrastr)
library(RColorBrewer)
library(patchwork)
library(ggpubr)
library(uwot)
source("~/R/project/script/buildReference_on_UMAP.R")
source("~/R/project/script/symphony_utils.R")

## Load data ##
readRDS("~/R/project/data/hypoMap.rds")
readRDS("~/R/project/data/Hypo_integrated.rds")

## Build reference ##
ref_exp_full <- hypoMap@assays[["RNA"]]@counts
ref_metadata <- data.frame(Cell = hypoMap@meta.data[["Cell_ID"]],
                           Cell_type = hypoMap@meta.data[["Author_Class_Curated"]],
                           Dataset = hypoMap@meta.data[["Dataset"]],
                           Technology = hypoMap@meta.data$Technology)

sc <- as.matrix(hypoMap@reductions$umap_scvi@cell.embeddings)
colnames(sc) <- c("UMAP1", "UMAP2")
set.seed(0)
reference_umap = buildReference_on_UMAP(
  ref_exp_full,
  ref_metadata,
  vars = NULL,               # variables to integrate over
  sc2d = sc,
  theta = c(0.8),
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = TRUE,            # can set to FALSE if want to run umap separately later
  do_normalize = T,          # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = NULL,    # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20,                    # number of PCs
  save_uwot_path = '~/R/project/data/testing_uwot_model_on_umap' 
)




## Map query ##
query_exp <- Hypo_integrated@assays[["RNA"]]@counts
Cluster <- Hypo_integrated@meta.data[["seurat_clusters"]]
Cell <- Hypo_integrated@meta.data[["Cell"]]
Cell_type <- Hypo_integrated@meta.data[["Celltype"]]
query_metadata <- data.frame(Cell = Cell,
                             Cluster = Cluster,
                             Cell_type = Cell_type)

query_umap <- mapQuery(query_exp,
                       query_metadata,
                       reference_umap,
                       do_normalize = TRUE,
                       do_umap = TRUE)

query_umap <- knnPredict(query_umap, 
                         reference_umap,
                         reference_umap$meta_data$Cell_type,
                         k = 5,
                         confidence = TRUE)


query_labels_u <- cbind(query_metadata, query_umap$umap)
umap_labels_ref = cbind(reference_umap$meta_data, reference_umap$umap$embedding)
cell_type_col <- c('Astrocytes'="#526E2DFF",
                   'Endothelial cells'="#f6bb86", 
                   'Interneurons'="#82491EFF",
                   'B cells'="#69C8ECFF",
                   'Microglial/Macrophages'="#24325FFF",
                   'Neurons'="#E89242FF", 
                   'Oligodendrocytes'="#917C5DFF", 
                   'Stroma cells'="#9698dc",
                   'Tanycytes'= "#4f8c9d")

query_labels_u$Cell_type <- factor(query_labels_u$Cell_type,
                                   levels = c("Astrocytes", "B cells", "Endothelial cells", "Interneurons", 
                                              "Microglial/Macrophages", "Oligodendrocytes", "Tanycytes", "Neurons"))
options(repr.plot.height = 3,
        repr.plot.width = 4)

umap_labels_ref <- umap_labels_ref %>%
  sample_frac(1L)

query_labels_u <- query_labels_u %>%
  sample_frac(1L)

## the plot of query annotated with cell type ##
p1 <- ggplot(data = query_labels_u, mapping = aes(x = UMAP1, y = UMAP2, color = Cell_type)) +
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
                                'Stroma cells'="#9698dc",
                                'Tanycytes'= "#4f8c9d"))


## the plot of query projection to the reference ##
p2 <- ggplot(NULL) +
  geom_point(data = umap_labels_ref, mapping = aes(x = UMAP1, y = UMAP2), col = "#CCCCCC", size = 0.5) +
  geom_point(data = query_labels_u, mapping = aes(x = UMAP1, y = UMAP2), col = "#87CEFA", size = 0.5) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  theme_bw() +
  theme(legend.position = "none") +
  theme_void()

## the plot of reference annotated with cell type ##

options(repr.plot.width = 4, repr.plot.height = 3)
p3 <- ggplot(data = umap_labels_ref, mapping = aes(x = UMAP1, y = UMAP2, color = Cell_type)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  theme_bw() +
  theme_classic()


## the plot of query projection to the reference annotated with cell type ##
p4 <- ggplot(NULL) +
  geom_point(data = umap_labels_ref, mapping = aes(x = UMAP1, y = UMAP2), col = "#C0C0C0", size = 0.5) +
  geom_point(data = query_labels_u, mapping = aes(x = UMAP1, y = UMAP2, col = Cell_type), size = 0.5) +
  scale_color_manual(values = c('Astrocytes'="#526E2DFF",
                                'Endothelial cells'="#f6bb86", 
                                'Interneurons'="#82491EFF",
                                'B cells'="#69C8ECFF",
                                'Microglial/Macrophages'="#24325FFF",
                                'Neurons'="#E89242FF", 
                                'Oligodendrocytes'="#917C5DFF", 
                                'Stroma cells'="#9698dc",
                                'Tanycytes'= "#4f8c9d")) +
  labs(x = "UMAP1",
       y = "UMAP2") +
  theme_bw() +
  theme_void()  +
  theme(legend.position = "none")

## Use Dimplot() to draw a plot of reference ##
pd <- DimPlot(hypoMap,
              reduction = "umap_scvi",
              label = TRUE,
              group.by = "Author_Class_Curated",
              label.size = 4,
              repel = TRUE)
pd + theme_void() + theme(legend.position = "none") + ggtitle("")


pd.nolable <- DimPlot(hypoMap,
                      reduction = "umap_scvi",
                      label = FALSE,
                      group.by = "Author_Class_Curated")
pd.nolable + theme_void() + theme(legend.position = "none")