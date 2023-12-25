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
# figureS2A---------------------------------
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

# figureS2B---------------------------------
sample.color <- c('m36m'="#D97176FF",'m39f'="#917C5DFF",'m44m'='#526E2DFF', 'm47f'='#82491EFF', 'm48f'='#4f8c9d', 'm50m'='#9698dc', 'm6f'='#24325FFF', 'm71m'='#f6bb86')
batch.color <- c('batch1'="#D97176FF",'batch2'="#917C5DFF",'batch3'="#526E2DFF",'batch4'="#f6bb86") 
sex.color <- c('female'="#f6bb86",'male'="#526E2DFF") 
group.color <- c('male_LFD'="#E89242FF",'male_HFD'="#917C5DFF",'female_LFD'='#526E2DFF', 'female_HFD'='#82491EFF')

DimPlot(Hypo_integrated, group.by = "orig.ident", cols = sample.color)
DimPlot(Hypo_integrated, group.by = "batch", cols = batch.color)
DimPlot(Hypo_integrated, group.by = "sex", cols = sex.color)
DimPlot(Hypo_integrated, group.by = "group", cols = group.color)

# figureS2C---------------------------------
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

# figureS2D---------------------------------

# figureS2E---------------------------------
##Calculate the DEGs between maternal dietary groups in major clusters
mouse.combined.cca.clean <- readRDS(file = "mouse.combined.cca.clean_recalculate_sct.rds")
mouse.combined.cca.clean$new_id <- paste(mouse.combined.cca.clean$Celltype,"_",mouse.combined.cca.clean$sex,"_",mouse.combined.cca.clean$diet,sep = "")
clusters <- names(table(mouse.combined.cca.clean$Celltype))
mouse.combined.cca.clean <- SetIdent(mouse.combined.cca.clean,
                                     value = mouse.combined.cca.clean@meta.data$new_id)

Idents(mouse.combined.cca.clean)
DefaultAssay(mouse.combined.cca.clean) <- "SCT"

#mouse.combined.cca.clean <- PrepSCTFindMarkers(mouse.combined.cca.clean, verbose = T)
#Minimum UMI unchanged. Skipping re-correction.


for (i in clusters){
  DE_LFM_vs_HFM_GSEA <- FindMarkers(mouse.combined.cca.clean, assay = "SCT", ident.1 = as.character(paste(i,"_male_HFD",sep = "")), ident.2  = as.character(paste(i,"_male_LFD",sep = "")),
                                   logfc.threshold = 0, min.pct = 0, min.cells.feature = 1,test.use = "MAST")
  write.csv(DE_LFM_vs_HFM_GSEA, file =paste("./final_DEG_SCTV2/Major/",i,"_LFM_vs_HFM_DEG_GSEA_SCT_MAST_20221009.csv",sep =""))
  DE_LFM_vs_HFM <- DE_LFM_vs_HFM_GSEA %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC) 
  write.csv(DE_LFM_vs_HFM, file =paste("./final_DEG_SCTV2/Major/",i,"_LFM_vs_HFM_DEG_SCT_MAST_20221009.csv",sep =""))
}


for (i in clusters){
  DE_LFF_vs_HFF_GSEA<- FindMarkers(mouse.combined.cca.clean, assay = "SCT", ident.1 = as.character(paste(i,"_female_HFD",sep = "")), ident.2  = as.character(paste(i,"_female_LFD",sep = "")),
                                   logfc.threshold = 0, min.pct = 0, min.cells.feature = 1,test.use = "MAST")
  write.csv(DE_LFF_vs_HFF_GSEA, file =paste("./final_DEG_SCTV2/Major/",i,"_LFF_vs_HFF_DEG_GSEA_SCT_MAST_20221009.csv",sep =""))
  DE_LFF_vs_HFF <- DE_LFF_vs_HFF_GSEA %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC) 
  write.csv(DE_LFF_vs_HFF, file =paste("./final_DEG_SCTV2/Major/",i,"_LFF_vs_HFF_DEG_SCT_MAST_20221009.csv",sep =""))
}

for (i in clusters){
  DE_HFF_vs_HFM_GSEA<- FindMarkers(mouse.combined.cca.clean, assay = "SCT", ident.1 = as.character(paste(i,"_male_HFD",sep = "")), ident.2  = as.character(paste(i,"_female_HFD",sep = "")),
                                   logfc.threshold = 0, min.pct = 0, min.cells.feature = 1,test.use = "MAST")
  write.csv(DE_HFF_vs_HFM_GSEA, file =paste("./final_DEG_SCTV2/Major/",i,"_HFF_vs_HFM_DEG_GSEA_SCT_MAST_20221009.csv",sep =""))
  DE_HFF_vs_HFM <- DE_HFF_vs_HFM_GSEA %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC) 
  write.csv(DE_HFF_vs_HFM, file =paste("./final_DEG_SCTV2/Major/",i,"_HFF_vs_HFM_DEG_SCT_MAST_20221009.csv",sep =""))
}

for (i in clusters){
  DE_LFF_vs_LFM_GSEA<- FindMarkers(mouse.combined.cca.clean, assay = "SCT", ident.1 = as.character(paste(i,"_male_LFD",sep = "")), ident.2  = as.character(paste(i,"_female_LFD",sep = "")),
                                   logfc.threshold = 0, min.pct = 0, min.cells.feature = 1,test.use = "MAST")
  write.csv(DE_LFF_vs_LFM_GSEA, file =paste("./final_DEG_SCTV2/Major/",i,"_LFF_vs_LFM_DEG_GSEA_SCT_MAST_20221009.csv",sep =""))
  DE_LFF_vs_LFM <- DE_LFF_vs_LFM_GSEA %>% filter(p_val_adj<0.05) %>% arrange(avg_log2FC) 
  write.csv(DE_LFF_vs_LFM, file =paste("./final_DEG_SCTV2/Major/",i,"_LFF_vs_LFM_DEG_SCT_MAST_20221009.csv",sep =""))
}


##create a function to select the top and bottom n genes
major_DEG <- read.csv(file = "./final_DEG_SCTV2/Major_cluster_DEG_results_organization_20231123.csv", header = T)
major_DEG$celltype <-substr(major_DEG$Comparisons, 1, nchar(major_DEG$Comparisons) - 11)
major_DEG$Comparisons <- as.character(major_DEG$Comparisons)
major_DEG$groups <- substr(major_DEG$Comparisons, nchar(major_DEG$Comparisons) - 9, nchar(major_DEG$Comparisons))

select_top_bottom_genes <- function(df, n = 10) {
  # List to store results
  results_list <- list()

  # Unique comparisons
  comparisons <- unique(df$groups)

  # Loop through each comparison
  for (comp in comparisons) {
    # Subset DataFrame for the current comparison
    subset_df <- df[df$groups == comp, ]

    # Select top and bottom genes based on LFC
    top_genes <- head(subset_df[order(subset_df$avg_log2FC, decreasing = TRUE), ], n)
    bottom_genes <- head(subset_df[order(subset_df$avg_log2FC), ], n)

    results_list[[comp]] <- rbind(top_genes,bottom_genes)
  }
  
  out <- bind_rows(results_list)
  
  return(out) 
}

##create a function for the dot plot
dot_plot_DEG <- function(celltype){
DEG <- major_DEG[which(major_DEG$celltype %in% c(celltype)),]
DEG <- select_top_bottom_genes(DEG, n = 5) #Apply the select_top_bottom_genes function 
genes <- DEG$DEGs
comp <- names(table(DEG$groups))
DEG_list <- list()
for (i in comp) {
  GSEA_DEG <- read.csv(file = paste0("/Volumes/T7_1/Project_snRNA-seq_seurat/final_DEG_SCTV2/Major/",celltype,"_",i,"_DEG_GSEA_SCT_MAST_20221009.csv")) 
GSEA_DEG <- GSEA_DEG %>% rename(DEGs = X)
GSEA_DEG <- GSEA_DEG[which(GSEA_DEG$DEGs %in% genes),]
GSEA_DEG["celltype"] <- celltype
GSEA_DEG["Comparisons"] <- i
DEG_list[[i]] <- GSEA_DEG
}

#organize results
  
result_plot <- bind_rows(DEG_list)
result_plot$padj_1_plot <- NA
result_plot$padj_1_plot <- as.numeric(ifelse(result_plot$p_val_adj >= 0.05, "0.05", result_plot$p_val_adj))
result_plot$Comparisons <- factor(result_plot$Comparisons, levels = c("fLFD vs fHFD","mLFD vs mHFD","fLFD vs mLFD","fHFD vs mHFD"))

#plot
  
ggplot(data = result_plot, aes(x = Comparisons, y = DEGs, 
                               color = avg_log2FC, size = desc(padj_1_plot))) + 
  geom_point() +
  theme_classic() +
  ylab("") + 
  xlab("") + 
  scale_color_gradientn(
    colours =  c("#4f8c9d","lightgrey","#a40000"),    
    values = scales::rescale(c(-3.5, 0, 1)), 
    limits = c(min(result_plot$avg_log2FC), 1), 
    breaks = c(-3, 0, 1),
    na.value = "gray80"  # Set color for NA values
  ) +
  theme_ipsum() +
  theme(
    axis.text.x = element_text(face = "bold", size = 14,angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold.italic", size = 14)
  ) +
  labs(size = "Padj") +
  scale_y_discrete(limits = rev) +
  guides(color = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.3, "cm")))



ggsave(filename = paste0(celltype,'_top5_DEG.pdf'), units = "in",device = cairo_pdf, 
       path ="/Volumes/T7_1/Project_snRNA-seq_seurat/final_DEG_SCTV2/plots" ,width = 6, height = 8, dpi = 300)

return(result_plot) 
}

dotplot_astrocytes <- dot_plot_DEG(celltype = "Astrocytes")
dotplot_neuron <- dot_plot_DEG(celltype = "Neuron")
dotplot_oligo <- dot_plot_DEG(celltype = "Oligodendrocytes")


# figureS2F---------------------------------
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
