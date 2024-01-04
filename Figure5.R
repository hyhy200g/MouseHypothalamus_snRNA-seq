## import package
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.2")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(sctransform)
library(harmony)
library(dplyr)
library(plyr)
library(viridis)
library(EnhancedVolcano)
library(scales)
library(monocle3)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

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
        cols = c("Astrocytes-2"="#4e659b",
                 "Astrocytes-1"="#fdcf9e"),
        pt.size = 1) +
  theme_classic() +
  theme_void()

# figure5B
###### runmonocle3 ######
#### prepare the input for monocle3 #### 
expression_matrix <- Astrocyte@assays$RNA@counts
cell_metadata <- Astrocyte@meta.data
rownames(cell_metadata)
gene_annotation <- as.data.frame(rownames(Astrocyte@assays$RNA@counts))
colnames(gene_annotation) <- "gene_short_name"
rownames(gene_annotation) <- gene_annotation$gene_short_name

####  Make the CDS object from monocle3 #### 
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 100,method = "PCA")
cds <- align_cds(cds)
plot_pc_variance_explained(cds)

####  dimensional reduction in monocle3  #### 
cds <- reduce_dimension(cds,reduction_method = "UMAP",umap.min_dist = 0.3,preprocess_method = "Aligned")
cds <- cluster_cells(cds,reduction_method = "UMAP",k = 30,cluster_method = "louvain")
plot_cells(cds,graph_label_size=3,cell_size=1.5,show_trajectory_graph = F)

#### change the embeddings to UMAP embedddings from seurat   #### 
newcds<- cds

umap.results <- Astrocyte@reductions$umap@cell.embeddings
head(umap.results)
umap <- umap.results[,c("UMAP_1","UMAP_2")]
umap
umap <- as.data.frame(umap)
umap

umap2 <- newcds@int_colData$reducedDims$UMAP
umap2
umap <- umap[match(rownames(umap2),rownames(umap)),]
head(umap)
newcds@int_colData$reducedDims$UMAP <- umap

####  check the UMAP embeddings by plotting a gene, now the cluster was not changed yet  #### 
plot_cells(newcds, genes=c("Gfap"),show_trajectory_graph = F)
plot_cells(newcds, genes=c("Slc7a10"),show_trajectory_graph = F)

####  change to the seurat clusters  #### 
orig.clusters <- Astrocyte@meta.data$Celltype
orig.clusters
names(orig.clusters) <- rownames(Astrocyte@meta.data)
orig.clusters
newcds@clusters$UMAP$clusters <- orig.clusters

plot_cells(newcds,graph_label_size=3,cell_size=1.5,show_trajectory_graph = F)

#### check level of the partitions  #### 
newcds@clusters$UMAP$partitions
plot_cells(newcds,
           color_cells_by = "partition")

#### trajectory analysis #### 
newcds <- learn_graph(newcds)

plot_cells(newcds,
           labels_per_group = F,
           label_groups_by_cluster=T,
           label_roots = F,
           label_leaves=F,
           label_branch_points=F,
           group_label_size=4,
           cell_size=1.5)

get_earliest_principal_node <- function(newcds,cluster) {
  cell_ids <- which(colData(cds)[, "Celltype"] == cluster)
  closest_vertex <- newcds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(newcds), ])
  root_pr_nodes <- igraph::V(principal_graph(newcds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,
  ]))))]
  root_pr_nodes
}

colData(cds)
root.nodes <- get_earliest_principal_node(newcds,cluster = "Astrocytes-1")
newcds <- order_cells(newcds,root_pr_nodes = root.nodes)


p <- plot_cells(newcds,
                color_cells_by = "pseudotime",
                show_trajectory_graph = F,
                label_roots = T,
                label_cell_groups=F,
                label_leaves=F,
                label_branch_points=F,
                graph_label_size=1.5,
                group_label_size=4,
                cell_size=0.8)
p

plot_cells(newcds,
           group_label_size = 5,
           color_cells_by = "Celltype",
           show_trajectory_graph = T,
           label_roots = T,
           label_cell_groups = T,
           label_groups_by_cluster = T,
           label_leaves = F,
           label_branch_points = F)

#### module analyses ####
set.seed(1)
pr_graph_test_res <- graph_test(newcds, neighbor_graph="principal_graph",cores=8, verbose = FALSE) # Group relative to pseudotime
# save(pr_graph_test_res,file = "Astrocytes_pr_graph_test_res.rdata")

pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
length(pr_deg_ids)
set.seed(2)
gene_module_df <- find_gene_modules(newcds[pr_deg_ids,], resolution=1e-2)
# save(gene_module_df,file = "Astrocytes_gene_module_df.rdata")

cell_group_df <- tibble::tibble(cell=row.names(colData(newcds)), 
                                cell_group=clusters(newcds)[colnames(newcds)])
agg_mat <- aggregate_gene_expression(newcds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("cluster", colnames(agg_mat))


p <- pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                        scale="column", clustering_method="ward.D2",
                        fontsize=6)
p

#########################
matrix <- as.matrix(newcds@assays@data$counts)
used.genes.pr.res <- pr_graph_test_res[rownames(pr_graph_test_res),]
dim(used.genes.pr.res)
summary(used.genes.pr.res)
used.genes.pr.res <- used.genes.pr.res[order(used.genes.pr.res$morans_I,decreasing = T),]
used.genes.pr.res

used.pr.res <- used.genes.pr.res[1:50,]
used.pr.res


genes <- rownames(used.pr.res)
length(genes)
genes
genes <- genes[!genes %in% c("Mt1","Mt2","mt-Cytb")]
genes


pt.matrix <- as.matrix(matrix[match(genes,rowData(newcds)[,1]),order(pseudotime(newcds))])
cells.id <- rownames(astrocyte_merged@meta.data)[astrocyte_merged@meta.data$Celltype %in% c("Astrocytes-1","Astrocytes-2")]
pt.matrix <- pt.matrix[,colnames(pt.matrix) %in% cells.id]
dim(pt.matrix)
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
dim(pt.matrix)


rownames(pt.matrix)
genes.order <- c("Aqp4","Pla2g3","Neat1","Inpp5j","Slc36a2",
                 "Ablim2","Myoc","Padi2","Miat","Gfap",
                 "C4b","Tagln3","Atoh8","Il6st","Slc38a1",
                 "Id1","Aspg","Prune2","Meg3","Nav2",
                 "Id4","Slc38a2","Id3",
                 
                 "Hes5",
                 "Slc1a3","Grk3","Gm3764","Gria2","Nwd1",
                 "Ptprz1","Slc7a10","Plcb1","Etnppl","Abhd14b",
                 "A330076C08Rik","Slc1a2","Atp1a2","Gpr37l1","Slc6a11",
                 "Pla2g7","Plpp3","Htra1","Cst3","Sparcl1","Apoe","Cd81","Malat1")

pt.matrix <- pt.matrix[match(genes.order, rownames(pt.matrix)), ]
dim(pt.matrix)
set.seed(1)
ht3 <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "RdYlBu"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  # km = 3,
  row_title_rot                = 0,
  cluster_rows                 = F,
  cluster_row_slices           = FALSE,
  cluster_columns              = F,
  show_row_dend = T)
ht3


pdf(file = "./Figure5B_modules.pdf",width = 6,height = 7)
ht3
dev.off()


# figure5C
FeaturePlot(ODC, features = "Pdgfra", cols = c("lightgrey", "#8B7D6B")) + theme_void()
FeaturePlot(ODC, features = "Plp1", cols = c("lightgrey", "#8B7D6B")) + theme_void()

# figure5D
## UMAP
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


##Pseudotime analysis using monocle3 
#### prepare the input for monocle3 ###
expression_matrix <- ODC@assays$RNA@counts
cell_metadata <- ODC@meta.data
rownames(cell_metadata)
gene_annotation <- as.data.frame(rownames(ODC@assays$RNA@counts))
colnames(gene_annotation) <- "gene_short_name"
rownames(gene_annotation) <- gene_annotation$gene_short_name


# Make the CDS object from monocle3
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds,num_dim = 100,method = "PCA")
cds <- align_cds(cds)
plot_pc_variance_explained(cds)
str(cds)



######## dimensional reduction #####
cds <- reduce_dimension(cds,reduction_method = "UMAP",umap.min_dist = 0.3,preprocess_method = "Aligned")
cds <- cluster_cells(cds,reduction_method = "UMAP",k = 30,cluster_method = "louvain")
plot_cells(cds,graph_label_size=3,cell_size=1.5,show_trajectory_graph = F)


newcds<- cds

### umap embeddings for Seurat
umap.results <- ODC@reductions$umap@cell.embeddings
head(umap.results)
umap <- umap.results[,c("UMAP_1","UMAP_2")]
umap
umap <- as.data.frame(umap)
umap

### umap embeddings from Monocle3
umap2 <- newcds@int_colData$reducedDims$UMAP
umap2
umap <- umap[match(rownames(umap2),rownames(umap)),]
head(umap)
newcds@int_colData$reducedDims$UMAP <- umap

####################### change the clusters to seurat clusters #####################
orig.clusters <- ODC@meta.data$Celltype
orig.clusters
names(orig.clusters) <- rownames(ODC@meta.data)
orig.clusters
newcds@clusters$UMAP$clusters <- orig.clusters

### check the clustering 
plot_cells(newcds,graph_label_size=3,cell_size=1.5,show_trajectory_graph = F)


### check level of the partitions 
newcds@clusters$UMAP$partitions
plot_cells(newcds,
           color_cells_by = "partition")


####################### trajectory analysis #####################
newcds <- learn_graph(newcds)

plot_cells(newcds,
           labels_per_group = F,
           label_groups_by_cluster=T,
           label_roots = F,
           label_leaves=F,
           label_branch_points=F,
           group_label_size=4,
           cell_size=1.5)


###################### identify the root for pseudotime analysis  - we used cluster 1 here #####################
get_earliest_principal_node <- function(newcds,cluster) {
  cell_ids <- which(colData(cds)[, "Celltype"] == cluster)
  closest_vertex <- newcds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(newcds), ])
  root_pr_nodes <- igraph::V(principal_graph(newcds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,
  ]))))]
  root_pr_nodes
}

colData(cds)
root.nodes <- get_earliest_principal_node(newcds,cluster = "ODC-1")
newcds <- order_cells(newcds,root_pr_nodes = root.nodes)

p <- plot_cells(newcds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = F,
                  label_roots = T,
                  label_cell_groups=F,
                  label_leaves=F,
                  label_branch_points=F,
                  graph_label_size=1.5,
                  group_label_size=4,
                  cell_size=0.5)
p


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
  scale_fill_manual(values = c("#E89242FF",
                      "#fdcf9e",
                      "lightgreen",
                      "#4f8c9d",
                      "#4e659b"))
