run_DAseq_comparison <- function(object,
                                 groupA = "", 
                                 groupB = "",
                                 labels.1 = "", 
                                 labels.2 = "",
                                 path = "",
                                 outputname = ""){
  

  # extract umap embeddings
  meta.dadta <- object@meta.data
  umap.results <- object@reductions$umap@cell.embeddings
  umap.results <- as.data.frame(umap.results)
  umap.results$Cell <- rownames(umap.results)
  head(umap.results)
  umap.results <- merge(umap.results,meta.data, by = "Cell")
  s.umap <- umap.results[,c("UMAP_1","UMAP_2")]
  dim(s.umap)
  
  umap.results <- subset(umap.results,group == groupA | group == groupB)
  
  # extract labellings
  s.labels <- umap.results$orig.ident
  unique(s.labels)
  
  # extract pca
  s.pca <- object@reductions$pca@cell.embeddings
  s.pca
  s.pca <- s.pca[rownames(s.pca) %in% umap.results$Cell,]
  s.pca <- s.pca[match(umap.results$Cell,rownames(s.pca)),]
  

  ### DA-seq analysis ###
  da_cells <- getDAcells(
    X = s.pca,
    cell.labels = s.labels,
    labels.1 = labels.1,
    labels.2 = labels.2,
    k.vector = seq(20,500,200),
    plot.embedding = s.umap
  )
  
  # da_cells$pred.plot
  da_cells$rand.plot
  
  updated_da_cells <- updateDAcells(
    X = da_cells, pred.thres = c(-0.5,0.5),
    plot.embedding = s.umap
  )

  #################################
  cells.up.idx<-updated_da_cells$da.up
  cells.down.idx<-updated_da_cells$da.down
  #################################
  cells.info.up<-umap.results[cells.up.idx,]
  cells.info.down<-umap.results[cells.down.idx,]



  pdf(paste0(path,outputname))
  plot(s.umap,pch=19,cex=0.2,col="#D9D9D9FF", xaxt="n",yaxt="n",xlab="",ylab="")
  points(s.umap[rownames(s.umap) %in% rownames(cells.info.up),],pch=19,cex=0.2,col=alpha("#1984c5", 0.2))
  points(s.umap[rownames(s.umap) %in% rownames(cells.info.down),],pch=19,cex=0.2,col=alpha("#c23728", 0.2))
  legend("topleft",c(groupA,groupB),pch=19,col=c("#c23728","#1984c5"))
  dev.off()
  
}
