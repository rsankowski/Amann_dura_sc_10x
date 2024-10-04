library(Seurat)
#library(SeuratDisk)
library(tidyverse)
library(assertthat)
library(tidyquant)
library(viridis)
library(clusterProfiler)

source(file.path("R","functions.R"))

load(file.path("results", "harmony_integrated_naive_sctransform.RData"))

## reorder clusters
order_clusters <- data.frame(seurat_clusters= naive$seurat_clusters, row.names = rownames(naive[[]])) %>%
  bind_cols(as.data.frame(t(naive[["RNA"]]$scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(naive) <- rev(order_clusters)
naive$seurat_clusters <- factor(naive$seurat_clusters, levels = levels(naive))

## define colors
cols <- colors_many[-c(14,21)]
names(cols) <- levels(naive)

## plot clusters
DimPlot(naive, label = T) +
  theme_void() +
  NoLegend() +
  scale_color_manual(values = cols)

ggsave(file.path("plots","umaps","umap_cluster_naive.pdf"))

## compartments
DimPlot(naive, group.by = "compartment") +
  theme_void() +
  #NoLegend() +
  scale_color_tq() +
  theme(legend.position = "bottom")

ggsave(file.path("plots","umaps","umap_compartment_naive.pdf"))

## marimekko
naive[[]] %>% 
  mosaicGG2(X = "seurat_clusters", FILL = "compartment") +
  scale_fill_tq()

ggsave(file.path("plots","others","marimekko_compartment_naive.pdf"))

## statistical 
naive[[]] %>% 
  hyper_test_n(var1 = "seurat_clusters",var2 = "compartment") %>% 
  write.csv(file.path("results","naive_stat_testing.csv"))

## comparison of macrophages
sum_stat <- naive[[]] %>% 
  filter(seurat_clusters %in% c("10","5","0","19","3","1","4","16")) %>% 
  group_by(compartment, seurat_clusters) %>% 
  summarise(count=n()) %>% 
  group_by(compartment) %>% 
  mutate(ratio=count/sum(count),
         ratio_cum=cumsum(ratio)) 

naive[[]] %>% 
  filter(seurat_clusters %in% c("10","5","0","19","3","1","4","16")) %>% 
  ggplot(aes(x=compartment,  fill=seurat_clusters)) +
  geom_bar(position = "fill", width = 0.5) +
  scale_fill_manual(values = cols[c("10","5","0","19","3","1","4","16")]) +
  theme_linedraw()

ggsave(file.path("plots","others","comparison_compartments_macrophages.pdf"))

naive[[]] %>% 
  filter(seurat_clusters %in% c("10","5","0","19","3","1","4","16")) %>% 
  droplevels() %>% 
  hyper_test_n(var1 = "seurat_clusters",var2 = "compartment") %>% 
  write.csv(file.path("results","naive_stat_testing_macrophages_per_compartment.csv"))

## find cluster markers
if (!file.exists(file.path("results","naive_cluster_markers.csv"))) {
  all_markers <- FindAllMarkers(naive)
  
  save(all_markers, file = file.path("results","naive_cluster_markers.RData"))
  write.csv(all_markers, file.path("results","naive_cluster_markers.csv"))
} else {
  load(file.path("results","naive_cluster_markers.RData"))
}

## heatmap
## plot top 10 markers
      top10 <- all_markers %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
      
      heat <- DoHeatmap(naive,features = top10$gene, group.colors = cols, raster = T)
      heat + scale_fill_viridis(option = "A")
      
      ggsave(file.path("plots","heatmaps","naive-top10-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)
      
## top50 genes
      top50 <- all_markers %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) 
      
      heat <- DoHeatmap(naive,features = top50$gene, group.colors = cols, raster = T)
      heat + scale_fill_viridis(option = "A")
      
      ggsave(file.path("plots","heatmaps","naive-top50-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)
      
      ## export table with top50 genes per cluster
      write.csv(top50, file.path("results","top50_cluster_markers_naive.csv"))

## plot individual genes
genes <- c("Itgax","H2-Ab1","Cd74","Ccr2","Adgre1","Fcgr1","Mrc1","Cx3cr1","Mertk","Aif1","Clec10a","Folr2","Lyve1")

walk(genes, function(x) {
  plt <- FeaturePlot(naive, x) +
    theme_void() +
    scale_color_gradientn(colors = c("darkblue","lightblue2","yellow","red2"))
  print(plt)
  ggsave(file.path("plots","umaps",paste0(x,"-naive.pdf")))
})

## plot signatures
walk(colnames(naive[[]])[10:37], function(x) {
  plt <- FeaturePlot(naive, x) +
    theme_void() +
    scale_color_gradientn(colors = c("darkblue","lightblue2","yellow","red2"))
  print(plt)
  ggsave(file.path("plots","umaps",paste0(x,"-naive.pdf")))
})


genes2 <- c("Clec10a","Mrc1","Folr2","Mgl2","Lyve1","Cd74","H2-Aa","H2-Eb1","H2-Ab1","Ccr2","Itgax")

walk(genes2, function(x) {
  plt <- VlnPlot(naive, x) +
    scale_fill_manual(values = cols)
  
  print(plt)
  ggsave(file.path("plots","others",paste0(x,"-violin-naive.pdf")))
})

## go term analysis:
naive_sub <- subset(naive, cells = colnames(naive)[naive$seurat_clusters %in% c("0", "1", "3", "4", "5", "10")])

naive_sub$seurat_clusters <- droplevels(naive_sub$seurat_clusters)

## find markers
markers_sub <- FindAllMarkers(naive_sub)

top50 <- markers_sub %>% 
  filter(avg_log2FC>0 & p_val_adj <.05) %>% 
  group_by(cluster) %>% 
  top_n(200, wt=avg_log2FC)

## gsea cluster profiler
## convert gene names
genes <- bitr(unique(top50$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = 'org.Mm.eg.db')
colnames(genes)[1] <- "gene"

## comparison between clusters
top50 <- top50 %>% 
  left_join(genes[!duplicated(genes$gene),]) %>% 
  na.omit()

## run diff go terms
diff_go_bp <- compareCluster(ENTREZID ~ cluster,
                             data=top50, 
                             fun = enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             ont = "BP")

clusterProfiler::dotplot(diff_go_bp)
ggsave(file.path("plots","others","go_terms_naive_data_200_genes.pdf"))

write.csv(top50, file.path("results","top200_genes_got_terms_naive.csv"))
