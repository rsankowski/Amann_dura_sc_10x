library(Seurat)
#library(SeuratDisk)
library(tidyverse)
library(assertthat)
library(tidyquant)
library(viridis)
library(clusterProfiler)

source(file.path("R","functions.R"))

load(file.path("results", "harmony_integrated_all_sctransform.RData"))

## reorder clusters
order_clusters <- data.frame(seurat_clusters= all$seurat_clusters, row.names = rownames(all[[]])) %>%
  bind_cols(as.data.frame(t(all[["RNA"]]$scale.data))) %>%
  group_by(seurat_clusters) %>%
  summarize_all(.funs=mean) %>%
  as.data.frame()

rownames(order_clusters) <- order_clusters$seurat_clusters
order_clusters <- order_clusters$seurat_clusters[hclust(dist(order_clusters[,-1]))$order]

#reorder the clusters
levels(all) <- rev(order_clusters)
all$seurat_clusters <- factor(all$seurat_clusters, levels = levels(all))

## plot clusters
DimPlot(all, label = T) +
  theme_void() +
  NoLegend() +
  scale_color_manual(values = c(colors_many, colors_fig))

ggsave(file.path("plots","umaps","umap_cluster_all.pdf"))

## compartments
DimPlot(all, group.by = "compartment") +
  theme_void() +
  #NoLegend() +
  scale_color_tq() +
  theme(legend.position = "bottom")

ggsave(file.path("plots","umaps","umap_compartment_all.pdf"))

## marimekko
all[[]] %>% 
  mosaicGG2(X = "seurat_clusters", FILL = "compartment") +
  scale_fill_tq()

ggsave(file.path("plots","others","marimekko_compartment_all.pdf"))

## statistical 
all[[]] %>% 
  hyper_test_n(var1 = "seurat_clusters",var2 = "compartment") %>% 
  write.csv(file.path("results","all_stat_testing.csv"))

## conditions
DimPlot(all, group.by = "condition") +
  theme_void() +
  #NoLegend() +
  scale_color_tq(palette="dark") +
  theme(legend.position = "bottom")

ggsave(file.path("plots","umaps","umap_coondition_all.pdf"))

## marimekko
all[[]] %>% 
  mosaicGG2(X = "seurat_clusters", FILL = "condition") +
  scale_fill_tq(palette="dark")

ggsave(file.path("plots","others","marimekko_condition_all.pdf"))

## statistical 
all[[]] %>% 
  hyper_test_n(var1 = "seurat_clusters",var2 = "condition") %>% 
  write.csv(file.path("results","all_stat_testing_condition.csv"))

## find cluster markers
if (!file.exists(file.path("results","all_cluster_markers.csv"))) {
  all_markers <- FindAllMarkers(all)
  
  save(all_markers, file = file.path("results","all_cluster_markers.RData"))
  write.csv(all_markers, file.path("results","all_cluster_markers.csv"))
} else {
  load(file.path("results","all_cluster_markers.RData"))
}

## heatmap
## plot top 10 markers
top10 <- all_markers %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top10$gene, group.colors = c(colors_many, colors_fig), raster = T)
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","all-top10-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

## top50 genes
top50 <- all_markers %>% filter(p_val_adj <0.05) %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) 

heat <- DoHeatmap(all,features = top50$gene, group.colors = c(colors_many, colors_fig), raster = T)
heat + scale_fill_viridis(option = "A")

ggsave(file.path("plots","heatmaps","all-top50-gene-heatmap-viridis-A.pdf"), width = 20, height = 16)

## export table with top50 genes per cluster
write.csv(top50, file.path("results","top50_cluster_markers_all.csv"))

## plot individual genes
genes <- c("Itgax","H2-Ab1","Cd74","Ccr2","Adgre1","Fcgr1","Mrc1","Cx3cr1","Mertk","Aif1","Clec10a","Folr2")

walk(genes, function(x) {
  plt <- FeaturePlot(all, x) +
    theme_void() +
    scale_color_gradientn(colors = c("darkblue","lightblue2","yellow","red2"))
  print(plt)
  ggsave(file.path("plots","umaps",paste0(x,"-all.pdf")))
})

## plot signatures
walk(colnames(all[[]])[10:37], function(x) {
  plt <- FeaturePlot(all, x) +
    theme_void() +
    scale_color_gradientn(colors = c("darkblue","lightblue2","yellow","red2"))
  print(plt)
  ggsave(file.path("plots","umaps",paste0(x,"-all.pdf")))
})

## go term analysis:
all_sub <- subset(all, cells = colnames(all)[all$seurat_clusters %in% c("2", "4", "7", "28")])

all_sub$seurat_clusters <- droplevels(all_sub$seurat_clusters)

## find markers
markers_sub <- FindAllMarkers(all_sub)

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

write.csv(top50 ,file.path("results","macrophages_all_top200_marker_genes.csv"))

## run diff go terms
diff_go_bp <- compareCluster(ENTREZID ~ cluster,
                             data=top50, 
                             fun = enrichGO,
                             OrgDb = 'org.Mm.eg.db',
                             ont = "BP")

clusterProfiler::dotplot(diff_go_bp)
ggsave(file.path("plots","others","go_terms_all_data_top200_genes.pdf"))

