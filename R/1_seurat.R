library(SingleCellExperiment)
library(readxl)
library(scDblFinder)
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(assertthat)
library(scCustomize)
library(harmony)
library(UCell)

source(file.path("R","functions.R"))

# Load data
## dataset 1
pth <- "/Users/romansankowski/Downloads/dura"
fls1 <- list.files(pth)
obj1 <- map(fls1, function(x) {
 tryCatch({
    .obj <- Read10X_h5(file.path(pth,x,"outs","filtered_feature_bc_matrix.h5"))
    .obj_seurat <- CreateSeuratObject(.obj, min.cells = 1, min.features = 500) %>% 
      AddMetaData(x, "sample") %>% 
      AddMetaData(PercentageFeatureSet(., pattern = "^MT-"),"percent.mt") %>% 
      subset(subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
    
    ## exclude doublets
    sce <- SingleCellExperiment(assays=list(counts=.obj_seurat@assays$RNA$counts))
    sce <- scDblFinder(sce)
    .obj_seurat <- .obj_seurat[,sce@colData$scDblFinder.class == "singlet"]
    
    .obj_seurat
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
})

## merge lists
all <- Merge_Seurat_List(obj1, add.cell.ids = 1:10)

cond <- case_when(
  all$sample %in% c("GL-1","GL-7", "GL-4", "GL-9") ~ "Naive",
  all$sample %in% c("GL-2","GL-5") ~ "Onset",
  all$sample %in% c("GL-3", "GL-6") ~ "Peak",
  T ~ "Chronic"
) 

names(cond) <- colnames(all)

comp <- case_when(
    all$sample %in% c("GL-1","GL-2","GL-3","GL-7","GL-8") ~ "Sinus",
    T ~ "Extrasinusuidal"
  )
names(comp) <- colnames(all)

all$condition <- cond
all$compartment <- comp

## normalize
all <- all %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = "sample") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 1.5)

## add cell module expression
## add gene modules
lst_signatures <- read_excel(file.path("data","cell_signatures.xlsx"), sheet = 1, skip = 1) %>% 
  as.list() %>% 
  map(na.omit)


## set the module scores
## join layers 
all <- JoinLayers(all)

## add module scores
homing_scores <- ScoreSignatures_UCell(all[["RNA"]]$counts, features=lst_signatures)
all <- AddMetaData(all, as.data.frame(homing_scores))

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all <- CellCycleScoring(all, s.features = s.genes, g2m.features = g2m.genes)

save(all, file=file.path("results", "harmony_integrated_all.RData"))

## normalize using SCTransform
all <- all %>% 
  SCTransform() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 1)

save(all, file=file.path("results", "harmony_integrated_all_sctransform.RData"))

## extract naive cells
naive <- subset(all, cells= colnames(all)[all$condition=="Naive"])

naive <- naive %>% 
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 5000) %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunHarmony(group.by.vars = "sample") %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters()

save(naive, file=file.path("results", "harmony_integrated_naive.RData"))

## normalize using SCTransform
naive <- naive %>% 
  SCTransform() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 1)

save(naive, file=file.path("results", "harmony_integrated_naive_sctransform.RData"))
