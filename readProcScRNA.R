library(tidyverse)
library(Seurat)
library(openxlsx)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


# Count files
data.dir <- "../../Input/DeepToxo/scRNAseq/RH-RNA/outs/filtered_feature_bc_matrix/"

scRNA <- Read10X(data.dir = data.dir, gene.column = 1) ## column = 1 uses gene id; 2 uses gene name


## IDs
prod.desc  <- read.xlsx('../../Input/toxo_genomics/genes/ProductDescription_ME49_65.xlsx')

feats <- c("nFeature_RNA","nCount_RNA")

S.O <- CreateSeuratObject(counts = scRNA)
VlnPlot(S.O, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O, expression = nFeature_RNA > 150 & nFeature_RNA < 1700)
selected_f <- rownames(S.O)[ Matrix::rowSums(S.O) > 5]
S.O  <- subset(S.O, features=selected_f, cells=selected_c)
dim(S.O@assays$RNA@data)


saveRDS(S.O, '../../Input/DeepToxo/rds/S.O.RNA.initial.rds')

#set.seed(100)
#S.O <- subset(x = S.O.intra, downsample = 8000)
S.O <- prep_S.O(S.O, res = 0.4)
DimPlot(S.O, reduction = 'pca')


# Transfer cell cycle phase labels from available data https://elifesciences.org/articles/54129

S.O.list <- list(intra = S.O.intra)

S.O.tg.boothroyd <- readRDS('../Input_sub/toxo_cdc/rds_ME49_59/S.O.tg_RH_boothroyd.rds')

S.Os <- mclapply(S.O.list, function(S.O){
  S.O <- prep_S.O(S.O, res = 0.4)
  anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

spps <- names(S.Os)

S.Os <- lapply(1:length(S.Os), function(i){
  S.Os[[i]]@meta.data$spp <- spps[i]
  S.Os[[i]]
})




