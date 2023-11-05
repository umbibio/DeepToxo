library(Seurat)
source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


S.O.tg.boothroyd <- readRDS('../../Input/DeepToxo/rds/S.O.tg_RH_boothroyd.rds')


# Differential gene expression (inferred cell cycle phase)

Idents(S.O.tg.boothroyd) <- 'phase'
DefaultAssay(S.O.tg.boothroyd) <- 'RNA'

cell.cycle.markers <- FindAllMarkers(object = S.O.tg.boothroyd, only.pos = T, min.pct = 0)
cell.cycle.markers$GeneID <- gsub('-', '_', cell.cycle.markers$gene)
cell.cycle.markers <- cell.cycle.markers %>% dplyr::filter(avg_log2FC > log2(2) & p_val_adj < 0.01)

TGGT1_ME49 <- read.xlsx('../../Input/toxo_genomics/ID_convert/TGGT1_ME49 Orthologs.xlsx')
cell.cycle.markers <- inner_join(cell.cycle.markers, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

gs <- lapply(unique(cell.cycle.markers$cluster), function(p){
  gsub('_', '-', cell.cycle.markers$TGME49[cell.cycle.markers$cluster == p])
})

names(gs) <- unique(cell.cycle.markers$cluster)

## Read in deep scRNA data
S.O <- readRDS('../../Input/DeepToxo/rds/S.O.RNA.initial.rds')

# This function is a modification of the following code. Kourosh Zarringhalam, 10/21/2023
# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
es.max <- function(S.O, gs, scaled = F){
  scRNAseqData <- S.O[['RNA']]@counts
  
  marker_stat <- sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity <- data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), 
                                                                             to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  
  
  
  # subselect genes only found in data
  names_gs_cp <- names(gs)
  gs <- lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) <- names_gs_cp
  cell_markers_genes_score <- marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(as.matrix(scRNAseqData)))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] <- Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z <- Z[unique(unlist(gs)), ]
  
  # combine scores
  es <- do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z <- Z[gs[[gss_]], j]
      sum_t1 <- (sum(gs_z) / sqrt(length(gs_z)))
      sum_t1
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max_ <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  return(es.max_)
}


es.max_ <- es.max(S.O, gs, scaled = F)
cell.ident <- apply(es.max_, 2, which.max)
cell.ids <- names(cell.ident)
cell.ident <- data.frame(phase = rownames(es.max_)[cell.ident])
rownames(cell.ident) <- cell.ids

S.O <- AddMetaData(S.O, cell.ident)
Idents(S.O) <- 'phase'
DimPlot(S.O, reduction = 'pca', label = T)

saveRDS(S.O, '../../Input/DeepToxo/rds/S.O.RNA.phase.rds')


## Do a marker analysis on dee scRNA-seq with new phases
Idents(S.O) <- 'phase'
DefaultAssay(S.O) <- 'RNA'

cell.cycle.markers <- FindAllMarkers(object = S.O, only.pos = T, min.pct = 0)
cell.cycle.markers$GeneID <- gsub('-', '_', cell.cycle.markers$gene)
cell.cycle.markers.sig <- cell.cycle.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
prod.desc <- read.xlsx('../../Input/toxo_genomics/genes/ProductDescription_ME49_65.xlsx')
cell.cycle.markers.sig <- left_join(cell.cycle.markers.sig, prod.desc, by = 'GeneID')
write.xlsx(cell.cycle.markers.sig, '../../Output/DeepToxo/tables/cell_cycle_markers.xlsx')
ss <- cell.cycle.markers.sig %>% group_by(cluster) %>% summarise(n())
ss
