# The function 'ccapca' solves the problem 
# max_{u,v} a (u' X'X u + v' Y'Y v)/2 + (1 - a) u' X'Y v
# subject to u'u = v'v = 1 
# If k > 1 components (u,v) are required, the components 
# u1,...,uk are required to be orthogonal 
# and same for the v1, ..., vk 
# (along with the unit norm constraints)

# For a = 0, this is simply the SVD of the cross-covariance between X and Y
# (covariance-based CCA)
# For a = 1, this is simply the PCA of X and of Y
# For 0 < a < 1, this is a compromise between covariance-based CCA and PCA


ccapca <- function(x, y, a = 0, k, center = TRUE)
{
  if (!is.matrix(x)) x <- as.matrix(x)	
  if (!is.matrix(y)) y <- as.matrix(y)	
  stopifnot(nrow(x) == nrow(y))
  stopifnot(a >= 0 && a <= 1)
  k <- as.integer(k)
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  if (center) {
    x <- scale(x, scale = FALSE)
    y <- scale(y, scale = FALSE)
  }
  
  if (a == 0) { 
    k <- min(k, n, p, q)
    if (n <= max(p,q)/2) {
      qrx <- qr(t(x))
      qry <- qr(t(y))
      rr <- tcrossprod(qr.R(qrx), qr.R(qry))
      svdrr <- svd(rr, nu = k, nv = k)
      u <- qr.Q(qrx) %*% svdrr$u
      v <- qr.Q(qry) %*% svdrr$v
    } else {
      svdxy <- svd(crossprod(x, y), nu = k, nv = k)
      u <- svdxy$u
      v <- svdxy$v
    }
  } else if (a == 0.5) {
    k <- min(k, n, p + q)
    svdxy <- svd(cbind(x,y), nu = 0, nv = k)
    u <- svdxy$v[1:p,]
    v <- svdxy$v[-(1:p),]
  } else if (a == 1) { 
    u <- svd(x, nu = 0, nv = min(k, n, p))$v
    v <- svd(y, nu = 0, nv = min(k, n, q))$v
  } else { # general case
    k <- min(k, p + q)
    mat <- matrix(0, p + q, p + q)
    mat[1:p, 1:p] <- (a/2) * crossprod(x)
    mat[1:p, -(1:p)] <- (1 - a) * crossprod(x, y)
    mat[-(1:p), 1:p] <- t(mat[1:p, -(1:p)])
    mat[-(1:p), -(1:p)] <- (a/2) * crossprod(y)
    eigmat <- eigen(mat, TRUE)
    u <- eigmat$vectors[1:p, 1:k]
    v <- eigmat$vectors[-(1:p), 1:k]
  }
  
  list(u = u, v = v)
}

######
library(tidyverse)
library(Seurat)
library(openxlsx)

S.O.rna <- readRDS('../../Input/DeepToxo/rds/S.O.RNA.phase.pt.rds')
S.O.atac <- readRDS('../../Input/DeepToxo/rds/S.O_ATAC_labels_pt.rds')
DefaultAssay(S.O.rna) <- 'RNA'
DefaultAssay(S.O.atac) <- 'ACTIVITY'

cyclic.genes <- read.xlsx('../../Output/DeepToxo/tables/cyclic_genes.xlsx')

cyclic.genes <- cyclic.genes %>% dplyr::filter(rna.expressed == 1, atac.expressed == 1, rna.cyclic == 1, atac.cyclic == 1)



## This is used as a shortcut to find matching cells in the two dataset
## reduction can be set to 
transfer.anchors <- FindTransferAnchors(reference = S.O.rna, query = S.O.atac, 
                                        features = gsub('_', '-', cyclic.genes$GeneID),
                                        reference.assay = "RNA", query.assay = "ACTIVITY",
                                        reduction = 'cca')

match.inds <- data.frame(transfer.anchors@anchors) %>% group_by(cell1) %>% 
  summarise(cell2 = cell2[which.max(score)], score = score[which.max(score)]) %>% ungroup() %>% group_by(cell2) %>% 
  summarise(cell1 = cell1[which.max(score)])

rcs <- gsub('_reference', '', transfer.anchors@reference.cells)
qcs <- gsub('_query', '', transfer.anchors@query.cells)

rcs <- rcs[match.inds$cell1] ## take the anchers out from scRNA
qcs <- qcs[match.inds$cell2] ## take the anchers out from scATAC

DefaultAssay(S.O.rna) <- 'RNA'
X <- FetchData(S.O.rna, vars = gsub('_', '-', cyclic.genes$GeneID), cells = rcs, slot = 'data')
DefaultAssay(S.O.atac) <- 'ACTIVITY'
Y <- FetchData(S.O.atac, vars = gsub('_', '-', cyclic.genes$GeneID), cells = qcs, slot = 'data')


L0 <- ccapca(X, Y, a = 0, 3) ## CCA only
L1 <- ccapca(X, Y, a = 1, 3) ## PCA only
L5 <- ccapca(X, Y, a = 0.5, 3) ## CCA only

X <- scale(X, scale = F, center = T)
X.tilde <- as.matrix(X) %*% L1$u
plot(X.tilde[,1], X.tilde[,2], cex = 0.4, pch = 20)


## Project the entire data
X.all <- FetchData(S.O.rna, vars = gsub('_', '-', cyclic.genes$GeneID), slot = 'data')
X.all <- scale(X.all, scale = F, center = T)
X.tilde <- as.matrix(X.all) %*% L5$u
cols <- S.O.rna@meta.data$phase[match(rownames(X.tilde), S.O.rna@meta.data$Sample)]
plot(X.tilde[,1], X.tilde[,2], pch=20, cex = 0.4, col = as.factor(cols))

Y.all <- FetchData(S.O.atac, vars = gsub('_', '-', cyclic.genes$GeneID), slot = 'data')
Y.all <- scale(Y.all, scale = F, center = T)
Y.tilde <- as.matrix(Y.all) %*% L0$v
cols <- S.O.atac@meta.data$phase[match(rownames(Y.tilde), S.O.atac@meta.data$Sample)]
plot(Y.tilde[,1], Y.tilde[,2], pch = 20, cex = 0.4, col = as.factor(cols))


## Try with the select cells in a new Seurat Object.
s.o.atac.sub <- subset(S.O.atac, cells = qcs)
S.O.tmp <- CreateSeuratObject(counts = s.o.atac.sub@assays$ACTIVITY@counts)
S.O.tmp <- prep_S.O(S.O.tmp)
DimPlot(object = S.O.tmp, reduction = 'pca') + NoLegend()
DimPlot(object = S.O.atac, reduction = 'pca') + NoLegend()

## Try with the select genes in a new Seurat Object.
s.o.atac.sub <- subset(S.O.atac, features =gsub('_', '-', cyclic.genes$GeneID))
S.O.tmp <- CreateSeuratObject(counts = s.o.atac.sub@assays$ACTIVITY@counts)
S.O.tmp <- prep_S.O(S.O.tmp)
DimPlot(object = S.O.tmp, reduction = 'pca') + NoLegend()
DimPlot(object = S.O.atac, reduction = 'pca') + NoLegend()
## Seems like ATAC PCA directions learned from new cells (cca mathched with scRNA) 
## are significantly different from ATAC PCA using all cells. When using slected cells and
## projecting all data, we see a circular pattern, even for ALL cells. 
## That means the cells selected, induce a lower dimention, already influenced by priodicity
## in RNA. 
