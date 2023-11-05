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
S.O.rna <- readRDS('../../Input/DeepToxo/rds/S.O.RNA.phase.pt.rds')
S.O.atac <- readRDS('../../Input/DeepToxo/rds/S.O_ATAC_labels_pt.rds')
DefaultAssay(S.O.rna) <- 'RNA'
DefaultAssay(S.O.atac) <- 'ACTIVITY'

cyclic.genes <- read.xlsx('../../Output/DeepToxo/tables/cyclic_genes.xlsx')

cyclic.genes <- cyclic.genes %>% dplyr::filter(rna.expressed == 1, atac.expressed == 1, rna.cyclic == 1, atac.cyclic == 1)


S.O.rna.cyclic    <- subset(S.O.rna, features = gsub('_', '-', cyclic.genes$GeneID))
S.O.atac.cyclic   <- subset(S.O.atac, features = gsub('_', '-', cyclic.genes$GeneID))

transfer.anchors <- FindTransferAnchors(reference = S.O.rna.cyclic, query = S.O.atac.cyclic, 
                                        features = VariableFeatures(object = S.O.rna.cyclic),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

DefaultAssay(S.O.rna) <- 'RNA'
X <- FetchData(S.O.rna,vars = gsub('_', '-', cyclic.genes$GeneID), slot = 'data')
x <- t(as.matrix(X))
DefaultAssay(S.O.atac) <- 'ACTIVITY'
Y <- FetchData(S.O.atac, vars = gsub('_', '-', cyclic.genes$GeneID), slot = 'data')
y <- t(as.matrix(Y))

library(CCA)

x <- scale(x, scale = FALSE)
y <- scale(y, scale = FALSE)

##down sample for testing
comm.rows <- sample(1:nrow(x), 100)

x.ds <- x[comm.rows, sample(1:ncol(x), 200)]
y.ds <- y[comm.rows, sample(1:ncol(y), 200)]
cc_results <- cancor(x.ds,y.ds)

cc1_x <- x.ds %*% matrix(cc_results$xcoef[, 1], ncol = 1)





require(ggplot2)
require(GGally)
require(CCA)
require(CCP)
mm <- read.csv("https://stats.idre.ucla.edu/stat/data/mmreg.csv")
colnames(mm) <- c("Control", "Concept", "Motivation", "Read", "Write", "Math", 
                  "Science", "Sex")
summary(mm)

xtabs(~Sex, data = mm)


psych <- mm[, 1:3]
acad <- mm[, 4:8]

ggpairs(psych)
matcor(psych, acad)

cc1 <- cc(psych, acad)

# display the canonical correlations
cc1$cor









