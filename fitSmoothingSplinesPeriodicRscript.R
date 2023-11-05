library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

options <- commandArgs(trailingOnly = TRUE)

## in.file <- '../../Input/DeepToxo/rds/sc_rna_genes_expr_pt.rds'
in.file    <- options[1]
num.cores  <- as.numeric(options[2])
out.file   <- options[3]


sc.pt <- readRDS(in.file)



# Common genes between the data sets

comm.genes <- unique(sc.pt$GeneID)


# Expression

lbx <- 1.0
sc.spline.fits <- mclapply(1:5, function(i){
  tmp <- sc.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  
  
  y <- tmp$y
  t <- tmp$x
  
  
  ### Enforcing preiodicity
  ntimes <- length(t)
  y.ext <- c(y, y, y)
  t.ext <- c(t, t+6, t+12) 
  w.ext <- rep(1, length(y.ext))
  w.ext[which(y.ext == 0)] <- 1/5
  sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, spar = lbx, w = w.ext)
  #sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, cv = T, w = w.ext)
  #sparx <- sc.rna.sp.ext$spar
  #if(sparx < lbx){
  #  sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, spar = lbx, w = w.ext)
  #  
  #}
  sc.rna.sp.ext <- predict(sc.rna.sp.ext, seq(6, 12, by = 1/3)) 
  sc.rna.sp <- sc.rna.sp.ext
  sc.rna.sp$x <- round(sc.rna.sp$x - 6, 4)
  #plot(t, y)
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'red')
  
  ######
  
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.spline.fits <- bind_rows(sc.spline.fits)

cat('\nprocessing done\n')
saveRDS(sc.spline.fits, out.file)
