library(openxlsx)
library(ggplot2)
library(tidyverse)
library(Signac)
library(Seurat)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)
library(Seurat)


source('./util_funcs.R')

# Read scRAN-Seq data and convert TGGT1 Ids to TGME49
S.O <- readRDS('../../Input/DeepToxo/rds/S.O.RNA.phase.rds')

# IDs
prod.desc  <- read.xlsx('../../Input/toxo_genomics/genes/ProductDescription_ME49_65.xlsx')
TGGT1_ME49 <- read.xlsx('../../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


## prepare genome

ME49.fasta <- readDNAStringSet("../../Input/toxo_genomics/genome/ToxoDB-65_TgondiiME49_Genome.fasta")
chrs <- names(ME49.fasta)[grep("TGME49_chr", names(ME49.fasta))]

chr.len <- data.frame(chr = gsub(" ", "", unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 1))),
                      len = as.numeric(gsub('length=', '', unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 4)))))

txdb <- makeTxDbFromGFF(file="../../Input/toxo_genomics/genome/ToxoDB-65_TgondiiME49.gff",
                        dataSource="Toxodb",
                        organism="Toxoplasma")

trans_biotypes <- select(txdb, keys=keys(txdb, "TXID"), 
                         columns = "TXTYPE", keytype =  "TXID")

genome(txdb) <- 'ME49'

tx_trans <- exonsBy(txdb, by = "tx", use.names = TRUE)
tx_names <- names(tx_trans)
num.exons <- lapply(tx_trans, function(x) length(x))
tx_names <- rep(tx_names, unlist(num.exons))
tx_trans <- unlist(tx_trans)
tx_trans$tx_id <- tx_names
tx_trans$gene_id <- gsub('-t.*', '', tx_trans$tx_id)
tx_trans$gene_name <- tx_trans$gene_id
tx_trans$type <- 'exon'
tx_trans$gene_biotype <- 'protein_coding'
tx_trans$exon_name <- tx_trans$exon_rank

tmp <- chr.len$len[match(names(seqlengths(txdb)), chr.len$chr)]
names(tmp) <- names(seqlengths(txdb))
seqlengths(tx_trans) <- tmp
seqlevels(tx_trans)
#inds <- c(5,6,1,2,3,7,8,10,11,9,4,12,13,14) 
inds <- grep('TGME49', names(tmp))
inds <- inds[sort(as.numeric(as.roman(gsub('b', '', (gsub('a', '', gsub('TGME49_chr', '', names(tmp)[inds])))))), index.return = T)$ix]

tx_trans <- tx_trans[tx_trans@seqnames %in% seqlevels(tx_trans)[inds]]

seqlevels(tx_trans) <- seqlevels(tx_trans)[inds]
isCircular(tx_trans) <- rep(F, length(isCircular(tx_trans)))




## read scATAC data from cellranger 
counts <- Read10X_h5('../../Input/DeepToxo/scATACseq/Tgondii_ATAC_10_27_22_Deep/outs/filtered_peak_bc_matrix.h5')
metadata <- read.csv(
  file = "../../Input/DeepToxo/scATACseq/Tgondii_ATAC_10_27_22_Deep/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

metadata.filt <- metadata
metadata.filt$Sample <- rownames(metadata.filt)
metadata.filt <- metadata.filt[metadata.filt$Sample %in% colnames(counts), ]
peak_anno <- read.table("../../Input/DeepToxo/scATACseq/Tgondii_ATAC_10_27_22_Deep/outs/filtered_peak_bc_matrix/peaks.bed", 
                        sep = '\t', header = F, quote = "")
colnames(peak_anno) <- c('Chr', 'strt', 'stp')
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = seqinfo(tx_trans),
  fragments = '../../Input/DeepToxo/scATACseq/Tgondii_ATAC_10_27_22_Deep/outs/fragments.tsv.gz',
  min.cells = 5,
  min.features = 100
)

Tg_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata.filt
)


Tg_ATAC[['peaks']]
granges(Tg_ATAC)
annotations <- tx_trans

Annotation(Tg_ATAC) <- annotations
Tg_ATAC <- NucleosomeSignal(object = Tg_ATAC)
Tg_ATAC <- TSSEnrichment(object = Tg_ATAC, fast = FALSE)


# add blacklist ratio and a fraction of reads in peaks
Tg_ATAC$pct_reads_in_peaks <- Tg_ATAC$peak_region_fragments / Tg_ATAC$passed_filters * 100
Tg_ATAC$blacklist_ratio <- Tg_ATAC$blacklist_region_fragments / Tg_ATAC$peak_region_fragments

Tg_ATAC$high.tss <- ifelse(Tg_ATAC$TSS.enrichment > 3, 'High', 'Low')
TSSPlot(Tg_ATAC, group.by = 'high.tss') + NoLegend()

Tg_ATAC$nucleosome_group <- ifelse(Tg_ATAC$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = Tg_ATAC, group.by = 'nucleosome_group') # takes a while

VlnPlot(
  object = Tg_ATAC,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


Tg_ATAC <- subset(
  x = Tg_ATAC,
  subset = peak_region_fragments > 200 &
    peak_region_fragments < 6000 &
    pct_reads_in_peaks > 40 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# stats
dim(Tg_ATAC@assays$peaks)



Tg_ATAC <- RunTFIDF(Tg_ATAC)
Tg_ATAC <- FindTopFeatures(Tg_ATAC, min.cutoff = 'q0')
Tg_ATAC <- RunSVD(Tg_ATAC)

DepthCor(Tg_ATAC)

## remove highly correlating components
Tg_ATAC <- RunUMAP(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,5)])
Tg_ATAC <- FindNeighbors(object = Tg_ATAC, reduction = 'lsi', dims = seq(1:30)[-c(1,3)])
Tg_ATAC <- FindClusters(object = Tg_ATAC, verbose = FALSE, algorithm = 3)


DimPlot(object = Tg_ATAC, label = TRUE, reduction = 'umap') + NoLegend()

## Estimate RNA levels from atac-seq 
# quantify gene activity
gene.activities <- GeneActivity(Tg_ATAC, extend.upstream = 600,
                                extend.downstream = 200)

# add gene activities as a new assay
Tg_ATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(Tg_ATAC) <- "ACTIVITY"
Tg_ATAC <- NormalizeData(Tg_ATAC)
Tg_ATAC <- ScaleData(Tg_ATAC, features = rownames(Tg_ATAC))
Tg_ATAC <- RunPCA(Tg_ATAC, features = rownames(Tg_ATAC))

# Identify anchors
S.O <- readRDS( '../../Input/DeepToxo/rds/S.O.RNA.phase.pt.rds')

transfer.anchors <- FindTransferAnchors(reference = S.O, query = Tg_ATAC, features = VariableFeatures(object = S.O),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions1 <- TransferData(anchorset = transfer.anchors, refdata = S.O$phase,
                                     weight.reduction = Tg_ATAC[["lsi"]], dims = 2:30)
celltype.predictions1 <- celltype.predictions1 %>% transmute(phase = predicted.id)
celltype.predictions2 <- TransferData(anchorset = transfer.anchors, refdata = as.character(S.O$pt.shifted.scaled),
                                     weight.reduction = Tg_ATAC[["lsi"]], dims = 2:30)


celltype.predictions2 <- celltype.predictions2 %>% transmute(pt.shifted.scaled = as.numeric(predicted.id))
celltype.predictions <- cbind(celltype.predictions1, celltype.predictions2)

Tg_ATAC <- AddMetaData(Tg_ATAC, metadata = celltype.predictions)

Idents(Tg_ATAC) <- 'phase'
DefaultAssay(Tg_ATAC) <- 'Activity'
DimPlot(Tg_ATAC, reduction = 'umap', label = TRUE) + NoLegend()


## Pseudo time
genes.act <- as.matrix(Tg_ATAC@assays$ACTIVITY@data)

genes.df <- data.frame(GeneID = rownames(genes.act),
                       genes.act) %>%
  pivot_longer(-c(GeneID), names_to = 'Sample', values_to = 'log2.expr')

genes.df$Sample <- gsub('\\.', '-', genes.df$Sample)

tmp <- data.frame(Sample = Tg_ATAC@meta.data$Sample, pt.shifted.scaled = Tg_ATAC@meta.data$pt.shifted.scaled)
genes.df <- inner_join(tmp, genes.df, by = 'Sample')

saveRDS(genes.df, '../../Input/DeepToxo/rds/sc_atac_genes_activity_pt.rds')

## Integrate the RNA&ATAC for visualization
genes.use <- VariableFeatures(S.O)
refdata <- GetAssayData(S.O, assay = "RNA", slot = "data")[genes.use, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = Tg_ATAC[["lsi"]],
                           dims = 2:30)
Tg_ATAC[["RNA"]] <- imputation
saveRDS(Tg_ATAC, '../../Input/DeepToxo/rds/S.O_ATAC_labels_pt.rds')

S.O$orig.ident <- 'scRNA'
Tg_ATAC$orig.ident <- 'scATAC'
coembed <- merge(x = S.O, y = Tg_ATAC)
coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

Idents(coembed) <- 'orig.ident'
DimPlot(coembed, reduction = 'umap', label = T)
saveRDS(coembed, '../../Input/DeepToxo/rds/S.O_RNA_ATAC_coembed.rds')


#atac_sub <- subset(S.O.integrated, ident = 'scATAC')
#rna_sub <- subset(S.O.integrated, ident = 'intra')


