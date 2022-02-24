rm(list=ls())
source(here::here("scripts/utils.R"))

# sceasy::convertFormat(here::here('data/chromaffin/onefilepercell_A1_unique_and_others_J2CH1.loom'), from="loom", to="sce",
# 											outFile=here::here('data/chromaffin/chromaffin.rds'))

sce.o <-  as(readRDS(here::here('data/chromaffin/chromaffin.rds')), "SingleCellExperiment")
colnames(sce.o) <- gsub("_unique.bam","", gsub(".*:","", sce.o$CellID))
rownames(sce.o) <- rowData(sce.o)$Gene
cell.colors <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/cell.colors.rds"))
sce.o$cell_type <- factor(cell.colors[match(colnames(sce.o), names(cell.colors))])

### keep spliced and unspliced assays
assays(sce.o) <- assays(sce.o)[c("spliced", "unspliced")]
assays(sce.o) <- lapply(assays(sce.o), function(x) as(x, "dgCMatrix"))

### filter genes with more than 20 counts
idx <- which((rowSums(assay(sce.o, "spliced")) > 20) & (rowSums(assay(sce.o, "unspliced")) > 20))
sce.o <- sce.o[idx, ]

### filter cells with more than 200 genes
idx <- which((colSums(assay(sce.o, "spliced") > 0) > 200) & (colSums(assay(sce.o, "unspliced") > 0) > 200))
sce.o <- sce.o[, idx]


### scale total counts
sce.o <- logNormCounts(sce.o, assay.type= "spliced")


### Do PCA on spliced counts using HVG
dec <- modelGeneVar(sce.o)
top.hvgs <- getTopHVGs(dec, n=2000)

set.seed(100)
sce.o <- runPCA(sce.o, subset_row = top.hvgs)
sce.o <- runUMAP(sce.o, dimred="PCA")

### Do PCA on GO cell cycle genes
library(org.Mm.eg.db)
go.v <- AnnotationDbi::select(org.Mm.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

sce.o <- runPCA(sce.o, subset_row = rownames(sce.o) %in% go.v, name = "PCA.cc")
sce.o <- runUMAP(sce.o, dimred="PCA.cc", name = "UMAP.cc")

library(velociraptor)
steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = 5L)))
reducedDims(steady.o) <- reducedDims(sce.o)
rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))

dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = 5L)))
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))


qsave(sce.o, file = here::here("data/chromaffin/sce.o.qs"))
qsave(steady.o, file = here::here("data/chromaffin/steady.o.qs"))
qsave(dynamical.o, file = here::here("data/chromaffin/dynamical.o.qs"))




