rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")
adata <- scv$datasets$pancreas()
sce.o <- zellkonverter::AnnData2SCE(adata)

### keep spliced and unspliced assays
assays(sce.o) <- assays(sce.o)[c("spliced", "unspliced")]
assays(sce.o) <- lapply(assays(sce.o), function(x) as(x, "dgCMatrix"))
reducedDims(sce.o) <- list()

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
metadata(sce.o)$top.hvgs <- top.hvgs

set.seed(600)
sce.o <- runPCA(sce.o, subset_row = top.hvgs)
plotPCA(sce.o)
### do UMAP
library(umap)
umap.o <- umap::umap(reducedDim(sce.o, "PCA")[, 1:30], min_dist = 0.5)
umap.m <- umap.o$layout[, c(2, 1)]

reducedDim(sce.o, "UMAP") <- umap.m
plotUMAP(sce.o)


### scramble u
tmp <- assay(sce.o, "unspliced")[, sample(seq_len(ncol(sce.o)))]
dimnames(tmp) <- dimnames(assay(sce.o, "unspliced"))
assay(sce.o, "unspliced") <- tmp


### Do PCA on GO cell cycle genes
# library(org.Hs.eg.db)
# go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
# 															columns = "SYMBOL")[, "SYMBOL"]
# 
# sce.o <- runPCA(sce.o, subset_row = rownames(sce.o) %in% go.v, name = "PCA.cc")
# sce.o <- runUMAP(sce.o, dimred="PCA.cc", name = "UMAP.cc")

library(velociraptor)
steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA")
reducedDims(steady.o) <- reducedDims(sce.o)
rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))
metadata(steady.o) <- c(metadata(steady.o), metadata(sce.o)[c(1, 2, 6)])

dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))
metadata(dynamical.o) <- c(metadata(dynamical.o), metadata(sce.o)[c(1, 2, 6)])

qsave(sce.o, file = here::here("data/Pancreas/sce_scramble.o.qs"))
qsave(steady.o, file = here::here("data/Pancreas/steady_scramble.o.qs"))
qsave(dynamical.o, file = here::here("data/Pancreas/dynamical_scramble.o.qs"))
qsave(umap.o, file = here::here("data/Pancreas/umap_scramble.o.qs"))

