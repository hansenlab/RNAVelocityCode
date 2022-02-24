rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")
sce.o <- qread(here::here("data", "neurosphere", "wt1.o.qs"))
rowData(sce.o)$Accession <- str_split_fixed(rowData(sce.o)$AccessionVersion, fixed("."), 2)[, 1]

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
umap.m <- umap.o$layout

reducedDim(sce.o, "UMAP") <- umap.m
plotUMAP(sce.o)



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


dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))


qsave(sce.o, file = here::here("data/neurosphere_wt1/sce.o.qs"))
qsave(steady.o, file = here::here("data/neurosphere_wt1/steady.o.qs"))
qsave(dynamical.o, file = here::here("data/neurosphere_wt1/dynamical.o.qs"))
qsave(umap.o, file = here::here("data/neurosphere_wt1/umap.o.qs"))









rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")
sce.o <- qread(here::here("data", "neurosphere", "wt2.o.qs"))
rowData(sce.o)$Accession <- str_split_fixed(rowData(sce.o)$AccessionVersion, fixed("."), 2)[, 1]

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
umap.m <- umap.o$layout

reducedDim(sce.o, "UMAP") <- umap.m
plotUMAP(sce.o)



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


dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))


qsave(sce.o, file = here::here("data/neurosphere_wt2/sce.o.qs"))
qsave(steady.o, file = here::here("data/neurosphere_wt2/steady.o.qs"))
qsave(dynamical.o, file = here::here("data/neurosphere_wt2/dynamical.o.qs"))
qsave(umap.o, file = here::here("data/neurosphere_wt2/umap.o.qs"))





### together 
rm(list=ls())
source(here::here("scripts/utils.R"))

sce1.o <- qread(here::here("data/neurosphere_wt1/sce.o.qs"))
sce2.o <- qread(here::here("data/neurosphere_wt2/sce.o.qs"))

genes.v <- intersect(rownames(sce1.o), rownames(sce2.o))

sce.o <- cbind(sce1.o[genes.v, ], sce2.o[genes.v, ])

total.m <- assay(sce.o, "spliced") + assay(sce.o, "unspliced")

library(sva)
adjusted.m <- ComBat(log2(total.m + 1 ), batch = sce.o$sample, 
										 BPPARAM = SerialParam())

adjusted.m <- (2 ^ adjusted.m) - 1
adjusted.m[adjusted.m < 0] <- 0
adjusted.m[adjusted.m > max(total.m)] <- max(total.m)


altExp(sce.o, "unadjusted") <- sce.o

### replace spliced and unspliced counts
assay(sce.o, "spliced") <- (assay(altExp(sce.o, "unadjusted"), "spliced") / total.m) * adjusted.m
assay(sce.o, "spliced")[total.m == 0] <- 0
assay(sce.o, "spliced") <- drop0(assay(sce.o, "spliced"))
assay(sce.o, "unspliced") <- (assay(altExp(sce.o, "unadjusted"), "unspliced") / total.m) * adjusted.m
assay(sce.o, "unspliced")[total.m == 0] <- 0
assay(sce.o, "unspliced") <- drop0(assay(sce.o, "unspliced"))


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
umap.m <- umap.o$layout

reducedDim(sce.o, "UMAP") <- umap.m
plotUMAP(sce.o)



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


dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
reducedDims(dynamical.o) <- reducedDims(sce.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))


qsave(sce.o, file = here::here("data/neurosphere/sce.o.qs"))
qsave(steady.o, file = here::here("data/neurosphere/steady.o.qs"))
qsave(dynamical.o, file = here::here("data/neurosphere/dynamical.o.qs"))
qsave(umap.o, file = here::here("data/neurosphere/umap.o.qs"))





