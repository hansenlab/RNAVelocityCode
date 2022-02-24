rm(list=ls())
source(here::here("scripts/utils.R"))


sceasy::convertFormat(here::here('data/DentateGyrus/DentateGyrus.loom'), from="loom", to="sce",
											outFile=here::here('data/dentategyrus_lamanno/DentateGyrus.rds'))

sce.o <-  as(readRDS(here::here('data/dentategyrus_lamanno/DentateGyrus.rds')), "SingleCellExperiment")

colnames(sce.o) <- sce.o$CellID
rownames(sce.o) <- rowData(sce.o)$Gene

tsne.df <- read_csv(here::here('data/dentategyrus_lamanno/tsne.csv'))[, 2:3] %>% as.data.frame()
colnames(tsne.df) <- c("tsne1", "tsne2")
rownames(tsne.df) <- colnames(sce.o)
reducedDim(sce.o, "TSNE") <- tsne.df


adata <- zellkonverter::SCE2AnnData(sce.o)

scv <- reticulate::import("scvelo")
scv$pp$remove_duplicate_cells(adata)

scv$pp$filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000L)

#### 500 neighbors used in the La Manno paper
scv$pp$neighbors(adata, n_neighbors = 30L)
scv$pp$moments(adata,  n_neighbors = 30L)

sce.o <- zellkonverter::AnnData2SCE(adata)
reducedDimNames(sce.o)[2] <- "PCA"
metadata(sce.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
qsave(sce.o, file = here::here("data/dentategyrus_lamanno/sce.o.qs"))


scv$tl$velocity(adata, mode = "deterministic")
scv$tl$velocity_graph(adata) 

steady.o <- zellkonverter::AnnData2SCE(adata)
reducedDimNames(steady.o)[2] <- "PCA"
metadata(steady.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
qsave(steady.o, file = here::here("data/dentategyrus_lamanno/steady.o.qs"))

scv$tl$recover_dynamics(adata, n_jobs = 5L)
scv$tl$velocity(adata, mode = "dynamical")
scv$tl$velocity_graph(adata)

dynamical.o <- zellkonverter::AnnData2SCE(adata)
reducedDimNames(dynamical.o)[2] <- "PCA"
metadata(dynamical.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
qsave(dynamical.o, file = here::here("data/dentategyrus_lamanno/dynamical.o.qs"))






# 
# 
# ### keep spliced and unspliced assays
# assays(sce.o) <- assays(sce.o)[c("spliced", "unspliced")]
# assays(sce.o) <- lapply(assays(sce.o), function(x) as(x, "dgCMatrix"))
# #reducedDims(sce.o) <- list()
# 
# 
# ### filter genes with more than 20 counts
# idx <- which((rowSums(assay(sce.o, "spliced")) > 20) & (rowSums(assay(sce.o, "unspliced")) > 20))
# sce.o <- sce.o[idx, ]
# 
# ### filter cells with more than 200 genes
# idx <- which((colSums(assay(sce.o, "spliced") > 0) > 200) & (colSums(assay(sce.o, "unspliced") > 0) > 200))
# sce.o <- sce.o[, idx]
# 
# ### scale total counts
# sce.o <- logNormCounts(sce.o, assay.type= "spliced")
# 
# ## Do PCA on spliced counts using HVG
# dec <- modelGeneVar(sce.o)
# top.hvgs <- getTopHVGs(dec, n=2000)
# metadata(sce.o)$top.hvgs <- top.hvgs
# 
# set.seed(600)
# sce.o <- runPCA(sce.o, subset_row = top.hvgs)
# plotPCA(sce.o)
# ### do UMAP
# library(umap)
# umap.o <- umap::umap(reducedDim(sce.o, "PCA")[, 1:30], min_dist = 0.5)
# umap.m <- umap.o$layout
# 
# reducedDim(sce.o, "UMAP") <- umap.m
# plotUMAP(sce.o)
# 
# metadata(sce.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
# 
# ### Do PCA on GO cell cycle genes
# # library(org.Hs.eg.db)
# # go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
# # 															columns = "SYMBOL")[, "SYMBOL"]
# # 
# # sce.o <- runPCA(sce.o, subset_row = rownames(sce.o) %in% go.v, name = "PCA.cc")
# # sce.o <- runUMAP(sce.o, dimred="PCA.cc", name = "UMAP.cc")
# 
# library(velociraptor)
# steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA")
# reducedDims(steady.o) <- reducedDims(sce.o)
# rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
# colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))
# metadata(steady.o) <- c(metadata(steady.o), metadata(sce.o))
# 
# dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
# reducedDims(dynamical.o) <- reducedDims(sce.o)
# rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
# colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))
# metadata(dynamical.o) <- c(metadata(dynamical.o), metadata(sce.o))
# 
# qsave(sce.o, file = here::here("data/dentategyrus_lamanno/sce.o.qs"))
# qsave(steady.o, file = here::here("data/dentategyrus_lamanno/steady.o.qs"))
# qsave(dynamical.o, file = here::here("data/dentategyrus_lamanno/dynamical.o.qs"))
# 
# 


