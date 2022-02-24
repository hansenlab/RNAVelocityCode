rm(list=ls())
source(here::here("scripts/utils.R"))

### process pbmc68k as the Bergen review paper 
sceasy::convertFormat(here::here('data/pbmc68k/pbmc68k.loom'), from="loom", to="sce",
											outFile=here::here('data/pbmc68k/pbmc68k.rds'))

sce.o <-  as(readRDS(here::here('data/pbmc68k/pbmc68k.rds')), "SingleCellExperiment")


colnames(sce.o) <- sce.o$obs_names
rownames(sce.o) <- rowData(sce.o)$Gene

tsne.df <- read_csv(here::here('data/pbmc68k/tsne.csv'))[, 2:3] %>% as.data.frame()
colnames(tsne.df) <- c("tsne1", "tsne2")
rownames(tsne.df) <- colnames(sce.o)
reducedDim(sce.o, "TSNE") <- tsne.df

reducedDim(sce.o, "TSNE")[, 1] <- - reducedDim(sce.o, "TSNE")[, 1]

adata <- zellkonverter::SCE2AnnData(sce.o)

scv <- reticulate::import("scvelo")
scv$pp$remove_duplicate_cells(adata)

scv$pp$filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000L)
scv$pp$moments(adata)

sce.o <- zellkonverter::AnnData2SCE(adata)
qsave(sce.o, file = here::here("data/pbmc68k/sce.o.qs"))


scv$tl$velocity(adata, mode = "deterministic")
scv$tl$velocity_graph(adata)

steady.o <- zellkonverter::AnnData2SCE(adata)
reducedDimNames(steady.o)[2] <- "PCA"
qsave(steady.o, file = here::here("data/pbmc68k/steady.o.qs"))

scv$tl$recover_dynamics(adata, n_jobs = 5L)
scv$tl$velocity(adata, mode = "dynamical")
scv$tl$velocity_graph(adata)

dynamical.o <- zellkonverter::AnnData2SCE(adata)
reducedDimNames(dynamical.o)[2] <- "PCA"
qsave(dynamical.o, file = here::here("data/pbmc68k/dynamical.o.qs"))




