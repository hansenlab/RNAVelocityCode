rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")

for (k in c(30L, 50L, 100L, 250L, 500L, 800L)) {
	sce.o <- qread(here::here("data/dentategyrus_lamanno/sce.o.qs"))
	reducedDimNames(sce.o)[2] <- "X_pca"
	adata <- zellkonverter::SCE2AnnData(sce.o)
	scv$pp$neighbors(adata, n_neighbors = k)
	scv$pp$moments(adata,  n_neighbors = k)
	
	sce.o <- zellkonverter::AnnData2SCE(adata)
	reducedDimNames(sce.o)[2] <- "PCA"
	metadata(sce.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
	qsave(sce.o, file = here::here(str_c("data/dentategyrus_lamanno/sce_k", k, ".o.qs")))
	rm(sce.o)
	
	scv$tl$velocity(adata, mode = "deterministic")
	scv$tl$velocity_graph(adata)
	
	steady.o <- zellkonverter::AnnData2SCE(adata)
	reducedDimNames(steady.o)[2] <- "PCA"
	metadata(steady.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
	qsave(steady.o, here::here(str_c("data/dentategyrus_lamanno/steady_k", k, ".o.qs")))
	rm(steady.o)
	
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	scv$tl$velocity_graph(adata)
	
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	reducedDimNames(dynamical.o)[2] <- "PCA"
	metadata(dynamical.o)$cluster_colors <- c('#f2991a', '#d94c1a', '#cc051a' ,'#cf6eb8', '#9c21b8', '#e6cc4c', '#b2d199', '#72daf2', '#5966d1', '#3b4cb2', '#0d1c82', '#3387b5', '#1a734c', '#4c5980')
	qsave(dynamical.o, file = here::here(str_c("data/dentategyrus_lamanno/dynamical_k", k, ".o.qs")))
	
	rm(adata)
	rm(dynamical.o)
	gc()
	print(k)
}



