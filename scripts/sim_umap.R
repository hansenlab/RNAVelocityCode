rm(list=ls())
source(here::here("scripts/utils.R"))

make_sim <- function(noise_level = 1, n_obs = 500L, n_vars = 10L) {
	### steady model
	scv <- reticulate::import("scvelo")
	adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
	ncomp <- as.integer( min(30L, n_vars - 1L))
	scv$pp$pca(adata, n_comps = ncomp)
	scv$tl$velocity(adata)
	# scv$tl$recover_dynamics(adata)
	# scv$tl$velocity(adata, mode = "dynamical")
	adata$var["velocity_genes"] <- TRUE
	scv$tl$velocity_graph(adata)
	scv$tl$velocity_embedding(adata, basis = "pca")
	sce.o <- zellkonverter::AnnData2SCE(adata)
	
	ncomp <- min(30, ncol(adata$varm["PCs"]))
	rotation.m <- adata$varm["PCs"][, seq_len(ncomp)]
	
	
	library(umap)
	umap.o <- umap::umap(reducedDim(sce.o, "X_pca")[, seq_len(ncomp)], min_dist = 0.5)
	reducedDim(sce.o, "UMAP") <- umap.o$layout
	# plotUMAP(sce.o, colour_by = "true_t")
	# plotReducedDim(sce.o, dimred = "X_pca",colour_by = "true_t")
	adata <- zellkonverter::SCE2AnnData(sce.o)

	
	### plot
	embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p1 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (steady model)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) +
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	# projection
	# sce.o <- logNormCounts(sce.o, assay.type= "spliced")
	data.m <- assay(sce.o, "X")
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
	splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
	
	new_umap.m <- predict(umap.o, new_pca.m)

	
	reducedDim(sce.o, "new.UMAP") <- new_umap.m
	
	### projection 
	embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
	
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p2 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (steady model projection)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	
	### dynamical model
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	adata$var["velocity_genes"] <- TRUE
	scv$tl$velocity_graph(adata)
	scv$tl$velocity_embedding(adata, basis = "pca")
	sce.o <- zellkonverter::AnnData2SCE(adata)

	
	### plot
	embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p3 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (dynamical model)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) +
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	# projection
	# sce.o <- logNormCounts(sce.o, assay.type= "spliced")
	data.m <- assay(sce.o, "X")
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
	splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
	new_umap.m <- predict(umap.o, new_pca.m)
	plot(new_umap.m)
	
	### projection 
	embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
	
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p4 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (dynamical model projection)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	### use true velocity
	assay(sce.o, "velocity") <- sweep(assay(sce.o, "unspliced"), MARGIN = 1, STATS = rowData(sce.o)$true_beta, FUN = "*") - 
		sweep(assay(sce.o, "spliced"), MARGIN = 1, STATS = rowData(sce.o)$true_gamma, FUN = "*")
	
	adata <- zellkonverter::SCE2AnnData(sce.o)
	scv$tl$velocity_graph(adata)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	
	
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
	splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
	
	new_umap.m <- predict(umap.o, new_pca.m)
	plot(new_umap.m)
	### plot
	embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
	idx <- which(!is.na(embedded[, 1]))
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP")[idx, ], embedded[idx, ])
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p5 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (true v)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	
	### projection 
	embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
	
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p6 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (true v projection)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	### true v 2
	
	v.m <- t(sapply(seq_len(nrow(sce.o)), function(i) {
		tmp.m <- dynamical_thet(alpha = rowData(sce.o)$true_alpha[i], beta = rowData(sce.o)$true_beta[i],
														gamma = rowData(sce.o)$true_gamma[i], scaling = rowData(sce.o)$true_scaling[i],
														t_ = rowData(sce.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = sce.o$true_t)
		return(tmp.m[, "vt"])
	}))
	dimnames(v.m) <- dimnames(sce.o)
	
	assay(sce.o, "velocity") <- v.m
	
	adata <- zellkonverter::SCE2AnnData(sce.o)
	scv$tl$velocity_graph(adata)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
	splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
	new_umap.m <- predict(umap.o, new_pca.m)

	
	### plot
	embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
	idx <- which(!is.na(embedded[, 1]))
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP")[idx, ], embedded[idx, ])
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p7 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (true v2)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	
	### projection 
	embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
	
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p8 <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "Vector fields (true v2 projection)", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 3.01) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2) 
	
	
	
	mp <- wrap_plots(list(p1, p2, p3, p4, p5, p6, p7, p8), ncol = 2)
	
	save_plot(here::here("figs", "sim_umap", str_c("nvars_", n_vars, "_nobs_", n_obs, "_noiselevel_",noise_level, ".pdf")), mp,
						base_height = 2, base_width = 2*1.4, nrow = 1.2 * 4, ncol = 2 * 1.2, device = cairo_pdf)
	
	
}

make_sim(1)
make_sim(2)
make_sim(3)
make_sim(4)
make_sim(5)
make_sim(6)
make_sim(10)
make_sim(20)

make_sim(1, n_obs = 800L,  n_vars = 300L)
make_sim(2, n_obs = 800L,  n_vars = 300L)
make_sim(3, n_obs = 800L,  n_vars = 300L)
make_sim(4, n_obs = 800L,  n_vars = 300L)
make_sim(5, n_obs = 800L,  n_vars = 300L)
make_sim(6, n_obs = 800L,  n_vars = 300L)
make_sim(10, n_obs = 800L,  n_vars = 300L)
make_sim(20, n_obs = 800L,  n_vars = 300L)
make_sim(30, n_obs = 800L,  n_vars = 300L)

