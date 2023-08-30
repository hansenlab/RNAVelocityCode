rm(list = ls())
source(here::here("scripts/utils.R"))

p.l <- lapply(c(1, 3, 5, 7), function(noise_level) {
	n_obs = 500L
	n_vars = 10L
	
	
	scv <- reticulate::import("scvelo")
	adata <-
		scv$datasets$simulation(
			n_obs = n_obs,
			n_vars = n_vars,
			t_max = 25,
			alpha = 5,
			beta = .3,
			gamma = .5,
			noise_level = noise_level
		)
	ncomp <- as.integer(min(30L, n_vars - 1L))
	
	
	scv$pp$pca(adata, n_comps = ncomp)
	scv$pp$neighbors(adata)
	scv$pp$moments(adata)
	
	rotation.m <- adata$varm["PCs"]
	
	
	
	### get dynamical
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	
	### get steady
	scv$tl$velocity(adata, mode = "deterministic")
	steady.o <- zellkonverter::AnnData2SCE(adata)
	
	
	###
	### PCA streamline plot
	adata <- zellkonverter::SCE2AnnData(steady.o)
	adata$var["velocity_genes"] <- TRUE
	scv$tl$velocity_graph(adata)
	steady.o <- zellkonverter::AnnData2SCE(adata)
	
	embedded <- embedVelocity(reducedDim(steady.o, "X_pca"), steady.o)
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(steady.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	
	
	### steady projection
	data.m <- assay(steady.o, "X")
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	
	splusv.m <- assay(steady.o, "velocity") + assay(steady.o, "X")
	# splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx,]) %*% rotation.m[idx,])
	
	### projection
	embedded <- new_pca.m - reducedDim(steady.o, "X_pca")
	dimnames(embedded) <- dimnames(reducedDim(steady.o, "X_pca"))
	
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(steady.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	
	### dynamical
	adata <- zellkonverter::SCE2AnnData(dynamical.o)
	adata$var["velocity_genes"] <- TRUE
	scv$tl$velocity_graph(adata)
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	
	embedded <-
		embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
	v_tran_lknn.m <- embedded[, 1:2]
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(dynamical.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	
	
	
	### dynamical projection
	data.m <- assay(dynamical.o, "X")
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	
	splusv.m <-
		assay(dynamical.o, "velocity") + assay(dynamical.o, "X")
	# splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx,]) %*% rotation.m[idx,])
	
	### projection
	embedded <- new_pca.m - reducedDim(dynamical.o, "X_pca")
	dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "X_pca"))
	v_proj_lknn.m <- embedded[, 1:2]
	
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(dynamical.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	### real velocity projection
	v.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
		tmp.m <-
			dynamical_thet(
				alpha = rowData(dynamical.o)$true_alpha[i],
				beta = rowData(dynamical.o)$true_beta[i],
				gamma = rowData(dynamical.o)$true_gamma[i],
				scaling = rowData(dynamical.o)$true_scaling[i],
				t_ = rowData(dynamical.o)$true_t_[i],
				u0 = 0,
				s0 = 0,
				tmax = 25,
				new_t = dynamical.o$true_t
			)
		return(tmp.m[, "vt"])
	}))
	dimnames(v.m) <- dimnames(dynamical.o)
	
	sce.o <- dynamical.o
	assay(sce.o, "velocity") <- v.m
	
	adata <- zellkonverter::SCE2AnnData(sce.o)
	scv$tl$velocity_graph(adata)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
	# splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx,]) %*% rotation.m[idx,])
	
	### plot
	embedded <- embedVelocity(reducedDim(sce.o, "X_pca"), sce.o)
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(sce.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	### projection
	embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
	grid.df <-
		gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)
	
	
	
	### use true k-NN
	adata <-
		scv$datasets$simulation(
			n_obs = n_obs,
			n_vars = n_vars,
			t_max = 25,
			alpha = 5,
			beta = .3,
			gamma = .5,
			noise_level = noise_level
		)
	ncomp <- as.integer(min(30L, n_vars - 1L))
	
	### get true s
	sce.o <- zellkonverter::AnnData2SCE(adata)
	s.m <- t(sapply(seq_len(nrow(sce.o)), function(i) {
		dynamical_thet(
			alpha = rowData(sce.o)$true_alpha[i],
			beta = rowData(sce.o)$true_beta[i],
			gamma = rowData(sce.o)$true_gamma[i],
			scaling = rowData(sce.o)$true_scaling[i],
			t_ = rowData(sce.o)$true_t_[i],
			u0 = 0,
			s0 = 0,
			tmax = 25,
			new_t = sce.o$true_t
		)[, 2]
	}))
	dimnames(s.m) <- dimnames(sce.o)
	assay(sce.o, "X") <- s.m
	adata <- zellkonverter::SCE2AnnData(sce.o)
	
	### learn k-NN graph using true S
	scv$pp$pca(adata, n_comps = ncomp)
	scv$pp$neighbors(adata)
	scv$pp$moments(adata)
	
	### switch back and rerun PCA
	sce.o <- zellkonverter::AnnData2SCE(adata)
	reducedDim(sce.o, "true_pca") <- reducedDim(sce.o, "X_pca")
	assay(sce.o, "X") <- assay(sce.o, "spliced")
	adata <- zellkonverter::SCE2AnnData(sce.o)
	
	scv$pp$pca(adata, n_comps = ncomp)
	rotation.m <- adata$varm["PCs"]
	
	
	
	### get dynamical
	scv$tl$recover_dynamics(adata)
	scv$tl$velocity(adata, mode = "dynamical")
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	
	
	### get steady
	scv$tl$velocity(adata, mode = "deterministic")
	steady.o <- zellkonverter::AnnData2SCE(adata)
	
	
	###
	### PCA streamline plot
	adata <- zellkonverter::SCE2AnnData(steady.o)
	adata$var["velocity_genes"] <- TRUE
	scv$tl$velocity_graph(adata)
	steady.o <- zellkonverter::AnnData2SCE(adata)
	
	embedded <- embedVelocity(reducedDim(steady.o, "X_pca"), steady.o)
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(steady.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	
	### steady projection
	data.m <- assay(steady.o, "X")
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	
	splusv.m <- assay(steady.o, "velocity") + assay(steady.o, "X")
	# splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx,]) %*% rotation.m[idx,])
	
	### projection
	embedded <- new_pca.m - reducedDim(steady.o, "X_pca")
	dimnames(embedded) <- dimnames(reducedDim(steady.o, "X_pca"))
	
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(steady.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	### dynamical
	adata <- zellkonverter::SCE2AnnData(dynamical.o)
	adata$var["velocity_genes"] <- TRUE
	scv$tl$velocity_graph(adata)
	dynamical.o <- zellkonverter::AnnData2SCE(adata)
	
	embedded <-
		embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
	v_tran_tknn.m <- embedded[, 1:2]
	
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(dynamical.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	### dynamical projection
	data.m <- assay(dynamical.o, "X")
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	
	splusv.m <-
		assay(dynamical.o, "velocity") + assay(dynamical.o, "X")
	# splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx,]) %*% rotation.m[idx,])
	
	### projection
	embedded <- new_pca.m - reducedDim(dynamical.o, "X_pca")
	dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "X_pca"))
	v_proj_tknn.m <- embedded[, 1:2]
	
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(dynamical.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	### real velocity projection
	v.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
		tmp.m <-
			dynamical_thet(
				alpha = rowData(dynamical.o)$true_alpha[i],
				beta = rowData(dynamical.o)$true_beta[i],
				gamma = rowData(dynamical.o)$true_gamma[i],
				scaling = rowData(dynamical.o)$true_scaling[i],
				t_ = rowData(dynamical.o)$true_t_[i],
				u0 = 0,
				s0 = 0,
				tmax = 25,
				new_t = dynamical.o$true_t
			)
		return(tmp.m[, "vt"])
	}))
	dimnames(v.m) <- dimnames(dynamical.o)
	
	sce.o <- dynamical.o
	assay(sce.o, "velocity") <- v.m
	
	adata <- zellkonverter::SCE2AnnData(sce.o)
	scv$tl$velocity_graph(adata)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
	# splusv.m[which(splusv.m < 0)] <- 0
	data.m <- splusv.m
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx,]) %*% rotation.m[idx,])
	
	### plot
	embedded <- embedVelocity(reducedDim(sce.o, "X_pca"), sce.o)
	idx <- which(!is.na(embedded[, 1]))
	grid.df <-
		gridVectors(reducedDim(sce.o, "X_pca")[idx,], embedded[idx,], resolution = 20)
	
	
	### projection
	embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
	v_proj_true.m <- embedded[, 1:2]
	
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
	grid.df <-
		gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)
	
	cos_1.v <-
		sapply(seq_len(nrow(v_proj_true.m)), function(i)
			lsa::cosine(v_proj_true.m[i,], v_tran_lknn.m[i,]))
	cos_2.v <-
		sapply(seq_len(nrow(v_proj_true.m)), function(i)
			lsa::cosine(v_proj_true.m[i,], v_tran_tknn.m[i,]))
	cos_3.v <-
		sapply(seq_len(nrow(v_proj_true.m)), function(i)
			lsa::cosine(v_proj_true.m[i,], v_proj_lknn.m[i,]))
	cos_4.v <-
		sapply(seq_len(nrow(v_proj_true.m)), function(i)
			lsa::cosine(v_proj_true.m[i,], v_proj_tknn.m[i,]))
	
	
	### violin plots for cos. sim
	tmp.df <-
		data.frame(tl = cos_1.v,
							 tt = cos_2.v,
							 pl = cos_3.v,
							 pt = cos_4.v) %>%
		pivot_longer(cols = 1:4,
								 names_to = "type",
								 values_to = "cos")
	tmp.df$type <-
		factor(tmp.df$type, levels = c("tl",  "pl", "tt", "pt"))
	
	cos_violin.p <-
		ggplot(tmp.df, aes(x = type, y = cos, color = type)) +
		geom_violin(
			alpha = 0.5,
			scale = "width",
			adjust = .8,
			size = 0.5
		) +
		geom_quasirandom(
			varwidth = FALSE,
			size = 0.1,
			alpha = 0.5,
			width = 0.5,
			stroke = 0.6
		) +
		scale_color_brewer(palette = "Set1", guide = "none") +
		scale_x_discrete(
			labels = c(
				"Trans.\nlearned k-NN",
				"Proj.\nlearned k-NN",
				"Trans.\ntrue k-NN",
				"Proj.\ntrue k-NN"
			)
		) +
		labs(
			y = "Cosine similarity",
			title = str_c(
				"Cosine similarity to the \ntrue vector fields\nnoise level:",
				noise_level
			)
		) +
		theme(
			axis.title.x = element_blank(),
			axis.line.x = element_blank(),
			axis.text.x = element_text(
				size = 6,
				vjust = 0.5,
				hjust = 0.5,
				angle = 30
			)
		)
	cos_violin.p
})

mp <- plot_grid(
	plotlist = p.l,
	nrow = 2,
	ncol = 2,
	label_size = 10,
	labels = "auto"
)

save_plot(
	here::here("figs", "main", "main.emb_cos_all.pdf"),
	mp,
	base_height = 2 * 1.2 / 1.4,
	base_width = 2 * 1.4 / 1.4,
	nrow = 2 * 1.2,
	ncol = 2,
	device = cairo_pdf
)
