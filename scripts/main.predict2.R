rm(list=ls())
source(here::here("scripts/utils.R"))

noise_level = 1
n_obs = 500L
n_vars = 10L

scv <- reticulate::import("scvelo")
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))
scv$pp$pca(adata, n_comps = ncomp)
rotation.m <- adata$varm["PCs"]


umap.o <- umap::umap(adata$obsm["X_pca"][, seq_len(ncomp)], min_dist = 0.5, random_state = 123)

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)

### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)



### plot gene 1
i <- 1
tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
												gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
												t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)

mums.p <- plot_dynamicalmums(dynamical.o, rownm = as.character(i - 1), color_by = "true_t", hue.colors = c("#440154", "#482576", "#414487",
																																																					 "#35608D", "#2A788E", "#21908C",
																																																					 "#22A884", "#43BF71", "#7AD151",
																																																					 "#BBDF27", "#FDE725"),
														 title = str_c("Gene 1 (LL:", sprintf("%.3f", rowData(dynamical.o[i, ])$fit_likelihood),")"), color.name = "t",
														 point.size = 5.01) +
	geom_abline(slope = rowData(steady.o[i, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

tmp.df <- data.frame(x = dynamical.o$true_t, s = tmp.m[, 2], v = tmp.m[, 4],
										 ms = assay(dynamical.o, "Ms")[i, ],
										 steady_v = assay(steady.o, "velocity")[i,],
										 dynamical_v = assay(dynamical.o, "velocity")[i,])
ms_t.p <- plotxy(tmp.df$x, tmp.df$ms, title = "Ms ~ t", x_lab = "t",
								 y_lab = "Ms", color_by = "true_t", point.size = 5.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

steady_t.p <- plotxy(tmp.df$x, tmp.df$steady_v, title = "Steady v.", x_lab = "t",
										 y_lab = "Steady v.", color_by = "true_t", point.size = 5.01,
										 colors.v = tmp.df$x,
										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = v), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none", plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_t.p <- plotxy(tmp.df$x, tmp.df$dynamical_v, title = "Dynamical v.", x_lab = "t",
												y_lab = "Dynamical v.", color_by = "true_t", point.size = 5.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = v), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(plot.title = element_blank(), legend.position = "none",
				plot.margin = unit(c(2, 2, 3, 3), "pt"))


steady_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True v.",
										 y_lab = "Steady v.", color_by = "true_t", point.size = 5.01,
										 colors.v = tmp.df$x,
										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	.theme_noframe +
	theme(plot.title = element_blank(), legend.position = "none",
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_v.p <- plotxy(tmp.df$v, tmp.df$dynamical_v, title = "Dynamical v.", x_lab = "True v.",
												y_lab = "Dynamical v.", color_by = "true_t", point.size = 5.01,
												colors.v = tmp.df$x,
												hue.colors = c("#440154", "#482576", "#414487","#35608D",
																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																			 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	.theme_noframe +
	theme(plot.title = element_blank(), legend.position = "none",
				plot.margin = unit(c(2, 2, 3, 3), "pt"))


panel_a.p <- plot_grid(mums.p, ms_t.p,
											 steady_t.p, dynamical_t.p,
											 steady_v.p, dynamical_v.p,
											 nrow = 1, ncol = 6, label_size = 10, labels = c("a"))


### PCA streamline plot
adata <- zellkonverter::SCE2AnnData(steady.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
steady.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(steady.o, "X_pca"), steady.o)
steady_vec.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 8,
																	 stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																	 arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(steady.o$true_t, na.rm = TRUE)) +
	labs(  title = "Steady cos.sim") +
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))



### steady projection
data.m <- assay(steady.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(steady.o, "velocity") * 10 + assay(steady.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])

### projection 
embedded <- new_pca.m - reducedDim(steady.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(steady.o, "X_pca"))

steady_proj.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 8,
																		stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																		arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(steady.o$true_t, na.rm = TRUE)) +
	labs(  title = "Steady proj.") +
	
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
dynamical_vec.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 8,
																			stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																			arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(dynamical.o$true_t, na.rm = TRUE)) +
	labs(  title = "Dynamical cos.sim") +
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))



### dynamical projection
data.m <- assay(dynamical.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(dynamical.o, "velocity") * 10 + assay(dynamical.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])

### projection 
embedded <- new_pca.m - reducedDim(dynamical.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "X_pca"))

dynamical_proj.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 8,
																			 stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																			 arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(dynamical.o$true_t, na.rm = TRUE)) +
	labs(  title = "Dynamical proj.") +
	
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))


panel_b.p <- plot_grid(steady_vec.p, steady_proj.p,
											 nrow = 1, ncol = 2, label_size = 10, labels = c("b"))
panel_c.p <- plot_grid(dynamical_vec.p, dynamical_proj.p,
											 nrow = 1, ncol = 2, label_size = 10, labels = c("c"))

### real velocity projection
v.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
	tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
													gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
													t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)
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
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])

### plot
embedded <- embedVelocity(reducedDim(sce.o, "X_pca"), sce.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
panel_d.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. cos.sim", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.015, "inches")), inherit.aes = FALSE, size = 0.27)  +
	annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1))



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

panel_e.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.015, "inches")), inherit.aes = FALSE, size = 0.27) +
	annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none")


panel_bc.p <- plot_grid(panel_b.p, panel_c.p,
												nrow = 2, ncol = 1, label_size = 10, labels = NULL)
panel_bcde.p <- plot_grid(panel_bc.p, panel_d.p, panel_e.p,
													nrow = 1, ncol = 3, label_size = 10, labels = c("", "d", "e"))





### umap

reducedDim(steady.o, "UMAP") <- umap.o$layout

embedded <- embedVelocity(reducedDim(steady.o, "UMAP"), steady.o)
dimnames(embedded) <- dimnames(reducedDim(steady.o, "UMAP"))
steady_vec_umap.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 13,
																				stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																				arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(steady.o$true_t, na.rm = TRUE)) +
	labs(  title = "Steady cos.sim") +
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))



### steady projection
data.m <- assay(steady.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(steady.o, "velocity") * 10  + assay(steady.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)

### projection 
embedded <- new_umap.m - reducedDim(steady.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(steady.o, "UMAP"))

steady_proj_umap.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 13,
																				 stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																				 arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(steady.o$true_t, na.rm = TRUE)) +
	labs(  title = "Steady proj.") +
	
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))


### dynamical

reducedDim(dynamical.o, "UMAP") <- umap.o$layout

embedded <- embedVelocity(reducedDim(dynamical.o, "UMAP"), dynamical.o)
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "UMAP"))
dynamical_vec_umap.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 13,
																					 stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																					 arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(dynamical.o$true_t, na.rm = TRUE)) +
	labs(  title = "Dynamical cos.sim") +
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))



### dynamical projection
data.m <- assay(dynamical.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(dynamical.o, "velocity") * 10 + assay(dynamical.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)


### projection 
embedded <- new_umap.m - reducedDim(dynamical.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "UMAP"))

dynamical_proj_umap.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 13,
																						stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																						arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(dynamical.o$true_t, na.rm = TRUE)) +
	labs(  title = "Dynamical proj.") +
	
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))


panel_f.p <- plot_grid(steady_vec_umap.p, steady_proj_umap.p,
											 nrow = 1, ncol = 2, label_size = 10, labels = c("f"))
panel_g.p <- plot_grid(dynamical_vec_umap.p, dynamical_proj_umap.p,
											 nrow = 1, ncol = 2, label_size = 10, labels = c("g"))




### real velocity projection
reducedDim(sce.o, "UMAP") <- umap.o$layout
splusv.m <- assay(sce.o, "velocity") * 10 + assay(sce.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)

### projection 
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "UMAP")[idx, ], embedded[idx, ], resolution = 30, scale = FALSE)
panel_h.p <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "True v. cos.sim", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "time", point.size = 5.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.22)  +
	annotate("rect", xmin= -3, xmax= 0.3, ymin= - 10.5 , ymax= - 8, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none")



### projection 
embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = 30, scale = FALSE)

panel_i.p <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "True v. projection", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 5.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.22) +
	annotate("rect", xmin= -3, xmax= 0.3, ymin= - 10.5 , ymax= - 8, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none")


panel_fg.p <- plot_grid(panel_f.p, panel_g.p,
												nrow = 2, ncol = 1, label_size = 10, labels = NULL)
panel_fghi.p <- plot_grid(panel_fg.p, panel_h.p, panel_i.p,
													nrow = 1, ncol = 3, label_size = 10, labels = c("", "h", "i"))

mp <- plot_grid(panel_a.p, panel_bcde.p, panel_fghi.p,
								rel_heights = c(1, 2, 2),
								nrow = 3, ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "main", "main.predict_2.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2.5, ncol = 3, device = cairo_pdf)











