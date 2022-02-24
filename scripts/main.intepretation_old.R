rm(list=ls())
source(here::here("scripts/utils.R"))

noise_level = 5
n_obs = 500L
n_vars = 10L

scv <- reticulate::import("scvelo")
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))


scv$pp$pca(adata, n_comps = ncomp)
rotation.m <- adata$varm["PCs"]

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)

### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)


###
i <- 8
tmp.m <- dynamical_thet(alpha = rowData(dynamical.o)$true_alpha[i], beta = rowData(dynamical.o)$true_beta[i],
												gamma = rowData(dynamical.o)$true_gamma[i], scaling = rowData(dynamical.o)$true_scaling[i],
												t_ = rowData(dynamical.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = dynamical.o$true_t)
tmp.df <- data.frame(x = dynamical.o$true_t, s = tmp.m[, 2], v = tmp.m[, 4], u =tmp.m[, 1],
										 est = assay(dynamical.o, "fit_t")[i, ],
										 ms = assay(dynamical.o, "Ms")[i, ],
										 mu = assay(dynamical.o, "Mu")[i, ],
										 steady_v = assay(steady.o, "velocity")[i,],
										 dynamical_v = assay(dynamical.o, "velocity")[i,])

mums.p <- plot_dynamicalmums(dynamical.o, rownm = as.character(i - 1), color_by = "true_t", hue.colors = c("#440154", "#482576", "#414487",
																																																					 "#35608D", "#2A788E", "#21908C",
																																																					 "#22A884", "#43BF71", "#7AD151",
																																																					 "#BBDF27", "#FDE725"),
														 title = str_c("Noise level ", noise_level, " (LL:", sprintf("%.3f", rowData(dynamical.o[i, ])$fit_likelihood),")"), color.name = "True t",
														 point.size = 5.01) +
	geom_abline(slope = rowData(steady.o[i, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
	theme(legend.position = c(1, 0.05),
				legend.justification = c(1, 0))


ms_t.p <- plotTLoess(x = tmp.df$x, y = tmp.df$ms, x_lab = "True t", x.pos = 0.6, y.pos = 0.9, point.size = 5.01,
											 color_var = tmp.df$x, y_lab = "Ms", title = str_c("Noise level ", noise_level),
											 scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
											 																		colors = c("#440154", "#482576", "#414487","#35608D",
											 																							 "#2A788E", "#21908C","#22A884", "#43BF71",
											 																							 "#7AD151", "#BBDF27", "#FDE725")), addR2 = FALSE) +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")

mu_t.p <- plotTLoess(x = tmp.df$x, y = tmp.df$mu, x_lab = "True t", x.pos = 0.6, y.pos = 1.1, point.size = 5.01,
												color_var = tmp.df$x, y_lab = "Mu", title = str_c("Noise level ", noise_level),
												scale_color = scale_color_gradientn(name = "t", limits = range(tmp.df$x), 
																														colors = c("#440154", "#482576", "#414487","#35608D",
																																			 "#2A788E", "#21908C","#22A884", "#43BF71",
																																			 "#7AD151", "#BBDF27", "#FDE725")), addR2 = FALSE) +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "solid", color = "black", size = 0.5) +
	theme(legend.position = "none")

steadyv_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True velocity",
											y_lab = "Steady v.", color_by = "true_t", point.size = 5.01,
											colors.v = tmp.df$x,
											hue.colors = c("#440154", "#482576", "#414487","#35608D",
																		 "#2A788E", "#21908C","#22A884", "#43BF71",
																		 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none")

dynamicalv_v.p <- plotxy(tmp.df$v, tmp.df$dynamical_v, title = "Dynamical v.", x_lab = "True velocity",
											y_lab = "Dynamical v.", color_by = "true_t", point.size = 5.01,
											colors.v = tmp.df$x,
											hue.colors = c("#440154", "#482576", "#414487","#35608D",
																		 "#2A788E", "#21908C","#22A884", "#43BF71",
																		 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none")

est_t.p <- plotxy(tmp.df$x, tmp.df$est, title = "Estimated t.", x_lab = "True t.",
												 y_lab = "Estimated t.", color_by = "true_t", point.size = 5.01,
												 colors.v = tmp.df$x,
												 hue.colors = c("#440154", "#482576", "#414487","#35608D",
												 							 "#2A788E", "#21908C","#22A884", "#43BF71",
												 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	theme(legend.position = "none")



### PCA streamline plot
adata <- zellkonverter::SCE2AnnData(steady.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
steady.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(steady.o, "X_pca"), steady.o)
steady_vec.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 15,
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
												 plot.margin = unit(c(1, 1, 2, 2), "pt"))



### steady projection
data.m <- assay(steady.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(steady.o, "velocity") + assay(steady.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])

### projection 
embedded <- new_pca.m - reducedDim(steady.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(steady.o, "X_pca"))

steady_proj.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 15,
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
												 plot.margin = unit(c(1, 1, 2, 2), "pt"))


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
dynamical_vec.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = 15,
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
												 plot.margin = unit(c(1, 1, 2, 2), "pt"))



### dynamical projection
data.m <- assay(dynamical.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(dynamical.o, "velocity") + assay(dynamical.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])

### projection 
embedded <- new_pca.m - reducedDim(dynamical.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "X_pca"))

dynamical_proj.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution =15,
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
												 plot.margin = unit(c(1, 1, 2, 2), "pt"))


panel_g.p <- plot_grid(steady_vec.p, steady_proj.p,
											 nrow = 1, ncol = 2, label_size = 10, labels = c("g"))
panel_h.p <- plot_grid(dynamical_vec.p, dynamical_proj.p,
											 nrow = 1, ncol = 2, label_size = 10, labels = c("h"))

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
panel_i.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. cos.sim", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.015, "inches")), inherit.aes = FALSE, size = 0.20)  +

	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt"))



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

panel_j.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.015, "inches")), inherit.aes = FALSE, size = 0.20) +
	# annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt"))


panel_gh.p <- plot_grid(panel_g.p, panel_h.p,
												nrow = 2, ncol = 1, label_size = 10, labels = NULL)
panel_ghij.p <- plot_grid(panel_gh.p, panel_i.p, panel_j.p,
													nrow = 1, ncol = 3, label_size = 10, labels = c("", "i", "j"))

panel_abcdef.p <- plot_grid(mums.p, ms_t.p, mu_t.p, 
														steadyv_v.p, dynamicalv_v.p, est_t.p,
												nrow = 2, ncol = 3, label_size = 10, labels = "auto")

mp <- plot_grid(panel_abcdef.p, panel_ghij.p,
								rel_heights = c(2, 1.4),
														nrow = 2, ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "main", "main.intepretation.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3, ncol = 3, device = cairo_pdf)







