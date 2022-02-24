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

set.seed(100)
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
														 point.size = 4.01) +
	geom_abline(slope = rowData(steady.o[i, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

tmp.df <- data.frame(x = dynamical.o$true_t, u = tmp.m[, 1], s = tmp.m[, 2], v = tmp.m[, 4],
										 ms = assay(dynamical.o, "Ms")[i, ],
										 mu = assay(dynamical.o, "Mu")[i, ],
										 spliced = assay(dynamical.o, "spliced")[i, ],
										 unspliced = assay(dynamical.o, "unspliced")[i, ],
										 steady_v = assay(steady.o, "velocity")[i,],
										 dynamical_v = assay(dynamical.o, "velocity")[i,])

ms_t.p <- plotxy(tmp.df$x, tmp.df$ms, title = "Ms ~ t", x_lab = "True t.",
								 y_lab = "Ms", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

mu_t.p <- plotxy(tmp.df$x, tmp.df$mu, title = "Mu ~ t", x_lab = "True t.",
								 y_lab = "Mu", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

spliced_t.p <- plotxy(tmp.df$x, tmp.df$spliced, title = "Spliced ~ t", x_lab = "True t.",
								 y_lab = "Spliced", color_by = "true_t", point.size = 4.01,
								 colors.v = tmp.df$x,
								 hue.colors = c("#440154", "#482576", "#414487","#35608D",
								 							 "#2A788E", "#21908C","#22A884", "#43BF71",
								 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = s), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

unspliced_t.p <- plotxy(tmp.df$x, tmp.df$unspliced, title = "Unspliced ~ t", x_lab = "True t.",
											y_lab = "Unspliced", color_by = "true_t", point.size = 4.01,
											colors.v = tmp.df$x,
											hue.colors = c("#440154", "#482576", "#414487","#35608D",
																		 "#2A788E", "#21908C","#22A884", "#43BF71",
																		 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = u), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))


steady_t.p <- plotxy(tmp.df$x, tmp.df$steady_v, title = "Steady v.", x_lab = "True t.",
										 y_lab = "Steady v.", color_by = "true_t", point.size = 4.01,
										 colors.v = tmp.df$x,
										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = v), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(legend.position = "none", plot.title = element_blank(),
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_t.p <- plotxy(tmp.df$x, tmp.df$dynamical_v, title = "Dynamical v.", x_lab = "True t.",
										 y_lab = "Dynamical v.", color_by = "true_t", point.size = 4.01,
										 colors.v = tmp.df$x,
										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	geom_path(data = tmp.df, aes(x = x, y = v), alpha = 0.7, linetype = "dashed", color = "black", size = 0.5) +
	.theme_noframe +
	theme(plot.title = element_blank(), legend.position = "none",
				plot.margin = unit(c(2, 2, 3, 3), "pt"))


steady_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True v.",
										 y_lab = "Steady v.", color_by = "true_t", point.size = 4.01,
										 colors.v = tmp.df$x,
										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	.theme_noframe +
	theme(plot.title = element_blank(), legend.position = "none",
				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_v.p <- plotxy(tmp.df$v, tmp.df$dynamical_v, title = "Dynamical v.", x_lab = "True v.",
										 y_lab = "Dynamical v.", color_by = "true_t", point.size = 4.01,
										 colors.v = tmp.df$x,
										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
	.theme_noframe +
	theme(plot.title = element_blank(), legend.position = "none",
				plot.margin = unit(c(2, 2, 3, 3), "pt"))


panel_b.p <- plot_grid(mums.p, spliced_t.p, unspliced_t.p,
											 ms_t.p, mu_t.p,
											 
										 dynamical_t.p,
										 dynamical_v.p,
													 nrow = 1, ncol = 7, label_size = 10, labels = c("b"))


### PCA streamline plot
adata <- zellkonverter::SCE2AnnData(steady.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
steady.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(steady.o, "X_pca"), steady.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
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

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_proj.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj.", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt"))




### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec.p <- plotEmbExp(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = c(1, 1),
				legend.justification = c(1, 1),
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

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_proj.p <- plotEmbExp(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj.", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt"))




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

### cell idx for cells in the red box
cell.idx <- which((reducedDim(sce.o, "X_pca")[, 1] >= 13.2) &
										(reducedDim(sce.o, "X_pca")[, 1] <= 16.7) &
										(reducedDim(sce.o, "X_pca")[, 2] >= -4.3) &
										(reducedDim(sce.o, "X_pca")[, 2] <= 0.9))

### plot
embedded <- embedVelocity(reducedDim(sce.o, "X_pca"), sce.o)
trans_embedded <- embedded
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
true_vec.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t,
												title = "True v. trans.pro", x_lab = "PCA-1", 
												y_lab = "PCA-2", color.name = "time", point.size = 4.01,
												point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	.theme_noframe +
	theme(legend.position = "none")

### make zoom fig
grid.idx <- which((grid.df$start.1 >= 13.2) &
										(grid.df$start.1 <= 16.7) &
										(grid.df$start.2 >= -4.3) &
										(grid.df$start.2 <= 0.9))

true_vec_small.p <- plotEmbExp(sce.o[, cell.idx], "X_pca", color.v = sce.o$true_t[cell.idx],
					 title = "True v. trans.pro", x_lab = "PCA-1", 
					 y_lab = "PCA-2", color.name = "time", point.size = 12.01,
					 point.alpha = 0.6) + 
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "time", limits = range(sce.o$true_t, na.rm = TRUE)) +
	geom_segment(data=grid.df[grid.idx, ], mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.06, "inches")), inherit.aes = FALSE, size = 0.5, alpha = 0.6)  +
	annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	xlim(c(12.0, 16.8)) + ylim(c(-4.4, -0.8)) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				plot.margin = unit(c(1, 1, 2, 2), "pt"))


### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t,
												title = "True v. projection", x_lab = "PCA-1",
												y_lab = "PCA-2", color.name = "t", point.size = 4.01,
												point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6) +
	annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	.theme_noframe +
	theme(legend.position = "none")

grid.idx <- which((grid.df$start.1 >= 13.2) &
										(grid.df$start.1 <= 16.7) &
										(grid.df$start.2 >= -4.3) &
										(grid.df$start.2 <= 0.9))

true_proj_small.p <- plotEmbExp(sce.o[, cell.idx], "X_pca", color.v = sce.o$true_t[cell.idx],
															 title = "True v. projection", x_lab = "PCA-1", 
															 y_lab = "PCA-2", color.name = "time", point.size = 12.01,
															 point.alpha = 0.6) + 
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "time", limits = range(sce.o$true_t, na.rm = TRUE)) +
	geom_segment(data=grid.df[grid.idx, ], mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.06, "inches")), inherit.aes = FALSE, size = 0.5, alpha = 0.6)  +
	annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	xlim(c(12.0, 16.8)) + ylim(c(-4.4, -0.8)) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				plot.margin = unit(c(1, 1, 2, 2), "pt"))






### compare the cells in red rectangle


tmp.df <- data.frame(x = sqrt((embedded[, 1]) ^2 + (embedded[, 2]) ^2), y = sqrt((trans_embedded[, 1]) ^2 + (trans_embedded[, 2]) ^2))
tmp.df$color <- "black"
tmp.df$color[cell.idx] <- "red"
tmp.df$x[which(is.na(tmp.df$x))] <- 0
tmp.df$y[which(is.na(tmp.df$y))] <- 0

length_com.p <- ggplot(data = tmp.df, aes(x = x, y = y, color = color)) +
	geom_scattermore(pointsize = 5.01,  alpha = 0.8) +
	scale_color_manual(values = c("black", "red"), name = "",
										 labels = c(str_c("Other cell (n=", ncol(sce.o) - length(cell.idx), ")"), str_c("Cells in the box (n=",length(cell.idx), ")"))) + 
	labs(x = "Vector length by projection", y = "Vector length by\ntansition probability", title = "Unscaled vector length per cell") +
	theme(legend.position = c(0, 1),
				legend.justification = c(0, 1),
				legend.title = element_blank())











### umap

reducedDim(steady.o, "UMAP") <- umap.o$layout

embedded <- embedVelocity(reducedDim(steady.o, "UMAP"), steady.o)
dimnames(embedded) <- dimnames(reducedDim(steady.o, "UMAP"))
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "UMAP")[idx, ], embedded[idx, ], resolution = 20)
steady_vec_umap.p <- plotEmbExp(steady.o, "UMAP", color.v = steady.o$true_t, title = "Steady trans.pro", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
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
new_umap.m <- predict(umap.o, new_pca.m)

### projection 
embedded <- new_umap.m - reducedDim(steady.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(steady.o, "UMAP"))

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "UMAP")[idx, ], embedded[idx, ], resolution = 20)
steady_proj_umap.p <- plotEmbExp(steady.o, "UMAP", color.v = steady.o$true_t, title = "Steady UMAP-transform", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt"))



### dynamical

reducedDim(dynamical.o, "UMAP") <- umap.o$layout

embedded <- embedVelocity(reducedDim(dynamical.o, "UMAP"), dynamical.o)
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "UMAP"))
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "UMAP")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec_umap.p <- plotEmbExp(dynamical.o, "UMAP", color.v = dynamical.o$true_t, title = "Dynamical trans.pro", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
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
new_umap.m <- predict(umap.o, new_pca.m)


### projection 
embedded <- new_umap.m - reducedDim(dynamical.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "UMAP"))

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "UMAP")[idx, ], embedded[idx, ], resolution = 20)
dynamical_proj_umap.p <- plotEmbExp(dynamical.o, "UMAP", color.v =dynamical.o$true_t, title = "Dynamical UMAP-transform", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "time", point.size = 4.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt"))






### real velocity projection
reducedDim(sce.o, "UMAP") <- umap.o$layout
splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)

cell.idx <- which((reducedDim(sce.o, "UMAP")[, 1] >= -3) &
										(reducedDim(sce.o, "UMAP")[, 1] <= 0.3) &
										(reducedDim(sce.o, "UMAP")[, 2] >= -10.5) &
										(reducedDim(sce.o, "UMAP")[, 2] <= -8))
										 
										 
### projection 
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "UMAP")[idx, ], embedded[idx, ], resolution = 30, scale = FALSE)
true_vec_umap.p <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "True v. trans.pro", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "time", point.size = 4.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6)  +
	annotate("rect", xmin= -3, xmax= 0.3, ymin= - 10.5 , ymax= - 8, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	.theme_noframe +
	theme(legend.position = "none")

### make zoom fig
grid.idx <- which((grid.df$start.1 >= -3) &
										(grid.df$start.1 <= 0.3) &
										(grid.df$start.2 >= -10.5) &
										(grid.df$start.2 <= -8))

true_vec_umap_small.p <- plotEmbExp(sce.o[, cell.idx], "UMAP", color.v = sce.o$true_t[cell.idx],
															 title = "True v. trans.pro", x_lab = "UMAP-1", 
															 y_lab = "UMAP-2", color.name = "time", point.size = 9.01,
															 point.alpha = 0.6) + 
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "time", limits = range(sce.o$true_t, na.rm = TRUE)) +
	geom_segment(data=grid.df[grid.idx, ], mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.06, "inches")), inherit.aes = FALSE, size = 0.5, alpha = 0.6)  +
	annotate("rect", xmin= -3, xmax= 0.3, ymin= - 10.5 , ymax= - 8, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	xlim(c(-3.1, 1.4)) + ylim(c(-10.6, -7.9)) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				plot.margin = unit(c(1, 1, 2, 2), "pt"))


### projection 
embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = 30, scale = FALSE)

true_proj_umap.p <- plotEmbExp(sce.o, "UMAP", color.v = sce.o$true_t, title = "True v. UMAP-transform", x_lab = "UMAP-1", y_lab = "UMAP-2", color.name = "t", point.size = 4.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.023, "inches")), inherit.aes = FALSE, size = 0.27, alpha = 0.6) +
	annotate("rect", xmin= -3, xmax= 0.3, ymin= - 10.5 , ymax= - 8, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	.theme_noframe +
	theme(legend.position = "none")

grid.idx <- which((grid.df$start.1 >= -3) &
										(grid.df$start.1 <= 0.3) &
										(grid.df$start.2 >= -10.5) &
										(grid.df$start.2 <= -8))

true_proj_umap_small.p <- plotEmbExp(sce.o[, cell.idx], "UMAP", color.v = sce.o$true_t[cell.idx],
																		title = "True v. projection", x_lab = "UMAP-1", 
																		y_lab = "UMAP-2", color.name = "time", point.size = 9.01,
																		point.alpha = 0.6) + 
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																		"#35608D", "#2A788E", "#21908C",
																		"#22A884", "#43BF71", "#7AD151",
																		"#BBDF27", "#FDE725"),name = "time", limits = range(sce.o$true_t, na.rm = TRUE)) +
	geom_segment(data=grid.df[grid.idx, ], mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.06, "inches")), inherit.aes = FALSE, size = 0.5, alpha = 0.6)  +
	annotate("rect", xmin= -3, xmax= 0.3, ymin= - 10.5 , ymax= - 8, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.5) + 
	xlim(c(-3.1, 1.4)) + ylim(c(-10.6, -7.9)) +
	.theme_noframe +
	theme(legend.position = "none",
				plot.title = element_blank(),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				plot.margin = unit(c(1, 1, 2, 2), "pt"))






small.p <- ggdraw(ggplot() + theme_nothing()) +
	draw_plot(true_vec_small.p, 0, 1, .55, .6, hjust = 0, vjust = 1) +
	draw_plot(true_proj_small.p, 1, 0, .55, .6, hjust = 1, vjust = 0)

small_umap.p <- ggdraw(ggplot() + theme_nothing()) +
	draw_plot(true_vec_umap_small.p, 0, 1, 1, .52, hjust = 0, vjust = 1) +
	draw_plot(true_proj_umap_small.p, 1, 0, 1, .52, hjust = 1, vjust = 0)







### pancreas umap
sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))
umap.o <- qread(here::here("data/Pancreas/umap.o.qs"))

res <- 15
### plot
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
theta.v <- as.numeric(circular::coord2rad(embedded))
names(theta.v) <- colnames(sce.o)

pancreas_full.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																			stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																			arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas full data")) +
	.theme_noframe + theme(legend.position = "none")




### subset
### remove Beta
idx <- which(sce.o$clusters != "Beta")
new.o <- sce.o[, idx]
new.o$clusters <- droplevels(new.o$clusters)

adata <- zellkonverter::SCE2AnnData(new.o)
scv <- reticulate::import("scvelo")
scv$pp$neighbors(adata, use_rep = "PCA")
scv$tl$velocity_graph(adata)
new.o <- zellkonverter::AnnData2SCE(adata)

embedded_subset <- embedVelocity(reducedDim(new.o, "UMAP"), new.o)
theta_subset.v <- as.numeric(circular::coord2rad(embedded_subset))
names(theta_subset.v) <- colnames(new.o)


pancreas_subset.p <- plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.5, color= "#4DAF4A", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas without Beta cells")) +
	.theme_noframe + theme(legend.position = "none")

pancreas_legend.p <- get_legend(plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																									 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																									 arrow.angle = 18, arrow.length = 0.35) +
																	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Celltype", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
																	guides(color = guide_legend(override.aes = list(size = 1, alpha = 0.9))) +
																	.theme_noframe + theme(legend.key.size = unit(6, "pt")))

### compare direction
cells.v <- colnames(sce.o)[which((sce.o$clusters == "Pre-endocrine") & 
																 	(reducedDim(sce.o, "UMAP")[, 1] >= 3) & 
																 	(reducedDim(sce.o, "UMAP")[, 1] <= 5 ) &
																 	(reducedDim(sce.o, "UMAP")[, 2] >= - 0.5 ) &
																 	(reducedDim(sce.o, "UMAP")[, 2] <= 1 ))]

data.df <- rbind(data.frame(theta = theta.v[cells.v], type = "full"),
								 data.frame(theta = theta_subset.v[cells.v], type = "sub"))
data.df$type <- factor(data.df$type)

library(circular)
bw <- 30
d1.o <- density.circular(circular(theta.v[cells.v]), bw = bw)
d2.o <- density.circular(circular(theta_subset.v[cells.v]), bw = bw)

tmp.df <- rbind(data.frame(x = as.numeric(d1.o$x), y = d1.o$y, type = "full"),
								data.frame(x = as.numeric(d2.o$x), y = d2.o$y, type = "sub"))
maxx1.v <- 	as.numeric(d1.o$x)[which.max(d1.o$y)]
maxx2.v <- 	as.numeric(d2.o$x)[which.max(d2.o$y)]
max.v <- max(tmp.df$y)

watson.two.test(x = circular(theta.v[cells.v]),y = circular(theta_subset.v[cells.v]))  # P-value < 0.001 

dir_polar.p <- ggplot(tmp.df, aes(x = x, y = y * 0.5 + 2)) +
	geom_segment(aes(x = 0, y = 2, xend = 2 * pi, yend = 2), size = 0.5, linetype = "dashed", alpha = 0.3) + 
	geom_path(aes(color = type), size = 0.5, alpha = 0.8) +
	geom_segment(aes(x = 0, y = 3.5, xend = 2 * pi, yend = 3.5), size = 0.5, linetype = "dashed", alpha = 0.3) + 
	coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "on") +
	labs(title = str_c("Density of directions (n=", length(cells.v), ")"),
			 # x = str_c("Watson wheeler Pval:", format.pval(watson.two.test(x = circular(theta.v[cells.v]),y = circular(theta_subset.v[cells.v]))$p.v, digits = 3))) +
	x = str_c("Watson two test Pval < 0.001")) +
	scale_color_manual(values = c("#E41A1C", "#4DAF4A"), name = "",  guide ="none") +
	ylim(c(0, 3.5)) + 
	annotate("segment", x = maxx1.v, xend = maxx1.v,
					 y = 0, yend = 2, alpha = 0.8,
					 colour = "#E41A1C", size = 0.8, arrow = arrow(length = unit(0.1, "inches"))) +
	annotate("segment", x = maxx2.v, xend = maxx2.v,
					 y = 0, yend = 2, alpha = 0.8,
					 colour = "#4DAF4A", size = 0.8, arrow = arrow(length = unit(0.1, "inches"))) +
	.theme_noframe + theme(legend.position = "none",
												 plot.margin = unit(c(1, 1, 2, 2), "pt"),
												 
												 axis.title.y = element_blank())

pancreas_full_zoom.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																			stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2, point.size = 5.01,
																			arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas full data")) +
	.theme_noframe + theme(legend.position = "none") +
	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) )) +
	annotate("segment", x = 4, xend = 4 + 5 * cos(maxx1.v),
					 y = 0.25, yend = 0.25 + 5 * sin(maxx1.v), alpha = 0.8,
					 colour = "#E41A1C", size = 1.5, arrow = arrow(length = unit(0.2, "inches")))

pancreas_subset_zoom.p <- plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																				point.size = 5.01,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.5, color= "#4DAF4A", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas without Beta cells")) +
	.theme_noframe + theme(legend.position = "none") +
	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) )) +
	annotate("segment", x = 4, xend = 4 + 5 * cos(maxx2.v),
					 y = 0.25, yend = 0.25 + 5 * sin(maxx2.v), alpha = 0.8,
					 colour = "#4DAF4A", size = 1.5, arrow = arrow(length = unit(0.2, "inches")))




### put panels together
panel_a.p <- plot_grid(pancreas_full_zoom.p,  dir_polar.p, pancreas_subset_zoom.p, pancreas_legend.p,
											 rel_widths = c(1, 0.7, 1, 0.3),
											 nrow = 1, ncol = 4, label_size = 10, labels = "a")

panel_cd.p <- plot_grid(true_vec_umap.p, true_proj_umap.p, small_umap.p, 
												nrow = 1, ncol = 3, label_size = 10, labels = c("c", "", "d"))


panel_efg.p <- plot_grid(true_vec.p, true_proj.p, small.p,
													 ggplot() + theme_nothing(), length_com.p, ggplot() + theme_nothing(),
													nrow = 2, ncol = 3, label_size = 10, labels = c("e", "", "f", "", "g", ""))



mp <- plot_grid(panel_a.p, panel_b.p, panel_cd.p, panel_efg.p, 
								rel_heights = c(1, 0.5, 1, 2),
										 nrow = 4, ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "main", "main.predict.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4.5, ncol = 3, device = cairo_pdf)
	
	
	
	
	

### companion SFig
mp <- plot_grid(
								dynamical_vec.p, dynamical_proj.p,
								dynamical_vec_umap.p, dynamical_proj_umap.p, 
								steady_vec.p, steady_proj.p,
								steady_vec_umap.p, steady_proj_umap.p,
								
													nrow = 4, ncol = 2, label_size = 10, labels = c("a", " ", "b", " ", "c", "", "d", ""))
save_plot(here::here("figs", "sfigs", "sfig.predict.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4, ncol = 2, device = cairo_pdf)






