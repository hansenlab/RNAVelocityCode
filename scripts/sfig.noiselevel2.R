rm(list=ls())
source(here::here("scripts/utils.R"))

noise_level = 2
n_obs = 500L
n_vars = 10L

scv <- reticulate::import("scvelo")
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))


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
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (learned KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



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
steady_proj.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj. (learned KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec.p <- plotEmbExp(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (learned KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



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
dynamical_proj.p <- plotEmbExp(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj. (learned KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)


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
true_vec.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. trans.pro (learned KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection (learned KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2) +
	# annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



### use true KNN
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))

### get true s
sce.o <- zellkonverter::AnnData2SCE(adata)
s.m <- t(sapply(seq_len(nrow(sce.o)), function(i) {
	dynamical_thet(alpha = rowData(sce.o)$true_alpha[i], beta = rowData(sce.o)$true_beta[i],
								 gamma = rowData(sce.o)$true_gamma[i], scaling = rowData(sce.o)$true_scaling[i],
								 t_ = rowData(sce.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = sce.o$true_t)[, 2]
}))
dimnames(s.m) <- dimnames(sce.o)
assay(sce.o, "X") <- s.m
adata <- zellkonverter::SCE2AnnData(sce.o)

### learn KNN graph using true S
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
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec_trueknn.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (true KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



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
steady_proj_trueknn.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj. (true KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec_trueknn.p <- plotEmbExp(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (true KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



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
dynamical_proj_trueknn.p <- plotEmbExp(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj. (true KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)


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
true_vec_trueknn.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. trans.pro (true KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj_trueknn.p <- plotEmbExp(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection (true KNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2) +
	# annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = 15, y = 10, xend = 0, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = 1.5,
		inherit.aes = FALSE
	)




mp <- plot_grid(steady_vec.p, steady_proj.p, steady_vec_trueknn.p, steady_proj_trueknn.p,
								dynamical_vec.p, dynamical_proj.p, dynamical_vec_trueknn.p, dynamical_proj_trueknn.p,
								true_vec.p, true_proj.p, true_vec_trueknn.p, true_proj_trueknn.p,
								nrow = 3, ncol = 4, label_size = 10, labels = "auto")

save_plot(here::here("figs", "sfigs", "sfig.noiselevel2.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4, ncol = 4, device = cairo_pdf)








