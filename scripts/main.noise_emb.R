rm(list=ls())
source(here::here("scripts/utils.R"))

noise_level = 5
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
steady_vec.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
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
steady_proj.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj. (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
v_tran_lknn.m <- embedded[, 1:2]
idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)

### points with auxiliary line
aux.p <- plotEmbExp(dynamical.o, "X_pca", color.v = dynamical.o$true_t,
					 title = "Simulation with noise level 5", x_lab = "PCA-1",
					 y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_vline(xintercept = -3.5, color = "red", linetype = "dashed", size = 0.6) +
	geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.6) +
	# geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
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
v_proj_lknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_proj.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj. (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
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
true_vec.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. trans.pro (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6) +
	# annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)




### use true kNN
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

### learn kNN graph using true S
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
steady_vec_trueknn.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
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
steady_proj_trueknn.p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj. (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
v_tran_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec_trueknn.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
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
v_proj_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_proj_trueknn.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj. (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
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
true_vec_trueknn.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. trans.pro (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6)  +
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
v_proj_true.m <- embedded[, 1:2]

dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj_trueknn.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.05, "inches")), inherit.aes = FALSE, size = 0.3, alpha = 0.6) +
	# annotate("rect", xmin= 13.2, xmax=16.7, ymin= -4.3 , ymax=-0.9, alpha = 0.2, color="red", fill= NA, linetype = "dashed", size = 0.3) + 
	.theme_noframe +
	theme(legend.position = "none",
				plot.margin = unit(c(1, 1, 2, 2), "pt")) +
	geom_curve(
		data = data.frame(),
		aes(x = -11, y = 8, xend = -10, yend = -10),
		arrow = arrow(length=unit(0.2, "inches")),
		colour = "#EC7014",
		alpha = 0.6,
		size = 1.2,
		curvature = - 3,
		inherit.aes = FALSE
	)

cos_1.v <- sapply(seq_len(nrow(v_proj_true.m)), function(i) lsa::cosine(v_proj_true.m[i, ], v_tran_lknn.m[i, ]))
cos_2.v <- sapply(seq_len(nrow(v_proj_true.m)), function(i) lsa::cosine(v_proj_true.m[i, ], v_tran_tknn.m[i, ]))
cos_3.v <- sapply(seq_len(nrow(v_proj_true.m)), function(i) lsa::cosine(v_proj_true.m[i, ], v_proj_lknn.m[i, ]))
cos_4.v <- sapply(seq_len(nrow(v_proj_true.m)), function(i) lsa::cosine(v_proj_true.m[i, ], v_proj_tknn.m[i, ]))


### violin plots for cos. sim
tmp.df <- data.frame(tl = cos_1.v, tt = cos_2.v, pl = cos_3.v, pt = cos_4.v) %>%
	pivot_longer(cols = 1:4, names_to = "type", values_to = "cos")
tmp.df$type <- factor(tmp.df$type, levels = c("tl",  "pl", "tt", "pt"))

cos_violin.p <- ggplot(tmp.df, aes(x = type, y = cos, color = type)) +
	geom_violin(alpha = 0.5,  scale= "width", adjust = .8, size = 0.5) +
	geom_quasirandom(varwidth = FALSE, size = 0.1, alpha = 0.5, width = 0.5, stroke = 0.6) +
	scale_color_brewer(palette = "Set1", guide = "none") +
	scale_x_discrete(labels = c("(c) trans.\nlearned kNN",
															"(d) proj.\nlearned kNN",
															"(e) trans.\ntrue kNN",
															"(f) proj.\ntrue kNN")) +
	labs(y = "Cosine similarity", title = "Consine similarity to the \ntrue vector fields in (b)") +
	theme(axis.title.x = element_blank(),
				axis.line.x = element_blank(),
				axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5, angle = 30))


### PBMC 
pbmc68k.o <- qread(here::here("data/pbmc68k/dynamical.o.qs"))
embedded <- embedVelocity(reducedDim(pbmc68k.o, "TSNE"), pbmc68k.o)

pbmc68k.p <- plotVelocityStream(pbmc68k.o, embedded[, 1:2], use.dimred = "TSNE", color_by = "celltype", color.alpha = 0.8, grid.resolution = 20,
																stream.L = 4, stream.min.L = 1, stream.res = 10, stream.width = 0.4, point.size = 2.01,
																arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = brewer.pal(nlevels(factor(pbmc68k.o$celltype)), "Set3"), name = "Cell type",
										 labels = levels(factor(pbmc68k.o$celltype)), limits = levels(factor(pbmc68k.o$celltype))) +
	guides(color=guide_legend(ncol = 1)) + 
	labs( y = "TSNE-2", x = "TSNE-1") +
	ggtitle(str_c("PBMC68k (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = "right",
				legend.justification = c(0.5, 0.8),
				legend.key.size = unit(7, "pt"))

col1.p <- plot_grid(aux.p, dynamical_vec.p, dynamical_vec_trueknn.p,
										nrow = 3, ncol = 1, label_size = 10, labels = c("a", "c", "e"))
col2.p <- plot_grid(true_proj_trueknn.p + ggtitle("True v. projection"), dynamical_proj.p,  dynamical_proj_trueknn.p,
										nrow = 3, ncol = 1, label_size = 10, labels = c("b", "d", "f"))
col3.p <- plot_grid(cos_violin.p, pbmc68k.p + theme(legend.position = "none"), get_legend(pbmc68k.p),
										nrow = 3, ncol = 1, label_size = 10, labels = c("g", "h", ""))

# row1.p <- plot_grid(pbmc68k.p, aux.p, true_proj_trueknn.p + ggtitle("True v. projection"),
# 									 nrow = 1, ncol = 3, label_size = 10, labels = "auto")
# row23_1.p <- plot_grid(dynamical_vec.p, dynamical_proj.p,
# 										nrow = 2, ncol = 1, label_size = 10, labels = c("d", "f"))
# row23_2.p <- plot_grid(dynamical_vec_trueknn.p, dynamical_proj_trueknn.p,
# 											 nrow = 2, ncol = 1, label_size = 10, labels = c("e", "g"))
# row23_3.p <- plot_grid(ggplot() + theme_nothing(),  cos_violin.p, ggplot() + theme_nothing(),
# 											 nrow = 3, ncol = 1, rel_heights = c(0.5, 1, 0.5), label_size = 10, labels = c("", "h", ""))
# row23.p <- plot_grid(row23_1.p, row23_2.p, row23_3.p,
# 										 nrow = 1, ncol = 3, label_size = 10, labels = NULL)

mp <- plot_grid(col1.p, col2.p, col3.p,
								nrow = 1, ncol = 3, label_size = 10, labels = NULL)

save_plot(here::here("figs", "main", "main.noise_emb.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4, ncol = 3, device = cairo_pdf)






