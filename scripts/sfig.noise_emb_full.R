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
steady_tran_lknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec.p <- plotEmbExp_nopoint(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
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
steady_proj_lknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_proj.p <- plotEmbExp_nopoint(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj. (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
dynamical_tran_lknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
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
dynamical_proj_lknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_proj.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj. (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
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
true_tran_lknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
true_vec.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. trans.pro (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
true_proj_lknn.m <- embedded[, 1:2]

grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection (learned kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2) +
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
		curvature = - 2.5,
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
steady_tran_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec_trueknn.p <- plotEmbExp_nopoint(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
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
steady_proj_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_proj_trueknn.p <- plotEmbExp_nopoint(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady proj. (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
dynamical_tran_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec_trueknn.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
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
dynamical_proj_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_proj_trueknn.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v =dynamical.o$true_t, title = "Dynamical proj. (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
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
true_tran_tknn.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
true_vec_trueknn.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. trans.pro (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)



### projection 
embedded <- new_pca.m - reducedDim(sce.o, "X_pca")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "X_pca"))
true_proj_tknn.m <- embedded[, 1:2]

grid.df <- gridVectors(reducedDim(sce.o, "X_pca"), embedded, resolution = 20)

true_proj_trueknn.p <- plotEmbExp_nopoint(sce.o, "X_pca", color.v = sce.o$true_t, title = "True v. projection (true kNN)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "t", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2) +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)



### add 2 knn scenarios
### do smoothing with learned knn
adata <- scv$datasets$simulation(n_obs=n_obs, n_vars=n_vars, t_max=25, alpha=5, beta=.3, gamma=.5, noise_level=noise_level)
ncomp <- as.integer( min(30L, n_vars - 1L))

scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)

### emb is subsettable in python, but the operation is not allowed in R. So we use sce 

### get true s
sce.o <- zellkonverter::AnnData2SCE(adata)
s.m <- t(sapply(seq_len(nrow(sce.o)), function(i) {
	dynamical_thet(alpha = rowData(sce.o)$true_alpha[i], beta = rowData(sce.o)$true_beta[i],
								 gamma = rowData(sce.o)$true_gamma[i], scaling = rowData(sce.o)$true_scaling[i],
								 t_ = rowData(sce.o)$true_t_[i], u0 = 0, s0 = 0, tmax = 25, new_t = sce.o$true_t)[, 2]
}))
dimnames(s.m) <- dimnames(sce.o)
assay(sce.o, "X") <- s.m
reducedDim(sce.o, "noise_pca") <- reducedDim(sce.o, "X_pca")

adata <- zellkonverter::SCE2AnnData(sce.o)

### learn kNN graph using true S
scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)

sce.o <- zellkonverter::AnnData2SCE(adata)
reducedDim(sce.o, "true_pca") <- reducedDim(sce.o, "X_pca")
reducedDim(sce.o, "X_pca") <- reducedDim(sce.o, "noise_pca")
adata <- zellkonverter::SCE2AnnData(sce.o)

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
steady_tran_lknn2.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec2.p <- plotEmbExp_nopoint(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (learned+true)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)

### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
dynamical_tran_lknn2.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec2.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (learned+true)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)


### smoothing by true knn and trans.pro by learned knn
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
scv$pp$neighbors(adata)
rotation.m <- adata$varm["PCs"]
sce.o <- zellkonverter::AnnData2SCE(adata)

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
steady_tran_tknn2.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
steady_vec_trueknn2.p <- plotEmbExp_nopoint(steady.o, "X_pca", color.v = steady.o$true_t, title = "Steady trans.pro (true+learned)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)


### dynamical
adata <- zellkonverter::SCE2AnnData(dynamical.o)
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
dynamical.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(dynamical.o, "X_pca"), dynamical.o)
dynamical_tran_tknn2.m <- embedded[, 1:2]

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(dynamical.o, "X_pca")[idx, ], embedded[idx, ], resolution = 20)
dynamical_vec_trueknn2.p <- plotEmbExp_nopoint(dynamical.o, "X_pca", color.v = dynamical.o$true_t, title = "Dynamical trans.pro (true+learned)", x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
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
		curvature = - 2.5,
		inherit.aes = FALSE
	)




cos.m <- sapply(list(steady_tran_lknn.m, steady_proj_lknn.m,
										 steady_tran_tknn.m, steady_proj_tknn.m,
										 steady_tran_lknn2.m, steady_tran_tknn2.m,
						dynamical_tran_lknn.m, dynamical_proj_lknn.m,
						dynamical_tran_tknn.m, dynamical_proj_tknn.m,
						dynamical_tran_lknn2.m,  dynamical_tran_tknn2.m,
						true_tran_lknn.m, 
						true_tran_tknn.m), function(x) sapply(seq_len(nrow(true_proj_tknn.m)), function(i) lsa::cosine(true_proj_tknn.m[i, ], x[i, ])))
colnames(cos.m) <- str_c("t", 1:14)
tmp.df <- data.frame(cos.m) %>%
	pivot_longer(cols = 1:14, names_to = "type", values_to = "cos")
tmp.df$type <- factor(tmp.df$type)

cos_violin.p <- ggplot(tmp.df, aes(x = type, y = cos)) +
	geom_violin(alpha = 0.5,  scale= "width", adjust = .8, size = 0.5) +
	geom_quasirandom(varwidth = FALSE, size = 0.08, alpha = 0.5, width = 0.4, stroke = 0.3) +
	# scale_color_brewer(palette = "Set1", guide = "none") +
	scale_x_discrete(labels = c("(a) steady trans.\nlearned kNN",
															"(b) steady proj.\nlearned kNN",
															"(c) steady trans.\ntrue kNN",
															"(d) steady proj.\ntrue kNN",
															"(e) steady trans.\nlearned+true",
															"(f) steady trans.\ntrue+learned",
															"(g) dynamical trans.\nlearned kNN",
															"(h) dynamical proj.\nlearned kNN",
															"(i) dynamical trans.\ntrue kNN",
															"(j) dynamical proj.\ntrue kNN",
															"(k) dynamical trans.\nlearned+true",
															"(l) dynamical trans.\ntrue+learned",
															"(m) true.v trans.\nlearned kNN",
															"(o) true.v trans.\ntrue kNN"), limits = str_c("t", 1:14)) +
	labs(y = "Cosine similarity", title = "Consine similarity to the true vector fields in (n) or (p)") +
	theme(axis.title.x = element_blank(),
				axis.line.x = element_blank(),
				axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5, angle = 30))







mp <- plot_grid(steady_vec.p, steady_proj.p, steady_vec_trueknn.p, steady_proj_trueknn.p, 
								steady_vec2.p, steady_vec_trueknn2.p,
								dynamical_vec.p, dynamical_proj.p, dynamical_vec_trueknn.p, dynamical_proj_trueknn.p,
								dynamical_vec2.p, dynamical_vec_trueknn2.p,
								true_vec.p, true_proj.p, true_vec_trueknn.p, true_proj_trueknn.p,
								nrow = 4, ncol = 4, label_size = 10, labels = "auto")
mp <- plot_grid(mp, cos_violin.p,
								# plot_grid(ggplot() + theme_nothing(), cos_violin.p, ggplot() + theme_nothing(),
								# 					nrow = 1, ncol = 3, rel_widths = c(0.5, 2, 0.5), labels = NULL),
								nrow = 2, ncol = 1, rel_heights = c(4, 1), label_size = 10, labels = c("", "q"))

save_plot(here::here("figs", "sfigs", "sfig.noise_emb_full.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 4 / 3 * 5, ncol = 4, device = cairo_pdf)








