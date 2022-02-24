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




### umap
reducedDim(steady.o, "UMAP") <- umap.o$layout

epsilon <- 0.001

### steady projection
data.m <- assay(steady.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(steady.o, "velocity") * epsilon  + assay(steady.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))

new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)
diff.m <- new_umap.m - reducedDim(steady.o, "UMAP")



### projection 
embedded <- diff.m / epsilon
dimnames(embedded) <- dimnames(reducedDim(steady.o, "UMAP"))

a.p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 18,
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

### dynamical projection
data.m <- assay(dynamical.o, "X")
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

splusv.m <- assay(dynamical.o, "velocity") * epsilon + assay(dynamical.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))

new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)
diff.m <- new_umap.m - reducedDim(dynamical.o, "UMAP")


### projection 
embedded <- diff.m  / epsilon
dimnames(embedded) <- dimnames(reducedDim(dynamical.o, "UMAP"))

b.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 18,
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



### real velocity projection
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

reducedDim(sce.o, "UMAP") <- umap.o$layout
splusv.m <- assay(sce.o, "velocity") * epsilon + assay(sce.o, "X")
# splusv.m[which(splusv.m < 0)] <- 0
data.m <- splusv.m

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))

new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
new_umap.m <- predict(umap.o, new_pca.m)
diff.m <- new_umap.m - reducedDim(sce.o, "UMAP")


### projection 
embedded <- diff.m / epsilon

idx <- which(!is.na(embedded[, 1]))
grid.df <- gridVectors(reducedDim(sce.o, "UMAP")[idx, ], embedded[idx, ], resolution = 30, scale = FALSE)
c.p <- plotVelocityStream(dynamical.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "true_t", color.alpha = 0.7, grid.resolution = 18,
																stream.L = 4, stream.min.L = 1, stream.res = 3, stream.width = 0.15,
																arrow.angle = 18, arrow.length = 0.2, point.size = 5.01) +
	scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																		"#2A788E", "#21908C","#22A884", "#43BF71",
																		"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(dynamical.o$true_t, na.rm = TRUE)) +
	labs(  title = "True v. proj.") +
	
	.theme_noframe + theme(legend.position = "none",
												 plot.title = element_text(size = 6),
												 axis.title.y = element_blank(),
												 axis.title.x = element_blank(),
												 plot.margin = unit(c(2, 2, 3, 3), "pt"))


mp <- plot_grid(a.p, b.p, c.p,
								nrow = 2, ncol = 2, label_size = 10, labels = NULL)

save_plot(here::here("figs",  str_c("epsilon_", epsilon, ".pdf")), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2, device = cairo_pdf)


