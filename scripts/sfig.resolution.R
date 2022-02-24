rm(list=ls())
source(here::here("scripts/utils.R"))



### simulation data 
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
adata$var["velocity_genes"] <- TRUE
scv$tl$velocity_graph(adata)
steady.o <- zellkonverter::AnnData2SCE(adata)


embedded <- embedVelocity(reducedDim(steady.o, "X_pca"), steady.o)

row1.lp <- lapply(c(10, 15, 20, 25), function(res) {
	p <- plotVelocityStream(steady.o, embedded[, 1:2], use.dimred = "X_pca", color_by = "true_t", color.alpha = 0.7, grid.resolution = res,
													stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
													arrow.angle = 18, arrow.length = 0.45, point.size = 5.01) +
		scale_color_gradientn(colours = c("#440154", "#482576", "#414487","#35608D",
																			"#2A788E", "#21908C","#22A884", "#43BF71",
																			"#7AD151", "#BBDF27", "#FDE725"),name = "t", limits = range(steady.o$true_t, na.rm = TRUE)) +
		labs( y = "PCA-2", x = "PCA-1") +
		ggtitle(str_c("Streamline (res=", res, ")")) +
		.theme_noframe + theme(legend.position = "none")

	
	return(p)
})

row2.lp <- lapply(c(10, 15, 20, 25), function(res) {
	idx <- which(!is.na(embedded[, 1]))
	grid.df <- gridVectors(reducedDim(steady.o, "X_pca")[idx, ], embedded[idx, ], resolution = res)
	p <- plotEmbExp(steady.o, "X_pca", color.v = steady.o$true_t, title = str_c("Grid method (res=", res, ")"), x_lab = "PCA-1", y_lab = "PCA-2", color.name = "time", point.size = 5.01, point.alpha = 0.6) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.2) +
		.theme_noframe + theme(legend.position = "none") 
	
	
	
	return(p)
})




## FUCCI data
sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))
### plot
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

row3.lp <- lapply(c(10, 15, 20, 25), function(res) {
	p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "fucci_time", color.alpha = 0.7, grid.resolution = res,
													stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
													arrow.angle = 18, arrow.length = 0.45) +
		scale_color_gradientn(colours = c("#2E22EA","#9E3DFB","#F86BE2",
																			"#FCCE7B","#C4E416","#4BBA0F",
																			"#447D87","#2C24E9"),name = "FUCCI", limits = range(sce.o$fucci_time, na.rm = TRUE)) +
		labs( y = "UMAP-2", x = "UMAP-1") +
		ggtitle(str_c("Streamline (res=", res, ")")) +
		.theme_noframe + theme(legend.position = "none")
	
	return(p)
})

row4.lp <- lapply(c(10, 20, 30, 40), function(res) {
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = res)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p <- plotEmbFucci(sce.o, "UMAP", title = str_c("Grid method (res=", res, ")"), x_lab = "UMAP-1", y_lab = "UMAP-2") + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.2) +
		.theme_noframe + theme(legend.position = "none") 
	return(p)
})



### pancreas data
sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))

### plot
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

row5.lp <- lapply(c(10, 15, 20, 25), function(res) {
	p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
													stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
													arrow.angle = 18, arrow.length = 0.45) +
		scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
		labs( y = "UMAP-2", x = "UMAP-1") +
		ggtitle(str_c("Streamline (res=", res, ")")) +
		.theme_noframe + theme(legend.position = "none")
	return(p)
})

row6.lp <- lapply(c(10, 20, 30, 40), function(res) {
	grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = res)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	p <- plotEmbScat(sce.o, "UMAP", title = str_c("Grid method (res=", res, ")"), x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
									 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.2) +
		.theme_noframe + theme(legend.position = "none") 
	return(p)
})

mp <- plot_grid( plotlist = c(row1.lp, row2.lp, row3.lp, row4.lp, row5.lp, row6.lp),
								nrow = 6, ncol = 4, label_size = 10, labels = as.vector(rbind(letters[1:6], matrix(rep("", 3 * 6), ncol = 6))))
save_plot(here::here("figs", "sfigs", "sfig.resolution.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 6, ncol = 4, device = cairo_pdf)


