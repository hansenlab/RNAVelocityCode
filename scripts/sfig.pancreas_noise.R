rm(list=ls())
source(here::here("scripts/utils.R"))

scv <- reticulate::import("scvelo")
sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))

embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

pancreas.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas (dynamical)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))


grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
pancreas_zoom.p <- plotEmbScat(sce.o, "UMAP", title = "Zoom-in of right part", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
																			 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.4) +
	.theme_noframe + theme(legend.position = "none") +
	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) ))


#### add noise
velocity.m <- assay(sce.o, "velocity")
add_noise <- function(velocity.m, alpha) {
	n <- ncol(velocity.m)
	row.idx <- which(!rowAnyNAs(velocity.m))
	sd.v <- rowSds(velocity.m[row.idx, ])
	noise.m <- t(sapply(seq_along(sd.v), function(i) rnorm(n = n, mean = 0, sd = alpha * sd.v[i])))
	velocity.m[row.idx, ] <- velocity.m[row.idx, ] + noise.m
	return(velocity.m)
}

noise.lp <- lapply(c(0.01, 0.1, 0.5, 1, 5), function(alpha) {
	set.seed(alpha + 100)
	assay(sce.o, "velocity") <- add_noise(velocity.m, alpha)
	adata <- zellkonverter::SCE2AnnData(sce.o)
	scv$tl$velocity_graph(adata)
	sce.o <- zellkonverter::AnnData2SCE(adata)
	
	embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
	
  p1 <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																		stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																		arrow.angle = 18, arrow.length = 0.35) +
		scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
		labs( y = "UMAP-2", x = "UMAP-1") +
		ggtitle(str_c("Noise alpha: ", alpha)) +
		.theme_noframe + 
		theme(legend.position = "none",
					legend.justification = c(1, 0),
					legend.key.size = unit(6, "pt"))
  
  grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = 30)
  colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
  p2 <- plotEmbScat(sce.o, "UMAP", title = str_c("Noise alpha: ", alpha), x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
  															 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
  	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.4) +
  	.theme_noframe + theme(legend.position = "none") +
  	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
  									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) ))
	return(list(p1 = p1, p2 = p2))
})


mp <- plot_grid(pancreas.p, pancreas_zoom.p, noise.lp[[1]]$p1, noise.lp[[1]]$p2,
								noise.lp[[2]]$p1, noise.lp[[2]]$p2, noise.lp[[3]]$p1, noise.lp[[3]]$p2,
								noise.lp[[4]]$p1, noise.lp[[4]]$p2, noise.lp[[5]]$p1, noise.lp[[5]]$p2,
								 byrow = TRUE, 
								 nrow = 3, ncol = 4, label_size = 10, labels = "auto")
save_plot(here::here("figs", "sfigs", "sfig.pancreas_noise.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3*1.5, ncol = 4*1.5, device = cairo_pdf)


