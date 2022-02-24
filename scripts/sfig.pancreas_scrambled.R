rm(list=ls())
source(here::here("scripts/utils.R"))


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


sce.o <- qread(here::here("data/Pancreas/steady.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

pancreas_steady.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas (steady)")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))




sce.o <- qread(here::here("data/Pancreas/dynamical_scramble.o.qs"))

embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

pancreas_scramble.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas (dynamical) U scrambled")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/Pancreas/steady_scramble.o.qs"))
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

pancreas_steady_scramble.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 15,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.4,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas (steady) U scrambled")) +
	.theme_noframe + 
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0),
				legend.key.size = unit(6, "pt"))


mp <- plot_grid(pancreas.p, pancreas_scramble.p, 
								pancreas_steady.p, pancreas_steady_scramble.p,
								nrow = 2, ncol = 2, label_size = 10, labels = "auto")
save_plot(here::here("figs", "sfigs", "sfig.pancreas_scrambled.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2*1.5, ncol = 2*1.5, device = cairo_pdf)







