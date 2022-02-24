rm(list=ls())
source(here::here("scripts/utils.R"))

dynamical.l <- lapply(c(30L, 50L, 100L, 250L, 500L, 800L), function(k) {
	sce.o <- qread(here::here(str_c("data/dentategyrus_lamanno/dynamical_k", k, ".o.qs")))
	
	embedded <- embedVelocity(reducedDim(sce.o, "TSNE"), sce.o)
	grid.df <- gridVectors(reducedDim(sce.o, "TSNE"), embedded, resolution = 40)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	
	sce.o$ClusterName <- factor(sce.o$ClusterName, levels = c("RadialGlia", "RadialGlia2", "ImmAstro", 
																														"GlialProg", "OPC", "nIPC", "Nbl1", "Nbl2",
																														"ImmGranule1", "ImmGranule2", "Granule", "CA", 
																														"CA1-Sub", "CA2-3-4"))
	
	p <- plotEmbScat(sce.o, "TSNE", title = str_c("Dynamical k=", k), x_lab = "TSNE-1", y_lab = "TSNE-2",
							color_by = "ClusterName", labels.v = levels(sce.o$clusters),  point.size = 2.01,
							color.name = "Cell_type", colors.v = metadata(sce.o)$cluster_colors) + 
		scale_color_manual(values = metadata(sce.o)$cluster_colors, name = "Cell type", labels = levels(sce.o$ClusterName), limits = levels(sce.o$ClusterName)) +
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.3) +
		.theme_noframe + 
		theme(legend.position = c(1, 1),
					legend.justification = c(1, 1),
					legend.key.size = unit(6, "pt"))
	
	if (k > 30) p <- p + 	theme(legend.position = "none")
	
	return(p)
})

steady.l <- lapply(c(30L, 50L, 100L, 250L, 500L, 800L), function(k) {
	sce.o <- qread(here::here(str_c("data/dentategyrus_lamanno/steady_k", k, ".o.qs")))
	
	embedded <- embedVelocity(reducedDim(sce.o, "TSNE"), sce.o)
	grid.df <- gridVectors(reducedDim(sce.o, "TSNE"), embedded, resolution = 40)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	
	sce.o$ClusterName <- factor(sce.o$ClusterName, levels = c("RadialGlia", "RadialGlia2", "ImmAstro", 
																														"GlialProg", "OPC", "nIPC", "Nbl1", "Nbl2",
																														"ImmGranule1", "ImmGranule2", "Granule", "CA", 
																														"CA1-Sub", "CA2-3-4"))
	
	p <- plotEmbScat(sce.o, "TSNE", title = str_c("Steady k=", k), x_lab = "TSNE-1", y_lab = "TSNE-2",
									 color_by = "ClusterName", labels.v = levels(sce.o$clusters),  point.size = 2.01,
									 color.name = "Cell_type", colors.v = metadata(sce.o)$cluster_colors) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.3) +
		.theme_noframe + 
		theme(legend.position = "none")
	
	return(p)
})


mp <- plot_grid( plotlist = c(dynamical.l, steady.l),
								 nrow = 4, ncol = 3, label_size = 10, labels = c("a", rep("", 5), 'b'))
save_plot(here::here("figs", "sfigs", "sfig.dentategyrus_lamanno_k.pdf"), mp,
					base_height = 2 *1.2 , base_width = 2*1.4 , nrow = 4 * 1.5, ncol = 3 * 1.5, device = cairo_pdf)
