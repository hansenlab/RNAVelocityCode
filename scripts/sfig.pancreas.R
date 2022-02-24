rm(list=ls())
source(here::here("scripts/utils.R"))

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


pancreas_beta.p <- plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.5, color= "#4DAF4A", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas without Beta cells")) +
	.theme_noframe + theme(legend.position = "none")



### remove Pre-endocrine
idx <- which(sce.o$clusters != "Pre-endocrine")
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

grid.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded_subset, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
pancreas_preendocrine.p <- plotEmbScat(new.o, "UMAP", title = "Pancreas without Pre-endocrine cells", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
									color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2) +
	annotate("rect", xmin= 3, xmax= 7, ymin= - 1 , ymax= 1.5, alpha = 0.5, color= "black", fill= NA, linetype = "solid", size = 0.7) + 
	annotate("rect", xmin= 1, xmax= 4, ymin= - 7 , ymax= -4.5, alpha = 0.5, color= "black", fill= NA, linetype = "solid", size = 0.7) + 
	.theme_noframe + theme(legend.position = "none")



# Ngn3 high EP
idx <- which(sce.o$clusters != "Ngn3 high EP")
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

grid.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded_subset, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
pancreas_ngn3high.p <- plotEmbScat(new.o, "UMAP", title = "Pancreas without Ngn3 high EP cells", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
																			 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2) +
	annotate("rect", xmin= 0, xmax= 4.5, ymin= - 7 , ymax= -3.5, alpha = 0.5, color= "black", fill= NA, linetype = "solid", size = 0.7) + 
	.theme_noframe + theme(legend.position = "none")



pancreas_legend.p <- get_legend(plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																									 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																									 arrow.angle = 18, arrow.length = 0.35) +
																	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Celltype", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
																	guides(color = guide_legend(override.aes = list(size = 1, alpha = 0.9))) +
																	.theme_noframe + theme(legend.key.size = unit(6, "pt")))


mp <- plot_grid(pancreas_full.p, pancreas_beta.p,  pancreas_preendocrine.p, pancreas_ngn3high.p, 

								nrow = 2, ncol = 2, label_size = 10, labels = "auto")

save_plot(here::here("figs", "sfigs", "sfig.pancreas.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3, ncol = 3, device = cairo_pdf)










