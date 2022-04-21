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

grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")

# pancreas_full.p <- plotEmbScat(sce.o, "UMAP", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
# 															 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) +
# 	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)  +
# 	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
# 	labs( y = "UMAP-2", x = "UMAP-1") +
# 	ggtitle(str_c("Pancreas full data")) +
# 	.theme_noframe + theme(legend.position = "none")



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
scv$tl$velocity_graph(adata, n_recurse_neighbors = 5)
new.o <- zellkonverter::AnnData2SCE(adata)

embedded_subset <- embedVelocity(reducedDim(new.o, "UMAP"), new.o)
theta_subset.v <- as.numeric(circular::coord2rad(embedded_subset))
names(theta_subset.v) <- colnames(new.o)

grid_subset.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded_subset, resolution = 30)
colnames(grid_subset.df) <- c("start.1", "start.2", "end.1", "end.2")

pancreas_subset.p <- plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.5, color= "#4DAF4A", fill= NA, linetype = "solid", size = 0.7) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas without Beta cells")) +
	.theme_noframe + theme(legend.position = "none")

# pancreas_subset.p <- plotEmbScat(new.o, "UMAP", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
# 																 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) +
# 	geom_segment(data=grid_subset.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)  +
# 	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
# 	labs( y = "UMAP-2", x = "UMAP-1") +
# 	ggtitle(str_c("Pancreas without Beta cells")) +
# 	.theme_noframe + theme(legend.position = "none")


pancreas_legend.p <- get_legend(plotEmbScat(new.o, "UMAP", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
																						color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) +
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
	labs(title = str_c("Density of directions\n(n=", length(cells.v), ")"),
			 # x = str_c("Watson wheeler Pval:", format.pval(watson.two.test(x = circular(theta.v[cells.v]),y = circular(theta_subset.v[cells.v]))$p.v, digits = 3))) +
			 x = str_c("Watson two test Pval < 0.001")) +
	scale_color_manual(values = c("#E41A1C", "turquoise1"), name = "",  guide ="none") +
	ylim(c(0, 3.5)) + 
	annotate("segment", x = maxx1.v, xend = maxx1.v,
					 y = 0, yend = 2, alpha = 0.8,
					 colour = "#E41A1C", size = 0.8, arrow = arrow(length = unit(0.1, "inches"))) +
	annotate("segment", x = maxx2.v, xend = maxx2.v,
					 y = 0, yend = 2, alpha = 0.8,
					 colour = "turquoise1", size = 0.8, arrow = arrow(length = unit(0.1, "inches"))) +
	.theme_noframe + theme(legend.position = "none",
												 plot.margin = unit(c(1, 1, 2, 2), "pt"),
												 
												 axis.title.y = element_blank())

pancreas_full_zoom.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																					 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2, point.size = 5.01,
																					 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.3, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas full data\n(right part zoom-in)")) +
	.theme_noframe + theme(legend.position = "none") +
	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) )) #+
# annotate("segment", x = 4, xend = 4 + 5 * cos(maxx1.v),
# 				 y = 0.25, yend = 0.25 + 5 * sin(maxx1.v), alpha = 0.8,
# 				 colour = "#E41A1C", size = 1.5, arrow = arrow(length = unit(0.2, "inches")))
# pancreas_full_zoom.p <- plotEmbScat(sce.o, "UMAP", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
# 						color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) +
# 	annotate("rect", xmin= 3, xmax= 5, ymin= - 0 , ymax= 1, alpha = 0.15, color= "#E41A1C", fill= NA, linetype = "solid", size = 0.7) + 
# 	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.3)  +
# 	labs( y = "UMAP-2", x = "UMAP-1") +
# 	ggtitle(str_c("Pancreas full data")) +
# 	.theme_noframe + theme(legend.position = "none") +
# 	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
# 									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) )) +
# 	annotate("segment", x = 4, xend = 4 + 5 * cos(maxx1.v),
# 					 					 y = 0.25, yend = 0.25 + 5 * sin(maxx1.v), alpha = 0.5,
# 					 					 colour = "#E41A1C", size = 1.5, arrow = arrow(length = unit(0.2, "inches")))



pancreas_subset_zoom.p <- plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																						 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																						 point.size = 5.01,
																						 arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.5, color= "turquoise1", fill= NA, linetype = "solid", size = 0.7) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas without Beta cells\n(right part zoom-in)")) +
	.theme_noframe + theme(legend.position = "none") +
	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) )) #+
# annotate("segment", x = 4, xend = 4 + 5 * cos(maxx2.v),
# 			 y = 0.25, yend = 0.25 + 5 * sin(maxx2.v), alpha = 0.8,
# 			 colour = "#4DAF4A", size = 1.5, arrow = arrow(length = unit(0.2, "inches")))


# pancreas_subset_zoom.p <- plotEmbScat(new.o, "UMAP", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
# 						color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) +
# 	annotate("rect", xmin= 3, xmax= 5, ymin= - 0 , ymax= 1, alpha = 0.15, color= "#4DAF4A", fill= NA, linetype = "solid", size = 0.7) + 
# 	geom_segment(data=grid_subset.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.3)  +
# 	
# 	labs( y = "UMAP-2", x = "UMAP-1") +
# 	ggtitle(str_c("Pancreas without Beta cells")) +
# 	.theme_noframe + theme(legend.position = "none") +
# 	coord_cartesian(xlim = c(-2, max(reducedDim(sce.o, "UMAP")[, 1])),
# 									ylim=c(-2.5, max(reducedDim(sce.o, "UMAP")[, 2]) )) +
# 	annotate("segment", x = 4, xend = 4 + 5 * cos(maxx2.v),
# 					 y = 0.25, yend = 0.25 + 5 * sin(maxx2.v), alpha = 0.5,
# 					 colour = "#4DAF4A", size = 1.5, arrow = arrow(length = unit(0.2, "inches")))

### put panels together
row1.p <- plot_grid(pancreas_full.p, pancreas_legend.p, dir_polar.p, 
										rel_widths = c(1, 0.3, 0.7),
										nrow = 1, ncol = 3, label_size = 10, labels = c("a", "", "d"))

row2.p <- plot_grid(pancreas_full_zoom.p,  pancreas_subset_zoom.p, 
											 nrow = 1, ncol = 2, label_size = 10, labels = c("b", "c"))



mp <- plot_grid(row1.p, row2.p,
								nrow = 2, ncol = 1, label_size = 10, labels = NULL)

save_plot(here::here("figs", "main", "main.pancreas.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2, device = cairo_pdf)


