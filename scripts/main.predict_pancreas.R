
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


pancreas_subset.p <- plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																				stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																				arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
	annotate("rect", xmin= 3, xmax= 5, ymin= - 0.5 , ymax= 1, alpha = 0.5, color= "#4DAF4A", fill= NA, linetype = "solid", size = 0.7) + 
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle(str_c("Pancreas without Beta cells")) +
	.theme_noframe + theme(legend.position = "none")

pancreas_legend.p <- get_legend(plotVelocityStream(new.o, embedded_subset[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = res,
																									 stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
																									 arrow.angle = 18, arrow.length = 0.35) +
																	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Celltype", labels = levels(sce.o$clusters), limits = levels(sce.o$clusters)) +
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

dir_polar.p <- ggplot(tmp.df, aes(x = x, y = y * 0.5 + 2)) +
	geom_segment(aes(x = 0, y = 2, xend = 2 * pi, yend = 2), size = 0.5, linetype = "dashed", alpha = 0.3) + 
	geom_path(aes(color = type), size = 0.5, alpha = 0.8) +
	geom_segment(aes(x = 0, y = 3.5, xend = 2 * pi, yend = 3.5), size = 0.5, linetype = "dashed", alpha = 0.3) + 
	coord_polar(theta = "x", start = -pi / 2, direction = -1, clip = "on") +
	labs(title = str_c("Density of directions (n=", length(cells.v), ")"),
			 x = str_c("Watson wheeler Pval:", format.pval(watson.wheeler.test(x = circular(data.df$theta), group = data.df$type)$p.v, digits = 3))) +
	scale_color_manual(values = c("#E41A1C", "#4DAF4A"), name = "",  guide ="none") +
	ylim(c(0, 3.5)) + 
	annotate("segment", x = maxx1.v, xend = maxx1.v,
					 y = 0, yend = 3.5, alpha = 0.8,
					 colour = "#E41A1C", size = 0.8, arrow = arrow(length = unit(0.1, "inches"))) +
	annotate("segment", x = maxx2.v, xend = maxx2.v,
					 y = 0, yend = 3.5, alpha = 0.8,
					 colour = "#4DAF4A", size = 0.8, arrow = arrow(length = unit(0.1, "inches"))) +
	.theme_noframe + theme(legend.position = "none",
												 plot.margin = unit(c(1, 1, 2, 2), "pt"),
												 
												 axis.title.y = element_blank())

