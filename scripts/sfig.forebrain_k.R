rm(list=ls())
source(here::here("scripts/utils.R"))


scv <- reticulate::import("scvelo")


sce.o <- qread(here::here("data/forebrain/steady.o.qs"))
sce.o$Clusters <- factor(sce.o$Clusters)

embedded <- embedVelocity(reducedDim(sce.o, "PCA")[, 1:2], sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "PCA")[, 1:2], embedded, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
steady_k30.p <- plotEmbScat(sce.o, "PCA", title = str_c("Forebrain (steady k=30)"), x_lab = "PCA-1", y_lab = "PCA-2", color_by = "Clusters", labels.v = levels(sce.o$Clusters),
																			 color.name = "Cell_type" ) + 
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.3) +
	.theme_noframe +
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/forebrain/steady2.o.qs"))
sce.o$Clusters <- factor(sce.o$Clusters)

embedded <- embedVelocity(reducedDim(sce.o, "PCA")[, 1:2], sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "PCA")[, 1:2], embedded, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
steady_k550.p <- plotEmbScat(sce.o, "PCA", title = str_c("Forebrain (steady k=550)"), x_lab = "PCA-1", y_lab = "PCA-2", color_by = "Clusters", labels.v = levels(sce.o$Clusters),
														color.name = "Cell_type" ) + 
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.3) +
	.theme_noframe +
	theme(legend.position = "none",
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))



sce.o <- qread(here::here("data/forebrain/dynamical.o.qs"))
sce.o$Clusters <- factor(sce.o$Clusters)

embedded <- embedVelocity(reducedDim(sce.o, "PCA")[, 1:2], sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "PCA")[, 1:2], embedded, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
dynamical_k30.p <- plotEmbScat(sce.o, "PCA", title = str_c("Forebrain (dynamical k=30)"), x_lab = "PCA-1", y_lab = "PCA-2", color_by = "Clusters", labels.v = levels(sce.o$Clusters),
															 color.name = "Cell_type" ) + 
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.3) +
	.theme_noframe +
	theme(legend.position = "none",
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))


sce.o <- qread(here::here("data/forebrain/dynamical2.o.qs"))
sce.o$Clusters <- factor(sce.o$Clusters)

embedded <- embedVelocity(reducedDim(sce.o, "PCA")[, 1:2], sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "PCA")[, 1:2], embedded, resolution = 30)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
dynamical_k550.p <- plotEmbScat(sce.o, "PCA", title = str_c("Forebrain (dynamical k=550)"), x_lab = "PCA-1", y_lab = "PCA-2", color_by = "Clusters", labels.v = levels(sce.o$Clusters),
																color.name = "Cell_type" ) + 
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.03, "inches")), inherit.aes = FALSE, size = 0.3) +
	.theme_noframe +
	theme(legend.position = "none",
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))



mp <- plot_grid( steady_k30.p, steady_k550.p,
								 dynamical_k30.p, dynamical_k550.p,
								 byrow = TRUE, 
								 nrow = 2, ncol = 2, label_size = 10, labels = "auto")
save_plot(here::here("figs", "sfigs", "sfig.forebrain_k.pdf"), mp,
					base_height = 2 *1.2 , base_width = 2*1.4 , nrow = 2 * 1.5, ncol = 2 * 1.5, device = cairo_pdf)



sce.o <- qread(here::here("data/forebrain/sce.o.qs"))
steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = 550L)))
reducedDims(steady.o) <- reducedDims(sce.o)
rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))
metadata(steady.o) <- c(metadata(steady.o), metadata(sce.o))
sce.o <- steady.o

adata <- zellkonverter::SCE2AnnData(sce.o)
scv$pp$neighbors(adata, use_rep = 'PCA', n_neighbors = 150L)
scv$tl$velocity_graph(adata)
sce.o <- zellkonverter::AnnData2SCE(adata)

embedded <- embedVelocity(reducedDim(sce.o, "PCA"), sce.o)
sce.o$Clusters <- factor(sce.o$Clusters)
forebrain2.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "PCA", color_by = "Clusters", color.alpha = 0.7, grid.resolution = 20,
																	stream.L = 4, stream.min.L = 1, stream.res = 10, stream.width = 0.4,
																	arrow.angle = 18, arrow.length = 0.35) +
	scale_color_manual(values = brewer.pal(nlevels(factor(sce.o$Clusters)), "Set3"), name = "Cell type",
										 labels = c("Radial glia", "Radial glia",
										 					 "Neuroblast", "Neuroblast",
										 					 "Immature neuron", #"Immature neuron",
										 					 "Neuron", "Neuron"), limits = levels(sce.o$Clusters)) +
	labs( y = "PCA-2", x = "PCA-1") +
	ggtitle(str_c("Forebrain (steady, bigK)")) +
	.theme_noframe + 
	theme(legend.position = c(0, 0),
				legend.justification = c(0, 0),
				legend.key.size = unit(6, "pt"))
