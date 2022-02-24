rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))
umap.o <- qread(here::here("data/Pancreas/umap.o.qs"))

### plot
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)

grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded, resolution = 40)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p1 <- plotEmbScat(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
									color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


### remove Ngn3 high EP
idx <- which(sce.o$clusters != "Ngn3 high EP")
new.o <- sce.o[, idx]
new.o$clusters <- droplevels(new.o$clusters)

adata <- zellkonverter::SCE2AnnData(new.o)
scv <- reticulate::import("scvelo")
scv$pp$neighbors(adata, use_rep = "PCA")
scv$tl$velocity_graph(adata)
new.o <- zellkonverter::AnnData2SCE(adata)


embedded <- embedVelocity(reducedDim(new.o, "UMAP"), new.o)


grid.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded, resolution = 40)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p2 <- plotEmbScat(new.o, "UMAP", title = "Vector fields subset", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
									color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)



### remove Pre-endocrine
idx <- which(sce.o$clusters != "Pre-endocrine")
new.o <- sce.o[, idx]
new.o$clusters <- droplevels(new.o$clusters)

adata <- zellkonverter::SCE2AnnData(new.o)
scv <- reticulate::import("scvelo")
scv$pp$neighbors(adata, use_rep = "PCA")
scv$tl$velocity_graph(adata)
new.o <- zellkonverter::AnnData2SCE(adata)


embedded <- embedVelocity(reducedDim(new.o, "UMAP"), new.o)

grid.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded, resolution = 40)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p3 <- plotEmbScat(new.o, "UMAP", title = "Vector fields subset", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
									color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


### remove Beta
idx <- which(sce.o$clusters != "Beta")
new.o <- sce.o[, idx]
new.o$clusters <- droplevels(new.o$clusters)

adata <- zellkonverter::SCE2AnnData(new.o)
scv <- reticulate::import("scvelo")
scv$pp$neighbors(adata, use_rep = "PCA")
scv$tl$velocity_graph(adata)
new.o <- zellkonverter::AnnData2SCE(adata)


embedded <- embedVelocity(reducedDim(new.o, "UMAP"), new.o)


grid.df <- gridVectors(reducedDim(new.o, "UMAP"), embedded, resolution = 40)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p4 <- plotEmbScat(new.o, "UMAP", title = "Vector fields subset", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
									color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)







mp <- wrap_plots(list(p1, p2, p3, p4), ncol = 2)

save_plot(here::here("figs", "pancreas", str_c("umap_subset.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 1.5 * 2, ncol = 2 * 1.5, device = cairo_pdf)







