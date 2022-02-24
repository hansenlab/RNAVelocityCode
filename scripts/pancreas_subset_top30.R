rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))

### remove Pre-endocrine
idx <- which(sce.o$clusters != "Pre-endocrine")
new.o <- sce.o[, idx]
new.o$clusters <- droplevels(new.o$clusters)
reducedDims(new.o) <- list()

adata <- zellkonverter::SCE2AnnData(new.o)
scv <- reticulate::import("scvelo")
scv$pp$pca(adata )
scv$pp$neighbors(adata, use_rep = "X_pca")
scv$pp$moments(adata, use_rep = "X_pca")
scv$tl$umap(adata)
scv$tl$differential_kinetic_test(adata, var_names="Gnas", groupby='clusters')
new.o <- zellkonverter::AnnData2SCE(adata)

kwargs = dict(s=60, linewidth=3, add_linfit=True, frameon='artist')
scv$pl$scatter(adata, basis= c("umap", "Gnas"), add_outline='fit_diff_kinetics', dpi=80, linewidth=3, add_linfit=TRUE)



reducedDimNames(new.o) <- c("PCA", "UMAP")

dynamical.o <- scvelo(new.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
reducedDims(dynamical.o) <- reducedDims(new.o)
rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(new.o))
colData(dynamical.o) <- cbind(colData(dynamical.o), colData(new.o))
metadata(dynamical.o) <- c(metadata(dynamical.o), metadata(new.o))

sce_old.o <- sce.o
sce.o <- dynamical.o

### embeddings vectors
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)

update_geom_defaults("point",list(size = 0.2))
vec.p <- plotVelocityStream(sce.o, embedded[, 1:2], use.dimred = "UMAP", color_by = "clusters", color.alpha = 0.7, grid.resolution = 20,
														stream.L = 4, stream.min.L = 0.1, stream.res = 8, stream.width = 0.2,
														arrow.angle = 18, arrow.length = 0.2) +
	scale_color_manual(values = metadata(sce.o)$clusters_colors, name = "Cell_type", labels = levels(sce.o$clusters),
										 limits = levels(sce.o$clusters)) +
	labs( y = "UMAP-2", x = "UMAP-1") +
	ggtitle("Pancreas (dynamical model)") +
	.theme_noframe + theme(legend.position = "right")


# plotEmbScat(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "clusters", labels.v = levels(sce.o$clusters),
# 									 color.name = "Cell_type", colors.v = metadata(sce.o)$clusters_colors) + 
# geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


save_plot(here::here("figs", "pancreas", str_c("vec2.pdf")), vec.p,
					base_height = 2 / 1.4 * 1.2, base_width = 2*1.4 / 1.4, nrow = 1, ncol = 1 * 1.3, device = cairo_pdf)


### plot top 30 genes
top30.df <- rowData(sce.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top30.lp <- lapply(seq_len(30), function(i) {
	g <- rownames(top30.df)[i]
	p <- plot_dynamicalmums2(sce.o, rownm = g, color_by = "clusters", colors.v = metadata(sce.o)$clusters_colors[-4], point.size = 2.01,
													 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
		theme(legend.position = "none")
	return(p)
	
})

mp <- wrap_plots(top30.lp, ncol = 5)

save_plot(here::here("figs", "pancreas", str_c("top30_2.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 5, device = cairo_pdf)

qsave(sce.o, file = here::here("data/Pancreas/sce_subset.o.qs"))






