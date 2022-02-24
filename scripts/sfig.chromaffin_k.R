rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/chromaffin/sce.o.qs"))
colors.v <- levels(factor(sce.o$cell_type))
sce.o$cell_type <- fct_recode(factor(sce.o$cell_type), "SCPs" = "#0066FFFF", "chromaffin" = "#00FF66FF",
															"sympathoblasts" = "#CC00FFFF", "non-labelled" = "#CCFF00FF", "bridge" = "#FF0000FF")

library(velociraptor)
steady.lp <- lapply(c(5L, 10L, 20L, 30L), function(k) {
	steady.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = k)))
	reducedDims(steady.o) <- reducedDims(sce.o)
	rowData(steady.o) <- cbind(rowData(steady.o), rowData(sce.o))
	colData(steady.o) <- cbind(colData(steady.o), colData(sce.o))
	
	sce.o <- steady.o 
	rotation.m <- attr(reducedDim(sce.o, "PCA"), "rotation")[, 1:2]
	reducedDim(sce.o, "PCA") <- reducedDim(sce.o, "PCA")[, 1:2]
	sce.o <- logNormCounts(sce.o, assay.type= "spliced")
	data.m <- assay(sce.o, "logcounts")[match(rownames(rotation.m), rownames(sce.o)), ]
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	common_gene.v <- intersect(rownames(data.m), rownames(sce.o)[which(!rowAnyNAs(assay(sce.o, "velocity")))])
	
	### add velocity
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "spliced")
	splusv.m[which(splusv.m < 0)] <- 0
	data.m <- log2(splusv.m[match(rownames(rotation.m), rownames(sce.o)), ] + 1)
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
	
	embedded <- new_pca.m - reducedDim(sce.o, "PCA")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "PCA"))
	
	grid.df <- gridVectors(reducedDim(sce.o, "PCA"), embedded, resolution = 25)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	
	chromaffin_steady.p <- plotEmbScat(sce.o, "PCA", title = str_c("Steady (k=", k, ")"), x_lab = "PCA-1", y_lab = "PCA-2", color_by = "cell_type", labels.v = levels(sce.o$cell_type),
																				color.name = "Cell_type", colors.v = colors.v, point.size = 5.01, point.alpha = 0.7) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
		.theme_noframe +
		theme(legend.position = c(1, 1), legend.justification = c(1, 1))
	return(chromaffin_steady.p)
})

dynamical.lp <- lapply(c(5L, 10L, 20L, 30L), function(k) {
	dynamical.o <- scvelo(sce.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA", scvelo.params = list(moments = list(n_neighbors = k)))
	reducedDims(dynamical.o) <- reducedDims(sce.o)
	rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(sce.o))
	colData(dynamical.o) <- cbind(colData(dynamical.o), colData(sce.o))
	
	sce.o <- dynamical.o
	rotation.m <- attr(reducedDim(sce.o, "PCA"), "rotation")[, 1:2]
	reducedDim(sce.o, "PCA") <- reducedDim(sce.o, "PCA")[, 1:2]
	sce.o <- logNormCounts(sce.o, assay.type= "spliced")
	data.m <- assay(sce.o, "logcounts")[match(rownames(rotation.m), rownames(sce.o)), ]
	row_m.v <- rowMeans(data.m)
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	common_gene.v <- intersect(rownames(data.m), rownames(sce.o)[which(!rowAnyNAs(assay(sce.o, "velocity")))])
	
	### add velocity
	splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "spliced")
	splusv.m[which(splusv.m < 0)] <- 0
	data.m <- log2(splusv.m[match(rownames(rotation.m), rownames(sce.o)), ] + 1)
	
	data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
	idx <- which(!rowAnyNAs(as.matrix(data.m)))
	new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])
	
	embedded <- new_pca.m - reducedDim(sce.o, "PCA")
	dimnames(embedded) <- dimnames(reducedDim(sce.o, "PCA"))
	
	grid.df <- gridVectors(reducedDim(sce.o, "PCA"), embedded, resolution = 25)
	colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
	
	chromaffin_dynamical.p <- plotEmbScat(sce.o, "PCA", title = str_c("Dynamical (k=", k, ")"), x_lab = "PCA-1", y_lab = "PCA-2", color_by = "cell_type", labels.v = levels(sce.o$cell_type),
																		 color.name = "Cell_type", colors.v = colors.v, point.size = 5.01, point.alpha = 0.7) + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.04, "inches")), inherit.aes = FALSE, size = 0.2)  +
		.theme_noframe +
		theme(legend.position = c(1, 1), legend.justification = c(1, 1))
	return(chromaffin_dynamical.p)
})



mp <- plot_grid(plotlist = c(steady.lp, dynamical.lp),
								nrow = 4, ncol = 2, label_size = 10, labels = NULL, byrow = FALSE)

save_plot(here::here("figs", "sfigs", "sfig.chromaffin_k.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 5, ncol = 3, device = cairo_pdf)









