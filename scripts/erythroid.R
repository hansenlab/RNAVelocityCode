rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/erythroid/dynamical.o.qs"))


### plot top 30 genes
top30.df <- rowData(sce.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top30.lp <- lapply(seq_len(30), function(i) {
	g <- rownames(top30.df)[i]
	p <- plot_dynamicalmums2(sce.o, rownm = g, color_by = "celltype", colors.v = metadata(sce.o)$celltype_colors, point.size = 2.01,
													 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
		theme(legend.position = "none")
	return(p)
	
})

mp <- wrap_plots(top30.lp, ncol = 5)

save_plot(here::here("figs", "erythroid", str_c("top30.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 5, device = cairo_pdf)



### 
### remove Erythroid1
idx <- which(sce.o$celltype != "Erythroid1")
new.o <- sce.o[, idx]
# new.o$celltype <- droplevels(new.o$celltype)
reducedDims(new.o) <- list()

adata <- zellkonverter::SCE2AnnData(new.o)
scv <- reticulate::import("scvelo")
scv$pp$pca(adata )
scv$pp$neighbors(adata, use_rep = "X_pca")
scv$pp$moments(adata, use_rep = "X_pca")
scv$tl$umap(adata)
scv$tl$differential_kinetic_test(adata, groupby='celltype')
new.o <- zellkonverter::AnnData2SCE(adata)

reducedDimNames(new.o) <- c("PCA", "UMAP")

# dynamical.o <- scvelo(new.o, assay.X="spliced", mode =  "dynamical", use.dimred = "PCA")
# reducedDims(dynamical.o) <- reducedDims(new.o)
# rowData(dynamical.o) <- cbind(rowData(dynamical.o), rowData(new.o))
# colData(dynamical.o) <- cbind(colData(dynamical.o), colData(new.o))
# metadata(dynamical.o) <- c(metadata(dynamical.o), metadata(new.o))

sce_old.o <- sce.o
sce.o <- new.o



### plot top 30 genes
top30.df <- rowData(sce.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top30.lp <- lapply(seq_len(30), function(i) {
	g <- rownames(top30.df)[i]
	p <- plot_dynamicalmums2(sce.o, rownm = g, color_by = "celltype", colors.v = metadata(sce_old.o)$celltype_colors, point.size = 2.01,
													 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
		theme(legend.position = "none")
	return(p)
	
})

mp <- wrap_plots(top30.lp, ncol = 5)

save_plot(here::here("figs", "erythroid", str_c("top30_2.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 5, device = cairo_pdf)

qsave(sce.o, file = here::here("data/erythroid/sce_subset.o.qs"))



## figure
p1 <- plotEmbScat(sce_old.o, dimred = "UMAP", color_by = "celltype",
						colors.v = metadata(sce_old.o)$celltype_colors,
						labels.v = levels(sce_old.o$celltype),
						color.name = "Celltype", x_lab = "UMAP-1", y_lab = "UMAP-2", 
						title = "Gastrulation erythroid maturation", point.size = 4.01) +
	.theme_noframe + theme(legend.position = c(0.1, 0.2),
												 legend.justification = c(0, 0))
p2 <- plotEmbScat(sce.o, dimred = "UMAP", color_by = "celltype",
									colors.v = metadata(sce.o)$celltype_colors,
									labels.v = levels(sce.o$celltype),
									color.name = "Celltype", x_lab = "UMAP-1", y_lab = "UMAP-2", 
									title = "Without Erythroid1", point.size = 4.01) +
	.theme_noframe + theme(legend.position = "none",
												 legend.justification = c(0, 0))

g <- "Smim1"
p3 <- plot_dynamicalmums2(sce_old.o, rownm = g, color_by = "celltype", colors.v = metadata(sce_old.o)$celltype_colors, point.size = 4.01,
													title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
	theme(legend.position = "none")

p4 <- plot_dynamicalmums2(sce.o, rownm = g, color_by = "celltype", colors.v = metadata(sce_old.o)$celltype_colors[-3], point.size = 4.01,
													title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
	theme(legend.position = "none")

mp <- plot_grid(p1, p2, p3, p4,
								nrow = 2, ncol = 2, label_size = 10, labels = c("a", "b"))

save_plot(here::here("figs", "erythroid", "erythroid_smim1.pdf"), mp,
					base_height = 2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2*1.2, ncol = 2, device = cairo_pdf)

