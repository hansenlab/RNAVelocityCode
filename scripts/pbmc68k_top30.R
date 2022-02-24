rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/pbmc68k/dynamical.o.qs"))

library(org.Hs.eg.db)
go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

cc.df <- rowData(sce.o)[rownames(sce.o) %in% go.v, ] %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))


### embeddings vectors
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)


vec.p <- plotEmbScat(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2", color_by = "celltype", labels.v = levels(factor(sce.o$celltype)),
						color.name = "Cell_type", colors.v = brewer.pal(nlevels(factor(sce.o$celltype)), "Set3")) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.015, "inches")), inherit.aes = FALSE, size = 0.1)


save_plot(here::here("figs", "pbmc68k", str_c("vec.pdf")), vec.p,
					base_height = 2, base_width = 2*1.4, nrow = 1, ncol = 1*1.5, device = cairo_pdf)


### plot top 30 genes
top30.df <- rowData(sce.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top30.lp <- lapply(seq_len(30), function(i) {
	g <- rownames(top30.df)[i]
	p <- plot_dynamicalmums2(sce.o, gene = g, color_by = "celltype", colors.v = brewer.pal(nlevels(factor(sce.o$celltype)), "Set3"), point.size = 4.01,
													 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),")")) + 
		theme(legend.position = "none")
	return(p)
	
})

mp <- wrap_plots(top30.lp, ncol = 5)

save_plot(here::here("figs", "pbmc68k", str_c("top30.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 6, ncol = 5, device = cairo_pdf)



