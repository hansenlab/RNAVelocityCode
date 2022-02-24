rm(list=ls())
source(here::here("scripts/utils.R"))


dynamical.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))
steady.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums(dynamical.o, gene = g, color_by = "fucci_time",
															 title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
fucci.p <- plot_grid(plotlist = top10.lp,
								nrow = 1, ncol = 10, label_size = 10, labels = NULL)


dynamical.o <- qread(here::here("data/Pancreas/dynamical.o.qs"))
steady.o <- qread(here::here("data/Pancreas/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "clusters", colors.v = metadata(dynamical.o)$clusters_colors,
															 title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
pancreas.p <- plot_grid(plotlist = top10.lp,
										 nrow = 1, ncol = 10, label_size = 10, labels = NULL)



dynamical.o <-qread(here::here("data/chromaffin/dynamical.o.qs"))
steady.o <-qread(here::here("data/chromaffin/steady.o.qs"))

colors.v <- levels(factor(dynamical.o$cell_type))
dynamical.o$cell_type <- fct_recode(factor(dynamical.o$cell_type), "SCPs" = "#0066FFFF", "chromaffin" = "#00FF66FF",
															"sympathoblasts" = "#CC00FFFF", "non-labelled" = "#CCFF00FF", "bridge" = "#FF0000FF")

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "cell_type", colors.v = colors.v, point.size = 5.01,
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
chromaffin.p <- plot_grid(plotlist = top10.lp,
												nrow = 1, ncol = 10, label_size = 10, labels = NULL)


dynamical.o <- qread(here::here("data/bonemarrow/dynamical.o.qs"))
steady.o <- qread(here::here("data/bonemarrow/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "clusters", colors.v = metadata(dynamical.o)$clusters_colors,
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
Bonemarrow.p <- plot_grid(plotlist = top10.lp,
													nrow = 1, ncol = 10, label_size = 10, labels = NULL)



dynamical.o <- qread(here::here("data/dentategyrus_lamanno/dynamical.o.qs"))
steady.o <- qread(here::here("data/dentategyrus_lamanno/steady.o.qs"))

dynamical.o$ClusterName <- factor(dynamical.o$ClusterName, levels = c("RadialGlia", "RadialGlia2", "ImmAstro", 
																													"GlialProg", "OPC", "nIPC", "Nbl1", "Nbl2",
																													"ImmGranule1", "ImmGranule2", "Granule", "CA", 
																													"CA1-Sub", "CA2-3-4"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "ClusterName", colors.v = metadata(dynamical.o)$cluster_colors, point.size = 2.01,
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
Dentategyrus_lamanno.p <- plot_grid(plotlist = top10.lp,
													nrow = 1, ncol = 10, label_size = 10, labels = NULL)



dynamical.o <- qread(here::here("data/dentategyrus/dynamical.o.qs"))
steady.o <- qread(here::here("data/dentategyrus/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "clusters", colors.v = metadata(dynamical.o)$clusters_colors, point.size = 3.01,
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
Dentategyrus_hochgerner.p <- plot_grid(plotlist = top10.lp,
																		nrow = 1, ncol = 10, label_size = 10, labels = NULL)




dynamical.o <- qread(here::here("data/erythroid/dynamical.o.qs"))
steady.o <- qread(here::here("data/erythroid/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "celltype", colors.v = metadata(dynamical.o)$celltype_colors, 
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
erythroid.p <- plot_grid(plotlist = top10.lp,
																			 nrow = 1, ncol = 10, label_size = 10, labels = NULL)




dynamical.o <- qread(here::here("data/forebrain/dynamical.o.qs"))
steady.o <- qread(here::here("data/forebrain/steady.o.qs"))
dynamical.o$Clusters <- factor(dynamical.o$Clusters)

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "Clusters", colors.v = brewer.pal(nlevels(factor(dynamical.o$Clusters)), "Set3"), 
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
forebrain.p <- plot_grid(plotlist = top10.lp,
												 nrow = 1, ncol = 10, label_size = 10, labels = NULL)



dynamical.o <- qread(here::here("data/gastrulation_e75/dynamical.o.qs"))
steady.o <- qread(here::here("data/gastrulation_e75/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "celltype", colors.v = names(table(dynamical.o$colour)),  point.size = 2.01,
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
gastrulation_e75.p <- plot_grid(plotlist = top10.lp,
												 nrow = 1, ncol = 10, label_size = 10, labels = NULL)



dynamical.o <- qread(here::here("data/pbmc68k/dynamical.o.qs"))
steady.o <- qread(here::here("data/pbmc68k/steady.o.qs"))

genes.df <- rowData(dynamical.o) %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

top10.lp <- lapply(rownames(genes.df)[seq_len(10)], function(g) {
	mums.p <- plot_dynamicalmums2(dynamical.o, rownm = g, color_by = "celltype", colors.v = brewer.pal(nlevels(factor(dynamical.o$celltype)), "Set3"),  point.size = 4.01,
																title = as.expression(substitute(italic(g)~" LL:"~ll, list(g=g, ll = sprintf("%.3f", rowData(dynamical.o[g, ])$fit_likelihood)))), color.name = "FUCCI") +
		geom_abline(slope = rowData(steady.o[g, ])$velocity_gamma, intercept = 0, linetype = "dashed") +
		theme(legend.position = "none")
	return(mums.p)
})
pbmc68k.p <- plot_grid(plotlist = top10.lp,
																nrow = 1, ncol = 10, label_size = 10, labels = NULL)


mp <- plot_grid( forebrain.p, chromaffin.p, fucci.p, Bonemarrow.p, pancreas.p,
								 Dentategyrus_hochgerner.p, Dentategyrus_lamanno.p,
								  erythroid.p,  gastrulation_e75.p, pbmc68k.p,
								 nrow = 10, ncol = 1, label_size = 10, labels = "auto")
save_plot(here::here("figs", "sfigs", "sfig.realdata_top10.pdf"), mp,
					base_height = 2 *1.2 , base_width = 2*1.4 , nrow = 10, ncol = 10, device = cairo_pdf)
