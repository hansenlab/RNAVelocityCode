rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

library(org.Hs.eg.db)
go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

cc.df <- rowData(sce.o)[rownames(sce.o) %in% go.v, ] %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))

"TOP2A" %in% rownames(cc.df)

rowData(sce.o)["CDK1", ]

### embeddings vectors
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)



### plot
g <- "TOP2A"
us.p <- plot_us(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
mums.p <- plot_dynamicalmums(sce.o, gene = g, color_by = "fucci_time",
														 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),
														 							", v.gene:", ifelse(rowData(sce.o[g, ])$velocity_genes, "T", "F"),")"), color.name = "FUCCI")
exp_ms.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Ms")[g, ], color.name = "Ms", title = "Ms", x_lab = "UMAP-1", y_lab = "UMAP-2")
exp_mu.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Mu")[g, ], color.name = "Mu", title = "Mu", x_lab = "UMAP-1", y_lab = "UMAP-2")
velo.p <- plotEmbVelo(sce.o, dimred = "UMAP", color.v = assay(sce.o, "velocity")[g, ], color.name = "Velocity", title = "Velocity", x_lab = "UMAP-1", y_lab = "UMAP-2")
vec.p <- plotEmbFucci(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)
loess_ms.p <- plotFucciLoess(sce.o, assay(sce.o, "Ms")[g, ], y_lab = "Ms", title = "Ms ~ FUCCI")
loess_s.p <- plotFucciLoess(sce.o, assay(sce.o, "spliced")[g, ], y_lab = "spliced", title = "spliced ~ FUCCI")
loess_u.p <- plotFucciLoess(sce.o, assay(sce.o, "unspliced")[g, ], y_lab = "unspliced", title = "unspliced ~ FUCCI")
loess_mu.p <- plotFucciLoess(sce.o, assay(sce.o, "Mu")[g, ], y_lab = "Mu", title = "Mu ~ FUCCI")
loess_velo.p <-  plotFucciLoess(sce.o, assay(sce.o, "velocity")[g, ], y_lab = "Velocity", title = "Velocity ~ FUCCI")
vt_ft.p <- plotxy(sce.o$fucci_time, assay(sce.o, "fit_t")[g, ], title = "Velocity latent time ~ FUCCI", x_lab = "FUCCI pseudotime", y_lab = "Velocity latent time", color_by = "fucci_time",color.name = "FUCCI")

mp <- wrap_plots(list(mums.p, exp_ms.p, exp_mu.p,
											us.p, velo.p, vec.p,
											loess_s.p, loess_ms.p, loess_velo.p,
											loess_u.p, loess_mu.p, vt_ft.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c(g, "_dynamical.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)



g <- "CDK1"
us.p <- plot_us(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
mums.p <- plot_dynamicalmums(sce.o, gene = g, color_by = "fucci_time",
														 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),
														 							", v.gene:", ifelse(rowData(sce.o[g, ])$velocity_genes, "T", "F"),")"), color.name = "FUCCI")
exp_ms.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Ms")[g, ], color.name = "Ms", title = "Ms", x_lab = "UMAP-1", y_lab = "UMAP-2")
exp_mu.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Mu")[g, ], color.name = "Mu", title = "Mu", x_lab = "UMAP-1", y_lab = "UMAP-2")
velo.p <- plotEmbVelo(sce.o, dimred = "UMAP", color.v = assay(sce.o, "velocity")[g, ], color.name = "Velocity", title = "Velocity", x_lab = "UMAP-1", y_lab = "UMAP-2")
vec.p <- plotEmbFucci(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)
loess_ms.p <- plotFucciLoess(sce.o, assay(sce.o, "Ms")[g, ], y_lab = "Ms", title = "Ms ~ FUCCI")
loess_s.p <- plotFucciLoess(sce.o, assay(sce.o, "spliced")[g, ], y_lab = "spliced", title = "spliced ~ FUCCI")
loess_u.p <- plotFucciLoess(sce.o, assay(sce.o, "unspliced")[g, ], y_lab = "unspliced", title = "unspliced ~ FUCCI")
loess_mu.p <- plotFucciLoess(sce.o, assay(sce.o, "Mu")[g, ], y_lab = "Mu", title = "Mu ~ FUCCI")
loess_velo.p <-  plotFucciLoess(sce.o, assay(sce.o, "velocity")[g, ], y_lab = "Velocity", title = "Velocity ~ FUCCI")
vt_ft.p <- plotxy(sce.o$fucci_time, assay(sce.o, "fit_t")[g, ], title = "Velocity latent time ~ FUCCI", x_lab = "FUCCI pseudotime", y_lab = "Velocity latent time", color_by = "fucci_time",color.name = "FUCCI")

mp <- wrap_plots(list(mums.p, exp_ms.p, exp_mu.p,
											us.p, velo.p, vec.p,
											loess_s.p, loess_ms.p, loess_velo.p,
											loess_u.p, loess_mu.p, vt_ft.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c(g, "_dynamical.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)

g <- "UBE2C"
us.p <- plot_us(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
mums.p <- plot_dynamicalmums(sce.o, gene = g, color_by = "fucci_time",
														 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),
														 							", v.gene:", ifelse(rowData(sce.o[g, ])$velocity_genes, "T", "F"),")"), color.name = "FUCCI")
exp_ms.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Ms")[g, ], color.name = "Ms", title = "Ms", x_lab = "UMAP-1", y_lab = "UMAP-2")
exp_mu.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Mu")[g, ], color.name = "Mu", title = "Mu", x_lab = "UMAP-1", y_lab = "UMAP-2")
velo.p <- plotEmbVelo(sce.o, dimred = "UMAP", color.v = assay(sce.o, "velocity")[g, ], color.name = "Velocity", title = "Velocity", x_lab = "UMAP-1", y_lab = "UMAP-2")
vec.p <- plotEmbFucci(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)
loess_ms.p <- plotFucciLoess(sce.o, assay(sce.o, "Ms")[g, ], y_lab = "Ms", title = "Ms ~ FUCCI")
loess_s.p <- plotFucciLoess(sce.o, assay(sce.o, "spliced")[g, ], y_lab = "spliced", title = "spliced ~ FUCCI")
loess_u.p <- plotFucciLoess(sce.o, assay(sce.o, "unspliced")[g, ], y_lab = "unspliced", title = "unspliced ~ FUCCI")
loess_mu.p <- plotFucciLoess(sce.o, assay(sce.o, "Mu")[g, ], y_lab = "Mu", title = "Mu ~ FUCCI")
loess_velo.p <-  plotFucciLoess(sce.o, assay(sce.o, "velocity")[g, ], y_lab = "Velocity", title = "Velocity ~ FUCCI")
vt_ft.p <- plotxy(sce.o$fucci_time, assay(sce.o, "fit_t")[g, ], title = "Velocity latent time ~ FUCCI", x_lab = "FUCCI pseudotime", y_lab = "Velocity latent time", color_by = "fucci_time",color.name = "FUCCI")

mp <- wrap_plots(list(mums.p, exp_ms.p, exp_mu.p,
											us.p, velo.p, vec.p,
											loess_s.p, loess_ms.p, loess_velo.p,
											loess_u.p, loess_mu.p, vt_ft.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c(g, "_dynamical.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)



for (i in seq_len(10)) {
	g <- rownames(cc.df)[i]
	us.p <- plot_us(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
	mums.p <- plot_dynamicalmums(sce.o, gene = g, color_by = "fucci_time",
															 title = str_c(g, " (LL:", sprintf("%.3f", rowData(sce.o[g, ])$fit_likelihood),
															 							", v.gene:", ifelse(rowData(sce.o[g, ])$velocity_genes, "T", "F"),")"), color.name = "FUCCI")
	exp_ms.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Ms")[g, ], color.name = "Ms", title = "Ms", x_lab = "UMAP-1", y_lab = "UMAP-2")
	exp_mu.p <- plotEmbExp(sce.o, dimred = "UMAP", color.v = assay(sce.o, "Mu")[g, ], color.name = "Mu", title = "Mu", x_lab = "UMAP-1", y_lab = "UMAP-2")
	velo.p <- plotEmbVelo(sce.o, dimred = "UMAP", color.v = assay(sce.o, "velocity")[g, ], color.name = "Velocity", title = "Velocity", x_lab = "UMAP-1", y_lab = "UMAP-2")
	vec.p <- plotEmbFucci(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
		geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)
	loess_ms.p <- plotFucciLoess(sce.o, assay(sce.o, "Ms")[g, ], y_lab = "Ms", title = "Ms ~ FUCCI")
	loess_s.p <- plotFucciLoess(sce.o, assay(sce.o, "spliced")[g, ], y_lab = "spliced", title = "spliced ~ FUCCI")
	loess_u.p <- plotFucciLoess(sce.o, assay(sce.o, "unspliced")[g, ], y_lab = "unspliced", title = "unspliced ~ FUCCI")
	loess_mu.p <- plotFucciLoess(sce.o, assay(sce.o, "Mu")[g, ], y_lab = "Mu", title = "Mu ~ FUCCI")
	loess_velo.p <-  plotFucciLoess(sce.o, assay(sce.o, "velocity")[g, ], y_lab = "Velocity", title = "Velocity ~ FUCCI")
	vt_ft.p <- plotxy(sce.o$fucci_time, assay(sce.o, "fit_t")[g, ], title = "Velocity latent time ~ FUCCI", x_lab = "FUCCI pseudotime", y_lab = "Velocity latent time", color_by = "fucci_time",color.name = "FUCCI")
	
	mp <- wrap_plots(list(mums.p, exp_ms.p, exp_mu.p,
												us.p, velo.p, vec.p,
												loess_s.p, loess_ms.p, loess_velo.p,
												loess_u.p, loess_mu.p, vt_ft.p), ncol = 3)
	
	save_plot(here::here("figs", "FUCCI", str_c(i, "_", g, "_dynamical.pdf")), mp,
						base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)
	
}


