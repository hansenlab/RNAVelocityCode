rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))

go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

cc.df <- rowData(sce.o)[rownames(sce.o) %in% go.v, ] %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(velocity_r2))

"TOP2A" %in% rownames(cc.df)

### embeddings vectors
embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)


### TOP2A
g <- "TOP2A"
us.p <- plot_us(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
mums.p <- plot_steadymums(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
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
velo_ms.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(sce.o, "velocity")[g, ], title = "Velocity ~ Ms", x_lab = "Ms", y_lab = "Velocity", color_by = "fucci_time",color.name = "FUCCI")

mp <- wrap_plots(list(mums.p, exp_ms.p, exp_mu.p,
											us.p, velo.p, vec.p,
											loess_s.p, loess_ms.p, loess_velo.p,
											loess_u.p, loess_mu.p, velo_ms.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c(g, ".pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)




### CDK1
g <- "CDK1"
us.p <- plot_us(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
mums.p <- plot_steadymums(sce.o, gene = g, color_by = "fucci_time", title = g, color.name = "FUCCI")
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
velo_ms.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(sce.o, "velocity")[g, ], title = "Velocity ~ Ms", x_lab = "Ms", y_lab = "Velocity", color_by = "fucci_time",color.name = "FUCCI")

mp <- wrap_plots(list(mums.p, exp_ms.p, exp_mu.p,
											us.p, velo.p, vec.p,
											loess_s.p, loess_ms.p, loess_velo.p,
											loess_u.p, loess_mu.p, velo_ms.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c(g, ".pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 3, device = cairo_pdf)






### cdk1 top2a
u.p <- plotxy(assay(sce.o, "unspliced")["TOP2A", ], assay(sce.o, "unspliced")["CDK1", ], title = "Unspliced", x_lab = "TOP2A", y_lab = "CDK1", color_by = "fucci_time",color.name = "FUCCI")
s.p <- plotxy(assay(sce.o, "spliced")["TOP2A", ], assay(sce.o, "spliced")["CDK1", ], title = "Spliced", x_lab = "TOP2A", y_lab = "CDK1", color_by = "fucci_time",color.name = "FUCCI")
mu.p <- plotxy(assay(sce.o, "Mu")["TOP2A", ], assay(sce.o, "Mu")["CDK1", ], title = "Mu", x_lab = "TOP2A", y_lab = "CDK1", color_by = "fucci_time",color.name = "FUCCI")
ms.p <- plotxy(assay(sce.o, "Ms")["TOP2A", ], assay(sce.o, "Ms")["CDK1", ], title = "Ms", x_lab = "TOP2A", y_lab = "CDK1", color_by = "fucci_time",color.name = "FUCCI")
t.v <- runif(100, 0, 2 * pi)
sim.p <- plotxy(cos(t.v), cos(t.v + 0.5), title = "Simulation of two genes in cell cycle", x_lab = "cos(t)", y_lab = "cos(t + 0.5)")


mp <- wrap_plots(list(mu.p, ms.p, sim.p,
											u.p, s.p), ncol = 3)

save_plot(here::here("figs", "FUCCI", str_c( "CDK1_TOP2A.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 2, ncol = 3, device = cairo_pdf)









