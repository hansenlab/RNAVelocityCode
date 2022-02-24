rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

top2a_t_fucci.p <- plotFucciLoess(sce.o, assay(sce.o, "fit_t")["TOP2A", ], y_lab = "Velocity latent time",
												 title = expression(italic(TOP2A)~' latent t. over FUCCI'), 
												 x.pos = 0.6, y.pos = 0.9) +
	theme(legend.position = "none")

cdk1_t_fucci.p <- plotFucciLoess(sce.o, assay(sce.o, "fit_t")["CDK1", ], y_lab = "Velocity latent time",
																	title = expression(italic(CDK1)~' latent t. over FUCCI'), 
																	x.pos = 0.4, y.pos = 0.5) +
	theme(legend.position = c(1, 0.03),
				legend.justification = c(1, 0))

ms_com.p <- plotxy(assay(sce.o, "Ms")["TOP2A", ], assay(sce.o, "Ms")["CDK1", ],
									 title = "Comparisons of Ms", x_lab = "TOP2A", y_lab = "CDK1", color_by = "fucci_time",color.name = "FUCCI") +
	theme(legend.position = "none")

mu_com.p <- plotxy(assay(sce.o, "Mu")["TOP2A", ], assay(sce.o, "Mu")["CDK1", ],
									 title = "Comparisons of Mu", x_lab = "TOP2A", y_lab = "CDK1", color_by = "fucci_time",color.name = "FUCCI") +
	theme(legend.position = "none")


t.v <- runif(100, 0, 2 * pi)
ms_sim.p <- plotxy(cos(t.v), cos(t.v + 0.5), title = "Sim.", x_lab = "cos(t)", y_lab = "cos(t + 0.5)", point.size = 6.01) +
	theme_half_open(6) +
	theme(axis.ticks = element_blank(),
			axis.text.x = element_blank(),
			axis.text.y = element_blank(),
			plot.title = element_text(face = "plain", size = 6, hjust = 0.5))

ms_com.p <- ggdraw(ms_com.p) +
	draw_plot(ms_sim.p, 0.2, 0.55, .35, .35, hjust = 0, vjust = 0, halign = 0, valign = 0)




top2a_ms_t.p <- plotTLoess(x = assay(sce.o, "fit_t")["TOP2A", ], y = assay(sce.o, "Ms")["TOP2A", ],
					 color_var = sce.o$fucci_time, y_lab = "Ms", title = expression(italic(TOP2A)~' Ms over t'),
					 scale_color = scale_color_gradientn(name = "FUCCI", limits = range(sce.o$fucci_time), 
					 																										colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"))) +
	theme(legend.position = "none")

cdk1_ms_t.p <- plotTLoess(x = assay(sce.o, "fit_t")["CDK1", ], y = assay(sce.o, "Ms")["CDK1", ],
													 color_var = sce.o$fucci_time, y_lab = "Ms", title = expression(italic(CDK1)~' Ms over t'),
													 scale_color = scale_color_gradientn(name = "FUCCI", limits = range(sce.o$fucci_time), 
													 																		colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")),
													x.pos = 0.2, y.pos = 0.8) +
	theme(legend.position = "none")

top2a_mu_t.p <- plotTLoess(x = assay(sce.o, "fit_t")["TOP2A", ], y = assay(sce.o, "Mu")["TOP2A", ],
													 color_var = sce.o$fucci_time, y_lab = "Mu", title = expression(italic(TOP2A)~' Mu over t'),
													 scale_color = scale_color_gradientn(name = "FUCCI", limits = range(sce.o$fucci_time), 
													 																		colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"))) +
	theme(legend.position = "none")

cdk1_mu_t.p <- plotTLoess(x = assay(sce.o, "fit_t")["CDK1", ], y = assay(sce.o, "Mu")["CDK1", ],
													color_var = sce.o$fucci_time, y_lab = "Mu", title = expression(italic(CDK1)~' Mu over t'),
													scale_color = scale_color_gradientn(name = "FUCCI", limits = range(sce.o$fucci_time), 
																															colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")),
													x.pos = 0.2, y.pos = 0.8) +
	theme(legend.position = "none")


t_com.p <- plotxy(assay(sce.o, "fit_t")["TOP2A", ], assay(sce.o, "fit_t")["CDK1", ],
									 title = "Comparisons of latent time", x_lab = expression(italic(TOP2A)~' latent time'),
									y_lab = expression(italic(CDK1)~' latent time'), color_by = "fucci_time",color.name = "FUCCI") +
	theme(legend.position = "none")


mp <- plot_grid(top2a_t_fucci.p, cdk1_t_fucci.p, t_com.p, 
								top2a_ms_t.p, cdk1_ms_t.p, ms_com.p,
								top2a_mu_t.p, cdk1_mu_t.p, mu_com.p,
								nrow = 3, ncol = 3, label_size = 10, labels = c("auto"))

save_plot(here::here("figs", "main", "main.sync.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 3, ncol = 3, device = cairo_pdf)



