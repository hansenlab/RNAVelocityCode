rm(list=ls())
source(here::here("scripts/utils.R"))

steady.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))
dynamical.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))

### flip UMAP-1 
reducedDim(steady.o, "UMAP")[, 1] <-  - reducedDim(steady.o, "UMAP")[, 1]
reducedDim(dynamical.o, "UMAP")[, 1] <-  - reducedDim(dynamical.o, "UMAP")[, 1]


### use tricycle position instead of fucci time
steady.o <- tricycle::estimate_cycle_position(steady.o, gname.type = "SYMBOL", species = "human", exprs_values = "X")
dynamical.o <- tricycle::estimate_cycle_position(dynamical.o, gname.type = "SYMBOL", species = "human", exprs_values = "X")

library(org.Hs.eg.db)
go.v <- AnnotationDbi::select(org.Hs.eg.db, keytype="GOALL", keys="GO:0007049", 
															columns = "SYMBOL")[, "SYMBOL"]

cc.df <- rowData(dynamical.o)[rownames(dynamical.o) %in% go.v, ] %>% as.data.frame() %>%
	filter(`velocity_genes`) %>% arrange(desc(fit_likelihood))


### another Sfig for top 50 cc genes
lp <- do.call(c, lapply(lapply(1:5, function(i) seq(from = (i - 1) * 10 + 1, to = i * 10)), function(idx) {
	out <- lapply(rownames(cc.df)[idx], function(g) {
		## ms t
		loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$tricyclePosition, y = assay(dynamical.o, "Ms")[g, ])
		cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] 
		cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] 
		
		if (cutoff1.v < cutoff2.v) {
			g_color.v <- ifelse((dynamical.o$tricyclePosition > cutoff1.v) & (dynamical.o$tricyclePosition < cutoff2.v), "up", "down")
		} else {
			g_color.v <- ifelse((dynamical.o$tricyclePosition < cutoff2.v) | (dynamical.o$tricyclePosition > cutoff1.v), "up", "down")
		}
		
		ms_t.p <- plotFucciLoess3(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v,
															colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
															title = as.expression(substitute(italic(g)~" Ms over \u03b8", list(g=g)))) +
			theme(legend.position = "none") 
		
		message(g)
		return(ms_t.p)
	})
	return(out)
}))

mp <- plot_grid(plotlist = lp,
								nrow = 10, ncol = 5, label_size = 10, labels = NULL)
save_plot(here::here("figs", "sfigs", "sfig.fucci_top50.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 10, ncol = 5, device = cairo_pdf)

lp <- do.call(c, lapply(lapply(6:10, function(i) seq(from = (i - 1) * 10 + 1, to = i * 10)), function(idx) {
	out <- lapply(rownames(cc.df)[idx], function(g) {
		## ms t
		loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$tricyclePosition, y = assay(dynamical.o, "Ms")[g, ])
		cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] 
		cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] 
		
		if (cutoff1.v < cutoff2.v) {
			g_color.v <- ifelse((dynamical.o$tricyclePosition > cutoff1.v) & (dynamical.o$tricyclePosition < cutoff2.v), "up", "down")
		} else {
			g_color.v <- ifelse((dynamical.o$tricyclePosition < cutoff2.v) | (dynamical.o$tricyclePosition > cutoff1.v), "up", "down")
		}
		
		ms_t.p <- plotFucciLoess3(dynamical.o, assay(dynamical.o, "Ms")[g, ], y_lab = "Ms", color_var = g_color.v,
															colors.v = c("#80B1D3", "#FB8072"), labels.v = c("Down", "Up"), color.name = "Direction",
															title = as.expression(substitute(italic(g)~" Ms over \u03b8", list(g=g)))) +
			theme(legend.position = "none") 
		
		message(g)
		return(ms_t.p)
	})
	return(out)
}))

mp <- plot_grid(plotlist = lp,
								nrow = 10, ncol = 5, label_size = 10, labels = NULL)
save_plot(here::here("figs", "sfigs", "sfig.fucci_top50_100.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 10, ncol = 5, device = cairo_pdf)



