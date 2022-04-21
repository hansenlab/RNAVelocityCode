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


### get R2 for the top 100 genes
top100.m <- t(sapply(rownames(cc.df), function(g) {
	
	## ms t
	loess.l <- tricycle::fit_periodic_loess(theta.v = dynamical.o$tricyclePosition, y = assay(dynamical.o, "Ms")[g, ])
	cutoff1.v <- loess.l$pred.df$x[which.min(loess.l$pred.df$y)] 
	cutoff2.v <- loess.l$pred.df$x[which.max(loess.l$pred.df$y)] 
	
	if (cutoff1.v < cutoff2.v) {
		g_color.v <- ifelse((dynamical.o$tricyclePosition > cutoff1.v) & (dynamical.o$tricyclePosition < cutoff2.v), "up", "down")
	} else {
		g_color.v <- ifelse((dynamical.o$tricyclePosition < cutoff2.v) | (dynamical.o$tricyclePosition > cutoff1.v), "up", "down")
	}
	
	## inconsistent perct
	dynamical_incons.v <- mean((as.numeric(factor(g_color.v)) * 2 - 3) * sign(assay(dynamical.o, "velocity")[g, ]) == -1)
	steady_incons.v <- mean((as.numeric(factor(g_color.v)) * 2 - 3) * sign(assay(steady.o, "velocity")[g, ]) == -1)
	
	return(c(loess.l$rsquared, dynamical_incons.v, steady_incons.v))
}))

tmp.df <- data.frame(top100.m)
names(tmp.df) <- c("r2", "dynamical", "steady")
tmp.df %<>% filter(r2 > 0.5)

a.p <- plotxy(tmp.df$r2, tmp.df$dynamical * 100, title = "Dynamical model", x_lab = expression("Cell cycle position"~italic(R)^2),
			 y_lab = "Inconsistent percentage (%)", color_by = NULL, point.size = 4.01) +
	theme(legend.position = "none")

b.p <- plotxy(tmp.df$r2, tmp.df$steady * 100, title = "Steady model", x_lab = expression("Cell cycle position"~italic(R)^2),
							y_lab = "Inconsistent percentage (%)", color_by = NULL, point.size = 4.01) +
	theme(legend.position = "none")

c.p <- plotxy(tmp.df$dynamical * 100, tmp.df$steady * 100, title = "Comparison between models", x_lab = "Dynamical incons. perct. (%)",
							y_lab = "Steady incons. perct. (%)", color_by = NULL, point.size = 4.01) +
	geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
	theme(legend.position = "none")

d.p <- ggplot(tmp.df[, 2:3] %>% pivot_longer(
	cols = 1:2,
	names_to = "type",
	values_to = "perc"), aes(x = type, y = perc * 100)) +
	geom_boxplot(outlier.shape = 1, outlier.size = 0.2, size = 0.2, width = 0.5) +
	scale_x_discrete(name = "Model", labels = c("Dynamical", "Steady")) +
	labs( x = "Model", 
				y =  "Inconsistent percentage (%)", 
				title = "Comparison between models") 


mp <- plot_grid(a.p, b.p, c.p, d.p,
								nrow = 2, ncol = 2, label_size = 10, labels = "auto")
save_plot(here::here("figs", "sfigs", "sfig.fucci_incons_comp.pdf"), mp,
					base_height = 2 *1.2 / 1.4, base_width = 2*1.4 / 1.4, nrow = 2, ncol = 2, device = cairo_pdf)




