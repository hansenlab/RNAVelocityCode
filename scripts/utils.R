options(stringsAsFactors = F)

library("magrittr")
library("matrixStats")
library("Matrix")
library("RColorBrewer")
library("scales")
library("tidyverse")
library("cowplot")
library("ggbeeswarm")
library("viridis")
library("colorspace")
library("zeallot")
library("sparseMatrixStats")
library("qs")
library("ggridges")
library("scattermore")
library("SingleCellExperiment")
library("scater")
library("BiocParallel")
library(scuttle)
library(scran)
library(patchwork)
library(velociraptor)


### set ggplot theme
.new_theme <- theme_half_open(8) + 
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		#axis.line.x = element_blank(),
		#axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
		axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
		axis.text.y = element_text(size = 6),
		axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5),
		axis.title.y = element_text(size = 7),
		axis.title.x = element_text(size = 7),
		plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
		plot.margin = unit(c(3, 3, 4, 4), "pt"),
		legend.background = element_blank(),
		legend.position = "right",
		legend.justification = c(0.5, 0.5),
		legend.key.size = unit(7, "pt"),
		legend.key = element_blank(),
		legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
		legend.title = element_text(size = 5),
		legend.title.align = 0.5,
		legend.text.align = 0.5,
		legend.text = element_text(size = 5),
		strip.text = element_text(color = "white"),
		strip.background = element_rect(fill = "black", color = "black"))

.theme_frame <- theme_bw(base_size = 8) +
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		#axis.line.x = element_blank(),
		#axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
		axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
		axis.text.y = element_text(size = 6),
		axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5),
		axis.title.y = element_text(size = 7),
		axis.title.x = element_text(size = 7),
		plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
		plot.margin = unit(c(3, 3, 4, 4), "pt"),
		legend.background = element_blank(),
		legend.position = "right",
		legend.justification = c(0.5, 0.5),
		legend.key.size = unit(7, "pt"),
		legend.key = element_blank(),
		legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
		legend.title = element_text(size = 5),
		legend.title.align = 0.5,
		legend.text.align = 0.5,
		legend.text = element_text(size = 5),
		strip.text = element_text(color = "white"),
		strip.background = element_rect(fill = "black", color = "black"))


update_geom_defaults("point",list(size = 0.3))

.theme_noframe <- ggplot2::theme_minimal(base_size = 8) +
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		#axis.line.x = element_blank(),
		#axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
		axis.ticks = element_blank(),
		axis.text = element_blank(),
		axis.title.y = element_text(size = 7),
		axis.title.x = element_text(size = 7),
		plot.title = element_text(face = "plain", size = 8, hjust = 0.5),
		plot.margin = unit(c(3, 3, 4, 4), "pt"),
		legend.background = element_blank(),
		legend.position = "right",
		legend.justification = c(0.5, 0.5),
		legend.key.size = unit(7, "pt"),
		legend.key = element_blank(),
		legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
		legend.title = element_text(size = 5),
		legend.title.align = 0.5,
		legend.text.align = 0.5,
		legend.text = element_text(size = 5),
		strip.text = element_text(color = "white"),
		strip.background = element_rect(fill = "black", color = "black"))

theme_set(.new_theme)

theme_main <- theme_half_open(8) + 
	theme(
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		#axis.line.x = element_blank(),
		#axis.line.y = element_line(colour="black", linetype = 1, size= 0.25),
		axis.ticks = element_line(colour="black", linetype = 1, size= 0.25),
		axis.text.y = element_text(size = 6),
		axis.text.x = element_text(size = 6,  vjust = 0.5, hjust = 0.5),
		axis.title.y = element_text(size = 8),
		axis.title.x = element_text(size = 8),
		plot.title = element_text(face = "plain", size = 9, hjust = 0.5),
		plot.margin = unit(c(3, 3, 4, 4), "pt"),
		legend.background = element_blank(),
		legend.position = "right",
		legend.justification = c(0.5, 0.5),
		legend.key.size = unit(7, "pt"),
		legend.key = element_blank(),
		legend.margin = margin(t = 1, r = 2, b = 2, l = 2, unit = "pt"),
		legend.title = element_text(size = 5),
		legend.title.align = 0.5,
		legend.text.align = 0.5,
		legend.text = element_text(size = 5),
		strip.text = element_text(color = "white"),
		strip.background = element_rect(fill = "black", color = "black"))



### help function for annotate
.percent_range <- function(value.v, percent=0.9) {
	value.v <- na.omit(value.v)
	min(value.v) + percent*diff(range(value.v))
}




plot_us <- function(sce.o, rownm = NULL, gene = NULL, point.size = 3.01,
													line.size = 0.7,
													color_by = NULL, x_lab = "Spliced",
													y_lab = "Unspliced", color.name = NULL,
													title = NULL,
													hue.colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"),
													plot.legend = FALSE) {
	if (is.null(gene) & is.null(rownm)) stop("rownm or gene must be given.")
	if (is.null(gene)) gene <- rownm
	if (is.null(rownm)) rownm <- rownames(sce.o)[match(gene, rowData(sce.o)$Gene)]
	row.idx <- which(rownames(sce.o) == rownm)
	
	
	tmp.df <- data.frame(x = assay(sce.o, "spliced")[row.idx,],
											 y = assay(sce.o, "unspliced")[row.idx,],
											 color = colData(sce.o)[, color_by])
	if (is.null(color_by)) {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			
			labs( y = y_lab, x = x_lab, title = title)
	} else {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			scale_color_gradientn(name = color.name, limits = range(tmp.df$color), #labels = c(0, rep("", 498), 1),
														#breaks = seq(from = 0, to = 1, length.out = 500) ,
														colors = hue.colors) + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	return(scat.p)
}
plot_steadymums <- function(sce.o, rownm = NULL, gene = NULL, point.size = 3.01,
														line.size = 0.7,
														color_by = NULL, x_lab = "Ms",
														y_lab = "Mu", color.name = NULL,
														title = NULL,
														hue.colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"),
														hue.n = 500,
														plot.legend = FALSE) {
	if (is.null(gene) & is.null(rownm)) stop("rownm or gene must be given.")
	if (is.null(gene)) gene <- rownm
	if (is.null(rownm)) rownm <- rownames(sce.o)[match(gene, rowData(sce.o)$Gene)]
	row.idx <- which(rownames(sce.o) == rownm)
	gamma <- rowData(sce.o)$velocity_gamma[row.idx]
	
	tmp.df <- data.frame(x = assay(sce.o, "Ms")[row.idx,],
											 y = assay(sce.o, "Mu")[row.idx,],
											 color = colData(sce.o)[, color_by])
	if (is.null(color_by)) {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_abline(slope = gamma, intercept = 0, linetype = "dashed", size = line.size) +
			labs( y = y_lab, x = x_lab, title = title)
	} else {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_abline(slope = gamma, intercept = 0, linetype = "dashed", size = line.size) +
			scale_color_gradientn(name = color.name, limits = range(tmp.df$color), #labels = c(0, rep("", 498), 1),
														#breaks = seq(from = 0, to = 1, length.out = 500) ,
														colors = hue.colors) + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	return(scat.p)
}

plot_steadymums2 <- function(sce.o, rownm = NULL, gene = NULL, point.size = 3.01,
														line.size = 0.7,
														color_by = NULL, x_lab = "Ms",
														y_lab = "Mu", color.name = NULL,
														title = NULL,
														colors.v = NULL, labels.v = NULL,
														plot.legend = FALSE) {
	if (is.null(gene) & is.null(rownm)) stop("rownm or gene must be given.")
	if (is.null(gene)) gene <- rownm
	if (is.null(rownm)) rownm <- rownames(sce.o)[match(gene, rowData(sce.o)$Gene)]
	row.idx <- which(rownames(sce.o) == rownm)
	gamma <- rowData(sce.o)$velocity_gamma[row.idx]
	
	tmp.df <- data.frame(x = assay(sce.o, "Ms")[row.idx,],
											 y = assay(sce.o, "Mu")[row.idx,])
	if (is.null(color_by)) {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_abline(slope = gamma, intercept = 0, linetype = "dashed", size = line.size) +
			labs( y = y_lab, x = x_lab, title = title)
	} else {
		tmp.df$color <- colData(sce.o)[, color_by]
		if (is.null(labels.v)) labels.v <- levels(factor(tmp.df$color))
		scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_abline(slope = gamma, intercept = 0, linetype = "dashed", size = line.size) +
			scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = levels(factor(tmp.df$color))) +  
			labs( y = y_lab, x = x_lab, title = title)
	}
	return(scat.p)
}

plot_us2 <- function(sce.o, rownm = NULL, gene = NULL, point.size = 3.01,
										line.size = 0.7,
										color_by = NULL, x_lab = "Spliced",
										y_lab = "Unspliced", color.name = NULL,
										title = NULL,
										colors.v = NULL, labels.v = NULL,
										plot.legend = FALSE) {
	if (is.null(gene) & is.null(rownm)) stop("rownm or gene must be given.")
	if (is.null(gene)) gene <- rownm
	if (is.null(rownm)) rownm <- rownames(sce.o)[match(gene, rowData(sce.o)$Gene)]
	row.idx <- which(rownames(sce.o) == rownm)
	
	
	tmp.df <- data.frame(x = assay(sce.o, "spliced")[row.idx,],
											 y = assay(sce.o, "unspliced")[row.idx,],
											 color = colData(sce.o)[, color_by])
	if (is.null(color_by)) {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			
			labs( y = y_lab, x = x_lab, title = title)
	} else {
		tmp.df$color = colData(sce.o)[, color_by]
		if (is.null(labels.v)) labels.v <- levels(factor(tmp.df$color))
		scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = levels(factor(tmp.df$color)))  + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	return(scat.p)
}

dynamical_phase <- function(alpha, beta, gamma, scaling, t_, u0, s0, t, n.point = 1000) {
	### functions
	.updateEnvVar <- function(l, varnames = NULL, env = parent.frame()) {
		if (is.null(varnames)) varnames <- names(l)
		for (i in seq_along(l)) {
			assign(x = varnames[i], value = l[[i]], envir = env)
		}
	}
	
	.mRNA <- function(tau, u0, s0, alpha, beta, gamma) {
		expu <- exp(- beta * tau)
		exps <- exp(- gamma * tau)
		return(list(u = u0 * expu + alpha / beta * (1 - expu), s = s0 * exps + alpha / gamma * (1 - exps) + (alpha - u0 * beta)  * (exps - expu) * .invr(gamma - beta)))
	}
	
	.invr <- function(X) {
		o <- X != 0
		class(o) <- "numeric"
		X_invr <- 1 / X * o
		return(X_invr)
	}
	
	.vectorize <- function(t, t_, alpha, beta, gamma, alpha_ = 0, u0 = 0, s0 = 0) {
		o <- t < t_
		class(o) <- "numeric"
		tau = t * o + (t - t_) * (1 - o)
		u0_ = .unspliced(t_, u0, alpha, beta)
		s0_ = .spliced(t_, s0, u0, alpha, beta, gamma)
		# vectorize u0, s0 and alpha
		u0 = u0 * o + u0_ * (1 - o)
		s0 = s0 * o + s0_ * (1 - o)
		alpha = alpha * o + alpha_ * (1 - o)
		return(list(tau = tau, alpha = alpha, u0 = u0, s0 = s0))
	}
	
	.unspliced <- function(tau, u0, alpha, beta) {
		expu <- exp(-beta * tau)
		return(u0 * expu + alpha / beta * (1 - expu))
	}
	
	.spliced <- function(tau, s0, u0, alpha, beta, gamma) {
		co <- (alpha - u0 * beta) * .invr(gamma - beta)
		expu <- exp(-beta * tau)
		exps <- exp(-gamma * tau)
		return(s0 * exps + alpha / gamma * (1 - exps) + co * (exps - expu))
	}
	
	beta <- beta * scaling
	u0_offset <- u0
	s0_offset <- s0
	
	u0_ = .unspliced(t_, 0, alpha, beta) #steady u
	tmax <- max(t)
	# make new t for phase
	new_t <- c(seq(from = 0, to = t_, length.out = floor(n.point / 2)), seq(from = t_, to = tmax, length.out = floor(n.point / 2)))
	
	.updateEnvVar(l = .vectorize(sort(new_t), t_, alpha, beta, gamma), env = environment())
	.updateEnvVar(l = .mRNA(tau, u0, s0, alpha, beta, gamma), env = environment())
	
	ut <- u * scaling + u0_offset
	st <- s + s0_offset
	return(data.frame(ut = ut, st = st))
}

plot_dynamicalmums <- function(sce.o, rownm = NULL, gene = NULL, point.size = 3.01,
														line.size = 0.7,
														color_by = NULL, x_lab = "Ms",
														y_lab = "Mu", color.name = NULL,
														title = NULL,
														hue.colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"),
														hue.n = 500,
														plot.legend = FALSE) {
	if (is.null(gene) & is.null(rownm)) stop("rownm or gene must be given.")
	if (is.null(gene)) gene <- rownm
	if (is.null(rownm)) rownm <- rownames(sce.o)[match(gene, rowData(sce.o)$Gene)]
	row.idx <- which(rownames(sce.o) == rownm)
	gamma <- rowData(sce.o)$fit_gamma[row.idx]
	
	rowData.df <- rowData(sce.o)
	tmp.df <- data.frame(x = assay(sce.o, "Ms")[row.idx,],
											 y = assay(sce.o, "Mu")[row.idx,],
											 color = colData(sce.o)[, color_by],
											 t = assay(sce.o, "fit_t")[row.idx, ])
	
	
	phase.df <- dynamical_phase(alpha = rowData.df$fit_alpha[row.idx], beta = rowData.df$fit_beta[row.idx],
															gamma = rowData.df$fit_gamma[row.idx], scaling = rowData.df$fit_scaling[row.idx], t_ = rowData.df$fit_t_[row.idx],
															u0 = rowData.df$fit_u0[row.idx], s0 = rowData.df$fit_s0[row.idx], t = tmp.df$t)
	
	ymax.v <- rowData.df$fit_gamma[row.idx] / rowData.df$fit_beta[row.idx] * (max(phase.df$st) - min(phase.df$st)) + min(phase.df$ut)
	
	
	if (is.null(color_by)) {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_path(data = phase.df, aes(x = st, y = ut), color = "#800080") +
			geom_segment(aes(x = min(phase.df$st), xend = max(phase.df$st), y = min(phase.df$ut), yend = ymax.v), linetype = "dashed", color = '#800080') +
			labs( y = y_lab, x = x_lab, title = title)
	} else {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_path(data = phase.df, aes(x = st, y = ut), color = "#800080") +
			geom_segment(aes(x = min(phase.df$st), xend = max(phase.df$st), y = min(phase.df$ut), yend = ymax.v), linetype = "dashed", color = '#800080') +
			scale_color_gradientn(name = color.name, limits = range(tmp.df$color), #labels = c(0, rep("", 498), 1),
														#breaks = seq(from = 0, to = 1, length.out = 500) ,
														colors = hue.colors) + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	return(scat.p)
}

plot_dynamicalmums2 <- function(sce.o, rownm = NULL, gene = NULL, point.size = 3.01,
															 line.size = 0.7,
															 color_by = NULL, x_lab = "Ms",
															 y_lab = "Mu", color.name = NULL,
															 title = NULL,
																colors.v = NULL, labels.v = NULL,
															 plot.legend = FALSE) {
	if (is.null(gene) & is.null(rownm)) stop("rownm or gene must be given.")
	if (is.null(gene)) gene <- rownm
	if (is.null(rownm)) rownm <- rownames(sce.o)[match(gene, rowData(sce.o)$Gene)]
	row.idx <- which(rownames(sce.o) == rownm)[1]
	gamma <- rowData(sce.o)$fit_gamma[row.idx]
	
	rowData.df <- rowData(sce.o)
	tmp.df <- data.frame(x = assay(sce.o, "Ms")[row.idx,],
											 y = assay(sce.o, "Mu")[row.idx,],
											 t = assay(sce.o, "fit_t")[row.idx, ])
	
	phase.df <- dynamical_phase(alpha = rowData.df$fit_alpha[row.idx], beta = rowData.df$fit_beta[row.idx],
															gamma = rowData.df$fit_gamma[row.idx], scaling = rowData.df$fit_scaling[row.idx], t_ = rowData.df$fit_t_[row.idx],
															u0 = rowData.df$fit_u0[row.idx], s0 = rowData.df$fit_s0[row.idx], t = tmp.df$t)
	
	ymax.v <- rowData.df$fit_gamma[row.idx] / rowData.df$fit_beta[row.idx] * (max(phase.df$st) - min(phase.df$st)) + min(phase.df$ut)
	
	
	if (is.null(color_by)) {
		scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_path(data = phase.df, aes(x = st, y = ut), color = "#800080") +
			geom_segment(aes(x = min(phase.df$st), xend = max(phase.df$st), y = min(phase.df$ut), yend = ymax.v), linetype = "dashed", color = '#800080') +
			labs( y = y_lab, x = x_lab, title = title)
	} else {
		tmp.df$color <- colData(sce.o)[, color_by]
		if (is.null(labels.v)) labels.v <- levels(factor(tmp.df$color))
		scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
			geom_scattermore(pointsize = point.size,  alpha = 0.8) +
			geom_path(data = phase.df, aes(x = st, y = ut), color = "#800080") +
			geom_segment(aes(x = min(phase.df$st), xend = max(phase.df$st), y = min(phase.df$ut), yend = ymax.v), linetype = "dashed", color = '#800080') +
			scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = levels(factor(tmp.df$color))) + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	return(scat.p)
}


plotxy <- function(x, y, title = NULL, point.size = 3.01, colors.v = NULL,
									 hue.colors = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"),
													 x_lab = "",  y_lab = "", color_by = NULL, color.name = NULL) {
	
	tmp.df <- data.frame(x = x, y= y)
	if (is.null(color_by)) {
		p <- ggplot(data = tmp.df , aes(x = x, y = y)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			labs( y = y_lab , x = x_lab, title = title) 
	} else {
		if (is.null(colors.v)) colors.v <- colData(sce.o)[, color_by]
		tmp.df$color <- colors.v
		p <- ggplot(data = tmp.df , aes(x = x, y = y, color = color)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			scale_color_gradientn(name = color.name, limits = range(tmp.df$color), #labels = c(0, rep("", 498), 1),
														#breaks = seq(from = 0, to = 1, length.out = 500) ,
														colors = hue.colors) + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	

	return(p)
}
plotxy2 <- function(x, y, title = NULL, point.size = 3.01, color_var = NULL,
									 colors.v = NULL, labels.v = NULL,
									 x_lab = "",  y_lab = "",  color.name = NULL) {
	
	tmp.df <- data.frame(x = x, y= y)
	if (is.null(color_var)) {
		p <- ggplot(data = tmp.df , aes(x = x, y = y)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			labs( y = y_lab , x = x_lab, title = title) 
	} else {

		tmp.df$color <- color_var
		p <- ggplot(data = tmp.df , aes(x = x, y = y, color = color)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = levels(factor(tmp.df$color))) + 
			labs( y = y_lab, x = x_lab, title = title)
	}
	
	
	return(p)
}


plotFucciLoess <- function(sce.o, y, title = NULL, point.size = 3.01, 
													 color.name = "FUCCI",
													 #col.outname = NULL, 
													 x_lab = "FUCCI pseudotime",  colors.v = NULL,
													 hue.n = 500, log2.trans = FALSE, y_lab = NULL, addR2 = TRUE, r2size = 2.5,
													 x.pos = 0, y.pos = 1) {

	if (is.null(colors.v)) colors.v <- c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")
	if (log2.trans) y <- log2(y + 1)
	#if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))
	nna.idx <- which(!is.na(sce.o$fucci_time))
	sce.o <- sce.o[, nna.idx]
	y <- y[nna.idx]
	theta.v <- sce.o$fucci_time * 2 * pi
	
	tmp.df <- data.frame(x = sce.o$fucci_time, theta = theta.v, color = theta.v, y = y)
	scale_color <- scale_color_gradientn(name = color.name, limits = range(0, 1), 
																			 colors = colors.v)
	
	loess.l <- tricycle::fit_periodic_loess(theta.v = theta.v, y = y)
	
	p <- ggplot(data = tmp.df , aes(x = x , y = y, color = x)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = 0.8) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = 0.8) +
		geom_path(data = loess.l$pred.df, aes(x = x / (2 * pi) , y = y), linetype = "dashed", color = "black", size = 0.7, alpha = 0.6, inherit.aes = FALSE) +
		scale_color +
		labs( y = y_lab , x = x_lab, title = title) +
		scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) , labels = c(0, str_c(seq(0.25, 1, 0.25))), limits = c(0, 1)) 
	if (addR2) p <- p + annotate(geom = "text", x = .percent_range(tmp.df$x, x.pos), y = .percent_range(tmp.df$y, y.pos), size = r2size, hjust = 0, vjust = 1,
															 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)
	
	return(p)
}


plotFucciLoess2 <- function(sce.o, y, title = NULL, point.size = 3.01, 
													 color.name = NULL, color_var = NULL,
														labels.v = NULL,
													 #col.outname = NULL, 
													 x_lab = "FUCCI pseudotime",  colors.v = NULL,
													 hue.n = 500, log2.trans = FALSE, y_lab = NULL, addR2 = TRUE, r2size = 2.5,
													 x.pos = 0, y.pos = 1) {
	
	if (log2.trans) y <- log2(y + 1)
	#if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))
	nna.idx <- which(!is.na(sce.o$fucci_time))
	sce.o <- sce.o[, nna.idx]
	y <- y[ nna.idx]
	theta.v <- sce.o$fucci_time * 2 * pi
	
	tmp.df <- data.frame(x = sce.o$fucci_time, theta = theta.v, y = y)
	# scale_color <- scale_color_gradientn(name = color.name, limits = range(0, 1), 
	# 																		 colors = colors.v)
	loess.l <- tricycle::fit_periodic_loess(theta.v = theta.v, y = y)
	
	if (is.null(color_var)) {
		p <- ggplot(data = tmp.df , aes(x = x , y = y)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			geom_path(data = loess.l$pred.df, aes(x = x / (2 * pi) , y = y), linetype = "dashed", color = "black", size = 0.7, alpha = 0.6, inherit.aes = FALSE) +
			labs( y = y_lab , x = x_lab, title = title) +
			scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) , labels = c(0, str_c(seq(0.25, 1, 0.25))), limits = c(0, 1)) 
	} else {
		tmp.df$color <- color_var
		p <- ggplot(data = tmp.df , aes(x = x , y = y, color = color)) +
			geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size, alpha = 0.8) +
			geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size, alpha = 0.8) +
			geom_path(data = loess.l$pred.df, aes(x = x / (2 * pi) , y = y), linetype = "dashed", color = "black", size = 0.7, alpha = 0.6, inherit.aes = FALSE) +
			scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = levels(factor(tmp.df$color))) + 
			labs( y = y_lab , x = x_lab, title = title) +
			scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1) , labels = c(0, str_c(seq(0.25, 1, 0.25))), limits = c(0, 1)) 
	}
	

	if (addR2) p <- p + annotate(geom = "text", x = .percent_range(tmp.df$x, x.pos), y = .percent_range(tmp.df$y, y.pos), size = r2size, hjust = 0, vjust = 1,
															 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)
	
	return(p)
}

plotEmbExp <- function(sce.o, dimred, 
															 color.v,
															 color.name,
											 point.alpha = 0.8,
											 point.size = 3.01,
															 title = NULL,
															 x_lab = NULL,
															 y_lab = NULL) {
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = color.v)
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = point.size,  alpha = point.alpha) +
		scale_color_gradientn(colours = c("#440154", "#482576", "#414487",
																			"#35608D", "#2A788E", "#21908C",
																			"#22A884", "#43BF71", "#7AD151",
																			"#BBDF27", "#FDE725"),name = color.name, limits = range(tmp.df$color, na.rm = TRUE)) + 
		labs( y = y_lab, x = x_lab) +
		ggtitle(title) +
		xlim(x_lim) + ylim(y_lim)
	return(scat.p)
}

plotEmbExp_nopoint <- function(sce.o, dimred, 
											 color.v,
											 color.name,
											 xintercept = -3.5,
											 yintercept = 0,
											 point.alpha = 0.8,
											 point.size = 3.01,
											 title = NULL,
											 x_lab = NULL,
											 y_lab = NULL) {
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = color.v)
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y)) +
		geom_vline(xintercept = xintercept, color = "red", linetype = "dashed", size = 0.6) +
		geom_hline(yintercept = yintercept, color = "red", linetype = "dashed", size = 0.6)  + 
		labs( y = y_lab, x = x_lab) +
		ggtitle(title) +
		xlim(x_lim) + ylim(y_lim)
	return(scat.p)
}

plotEmbVelo <- function(sce.o, dimred, 
											 color.v,
											 color.name,
											 point.size = 3.01,
												point.alpha = 1,
											 title = NULL,
											 x_lab = NULL,
											 y_lab = NULL) {
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = color.v)
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = point.size,  alpha = point.alpha) +
		scale_color_gradientn(colours = c("#A50026", "#D73027", "#F46D43",
																			"#FDAE61", "#FEE08B", "#FFFFBF",
																			"#D9EF8B", "#A6D96A", "#66BD63",
																			"#1A9850", "#006837"), limits = max(abs(tmp.df$color))*c(-1, 1), name = color.name) +
		labs( y = y_lab, x = x_lab) +
		ggtitle(title) +
		xlim(x_lim) + ylim(y_lim)
	return(scat.p)
}

plotEmbFucci <- function(sce.o, dimred, 
											 color.v = c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9"),
											 color.name = "FUCCI",
											 point.size = 3.01,
											 title = NULL,
											 x_lab = NULL,
											 y_lab = NULL) {
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = sce.o$fucci_time)
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(pointsize = point.size,  alpha = 0.8) +
		scale_color_gradientn(name = color.name, limits = range(0, 1), #labels = c(0, rep("", 498), 1),
													#breaks = seq(from = 0, to = 1, length.out = 500) ,
													colors = color.v) + 
		labs( y = y_lab, x = x_lab) +
		ggtitle(title) +
		xlim(x_lim) + ylim(y_lim)
	return(scat.p)
}



plotEmbScat <- function(sce.o, dimred, 
												color_by,
												colors.v,
												labels.v,
												color.name,
												x_lab = NULL,
												y_lab = NULL,
												facet_by = NULL,
												facet_labels = NULL,
												point.size = 3.01,
												point.alpha = 0.8,
												title = NULL) {
	
	emb.m <- reducedDim(sce.o, dimred)
	
	x_lim <- range(emb.m[, 1]) + c(diff(range(emb.m[, 1])) * c(-0.05, 0.05))
	y_lim <- range(emb.m[, 2]) + c(diff(range(emb.m[, 2])) * c(-0.05, 0.05))
	
	tmp.df <- data.frame(x = emb.m[, 1],
											 y = emb.m[, 2],
											 color = colData(sce.o)[, color_by])
	if (any(is.na(tmp.df$color))) {
		tmp.df$color <- fct_explicit_na(tmp.df$color, na_level = "NA") %>% fct_relevel("NA", after = Inf)
		scale_color <- scale_color_manual(values = c(colors.v, "grey"), name = color.name, labels =  c(labels.v, "NA"), limits =   c(labels.v, "NA"))
	} else {
		scale_color <- scale_color_manual(values = colors.v, name = color.name, labels = labels.v, limits = labels.v)
	}
	
	
	scat.p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` == "NA"), pointsize = point.size,  alpha = point.alpha) +
		geom_scattermore(data = tmp.df %>% dplyr::filter(`color` != "NA"), pointsize = point.size,  alpha = point.alpha) +
		scale_color +
		guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
		labs( y = y_lab, x = x_lab, title = title) +
		xlim(x_lim) + ylim(y_lim)
	
	if (!is.null(facet_by)) {
		if (is.null(facet_labels)) facet_labels <- levels(factor(colData(sce.o)[, facet_by]))
		tmp.df$facet <- factor(colData(sce.o)[, facet_by])
		lp <- lapply(seq_len(nlevels(factor(tmp.df$facet))), function(idx) {
			p <- ggplot(tmp.df, aes(x = x, y = y, color = color)) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` != levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = 0.4, color = "gray90", show.legend = FALSE) +
				geom_scattermore(data = tmp.df %>% dplyr::filter(`facet` == levels(factor(tmp.df$facet))[idx]), pointsize = point.size,  alpha = point.alpha) +
				scale_color +
				guides(color = guide_legend(override.aes = list(alpha = 1, size = 1))) + 
				labs( y = y_lab, x = x_lab, title = str_c(facet_labels[idx], " (n=", sum(tmp.df$facet == levels(factor(tmp.df$facet))[idx], na.rm = TRUE), ")")) +
				xlim(x_lim) + ylim(y_lim)
			return(p)
		})
		return(c(list(scat.p), lp))
	}
	return(scat.p)
}



plot_comV <- function(steady.o, dynamical.o, gene, point.size = 3.01, title = NULL, r2size = 2, x.pos = 0.6, y.pos = 0.3) {
	if (is.null(title)) title <- gene
	tmp.df <- data.frame(x = assay(steady.o, "velocity")[which(rownames(steady.o) == gene),],
											 y = assay(dynamical.o, "velocity")[which(rownames(dynamical.o) == gene),])
	tmp.df$color <- ifelse(sign(tmp.df$x * tmp.df$y) == -1, "non", "up")
	tmp.df$color[(tmp.df$x < 0) & (tmp.df$y < 0)] <- "down"
	pcc.v <- cor(tmp.df$x, tmp.df$y)
	inconsistency.v <- mean(sign(tmp.df$x * tmp.df$y) == -1)
	
	p <- ggplot(data = tmp.df, aes(x = x, y = y, color = color)) + 
		geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, size = 0.5) +
		geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5, size = 0.5) +
		geom_scattermore(pointsize = point.size,  alpha = 0.8) +
		scale_color_manual(values = c("#80B1D3", "black","#FB8072"), name = "Direction", labels = c("Down", "Incons.", "Up")) +
		labs(x = "Velocity - steady model", y = "Velocity - dynamical model", title = title) +
		annotate(geom = "text", x = .percent_range(tmp.df$x, x.pos), y = .percent_range(tmp.df$y, y.pos - 0.1), size = r2size, hjust = 0, vjust = 1,
						 label = str_c("PCC: ", sprintf("%.3f", pcc.v)), parse = FALSE) +
		annotate(geom = "text", x = .percent_range(tmp.df$x, x.pos), y = .percent_range(tmp.df$y, y.pos - 0.2), size = r2size, hjust = 0, vjust = 1,
						 label = str_c("Incons.: ", sprintf("%0.1f", inconsistency.v * 100), "%"), parse = FALSE)
	
	
	return(p)
	
}

dynamical_thet <- function(alpha, beta, gamma, scaling, t_, u0, s0, tmax, new_t) {
	.updateEnvVar <- function(l, varnames = NULL, env = parent.frame()) {
		if (is.null(varnames)) varnames <- names(l)
		for (i in seq_along(l)) {
			assign(x = varnames[i], value = l[[i]], envir = env)
		}
	}
	
	.mRNA <- function(tau, u0, s0, alpha, beta, gamma) {
		expu <- exp(- beta * tau)
		exps <- exp(- gamma * tau)
		return(list(u = u0 * expu + alpha / beta * (1 - expu), s = s0 * exps + alpha / gamma * (1 - exps) + (alpha - u0 * beta)  * (exps - expu) * .invr(gamma - beta)))
	}
	
	.invr <- function(X) {
		o <- X != 0
		class(o) <- "numeric"
		X_invr <- 1 / X * o
		return(X_invr)
	}
	
	.vectorize <- function(t, t_, alpha, beta, gamma, alpha_ = 0, u0 = 0, s0 = 0) {
		o <- t < t_
		class(o) <- "numeric"
		tau = t * o + (t - t_) * (1 - o)
		u0_ = .unspliced(t_, u0, alpha, beta)
		s0_ = .spliced(t_, s0, u0, alpha, beta, gamma)
		# vectorize u0, s0 and alpha
		u0 = u0 * o + u0_ * (1 - o)
		s0 = s0 * o + s0_ * (1 - o)
		alpha = alpha * o + alpha_ * (1 - o)
		return(list(tau = tau, alpha = alpha, u0 = u0, s0 = s0))
	}
	
	.unspliced <- function(tau, u0, alpha, beta) {
		expu <- exp(-beta * tau)
		return(u0 * expu + alpha / beta * (1 - expu))
	}
	
	.spliced <- function(tau, s0, u0, alpha, beta, gamma) {
		co <- (alpha - u0 * beta) * .invr(gamma - beta)
		expu <- exp(-beta * tau)
		exps <- exp(-gamma * tau)
		return(s0 * exps + alpha / gamma * (1 - exps) + co * (exps - expu))
	}
	
	beta <- beta * scaling
	u0_offset <- u0
	s0_offset <- s0
	
	u0_ = .unspliced(t_, 0, alpha, beta) #steady u
	
	.updateEnvVar(l = .vectorize(new_t , t_, alpha, beta, gamma), env = environment())
	.updateEnvVar(l = .mRNA(tau, u0, s0, alpha, beta, gamma), env = environment())
	
	ut <- u * scaling + u0_offset
	st <- s + s0_offset
	vt <- (ut * beta)/scaling - st * gamma  # ds/dt
	wt <- (alpha - beta * ut) * scaling  # du/dt
	
	# su_t.df <- data.frame(ut = ut, st = st, t = new_t) %>% gather(key = "type", value = "value", 1:2)
	# vw_t.df <- data.frame(vt = vt, wt = wt, t = new_t) %>% gather(key = "type", value = "value", 1:2)
	return(data.frame(ut = ut, st = st, new_t = new_t, vt = vt, wt = wt))
}


plotVelocityStream <- function(sce, embedded, use.dimred = 1, point.size = 3.01,
															 color_by = "#444444", color.alpha = 0.2,
															 grid.resolution = 60, scale = TRUE,
															 stream.L = 10, stream.min.L = 0, stream.res = 4,
															 stream.width = 8,
															 arrow.angle = 8, arrow.length = 0.8) {
	if (!identical(ncol(sce), nrow(embedded))) {
		stop("'sce' and 'embedded' do not have consistent dimensions.")
	}
	if (is.numeric(use.dimred)) {
		stopifnot(exprs = {
			identical(length(use.dimred), 1L)
			use.dimred <= length(reducedDims(sce))
		})
		use.dimred <- reducedDimNames(sce)[use.dimred]
	}
	else if (is.character(use.dimred)) {
		stopifnot(exprs = {
			length(use.dimred) == 1L
			use.dimred %in% reducedDimNames(sce)
		})
	}
	else {
		stop("'use.dimred' is not a valid value for use in reducedDim(sce, use.dimred)")
	}
	if (!requireNamespace("ggplot2")) {
		stop("'plotVelocityStream' requires the package 'ggplot2'.")
	}
	
	# get coordinates in reduced dimensional space
	xy <- reducedDim(sce, use.dimred)[, 1:2]
	colnames(xy) <-c("V1", "V2")
	
	# plot it using ggplot2 and metR::geom_streamline
	plotdat1 <- data.frame(xy)
	colnames(plotdat1) <- c("x", "y")
	if (is.character(color_by) && length(color_by) == 1L && color_by %in% colnames(colData(sce))) {
		plotdat1 <- cbind(plotdat1, col = colData(sce)[, color_by])
		colByFeat <- TRUE
	} else {
		colByFeat <- FALSE
	}
	
	idx <- which(!is.na(embedded[, 1]))
	sce <- sce[, idx]
	xy <- xy[idx, ]
	embedded <- embedded[idx, ]
	colnames(embedded) <- c("V1", "V2")
	# summarize velocities in a grid
	gr <- gridVectors(x = xy, embedded = embedded,
										resolution = grid.resolution, scale = scale,
										as.data.frame = FALSE,
										return.intermediates = TRUE)
	
	# now make it a regular grid needed for metR::geom_streamline
	xbreaks <- seq(gr$limits[1,1], gr$limits[2,1], by = gr$delta[1])
	ybreaks <- seq(gr$limits[1,2], gr$limits[2,2], by = gr$delta[2])
	plotdat2 <- expand.grid(x = xbreaks + gr$delta[1] / 2,
													y = ybreaks + gr$delta[2] / 2,
													dx = 0, dy = 0)
	allcategories <- DataFrame(expand.grid(V1 = seq(0, grid.resolution),
																				 V2 = seq(0, grid.resolution)))
	ivec <- match(gr$categories[sort(unique(gr$grp)), ], allcategories)
	plotdat2[ivec, c("dx", "dy")] <- gr$vec

	

	p <- ggplot2::ggplot(plotdat1, ggplot2::aes(x = !!ggplot2::sym("x"), y = !!ggplot2::sym("y"))) +
		ggplot2::labs(x = paste(use.dimred, "1"), y = paste(use.dimred, "2"))
	if (!colByFeat) {
		colMatrix <- grDevices::col2rgb(col = color_by, alpha = TRUE)
		if (any(colMatrix[4, ] != 255)) {
			warning("ignoring 'color.alpha' as 'color_by' already specifies alpha channels")
			color.alpha <- colMatrix[4, ] / 255
		}
		p <- p + geom_scattermore(color = color_by, pointsize = point.size,  alpha = color.alpha) 
	} else {
		p <- p + geom_scattermore(ggplot2::aes(color = !!ggplot2::sym("col")), alpha = color.alpha, pointsize = point.size) 
	}
	p <- p +
		metR::geom_streamline(mapping = ggplot2::aes(x = !!ggplot2::sym("x"),
																								 y = !!ggplot2::sym("y"),
																								 dx = !!ggplot2::sym("dx"),
																								 dy = !!ggplot2::sym("dy"),
																								 size = stream.width * !!ggplot2::sym("..step..")),
													data = plotdat2, size = 0.3, jitter = 2,
													L = stream.L, min.L = stream.min.L, 
													res = stream.res, arrow.angle = arrow.angle,
													arrow.length = arrow.length, inherit.aes = FALSE) +
		ggplot2::theme_minimal() +
		ggplot2::theme(axis.text = ggplot2::element_blank(),
									 panel.grid.major = ggplot2::element_blank(),
									 panel.grid.minor = ggplot2::element_blank())
	
	return(p)
}



fit_t_loess <- function(x, y, span = 0.3, length.out = 200, predict = TRUE, ...) {
	loess.o <- loess(y ~ x, span = span, ...)
	fitted.v <- loess.o$fitted
	residual.v <- loess.o$residuals
	ss.total <- sum(scale(y, scale = FALSE) ^ 2)
	rsquared <- 1 - sum(residual.v ^ 2) / ss.total
	
	if (!predict) {
		return(list(fitted = fitted.v, residual = residual.v, loess.o = loess.o, rsquared = rsquared))
	}
	pred.x <- seq(from = min(x), to = max(x), length.out = length.out)
	pred.y <- predict(loess.o, newdata = data.frame(x = pred.x))
	pred.df <-  data.frame(x = pred.x, y = pred.y)
	return(list(fitted = fitted.v, residual = residual.v, pred.df = pred.df, loess.o = loess.o, rsquared = rsquared))
}

plotTLoess <- function(x, y, title = NULL, point.size = 3.01, 
											 color_var = NULL, scale_color = NULL,
											 color.name = "FUCCI",
											 #col.outname = NULL, 
											 x_lab = "Velocity latent time",  colors.v = NULL,
											 hue.n = 500, log2.trans = FALSE, y_lab = NULL, addR2 = TRUE, r2size = 2.5,
											 x.pos = 0, y.pos = 1) {
	
	if (is.null(colors.v)) colors.v <- c("#2E22EA","#9E3DFB","#F86BE2","#FCCE7B","#C4E416","#4BBA0F","#447D87","#2C24E9")
	if (log2.trans) y <- log2(y + 1)
	#if (is.null(y_lab)) y_lab <- bquote(paste('log'['2'],'(expression of ', .(col.outname), ")"))
	na.idx <- which(!is.na(x))
	x <- x[na.idx]
	y <- y[na.idx]
	
	tmp.df <- data.frame(x = x, y = y)
	loess.l <- fit_t_loess(x = x, y = y)
	
	if (is.null(color_var)) {
		p <- ggplot(data = tmp.df , aes(x = x , y = y)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			geom_path(data = loess.l$pred.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.7, alpha = 0.6, inherit.aes = FALSE) +
			labs( y = y_lab , x = x_lab, title = title) 
	} else {
		tmp.df$color <- color_var
		if (is.null(scale_color)) stop("Need scale_color")
		p <- ggplot(data = tmp.df , aes(x = x , y = y, color = color)) +
			geom_scattermore(data = tmp.df, pointsize = point.size, alpha = 0.8) +
			geom_path(data = loess.l$pred.df, aes(x = x , y = y), linetype = "dashed", color = "black", size = 0.7, alpha = 0.6, inherit.aes = FALSE) +
			scale_color +
			labs( y = y_lab , x = x_lab, title = title) 
	}
	
	if (addR2) p <- p + annotate(geom = "text", x = .percent_range(tmp.df$x, x.pos), y = .percent_range(tmp.df$y, y.pos), size = r2size, hjust = 0, vjust = 1,
															 label = as.character(as.expression(substitute(italic(R)^2~"="~rsquared, list(rsquared = format(loess.l$rsquared, digits = 3))))), parse = TRUE)
	
	return(p)
}






runSeuratCC <- function(sce.o, gname, assay.type = 'spliced', species = c("mouse", "human"), id.type = c("ensembl", "name")) {
	require(Seurat)
	s.genes <- cc.genes$s.genes
	g2m.genes <- cc.genes$g2m.genes
	species <- match.arg(species)
	if (species == "mouse") {
		AnnotationDb <- org.Mm.eg.db::org.Mm.eg.db
	} else {
		AnnotationDb <- org.Hs.eg.db::org.Hs.eg.db
	}
	
	if (id.type == "ensembl") {
		gname <- suppressMessages(AnnotationDbi::mapIds(AnnotationDb, keys = x$first, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
	}
	counts.m <- assay(sce.o, assay.type)
	rownames(counts.m) <- toupper(gname)
	
	### remove NA gene names
	idx <- which(!is.na(rownames(counts.m)))
	o <- CreateSeuratObject(counts = counts.m[idx, ])
	o <- NormalizeData(o)
	o <- FindVariableFeatures(o, selection.method = "vst")
	o<- ScaleData(o, features = rownames(o))
	o <- CellCycleScoring(o, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
	return(factor(o$Phase, levels = c("G1", "S", "G2M")))
}

seuratIntegrate <- function(count.m, batch, nfeatures = 2000) {
	require(Seurat)
	seurat.o <- CreateSeuratObject(counts = count.m)
	seurat.o[["batch"]] <- batch
	
	seurat.list <- SplitObject(seurat.o, split.by = "batch")
	
	for (i in 1:length(seurat.list)) {
		seurat.list[[i]] <- NormalizeData(seurat.list[[i]], verbose = TRUE)
		seurat.list[[i]] <- FindVariableFeatures(seurat.list[[i]], selection.method = "vst", 
																						 nfeatures = nfeatures, verbose = TRUE)
	}
	
	seurat.anchors <- FindIntegrationAnchors(object.list = seurat.list, dims = 1:30)
	seurat.integrated <- IntegrateData(anchorset = seurat.anchors, dims = 1:30)
	
	corrected.m <- seurat.integrated@assays$integrated@data
	
	return(corrected.m)
}



nrmse <- function(true.v, estimation.v) {
	Metrics::rmse(true.v, estimation.v) / sd(true.v, na.rm = TRUE)
}
