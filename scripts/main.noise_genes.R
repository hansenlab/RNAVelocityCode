rm(list = ls())
source(here::here("scripts/utils.R"))

noise_level = 3
n_obs = 500L
n_vars = 10L

scv <- reticulate::import("scvelo")
adata <-
	scv$datasets$simulation(
		n_obs = n_obs,
		n_vars = n_vars,
		t_max = 25,
		alpha = 5,
		beta = .3,
		gamma = .5,
		noise_level = noise_level
	)
ncomp <- as.integer(min(30L, n_vars - 1L))


scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)
rotation.m <- adata$varm["PCs"]



### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)


###
i <- 1
tmp.m <-
	dynamical_thet(
		alpha = rowData(dynamical.o)$true_alpha[i],
		beta = rowData(dynamical.o)$true_beta[i],
		gamma = rowData(dynamical.o)$true_gamma[i],
		scaling = rowData(dynamical.o)$true_scaling[i],
		t_ = rowData(dynamical.o)$true_t_[i],
		u0 = 0,
		s0 = 0,
		tmax = 25,
		new_t = dynamical.o$true_t
	)
tmp.df <-
	data.frame(
		x = dynamical.o$true_t,
		s = tmp.m[, 2],
		v = tmp.m[, 4],
		u = tmp.m[, 1],
		est = assay(dynamical.o, "fit_t")[i,],
		ms = assay(dynamical.o, "Ms")[i,],
		mu = assay(dynamical.o, "Mu")[i,],
		spliced = assay(dynamical.o, "spliced")[i,],
		unspliced = assay(dynamical.o, "unspliced")[i,],
		steady_v = assay(steady.o, "velocity")[i, ],
		dynamical_v = assay(dynamical.o, "velocity")[i, ]
	)

mums.p <-
	plot_dynamicalmums(
		dynamical.o,
		rownm = as.character(i - 1),
		color_by = "true_t",
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		title = str_c("Using learned k-NN"),
		color.name = "True t",
		point.size = 5.01
	) +
	geom_abline(
		slope = rowData(steady.o[i,])$velocity_gamma,
		intercept = 0,
		linetype = "dashed"
	) +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0)) +
	xlab(TeX(r'(Ms$_g$)')) + ylab(TeX(r'(Mu$_g$)'))


ms_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$ms,
		title = "Using learned k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Ms$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = s),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")

mu_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$mu,
		title = "Using learned k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Mu$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = u),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")

spliced_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$spliced,
		title = "Using learned k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Spliced$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = s),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(
		legend.position = c(1, 1),
		legend.justification = c(1, 1),
		legend.key.size = unit(6, "pt")
	)

unspliced_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$unspliced,
		title = "Using learned k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Unspliced$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = u),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")



dynamical_t.p <-
	plotxy_zeoroy(
		tmp.df$x,
		tmp.df$dynamical_v,
		title = "Using learned k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Dynamical v$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = v),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")


# steady_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True v.",
# 										 y_lab = "Steady v.", color_by = "true_t", point.size = 4.01,
# 										 colors.v = tmp.df$x,
# 										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
# 										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
# 										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
# 	.theme_noframe +
# 	theme(plot.title = element_blank(), legend.position = "none",
# 				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_v.p <-
	plotxy_zeoroy(
		tmp.df$v,
		tmp.df$dynamical_v,
		title = "Using learned k-NN",
		x_lab = TeX(r'(True v$_g$)'),
		y_lab = TeX(r'(Dynamical v$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none") +
	annotate(
		geom = "text",
		x = .percent_range(tmp.df$v, 0.6),
		y = .percent_range(tmp.df$dynamical_v, 0.15),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("PCC: ", sprintf(
			"%.3f", cor(tmp.df$v, tmp.df$dynamical_v)
		)),
		parse = FALSE
	) +
	annotate(
		geom = "text",
		x = .percent_range(tmp.df$v, 0.6),
		y = .percent_range(tmp.df$dynamical_v, 0.05),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("NRMSE: ", sprintf(
			"%.3f", nrmse(tmp.df$v, tmp.df$dynamical_v)
		)),
		parse = FALSE
	)

est_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$est,
		title = "Using learned k-NN",
		x_lab = "True time",
		y_lab = "Estimated time",
		color_by = "true_t",
		point.size = 5.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none")




panel_a_h.p <-
	plot_grid(
		spliced_t.p,
		ms_t.p,
		unspliced_t.p,
		mu_t.p,
		mums.p,
		dynamical_t.p,
		dynamical_v.p,
		est_t.p,
		nrow = 2,
		ncol = 4,
		label_size = 10,
		labels = "auto"
	)

panel_efg.p <- plot_grid(
	mums.p,
	dynamical_t.p,
	dynamical_v.p,
	nrow = 1,
	ncol = 3,
	label_size = 10,
	labels = "auto"
)

### use true k-NN
adata <-
	scv$datasets$simulation(
		n_obs = n_obs,
		n_vars = n_vars,
		t_max = 25,
		alpha = 5,
		beta = .3,
		gamma = .5,
		noise_level = noise_level
	)
ncomp <- as.integer(min(30L, n_vars - 1L))

### get true s
sce.o <- zellkonverter::AnnData2SCE(adata)
s.m <- t(sapply(seq_len(nrow(sce.o)), function(i) {
	dynamical_thet(
		alpha = rowData(sce.o)$true_alpha[i],
		beta = rowData(sce.o)$true_beta[i],
		gamma = rowData(sce.o)$true_gamma[i],
		scaling = rowData(sce.o)$true_scaling[i],
		t_ = rowData(sce.o)$true_t_[i],
		u0 = 0,
		s0 = 0,
		tmax = 25,
		new_t = sce.o$true_t
	)[, 2]
}))
dimnames(s.m) <- dimnames(sce.o)
assay(sce.o, "X") <- s.m
adata <- zellkonverter::SCE2AnnData(sce.o)

### learn k-NN graph using true S
scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)

### switch back and rerun PCA
sce.o <- zellkonverter::AnnData2SCE(adata)
reducedDim(sce.o, "true_pca") <- reducedDim(sce.o, "X_pca")
assay(sce.o, "X") <- assay(sce.o, "spliced")
adata <- zellkonverter::SCE2AnnData(sce.o)

scv$pp$pca(adata, n_comps = ncomp)
rotation.m <- adata$varm["PCs"]



### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)


###
i <- 1
tmp.m <-
	dynamical_thet(
		alpha = rowData(dynamical.o)$true_alpha[i],
		beta = rowData(dynamical.o)$true_beta[i],
		gamma = rowData(dynamical.o)$true_gamma[i],
		scaling = rowData(dynamical.o)$true_scaling[i],
		t_ = rowData(dynamical.o)$true_t_[i],
		u0 = 0,
		s0 = 0,
		tmax = 25,
		new_t = dynamical.o$true_t
	)
tmp.df <-
	data.frame(
		x = dynamical.o$true_t,
		s = tmp.m[, 2],
		v = tmp.m[, 4],
		u = tmp.m[, 1],
		est = assay(dynamical.o, "fit_t")[i,],
		ms = assay(dynamical.o, "Ms")[i,],
		mu = assay(dynamical.o, "Mu")[i,],
		spliced = assay(dynamical.o, "spliced")[i,],
		unspliced = assay(dynamical.o, "unspliced")[i,],
		steady_v = assay(steady.o, "velocity")[i, ],
		dynamical_v = assay(dynamical.o, "velocity")[i, ]
	)

mums.p <-
	plot_dynamicalmums(
		dynamical.o,
		rownm = as.character(i - 1),
		color_by = "true_t",
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		title = str_c("Using true k-NN"),
		color.name = "True t",
		point.size = 5.01
	) +
	geom_abline(
		slope = rowData(steady.o[i,])$velocity_gamma,
		intercept = 0,
		linetype = "dashed"
	) +
	theme(legend.position = "none") +
	xlab(TeX(r'(Ms$_g$)')) + ylab(TeX(r'(Mu$_g$)'))


ms_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$ms,
		title = "Using true k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Ms$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = s),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")

mu_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$mu,
		title = "Using true k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Mu$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = u),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")

spliced_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$spliced,
		title = "Using true k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Spliced$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = s),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(
		legend.position = c(1, 1),
		legend.justification = c(1, 1),
		legend.key.size = unit(6, "pt")
	)

unspliced_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$unspliced,
		title = "Using true k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Unspliced$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = u),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")



dynamical_t.p <-
	plotxy_zeoroy(
		tmp.df$x,
		tmp.df$dynamical_v,
		title = "Using true k-NN",
		x_lab = "True time",
		y_lab = TeX(r'(Dynamical v$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	geom_path(
		data = tmp.df,
		aes(x = x, y = v),
		alpha = 0.7,
		linetype = "solid",
		color = "black",
		size = 0.5
	) +
	theme(legend.position = "none")


# steady_v.p <- plotxy(tmp.df$v, tmp.df$steady_v, title = "Steady v.", x_lab = "True v.",
# 										 y_lab = "Steady v.", color_by = "true_t", point.size = 4.01,
# 										 colors.v = tmp.df$x,
# 										 hue.colors = c("#440154", "#482576", "#414487","#35608D",
# 										 							 "#2A788E", "#21908C","#22A884", "#43BF71",
# 										 							 "#7AD151", "#BBDF27", "#FDE725"), color.name = "t") +
# 	.theme_noframe +
# 	theme(plot.title = element_blank(), legend.position = "none",
# 				plot.margin = unit(c(2, 2, 3, 3), "pt"))

dynamical_v.p <-
	plotxy_zeoroy(
		tmp.df$v,
		tmp.df$dynamical_v,
		title = "Using true k-NN",
		x_lab = TeX(r'(True v$_g$)'),
		y_lab = TeX(r'(Dynamical v$_g$)'),
		color_by = "true_t",
		point.size = 4.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none") +
	annotate(
		geom = "text",
		x = .percent_range(tmp.df$v, 0.6),
		y = .percent_range(tmp.df$dynamical_v, 0.15),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("PCC: ", sprintf(
			"%.3f", cor(tmp.df$v, tmp.df$dynamical_v)
		)),
		parse = FALSE
	) +
	annotate(
		geom = "text",
		x = .percent_range(tmp.df$v, 0.6),
		y = .percent_range(tmp.df$dynamical_v, 0.05),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("NRMSE: ", sprintf(
			"%.3f", nrmse(tmp.df$v, tmp.df$dynamical_v)
		)),
		parse = FALSE
	)


est_t.p <-
	plotxy(
		tmp.df$x,
		tmp.df$est,
		title = "Using true k-NN",
		x_lab = "True time",
		y_lab = "Estimated time",
		color_by = "true_t",
		point.size = 5.01,
		colors.v = tmp.df$x,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none")




panel_i_p.p <-
	plot_grid(
		spliced_t.p,
		ms_t.p,
		unspliced_t.p,
		mu_t.p,
		mums.p,
		dynamical_t.p,
		dynamical_v.p,
		est_t.p,
		nrow = 2,
		ncol = 4,
		label_size = 10,
		labels = letters[9:16]
	)

panel_mno.p <- plot_grid(
	mums.p,
	dynamical_t.p,
	dynamical_v.p,
	nrow = 1,
	ncol = 3,
	label_size = 10,
	labels = letters[4:6]
)


mp <- plot_grid(
	panel_a_h.p,
	panel_i_p.p,
	nrow = 2,
	ncol = 1,
	label_size = 10,
	labels = NULL
)

save_plot(
	here::here("figs", "sfigs", "sfig.noise_genes.pdf"),
	mp,
	base_height = 2 * 1.2 / 1.4,
	base_width = 2 * 1.4 / 1.4,
	nrow = 4 * 0.8,
	ncol = 3,
	device = cairo_pdf
)


### add 4 panels for comparisons of high dim. vector length
noise_level = 3
n_obs = 500L
n_vars = 10L


adata <-
	scv$datasets$simulation(
		n_obs = n_obs,
		n_vars = n_vars,
		t_max = 25,
		alpha = 5,
		beta = .3,
		gamma = .5,
		noise_level = noise_level
	)
sce.o <- zellkonverter::AnnData2SCE(adata)
ncomp <- as.integer(min(30L, n_vars - 1L))
scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)


### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)

### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)


true_v.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
	tmp.m <-
		dynamical_thet(
			alpha = rowData(dynamical.o)$true_alpha[i],
			beta = rowData(dynamical.o)$true_beta[i],
			gamma = rowData(dynamical.o)$true_gamma[i],
			scaling = rowData(dynamical.o)$true_scaling[i],
			t_ = rowData(dynamical.o)$true_t_[i],
			u0 = 0,
			s0 = 0,
			tmax = 25,
			new_t = dynamical.o$true_t
		)
	return(tmp.m[, "vt"])
}))
dimnames(true_v.m) <- dimnames(dynamical.o)

true_s.m <- t(sapply(seq_len(nrow(dynamical.o)), function(i) {
	tmp.m <-
		dynamical_thet(
			alpha = rowData(dynamical.o)$true_alpha[i],
			beta = rowData(dynamical.o)$true_beta[i],
			gamma = rowData(dynamical.o)$true_gamma[i],
			scaling = rowData(dynamical.o)$true_scaling[i],
			t_ = rowData(dynamical.o)$true_t_[i],
			u0 = 0,
			s0 = 0,
			tmax = 25,
			new_t = dynamical.o$true_t
		)
	return(tmp.m[, "st"])
}))
dimnames(true_s.m) <- dimnames(dynamical.o)



row.idx <- which(!rowAnyNAs(assay(steady.o, "velocity")))
steady_len_est.v <-
	sqrt(colSums(assay(steady.o, "velocity")[row.idx,] ^ 2))
len_true.v <- sqrt(colSums(true_v.m[row.idx,] ^ 2))

row.idx <- which(!rowAnyNAs(assay(dynamical.o, "velocity")))
dynamical_len_est.v <-
	sqrt(colSums(assay(dynamical.o, "velocity")[row.idx,] ^ 2))


### use true k-NN
assay(sce.o, "X") <- true_s.m
adata <- zellkonverter::SCE2AnnData(sce.o)

### learn k-NN graph using true S
scv$pp$pca(adata, n_comps = ncomp)
scv$pp$neighbors(adata)
scv$pp$moments(adata)

### switch back
sce.o <- zellkonverter::AnnData2SCE(adata)
reducedDim(sce.o, "true_pca") <- reducedDim(sce.o, "X_pca")
assay(sce.o, "X") <- assay(sce.o, "spliced")
adata <- zellkonverter::SCE2AnnData(sce.o)


### get steady
scv$tl$velocity(adata, mode = "deterministic")
steady.o <- zellkonverter::AnnData2SCE(adata)

### get dynamical
scv$tl$recover_dynamics(adata)
scv$tl$velocity(adata, mode = "dynamical")
dynamical.o <- zellkonverter::AnnData2SCE(adata)


row.idx <- which(!rowAnyNAs(assay(steady.o, "velocity")))
steady_len_est_trueknn.v <-
	sqrt(colSums(assay(steady.o, "velocity")[row.idx,] ^ 2))

row.idx <- which(!rowAnyNAs(assay(dynamical.o, "velocity")))
dynamical_len_est_trueknn.v <-
	sqrt(colSums(assay(dynamical.o, "velocity")[row.idx,] ^ 2))

### 4 panels showing example
steady_leanred.p <-
	plotxy(
		len_true.v,
		steady_len_est.v,
		title = "Speed from steady model using learned k-NN",
		x_lab = Tex(r'(True speed)'),
		y_lab = "Estimated speed from steady model",
		color_by = "true_t",
		point.size = 5.01,
		colors.v = steady.o$true_t,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = c(1, 0),
				legend.justification = c(1, 0)) +
	annotate(
		geom = "text",
		x = .percent_range(len_true.v, 0.05),
		y = .percent_range(steady_len_est.v, 0.95),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("PCC: ", sprintf(
			"%.3f", cor(len_true.v, steady_len_est.v)
		)),
		parse = FALSE
	)

steady_true.p <-
	plotxy(
		len_true.v,
		steady_len_est_trueknn.v,
		title = "Speed from steady model using using true k-NN",
		x_lab = "True speed",
		y_lab = "Estimated speed from steady model",
		color_by = "true_t",
		point.size = 5.01,
		colors.v = steady.o$true_t,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none") +
	annotate(
		geom = "text",
		x = .percent_range(len_true.v, 0.05),
		y = .percent_range(steady_len_est.v, 0.95),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("PCC: ", sprintf(
			"%.3f", cor(len_true.v, steady_len_est_trueknn.v)
		)),
		parse = FALSE
	)



dynamical_leanred.p <-
	plotxy(
		len_true.v,
		dynamical_len_est.v,
		title = "Speed from dynamical model using learned k-NN",
		x_lab = "True speed",
		y_lab = "Estimated speed from dynamical model",
		color_by = "true_t",
		point.size = 5.01,
		colors.v = dynamical.o$true_t,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none") +
	annotate(
		geom = "text",
		x = .percent_range(len_true.v, 0.05),
		y = .percent_range(dynamical_len_est.v, 0.95),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("PCC: ", sprintf(
			"%.3f", cor(len_true.v, dynamical_len_est.v)
		)),
		parse = FALSE
	)

dynamical_true.p <-
	plotxy(
		len_true.v,
		dynamical_len_est_trueknn.v,
		title = "Speed from dynamical model using true k-NN",
		x_lab = "True speed",
		y_lab = "Estimated speed from dynamical model",
		color_by = "true_t",
		point.size = 5.01,
		colors.v = dynamical.o$true_t,
		hue.colors = c(
			"#440154",
			"#482576",
			"#414487",
			"#35608D",
			"#2A788E",
			"#21908C",
			"#22A884",
			"#43BF71",
			"#7AD151",
			"#BBDF27",
			"#FDE725"
		),
		color.name = "t"
	) +
	theme(legend.position = "none") +
	annotate(
		geom = "text",
		x = .percent_range(len_true.v, 0.05),
		y = .percent_range(dynamical_len_est.v, 0.95),
		size = 2,
		hjust = 0,
		vjust = 1,
		label = str_c("PCC: ", sprintf(
			"%.3f", cor(len_true.v, dynamical_len_est_trueknn.v)
		)),
		parse = FALSE
	)


row3.p <-
	plot_grid(
		dynamical_leanred.p + ggtitle("Using learned k-NN") + ylab("Estimated speed") + xlab("True speed"),
		dynamical_true.p + ggtitle("Using true k-NN") + ylab("Estimated speed") + xlab("True speed"),
		nrow = 1,
		ncol = 3,
		label_size = 10,
		labels = c("g", "h")
	)





mp <- plot_grid(
	panel_efg.p,
	panel_mno.p,
	row3.p,
	nrow = 3,
	ncol = 1,
	label_size = 10,
	labels = NULL
)

save_plot(
	here::here("figs", "main", "main.noise_genes.pdf"),
	mp,
	base_height = 2 * 1.2 / 1.4,
	base_width = 2 * 1.4 / 1.4,
	nrow = 3 * 0.9,
	ncol = 3,
	device = cairo_pdf
)




### Sfig for steady model

mp <- plot_grid(
	steady_leanred.p,
	steady_true.p,
	dynamical_leanred.p,
	dynamical_true.p,
	
	nrow = 2,
	ncol = 2,
	label_size = 10,
	labels = c("auto")
)

save_plot(
	here::here("figs", "sfigs", "sfig.length_com.pdf"),
	mp,
	base_height = 2 * 1.2 ,
	base_width = 2 * 1.4 / 1,
	nrow = 2,
	ncol = 2,
	device = cairo_pdf
)
