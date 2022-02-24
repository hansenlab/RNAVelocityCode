rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/chromaffin/dynamical.o.qs"))
colors.v <- levels(factor(sce.o$cell_type))
sce.o$cell_type <- fct_recode(factor(sce.o$cell_type), "SCPs" = "#0066FFFF", "chromaffin" = "#00FF66FF",
					 "sympathoblasts" = "#CC00FFFF", "non-labelled" = "#CCFF00FF", "bridge" = "#FF0000FF")


rotation.m <- attr(reducedDim(sce.o, "PCA"), "rotation")[, 1:2]
reducedDim(sce.o, "PCA") <- reducedDim(sce.o, "PCA")[, 1:2]
sce.o <- logNormCounts(sce.o, assay.type= "spliced")
data.m <- assay(sce.o, "logcounts")[match(rownames(rotation.m), rownames(sce.o)), ]
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

common_gene.v <- intersect(rownames(data.m), rownames(sce.o)[which(!rowAnyNAs(assay(sce.o, "velocity")))])
new_pca.m <- as.matrix(t(data.m[common_gene.v, ]) %*% rotation.m[common_gene.v, ])


reducedDim(sce.o, "new.PCA") <- new_pca.m



### add velocity
splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "spliced")
splusv.m[which(splusv.m < 0)] <- 0
data.m <- log2(splusv.m[match(rownames(rotation.m), rownames(sce.o)), ] + 1)

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])


plot(new_pca.m)

### plot
embedded <- embedVelocity(reducedDim(sce.o, "PCA"), sce.o)


grid.df <- gridVectors(reducedDim(sce.o, "PCA"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p1 <- plotEmbScat(sce.o, "PCA", title = "Vector fields", x_lab = "PCA-1", y_lab = "PCA-2", color_by = "cell_type", labels.v = levels(sce.o$cell_type),
									color.name = "Cell_type", colors.v = colors.v, point.size = 3.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


### projection 
embedded <- new_pca.m - reducedDim(sce.o, "PCA")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "PCA"))

grid.df <- gridVectors(reducedDim(sce.o, "PCA"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p2 <- plotEmbScat(sce.o, "PCA", title = "Vector fields (umap.predict)", x_lab = "PCA-1", y_lab = "PCA-2", color_by = "cell_type", labels.v = levels(sce.o$cell_type),
									color.name = "Cell_type", colors.v = colors.v, point.size = 3.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)




embedded <- embedVelocity(reducedDim(sce.o, "new.PCA"), sce.o)

grid.df <- gridVectors(reducedDim(sce.o, "new.PCA"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p3 <- plotEmbScat(sce.o, "new.PCA", title = "Vector fields on subset data ", x_lab = "PCA-1", y_lab = "PCA-2", color_by = "cell_type", labels.v = levels(sce.o$cell_type),
									color.name = "Cell_type", colors.v = colors.v, point.size = 3.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)


### projection 
embedded <- new_pca.m - reducedDim(sce.o, "new.PCA")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "new.PCA"))


grid.df <- gridVectors(reducedDim(sce.o, "new.PCA"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p4 <- plotEmbScat(sce.o, "new.PCA", title = "Vector fields (project on subset data)", x_lab = "PCA-1", y_lab = "PCA-2", color_by = "cell_type", labels.v = levels(sce.o$cell_type),
									color.name = "Cell_type", colors.v = colors.v, point.size = 3.01) + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.02, "inches")), inherit.aes = FALSE, size = 0.2)





mp <- wrap_plots(list(p1, p2, p3, p4), ncol = 2)

save_plot(here::here("figs", "chromaffin", str_c("pca_predict.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow =  2*1.2, ncol = 2*1.2 , device = cairo_pdf)




