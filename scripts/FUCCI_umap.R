rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/dynamical.o.qs"))
umap.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/umap.o.qs"))
rotation.m <- attr(reducedDim(sce.o, "PCA"), "rotation")[, 1:30]
sce.o <- logNormCounts(sce.o, assay.type= "spliced")
data.m <- assay(sce.o, "logcounts")[match(rownames(rotation.m), rownames(sce.o)), ]
row_m.v <- rowMeans(data.m)
data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")

common_gene.v <- intersect(rownames(data.m), rownames(sce.o)[which(!rowAnyNAs(assay(sce.o, "velocity")))])
new_pca.m <- as.matrix(t(data.m[common_gene.v, ]) %*% rotation.m[common_gene.v, ])

library(umap)
new_umap.m <- predict(umap.o, new_pca.m)[, c(2, 1)]
plot(new_umap.m)

reducedDim(sce.o, "new.UMAP") <- new_umap.m



### add velocity

splusv.m <- assay(sce.o, "velocity") + assay(sce.o, "spliced")
splusv.m[which(splusv.m < 0)] <- 0
data.m <- log2(splusv.m[match(rownames(rotation.m), rownames(sce.o)), ] + 1)

data.m <- sweep(data.m, MARGIN = 1, row_m.v, FUN = "-")
idx <- which(!rowAnyNAs(as.matrix(data.m)))
new_pca.m <- as.matrix(t(data.m[idx, ]) %*% rotation.m[idx, ])


new_umap.m <- predict(umap.o, new_pca.m)[, c(2, 1)]
plot(new_umap.m)

### plot

embedded <- embedVelocity(reducedDim(sce.o, "UMAP"), sce.o)
plotVelocityStream(sce.o, embedded = embedded, use.dimred = "UMAP")

grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p1 <- plotEmbFucci(sce.o, "UMAP", title = "Vector fields", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)


### projection 
embedded <- new_umap.m - reducedDim(sce.o, "UMAP")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "UMAP"))
plotVelocityStream(sce.o, embedded = embedded, use.dimred = "UMAP")

grid.df <- gridVectors(reducedDim(sce.o, "UMAP"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p2 <-  plotEmbFucci(sce.o, "UMAP", title = "Vector fields (umap.predict)", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)





embedded <- embedVelocity(reducedDim(sce.o, "new.UMAP"), sce.o)
plotVelocityStream(sce.o, embedded = embedded, use.dimred = "new.UMAP")

grid.df <- gridVectors(reducedDim(sce.o, "new.UMAP"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p3 <- plotEmbFucci(sce.o, "new.UMAP", title = "Vector fields on subset data", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)


### projection 
embedded <- new_umap.m - reducedDim(sce.o, "new.UMAP")
dimnames(embedded) <- dimnames(reducedDim(sce.o, "new.UMAP"))
plotVelocityStream(sce.o, embedded = embedded, use.dimred = "new.UMAP")

grid.df <- gridVectors(reducedDim(sce.o, "new.UMAP"), embedded)
colnames(grid.df) <- c("start.1", "start.2", "end.1", "end.2")
p4 <- plotEmbFucci(sce.o, "new.UMAP", title = "Vector fields (umap.predict on subset data)", x_lab = "UMAP-1", y_lab = "UMAP-2") + 
	geom_segment(data=grid.df, mapping=aes(x=start.1, y=start.2, xend=end.1, yend=end.2), arrow = arrow(length=unit(0.01, "inches")), inherit.aes = FALSE, size = 0.2)






mp <- wrap_plots(list(p1, p2, p3, p4), ncol = 2)

save_plot(here::here("figs", "FUCCI", str_c("umap_predict.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 1.5 * 2, ncol = 2 * 1.5, device = cairo_pdf)








