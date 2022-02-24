
rm(list=ls())
source(here::here("scripts/utils.R"))

sce.o <- qread(here::here("data/sourcedata/Mahdessian_v1.2/steady.o.qs"))

### run PCA on unspliced counts
sce.o <- logNormCounts(sce.o, assay.type= "unspliced")

### Do PCA on spliced counts using HVG
dec <- modelGeneVar(sce.o)
top.hvgs <- getTopHVGs(dec, n=2000)

set.seed(100)
sce.o <- runPCA(sce.o, subset_row = top.hvgs, name = "PCA.u")

library(velociraptor)
un_knn.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA.u")
plotPCA(sce.o, colour_by="velocity_pseudotime")
plotReducedDim(un_knn.o, dimred = "X_pca",colour_by="velocity_pseudotime")

### replace with new Mu
# assay(sce.o, "Mu") <- assay(un_knn.o, "Mu")

g <- "TOP2A"
p1.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(sce.o, "Mu")[g, ], title = g, x_lab = "Ms", y_lab = "Mu (KNN by s)", color_by = "fucci_time",color.name = "FUCCI")
p2.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(un_knn.o, "Mu")[g, ], title = g, x_lab = "Ms", y_lab = "Mu (KNN by u)", color_by = "fucci_time",color.name = "FUCCI")
p3.p <- plotxy(assay(sce.o, "Mu")[g, ], assay(un_knn.o, "Mu")[g, ], title = g, x_lab = "Mu (KNN by s)", y_lab = "Mu (KNN by u)", color_by = "fucci_time",color.name = "FUCCI")

g <- "CDK1"
p4.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(sce.o, "Mu")[g, ], title = g, x_lab = "Ms", y_lab = "Mu (KNN by s)", color_by = "fucci_time",color.name = "FUCCI")
p5.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(un_knn.o, "Mu")[g, ], title = g, x_lab = "Ms", y_lab = "Mu (KNN by u)", color_by = "fucci_time",color.name = "FUCCI")
p6.p <- plotxy(assay(sce.o, "Mu")[g, ], assay(un_knn.o, "Mu")[g, ], title = g, x_lab = "Mu (KNN by s)", y_lab = "Mu (KNN by u)", color_by = "fucci_time",color.name = "FUCCI")


mp <- wrap_plots(list(p1.p, p4.p, p2.p, p5.p, p3.p, p6.p), ncol = 2)

save_plot(here::here("figs", "FUCCI", str_c( "KNN_by_u.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 3, ncol = 2, device = cairo_pdf)








###  use 60 neighbors

knn60.o <- scvelo(sce.o, assay.X="spliced", mode =  "deterministic", use.dimred = "PCA", scvelo.params= list(moments = list(n_neighbors=60L)))

g <- "TOP2A"
p1.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(sce.o, "Mu")[g, ], title = g, x_lab = "Ms (KNN n=30)", y_lab = "Mu (KNN n=30)", color_by = "fucci_time",color.name = "FUCCI")
p2.p <- plotxy(assay(knn60.o, "Ms")[g, ], assay(knn60.o, "Mu")[g, ], title = g, x_lab = "Ms (KNN n=60)", y_lab = "Mu (KNN n=60)", color_by = "fucci_time",color.name = "FUCCI")
p3.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(knn60.o, "Ms")[g, ], title = g, x_lab = "Ms (KNN n=30)", y_lab = "Ms (KNN n=60)", color_by = "fucci_time",color.name = "FUCCI")
p4.p <- plotxy(assay(sce.o, "Mu")[g, ], assay(knn60.o, "Mu")[g, ], title = g, x_lab = "Mu (KNN n=30)", y_lab = "Mu (KNN n=60)", color_by = "fucci_time",color.name = "FUCCI")


g <- "CDK1"
p5.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(sce.o, "Mu")[g, ], title = g, x_lab = "Ms (KNN n=30)", y_lab = "Mu (KNN n=30)", color_by = "fucci_time",color.name = "FUCCI")
p6.p <- plotxy(assay(knn60.o, "Ms")[g, ], assay(knn60.o, "Mu")[g, ], title = g, x_lab = "Ms (KNN n=60)", y_lab = "Mu (KNN n=60)", color_by = "fucci_time",color.name = "FUCCI")
p7.p <- plotxy(assay(sce.o, "Ms")[g, ], assay(knn60.o, "Ms")[g, ], title = g, x_lab = "Ms (KNN n=30)", y_lab = "Ms (KNN n=60)", color_by = "fucci_time",color.name = "FUCCI")
p8.p <- plotxy(assay(sce.o, "Mu")[g, ], assay(knn60.o, "Mu")[g, ], title = g, x_lab = "Mu (KNN n=30)", y_lab = "Mu (KNN n=60)", color_by = "fucci_time",color.name = "FUCCI")

mp <- wrap_plots(list(p1.p, p5.p, p2.p, p6.p, p3.p, p7.p, p4.p, p8.p), ncol = 2)

save_plot(here::here("figs", "FUCCI", str_c( "KNN_by_60.pdf")), mp,
					base_height = 2, base_width = 2*1.4, nrow = 4, ncol = 2, device = cairo_pdf)






