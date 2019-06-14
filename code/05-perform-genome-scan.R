## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("05-")


## ----load_dependencies, include=FALSE------------------------------------
library(qtl2)
iron <- read_cross2(file = system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
Xcovar <- get_x_covar(iron)


## ----scan1---------------------------------------------------------------
out <- scan1(genoprobs = pr, pheno = iron$pheno, Xcovar=Xcovar)


## ----scan1_multicore, eval=FALSE-----------------------------------------
## out <- scan1(genoprobs = pr, pheno = iron$pheno, Xcovar=Xcovar, cores=4)


## ----head_scan-----------------------------------------------------------
head(out, n=10)


## ----plot_lod, eval=FALSE------------------------------------------------
## plot_scan1(out, map = map, lodcolumn = "liver")

