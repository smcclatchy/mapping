## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("07-")
getwd()


## ---- include=FALSE------------------------------------------------------
library(qtl2)
iron <- read_cross2(file = system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
Xcovar <- get_x_covar(iron)
out <- scan1(genoprobs = pr, pheno = iron$pheno, Xcovar=Xcovar)


## ----find_peaks----------------------------------------------------------
operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000)
thr = summary(operm)
find_peaks(scan1_output = out, map = map, threshold = thr, prob = 0.95, expand2markers = FALSE)


## ----find_multiple_peaks-------------------------------------------------
find_peaks(scan1_output = out, map = map, threshold = thr, peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)

