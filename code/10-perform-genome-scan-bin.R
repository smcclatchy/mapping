## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("10-")


## ----threshold_phenotypes, eval=FALSE------------------------------------
## bin_pheno <- apply(iron$pheno, 2, function(a) as.numeric(a > median(a)))
## rownames(bin_pheno) <- rownames(iron$pheno)


## ----binary_trait_scan, eval=FALSE---------------------------------------
## out_bin <- scan1(pr, bin_pheno, Xcovar=Xcovar, model="binary")


## ----plot_bin_scan, eval=FALSE-------------------------------------------
## par(mar=c(5.1, 4.1, 1.1, 1.1))
## ymx <- maxlod(out_bin)
## plot(out_bin, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
## plot(out_bin, map, lodcolumn=2, col="violetred", add=TRUE)
## legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out_bin), bg="gray90")


## ----find_peaks_bin_scan, eval=FALSE-------------------------------------
## find_peaks(out_bin, map, threshold=3.5, drop=1.5)

