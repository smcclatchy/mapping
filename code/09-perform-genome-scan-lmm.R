## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("09-")


## ----scan1_pg, eval=FALSE------------------------------------------------
## out_pg <- scan1(pr, iron$pheno, kinship=kinship, Xcovar=Xcovar)


## ----scan1_pg_multicore, eval=FALSE--------------------------------------
## out_pg <- scan1(pr, iron$pheno, kinship, Xcovar=Xcovar, cores=4)


## ----calc_kinship_loco, eval=FALSE---------------------------------------
## kinship_loco <- calc_kinship(pr, "loco")


## ----scan1_pg_loco, eval=FALSE-------------------------------------------
## out_pg_loco <- scan1(pr, iron$pheno, kinship_loco, Xcovar=Xcovar)


## ---- eval=FALSE---------------------------------------------------------
## plot_scan1(out_pg_loco, map = map, lodcolumn = "liver", col = "black")
## plot_scan1(out_pg, map = map, lodcolumn = "liver", col = "blue", add = TRUE)
## plot_scan1(out, map = map, lodcolumn = "liver", add = TRUE, col = "green")


## ------------------------------------------------------------------------
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/B6BTBR/b6btbr.zip")
b6btbr <- read_cross2(file)

