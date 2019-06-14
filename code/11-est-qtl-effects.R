## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("11-")


## ----load_dependencies, echo=FALSE---------------------------------------
library(qtl2)
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
kinship <- calc_kinship(pr)
kinship_loco <- calc_kinship(pr, "loco")
bin_pheno <- apply(iron$pheno, 2, function(a) as.numeric(a > median(a)))
rownames(bin_pheno) <- rownames(iron$pheno)


## ----est_effects_liver_c2------------------------------------------------
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])


## ----effects_matrix_dim--------------------------------------------------
dim(c2eff)
head(c2eff)


## ----plot_effects_liver_c2-----------------------------------------------
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["2"], columns=1:3, col=col)
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


## ----est_effects_liver_c2_contr------------------------------------------
c2effB <- scan1coef(pr[,"2"], iron$pheno[,"liver"],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))



## ----add_dom_matrix_dim--------------------------------------------------
dim(c2effB)
head(c2effB)


## ----plot_effects_liver_c2_contr-----------------------------------------
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB, map["2"], columns=2:3, col=col)
last_coef <- unclass(c2effB)[nrow(c2effB),2:3] # last two coefficients
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


## ----est_effects_pg_liver_c2---------------------------------------------
c2eff_pg <- scan1coef(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]])
dim(c2eff_pg)
head(c2eff_pg)


## ----plot_effects_pg_liver_c2--------------------------------------------
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff_pg, map["2"], columns=1:3, col=col, ylab="Phenotype average")
last_coef <- unclass(c2eff_pg)[nrow(c2eff_pg),]
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


## ----est_effects_pg_liver_c2_contr---------------------------------------
c2effB_pg <- scan1coef(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]],
                       contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))


## ----plot_effects_pg_liver_c2_contr--------------------------------------
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB_pg, map["2"], columns=2:3, col=col)
last_coef <- unclass(c2effB_pg)[nrow(c2effB_pg),2:3]
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


## ----scan1blup-----------------------------------------------------------
c2blup <- scan1blup(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]])


## ----plot_scan1blup------------------------------------------------------
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
ylim <- range(c(c2blup, c2eff))+c(-1,1)
plot(c2eff, map["2"], columns=1:3, col=col, ylab="Phenotype average", ylim=ylim,
     xlab="Chr 2 position")
plot(c2blup, map["2"], columns=1:3, col=col, add=TRUE, lty=2)
last_coef <- unclass(c2eff)[nrow(c2eff),]
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


## ----scan1coef_binary----------------------------------------------------
c2eff_bin <- scan1coef(pr[,"2"], bin_pheno[,"liver"], model="binary")


## ----plot_effects_binary-------------------------------------------------
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff_bin, map["2"], columns=1:3, col=col)
last_coef <- unclass(c2eff_bin)[nrow(c2eff_bin),] # pull out last coefficients
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])


## ----get_inferred_genotypes----------------------------------------------
g <- maxmarg(pr, map, chr=2, pos=28.6, return_char=TRUE)


## ----plot_pheno_geno-----------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, iron$pheno[,"liver"], ylab="Liver phenotype")


## ----plot_pheno_geno_se--------------------------------------------------
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, iron$pheno[,"liver"], SEmult=2, swap_axes=TRUE, xlab="Liver phenotype")

