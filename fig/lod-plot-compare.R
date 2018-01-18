library(qtl2)
iron <- read_cross2(file = system.file("extdata", 
                                       "iron.zip", 
                                       package="qtl2geno") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)
Xcovar <- get_x_covar(iron)
out <- scan1(genoprobs = pr, pheno = iron$pheno, Xcovar=Xcovar)
out_pg <- scan1(pr, iron$pheno, kinship, Xcovar=Xcovar, cores=4)
out_pg_loco <- scan1(pr, iron$pheno, 
                     kinship_loco, Xcovar=Xcovar)

color <- c("slateblue", "violetred", "green3")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(maxlod(out), maxlod(out_pg), maxlod(out_pg_loco))

png(filename = "lod-plot-compare-liver.png", width = 600, height = 600)
plot(out, map, lodcolumn=1, 
     col=color[1], 
     main=colnames(iron$pheno)[1],
     ylim=c(0, ymx*1.02))
plot(out_pg, map, lodcolumn=1, col=color[2], add=TRUE)
plot(out_pg_loco, map, lodcolumn=1,
     col=color[3], add=TRUE, lty=2)
legend("topleft", lwd=2, col=color, 
       c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
dev.off()

png(filename = "lod-plot-compare-spleen.png", width = 600, height = 600)
plot(out, map, lodcolumn=2, 
     col=color[1], 
     main=colnames(iron$pheno)[2],
     ylim=c(0, ymx*1.02))
plot(out_pg, map, lodcolumn=2, col=color[2], add=TRUE)
plot(out_pg_loco, map, lodcolumn=2,
     col=color[3], add=TRUE, lty=2)
legend("topleft", lwd=2, col=color, 
       c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
dev.off()