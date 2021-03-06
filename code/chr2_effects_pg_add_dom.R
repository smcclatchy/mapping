# load library, data, and create dependencies
library(qtl2)
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)

# chromosome 2 coefficients - additive and dominant effects
c2effB_pg <- scan1coef(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]],
                       contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))

# output figure
png(filename = "./fig/chr2_effects_pg_add_dom.png")
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB_pg, map["2"], columns=2:3, col=col)
last_coef <- unclass(c2effB_pg)[nrow(c2effB_pg),2:3]
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
dev.off()