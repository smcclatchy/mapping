# load library, data, and create dependencies
library(qtl2)
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)

# chromosome 2 contrasts with estimates of mu, a, and d
c2effB <- scan1coef(pr[,"2"], iron$pheno[,"liver"],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))

# output figure
png(filename = "./fig/chr2_effects_contrasts.png")
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB, map["2"], columns=2:3, col=col)
last_coef <- unclass(c2effB)[nrow(c2effB),2:3] # last two coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
dev.off()