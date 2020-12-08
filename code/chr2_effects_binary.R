# load library, data, and create dependencies
library(qtl2)
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map=iron$gmap, step=1)
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)

# chromosome 2 coefficients for binary model
c2eff_bin <- scan1coef(pr[,"2"], bin_pheno[,"liver"], model="binary")

# output figure
png(filename = "./fig/chr2_effects_binary.png")
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff_bin, map["2"], columns=1:3, col=col)
last_coef <- unclass(c2eff_bin)[nrow(c2eff_bin),] # pull out last coefficients
for(i in seq(along=last_coef))
  axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
dev.off()