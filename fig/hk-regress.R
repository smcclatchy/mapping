x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]

plot(y ~ x, type="n", xlab="Genotype", ylab="Phenotype",
     xlim=c(0.5,2.5), xaxs="i", xaxt="n")
abline(v=1:2, lty=2, col="gray50")
points(x+u, y)
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments((1:2)-0.15, me, (1:2)+0.15, me, lwd=3, col="green2")
fakey <- rnorm(n=125, mean = mean(y), sd=sd(y))
points(1.5+(4*u[1:125]), fakey, col = "gray70")
segments((1.5)-0.15, mean(y), (1.5)+0.15, mean(y), lwd=3, col="green2")
mtext("P(AA) = 1\nP(AB) = 0", at = 1, cex = 0.6)
mtext("P(AA) = 0\nP(AB) = 1", at = 2, cex = 0.6)
mtext("P(AA) = P(AB)\n      = 0.5", at = 1.5, cex = 0.6)
segments(x0=1.5, y0= 120.5, y1=130, lty=2, col="gray50")
points(x = 1.4, y = 120, col = "red")
abline(lm(y~x)$coef, col="blue", lwd=2)

png(filename = "hk-regress.png", width = 610, height = 455)
dev.off()