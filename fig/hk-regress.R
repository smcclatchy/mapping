x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]

png(filename = "hk-regress.png", width = 600, height = 600)
plot(y ~ x, type="n", xlab="Genotype", ylab="Blood Pressure",
     xlim=c(0.5,2.5), xaxs="i", xaxt="n", main = "Haley-Knott Regression")
abline(v=1:2, lty=2, col="gray50")
points(x+0.3*u, y)
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments((1:2)-0.15, me, (1:2)+0.15, me, lwd=3, col="orange1")
fakey <- rnorm(n=125, mean = mean(y), sd=sd(y))
points(1.5+(4*u[1:125]), fakey, col = "gray70")
segments((1.5)-0.15, mean(y), (1.5)+0.15, mean(y), lwd=3, col="orange1")
mtext("P(AA) = 1\nP(AB) = 0", side = 1, line = -1, at = .85, cex = 0.7)
mtext("P(AA) = 0\nP(AB) = 1", side = 1, line = -1, at = 2.15, cex = 0.7)
mtext("P(AA) = P(AB)\n      = 0.5", side = 1, line = 0.4, at = 1.5, cex = 0.7)
abline(v=1.5, lty=2, col="gray50")
segments(x0 = seq(from = 1, to = 2, by = 0.1), y0= 80, y1=83)
points(x = 1.4, y = 85, col = "green", lwd=2)
abline(lm(y~x)$coef, col="blue", lwd=2)
dev.off()