library(qtl)
data(hyper)
x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]

png(filename = "model-equation.png", width = 600, height = 600)
plot(y ~ x, type="n", xlab="Genotype", ylab="Blood Pressure",
     xlim=c(0, 3), xaxs="i", xaxt="n", main = "Null and Alternative\nHypotheses")
abline(v=1:2, lty=2, col="gray90")
points(x + 0.3 * u, y, col = "gray90")
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments(x0 = 0, y0 = mean(y), x1 = 3, y1 = mean(y), lwd=2, col="purple4")
segments((1:2)-0.15, mean(y), (1:2)+0.15, mean(y), lwd=3, col="green2")
abline(lm(y~x)$coef, col="sienna1", lwd=2)
segments((1:2)-0.15, me, (1:2)+0.15, me, lwd=3, col="orange1")
text(x = 0.4, y = 99.5, labels = "y = 101.6 + error", cex=0.9)
text(x = 0.5, y = 111, labels = "y = 110.2 + -5.8X + error", cex=0.9)
points(1, 119.58, col = "red")
segments(x0 = 1, y0 = 104.4, y1 = 119.5, lwd=2, col="red")
text(x = 1.15, y = 117, labels = "error or\n residual", cex=0.9)
dev.off()