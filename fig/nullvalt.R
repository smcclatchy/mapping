x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]

png(filename = "nullvalt.png", width = 600, height = 600)
plot(y ~ x, type="n", xlab="Genotype", ylab="Blood Pressure",
     xlim=c(0.5,2.5), xaxs="i", xaxt="n", main = "Null and Alternative\nHypotheses")
abline(v=1:2, lty=2, col="gray90")
points(x + 0.3 * u, y, col = "gray90")
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments(x0 = 0.5, y0 = mean(y), x1 = 2.5, y1 = mean(y), lwd=2, col="purple4")
segments((1:2)-0.15, mean(y), (1:2)+0.15, mean(y), lwd=3, col="green2")
abline(lm(y~x)$coef, col="sienna1", lwd=2)
segments((1:2)-0.15, me, (1:2)+0.15, me, lwd=3, col="orange1")
text(x = 0.8, y = 98, labels = "no QTL exists", cex=0.7)
text(x = 0.8, y = 110, labels = "QTL exists\nhere", cex=0.7)
text(x = 2.35, y = 99.5, labels = "mean\nphenotype", cex=0.7)
dev.off()