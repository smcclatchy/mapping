x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]

png(filename = "residual.png", width = 600, height = 600)
plot(y ~ x, type="n", xlab="Genotype", ylab="Blood Pressure",
     xlim=c(0, 3), xaxs="i", xaxt="n", main = "Line of Best Fit")
abline(v=1:2, lty=2, col="gray90")
points(x + 0.3 * u, y, col = "gray90")
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments((1:2)-0.15, me, (1:2)+0.15, me, lwd=3, col="orange1")
abline(lm(y~x)$coef, col="blue", lwd=2)
points(1, 119.58, col = "red")
segments(x0 = 1, y0 = 104.4, y1 = 119.5, lwd=2, col="red")
text(x = 0.8, y = 112, labels = "error or\n residual")
dev.off()