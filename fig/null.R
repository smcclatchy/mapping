x <- pull.geno(hyper)[,"D4Mit214"]
u <- runif(length(x), -0.075, 0.075)
y <- hyper$pheno[,1]

png(filename = "null.png", width = 600, height = 600)
plot(y ~ x, type="n", xlab="Genotype", ylab="Blood Pressure",
     xlim=c(0.5,2.5), xaxs="i", xaxt="n", main = "Null Hypothesis")
abline(v=1:2, lty=2, col="gray90")
points(x + 0.3 * u, y, col = "gray90")
axis(side=1, at=1:2, labels=c("AA","AB"))
me <- tapply(y, x, mean)
segments(x0 = 0.5, y0 = mean(y), x1 = 2.5, y1 = mean(y), lwd=2, col="blue")
segments((1:2)-0.15, mean(y), (1:2)+0.15, mean(y), lwd=3, col="green2")
points(1, 119.58, col = "red")
segments(x0 = 1, y0 = mean(y), y1 = 119.5, lwd=2, col="red")
segments(x0 = 1, y0 = mean(y), x1 = 1.66, y1 = mean(y), lwd=2, col="red")
segments(x0 = 1, y0 = 119.5, x1 = 1.66, y1 = 119.5, lwd=2, col="red")
segments(x0 = 1.66, y0 = mean(y), y1 = 119.5, lwd=2, col="red")
text(x = 1.3, y = 110, labels = "residual\nsquared")
text(x = 1.5, y = 90, labels = "no QTL exists")
dev.off()