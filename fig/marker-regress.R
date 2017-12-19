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
abline(lm(y~x)$coef, col="blue", lwd=2)
png(filename = "marker-regress.png", width = 610, height = 600)
dev.off()