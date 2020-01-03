library(qtl2)
iron <- read_cross2(file = system.file("extdata",
                                       "iron.zip",
                                       package = "qtl2"))
g <- maxmarg(pr, map, chr=2, pos=56.8, return_char=TRUE)
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])

png(filename = "residual.png", width = 600, height = 600)

plot_pxg(g, iron$pheno[,"liver"], ylab="Liver phenotype", 
         main = "Line of Best Fit",
         seg_col = "#d7191c", 
         seg_lwd = 4,
         bg = "#f7f7f7",
         lwd = 1,
         col = "#969696")

abline(a = (c2eff["D2Mit17",4] + (2 * c2eff["D2Mit17",3])), 
       b = c2eff["D2Mit17", 1],
       col = "#d7191c")
points(2, 157, lwd = 3)
segments(x0 = 2, y0 = 95, y1 = 155,
         lwd=2)
text(x = 1.7, y = 160, labels = "error or\n residual", cex=0.8)
segments(x0 = 1.8, y0 = 159, x1 = 2, y1 = 151, 
         lwd=1)
dev.off()