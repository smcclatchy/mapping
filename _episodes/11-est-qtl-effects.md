---
title: "Estimated QTL effects"
teaching: 10
exercises: 20
questions:
- "How do I find the estimated effects of a QTL on a phenotype?"
objectives:
- Obtain estimated QTL effects.
- Plot estimated QTL effects.
keypoints:
- "."
- "."
source: Rmd
---




~~~
Error in insert_pseudomarkers(map = iron$gmap, step = 1): object 'iron' not found
~~~
{: .error}



~~~
Error in "cross2" %in% class(x): object 'iron' not found
~~~
{: .error}



~~~
Error in calc_kinship(pr): object 'pr' not found
~~~
{: .error}

The `scan1()` function returns only LOD scores. To obtain estimated QTL effects, use the function `scan1coef()`. This function takes a single phenotype and the genotype probabilities for a single chromosome and returns a matrix with the estimated coefficients at each putative QTL location along the chromosome.

For example, to get the estimated effects on chromosome 2 for the liver phenotype, we'd do the following:


~~~
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])
~~~
{: .r}



~~~
Error in scan1coef(pr[, "2"], iron$pheno[, "liver"]): object 'pr' not found
~~~
{: .error}

The result is a matrix,  positions &times;  genotypes. To plot the effects, use the function `plot_coef()`. There is again an S3 method function `plot.scan1coef()`, so one can just type `plot()`, but in either case you need to provide the map for the relevant chromosome. Use the argument columns to indicate which coefficient columns to plot.


~~~
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff, map["2"], columns=1:3, col=col)
~~~
{: .r}



~~~
Error in plot(c2eff, map["2"], columns = 1:3, col = col): object 'c2eff' not found
~~~
{: .error}



~~~
last_coef <- unclass(c2eff)[nrow(c2eff),] # pull out last coefficients
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2eff' not found
~~~
{: .error}



~~~
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
~~~
{: .r}



~~~
Error in seq(along = last_coef): object 'last_coef' not found
~~~
{: .error}

The default is to provide phenotype averages for each genotype group. If instead you want additive and dominance effects, you can provide a square matrix of _contrasts_, as follows:


~~~
c2effB <- scan1coef(pr[,"2"], iron$pheno[,"liver"],
                    contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))
~~~
{: .r}



~~~
Error in scan1coef(pr[, "2"], iron$pheno[, "liver"], contrasts = cbind(mu = c(1, : object 'pr' not found
~~~
{: .error}

The result will then contain the estimates of `mu`, `a`, and `d`. Here's a plot of the additive and dominance effects, which are in the
second and third columns.


~~~
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB, map["2"], columns=2:3, col=col)
~~~
{: .r}



~~~
Error in plot(c2effB, map["2"], columns = 2:3, col = col): object 'c2effB' not found
~~~
{: .error}



~~~
last_coef <- unclass(c2effB)[nrow(c2effB),2:3] # last two coefficients
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2effB' not found
~~~
{: .error}



~~~
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
~~~
{: .r}



~~~
Error in seq(along = last_coef): object 'last_coef' not found
~~~
{: .error}

If you provide a kinship matrix to `scan1coef()`, it fits a linear mixed model (LMM) to account for a residual polygenic effect. Here let's use the kinship matrix from the LOCO method.


~~~
c2eff_pg <- scan1coef(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]])
~~~
{: .r}



~~~
Error in scan1coef(pr[, "2"], iron$pheno[, "liver"], kinship_loco[["2"]]): object 'pr' not found
~~~
{: .error}

Here's a plot of the estimates.


~~~
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff_pg, map["2"], columns=1:3, col=col, ylab="Phenotype average")
~~~
{: .r}



~~~
Error in plot(c2eff_pg, map["2"], columns = 1:3, col = col, ylab = "Phenotype average"): object 'c2eff_pg' not found
~~~
{: .error}



~~~
last_coef <- unclass(c2eff_pg)[nrow(c2eff_pg),]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2eff_pg' not found
~~~
{: .error}



~~~
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
~~~
{: .r}



~~~
Error in seq(along = last_coef): object 'last_coef' not found
~~~
{: .error}

You can also get estimated additive and dominance effects, using a matrix of contrasts.


~~~
c2effB_pg <- scan1coef(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]],
                       contrasts=cbind(mu=c(1,1,1), a=c(-1, 0, 1), d=c(-0.5, 1, -0.5)))
~~~
{: .r}



~~~
Error in scan1coef(pr[, "2"], iron$pheno[, "liver"], kinship_loco[["2"]], : object 'pr' not found
~~~
{: .error}

Here's a plot of the results.


~~~
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2effB_pg, map["2"], columns=2:3, col=col)
~~~
{: .r}



~~~
Error in plot(c2effB_pg, map["2"], columns = 2:3, col = col): object 'c2effB_pg' not found
~~~
{: .error}



~~~
last_coef <- unclass(c2effB_pg)[nrow(c2effB_pg),2:3]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2effB_pg' not found
~~~
{: .error}



~~~
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
~~~
{: .r}



~~~
Error in seq(along = last_coef): object 'last_coef' not found
~~~
{: .error}

Another option for estimating the QTL effects is to treat them as random effects and calculate Best Linear Unbiased Predictors (BLUPs). This is particularly valuable for multi-parent populations such as the Collaborative Cross and Diversity Outbred mice, where the large number of possible genotypes at a QTL lead to considerable variability in the effect estimates. To calculate BLUPs, use `scan1blup()`; it takes the same arguments as `scan1coef()`, including
the option of a kinship matrix to account for a residual polygenic effect.


~~~
c2blup <- scan1blup(pr[,"2"], iron$pheno[,"liver"], kinship_loco[["2"]])
~~~
{: .r}



~~~
Error in scan1blup(pr[, "2"], iron$pheno[, "liver"], kinship_loco[["2"]]): object 'pr' not found
~~~
{: .error}

Here is a plot of the BLUPs (as dashed curves) alongside the standard estimates. Note that the BLUPs are centered at 0, while the coefficient estimates are centered at the phenotype average.


~~~
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
ylim <- range(c(c2blup, c2eff))+c(-1,1)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2blup' not found
~~~
{: .error}



~~~
plot(c2eff, map["2"], columns=1:3, col=col, ylab="Phenotype average", ylim=ylim,
     xlab="Chr 2 position")
~~~
{: .r}



~~~
Error in plot(c2eff, map["2"], columns = 1:3, col = col, ylab = "Phenotype average", : object 'c2eff' not found
~~~
{: .error}



~~~
plot(c2blup, map["2"], columns=1:3, col=col, add=TRUE, lty=2)
~~~
{: .r}



~~~
Error in plot(c2blup, map["2"], columns = 1:3, col = col, add = TRUE, : object 'c2blup' not found
~~~
{: .error}



~~~
last_coef <- unclass(c2eff)[nrow(c2eff),]
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2eff' not found
~~~
{: .error}



~~~
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
~~~
{: .r}



~~~
Error in seq(along = last_coef): object 'last_coef' not found
~~~
{: .error}

The `scan1coef` function can also provide estimated QTL effects for binary traits, with `model="binary"`. (However, `scan1blup` has not yet been implemented for binary traits.)


~~~
c2eff_bin <- scan1coef(pr[,"2"], bin_pheno[,"liver"], model="binary")
~~~
{: .r}



~~~
Error in scan1coef(pr[, "2"], bin_pheno[, "liver"], model = "binary"): object 'pr' not found
~~~
{: .error}

Here's a plot of the effects. They're a bit tricky to interpret, as they are basically log odds ratios.


~~~
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
col <- c("slateblue", "violetred", "green3")
plot(c2eff_bin, map["2"], columns=1:3, col=col)
~~~
{: .r}



~~~
Error in plot(c2eff_bin, map["2"], columns = 1:3, col = col): object 'c2eff_bin' not found
~~~
{: .error}



~~~
last_coef <- unclass(c2eff_bin)[nrow(c2eff_bin),] # pull out last coefficients
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'c2eff_bin' not found
~~~
{: .error}



~~~
for(i in seq(along=last_coef))
    axis(side=4, at=last_coef[i], names(last_coef)[i], tick=FALSE, col.axis=col[i])
~~~
{: .r}



~~~
Error in seq(along = last_coef): object 'last_coef' not found
~~~
{: .error}

Finally, to plot the raw phenotypes against the genotypes at a single putative QTL position, you can use the function plot_pxg(). This takes a vector of genotypes as produced by the maxmarg() function, which picks the most likely genotype from a set of genotype probabilities, provided it is greater than some specified value (the argument minprob). Note that the “marg” in “maxmarg” stands for “marginal”, as this function is selecting the genotype at each position that has maximum marginal probability.

For example, we could get inferred genotypes at the chr 2 QTL for the liver phenotype (at 28.6 cM) as follows:


~~~
g <- maxmarg(pr, map, chr=2, pos=28.6, return_char=TRUE)
~~~
{: .r}



~~~
Error in maxmarg(pr, map, chr = 2, pos = 28.6, return_char = TRUE): object 'pr' not found
~~~
{: .error}

We use `return_char=TRUE` to have `maxmarg()` return a vector of character strings with the genotype labels.

We then plot the liver phenotype against these genotypes as follows:


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, iron$pheno[,"liver"], ylab="Liver phenotype")
~~~
{: .r}



~~~
Error in plot_pxg(g, iron$pheno[, "liver"], ylab = "Liver phenotype"): object 'g' not found
~~~
{: .error}

We can use swap_axes=TRUE to have the phenotype on the x-axis. And we can use SEmult=2 to include mean ± 2 SE intervals.


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_pxg(g, iron$pheno[,"liver"], SEmult=2, swap_axes=TRUE, xlab="Liver phenotype")
~~~
{: .r}



~~~
Error in plot_pxg(g, iron$pheno[, "liver"], SEmult = 2, swap_axes = TRUE, : object 'g' not found
~~~
{: .error}

