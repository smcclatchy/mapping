---
title: "Performing a genome scan with a linear mixed model"
teaching: 0
exercises: 0
questions:
- "?"
objectives:
- 
- 
keypoints:
- "."
- "."
source: Rmd
---



## Performing a genome scan with a linear mixed model

To perform a genome scan using a linear mixed model, accounting for
relationships among individuals using a random polygenic effect, you also use
the function `scan1`; you just need to provide the argument `kinship`,
a kinship matrix (or, for the LOCO method, a list of kinship
matrices).


~~~
out_pg <- scan1(pr, iron$pheno, kinship, Xcovar=Xcovar)
~~~
{: .r}



~~~
Error in scan1(pr, iron$pheno, kinship, Xcovar = Xcovar): object 'kinship' not found
~~~
{: .error}

Again, on a multi-core machine, you can get some speed-up using the
`cores` argument.


~~~
out_pg <- scan1(pr, iron$pheno, kinship, Xcovar=Xcovar, cores=4)
~~~
{: .r}

For the LOCO (leave one chromosome out) method, provide the list of
kinship matrices as obtained from `calc_kinship()` with
`method="loco"`.


~~~
out_pg_loco <- scan1(pr, iron$pheno, kinship_loco, Xcovar=Xcovar)
~~~
{: .r}



~~~
Error in scan1(pr, iron$pheno, kinship_loco, Xcovar = Xcovar): object 'kinship_loco' not found
~~~
{: .error}

To plot the results, we again use `plot_scan1()` from the
[qtl2plot](https://github.com/rqtl/qtl2plot) package, or just type `plot()`.

Here is a plot of the LOD scores, by Haley-Knott regression and the linear
mixed model using either the standard kinship matrix or the LOCO
method.


~~~
color <- c("slateblue", "violetred", "green3")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(maxlod(out), maxlod(out_pg), maxlod(out_pg_loco))
~~~
{: .r}



~~~
Error in "scan1coef" %in% class(scan1_output): object 'out' not found
~~~
{: .error}



~~~
for(i in 1:2) {
    plot(out, map, lodcolumn=i, col=color[1], main=colnames(iron$pheno)[i],
              ylim=c(0, ymx*1.02))
    plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
    plot(out_pg_loco, map, lodcolumn=i, col=color[3], add=TRUE, lty=2)
    legend("topleft", lwd=2, col=color, c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
}
~~~
{: .r}



~~~
Error in plot(out, map, lodcolumn = i, col = color[1], main = colnames(iron$pheno)[i], : object 'out' not found
~~~
{: .error}

For the liver phenotype (top panel), the three methods give quite
different results. The linear mixed model with an overall kinship
matrix gives much lower LOD scores than the other two methods.  On
chromosomes with some evidence of a QTL, the LOCO method gives higher
LOD scores than Haley-Knott, except on chromosome 16 where it gives
lower LOD scores.

For the spleen phenotype (bottom panel), the linear mixed model with an
overall kinship matrix again gives much lower LOD scores than the
other two methods. However, in this case Haley-Knott regression and
the LOCO method give quite similar results.
