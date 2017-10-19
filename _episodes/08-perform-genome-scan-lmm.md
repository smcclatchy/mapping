---
title: "Performing a genome scan with a linear mixed model"
teaching: 20
exercises: 10
questions:
- "How do I use a linear mixed model in a genome scan?"
- "How do different mapping and kinship calculation methods differ?"
objectives:
- Create a genome scan with a linear mixed model.
- Compare LOD plots for Haley-Knott regression and linear mixed model methods.
- Compare LOD plots for the standard kinship matrix with the leave-one-chromosome-out (LOCO) method.
keypoints:
- "To perform a genome scan with a linear mixed model, supply a kinship matrix."
- "Different mapping and kinship calculation methods give different results."
source: Rmd
---



Why would a use a linear mixed model?
explanatory text and a graphic
use an example

When experimental design uses multiple groups or batches, a linear mixed model is in order. For example, 2 strains of cucumber on 4 different plots at each of 5 different farms. group effect or batch effect; 

To perform a genome scan using a linear mixed model, accounting for relationships among individuals using a random polygenic effect, you also use the function `scan1`; you just need to provide the argument `kinship`, a kinship matrix (or, for the LOCO method, a list of kinship matrices).


~~~
out_pg <- scan1(pr, iron$pheno, kinship, Xcovar=Xcovar)
~~~
{: .r}

Again, on a multi-core machine, you can get some speed-up using the `cores` argument.


~~~
out_pg <- scan1(pr, iron$pheno, kinship, Xcovar=Xcovar, cores=4)
~~~
{: .r}

For the LOCO (leave one chromosome out) method, provide the list of kinship matrices as obtained from `calc_kinship()` with `method="loco"`.


~~~
out_pg_loco <- scan1(pr, iron$pheno, kinship_loco, Xcovar=Xcovar)
~~~
{: .r}

To plot the results, we again use `plot_scan1()` from the [qtl2plot](https://github.com/rqtl/qtl2plot) package, or just type `plot()`.

Here is a plot of the LOD scores, by Haley-Knott regression and the linear mixed model using either the standard kinship matrix or the LOCO method.


~~~
color <- c("slateblue", "violetred", "green3")
par(mar=c(4.1, 4.1, 1.6, 1.1))
ymx <- max(maxlod(out), maxlod(out_pg), maxlod(out_pg_loco))
for(i in 1:2) {
    plot(out, map, lodcolumn=i, col=color[1], main=colnames(iron$pheno)[i],
              ylim=c(0, ymx*1.02))
    plot(out_pg, map, lodcolumn=i, col=color[2], add=TRUE)
    plot(out_pg_loco, map, lodcolumn=i, col=color[3], add=TRUE, lty=2)
    legend("topleft", lwd=2, col=color, c("H-K", "LMM", "LOCO"), bg="gray90", lty=c(1,1,2))
}
~~~
{: .r}

For the liver phenotype (top panel), the three methods give quite different results. The linear mixed model with an overall kinship matrix gives much lower LOD scores than the other two methods.  On chromosomes with some evidence of a QTL, the LOCO method gives higher LOD scores than Haley-Knott, except on chromosome 16 where it gives lower LOD scores.

For the spleen phenotype (bottom panel), the linear mixed model with an overall kinship matrix again gives much lower LOD scores than the other two methods. However, in this case Haley-Knott regression and
the LOCO method give quite similar results.
