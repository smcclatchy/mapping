---
title: "Performing a genome scan"
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



## Performing a genome scan

To perform a genome scan by Haley-Knott regression
([Haley and Knott 1992](https://www.ncbi.nlm.nih.gov/pubmed/16718932)),
use the function `scan1()` in
[qtl2scan](https://github.com/rqtl/qtl2scan). (The functions above
were all from [qtl2geno](https://github.com/rqtl/qtl2geno); from here
forward, the functions are all from
[qtl2scan](https://github.com/rqtl/qtl2scan).)

`scan1()` takes as input the genotype probabilities, a matrix of
phenotypes, and then optional additive and interactive covariates, and
the special X chromosome covariates. Another option is to provide a
vector of weights.


~~~
library(qtl2scan)
out <- scan1(pr, iron$pheno, Xcovar=Xcovar)
~~~
{: .r}



~~~
Error in check4names(pheno, addcovar, Xcovar, intcovar): object 'iron' not found
~~~
{: .error}

On a multi-core machine, you can get some speed-up via the `cores`
argument, as with `calc_genoprob()` and `calc_kinship()`.


~~~
out <- scan1(pr, iron$pheno, Xcovar=Xcovar, cores=4)
~~~
{: .r}

The output of `scan1()` is a matrix of LOD scores, positions &times;
phenotypes. (Well, actually, the output is a list including
`"lod"` which is this matrix of LOD scores, but also including the
genetic map and other details.)

The function `plot_scan1()` in the
[qtl2plot](https://github.com/rqtl/qtl2plot) package can be used to
plot the LOD curves. You can just write `plot()`, as there's an S3
method `plot.scan1()` and the output of `scan1()` has class `"scan1"`.
Use the argument `lodcolumn` to indicate which column to plot. Unlike the
`plot.scanone()` function in [R/qtl](http://rqtl.org), the
`plot_scan1()` function will only plot one set of LOD curves at a
time, and you need to provide the marker/pseudomarker map (created by
`insert_pseudomarkers()`).


~~~
library(qtl2plot)
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out) # overall maximum LOD score
~~~
{: .r}



~~~
Error in "scan1coef" %in% class(scan1_output): object 'out' not found
~~~
{: .error}



~~~
plot(out, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
~~~
{: .r}



~~~
Error in plot(out, map, lodcolumn = 1, col = "slateblue", ylim = c(0, : object 'out' not found
~~~
{: .error}



~~~
plot(out, map, lodcolumn=2, col="violetred", add=TRUE)
~~~
{: .r}



~~~
Error in plot(out, map, lodcolumn = 2, col = "violetred", add = TRUE): object 'out' not found
~~~
{: .error}



~~~
legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out), bg="gray90")
~~~
{: .r}



~~~
Error in is.data.frame(x): object 'out' not found
~~~
{: .error}

