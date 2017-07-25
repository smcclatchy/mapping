---
title: "Performing a genome scan with binary traits"
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



## Performing a genome scan with binary traits

The genome scans above were performed assuming that the residual
variation followed a normal distribution. This will often provide
reasonable results even if the residuals are not normal, but an
important special case is that of a binary trait, with values 0 and 1,
which is best treated differently. The `scan1` function can perform a
genome scan with binary traits by logistic regression, using the
argument `model="binary"`. (The default value for the `model` argument
is `"normal"`.) At present, we _can not_ account for relationships
among individuals in this analysis.

Let's first turn our two phenotypes into binary traits by thresholding
at the median. One would generally _not_ do this in practice;
this is just for illustration.


~~~
bin_pheno <- apply(iron$pheno, 2, function(a) as.numeric(a > median(a)))
~~~
{: .r}



~~~
Error in apply(iron$pheno, 2, function(a) as.numeric(a > median(a))): object 'iron' not found
~~~
{: .error}



~~~
rownames(bin_pheno) <- rownames(iron$pheno)
~~~
{: .r}



~~~
Error in rownames(iron$pheno): object 'iron' not found
~~~
{: .error}

We now perform the genome scan as before, including `model="binary"`
to indicates that the phenotypes are binary traits with values 0 and
1.


~~~
out_bin <- scan1(pr, bin_pheno, Xcovar=Xcovar, model="binary")
~~~
{: .r}



~~~
Error in scan1(pr, bin_pheno, Xcovar = Xcovar, model = "binary"): could not find function "scan1"
~~~
{: .error}

Here is a plot of the two LOD curves.


~~~
par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out_bin)
~~~
{: .r}



~~~
Error in maxlod(out_bin): could not find function "maxlod"
~~~
{: .error}



~~~
plot(out_bin, map, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
~~~
{: .r}



~~~
Error in plot(out_bin, map, lodcolumn = 1, col = "slateblue", ylim = c(0, : object 'out_bin' not found
~~~
{: .error}



~~~
plot(out_bin, map, lodcolumn=2, col="violetred", add=TRUE)
~~~
{: .r}



~~~
Error in plot(out_bin, map, lodcolumn = 2, col = "violetred", add = TRUE): object 'out_bin' not found
~~~
{: .error}



~~~
legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(out_bin), bg="gray90")
~~~
{: .r}



~~~
Error in is.data.frame(x): object 'out_bin' not found
~~~
{: .error}

We can use `find_peaks` as before.


~~~
find_peaks(out_bin, map, threshold=3.5, drop=1.5)
~~~
{: .r}



~~~
Error in find_peaks(out_bin, map, threshold = 3.5, drop = 1.5): could not find function "find_peaks"
~~~
{: .error}
