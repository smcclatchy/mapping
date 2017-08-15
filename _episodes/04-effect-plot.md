---
title: "QTL Effects"
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




We may obtain plots indicating the estimated effects of the QTL via `plotPXG`, which creates a dot plot, or
`effectplot`, which plots the average phenotype for each genotype group.For `plotPXG`, we must first identify the marker closest to the QTL peak. Use
`find.marker`.


~~~
max(out.hk)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'out.hk' not found
~~~
{: .error}



~~~
mar <- find.marker(, chr=7, pos=47.7)
~~~
{: .r}



~~~
Error in find.marker(, chr = 7, pos = 47.7): argument "cross" is missing, with no default
~~~
{: .error}



~~~
plotPXG(hyper, marker=mar)
~~~
{: .r}



~~~
Error in b %in% colnames(a$data): object 'mar' not found
~~~
{: .error}

Note that red dots correspond to inferred genotypes (based on a single imputation).

The function `effectplot` uses the multiple imputation results from `sim.geno`.


~~~
effectplot(hyper, mname1=mar)
~~~
{: .r}



~~~
Warning in effectplot(hyper, mname1 = mar): -Running sim.geno.
~~~
{: .error}



~~~
Error in grep(pmalt.pattern, mname): object 'mar' not found
~~~
{: .error}

We may use `effectplot` at a position on the “grid” between markers, using "7@47.7" to indicate the position at 47.7 cM on chr 7.


~~~
effectplot(hyper, mname1="7@47.7")
~~~
{: .r}



~~~
Warning in effectplot(hyper, mname1 = "7@47.7"): -Running sim.geno.
~~~
{: .error}

<img src="../fig/rmd-04-unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" style="display: block; margin: auto;" />

Similar plots may be obtained for the locus on chr 15.


~~~
max(out.hk, chr=15)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'out.hk' not found
~~~
{: .error}



~~~
mar2 <- find.marker(hyper, chr=15, pos=12)
plotPXG(hyper, marker=mar2)
~~~
{: .r}

<img src="../fig/rmd-04-unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

~~~
effectplot(hyper, mname1="15@12")
~~~
{: .r}



~~~
Warning in effectplot(hyper, mname1 = "15@12"): -Running sim.geno.
~~~
{: .error}

<img src="../fig/rmd-04-unnamed-chunk-5-2.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" style="display: block; margin: auto;" />

We may plot the joint effects of the two loci via
`plotPXG` as follows:


~~~
plotPXG(hyper, marker=c(mar, mar2))
~~~
{: .r}



~~~
Error in b %in% colnames(a$data): object 'mar' not found
~~~
{: .error}



~~~
plotPXG(hyper, marker=c(mar2, mar))
~~~
{: .r}



~~~
Error in b %in% colnames(a$data): object 'mar' not found
~~~
{: .error}

The function `effectplot` gives more readable figures in this case; it’s often useful to look at it in both ways.


~~~
effectplot(hyper, mname1="7@47.7", mname2="15@12")
~~~
{: .r}



~~~
Warning in effectplot(hyper, mname1 = "7@47.7", mname2 = "15@12"): -Running
sim.geno.
~~~
{: .error}

<img src="../fig/rmd-04-unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />

~~~
effectplot(hyper, mname2="7@47.7", mname1="15@12")
~~~
{: .r}



~~~
Warning in effectplot(hyper, mname2 = "7@47.7", mname1 = "15@12"): -Running
sim.geno.
~~~
{: .error}

<img src="../fig/rmd-04-unnamed-chunk-7-2.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" style="display: block; margin: auto;" />

The two loci do not appear to interact.
