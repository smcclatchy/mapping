---
title: "Special covariates for the X chromosome"
teaching: 15
exercises: 15
questions:
- "How do I find the chromosome X covariates for a cross?"
objectives:
- Get the X covariates for a cross.
keypoints:
- "The X chromosome requires special treatment in QTL mapping."
- "Special covariates such as sex should be included to avoid spurious evidence of linkage."
source: Rmd
---



The X chromosome must be treated differently than the autosomes in the mapping of quantitative trait loci (QTL). If the X chromosome is treated like an autosome, a sex difference in a phenotype, such as weight or height, can lead to spurious linkage on the X chromosome. The X chromosome varies depending on the sex of the animal and the direction of the cross, so accounting for these covariates is important under the null hypothesis of no QTL, to avoid spurious evidence of linkage. (See [Broman et al. (2006) Genetics 174:2151-2158](http://www.genetics.org/content/174/4/2151.long).)

The particular X chromosome covariates depends on the cross, and can be obtained with the [qtl2geno](https://github.com/rqtl/qtl2geno) function `get_x_covar()`. In the iron data, sex and cross direction are indicated with 0s and 1s.


~~~
Xcovar <- get_x_covar(iron)
head(Xcovar)
~~~
{: .r}
