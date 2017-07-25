---
title: "Special covariates for the X chromosome"
teaching: 0
exercises: 0
questions:
- "How do I find the chromosome X covariates for a cross?"
objectives:
- Get the X covariates for a cross.
keypoints:
- "The X chromosome requires special treatment in QTL mapping."
- "Special covariates such as sex should be included to avoid spurious evidence of linkage."
source: Rmd
---



In a QTL scan of the X chromosome, special covariates (such as sex)
may need to be included under the null hypothesis of no QTL, to avoid
spurious evidence of linkage. (See
[Broman et al. (2006) Genetics 174:2151-2158](http://www.genetics.org/content/174/4/2151.long).)

The particular X chromosome covariates depends on the cross, and can
be obtained with the [qtl2geno](https://github.com/rqtl/qtl2geno)
function `get_x_covar()`.


~~~
Xcovar <- get_x_covar(iron)
~~~
{: .r}



~~~
Error in "cross2" %in% class(x): object 'iron' not found
~~~
{: .error}
