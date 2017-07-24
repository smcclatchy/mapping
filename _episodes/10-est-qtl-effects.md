---
title: "Estimated QTL effects"
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




## Estimated QTL effects

The `scan1()` function return only LOD scores. To
obtain estimated QTL effects, use the function `scan1coef()`.
This function takes a single phenotype and the
genotype probabilities for a single chromosome and returns a matrix
with the estimated coefficients at each putative QTL location along
the chromosome.

For example, to get the estimated effects on chromosome 2 for the
liver phenotype, we'd do the following:


~~~
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])
~~~
{: .r}



~~~
Error in check4names(pheno, addcovar, NULL, intcovar, nullcovar): object 'iron' not found
~~~
{: .error}






















