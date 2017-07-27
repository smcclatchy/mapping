---
title: "Performing a permutation test"
teaching: 0
exercises: 0
questions:
- "How can I evaluate the statistical significance of genome scan results?"
objectives:
- Run a permutation test to establish LOD score thresholds.
keypoints:
- "A permutation test establishes the statistical  significance of a genome scan."
source: Rmd
---




## Performing a permutation test

To perform a permutation test to establish the statistical significance of
the results of a genome scan, use the function `scan1perm()`. (In
[R/qtl](http://rqtl.org), a single function, `scanone()`, was used for
both performing a genome scan and for getting permutation-based
significance thresholds, but in [R/qtl2](http://kbroman.org/qtl2),
we've decided to make two separate functions).

The `scan1perm()` function takes the same arguments as `scan1()`, plus
additional arguments to control the permutations:

- `n_perm` is the number of permutation replicates.
- `perm_Xsp` controls whether to perform autosome/X chromosome
  specific permutations (with `perm_Xsp=TRUE`) or not (the default is
  to not).
- `perm_strata` is a vector that defines the strata for a stratified permutation
  test.
- `chr_lengths` is a vector of chromosome lengths, used in the case
  that `perm_Xsp=TRUE`.

As with `scan1()`, you may provide a kinship matrix (or vector of
kinship matrices, for the "leave one chromosome out" (loco) approach),
in order to fit linear mixed models to account for accounting for the
relationships among individuals (in other words, including a random
polygenic effect). If `kinship` is unspecified, the function performs
ordinary Haley-Knott regression.



To perform a permutation test with the `iron` data, we do the following:


~~~
operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000)
~~~
{: .r}

Note the need to specify special covariates for the X chromosome (via
`Xcovar`), to be included under the null hypothesis of no QTL.
And note that when these are provided, the default is to perform a
stratified permutation test, using strata defined by the rows in
`Xcovar`. In general, when the X chromosome is considered, one will
wish to stratify at least by sex.

Also note that, as with `scan1()`, you can speed up the calculations on a multi-core
machine by specifying the argument `cores`. With `cores=0`, the number
of available cores will be detected via `parallel::detectCores()`.
Otherwise, specify the number of cores as a positive integer.
For large datasets, be mindful of the amount of memory that will be
needed; you may need to use fewer than the maximum number of cores,
to avoid going beyond the available memory.


~~~
operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000, cores=0)
~~~
{: .r}

To get estimated significance thresholds, use the function
`summary()`.




~~~
summary(operm)
~~~
{: .r}


~~~
        liver   spleen
0.05 3.461161 3.463352
attr(,"class")
[1] "summary.scan1perm" "matrix"           
attr(,"n_perm")
     liver spleen
[1,]  1000   1000
~~~
{: .output}

The default is to return the 5% significance thresholds. Thresholds for
other (or for multiple) significance levels can be obtained via the
`alpha` argument.




~~~
summary(operm, alpha=c(0.2, 0.05))
~~~
{: .r}


~~~
        liver   spleen
0.2  2.625936 2.636614
0.05 3.461161 3.463352
attr(,"class")
[1] "summary.scan1perm" "matrix"           
attr(,"n_perm")
     liver spleen
[1,]  1000   1000
~~~
{: .output}

To obtain autosome/X chromosome-specific significance thresholds,
specify `perm_Xsp=TRUE`. In this case, you need to provide chromosome
lengths, which may be obtained with the function `chr_lengths()`.




~~~
operm2 <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000,
                    perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
~~~
{: .r}

Separate permutations are performed for the autosomes and X
chromosome, and considerably more permutation replicates are needed
for the X chromosome. The computations take about twice as much time.
See [Broman et al. (2006) Genetics
174:2151-2158](https://www.ncbi.nlm.nih.gov/pubmed/17028340).

The significance thresholds are again derived via `summary()`:




~~~
summary(operm2, alpha=c(0.2, 0.05))
~~~
{: .r}


~~~
$A
        liver   spleen
0.2  2.654181 2.543013
0.05 3.424862 3.224816

$X
        liver   spleen
0.2  3.097920 4.017580
0.05 3.896456 5.179283

attr(,"class")
[1] "summary.scan1perm" "list"             
attr(,"n_perm")
  liver spleen
A  1000   1000
X 28243  28243
~~~
{: .output}

Permutations for a genome scan with a linear mixed model-based are
performed by specifying the `kinship` argument. We can
use the "leave one chromosome out" (loco) method by providing
`kinship_loco`, the list of kinship matrices calculated above with
`calc_kinship()`.




~~~
operm3 <- scan1perm(pr, iron$pheno, kinship_loco, Xcovar=Xcovar, n_perm=1000,
                    perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
~~~
{: .r}

Here are the estimated significance thresholds:




~~~
summary(operm3, alpha=c(0.2, 0.05))
~~~
{: .r}


~~~
$A
        liver   spleen
0.2  2.641581 2.623347
0.05 3.286362 3.285426

$X
        liver   spleen
0.2  3.138534 4.365114
0.05 3.816997 5.496583

attr(,"class")
[1] "summary.scan1perm" "list"             
attr(,"n_perm")
  liver spleen
A  1000   1000
X 28243  28243
~~~
{: .output}

As with `scan1`, we can use `scan1perm` with binary traits, using the
argument `model="binary"`. Again, this can't be used with a kinship
matrix, but all of the other arguments can be applied.


~~~
operm_bin <- scan1perm(pr, bin_pheno, Xcovar=Xcovar, model="binary",
                       n_perm=1000, perm_Xsp=TRUE, chr_lengths=chr_lengths(map))
~~~
{: .r}

Here are the estimated 5% and 20% significance thresholds.




~~~
summary(operm_bin, alpha=c(0.2, 0.05))
~~~
{: .r}


~~~
$A
        liver   spleen
0.2  2.598312 2.629139
0.05 3.329495 3.409336

$X
        liver   spleen
0.2  3.161158 3.061490
0.05 3.858628 3.766575

attr(,"class")
[1] "summary.scan1perm" "list"             
attr(,"n_perm")
  liver spleen
A  1000   1000
X 28243  28243
~~~
{: .output}
