---
title: "SNP association"
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



## SNP association

For multi-parent crosses, it can be useful to collapse the genotype or
allele probabilities according to the founder genotypes of the various
SNPs in the region of a QTL.

### QTL analysis in Diversity Outbred mice

To illustrate this sort of SNP association analysis, we'll consider some
Diversity Outbred mouse data. The Diversity Outcross (DO) mice are an
advanced intercross population derived from the same eight founder strains
as the Collaborative Cross (CC). See
[Svenson et al. (2012)](https://www.ncbi.nlm.nih.gov/pubmed/22345611)
and
[Gatti et al. (2014)](https://www.ncbi.nlm.nih.gov/pubmed/25237114).

We'll consider a subset of the data from
[Recla et al. (2014)](https://www.ncbi.nlm.nih.gov/pubmed/24700285),
available as part of the
[qtl2data github repository](https://github.com/rqtl/qtl2data). (The
full data are in
[`DO_Recla`](https://github.com/rqtl/qtl2data/tree/master/DO_Recla);
the directory
[`DOex`](https://github.com/rqtl/qtl2data/tree/master/DOex) contains a
reduced set, with just three chromosomes, one phenotype
(`OF_immobile_pct`, percent immobile in the open field test), and a
reduced set of markers.

You can download the data from a single zip file, as follows:


~~~
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex.zip")
DOex <- read_cross2(file)
~~~
{: .r}

Let's quickly whip through a basic analysis.

We first calculate genotype probabilities and convert them to allele
probabilities. We'll just use marker locations and not insert any
pseudomarkers.


~~~
pr <- calc_genoprob(DOex, error_prob=0.002)
apr <- genoprob_to_alleleprob(pr)
~~~
{: .r}



We calculate kinship matrices (using the "loco" method, though with the
caveat that here we are only considering genotypes on three chromosomes).


~~~
k <- calc_kinship(apr, "loco")
~~~
{: .r}

We create a numeric covariate for sex; be sure to include the individual
IDs as names.


~~~
sex <- (DOex$covar$Sex == "male")*1
names(sex) <- rownames(DOex$covar)
~~~
{: .r}

We perform a genome scan with a linear mixed model (adjusting for
a residual polygenic effect), with sex as an additive covariate.


~~~
out <- scan1(apr, DOex$pheno, k, sex)
~~~
{: .r}



~~~
Error in scan1(apr, DOex$pheno, k, sex): could not find function "scan1"
~~~
{: .error}

Here's a plot of the results.


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out, DOex$gmap)
~~~
{: .r}



~~~
Error in plot(out, DOex$gmap): object 'out' not found
~~~
{: .error}

There's a strong peak on chromosome 2. Let's look at the QTL effects.
We estimate them with `scan1coef()`. We need to subset the allele
probabilities and the list of kinship matrices.


~~~
coef_c2 <- scan1coef(apr[,"2"], DOex$pheno, k[["2"]], sex)
~~~
{: .r}



~~~
Error in scan1coef(apr[, "2"], DOex$pheno, k[["2"]], sex): could not find function "scan1coef"
~~~
{: .error}

For the DO, with 8 QTL alleles, we can use the function `plot_coefCC`
in the [R/qtl2plot](https://github.com/rqtl/qtl2plot) package, which
plots the 8 allele effects in the "official" Collaborative Cross (CC)
colors. (Well, actually _slightly_ modified colors, because I think the
official colors are kind of ugly.) The strong locus seems to be mostly
due to the NZO allele. Note that `CCcolors` is a vector of colors
included in the qtl2plot package; there's also a `CCorigcolors` object
with the _official_ colors.


~~~
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_coefCC(coef_c2, DOex$gmap["2"], bgcolor="gray95")
~~~
{: .r}



~~~
Error in plot_coefCC(coef_c2, DOex$gmap["2"], bgcolor = "gray95"): could not find function "plot_coefCC"
~~~
{: .error}



~~~
legend("bottomleft", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")
~~~
{: .r}



~~~
Error in legend("bottomleft", col = CCcolors, names(CCcolors), ncol = 2, : object 'CCcolors' not found
~~~
{: .error}






















