---
title: "Introduction"
teaching: 15
exercises: 0
questions:
- "What is quantitative trait mapping?"
- "Where can I find sample data for mapping with the qtl and qtl2 packages?"
objectives:
- "Be able to find sample data for qtl mapping."
keypoints:
- "Published and public data already formatted for QTL mapping are available on the web."
- "These data can be used as a model for formatting your own QTL data."
---

Quantitative trait mapping is used in biomedical, agricultural, and evolutionary studies
to find causal genes for quantitative traits, to aid crop and breed selection in agriculture,
and to shed light on natural selection. The goal of quantitative trait locus (QTL) analysis
is to identify genomic regions linked to a phenotype, to map these regions precisely,
and to define the effects, number, and interactions of QTL. QTL analysis can 
be performed in natural populations and in experimental crosses, and can be studied in 
humans and non-human species. Human studies, however, are very expensive, lack environmental
control, and can be confounded by population structure such that associations between
genotype and phenotype  are not necessarily causal.

QTL analysis in experimental crosses requires two or more strains that differ genetically
with regard to a phenotype of interest. Genetic markers, such as SNPs or microsatellites,
distinguish between parental strains in the experimental cross. Markers that are genetically
linked to a phenotype will segregate more often with phenotype values (high or low values, for example),
while unlinked markers will not be significantly associated with the phenotype. The markers 
themselves might be associated with the phenotype but are not causal. Rather, markers may 
be associated with the phenotype through linkage to nearby QTL. They serve as signposts
indicating the neighborhood of a QTL that influences a phenotype. Covariates such as sex
or diet can also influence the phenotype.

R/qtl is a package for mapping quantitative trait loci (QTL) in experimental crosses.
R/qtl methods employ hidden Markov models specifically for handling missing genotype
data. Hidden Markov model algorithms are implemented in R/qtl for 
backcrosses, intercrosses, and phase-known four-way crosses. You can use
R/qtl to estimate genetic maps, identify genotyping errors, and perform 
single-QTL genome scans and two-QTL, two-dimensional genome scans, 
by interval mapping (with the EM algorithm), Haley-Knott regression, 
and multiple imputation. You can also use covariates such as sex, age or treatment
in genome scans.

[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) is a reimplementation of the QTL analysis software
[R/qtl](http://rqtl.org), to better handle high-dimensional data
and complex cross designs such as the Diversity Outbred. Typically R/qtl2 will 
be employed in "batch" (for example, on a cluster) rather than interactively. And
so the software is split into three parts:
[qtl2geno](https://github.com/rqtl/qtl2geno) for genotype probability
calculations, [qtl2scan](https://github.com/rqtl/qtl2scan) for QTL
scans, and [qtl2plot](https://github.com/rqtl/qtl2plot) for data
visualization.

A further package, [qtl2convert](https://github.com/rqtl/qtl2convert),
contains functions for converting data among the R/qtl2,
[DOQTL](https://www.bioconductor.org/packages/release/bioc/html/DOQTL.html),
and [R/qtl](http://rqtl.org) formats, for example to convert genotype
probabilities produced by DOQTL to the format needed by qtl2scan, or
to convert qtl2scan results to the format produced by `scanone` in
R/qtl, so that they may be graphed with the R/qtl functions.

This lesson will focus on the R/qtl2 package in R. To learn how to use
R/qtl, see Karl Broman's [tutorials](http://rqtl.org/tutorials).  

## Sample data sets

The R/qtl2 web site includes
[sample data files](http://kbroman.org/qtl2/pages/sampledata.html) in
the new format. Zipped versions of these datasets are included with
the [qtl2geno](https://github.com/rqtl/qtl2geno) package and can be
loaded into R using the `read_cross2()` function.

In the [qtl2geno package source](https://github.com/rqtl/qtl2geno),
the sample zip files are located in
[`qtl2geno/inst/extdata`](https://github.com/rqtl/qtl2geno/tree/master/inst/extdata).
In the installed version of the package, they are in
`qtl2geno/extdata`, within whatever directory your R packages were
installed. The R function `system.file()` can be used to construct the
path to these files.

For example, one of the sample data sets concerns a gravitropism
phenotype in a set of Arabidopsis recombinant inbred lines (RIL), from
[Moore et al. (2013) Genetics 195:1077-1086](http://www.genetics.org/content/195/3/1077.abstract).
The data are in `qtl2geno/extdata/grav2.zip`, which can be loaded as
follows:

~~~
library(qtl2geno)
grav2 <- read_cross2( system.file("extdata", "grav2.zip", package="qtl2geno") )
~~~
{: r}

Additional sample data sets, including data on Diversity Outbred (DO)
mice, are available at <https://github.com/rqtl/qtl2data>.


To cite R/qtl in publications, use
Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping
in experimental crosses. Bioinformatics 19:889-89


