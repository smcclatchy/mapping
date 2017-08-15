---
title: "Introduction"
teaching: 0
exercises: 0
questions:
- "Key question"
objectives:
- "First objective."
keypoints:
- "First key point."
---

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


