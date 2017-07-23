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

[R/qtl2](http://kbroman.org/qtl2) (aka qtl2) is a reimplementation of the QTL analysis software
[R/qtl](http://rqtl.org), to better handle high-dimensional data
and complex cross designs.

We expect that basic analyses with R/qtl2 will generally be performed
in "batch" (for example, on a cluster) rather than interactively. And
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

