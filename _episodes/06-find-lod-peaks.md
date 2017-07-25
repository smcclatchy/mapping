---
title: "Finding LOD peaks"
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



## Finding LOD peaks

The function `find_peaks()` in the
[qtl2scan](https://github.com/rqtl/qtl2scan) package can be used to
identify a set of LOD peaks that exceed some threshold. It can also
provide LOD support or Bayes credible intervals, by using the
arguments `drop` (the amount to drop in the LOD support intervals)
or `prob` (the nominal coverage for the Bayes credible intervals).

You need to provide both the `scan1()` output as well as the
marker/pseudomarker map.


~~~
find_peaks(out, map, threshold=4, drop=1.5)
~~~
{: .r}



~~~
Error in find_peaks(out, map, threshold = 4, drop = 1.5): could not find function "find_peaks"
~~~
{: .error}

The `find_peaks()` function can also pick out multiple peaks on a
chromosome: each peak must exceed the chosen threshold, and the
argument `peakdrop` indicates the amount that the LOD curve must drop
between the lowest of two adjacent peaks.  Use this feature with
caution.


~~~
find_peaks(out, map, threshold=4, peakdrop=1.8, drop=1.5)
~~~
{: .r}



~~~
Error in find_peaks(out, map, threshold = 4, peakdrop = 1.8, drop = 1.5): could not find function "find_peaks"
~~~
{: .error}

The functions `lod_int()` and `bayes_int()` can be used to derive the
LOD support or Bayes credible intervals for QTL, for a specific
chromosome and LOD score column. For example, to obtain the Bayes
interval for the locus on chromosome 9 for the second phenotype
("spleen"):


~~~
bayes_int(out, map, lodcolumn=2, chr=9, prob=0.95)
~~~
{: .r}



~~~
Error in bayes_int(out, map, lodcolumn = 2, chr = 9, prob = 0.95): could not find function "bayes_int"
~~~
{: .error}

Both `lod_int()` and `bayes_int()` take a `peakdrop` argument, if you
wish to try to identify multiple peaks on a chromosome. Again, use
this feature with caution.


~~~
lod_int(out, map, lodcolumn=1, chr=7, peakdrop=1.8, drop=1.5)
~~~
{: .r}



~~~
Error in lod_int(out, map, lodcolumn = 1, chr = 7, peakdrop = 1.8, drop = 1.5): could not find function "lod_int"
~~~
{: .error}

Each row is a different peak; the columns are the lower interval endpoint, the
estimated QTL position, and the upper interval endpoint.
