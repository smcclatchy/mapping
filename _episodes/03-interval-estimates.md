---
title: "Interval Estimates of QTL Location"
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



For the blood pressure phenotype, we’ve seen good evidence for QTL on chromosomes 7 and 15. Interval estimates of the location of QTL are commonly obtained via 1.5-LOD support intervals, which may be calculated via the function `lodint`. Alternatively, an approximate Bayes credible interval may be obtained with `bayesint`.

To obtain the 1.5-LOD support interval and 95% Bayes interval for the QTL on chromosome 7, type:


~~~
lodint(out.hk, chr=7)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}



~~~
bayesint(out.hk, chr=7)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}

The first and last rows define the ends of the intervals; the middle row is the estimated QTL location.

It is sometimes useful to identify the closest flanking markers; use
`expandtomarkers=TRUE`:


~~~
lodint(out.hk, chr=7, expandtomarkers=TRUE)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}



~~~
bayesint(out.hk, chr=7, expandtomarkers=TRUE)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}

We can calculate the 2-LOD support interval and the 99% Bayes interval as follows.


~~~
lodint(out.hk, chr=7, drop=2)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}



~~~
bayesint(out.hk, chr=7, prob=0.99)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}

The intervals for the chr 15 locus may be calculated as follows.


~~~
lodint(out.hk, chr=15)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}



~~~
bayesint(out.hk, chr=15)
~~~
{: .r}



~~~
Error in "scanone" %in% class(results): object 'out.hk' not found
~~~
{: .error}
