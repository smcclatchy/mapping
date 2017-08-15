---
title: "Interval estimates of QTL location"
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



For the blood pressure phenotype, we’ve seen good evidence for QTL on chromosomes 1 and 4. Interval estimates of the location of QTL are commonly obtained via 1.5-LOD support intervals, which may be calculated via the function `lodint`. Alternatively, an approximate Bayes credible interval may be obtained with `bayesint`.



To obtain the 1.5-LOD support interval and 95% Bayes interval for the QTL on chromosome 4, type:


~~~
lodint(out.hk, chr=4)
~~~
{: .r}



~~~
         chr  pos      lod
c4.loc19   4 19.0 6.579917
D4Mit164   4 29.5 8.086160
D4Mit178   4 30.6 6.467713
~~~
{: .output}



~~~
bayesint(out.hk, chr=4)
~~~
{: .r}



~~~
         chr  pos      lod
c4.loc17   4 17.0 6.225904
D4Mit164   4 29.5 8.086160
c4.loc31   4 31.0 6.368718
~~~
{: .output}

The first and last rows define the ends of the intervals; the middle row is the estimated QTL location.

It is sometimes useful to identify the closest flanking markers; use `expandtomarkers=TRUE`:


~~~
lodint(out.hk, chr=4, expandtomarkers=TRUE)
~~~
{: .r}



~~~
         chr  pos      lod
D4Mit286   4 18.6 6.493607
D4Mit164   4 29.5 8.086160
D4Mit178   4 30.6 6.467713
~~~
{: .output}



~~~
bayesint(out.hk, chr=4, expandtomarkers=TRUE)
~~~
{: .r}



~~~
         chr  pos      lod
D4Mit108   4 16.4 5.508811
D4Mit164   4 29.5 8.086160
D4Mit80    4 31.7 5.203672
~~~
{: .output}

We can calculate the 2-LOD support interval and the 99% Bayes interval as follows.


~~~
lodint(out.hk, chr=4, drop=2)
~~~
{: .r}



~~~
         chr  pos      lod
D4Mit108   4 16.4 5.508811
D4Mit164   4 29.5 8.086160
D4Mit80    4 31.7 5.203672
~~~
{: .output}



~~~
bayesint(out.hk, chr=4, prob=0.99)
~~~
{: .r}



~~~
         chr  pos      lod
c4.loc13   4 13.0 5.606094
D4Mit164   4 29.5 8.086160
c4.loc31   4 31.0 6.368718
~~~
{: .output}

The intervals for the chr 1 locus may be calculated as follows.


~~~
lodint(out.hk, chr=1)
~~~
{: .r}



~~~
         chr  pos      lod
c1.loc32   1 35.3 1.846423
c1.loc45   1 48.3 3.548338
c1.loc82   1 85.3 2.034591
~~~
{: .output}



~~~
bayesint(out.hk, chr=1)
~~~
{: .r}



~~~
         chr  pos      lod
c1.loc33   1 36.3 2.447204
c1.loc45   1 48.3 3.548338
c1.loc80   1 83.3 2.702756
~~~
{: .output}
