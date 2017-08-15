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



For the blood pressure phenotype, we’ve seen good evidence for QTL on chromosomes 7 and 15. Interval estimates of the location of QTL are commonly obtained via 1.5-LOD support intervals, which may be calculated via the function `lodint`. Alternatively, an approximate Bayes credible interval may be obtained with `bayesint`.



To obtain the 1.5-LOD support interval and 95% Bayes interval for the QTL on chromosome 7, type:


~~~
lodint(out.hk, chr=7)
~~~
{: .r}



~~~
         chr  pos        lod
D7Mit306   7  1.1 0.22879044
D7Mit297   7 26.2 0.47251936
D7Nds4     7 55.6 0.08926882
~~~
{: .output}



~~~
bayesint(out.hk, chr=7)
~~~
{: .r}



~~~
         chr  pos        lod
D7Mit306   7  1.1 0.22879044
D7Mit297   7 26.2 0.47251936
D7Nds4     7 55.6 0.08926882
~~~
{: .output}

The first and last rows define the ends of the intervals; the middle row is the estimated QTL location.

It is sometimes useful to identify the closest flanking markers; use `expandtomarkers=TRUE`:


~~~
lodint(out.hk, chr=7, expandtomarkers=TRUE)
~~~
{: .r}



~~~
         chr  pos        lod
D7Mit306   7  1.1 0.22879044
D7Mit297   7 26.2 0.47251936
D7Nds4     7 55.6 0.08926882
~~~
{: .output}



~~~
bayesint(out.hk, chr=7, expandtomarkers=TRUE)
~~~
{: .r}



~~~
         chr  pos        lod
D7Mit306   7  1.1 0.22879044
D7Mit297   7 26.2 0.47251936
D7Nds4     7 55.6 0.08926882
~~~
{: .output}

We can calculate the 2-LOD support interval and the 99% Bayes interval as follows.


~~~
lodint(out.hk, chr=7, drop=2)
~~~
{: .r}



~~~
         chr  pos        lod
D7Mit306   7  1.1 0.22879044
D7Mit297   7 26.2 0.47251936
D7Nds4     7 55.6 0.08926882
~~~
{: .output}



~~~
bayesint(out.hk, chr=7, prob=0.99)
~~~
{: .r}



~~~
         chr  pos        lod
D7Mit306   7  1.1 0.22879044
D7Mit297   7 26.2 0.47251936
D7Nds4     7 55.6 0.08926882
~~~
{: .output}

The intervals for the chr 15 locus may be calculated as follows.


~~~
lodint(out.hk, chr=15)
~~~
{: .r}



~~~
          chr  pos      lod
D15Mit11   15  5.5 1.114985
c15.loc13  15 18.5 1.750730
D15Mit79   15 63.4 1.722788
~~~
{: .output}



~~~
bayesint(out.hk, chr=15)
~~~
{: .r}



~~~
          chr  pos      lod
c15.loc6   15 11.5 1.178943
c15.loc13  15 18.5 1.750730
D15Mit79   15 63.4 1.722788
~~~
{: .output}
