---
title: "Finding LOD peaks"
teaching: 30
exercises: 30
questions:
- "How do I locate LOD peaks above a certain threshold value?"
objectives:
- Locate LOD peaks above a threshold value throughout the genome.
- Identify the Bayesian credible interval for a QTL peak.
keypoints:
- "LOD peaks and support intervals can be identified with find_peaks()."
source: Rmd
---





Once we have LOD scores from a genome scan and a significance threshold, we can look for QTL peaks associated with the phenotype. High LOD scores indicate the neighborhood of a QTL but don't give its precise position. To find the exact position of a QTL, we define an interval that is likely to hold the QTL.

We'll use the Bayesian credible interval, which is a method for defining QTL intervals. It describes the probability that the interval contains the true value. Credible intervals make a probabilistic statement about the true value, for example, a 95% credible interval states that there is a 95% chance that the true value lies within the interval.

To find peaks above a given threshold LOD value, use the function `find_peaks()`. It can also provide Bayesian credible intervals by using the argument `prob` (the nominal coverage for the Bayes credible intervals). Set the argument `expand2markers = FALSE` to keep from expanding the interval out to typed markers, or exclude this argument if you'd like to include flanking markers.

You need to provide both the `scan1()` output, the marker/pseudomarker map and a threshold. We will use the 95% threshold from the permutations in the previous lesson.


~~~
operm <- scan1perm(pr, iron$pheno, Xcovar=Xcovar, n_perm=1000)
thr = summary(operm)
find_peaks(scan1_output = out, map = map, threshold = thr, prob = 0.95, expand2markers = FALSE)
~~~
{: .r}



~~~
  lodindex lodcolumn chr  pos       lod ci_lo ci_hi
1        1     liver   2 56.8  4.957564  54.3  70.3
2        1     liver   7 50.1  4.050766  17.1  53.6
3        1     liver   8 40.0  3.802511  34.0  55.0
4        1     liver  16 28.6  7.681569  21.6  32.6
5        2    spleen   8 13.6  4.302919   5.0  23.0
6        2    spleen   9 56.6 12.063226  54.6  58.6
~~~
{: .output}

The `find_peaks()` function can also pick out multiple peaks on a chromosome: each peak must exceed the chosen threshold, and the argument `peakdrop` indicates the amount that the LOD curve must drop between the lowest of two adjacent peaks.  Use this feature with caution.


~~~
find_peaks(scan1_output = out, map = map, threshold = thr, peakdrop = 1.8, prob = 0.95, expand2markers = FALSE)
~~~
{: .r}



~~~
  lodindex lodcolumn chr  pos       lod ci_lo ci_hi
1        1     liver   2 56.8  4.957564  54.3  70.3
2        1     liver   7 25.1  4.040021  15.1  27.1
3        1     liver   7 50.1  4.050766  41.1  53.6
4        1     liver   8 40.0  3.802511  34.0  55.0
5        1     liver  16 28.6  7.681569  21.6  32.6
6        2    spleen   8 13.6  4.302919   5.0  23.0
7        2    spleen   9 56.6 12.063226  54.6  58.6
~~~
{: .output}

Each row shows a different peak; the columns show the peak location, LOD score and the lower and upper interval endpoints.

> ## Challenge 1
> Find peaks in the genome scan object called `out` that meet a threshold of 3 and are in the interval described by a 2 point LOD drop from the peak. How many peaks meet the LOD threshold of 3 and lie within the interval defined by a 2 point LOD drop from the maximum peaks on each chromosome?
>
> > ## Solution to Challenge 1
> > `find_peaks(out, map, threshold=3, drop=2)` produces 7 peaks on 6 different chromosomes that meet a LOD threshold of 3 and are within the interval defined by a 2-LOD drop from the maximum peak on each chromosome.
> {: .solution}
{: .challenge}


> ## Challenge 2
> 1). Calculate the 90% Bayes credible interval on chromosome 16 for the liver phenotype (lodcolumn = 1).
What is the range of this interval that has a 90% chance of containing the true QTL position?  
2). Calculate the 95% Bayes credible interval for the same chromosome and phenotype. How does the interval change as you increase the probability? Why?
>
> > ## Solution to Challenge 2
> >
> > 1). `find_peaks(out, map, lodcolumn=1, chr=16, prob=0.90, expand2markers = FALSE)` produces a range from 25.1 to 40.4.  
> > 2). `find_peaks(out, map, lodcolumn=1, chr=16, prob=0.95, expand2markers = FALSE)` produces a range from 6.6 to 40.4, which is much broader than that of a 90% probability. The interval widens because the probability that the interval contains the true QTL position has increased. 
> {: .solution}
{: .challenge}
