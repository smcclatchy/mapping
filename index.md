---
layout: lesson
root: .
---

<head>
<!-- Global site tag (gtag.js) - Google Analytics -->
<script async src="https://www.googletagmanager.com/gtag/js?id=UA-146162546-1"></script>
<script>
window.dataLayer = window.dataLayer || [];
function gtag(){dataLayer.push(arguments);}
gtag('js', new Date());
gtag('config', 'UA-146162546-1');
</script>
</head>

This lesson introduces genetic mapping using qtl2, a R package for analyzing quantitative phenotypes and genetic data from complex crosses like the Diversity Outbred (DO). Genetic mapping with qtl2 allows researchers in fields as diverse as medicine, evolution, and agriculture to identify specific chromosomal regions that contribute to variation in phenotypes (quantitative trait loci or QTL). The goal is to identify the action, interaction, number, and precise location of these regions.

Participants will learn to
* calculate genotype and allele probabilities
* perform a genome scan and plot the results
* evaluate statistical significance of results
* find estimated effects of a QTL on a phenotype
* account for relationships among individuals by using a kinship matrix
* perform SNP association analysis

The lesson concludes with a complete analytical workflow from a study of DO mice.The lesson is adapted from [Karl Broman's](http://kbroman.org/) [software](http://kbroman.org/pages/software.html), tutorials, and book co-authored with [Saunak Sen](http://senresearch.org/), [A Guide to QTL Mapping with R/qtl](http://www.rqtl.org/book/).

> ## Prerequisites
> Understand fundamental genetic principles  
> Know how to access files not in the working directory by specifying the path    
> Know how to install a R package  
> Know how to assign a value to a variable  
> Know how to apply a built-in function  
{: .prereq}
