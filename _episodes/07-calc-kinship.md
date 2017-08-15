---
title: "Calculating A Kinship Matrix"
teaching: 0
exercises: 0
questions:
- "How to I calculate kinship between individuals?"
objectives:
- Create a kinship matrix for individuals.
keypoints:
- "Kinship matrices account for relationships among individuals."
- "Kinship is calculated as the proportion of shared alleles between individuals."
- "Kinship calculation is a precursor to a genome scan via linear mixed model."
source: Rmd
---



If you wish to perform a genome scan by a linear mixed model, accounting for the relationships among individuals (in other words, including a random polygenic effect), you'll need to calculate a
kinship matrix for the individuals. This is accomplished with the `calc_kinship()` function in
[qtl2geno](https://github.com/rqtl/qtl2geno).
It takes the genotype probabilities as input.


~~~
kinship <- calc_kinship(pr)
~~~
{: .r}

By default, the genotype probabilities are converted to allele probabilities, and the kinship matrix is calculated as the proportion of shared alleles. To use genotype probabilities instead, use `use_allele_probs=FALSE`. Further, by default we omit the X chromosome and only use the autosomes. To include the X chromosome, use `omit_x=FALSE`.

In calculating the kinship matrix, one may wish to eliminate the effect of varying marker density across
the genome, and only use the probabilities along the grid of pseudomarkers (defined by the `step`
argument to `insert_pseudomarkers()`. To do so, we need to first use `calc_grid()` to determine the grid of pseudomarkers, and then `probs_to_grid()` to probabilities for positions that are not on the
grid.


~~~
grid <- calc_grid(iron$gmap, step=1)
pr_grid <- probs_to_grid(pr, grid)
kinship_grid <- calc_kinship(pr_grid)
~~~
{: .r}

If, for your linear mixed model genome scan, you wish to use the "leave one chromosome out" (LOCO) method (scan each chromosome using a kinship matrix that is calculated using data from all other chromosomes), use `type="loco"` in the call to `calc_kinship()`.


~~~
kinship_loco <- calc_kinship(pr, "loco")
~~~
{: .r}

On a multi-core machine, you can get some speed-up via the `cores` argument, as with `calc_genoprob()`.


~~~
kinship_loco <- calc_kinship(pr, "loco", cores=4)
~~~
{: .r}
