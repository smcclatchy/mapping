---
title: "QTL analysis in Diversity Outbred Mice"
teaching: 30
exercises: 30
questions:
- "How do I bring together each step in the workflow?"
- "How is the workflow implemented in an actual study?"
objectives:
- Work through an actual QTL mapping workflow from a published study.
keypoints:
- "."
- "."
source: Rmd
---





This tutorial will take you through the process of mapping a QTL and searching for candidate genes.

The data comes from a toxicology study in which Diversity Outbred (DO) mice were exposed to benzene via inhalation for 6 hours a day, 5 days a week for 4 weeks  [(French, J. E., et al. (2015) <i>Environ Health Perspect</i> 123(3): 237-245)](http://ehp.niehs.nih.gov/1408202/). The study was conducted in two equally sized cohort of 300 male mice each, for a total of 600 mice. They were then sacrificed and reticulocytes (red blood cell precursors) were isolated from bone marrow. The number of micro-nucleated reticulocytes, a measure of DNA damage, was then measured in each mouse. The goal is to map gene(s) that influence the level of DNA damage in the bone marrow.

![](../fig/benzene_study_design.png)


### Load and explore the data

Make sure that you are in your scripts directory. If you're not sure where you are working right now, you can check your working directory with `getwd()`. If you are not in your scripts directory, use `setwd()` in the Console or Session -> Set Working Directory -> Choose Directory in the RStudio menu to set your working directory to the scripts directory.

Once you are in your scripts directory, create a new R script.

The data for this tutorial has been saved as an R binary file that contains several data objects.  Load it in now by running the following command.


~~~
load("../data/qtl2_demo.Rdata")
~~~
{: .r}



~~~
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file '../
data/qtl2_demo.Rdata', probable reason 'No such file or directory'
~~~
{: .error}



~~~
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
~~~
{: .error}



~~~
sessionInfo()
~~~
{: .r}

The call to `sessionInfo` provides information about the R version running on your machine and the R packages that are installed. This information can help you to troubleshoot. If you can't load the data, try downloading the data again by typing the following into your browser: ftp://ftp.jax.org/dgatti/qtl2_workshop/qtl2_demo.Rdata

We loaded in two data objects. Look in the Environment pane to see what was loaded.  You should see an object called `pheno` with 143 observations (in rows) of 5 variables (in columns), an object called `map` and an object called `probs`.

#### Phenotypes  
`pheno` is a data frame containing the phenotype data. Click on the triangle to the left of `pheno` in the Environment pane to view its contents.

**NOTE:** the sample IDs must be in the rownames of `pheno`. For more information about data file format, see the [Setup](../setup.md) instructions.

`pheno` contains the sample ID, sex, the study cohort, the concentration of benzene and the proportion of bone marrow reticulocytes that were micro-nucleated `(prop.bm.MN.RET)`.  Note that the sample IDs are also stored in the rownames of `pheno`. In order to save time for this tutorial, we will only map with 143 samples from the 100 ppm dosing group.

> ## Challenge 1 Look at the Data
> Determine the dimensions of `pheno`.  
> 1). How many rows and columns does it have?  
> 2). What are the names of the variables it contains?  
> 3). What does the distribution of the `prop.bm.MN.RET` column look like? Is it normally distributed?  
>
> > ## Solution to Challenge 1
> > `dim(pheno)`  
> > 1). 143 rows and 5 columns  
> > 2). Use `colnames(pheno)` or click the triangle to the left of `pheno`
> > in the Environment tab.  
> > 3). Use `hist(pheno$prop.bm.MN.RET)` to plot the distribution of the data. The data has a long right 
> > tail and is not normally distributed.  
> {: .solution}
{: .challenge}

Many statistical models, including the QTL mapping model in qtl2, expect that the incoming data will be normally distributed. You may use transformations such as log or square root to make your data more normally distributed. We will be mapping the proportion of micronucleated reticulocytess in bone marrow (`prop.bm.MN.RET`) and, since the data does not contain zeros or negative numbers, we will log transform the data.


~~~
pheno$log.MN.RET = log(pheno$prop.bm.MN.RET)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'pheno' not found
~~~
{: .error}

Now, let's make a histogram of the log-transformed data.


~~~
hist(pheno$log.MN.RET)
~~~
{: .r}



~~~
Error in hist(pheno$log.MN.RET): object 'pheno' not found
~~~
{: .error}

This looks much better!

Some researchers are concerned about the reproducibility of DO studies. The argument is that each DO mouse is unique and therefor can never be reproduced. But this misses the point of using the DO. While mice are the sampling unit, in QTL mapping we are sampling the founder alleles at each locus. And an average of 1/8th of the alleles should come from each founder at any given locus. Also, DO mice are a *population* of mice, not a single strain. While it is true that results in an individual DO mouse may not be reproducible, results at the population level should be reproducible. This is similar to the human population in that results from one individual may not represent all humans, but results at the population level should be reproducible.

The benzene inhalation study was conducted in two separate cohorts (termed "Study" in the `pheno` file). Below, we plot the proportion of micronucleated reticulocytess in bone marrow versus the benzene concentration for each study cohort. As you can see, while individual mice have varying micronucleated reticulocytess, there is a dose-dependent increase in micronucleated reticulocytess in both cohorts. This is an example of how results in the DO reproduce at the population level.

![](../fig/bm_mnret_by_dose_cohort.png)


#### The Marker Map  

The marker map for each chromosome is stored in the `map` object. Each list element is a numeric vector with each marker position in megabases (Mb). This is used to plot the LOD scores calculated at each marker during QTL mapping.

The markers are from a genotyping array called the Mouse Universal Genotyping Array (MUGA) and contains 7,856 SNP markers. Marker locations for the MUGA and other mouse arrays are available from [The Jackson Laboratory's FTP site](ftp://ftp.jax.org/MUGA).

Look at the structure of `map` in the Environment tab in RStudio. 

> ## Challenge 2 Data Dimensions II
> 1). Determine the length of `map`.  
> 2). How mamy markers are chromosome 1?  
>
> > ## Solution to Challenge 2  
> > 1). `length(map)`  
> > 2). `length(map[[1]])`  
> {: .solution}
{: .challenge} 

#### Genotype (allele) probabilities  
Each element of `probs` is a 3 dimensional array containing the founder allele dosages for each sample at each marker on one chromosome. These probabilities have been pre-calculated for you, so you can skip the step for [calculating genotype probabilities](https://smcclatchy.github.io/mapping/03-calc-genoprob/) and the optional step for calculating allele probabilities.

Next, we look at the dimensions of `probs` for chromosome 1:


~~~
dim(probs[[1]])
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'probs' not found
~~~
{: .error}

`probs` is a three dimensional array containing the proportion of each founder haplotype at each marker for each DO sample.  The 143 samples are in the first dimension, the 8 founders in the second and the markers along chromosome 1 are in the third dimension.
Let's return to the `probs` object. Look at the contents for the first 500 markers of one sample.

**NOTE:** the sample IDs must be in the rownames of `probs`.


~~~
image(1:500, 1:ncol(probs[[1]]), t(probs[[1]][1,8:1,1:500]), breaks = 0:100/100,
      col = grey(99:0/100), axes = F, xlab = "Markers", ylab = "Founders",
      main = "Founder Allele Contributions for Sample 1")
~~~
{: .r}



~~~
Error in ncol(probs[[1]]): object 'probs' not found
~~~
{: .error}



~~~
abline(h = 0:8 + 0.5, col = "grey70")
~~~
{: .r}



~~~
Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
~~~
{: .error}



~~~
usr = par("usr")
rect(usr[1], usr[3], usr[2], usr[4])
~~~
{: .r}



~~~
Error in rect(usr[1], usr[3], usr[2], usr[4]): plot.new has not been called yet
~~~
{: .error}



~~~
axis(side = 1, at = 0:5 * 100, labels = 0:5 * 100)
~~~
{: .r}



~~~
Error in axis(side = 1, at = 0:5 * 100, labels = 0:5 * 100): plot.new has not been called yet
~~~
{: .error}



~~~
axis(side = 2, at = 1:8, labels = LETTERS[8:1], las = 1, tick = F)
~~~
{: .r}



~~~
Error in axis(side = 2, at = 1:8, labels = LETTERS[8:1], las = 1, tick = F): plot.new has not been called yet
~~~
{: .error}

In the plot above, the founder contributions, which range between 0 and 1, are colored from white (= 0) to black (= 1.0). A value of ~0.5 is grey. The markers are on the X-axis and the eight founders (denoted by the letters A through H) on the Y-axis. Starting at the left, we see that this sample has genotype CD because both rows C and D are grey, indicating values of 0.5 for each one. Moving along the genome to the right, the genotype becomes DD where row D is black, then CD, AC, CH, CD, CH, etc. The values at each marker sum to 1.0.  


### [Calculating A Kinship Matrix](https://smcclatchy.github.io/mapping/04-calc-kinship/)
Next, we need to create a matrix that accounts for the kinship relationships between the mice. We do this by looking at the correlation between the founder haplotypes for each sample at each SNP. For each chromosome, we create a kinship matrix using all markers *except* the ones on the current chromosome using the loco (leave-one-chromosome-out) method. Simulations suggest that mapping using this approach increases the power to detect QTL.
           

~~~
K = calc_kinship(probs = probs, type = "loco", use_allele_probs = TRUE)
~~~
{: .r}



~~~
Error in calc_kinship(probs = probs, type = "loco", use_allele_probs = TRUE): object 'probs' not found
~~~
{: .error}

Kinship values between pairs of samples range between 0 (no relationship) and 1.0 (completely identical). Let's look at the kinship matrix.


~~~
image(1:nrow(K[[1]]), 1:ncol(K[[1]]), K[[1]][,ncol(K[[1]]):1], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship between samples", 
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
~~~
{: .r}



~~~
Error in nrow(K[[1]]): object 'K' not found
~~~
{: .error}

The figure above shows kinship between all pairs of samples. White ( = 1) indicates no kinship and red ( = 0) indicates full kinship. Orange values indicate varying levels of kinship between 0 and 1. The white diagonal of the matrix indicates that each sample is identical to itself. The lighter yellow blocks off of the diagonal may indicate siblings or cousins.

#### Covariates    
Next, we need to create additive covariates that will be used in the mapping model. We will use study cohort as a covariate in the mapping model. If we were mapping with *all* mice, we would also add benzene concentration to the model. This study contained only male mice, but in most cases, you would include sex as an additive covariate as well.


~~~
addcovar = model.matrix(~Study, data = pheno)[,-1]
~~~
{: .r}



~~~
Error in terms.formula(object, data = data): object 'pheno' not found
~~~
{: .error}

The code above copies the `rownames(pheno)` to `rownames(addcovar)` as a side-effect.

**NOTE:** the sample IDs must be in the rownames of `pheno`, `addcovar`, `genoprobs` and `K`. `qtl2` uses the sample IDs to align the samples between objects. For more information about data file format, see [Karl Broman's vignette on input file format](http://kbroman.org/qtl2/assets/vignettes/input_files.html).

### [Performing a genome scan](https://smcclatchy.github.io/mapping/06-perform-genome-scan/)  

Before we run the mapping function, let's look at the mapping model. At each marker on the genotyping array, we will fit a model that regresses the phenotype (micronucleated reticulocytess) on covariates and the founder allele proportions.  

![](../fig/equation1.png)

  where:  
  
<ul>
<li><i>y<sub>i</sub></i> is the phenotype for mouse <i>i</i>,</li>
  <li><i>&beta;<sub>s</sub></i> is the effect of study cohort,</li>
  <li><i>s<sub>i</sub></i> is the study cohort for mouse <i>i</i>,</li>
  <li><i>&beta;<sub>j</sub></i> is the effect of founder allele <i>j</i>,</li>
  <li><i>g<sub>ij</sub></i> is the probability that mouse <i>i</i> carries an allele from founder <i>j</i>,</li>
  <li><i>&lambda;<sub>i</sub></i> is an adjustment for kinship-induced correlated errors for mouse <i>i</i>,</li>
  <li><i>&epsilon;<sub>i</sub></i> is the residual error for mouse <i>i</i>.</li>
  </ul>  

Note that this model will give us an estimate of the effect of each founder allele at each marker. There are eight founder strains that contributed to the DO, so we will get eight founder allele effects.

There are almost 600 samples in this data set and it may take several minutes to map one trait. In order to save some time, we will map using only the samples in the 100 ppm concentration group.


~~~
c100 = which(pheno$Conc == 100)
~~~
{: .r}



~~~
Error in which(pheno$Conc == 100): object 'pheno' not found
~~~
{: .error}

In order to map the proportion of bone marrow reticulocytes that were micro-nucleated, you will use the [scan1](https://github.com/rqtl/qtl2/blob/master/R/plot_scan1.R) function. To see the arguments for [scan1](https://github.com/rqtl/qtl2/blob/master/R/plot_scan1.R), you can type `help(scan1)`. First, let's map the *untransformed* phenotype. (Recall that we log-transformed it above).



























