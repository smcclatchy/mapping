---
title: "QTL analysis in Diversity Outbred Mice"
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



DO QTL Mapping Tutorial
========================================================
Introduction
--------------------------------------------------------

This tutorial will take you through the process of mapping a QTL and searching for candidate genes. It assumes that you have already phenotyped and genotyped your mice, and have reconstructed their genomes using the DOQTL calc.genoprob() function. DOQTL runs in the R statistical package and this demo assumes some familiarity with that software. You should have been provided the following files for the demo:

* DOQTL_demo.Rdata: contianing the phenotype and genotype data.

The data comes from a toxicology study in which mice were exposed to benzene via inhalation for 6 hours a day, 5 days a week for 4 weeks. The study was conducted in two equally sized cohort of 300 male mice each, for a total of 600 mice. They were then sacrificed and reticulocytes (red blood cell precursors) were isolated from bone marrow. The number of micro-nucleated reticulocytes, a measure of DNA damage, was then measured in each mouse. The goal is to map gene(s) that influence the level of DNA damage in the bone marrow.

![](figure/benzene_study_design.png)

Loading the DOQTL package
--------------------------------------------------------

The following code loads DOQTL's dependencie and installs DOQTL.  You do not need to run it now, but this is what you would run to install DOQTL on your computer at home.  Note that the installation of dependencies takes a while becuase the mouse and human genomes are downloaded.

   source("http://bioconductor.org/biocLite.R")  
   biocLite("DOQTL")

First, we load DOQTL into the R environment.


~~~
library(DOQTL)
~~~
{: .r}

Loading the data
--------------------------------------------------------

The data for this tutorial has been saved as an R binary file that contains several data objects.  Load it in now by running the following command.


~~~
load("/data/DOQTL_demo.Rdata")
~~~
{: .r}



~~~
Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
'/data/DOQTL_demo.Rdata', probable reason 'No such file or directory'
~~~
{: .error}



~~~
Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
~~~
{: .error}

This loaded in two data objects. Look in the Global Environment panel to see what was loaded.  You should see an obejct called `pheno` with 143 rows and 5 columns and an object called `probs`.

`pheno` is a data frame containing the phenotype data. `probs` is a 3 dimensional array containing the founder allele dosages for each sample at each marker on the array.  Double-click on `pheno` in the Global Environment panel to view its contents.

**NOTE:** the sample IDs must be in the rownames of `pheno`.

It contains the sample ID, the study cohort, the dose of benzene and the proportion of bone marrow reticulocytes that were micro-nucleated (prop.bm.MN.RET).  Note that the sample IDs are also stored in the rownames of pheno. In order to save time for this tutorial, we will only map with 143 samples from the 100 ppm dosing group.

Next, we look at the contents of `probs`:


~~~
dim(probs)
~~~
{: .r}



~~~
Error in eval(expr, envir, enclos): object 'probs' not found
~~~
{: .error}

`probs` is a three dimensional array containing the proportion of each founder haplotype at each marker for each DO sample.  The 143 samples are in the first dimension, the 8 founders in the second and the markers along the mouse genome are in the third dimension. Let's look at the contents for the first 500 markers of one sample.

**NOTE:** the sample IDs must be in the rownames of `probs`.


~~~
image(1:500, 1:ncol(probs), t(probs[1,8:1,1:500]), breaks = 0:100/100,
      col = grey(99:0/100), axes = F, xlab = "Markers", ylab = "Founders",
      main = "Founder Allele Contributions for Sample 1")
~~~
{: .r}



~~~
Error in ncol(probs): object 'probs' not found
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

In the plot above, the founder contributions, which range between 0 and 1, are colored from white (= 0) to black (= 1.0). A value of ~0.5 is grey. The markers are on the X-axis and the eight founders (denoted by the letters A through H) on the Y-axis. Starting at the left, we see that this sample has genotype CD because both rows C and D are grey, indicating values o 0.5 for each one. Moving along the genome to the right, the genotype becomes DD where row D is black, then CD, AC, CH, CD, CH, etc. The value at each marker sum to 1.0.

QTL Mapping
-------------------------------------------------------

First, we need the locations of the markers on the genotyping array. The array is called the Mouse Universal Genotyping Array (MUGA) and contain 7,856 SNP markers. Their locations are on [The Jackson Laboratory's FTP site](ftp://ftp.jax.org/MUGA):


~~~
load(url("ftp://ftp.jax.org/MUGA/muga_snps.Rdata"))
~~~
{: .r}

Next, we need to create a matrix that accounts for the kinship relationships between the mice. We do this by looking at the correlation between the founder haplotypes for each sample at each SNP. For each chromosome, we create a kinship matrix using the all markers *except* the ones on the current chromosome. Simulations suggest that mapping using this approach increases the power to detect QTL.
           

~~~
K = kinship.probs(probs, snps = muga_snps, bychr = TRUE)
~~~
{: .r}



~~~
Error in kinship.probs(probs, snps = muga_snps, bychr = TRUE): could not find function "kinship.probs"
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



~~~
axis(side = 2, at = 20 * 0:7, labels = 20 * 7:0, las = 1)
~~~
{: .r}



~~~
Error in axis(side = 2, at = 20 * 0:7, labels = 20 * 7:0, las = 1): plot.new has not been called yet
~~~
{: .error}

The figure above shows kinship between all pairs of samples. White (= 1) indicates no kinship and red (= 0) indicates full kinship. Orange values indicate varying levels of kinship between 0 and 1. The white diagonal of the matrix indicates that each sample is identical to itself. The lighter yellow blocks off of the diagonal may indicate siblings or cousins.

Next, we need to create additive covariates that wil be used in the mapping model. We will use sex and study cohort as a covariate in the mapping model.  While all of the samples We must add the sample IDs to the rownames of the covariates becuase the 'scanone' function will match up sample IDs in all of the data.


~~~
addcovar = model.matrix(~Study, data = pheno)
~~~
{: .r}



~~~
Error in terms.formula(object, data = data): object 'pheno' not found
~~~
{: .error}



~~~
colnames(addcovar)[1] = "Sex"
~~~
{: .r}



~~~
Error in colnames(addcovar)[1] = "Sex": object 'addcovar' not found
~~~
{: .error}

The code above copies the rownames(pheno) to rownames(addcovar).

**NOTE:** the sample IDs must be in the rownames of `addcovar`.

In order to map prop.bm.MN.RET, you will use the scanone() function. To see the arguments for `scanone`, you can type 'help(scanone)'.





























