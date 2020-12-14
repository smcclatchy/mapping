# this script demonstrates how to use the qtl2 package
# for mapping quantitative traits
# Sue McClatchy
# 2020-12-14

# load the library after installing
library(qtl2)

# load the iron data
iron <- read_cross2(
  file = system.file(
    "extdata", "iron.zip", package="qtl2") )

# get an overview of the data
summary(iron)
names(iron)

# look at the first few entries in the genetic map
head(iron$gmap)

# create a new map by inserting pseudomarkers at 1 cM
# intervals between markers in the genetic map
map <- insert_pseudomarkers(map=iron$gmap, step=1)

# check the first 2 chromosomes in the new map
head(map, n=2)

# calculate genotype probabilities for each pseudomarker
# assuming an error rate of 0.002 (default is 0.0001)
pr <- calc_genoprob(
  cross=iron, map=map, error_prob=0.002)

# calc_genoprob() produces a 3-dimensional array 
# of genotype probabilities for each chromosome
names(pr)

# look at the 3 dimensions for chromosome 19
dimnames(pr$`19`)

# view the genotype probabilities on chromosome 19 
# for a genotyped marker and the two adjacent pseudomarkers 
# located at 1 cM intervals away for the first three individuals. 
# Compare the probabilities for each 
# pseudomarker genotype with those of the genotyped marker.
(pr$`19`)[1:3,,"D19Mit68"] # genotyped marker
(pr$`19`)[1:3,,"c19.loc4"] # pseudomarker 1 cM away
(pr$`19`)[1:3,,"c19.loc5"] # the next pseudomarker

# visualize the genotype probabilities for chromosome 19
# for individual 1
plot_genoprob(pr, map, ind = 1, chr = 19)
