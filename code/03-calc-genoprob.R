## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("03-")


## ----load_data-----------------------------------------------------------
library(qtl2)
iron <- read_cross2(file = system.file("extdata", "iron.zip", package="qtl2") )


## ----load_my_data, eval=FALSE, error=FALSE-------------------------------
## myQTLdata <- read_cross2(file = "/Users/myUserName/qtlProject/data/myqtldata.yaml" )


## ----load_my_zipdata, eval=FALSE, error=FALSE----------------------------
## myQTLdata <- read_cross2(file = "/Users/myUserName/qtlProject/data/myqtldata.zip" )


## ----summary_data--------------------------------------------------------
summary(iron)
names(iron)


## ----map_data------------------------------------------------------------
head(iron$gmap)


## ----insert_pseudomarkers------------------------------------------------
map <- insert_pseudomarkers(map=iron$gmap, step=1)


## ----view_map------------------------------------------------------------
head(map, n=2)


## ----calc_genoprob-------------------------------------------------------
pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002)


## ----list_chrs-----------------------------------------------------------
names(pr)


## ----view_array----------------------------------------------------------
dimnames(pr$`19`)


## ----view_genoprob-------------------------------------------------------
(pr$`19`)[1:3,,"D19Mit68"] # genotyped marker
(pr$`19`)[1:3,,"c19.loc4"] # pseudomarker 1 cM away
(pr$`19`)[1:3,,"c19.loc5"] # the next pseudomarker


## ----calc_genoprob_multicore, eval=FALSE---------------------------------
## pr <- calc_genoprob(cross=iron, map=map, error_prob=0.002, cores=4)


## ----allele_probs--------------------------------------------------------
apr <- genoprob_to_alleleprob(probs=pr)

