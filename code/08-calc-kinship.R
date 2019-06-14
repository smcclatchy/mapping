## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("08-")
library()


## ----load_data, echo=FALSE-----------------------------------------------
library(qtl2)
iron <- read_cross2(file = system.file("extdata", "iron.zip", package="qtl2") )
map <- insert_pseudomarkers(map = iron$gmap, step = 1)
pr <- calc_genoprob(cross = iron, map = map, error_prob = 0.002)


## ----table, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
|                              |subpop1|subpop2|overall pop
|:-----------------------------|:-----:|:-----:|:-----:|
| frequency                    |  0.5  |  0.5  |   1   |
| probability of AA genotype   |  0.1  |  0.9  |  0.5  |
| probability of disease       |  0.9  |  0.1  |  0.5  |
| probability of disease & AA  |  0.09 |  0.09 |  0.09 |
"
cat(tabl) 


## ----calc_kinship--------------------------------------------------------
kinship <- calc_kinship(probs = pr)


## ----view_kinship--------------------------------------------------------
kinship[1:5, 1:5]


## ----calc_kinship_grid, eval=FALSE---------------------------------------
## grid <- calc_grid(map = iron$gmap, step=1)
## pr_grid <- probs_to_grid(probs = pr, grid = grid)
## kinship_grid <- calc_kinship(probs = pr_grid)


## ----calc_kinship_loco_multicore, eval=FALSE, eval=FALSE-----------------
## kinship <- calc_kinship(pr, cores=4)

