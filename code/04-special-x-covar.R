## ---- include=FALSE------------------------------------------------------
source("../bin/chunk-options.R")
knitr_fig_path("04-")


## ----load_data, echo=FALSE-----------------------------------------------
library(qtl2)
iron <- read_cross2(file = system.file("extdata", "iron.zip", package="qtl2") )


## ----x_covar-------------------------------------------------------------
Xcovar <- get_x_covar(iron)
head(Xcovar)

