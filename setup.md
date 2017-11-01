---
layout: page
title: Setup
permalink: /setup/
---
## Installation

This lesson assumes you have the R and RStudio software installed on your computer.

R can be downloaded from [CRAN](https://cran.r-project.org/mirrors.html).

RStudio is an environment for developing using R.
It can be downloaded [here](https://www.rstudio.com/products/rstudio/download/).
You will need the Desktop version for your computer.

The [qtl2](https://github.com/rqtl/qtl2) package is
inspired by the
[tidyverse package](https://cran.r-project.org/package=tidyverse);
it is basically empty, but when you install it, the
[qtl2geno](https://github.com/rqtl/qtl2geno),
[qtl2scan](https://github.com/rqtl/qtl2scan),
[qtl2plot](https://github.com/rqtl/qtl2plot),
[qtl2db](https://github.com/rqtl/qtl2db), and
[qtl2convert](https://github.com/rqtl/qtl2convert) packages, plus a
bunch of dependencies, will be installed.

R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), and
is still under development. To retrieve the latest version, install R/qtl2 from its source on
[GitHub](https://github.com/rqtl). (But note that compiling the C++
code can be rather slow.)

On _Windows_, you'll need [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

On _Mac OS X_, you'll need the
[command-line developer tools](https://mac-how-to.gadgethacks.com/how-to/install-command-line-developer-tools-without-xcode-0168115/),
as well as [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS).

You then need to install the
[devtools](https://github.com/hadley/devtools) package, plus a set of
package dependencies: [yaml](https://cran.r-project.org/package=yaml),
[jsonlite](https://cran.r-project.org/package=jsonlite),
[data.table](https://cran.r-project.org/package=data.table),
and [RcppEigen](https://github.com/RcppCore/RcppEigen).
(Additional, secondary dependencies will also be installed.) Start RStudio, then copy and paste the following code into the R console in RStudio.

~~~
install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))
~~~
{: .r}

Finally, install R/qtl2 using `devtools::install_github()`. Copy and paste the following code into the R console in RStudio.

~~~
library(devtools)
install_github("rqtl/qtl2")
install_github("rqtl/qtl2convert")
install_github("rqtl/qtl2db")
~~~
{: .r}

## Data files and project organization

1. Make a new folder in your Desktop called `mapping`. Move into this new folder.

2. Create  a `data` folder to hold the data, a `scripts` folder to house your scripts, and a `results` folder to hold results. 

Alternatively, you can use the R console to run the following commands for steps 1 and 2.

~~~
setwd("~/Desktop")
dir.create("./mapping")
setwd("~/Desktop/mapping")
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
~~~
{: .r}

3. Please download the following large files **before the workshop**, and place them in your `data` folder. You will 
need these for the final lesson episode on QTL analysis in Diversity Outbred mice.

[DO QTL data](ftp://ftp.jax.org/dgatti/MDIBL_Aging2016/DOQTL_demo.Rdata)

[MUGA SNPs](ftp://ftp.jax.org/MUGA/muga_snps.Rdata)

[CC Variants](ftp://ftp.jax.org/dgatti/CC_SNP_DB/cc_variants.sqlite)

[Mouse Genes](ftp://ftp.jax.org/dgatti/CC_SNP_DB/mouse_genes.sqlite)
