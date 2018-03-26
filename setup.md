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

The [qtl2](https://github.com/rqtl/qtl2) package contains code for haplotype reconstruction, QTL mapping and plotting.

Option 1: R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), but it can be installed from a mini-CRAN at [rqtl.org](http://www.rqtl.org/). Make sure you have the latest version of [R (3.4.4)](https://cran.r-project.org).

Option 2: Alternatively, you can install R/qtl2 from its source on [GitHub](https://github.com/rqtl).
(But note that compiling the C++ code can be rather slow.)

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
install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen", "RSQLite", "qtl"))
~~~
{: .r}

Finally, install R/qtl2 using `devtools::install_github()`. Copy and paste the following code into the R console in RStudio.

~~~
library(devtools)
install_github("rqtl/qtl2")
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

3. Please download the following large files **before the workshop**, and place them in your `data` folder. You can download the files from the URLs below and move the files the same way that you would for downloading and moving any other kind of data.

- [cc_variants.sqlite doi:10.6084/m9.figshare.5280229.v2](https://figshare.com/articles/SQLite_database_of_variants_in_Collaborative_Cross_founder_mouse_strains/5280229/2), variants in the Collaborative Cross founders (3 GB)
- [mouse_genes.sqlite doi:10.6084/m9.figshare.5280238.v4](https://figshare.com/articles/SQLite_database_with_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5280238/4) full set of mouse gene annotations (677 MB)
- [mouse_genes_mgi.sqlite doi:10.6084/m9.figshare.5286019.v5](https://figshare.com/articles/SQLite_database_with_MGI_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5286019/5) just the MGI mouse gene annotations (11 MB)
- [DO QTL data](ftp://ftp.jax.org/dgatti/qtl2_workshop/qtl2_demo.Rdata) Dan's DO data from the benzene study

Alternatively, you can copy and paste the following into the R console to download the data.
~~~
setwd("~/Desktop/mapping")
download.file("https://ndownloader.figshare.com/files/9746485", "./data/cc_variants.sqlite")
download.file("https://ndownloader.figshare.com/files/9746458", "./data/mouse_genes.sqlite")
download.file("https://ndownloader.figshare.com/files/9746452", "./data/mouse_genes_mgi.sqlite")
download.file("ftp://ftp.jax.org/dgatti/qtl2_workshop/qtl2_demo.Rdata", "./data/qtl2_demo.Rdata")
~~~
{: .r}


 You will
need these for the final lesson episode on QTL analysis in Diversity Outbred mice.


