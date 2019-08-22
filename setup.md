---
layout: page
title: Setup
permalink: /setup/
---
## Installation

R is a programming language that is especially powerful for data exploration, visualization, and statistical analysis. To interact with R, we use RStudio. 

1. Install the latest version of R from [CRAN](https://cran.r-project.org/).

2. Install the latest version of RStudio [here](https://www.rstudio.com/products/rstudio/download/). Choose the free RStudio Desktop version for Windows, Mac, or Linux. 

3. Start RStudio. The [qtl2](https://github.com/rqtl/qtl2) package contains code for haplotype reconstruction, QTL mapping and plotting. Install qtl2 by running the following code in the R console.

~~~
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
~~~
{: .r}

Make sure that the installation was successful by loading the qtl2 library. You shouldn't get any error messages.

~~~
library(qtl2)
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


You will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice.


Make sure that both the SNP and gene files downloaded correctly by running the following code. If you get an error, check the file path carefully or download the files again. Make sure to change the file path to the location where you saved the file.


Check the SNP file.

~~~
snp_func = create_variant_query_func("~/Desktop/mapping/data/cc_variants.sqlite")
snps = snp_func(1, 10, 11)
dim(snps)
~~~
{: .r}


You should get a result that is 13150 rows by 16 columns.


Check the gene file.


~~~
gene_func = create_gene_query_func("~/Desktop/mapping/data/mouse_genes_mgi.sqlite")
genes = gene_func(1, 10, 11)
dim(genes)
~~~
{: .r}


You should get a result that is 18 rows by 15 columns.
