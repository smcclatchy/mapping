---
layout: page
title: Setup
permalink: /setup/
---
## Installation

R is a programming language that is especially powerful for data exploration, visualization, and statistical analysis. To interact with R, we use RStudio. 

1. Install the latest version of R from [CRAN](https://cran.r-project.org/).

2. Install the latest version of RStudio [here](https://www.rstudio.com/products/rstudio/download/). Choose the free RStudio Desktop version for Windows, Mac, or Linux. 

3. Start RStudio. The [qtl2](https://github.com/rqtl/qtl2) package contains code for haplotype reconstruction, QTL mapping and plotting. Install qtl2 by copying and pasting the following code in the R console.

~~~
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
~~~
{: .r}

Make sure that the installation was successful by loading the qtl2 library, either by copy-paste to the Console or by checking the box next to `qtl2` in the RStudio Packages tab. You shouldn't get any error messages.

~~~
library(qtl2)
~~~
{: .r}

## Data files and project organization

1. Create a new project in your Desktop called `mapping`. 
- Click the `File` menu button, then `New Project`.
- Click `New Directory`. 
- Click `New Project`.
- Type `mapping` as the directory name. Browse to your Desktop to create the project there.
- Click the `Create Project` button.

2. Use the `Files` tab to create  a `data` folder to hold the data, a `scripts` folder to house your scripts, and a `results` folder to hold results. 

Alternatively, you can use the R console to run the following commands for step 2.

~~~
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
~~~
{: .r}

3. Please download the following large files **before the workshop**, and place them in your `data` folder. You can download the files from the URLs below and move the files the same way that you would for downloading and moving any other kind of data.
 
+ [SQLite database of variants in Collaborative Cross founder mouse strains (v3)](https://figshare.com/articles/SQLite_database_of_variants_in_Collaborative_Cross_founder_mouse_strains/5280229/3): SNP, indel, and structural variants in the Collaborative Cross founders (4.91 GB)
+ [SQLite database with mouse gene annotations from Mouse Genome Informatics (v7)](https://figshare.com/articles/dataset/SQLite_database_with_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5280238): full set of mouse gene annotations from build 38 mm10 (1.32 GB)
+ [SQLite database with MGI mouse gene annotations from Mouse Genome Informatics (v8)](https://figshare.com/articles/dataset/SQLite_database_with_MGI_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5286019): like the previous, but including only non-duplicate gene records sourced from MGI (11.36 MB)
+ [DO QTL data](ftp://ftp.jax.org/dgatti/qtl2_workshop/qtl2_demo.Rdata) from benzene study described in [French, John E., et al. Env Health Perspectives (2015): 237-245.](https://ehp.niehs.nih.gov/doi/10.1289/ehp.1408202) (240.8 MB)

Alternatively, you can copy and paste the following into the R console to download the data.
~~~
setwd("~/Desktop/mapping")
download.file("https://ndownloader.figshare.com/files/18533342", "./data/cc_variants.sqlite")
download.file("https://ndownloader.figshare.com/files/17609261", "./data/mouse_genes.sqlite")
download.file("https://ndownloader.figshare.com/files/17609252", "./data/mouse_genes_mgi.sqlite")
download.file("ftp://ftp.jax.org/dgatti/qtl2_workshop/qtl2_demo.Rdata", "./data/qtl2_demo.Rdata")
~~~
{: .r}


You will need these for the final lesson episodes on SNP association mapping and QTL analysis in Diversity Outbred mice.


Make sure that both the SNP and gene files downloaded correctly by running the following code. If you get an error, check the file path (e.g. "~/Desktop/mapping/data/cc_variants.sqlite") carefully or download the files again. Make sure to change the file path to the location where you saved the file.


Check part of the SNP file. It is a very large file, so checking only a sample of the file should do.

~~~
# create a function to query the SNP file, then use this new function  
# to select SNPs on chromosome 1 from 10 to 11 Mbp
snp_func = create_variant_query_func("~/Desktop/mapping/data/cc_variants.sqlite") 
snps = snp_func(chr = 1, start = 10, stop = 11) 

# check the dimensions of this sample of the SNP file
dim(snps) 
~~~
{: .r}


You should get a result that is 13150 rows by 16 columns.


Check the gene file in the same way.


~~~
# create a function to query the gene file, then select genes in the same region as above
gene_func = create_gene_query_func("~/Desktop/mapping/data/mouse_genes_mgi.sqlite") 
genes = gene_func(chr = 1, start = 10, stop = 11) 
dim(genes) # check the dimensions
~~~
{: .r}


You should get a result that is 18 rows by 15 columns.
