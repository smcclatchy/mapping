---
layout: page
title: Setup
permalink: /setup/
---
## Installation

To install R/qtl, type

~~~
install.packages("qtl")
~~~
{: .r}

or use the Install button on the Packages tab in RStudio.

A manual, tutorials, and a sample chapter from Broman & Sen's book [A Guide to QTL Mapping with R/qtl](http://www.rqtl.org/book/) are available from the [R/qtl website](http://www.rqtl.org/).

R/qtl2 is not yet available on [CRAN](https://cran.r-project.org), but
it can be installed from a mini-CRAN at [rqtl.org](http://rqtl.org).

~~~
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
~~~
{: .r}

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

Alternatively, you can install R/qtl2 from its source on
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
(Additional, secondary dependencies will also be installed.)

~~~
install.packages(c("devtools", "yaml", "jsonlite", "data.table", "RcppEigen"))
~~~
{: .r}

Finally, install R/qtl2 using `devtools::install_github()`.

~~~
library(devtools)
install_github("rqtl/qtl2")
install_github("rqtl/qtl2convert")
install_github("rqtl/qtl2db")
~~~
{: .r}


## Data file format

The recommended file format for R/qtl is comma-separated values (csv).
In the simplest case all phenotype and genotype data and the genetic
map of the genotyped markers are contained in a single csv file. 
The first columns should be phenotypes and can include animal ID 
numbers, sex, or treatment. The following columns contain marker
genotypes. The input data file formats for [R/qtl](http://rqtl.org) cannot
handle complex crosses, and so for R/qtl2, we have defined a new
format for the data files. We'll describe it here briefly; for
details, see the separate
[vignette on the input file format](http://kbroman.org/qtl2/assets/vignettes/input_files.html).

QTL mapping data consists of a set of tables of data: marker
genotypes, phenotypes, marker maps, etc. In the new format, these
different tables are in separate comma-delimited (CSV) files. In each
file, the first column is a set of IDs for the rows, and the first row
is a set of IDs for the columns. For example, the phenotype data file
will have individual IDs in the first column and phenotype names in
the first row.

A few important changes in the tabular data:

- We will use not just the genetic marker map, but also a physical map
(if available).
- Previously, phenotypes and covariates were combined. In the new
format, we separate numeric phenotypes from the often
non-numeric covariates.
- We define a table of &ldquo;phenotype covariates.&rdquo; These are
[metadata](https://en.wikipedia.org/wiki/Metadata) describing the
phenotypes. For example, in the case of a phenotype measured over
time, one column in the phenotype covariate data could be the
time of measurement. For gene expression data, we would have columns
representing chromosome and physical position of genes, as well as
gene IDs.

In addition to the set of CSV files with the primary data, we need a
separate &ldquo;control&rdquo; file with various control parameters
(or metadata), including the names of all of the other data files and
the genotype codes used in the genotype data file. The control file is
in a specific format using either [YAML](http://www.yaml.org) or
[JSON](http://json.org); these are human-readable text files for
representing relatively complex data.

A big advantage of this control file scheme is that it greatly
simplifies the function for reading in the data. That function,
`read_cross2()`, has a _single_ argument: the name (with path) of the
control file. So you can read in data like this:

~~~
library(qtl2geno)
grav2 <- read_cross2("~/my_data/grav2.yaml")
~~~
{: .r}

_Note_: it can sometimes confusing which packages you need to load. If
you install the (largely empty) [qtl2](https://github.com/rqtl/qtl2)
package, you can type just

~~~
library(qtl2)
~~~
{: .r}

and the three main packages,
[qtl2geno](https://github.com/rqtl/qtl2geno),
[qtl2scan](https://github.com/rqtl/qtl2scan), and
[qtl2plot](https://github.com/rqtl/qtl2plot), will be loaded.

The large number of files is a bit cumbersome, so we've made it
possible to use a
[zip file](https://en.wikipedia.org/wiki/Zip_(file_format)) containing
all of the data files, and to read that zip file directly. There's even a
function for creating the zip file:

~~~
zip_datafiles("~/my_data/grav2.yaml")
~~~
{: .r}

This `zip_datafiles()` function will read the control file to identify
all of the relevant data files and then zip them up into a file with
the same name and location, but with the extension `.zip` rather than
`.yaml` or `.json`.

To read the data back in, we use the same `read_cross2()` function,
providing the name (and path) of the zip file rather than the control
file.

~~~
grav2 <- read_cross2("~/my_data/grav2.zip")
~~~
{: .r}

This can even be done with remote files.

~~~
grav2 <- read_cross2("http://kbroman.org/qtl2/assets/sampledata/grav2/grav2.zip")
~~~
{: .r}

Of course, the other advantage of the zip file is that it is
_compressed_ and so smaller than the combined set of CSV files.

The control file may be confusing for some users. To assist in its
construction, there's a function `write_control_file()` that takes the
large set of control parameters as input and then writes the YAML
control file in the appropriate format.


