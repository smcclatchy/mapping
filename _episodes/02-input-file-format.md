---
title: "Input File Format"
teaching: 15
exercises: 30
questions:
- "How are the data files formatted for qtl2?"
- "Which data files are required for qtl2?"
- "How do input files compare between qtl and qtl2?"
objectives:
- To understand the which input files are required for qtl2 and how they should be formatted.
- To compare input files between qtl and qtl2.
keypoints:
- "QTL mapping data consists of a set of tables of data: marker genotypes, phenotypes, marker maps, etc."
- "These different tables are in separate comma-delimited (CSV) files."
- "In each file, the first column is a set of IDs for the rows, and the first row is a set of IDs for the columns."
- "In addition to primary data, a separate file with control parameters (or metadata) in either [YAML](http://www.yaml.org) or [JSON](http://json.org) format is required."
source: Rmd
---



QTL mapping data consists of a set of tables of data: marker
genotypes, phenotypes, marker maps, etc. These different tables are in different comma-separated value (CSV) files. In each file, the first column is a set of IDs for the rows, and the first row is a set of IDs for the columns. For example, the phenotype data file will have individual IDs in the first column and phenotype names in the first row.

R/qtl2 accepts the following files:
1. genotypes
2. phenotypes
3. phenotype covariates (*i.e.* tissue type, time points)  
4. genetic map  
5. physical map (optional)  
6. control file (YAML or JSON format, not CSV)

We use both a genetic marker map and a physical map (if available). Numeric phenotypes are separate from the often non-numeric covariates. Phenotype covariates are [metadata](https://en.wikipedia.org/wiki/Metadata) describing the
phenotypes. For example, in the case of a phenotype measured over
time, one column in the phenotype covariate data could be the
time of measurement. For gene expression data, we would have columns
representing chromosome and physical position of genes, as well as
gene IDs.

In addition to the set of CSV files with the primary data, we need a
separate control file with various control parameters
(or metadata), including the names of all of the other data files and
the genotype codes used in the genotype data file. The control file is
in a specific format using either [YAML](http://www.yaml.org) or
[JSON](http://json.org); these are human-readable text files for
representing relatively complex data.

A big advantage of this control file scheme is that it greatly
simplifies the function for reading in the data. That function,
`read_cross2()`, has a _single_ argument: the name (with path) of the
control file.

For further details, see the separate [vignette on the input file format](http://kbroman.org/qtl2/assets/vignettes/input_files.html).

In this lesson, we'll work with data sets included in the `qtl2` package. Additional sample data sets, including data on Diversity Outbred (DO) mice, are available at <https://github.com/rqtl/qtl2data>.

> ## Challenge 1
> Which data files are required by `qtl2`?  
> Which ones are optional?  
> How should they be formatted?
>
> > ## Solution to Challenge 1
> >
> {: .solution}
{: .challenge}

> ## Challenge 2
> Go to <https://github.com/rqtl/qtl2data> to view additional sample data.
> 1). Find the Recla data and locate the phenotype data file. Open the file. What is in the first column? the first row?  
> 2). Locate the genotype data file and open it.  
> What is in the first column? the first row?  
> 3). Locate the control file (YAML or JSON format) and open it.  
> What kind of information does this file contain?
> 4). Locate the phenotype covariates file and open it.  
> What kind of information does this file contain?
>
> > ## Solution to Challenge 2
> >
> {: .solution}
{: .challenge}






