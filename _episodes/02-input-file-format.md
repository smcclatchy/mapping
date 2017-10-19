---
title: "Input File Format"
teaching: 0
exercises: 0
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



If you are accustomed to input files for [R/qtl](http://rqtl.org), you will find that there are some similarities between formats for R/qtl and R/qtl2, and some important differences. The input data file formats for R/qtl cannot handle complex crosses, and so for R/qtl2, there is a new format for the data files. The table below provides a quick comparison between data files and file formats for `qtl` and `qtl2`.

| Data                      | R/qtl     | R/qtl2                  |  
|:--------------------------|:----------|:------------------------|
| Format                    |  CSV      |  CSV                    | 
| Marker map                |  genetic  |  genetic, physical      | 
| Phenotypes and covariates |  separate |  combined               | 
| Control file              |  none     |  required (YAML or JSON)| 

R/qtl2 accepts the following files:
1. genotypes (CSV)  
2. phenotypes (CSV)  
3. phenotype covariates (*i.e.* tissue type, time points as CSV)  
4. genetic map (CSV)  
5. physical map (CSV; optional)  
6. control file (YAML or JSON)


In this lesson, we'll work with data sets included in the `qtl2` package. Additional sample data sets, including data on Diversity Outbred (DO)
mice, are available at <https://github.com/rqtl/qtl2data>.

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
> 1). Find the Recla data and locate the phenotype data file. Open the file.  > What is in the first column? the first row?  
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



