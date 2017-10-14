---
title: "Calculating Genotype Probabilities"
teaching: 0
exercises: 0
questions:
- "How do I calculate QTL at positions between genotyped markers?"
- "How do I calculate QTL genotype probabilities?"
- "How do I calculate allele probabilities?"
- "How can I speed up calculations if I have a large data set?"
objectives:
- To explain why the first step in QTL analysis is to calculate genotype probabilities.
- To insert pseudomarkers between genotyped markers.
- To calculate genotype probabilities.
keypoints:
- "The first step in QTL analysis is to calculate genotype probabilities."
- "Insert pseudomarkers to calculate QTL at positions between genotyped markers."
- "Calculate genotype probabilities between genotyped markers with calc_genoprob()."
source: Rmd
---



The first basic task in QTL analysis is to calculate conditional genotype probabilities, given the observed marker data, at each putative QTL position. For example, the first step would be to determine the probabilities for genotypes AA and AB at the locus indicated below.

![](../fig/unknown_genotype.png)

Hidden Markov models (HMM) are used to calculate
QTL genotype probabilities, to simulate from the joint genotype distribution and to calculate the most likely sequence of underlying genotypes (all conditional on the observed marker data). This is done in a quite general way, with possible allowance for the presence of genotyping errors. For convenience we assume no crossover interference.

The `calc_genoprob()` function in the [qtl2geno](https://github.com/rqtl/qtl2geno)
package calculates QTL genotype probabilities, conditional on the available marker data. These are needed for most of the QTL mapping functions. Unlike the corresponding function in
[R/qtl](http://rqtl.org), `calc.genoprob()`, the result is not inserted back into the input cross object, but is returned as a list of three-dimensional arrays (one per chromosome). Each 3d array of probabilities is arranged as individuals &times; genotypes &times; positions.

[image of 3d array w/probs, inds, etc]

If we wish to perform QTL calculations at positions between markers (so called "pseudomarkers"), we first need to insert such positions into the genetic map with the function `insert_pseudomarkers()`. Unlike [R/qtl], the map is kept separate from the genotype
probabilities.

We'll use the
[iron dataset](https://github.com/kbroman/qtl2/tree/gh-pages/assets/sampledata/iron)
from
[Grant et al. (2006) Hepatology 44:174-185](https://www.ncbi.nlm.nih.gov/pubmed/16799992)
(an intercross) as an example. We first load the data:


~~~
library(qtl2geno)
iron <- read_cross2( system.file("extdata", "iron.zip", package="qtl2geno") )
~~~
{: .r}

(_Note_: you can use `library(qtl2)` to load the
three main packages,
[qtl2geno](https://github.com/rqtl/qtl2geno),
[qtl2scan](https://github.com/rqtl/qtl2scan), and
[qtl2plot](https://github.com/rqtl/qtl2plot), all at once.)


~~~
summary(iron)
~~~
{: .r}



~~~
Object of class cross2 (crosstype "f2")

Total individuals            284
No. genotyped individuals    284
No. phenotyped individuals   284
No. with both geno & pheno   284

No. phenotypes                 2
No. covariates                 2
No. phenotype covariates       1

No. chromosomes               20
Total markers                 66

No. markers by chr:
 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  X 
 3  5  2  2  2  2  7  8  5  2  7  2  2  2  2  5  2  2  2  2 
~~~
{: .output}



~~~
names(iron)
~~~
{: .r}



~~~
 [1] "crosstype"  "geno"       "gmap"       "pmap"       "pheno"     
 [6] "covar"      "phenocovar" "is_x_chr"   "is_female"  "cross_info"
[11] "alleles"   
~~~
{: .output}

Have a look at the markers listed in the genetic map, `gmap`. Markers are listed by chromosome and described by cM position.


~~~
iron$gmap
~~~
{: .r}



~~~
$`1`
D1Mit18 D1Mit80 D1Mit17 
   27.3    51.4   110.4 

$`2`
D2Mit379  D2Mit75  D2Mit17 D2Mit304  D2Mit48 
    38.3     48.1     56.8     58.0     73.2 

$`3`
D3Mit22 D3Mit18 
   25.1    54.6 

$`4`
  D4Mit2 D4Mit352 
    10.9     53.6 

$`5`
D5Mit11 D5Mit30 
   17.5    62.3 

$`6`
D6Mit104  D6Mit15 
    41.5     66.7 

$`7`
D7Mit74 D7Mit25  D7Nds5 D7mit30 D7Mit31 D7Mit17 D7Mit71 
    1.1    13.1    23.0    28.4    31.7    37.2    53.6 

$`8`
D8Mit124   D8Mit4 D8Mit195  D8Mit31 D8Mit294  D8Mit40 D8Mit120  D8Mit36 
     0.0     13.6     17.3     32.7     39.1     45.5     69.9     75.3 

$`9`
 D9Mit42  D9Mit31  D9Mit10 D9Mit182  D9Mit17 
     6.6     33.9     43.7     53.6     61.2 

$`10`
D10Mit61 D10Mit70 
    24.0     57.9 

$`11`
 D11Mit20   D11Mit4  D11Mit36  D11Mit41 D11Mit288   D8Mit18 D11Mit101 
      0.0      16.4      28.2      32.3      36.2      42.2      56.9 

$`12`
 D12Mit88 D12Mit134 
     19.7      57.9 

$`13`
D13Mit10 D13Mit51 
    17.5     40.4 

$`14`
 D14Mit54 D14Mit195 
     19.7      52.5 

$`15`
 D15Mit22 D15Mit159 
     16.4      49.2 

$`16`
D16Mit131   D16Mit4  D16Mit30  D16Mit19  D16Mit70 
      6.6      25.1      30.6      40.4      51.4 

$`17`
D17Mit46 D17Mit93 
     3.3     39.3 

$`18`
 D18Mit20 D18Mit186 
      4.4      30.6 

$`19`
D19Mit68 D19Mit37 
     3.3     38.3 

$X
 DXMit16 DXMit186 
    29.5     57.9 

attr(,"is_x_chr")
    1     2     3     4     5     6     7     8     9    10    11    12 
FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
   13    14    15    16    17    18    19     X 
FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
~~~
{: .output}

We then use `insert_pseudomarkers()` to insert pseudomarkers into the
genetic map, which we grab from the `iron` object as `iron$gmap`:


~~~
map <- insert_pseudomarkers(map=iron$gmap, step=1)
~~~
{: .r}

Now have a look at the new `map`.


~~~
map
~~~
{: .r}



~~~
$`1`
  D1Mit18  c1.loc28  c1.loc29  c1.loc30  c1.loc31  c1.loc32  c1.loc33 
     27.3      28.3      29.3      30.3      31.3      32.3      33.3 
 c1.loc34  c1.loc35  c1.loc36  c1.loc37  c1.loc38  c1.loc39  c1.loc40 
     34.3      35.3      36.3      37.3      38.3      39.3      40.3 
 c1.loc41  c1.loc42  c1.loc43  c1.loc44  c1.loc45  c1.loc46  c1.loc47 
     41.3      42.3      43.3      44.3      45.3      46.3      47.3 
 c1.loc48  c1.loc49  c1.loc50  c1.loc51   D1Mit80  c1.loc52  c1.loc53 
     48.3      49.3      50.3      51.3      51.4      52.3      53.3 
 c1.loc54  c1.loc55  c1.loc56  c1.loc57  c1.loc58  c1.loc59  c1.loc60 
     54.3      55.3      56.3      57.3      58.3      59.3      60.3 
 c1.loc61  c1.loc62  c1.loc63  c1.loc64  c1.loc65  c1.loc66  c1.loc67 
     61.3      62.3      63.3      64.3      65.3      66.3      67.3 
 c1.loc68  c1.loc69  c1.loc70  c1.loc71  c1.loc72  c1.loc73  c1.loc74 
     68.3      69.3      70.3      71.3      72.3      73.3      74.3 
 c1.loc75  c1.loc76  c1.loc77  c1.loc78  c1.loc79  c1.loc80  c1.loc81 
     75.3      76.3      77.3      78.3      79.3      80.3      81.3 
 c1.loc82  c1.loc83  c1.loc84  c1.loc85  c1.loc86  c1.loc87  c1.loc88 
     82.3      83.3      84.3      85.3      86.3      87.3      88.3 
 c1.loc89  c1.loc90  c1.loc91  c1.loc92  c1.loc93  c1.loc94  c1.loc95 
     89.3      90.3      91.3      92.3      93.3      94.3      95.3 
 c1.loc96  c1.loc97  c1.loc98  c1.loc99 c1.loc100 c1.loc101 c1.loc102 
     96.3      97.3      98.3      99.3     100.3     101.3     102.3 
c1.loc103 c1.loc104 c1.loc105 c1.loc106 c1.loc107 c1.loc108 c1.loc109 
    103.3     104.3     105.3     106.3     107.3     108.3     109.3 
c1.loc110   D1Mit17 
    110.3     110.4 

$`2`
D2Mit379 c2.loc39 c2.loc40 c2.loc41 c2.loc42 c2.loc43 c2.loc44 c2.loc45 
    38.3     39.3     40.3     41.3     42.3     43.3     44.3     45.3 
c2.loc46 c2.loc47  D2Mit75 c2.loc48 c2.loc49 c2.loc50 c2.loc51 c2.loc52 
    46.3     47.3     48.1     48.3     49.3     50.3     51.3     52.3 
c2.loc53 c2.loc54 c2.loc55 c2.loc56  D2Mit17 c2.loc57 D2Mit304 c2.loc58 
    53.3     54.3     55.3     56.3     56.8     57.3     58.0     58.3 
c2.loc59 c2.loc60 c2.loc61 c2.loc62 c2.loc63 c2.loc64 c2.loc65 c2.loc66 
    59.3     60.3     61.3     62.3     63.3     64.3     65.3     66.3 
c2.loc67 c2.loc68 c2.loc69 c2.loc70 c2.loc71 c2.loc72  D2Mit48 
    67.3     68.3     69.3     70.3     71.3     72.3     73.2 

$`3`
 D3Mit22 c3.loc26 c3.loc27 c3.loc28 c3.loc29 c3.loc30 c3.loc31 c3.loc32 
    25.1     26.1     27.1     28.1     29.1     30.1     31.1     32.1 
c3.loc33 c3.loc34 c3.loc35 c3.loc36 c3.loc37 c3.loc38 c3.loc39 c3.loc40 
    33.1     34.1     35.1     36.1     37.1     38.1     39.1     40.1 
c3.loc41 c3.loc42 c3.loc43 c3.loc44 c3.loc45 c3.loc46 c3.loc47 c3.loc48 
    41.1     42.1     43.1     44.1     45.1     46.1     47.1     48.1 
c3.loc49 c3.loc50 c3.loc51 c3.loc52 c3.loc53 c3.loc54  D3Mit18 
    49.1     50.1     51.1     52.1     53.1     54.1     54.6 

$`4`
  D4Mit2 c4.loc12 c4.loc13 c4.loc14 c4.loc15 c4.loc16 c4.loc17 c4.loc18 
    10.9     11.9     12.9     13.9     14.9     15.9     16.9     17.9 
c4.loc19 c4.loc20 c4.loc21 c4.loc22 c4.loc23 c4.loc24 c4.loc25 c4.loc26 
    18.9     19.9     20.9     21.9     22.9     23.9     24.9     25.9 
c4.loc27 c4.loc28 c4.loc29 c4.loc30 c4.loc31 c4.loc32 c4.loc33 c4.loc34 
    26.9     27.9     28.9     29.9     30.9     31.9     32.9     33.9 
c4.loc35 c4.loc36 c4.loc37 c4.loc38 c4.loc39 c4.loc40 c4.loc41 c4.loc42 
    34.9     35.9     36.9     37.9     38.9     39.9     40.9     41.9 
c4.loc43 c4.loc44 c4.loc45 c4.loc46 c4.loc47 c4.loc48 c4.loc49 c4.loc50 
    42.9     43.9     44.9     45.9     46.9     47.9     48.9     49.9 
c4.loc51 c4.loc52 c4.loc53 D4Mit352 
    50.9     51.9     52.9     53.6 

$`5`
   D5Mit11 c5.loc18.5 c5.loc19.5 c5.loc20.5 c5.loc21.5 c5.loc22.5 
      17.5       18.5       19.5       20.5       21.5       22.5 
c5.loc23.5 c5.loc24.5 c5.loc25.5 c5.loc26.5 c5.loc27.5 c5.loc28.5 
      23.5       24.5       25.5       26.5       27.5       28.5 
c5.loc29.5 c5.loc30.5 c5.loc31.5 c5.loc32.5 c5.loc33.5 c5.loc34.5 
      29.5       30.5       31.5       32.5       33.5       34.5 
c5.loc35.5 c5.loc36.5 c5.loc37.5 c5.loc38.5 c5.loc39.5 c5.loc40.5 
      35.5       36.5       37.5       38.5       39.5       40.5 
c5.loc41.5 c5.loc42.5 c5.loc43.5 c5.loc44.5 c5.loc45.5 c5.loc46.5 
      41.5       42.5       43.5       44.5       45.5       46.5 
c5.loc47.5 c5.loc48.5 c5.loc49.5 c5.loc50.5 c5.loc51.5 c5.loc52.5 
      47.5       48.5       49.5       50.5       51.5       52.5 
c5.loc53.5 c5.loc54.5 c5.loc55.5 c5.loc56.5 c5.loc57.5 c5.loc58.5 
      53.5       54.5       55.5       56.5       57.5       58.5 
c5.loc59.5 c5.loc60.5 c5.loc61.5    D5Mit30 
      59.5       60.5       61.5       62.3 

$`6`
  D6Mit104 c6.loc42.5 c6.loc43.5 c6.loc44.5 c6.loc45.5 c6.loc46.5 
      41.5       42.5       43.5       44.5       45.5       46.5 
c6.loc47.5 c6.loc48.5 c6.loc49.5 c6.loc50.5 c6.loc51.5 c6.loc52.5 
      47.5       48.5       49.5       50.5       51.5       52.5 
c6.loc53.5 c6.loc54.5 c6.loc55.5 c6.loc56.5 c6.loc57.5 c6.loc58.5 
      53.5       54.5       55.5       56.5       57.5       58.5 
c6.loc59.5 c6.loc60.5 c6.loc61.5 c6.loc62.5 c6.loc63.5 c6.loc64.5 
      59.5       60.5       61.5       62.5       63.5       64.5 
c6.loc65.5 c6.loc66.5    D6Mit15 
      65.5       66.5       66.7 

$`7`
 D7Mit74  c7.loc2  c7.loc3  c7.loc4  c7.loc5  c7.loc6  c7.loc7  c7.loc8 
     1.1      2.1      3.1      4.1      5.1      6.1      7.1      8.1 
 c7.loc9 c7.loc10 c7.loc11 c7.loc12  D7Mit25 c7.loc14 c7.loc15 c7.loc16 
     9.1     10.1     11.1     12.1     13.1     14.1     15.1     16.1 
c7.loc17 c7.loc18 c7.loc19 c7.loc20 c7.loc21 c7.loc22   D7Nds5 c7.loc23 
    17.1     18.1     19.1     20.1     21.1     22.1     23.0     23.1 
c7.loc24 c7.loc25 c7.loc26 c7.loc27 c7.loc28  D7mit30 c7.loc29 c7.loc30 
    24.1     25.1     26.1     27.1     28.1     28.4     29.1     30.1 
c7.loc31  D7Mit31 c7.loc32 c7.loc33 c7.loc34 c7.loc35 c7.loc36 c7.loc37 
    31.1     31.7     32.1     33.1     34.1     35.1     36.1     37.1 
 D7Mit17 c7.loc38 c7.loc39 c7.loc40 c7.loc41 c7.loc42 c7.loc43 c7.loc44 
    37.2     38.1     39.1     40.1     41.1     42.1     43.1     44.1 
c7.loc45 c7.loc46 c7.loc47 c7.loc48 c7.loc49 c7.loc50 c7.loc51 c7.loc52 
    45.1     46.1     47.1     48.1     49.1     50.1     51.1     52.1 
c7.loc53  D7Mit71 
    53.1     53.6 

$`8`
D8Mit124  c8.loc1  c8.loc2  c8.loc3  c8.loc4  c8.loc5  c8.loc6  c8.loc7 
     0.0      1.0      2.0      3.0      4.0      5.0      6.0      7.0 
 c8.loc8  c8.loc9 c8.loc10 c8.loc11 c8.loc12 c8.loc13   D8Mit4 c8.loc14 
     8.0      9.0     10.0     11.0     12.0     13.0     13.6     14.0 
c8.loc15 c8.loc16 c8.loc17 D8Mit195 c8.loc18 c8.loc19 c8.loc20 c8.loc21 
    15.0     16.0     17.0     17.3     18.0     19.0     20.0     21.0 
c8.loc22 c8.loc23 c8.loc24 c8.loc25 c8.loc26 c8.loc27 c8.loc28 c8.loc29 
    22.0     23.0     24.0     25.0     26.0     27.0     28.0     29.0 
c8.loc30 c8.loc31 c8.loc32  D8Mit31 c8.loc33 c8.loc34 c8.loc35 c8.loc36 
    30.0     31.0     32.0     32.7     33.0     34.0     35.0     36.0 
c8.loc37 c8.loc38 c8.loc39 D8Mit294 c8.loc40 c8.loc41 c8.loc42 c8.loc43 
    37.0     38.0     39.0     39.1     40.0     41.0     42.0     43.0 
c8.loc44 c8.loc45  D8Mit40 c8.loc46 c8.loc47 c8.loc48 c8.loc49 c8.loc50 
    44.0     45.0     45.5     46.0     47.0     48.0     49.0     50.0 
c8.loc51 c8.loc52 c8.loc53 c8.loc54 c8.loc55 c8.loc56 c8.loc57 c8.loc58 
    51.0     52.0     53.0     54.0     55.0     56.0     57.0     58.0 
c8.loc59 c8.loc60 c8.loc61 c8.loc62 c8.loc63 c8.loc64 c8.loc65 c8.loc66 
    59.0     60.0     61.0     62.0     63.0     64.0     65.0     66.0 
c8.loc67 c8.loc68 c8.loc69 D8Mit120 c8.loc70 c8.loc71 c8.loc72 c8.loc73 
    67.0     68.0     69.0     69.9     70.0     71.0     72.0     73.0 
c8.loc74 c8.loc75  D8Mit36 
    74.0     75.0     75.3 

$`9`
 D9Mit42  c9.loc8  c9.loc9 c9.loc10 c9.loc11 c9.loc12 c9.loc13 c9.loc14 
     6.6      7.6      8.6      9.6     10.6     11.6     12.6     13.6 
c9.loc15 c9.loc16 c9.loc17 c9.loc18 c9.loc19 c9.loc20 c9.loc21 c9.loc22 
    14.6     15.6     16.6     17.6     18.6     19.6     20.6     21.6 
c9.loc23 c9.loc24 c9.loc25 c9.loc26 c9.loc27 c9.loc28 c9.loc29 c9.loc30 
    22.6     23.6     24.6     25.6     26.6     27.6     28.6     29.6 
c9.loc31 c9.loc32 c9.loc33 c9.loc34  D9Mit31 c9.loc35 c9.loc36 c9.loc37 
    30.6     31.6     32.6     33.6     33.9     34.6     35.6     36.6 
c9.loc38 c9.loc39 c9.loc40 c9.loc41 c9.loc42 c9.loc43 c9.loc44  D9Mit10 
    37.6     38.6     39.6     40.6     41.6     42.6     43.6     43.7 
c9.loc45 c9.loc46 c9.loc47 c9.loc48 c9.loc49 c9.loc50 c9.loc51 c9.loc52 
    44.6     45.6     46.6     47.6     48.6     49.6     50.6     51.6 
c9.loc53 D9Mit182 c9.loc55 c9.loc56 c9.loc57 c9.loc58 c9.loc59 c9.loc60 
    52.6     53.6     54.6     55.6     56.6     57.6     58.6     59.6 
c9.loc61  D9Mit17 
    60.6     61.2 

$`10`
 D10Mit61 c10.loc25 c10.loc26 c10.loc27 c10.loc28 c10.loc29 c10.loc30 
     24.0      25.0      26.0      27.0      28.0      29.0      30.0 
c10.loc31 c10.loc32 c10.loc33 c10.loc34 c10.loc35 c10.loc36 c10.loc37 
     31.0      32.0      33.0      34.0      35.0      36.0      37.0 
c10.loc38 c10.loc39 c10.loc40 c10.loc41 c10.loc42 c10.loc43 c10.loc44 
     38.0      39.0      40.0      41.0      42.0      43.0      44.0 
c10.loc45 c10.loc46 c10.loc47 c10.loc48 c10.loc49 c10.loc50 c10.loc51 
     45.0      46.0      47.0      48.0      49.0      50.0      51.0 
c10.loc52 c10.loc53 c10.loc54 c10.loc55 c10.loc56 c10.loc57  D10Mit70 
     52.0      53.0      54.0      55.0      56.0      57.0      57.9 

$`11`
 D11Mit20  c11.loc1  c11.loc2  c11.loc3  c11.loc4  c11.loc5  c11.loc6 
      0.0       1.0       2.0       3.0       4.0       5.0       6.0 
 c11.loc7  c11.loc8  c11.loc9 c11.loc10 c11.loc11 c11.loc12 c11.loc13 
      7.0       8.0       9.0      10.0      11.0      12.0      13.0 
c11.loc14 c11.loc15 c11.loc16   D11Mit4 c11.loc17 c11.loc18 c11.loc19 
     14.0      15.0      16.0      16.4      17.0      18.0      19.0 
c11.loc20 c11.loc21 c11.loc22 c11.loc23 c11.loc24 c11.loc25 c11.loc26 
     20.0      21.0      22.0      23.0      24.0      25.0      26.0 
c11.loc27 c11.loc28  D11Mit36 c11.loc29 c11.loc30 c11.loc31 c11.loc32 
     27.0      28.0      28.2      29.0      30.0      31.0      32.0 
 D11Mit41 c11.loc33 c11.loc34 c11.loc35 c11.loc36 D11Mit288 c11.loc37 
     32.3      33.0      34.0      35.0      36.0      36.2      37.0 
c11.loc38 c11.loc39 c11.loc40 c11.loc41 c11.loc42   D8Mit18 c11.loc43 
     38.0      39.0      40.0      41.0      42.0      42.2      43.0 
c11.loc44 c11.loc45 c11.loc46 c11.loc47 c11.loc48 c11.loc49 c11.loc50 
     44.0      45.0      46.0      47.0      48.0      49.0      50.0 
c11.loc51 c11.loc52 c11.loc53 c11.loc54 c11.loc55 c11.loc56 D11Mit101 
     51.0      52.0      53.0      54.0      55.0      56.0      56.9 

$`12`
 D12Mit88 c12.loc21 c12.loc22 c12.loc23 c12.loc24 c12.loc25 c12.loc26 
     19.7      20.7      21.7      22.7      23.7      24.7      25.7 
c12.loc27 c12.loc28 c12.loc29 c12.loc30 c12.loc31 c12.loc32 c12.loc33 
     26.7      27.7      28.7      29.7      30.7      31.7      32.7 
c12.loc34 c12.loc35 c12.loc36 c12.loc37 c12.loc38 c12.loc39 c12.loc40 
     33.7      34.7      35.7      36.7      37.7      38.7      39.7 
c12.loc41 c12.loc42 c12.loc43 c12.loc44 c12.loc45 c12.loc46 c12.loc47 
     40.7      41.7      42.7      43.7      44.7      45.7      46.7 
c12.loc48 c12.loc49 c12.loc50 c12.loc51 c12.loc52 c12.loc53 c12.loc54 
     47.7      48.7      49.7      50.7      51.7      52.7      53.7 
c12.loc55 c12.loc56 c12.loc57 c12.loc58 D12Mit134 
     54.7      55.7      56.7      57.7      57.9 

$`13`
   D13Mit10 c13.loc18.5 c13.loc19.5 c13.loc20.5 c13.loc21.5 c13.loc22.5 
       17.5        18.5        19.5        20.5        21.5        22.5 
c13.loc23.5 c13.loc24.5 c13.loc25.5 c13.loc26.5 c13.loc27.5 c13.loc28.5 
       23.5        24.5        25.5        26.5        27.5        28.5 
c13.loc29.5 c13.loc30.5 c13.loc31.5 c13.loc32.5 c13.loc33.5 c13.loc34.5 
       29.5        30.5        31.5        32.5        33.5        34.5 
c13.loc35.5 c13.loc36.5 c13.loc37.5 c13.loc38.5 c13.loc39.5    D13Mit51 
       35.5        36.5        37.5        38.5        39.5        40.4 

$`14`
 D14Mit54 c14.loc21 c14.loc22 c14.loc23 c14.loc24 c14.loc25 c14.loc26 
     19.7      20.7      21.7      22.7      23.7      24.7      25.7 
c14.loc27 c14.loc28 c14.loc29 c14.loc30 c14.loc31 c14.loc32 c14.loc33 
     26.7      27.7      28.7      29.7      30.7      31.7      32.7 
c14.loc34 c14.loc35 c14.loc36 c14.loc37 c14.loc38 c14.loc39 c14.loc40 
     33.7      34.7      35.7      36.7      37.7      38.7      39.7 
c14.loc41 c14.loc42 c14.loc43 c14.loc44 c14.loc45 c14.loc46 c14.loc47 
     40.7      41.7      42.7      43.7      44.7      45.7      46.7 
c14.loc48 c14.loc49 c14.loc50 c14.loc51 c14.loc52 D14Mit195 
     47.7      48.7      49.7      50.7      51.7      52.5 

$`15`
 D15Mit22 c15.loc17 c15.loc18 c15.loc19 c15.loc20 c15.loc21 c15.loc22 
     16.4      17.4      18.4      19.4      20.4      21.4      22.4 
c15.loc23 c15.loc24 c15.loc25 c15.loc26 c15.loc27 c15.loc28 c15.loc29 
     23.4      24.4      25.4      26.4      27.4      28.4      29.4 
c15.loc30 c15.loc31 c15.loc32 c15.loc33 c15.loc34 c15.loc35 c15.loc36 
     30.4      31.4      32.4      33.4      34.4      35.4      36.4 
c15.loc37 c15.loc38 c15.loc39 c15.loc40 c15.loc41 c15.loc42 c15.loc43 
     37.4      38.4      39.4      40.4      41.4      42.4      43.4 
c15.loc44 c15.loc45 c15.loc46 c15.loc47 c15.loc48 D15Mit159 
     44.4      45.4      46.4      47.4      48.4      49.2 

$`16`
D16Mit131  c16.loc8  c16.loc9 c16.loc10 c16.loc11 c16.loc12 c16.loc13 
      6.6       7.6       8.6       9.6      10.6      11.6      12.6 
c16.loc14 c16.loc15 c16.loc16 c16.loc17 c16.loc18 c16.loc19 c16.loc20 
     13.6      14.6      15.6      16.6      17.6      18.6      19.6 
c16.loc21 c16.loc22 c16.loc23 c16.loc24 c16.loc25   D16Mit4 c16.loc26 
     20.6      21.6      22.6      23.6      24.6      25.1      25.6 
c16.loc27 c16.loc28 c16.loc29 c16.loc30  D16Mit30 c16.loc32 c16.loc33 
     26.6      27.6      28.6      29.6      30.6      31.6      32.6 
c16.loc34 c16.loc35 c16.loc36 c16.loc37 c16.loc38 c16.loc39 c16.loc40 
     33.6      34.6      35.6      36.6      37.6      38.6      39.6 
 D16Mit19 c16.loc41 c16.loc42 c16.loc43 c16.loc44 c16.loc45 c16.loc46 
     40.4      40.6      41.6      42.6      43.6      44.6      45.6 
c16.loc47 c16.loc48 c16.loc49 c16.loc50 c16.loc51  D16Mit70 
     46.6      47.6      48.6      49.6      50.6      51.4 

$`17`
 D17Mit46  c17.loc4  c17.loc5  c17.loc6  c17.loc7  c17.loc8  c17.loc9 
      3.3       4.3       5.3       6.3       7.3       8.3       9.3 
c17.loc10 c17.loc11 c17.loc12 c17.loc13 c17.loc14 c17.loc15 c17.loc16 
     10.3      11.3      12.3      13.3      14.3      15.3      16.3 
c17.loc17 c17.loc18 c17.loc19 c17.loc20 c17.loc21 c17.loc22 c17.loc23 
     17.3      18.3      19.3      20.3      21.3      22.3      23.3 
c17.loc24 c17.loc25 c17.loc26 c17.loc27 c17.loc28 c17.loc29 c17.loc30 
     24.3      25.3      26.3      27.3      28.3      29.3      30.3 
c17.loc31 c17.loc32 c17.loc33 c17.loc34 c17.loc35 c17.loc36 c17.loc37 
     31.3      32.3      33.3      34.3      35.3      36.3      37.3 
c17.loc38  D17Mit93 
     38.3      39.3 

$`18`
 D18Mit20  c18.loc5  c18.loc6  c18.loc7  c18.loc8  c18.loc9 c18.loc10 
      4.4       5.4       6.4       7.4       8.4       9.4      10.4 
c18.loc11 c18.loc12 c18.loc13 c18.loc14 c18.loc15 c18.loc16 c18.loc17 
     11.4      12.4      13.4      14.4      15.4      16.4      17.4 
c18.loc18 c18.loc19 c18.loc20 c18.loc21 c18.loc22 c18.loc23 c18.loc24 
     18.4      19.4      20.4      21.4      22.4      23.4      24.4 
c18.loc25 c18.loc26 c18.loc27 c18.loc28 c18.loc29 c18.loc30 D18Mit186 
     25.4      26.4      27.4      28.4      29.4      30.4      30.6 

$`19`
 D19Mit68  c19.loc4  c19.loc5  c19.loc6  c19.loc7  c19.loc8  c19.loc9 
      3.3       4.3       5.3       6.3       7.3       8.3       9.3 
c19.loc10 c19.loc11 c19.loc12 c19.loc13 c19.loc14 c19.loc15 c19.loc16 
     10.3      11.3      12.3      13.3      14.3      15.3      16.3 
c19.loc17 c19.loc18 c19.loc19 c19.loc20 c19.loc21 c19.loc22 c19.loc23 
     17.3      18.3      19.3      20.3      21.3      22.3      23.3 
c19.loc24 c19.loc25 c19.loc26 c19.loc27 c19.loc28 c19.loc29 c19.loc30 
     24.3      25.3      26.3      27.3      28.3      29.3      30.3 
c19.loc31 c19.loc32 c19.loc33 c19.loc34 c19.loc35 c19.loc36 c19.loc37 
     31.3      32.3      33.3      34.3      35.3      36.3      37.3 
 D19Mit37 
     38.3 

$X
   DXMit16 cX.loc30.5 cX.loc31.5 cX.loc32.5 cX.loc33.5 cX.loc34.5 
      29.5       30.5       31.5       32.5       33.5       34.5 
cX.loc35.5 cX.loc36.5 cX.loc37.5 cX.loc38.5 cX.loc39.5 cX.loc40.5 
      35.5       36.5       37.5       38.5       39.5       40.5 
cX.loc41.5 cX.loc42.5 cX.loc43.5 cX.loc44.5 cX.loc45.5 cX.loc46.5 
      41.5       42.5       43.5       44.5       45.5       46.5 
cX.loc47.5 cX.loc48.5 cX.loc49.5 cX.loc50.5 cX.loc51.5 cX.loc52.5 
      47.5       48.5       49.5       50.5       51.5       52.5 
cX.loc53.5 cX.loc54.5 cX.loc55.5 cX.loc56.5 cX.loc57.5   DXMit186 
      53.5       54.5       55.5       56.5       57.5       57.9 

attr(,"is_x_chr")
    1     2     3     4     5     6     7     8     9    10    11    12 
FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE 
   13    14    15    16    17    18    19     X 
FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE 
~~~
{: .output}

Notice that pseudomarkers are now spaced at 1 cM intervals from genotyped markers. The argument `step=1` generated pseudomarkers at these intervals. 

Next we use `calc_genoprob()` to calculate the QTL genotype probabilities.


~~~
pr <- calc_genoprob(cross=iron, map=map, err=0.002)
~~~
{: .r}

View the first several rows of genotype probabilities for pseudomarker c1.loc28 on chromosome 1.


~~~
head((pr$`1`)[,,"c1.loc28"])
~~~
{: .r}



~~~
           SS         SB          BB
1 0.001066497 0.04841359 0.950519910
2 0.250000000 0.50000000 0.250000000
3 0.001066497 0.04841359 0.950519910
4 0.004495825 0.99099565 0.004508526
5 0.250000000 0.50000000 0.250000000
6 0.250000000 0.50000000 0.250000000
~~~
{: .output}

**A bit more advanced (optional)** To speed up the calculations with large datasets on a multi-core machine, you can use the argument `cores`. With `cores=0`, the number of available cores will be detected via `parallel::detectCores()`. Otherwise, specify the number of cores as a positive integer.


~~~
pr <- calc_genoprob(cross=iron, map=map, err=0.002, cores=4)
~~~
{: .r}

**(optional)** The genome scan functions use genotype probabilities as well as a matrix of phenotypes. If you wished to perform a genome scan via an additive allele model, you would first convert the genotype probabilities to allele probabilities, using the function `genoprob_to_alleleprob()`.


~~~
apr <- genoprob_to_alleleprob(probs=pr)
~~~
{: .r}


> ## Challenge 1
> Explain why calculating genotype probabilities is the first step in QTL analysis.
>
> > ## Solution to Challenge 1
> >
> {: .solution}
{: .challenge}


> ## Challenge 2
> Load the grav2.zip file into an object called `grav`. 
> Insert pseudomarkers and calculate genotype probabilities.
>
> > ## Solution to Challenge 2
> >
> > grav <- read_cross2( system.file("extdata", "grav2.zip", package="qtl2geno") )
> > summary(grav)
> > grav$gmap
> > gravmap <- insert_pseudomarkers(map = grav$gmap, step=2)
> > gravpr <- calc_genoprob(grav, map = gravmap)
> > head((gravpr$`5`)[,,"c5.loc4"])
> {: .solution}
{: .challenge}



