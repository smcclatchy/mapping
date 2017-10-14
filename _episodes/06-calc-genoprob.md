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

![](../fig/threeD_array.png)


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

Now look at a summary of the cross data and the names of each object within the data.


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

Have a look at the markers listed in the genetic map, `gmap`. Markers are listed by chromosome and described by cM position. View only the markers on the first several chromosomes.


~~~
head(iron$gmap)
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
~~~
{: .output}

We then use `insert_pseudomarkers()` to insert pseudomarkers into the
genetic map, which we grab from the `iron` object as `iron$gmap`:


~~~
map <- insert_pseudomarkers(map=iron$gmap, step=1)
~~~
{: .r}

Now have a look at the new object called `map`. View only the first several chromosomes.


~~~
head(map)
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
~~~
{: .output}

Notice that pseudomarkers are now spaced at 1 cM intervals from genotyped markers. The argument `step=1` generated pseudomarkers at these intervals. 

Next we use `calc_genoprob()` to calculate the QTL genotype probabilities.


~~~
pr <- calc_genoprob(cross=iron, map=map, err=0.002)
~~~
{: .r}

Recall that the result of `calc_genoprob`, `pr`, is a list of three-dimensional arrays (one per chromosome). 


~~~
names(pr)
~~~
{: .r}



~~~
 [1] "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
[15] "15" "16" "17" "18" "19" "X" 
~~~
{: .output}

Each 3d array of probabilities is arranged as individuals &times; genotypes &times; positions. Have a look at the names of each of the three dimensions for chromosome 1.


~~~
dimnames(pr$`1`)
~~~
{: .r}



~~~
[[1]]
  [1] "1"   "2"   "3"   "4"   "5"   "6"   "7"   "8"   "9"   "10"  "11" 
 [12] "12"  "13"  "14"  "15"  "16"  "17"  "18"  "19"  "20"  "21"  "22" 
 [23] "23"  "24"  "25"  "26"  "27"  "28"  "29"  "30"  "31"  "32"  "33" 
 [34] "34"  "35"  "36"  "37"  "38"  "39"  "40"  "41"  "42"  "43"  "44" 
 [45] "45"  "46"  "47"  "48"  "49"  "50"  "51"  "52"  "53"  "54"  "55" 
 [56] "56"  "57"  "58"  "59"  "60"  "61"  "62"  "63"  "64"  "65"  "66" 
 [67] "67"  "68"  "69"  "70"  "71"  "72"  "73"  "74"  "75"  "76"  "77" 
 [78] "78"  "79"  "80"  "81"  "82"  "83"  "84"  "85"  "86"  "87"  "88" 
 [89] "89"  "90"  "91"  "92"  "93"  "94"  "95"  "96"  "97"  "98"  "99" 
[100] "100" "101" "102" "103" "104" "105" "106" "107" "108" "109" "110"
[111] "111" "112" "113" "114" "115" "116" "117" "118" "119" "120" "121"
[122] "122" "123" "124" "125" "126" "127" "128" "129" "130" "131" "132"
[133] "133" "134" "135" "136" "137" "138" "139" "140" "141" "142" "143"
[144] "144" "145" "146" "147" "148" "149" "150" "151" "152" "153" "154"
[155] "155" "156" "157" "158" "159" "160" "161" "162" "163" "164" "165"
[166] "166" "167" "168" "169" "170" "171" "172" "173" "174" "175" "176"
[177] "177" "178" "179" "180" "181" "182" "183" "184" "185" "186" "187"
[188] "188" "189" "190" "191" "192" "193" "194" "195" "196" "197" "198"
[199] "199" "200" "201" "202" "203" "204" "205" "206" "207" "208" "209"
[210] "210" "211" "212" "213" "214" "215" "216" "217" "218" "219" "220"
[221] "221" "222" "223" "224" "225" "226" "227" "228" "229" "230" "231"
[232] "232" "233" "234" "235" "236" "237" "238" "239" "240" "241" "242"
[243] "243" "244" "245" "246" "247" "248" "249" "250" "251" "252" "253"
[254] "254" "255" "256" "257" "258" "259" "260" "261" "262" "263" "264"
[265] "265" "266" "267" "268" "269" "270" "271" "272" "273" "274" "275"
[276] "276" "277" "278" "279" "280" "281" "282" "283" "284"

[[2]]
[1] "SS" "SB" "BB"

[[3]]
 [1] "D1Mit18"   "c1.loc28"  "c1.loc29"  "c1.loc30"  "c1.loc31" 
 [6] "c1.loc32"  "c1.loc33"  "c1.loc34"  "c1.loc35"  "c1.loc36" 
[11] "c1.loc37"  "c1.loc38"  "c1.loc39"  "c1.loc40"  "c1.loc41" 
[16] "c1.loc42"  "c1.loc43"  "c1.loc44"  "c1.loc45"  "c1.loc46" 
[21] "c1.loc47"  "c1.loc48"  "c1.loc49"  "c1.loc50"  "c1.loc51" 
[26] "D1Mit80"   "c1.loc52"  "c1.loc53"  "c1.loc54"  "c1.loc55" 
[31] "c1.loc56"  "c1.loc57"  "c1.loc58"  "c1.loc59"  "c1.loc60" 
[36] "c1.loc61"  "c1.loc62"  "c1.loc63"  "c1.loc64"  "c1.loc65" 
[41] "c1.loc66"  "c1.loc67"  "c1.loc68"  "c1.loc69"  "c1.loc70" 
[46] "c1.loc71"  "c1.loc72"  "c1.loc73"  "c1.loc74"  "c1.loc75" 
[51] "c1.loc76"  "c1.loc77"  "c1.loc78"  "c1.loc79"  "c1.loc80" 
[56] "c1.loc81"  "c1.loc82"  "c1.loc83"  "c1.loc84"  "c1.loc85" 
[61] "c1.loc86"  "c1.loc87"  "c1.loc88"  "c1.loc89"  "c1.loc90" 
[66] "c1.loc91"  "c1.loc92"  "c1.loc93"  "c1.loc94"  "c1.loc95" 
[71] "c1.loc96"  "c1.loc97"  "c1.loc98"  "c1.loc99"  "c1.loc100"
[76] "c1.loc101" "c1.loc102" "c1.loc103" "c1.loc104" "c1.loc105"
[81] "c1.loc106" "c1.loc107" "c1.loc108" "c1.loc109" "c1.loc110"
[86] "D1Mit17"  
~~~
{: .output}

View the first several rows of genotype probabilities for the first genotyped marker on chromosome 1, and the two adjacent pseudomarkers located at 1 cM intervals away. Compare the probabilities for each pseudomarker genotype with those of the genotyped marker.


~~~
head((pr$`1`)[,,"D1Mit18"])
~~~
{: .r}



~~~
            SS          SB          BB
1 0.0009965527 0.004443869 0.994559578
2 0.2500000000 0.500000000 0.250000000
3 0.0009965527 0.004443869 0.994559578
4 0.0002244510 0.999550492 0.000225057
5 0.2500000000 0.500000000 0.250000000
6 0.2500000000 0.500000000 0.250000000
~~~
{: .output}



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



~~~
head((pr$`1`)[,,"c1.loc29"])
~~~
{: .r}



~~~
           SS         SB          BB
1 0.001307861 0.09166068 0.907031456
2 0.250000000 0.50000000 0.250000000
3 0.001307861 0.09166068 0.907031456
4 0.008374347 0.98322650 0.008399148
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
> View a summary of the `grav` data.
> View the genetic map for the `grav` data.
> Insert pseudomarkers and calculate genotype probabilities.
> View the first several rows of genotype probabilities for 
> any marker on chromosome 5.
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



