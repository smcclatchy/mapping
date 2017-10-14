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
package calculates QTL genotype probabilities, conditional on the available marker data. These are needed for most of the QTL mapping functions. Unlike the corresponding function in [R/qtl](https://rqtl.org), `calc.genoprob()`, the result is not inserted back into the input cross object, but is returned as a list of three-dimensional arrays (one per chromosome). Each 3d array of probabilities is arranged as individuals &times; genotypes &times; positions.

![](../fig/threeD_array.png)


If we wish to perform QTL calculations at positions between markers (so called "pseudomarkers"), we first need to insert such positions into the genetic map with the function `insert_pseudomarkers()`. Unlike [R/qtl](https://rqtl.org), the map is kept separate from the genotype
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

Now look at a summary of the cross data and the names of each variable within the data.


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

Now have a look at the new object called `map`. View only the first two chromosomes.


~~~
head(map, n=2)
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
head((pr$`1`)[,,"D1Mit18"]) # genotyped marker
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
head((pr$`1`)[,,"c1.loc28"]) # pseudomarker 1 cM away
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
head((pr$`1`)[,,"c1.loc29"]) # the next pseudomarker
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



