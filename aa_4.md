Análisis de agrupamiento (cluster analysis). <br> Parte 4: Especies indicadoras, especies con preferencia por hábitats
================
JR
15 de noviembre, 2020

``` r
knitr::opts_chunk$set(fig.width=12, fig.height=8)
```

Preámbulo
---------

### Cargar paquetes

``` r
library(indicspecies)
```

    ## Loading required package: permute

``` r
source('biodata/funciones.R')
```

### Cargar datos

``` r
load('biodata/Myrtaceae.Rdata')
mi_fam <- mc_myrtc
grupos_upgma_k2 <- readRDS('grupos_upgma_k2.RDS')
table(grupos_upgma_k2)
```

    ## grupos_upgma_k2
    ##  1  2 
    ## 48  2

``` r
grupos_ward_k4 <- readRDS('grupos_ward_k4.RDS')
table(grupos_ward_k4)
```

    ## grupos_ward_k4
    ##  1  2  3  4 
    ## 20 13  2 15

Análisis de especies indicadoras mediante IndVal
------------------------------------------------

### UPGMA

``` r
iva_upgma_k2 <- multipatt(
  x = mi_fam,
  cluster = grupos_upgma_k2,
  func = 'IndVal.g',
  max.order = 1,
  control = how(nperm = 999))
summary(iva_upgma_k2, indvalcomp = TRUE)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 7
    ##  Selected number of species: 1 
    ##  Number of species associated to 1 group: 1 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 2  #sps.  1 
    ##                         A      B  stat p.value   
    ## Chamguava schippii 0.9724 1.0000 0.986   0.002 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ##          Chamguava schippii       Eugenia coloradoensis 
    ##                         541                         609 
    ##        Eugenia galalonensis           Eugenia nesiotica 
    ##                        1975                         502 
    ##         Eugenia oerstediana           Myrcia gatunensis 
    ##                        1838                          56 
    ## Psidium friedrichsthalianum 
    ##                          58

``` r
(p_upgma_adj <- p.adjust(iva_upgma_k2$sign$p.value))
```

    ## [1] 0.014 1.000 1.000 1.000 1.000 0.792 1.000

``` r
(iva_upgma_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_upgma_k2,
  func = "IndVal.g",
  nboot = 1000))
```

    ## $stat
    ##                                     1         2
    ## Chamguava schippii          0.1149116 0.9861250
    ## Eugenia coloradoensis       0.7994359 0.6007514
    ## Eugenia galalonensis        0.6434162 0.7655165
    ## Eugenia nesiotica           0.6161531 0.7876264
    ## Eugenia oerstediana         0.7160733 0.6980251
    ## Myrcia gatunensis           0.4103119 0.8053873
    ## Psidium friedrichsthalianum 0.5931710 0.3849002
    ## 
    ## $lowerCI
    ##                                      1 2
    ## Chamguava schippii          0.06399835 0
    ## Eugenia coloradoensis       0.72843136 0
    ## Eugenia galalonensis        0.57273096 0
    ## Eugenia nesiotica           0.58169045 0
    ## Eugenia oerstediana         0.67665801 0
    ## Myrcia gatunensis           0.27355060 0
    ## Psidium friedrichsthalianum 0.43803974 0
    ## 
    ## $upperCI
    ##                                     1         2
    ## Chamguava schippii          0.7348469 0.9948789
    ## Eugenia coloradoensis       1.0000000 0.6850903
    ## Eugenia galalonensis        1.0000000 0.8195789
    ## Eugenia nesiotica           1.0000000 0.8130874
    ## Eugenia oerstediana         1.0000000 0.7362489
    ## Myrcia gatunensis           0.7483315 0.8976159
    ## Psidium friedrichsthalianum 0.7693093 0.7258662

Ward

``` r
iva_ward_k4 <- multipatt(
  x = mi_fam,
  cluster = grupos_ward_k4,
  func = 'IndVal.g',
  max.order = 2,
  control = how(nperm = 999))
summary(iva_ward_k4, indvalcomp = TRUE)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: IndVal.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 7
    ##  Selected number of species: 3 
    ##  Number of species associated to 1 group: 1 
    ##  Number of species associated to 2 groups: 2 
    ##  Number of species associated to 3 groups: 0 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 3  #sps.  1 
    ##                         A      B  stat p.value   
    ## Chamguava schippii 0.9248 1.0000 0.962   0.002 **
    ## 
    ##  Group 1+2  #sps.  1 
    ##                           A     B  stat p.value  
    ## Eugenia coloradoensis 0.639 1.000 0.799   0.033 *
    ## 
    ##  Group 3+4  #sps.  1 
    ##                          A      B  stat p.value  
    ## Eugenia oerstediana 0.6515 1.0000 0.807   0.023 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ##          Chamguava schippii       Eugenia coloradoensis 
    ##                         541                         609 
    ##        Eugenia galalonensis           Eugenia nesiotica 
    ##                        1975                         502 
    ##         Eugenia oerstediana           Myrcia gatunensis 
    ##                        1838                          56 
    ## Psidium friedrichsthalianum 
    ##                          58

``` r
(p_ward_adj <- p.adjust(iva_ward_k4$sign$p.value))
```

    ## [1] 0.014 0.165 0.284 0.465 0.138 0.465 0.465

``` r
(iva_ward_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_ward_k4,
  func = "IndVal.g",
  nboot = 1000))
```

    ## $stat
    ##                                     1          2         3         4
    ## Chamguava schippii          0.1306131 0.04362739 0.9616600 0.1376320
    ## Eugenia coloradoensis       0.4909648 0.63082020 0.3925476 0.4548909
    ## Eugenia galalonensis        0.5232772 0.41183028 0.5724933 0.4783598
    ## Eugenia nesiotica           0.4310238 0.50071682 0.5896282 0.4645856
    ## Eugenia oerstediana         0.3844816 0.44796082 0.4838530 0.6460583
    ## Myrcia gatunensis           0.2653687 0.45884405 0.6087975 0.2277912
    ## Psidium friedrichsthalianum 0.2297755 0.45986468 0.2408702 0.4763576
    ## 
    ## $lowerCI
    ##                                      1           2 3          4
    ## Chamguava schippii          0.07249781 0.009576414 0 0.02890765
    ## Eugenia coloradoensis       0.43868859 0.570911217 0 0.40284920
    ## Eugenia galalonensis        0.47082852 0.353497614 0 0.42592909
    ## Eugenia nesiotica           0.38458748 0.433735812 0 0.41370724
    ## Eugenia oerstediana         0.35075060 0.407200995 0 0.60780143
    ## Myrcia gatunensis           0.12182961 0.235341451 0 0.10263385
    ## Psidium friedrichsthalianum 0.09684981 0.239113239 0 0.27745013
    ## 
    ## $upperCI
    ##                                     1         2         3         4
    ## Chamguava schippii          0.5661128 0.2261610 0.9871316 0.5604640
    ## Eugenia coloradoensis       0.5638321 0.7122926 0.4696290 0.5244808
    ## Eugenia galalonensis        0.6607778 0.5178623 0.6442954 0.6090145
    ## Eugenia nesiotica           0.5569527 0.6431658 0.6253853 0.6041680
    ## Eugenia oerstediana         0.4564260 0.5248630 0.5132752 0.7532527
    ## Myrcia gatunensis           0.4479367 0.6891586 0.7519679 0.4155301
    ## Psidium friedrichsthalianum 0.3713326 0.6718252 0.5038250 0.6797699

Análisis de especies con preferencia por hábitat mediante el coeficiente de correlación biserial puntual
--------------------------------------------------------------------------------------------------------

### UPGMA

``` r
phi_upgma_k2 <- multipatt(
  mi_fam,
  grupos_upgma_k2,
  func = "r.g",
  max.order = 1,
  control = how(nperm = 999))
summary(phi_upgma_k2)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: r.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 7
    ##  Selected number of species: 1 
    ##  Number of species associated to 1 group: 1 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 2  #sps.  1 
    ##                    stat p.value   
    ## Chamguava schippii 0.74   0.003 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ##          Chamguava schippii       Eugenia coloradoensis 
    ##                         541                         609 
    ##        Eugenia galalonensis           Eugenia nesiotica 
    ##                        1975                         502 
    ##         Eugenia oerstediana           Myrcia gatunensis 
    ##                        1838                          56 
    ## Psidium friedrichsthalianum 
    ##                          58

``` r
(phi_upgma_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_upgma_k2,
  func = "r.g",
  nboot = 1000))
```

    ## $stat
    ##                                       1           2
    ## Chamguava schippii          -0.74006390  0.74006390
    ## Eugenia coloradoensis        0.45133802 -0.45133802
    ## Eugenia galalonensis        -0.40514421  0.40514421
    ## Eugenia nesiotica           -0.63963816  0.63963816
    ## Eugenia oerstediana          0.05738124 -0.05738124
    ## Myrcia gatunensis           -0.28411350  0.28411350
    ## Psidium friedrichsthalianum  0.25934454 -0.25934454
    ## 
    ## $lowerCI
    ##                                      1          2
    ## Chamguava schippii          -0.9995333 -0.3598624
    ## Eugenia coloradoensis        0.1332122 -0.8060062
    ## Eugenia galalonensis        -0.8627533 -0.8670791
    ## Eugenia nesiotica           -0.8041435 -0.8319311
    ## Eugenia oerstediana         -0.2073077 -0.7894229
    ## Myrcia gatunensis           -0.7936029 -0.4562658
    ## Psidium friedrichsthalianum -0.0767718 -0.4961394
    ## 
    ## $upperCI
    ##                                     1           2
    ## Chamguava schippii          0.3575464  0.99953039
    ## Eugenia coloradoensis       0.8051166 -0.14154961
    ## Eugenia galalonensis        0.8669008  0.86236788
    ## Eugenia nesiotica           0.8317054  0.80409400
    ## Eugenia oerstediana         0.7877167  0.20672747
    ## Myrcia gatunensis           0.4551856  0.79078058
    ## Psidium friedrichsthalianum 0.4956841  0.07590529

Ward

``` r
phi_ward_k4 <- multipatt(
  mi_fam,
  grupos_ward_k4,
  func = "r.g",
  max.order = 2,
  control = how(nperm = 999))
summary(phi_ward_k4)
```

    ## 
    ##  Multilevel pattern analysis
    ##  ---------------------------
    ## 
    ##  Association function: r.g
    ##  Significance level (alpha): 0.05
    ## 
    ##  Total number of species: 7
    ##  Selected number of species: 3 
    ##  Number of species associated to 1 group: 3 
    ##  Number of species associated to 2 groups: 0 
    ##  Number of species associated to 3 groups: 0 
    ## 
    ##  List of species associated to each combination: 
    ## 
    ##  Group 2  #sps.  1 
    ##                        stat p.value  
    ## Eugenia coloradoensis 0.568   0.048 *
    ## 
    ##  Group 3  #sps.  1 
    ##                    stat p.value    
    ## Chamguava schippii  0.8   0.001 ***
    ## 
    ##  Group 4  #sps.  1 
    ##                     stat p.value   
    ## Eugenia oerstediana 0.74    0.01 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
colSums(mi_fam)
```

    ##          Chamguava schippii       Eugenia coloradoensis 
    ##                         541                         609 
    ##        Eugenia galalonensis           Eugenia nesiotica 
    ##                        1975                         502 
    ##         Eugenia oerstediana           Myrcia gatunensis 
    ##                        1838                          56 
    ## Psidium friedrichsthalianum 
    ##                          58

``` r
(phi_ward_boot <- strassoc(
  X = mi_fam,
  cluster = grupos_ward_k4,
  func = "r.g",
  nboot = 1000))
```

    ## $stat
    ##                                       1            2          3          4
    ## Chamguava schippii          -0.26279504 -0.289181756  0.8003513 -0.2483745
    ## Eugenia coloradoensis       -0.03440221  0.568409120 -0.3685022 -0.1655047
    ## Eugenia galalonensis         0.11926988 -0.402568107  0.3893127 -0.1060145
    ## Eugenia nesiotica           -0.32086136  0.003584071  0.4879556 -0.1706783
    ## Eugenia oerstediana         -0.45139461 -0.217940113 -0.0701841  0.7395188
    ## Myrcia gatunensis           -0.12307329  0.153322624  0.2007726 -0.2310219
    ## Psidium friedrichsthalianum -0.15181823  0.218559027 -0.2051198  0.1383790
    ## 
    ## $lowerCI
    ##                                      1          2          3          4
    ## Chamguava schippii          -0.3484193 -0.4201525 -0.2417644 -0.3520465
    ## Eugenia coloradoensis       -0.2537072  0.3421685 -0.7021498 -0.3632016
    ## Eugenia galalonensis        -0.1464016 -0.6052462 -0.7606306 -0.3576504
    ## Eugenia nesiotica           -0.4962857 -0.3190691 -0.7330294 -0.3800154
    ## Eugenia oerstediana         -0.5660056 -0.3486656 -0.6841071  0.6018420
    ## Myrcia gatunensis           -0.3318752 -0.2620695 -0.3469129 -0.3962477
    ## Psidium friedrichsthalianum -0.3073850 -0.1098463 -0.3937399 -0.1418937
    ## 
    ## $upperCI
    ##                                       1           2           3
    ## Chamguava schippii           0.28736023 -0.11007549  0.99934101
    ## Eugenia coloradoensis        0.21093875  0.75761770 -0.12587235
    ## Eugenia galalonensis         0.53898671  0.07984603  0.77081277
    ## Eugenia nesiotica            0.15101776  0.46266715  0.67063020
    ## Eugenia oerstediana         -0.11121535  0.08514380  0.05271474
    ## Myrcia gatunensis            0.22949779  0.46450917  0.70341461
    ## Psidium friedrichsthalianum  0.02926989  0.54249618  0.01481957
    ##                                        4
    ## Chamguava schippii           0.355155615
    ## Eugenia coloradoensis        0.076118360
    ## Eugenia galalonensis         0.345110542
    ## Eugenia nesiotica            0.310075036
    ## Eugenia oerstediana          0.855865927
    ## Myrcia gatunensis           -0.002554268
    ## Psidium friedrichsthalianum  0.429332136
