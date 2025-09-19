# BLEND

BLEND is a computational tool to estimate cell type fractions from bulk
RNA-seq count data utilizing multiple references. BLEND individualizes
references for each bulk sample and employs a bag-of-words
representation for deconvolution.

**No cell type marker gene selection needed!**

**No data normalization/transformation needed!**

**No reference quality evaluation needed!**

<figure>
  <img src="BLEND logo.png" style="width:15%;" alt="BLEND logo" />
</figure>

## Installation

``` r
if (!"BLEND" %in% rownames(installed.packages())) {
  devtools::install_github('Penghuihuang2000/BLEND')
}
library(BLEND)
```

## Data input

Here, we provide an example on how to use BLEND. The bulk samples are a
small subset of the MSBB data in the manuscript. The reference list was
collected by Sutton et al. Original bulk sample names were replaced by
pseudonames.

``` r
load("BLEND_example.RData")
names(BLEND_example)
```

    ## [1] "bulk"      "reference" "cell_size"

Bulk data must be counts. Linear scale transformed reference is allowed.

``` r
dim(BLEND_example$bulk) # 6377 genes 30 bulk samples
```

    ## [1] 6377   30

``` r
BLEND_example$bulk[1:6,1:6] 
```

    ##                 sample_1 sample_2 sample_3 sample_4 sample_5 sample_6
    ## ENSG00000162929      626      363      638      992     1015      603
    ## ENSG00000162384      307      212      419      509      541      287
    ## ENSG00000110696     1146      790     1353     2318     2166     1216
    ## ENSG00000164972       55       44       78       80       44       62
    ## ENSG00000182307      517      410      402      528      460      571
    ## ENSG00000149179      873      515      527      984      748      822

References should be provided as a list, each element of which
represents the cell type-specific gene expression matrix of a cell type
from multiple references. Names of the reference list should be cell
type names. Column names of the matrix should be the reference names.
Note that different cell types are allowed to have different number of
references.

``` r
names(BLEND_example$reference)
```

    ## [1] "Astrocytes"       "Endothelia"       "Excitatory"       "Inhibitory"      
    ## [5] "Microglia"        "Oligodendrocytes" "OPCs"

``` r
dim(BLEND_example$reference$Astrocytes)
```

    ## [1] 6377    9

``` r
colnames(BLEND_example$reference$Astrocytes) # 9 references
```

    ## [1] "CA" "DM" "IP" "F5" "MM" "LK" "VL" "NG" "TS"

``` r
colnames(BLEND_example$reference$Excitatory) # 5 references
```

    ## [1] "CA" "LK" "VL" "NG" "TS"

BLEND is an RNA-based model. Thus, the estimated fractions are RNA
molecule fractions of cell types. If the cell fractions are needed, cell
size adjustment is required. Cell size vector quantifies the relative
RNA abundance of different cell types. It can be estimated from
scRNA-seq(sn) data or provided by other sources.

``` r
BLEND_example$cell_size
```

    ##      astro       vasc         ex        inh     immune      oligo        OPC 
    ## 0.10787251 0.05709066 0.36358833 0.19819637 0.04936413 0.07609443 0.14779357

## Deconvolution

Three estimation strategies are provided: mixSQP, Gibbs sampling, and EM-MAP
algorithm. They provide consistent estimates. The mixSQP algorithm is our latest algorithm and is
faster than the EM-MAP algorithm and the Gibbs sampler. Speed: mixSQP > EM-MAP >> Gibbs.

``` r
## Default estimation algorithm (the fastest)
## Run mixSQP for parameter estimation
## Use 5 cores for computation
time.mixSQP <- system.time(res.mixSQP <- BLEND(bulk = BLEND_example$bulk,
             phi = BLEND_example$reference,
             method = "mixSQP",
             ncore = 5))[3]
cat("Average running time per sample using one core: ",round((time.mixSQP*5)/30,1), "sec")
```

    ## Average running time per sample using one core:  0.6 sec


``` r
## Run EMMAP for parameter estimation
## Use 5 cores for computation
time.EMMAP <- system.time(res.EMMAP <- BLEND(bulk = BLEND_example$bulk,
             phi = BLEND_example$reference,
             method = "EMMAP",
             ncore = 5))[3]
cat("Average running time per sample using one core: ",round((time.EMMAP*5)/(30*60),1), "min")
```

    ## Average running time per sample using one core:  0.5 min

``` r
## Run Gibbs sampler for parameter estimation
## Use 5 cores for computation
#time.GIBBS <- system.time(res.GIBBS <- BLEND(bulk = BLEND_example$bulk,
#             phi = BLEND_example$reference,
#             method = "GIBBS",
#             ncore = 5))[3]
#cat("Average running time per sample using one core: ",round((time.GIBBS*30)/(30*60),1), "min")
```


Here, we retrieve estimation results.

``` r
## cellular fractions
round(res.mixSQP$cellular_fraction[1:6,], 3)
```

    ##          Astrocytes Endothelia Excitatory Inhibitory Microglia Oligodendrocytes
    ## sample_1      0.509      0.072      0.203      0.000     0.029            0.187
    ## sample_2      0.541      0.067      0.057      0.000     0.013            0.321
    ## sample_3      0.484      0.033      0.073      0.033     0.032            0.303
    ## sample_4      0.275      0.029      0.541      0.048     0.002            0.091
    ## sample_5      0.250      0.032      0.451      0.028     0.013            0.209
    ## sample_6      0.365      0.062      0.362      0.025     0.026            0.161
    ##           OPCs
    ## sample_1 0.000
    ## sample_2 0.000
    ## sample_3 0.041
    ## sample_4 0.013
    ## sample_5 0.016
    ## sample_6 0.000

``` r
## reference mixing proportions for excitatory neuron
## Interpretation: the underlying excitatory neuron cell type-specific gene expression is likely in between CA and LK references' expression 
round(res.mixSQP$reference_mix_prop$Excitatory[1:6,],3)
```

    ##             CA    LK VL    NG TS
    ## sample_1 0.539 0.461  0 0.000  0
    ## sample_2 0.699 0.301  0 0.000  0
    ## sample_3 1.000 0.000  0 0.000  0
    ## sample_4 0.634 0.260  0 0.106  0
    ## sample_5 0.748 0.252  0 0.000  0
    ## sample_6 0.476 0.421  0 0.103  0
