# BLEND

BLEND is a computational tool to estimate cell type fractions from bulk
RNA-seq count data utilizing multiple references. BLEND individualizes
references for each bulk sample and employs a bag-of-words
representation for deconvolution.

**No cell type marker gene selection needed!**

**No data normalization/transformation needed!**

**No reference quality evaluation needed!**

<figure>
  <img src="BLEND logo.pdf" style="width:35.0%" alt="BLEND logo" />
  <figcaption aria-hidden="true">BLEND logo</figcaption>
</figure>

## Installation

``` r
if (!"BLEND" %in% rownames(installed.packages())) {
  install_github('Penghuihuang2000/BLEND')
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

Bulk data must be counts.

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

Two estimation strategies are provided: Gibbs sampling and EM-MAP
algorithm. They provide consistent estimates. EM-MAP algorithm is times
faster than the Gibbs sampler. Because it estimates parameters directly
and is implemented using Rcpp.

``` r
## Run EMMAP for parameter estimation
## Use 30 cores for computation
time.EMMAP <- system.time(res.EMMAP <- BLEND(bulk = BLEND_example$bulk,
             phi = BLEND_example$reference,
             method = "EMMAP",
             ncore = 30))[3]
cat("\n Average running time per sample using one core: ",round((time.EMMAP*30)/(30*60),1), "min")
```

    ## 
    ##  Average running time per sample using one core:  2 min

``` r
## Run Gibbs sampler for parameter estimation
## Use 30 cores for computation
time.GIBBS <- system.time(res.GIBBS <- BLEND(bulk = BLEND_example$bulk,
             phi = BLEND_example$reference,
             method = "GIBBS",
             ncore = 30))[3]
```

    ## No group or design set. Assuming all samples belong to one group.

``` r
cat("\n Average running time per sample using one core: ",round((time.GIBBS*30)/(30*60),1), "min")
```

    ## 
    ##  Average running time per sample using one core:  13.8 min

Here, we explain the results using the first bulk sample as an example.

``` r
## cellular fractions
res.EMMAP[[1]]$`cellular frac`
```

    ##       Astrocytes       Endothelia       Excitatory       Inhibitory 
    ##     0.5081756057     0.0716156096     0.2021741995     0.0001397791 
    ##        Microglia Oligodendrocytes             OPCs 
    ##     0.0302405858     0.1871124131     0.0005418072

``` r
## Reference mixing proportions
res.EMMAP[[1]]$`ref prop`
```

    ## $Astrocytes
    ##           CA           DM           IP           F5           MM           LK 
    ## 4.500419e-03 1.482121e-05 4.847398e-01 2.771210e-01 1.242334e-01 4.035938e-05 
    ##           VL           NG           TS 
    ## 2.071604e-05 1.093285e-01 9.969070e-07 
    ## 
    ## $Endothelia
    ##           CA           DM           IP           MM           LK           VL 
    ## 4.992578e-05 2.126598e-05 7.366538e-01 3.597212e-05 2.630072e-01 1.885855e-04 
    ##           NG           TS 
    ## 3.551273e-05 7.726283e-06 
    ## 
    ## $Excitatory
    ##           CA           LK           VL           NG           TS 
    ## 5.356335e-01 4.641954e-01 3.221351e-05 1.360754e-04 2.817460e-06 
    ## 
    ## $Inhibitory
    ##          CA          LK          VL          NG          TS 
    ## 0.054266383 0.726236897 0.039835533 0.175864093 0.003797094 
    ## 
    ## $Microglia
    ##           CA           DM           IP           MM           LK           VL 
    ## 1.824513e-04 1.162601e-04 5.250361e-05 2.076964e-05 1.975170e-02 9.674254e-01 
    ##           NG           TS 
    ## 1.245092e-02 0.000000e+00 
    ## 
    ## $Oligodendrocytes
    ##           CA           DM           IP           MM           LK           VL 
    ## 8.304891e-02 3.497426e-05 1.970021e-01 1.578281e-01 5.618905e-01 7.749996e-05 
    ##           NG           TS 
    ## 1.149087e-04 2.959614e-06 
    ## 
    ## $OPCs
    ##          CA          DM          LK          VL          NG          TS 
    ## 0.087776815 0.013908403 0.874975991 0.008380233 0.013927592 0.001030967

``` r
## number of updates till convergence
res.EMMAP[[1]]$n.cvg
```

    ## [1] 638

Now, we extract cellular fractions estimated using EM-MAP and Gibbs
sampling. They provide consistent fraction estimates. Thus, in
applications, we recommend using EMMAP option.

``` r
frac.EMMAP <- rlist::list.rbind(lapply(res.EMMAP, function(x){x[[1]]}))
frac.GIBBS <- rlist::list.rbind(lapply(res.GIBBS, function(x){x[[1]]}))
mean(abs(frac.EMMAP - frac.GIBBS))
```

    ## [1] 0.0003759182

Use cell size adjustment to get cell fractions.

``` r
frac.cell <- t(apply(frac.EMMAP, 1, function(x){(x/BLEND_example$cell_size)/
    sum(x/BLEND_example$cell_size)}))
```
