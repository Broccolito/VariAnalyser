
# VariAnalyser

<!-- badges: start -->

<!-- badges: end -->

The goal of VariAnalyser is to functionally annotate list of variants from 
Genome-wide Association Studies (GWAS). Particularly, VariAnalyser considers 
CADD functional predictions, Gene Oncology (GO) and KEGG pathways and disease 
associations with genes.

## Installation

You can install the most recent version of VariAnalyser like so:

``` r
if(!require("devtools")){
    install.packages("devtools")
    library("devtools")
}
install_github("Broccolito/VariAnalyser")
```

## Example

To annotate the GWAS results and include pathway and disease information

``` r
library(VariAnalyser)
PathwayDiseaseAnnotator("data/example/FHS_EA_MRS_merged.txt")
```

To make a Manhattan plot

```r
ManhattanPlotter("data/example/FHS_EA_MRS_merged.txt")
```

