# TIGRESS: Trustful Inference of Gene REgulation using Stability Selection

This R package implements the TIGRESS method for gene regulatory network (GRN) 
inference from gene expression data.

TIGRESS formulates GRN inference as a sparse regression problem, and solves it 
with a LARS regression model combined with stability selection. Check reference 
below for more details.

## Installation
To install TIGRESS from R, type:

```{r}
install.packages('devtools')
library(devtools)
install_github("jpvert/tigress")
```

## Citation
A.-C. Haury, F. Mordelet, P. Vera-Licona and J.-P. Vert. [TIGRESS: trustful inference of gene regulation using stability selection](https://doi.org/10.1186/1752-0509-6-145). BMC systems biology 6(1), 145-153, 2012
