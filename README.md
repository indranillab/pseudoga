- [pseudoga](#pseudoga)
- [Installation and Usage](#usage-)
- [Output](#output)
- [Recommendations](#recommendations)
- [Analysis of Human midbrain cell development data](#analysis-of-human-midbrain-cell-development-data)
  * [Gene Selection](#gene-selection)
  * [Pseudotime estimation by PseudoGA](#pseudotime-estimation-by-pseudoga)
  * [Visualization of estimated pseudotime with respect to actual pseudotime](#visualization-of-estimated-pseudotime-with-respect-to-actual-pseudotime)
  * [Genes with highest correlation with the pseudotime](#genes-with-highest-correlation-with-the-pseudotime)

# pseudoga
Cell pseudotime reconstruction based on genetic algorithm

The package pseudoga can be used to perform pseudotime analysis on single
cell gene expression data. Given a homogeneous population of cells, the cells can be ordered
to form a trajectory. Given a heterogeneous population and cell cluster ids, the packages
can be used to find a tree structure based on pseudotime ordering of cells.

Input must be provided as SingleCellExperiment object with the expression matrix denoting rows as genes 
and columns as cells.


# Installation and Usage
The package can be installed from github using the following command:
```
library(devtools)
install_github("indranillab/pseudoga")
```
The following commands illustarte a simple example of how PseudoGA should be used.
```

library(pseudoga)

library(SingleCellExperiment)

counts <- matrix(rpois(10000, lambda = 10), ncol=100, nrow=100) 

sce <- SingleCellExperiment(list(counts=counts))

sce<-pseudoga(sce) #Usual PseudoGA

sce1<-pseudoga_parallel(sce) #PseudoGA based on subsampling 
```

# Output 

The object "Pseudotime" under "colData" contains inferred pseduotime by PseudoGA.

# Recommendations

For large number of cells, "pseudoga_parallel" is more suitable. One should check all the parameters carefully before applying these two functions.
For details about the parameters, type:
```
?pseudoga

?pseudoga_parallel
```
