# pseudoga
Cell pseudotime reconstruction based on genetic algorithm

The package pseudoga can be used to perform pseudotime analysis on single
cell gene expression data. Given a homogeneous population of cells, the cells can be ordered
to form a trajectory. Given a heterogeneous population and cell cluster ids, the packages
can be used to find a tree structure based on pseudotime ordering of cells.

Input must be provided as SingleCellExperiment object with the expression matrix denoting rows as genes 
and columns as cells.


# Usage:
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

# Analysis of Human midbrain cell development data

The data can be downloaded from <https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/manno_human.rds> . For a collection of similar datasets one can visit <https://hemberg-lab.github.io/scRNA.seq.datasets/> .







