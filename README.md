- [pseudoga](#pseudoga)
- [Installation and Usage](#usage-)
- [Output](#output)
- [Recommendations](#recommendations)


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


# Analysis of Human embryo development data

The data can be downloaded from <https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/yan.rds> . The dataset contains whole transcriptome profiles of human embryo cells at different developmental stages ranging from oocyte to late-blast stage. More details about the dataset can be found here: <https://www.nature.com/articles/nsmb.2660> . For a collection of similar datasets one can visit <https://hemberg-lab.github.io/scRNA.seq.datasets/> . The dataset is available in the form of a cell by gene expression matrix. The read counts have been estimated from raw reads in an single cell RNA-seq experiment. 

## Gene Selection

All genes may not change with pseudotime. So, gene filtering before pseudotime estimation should improve the accuracy the algorithm. We divide the cells into two clusters and find top 2000 differentally expressed genes between these two clusters. However, any other feature selection method can also be used to at this step. 

```
a<-readRDS("yan.rds")
library(pseudoga)
sce<-SingleCellExperiment(list(expression=assays(a)$logcounts))
sce1<-select_genes(sce,type="expression")
```

## Pseudotime estimation by PseudoGA

Next, we perform pseudotime estimation by PseudoGA. Since this can be considered a small daatset, usual PseudoGA should be used here.

```
sce2<-pseudoga(a1,type="expression",normalization="cpm")

```



