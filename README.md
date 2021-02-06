- [pseudoga](#pseudoga)
- [Installation and Usage](#installation-and-usage)
- [Output](#output)
- [Recommendations](#recommendations)
- [Analysis of Human embryo development data](#analysis-of-human-embryo-development-data)
  * [Gene Selection](#gene-selection)
  * [Pseudotime estimation by PseudoGA](#pseudotime-estimation-by-pseudoga)
  * [Comparison with developmental stage](#comparison-with-developmental-stage)
  * [Genes with highest correlation with the pseudotime](#genes-with-highest-correlation-with-the-pseudotime)
- [Analysis of Mouse brain development data](#analysis-of-mouse-brain-development-data)
  * [Gene Selection](#gene-selection-1)
  * [Pseudotime estimation by PseudoGA](#pseudotime-estimation-by-pseudoga-1)
  * [Visualization with principal components](#visualization-with-principal-components)
  * [Genes with highest correlation with the pseudotime](#genes-with-highest-correlation-with-the-pseudotime-1)



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
sce<-readRDS("yan.rds")
library(pseudoga)
assays(sce)$expression<-assays(sce)$logcounts
sce1<-select_genes(sce,type="expression")
```

## Pseudotime estimation by PseudoGA

Next, we perform pseudotime estimation by PseudoGA. Since this can be considered a small daatset, usual PseudoGA should be used here.

```
sce2<-pseudoga(sce1,type="expression",normalization="cpm")

```

## Comparison with developmental stage

Next we compare the estimated pseudotime by PseudoGA with the developmental stages of the cells. 

```
type<-colData(sce2)$cell_type2
type1<-rep(NA,length(type))
type1[(type=="oocyte")]<-1
type1[(type=="zygote")]<-2
type1[(type=="2cell")]<-3
type1[(type=="4cell")]<-4
type1[(type=="8cell")]<-5
type1[(type=="morula")]<-6
type1[(type=="lateblast")]<-7

colData(sce2)$cell_type3<-as.character(type1)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)

p<-ggplot(as.data.frame(colData(sce2)), aes(x = rank(Pseudotime,ties.method="random"), y = cell_type3, 
                                             colour = cell_type3)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic()+
  xlab("PseudoGA pseudotime ordering") + ylab("Developmental Stage") +
  ggtitle("Estimated Pseudotime and developmental stage")

p<-p+scale_y_discrete( labels = unique(colData(sce2)$cell_type2))
p<-p+  scale_color_discrete(name = "Stage", labels = unique(colData(sce2)$cell_type2))

ggsave("comparison_plot1.png")
```

![](https://github.com/pronoymondal/pseudogadata/blob/main/comparison_plot.png)

## Genes with highest correlation with the pseudotime
The following command shows names of top few genes that have highest linear rank correlation with the estimated pseudotime. 

```
cors1<-NULL
data<-assays(sce2)$expression
cors1<-cor(t(data),colData(sce2)$Pseudotime,method="spearman")
abscors1<-abs(cors1)
geneord1<-order(abscors1,decreasing=TRUE)
pseudogene1<-rownames(data)[geneord1]
head(pseudogene1)

```
The following commands plots first two genes with highest correlation with the estimated pseudotime:
  
```
par(mar = c(5, 5, 4, 2))
plot(rank(colData(sce2)$Pseudotime,ties.method="random"),data[geneord1[1],],col="red",main="Expression of ACCSL",xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)

dev.new()
par(mar = c(5, 5, 4, 2))
plot(rank(colData(sce2)$Pseudotime,ties.method="random"),data[geneord1[2],],col="red",main="Expression of C21orf7",xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)

```

![alt](https://github.com/pronoymondal/pseudogadata/blob/main/yan_gene1.png) | ![alt](https://github.com/pronoymondal/pseudogadata/blob/main/yan_gene2.png)


# Analysis of Mouse brain development data

The data can be downloaded from <https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/marques.rds> . This dataset contains transcriptome of 5053 mouse oligodendrocyte cells. Details of the dataset can be found at <https://science.sciencemag.org/content/352/6291/1326>

## Gene Selection
We select top 2000 genes that are differentially expressed between two clusters of cells generated from the datset.

```
sce<-readRDS("marques.rds")
assays(sce)$expression<-assays(sce)$logcounts
sce1<-select_genes(sce,numgenes=2000,type="expression")
```

## Pseudotime estimation by PseudoGA

We apply PseudoGA with subsampling. Parameters like normalization method, subsample size or number of replicates used in this estimation can be altered. The parameters given here are only for illustration.

```
sce4<-pseudoga_parallel(sce2,type="expression",normalization="cpm",subsample=300,repl=20)
```

## Visualization with principal components

First principal component is plotted against the estimated pseudotime showing continuum of different cell types.

```
prc<-prcomp(t(assays(sce4)$expression))

colData(sce4)$prcomp1<-prc$x[,1]
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)

set.seed(0)
p<-ggplot(as.data.frame(colData(sce4)), aes(x = rank(Pseudotime,ties.method="random"), y = prcomp1, 
                                            colour = cell_type1)) +
  geom_quasirandom(groupOnX = FALSE) +
  scale_color_tableau() + theme_classic()+
  xlab("PseudoGA pseudotime ordering") + ylab("Principal Component I") +
  ggtitle("Estimated Pseudotime and first principal component")

p<-p+scale_y_discrete( labels = unique(colData(sce2)$cell_type1))
p<-p+  scale_color_discrete(name = "Stage", labels = unique(colData(sce2)$cell_type1))

ggsave("prcomp_plot1.png")

```

![](https://github.com/pronoymondal/pseudogadata/blob/main/prcomp_plot1.png)

## Genes with highest correlation with the pseudotime

To find the genes with highest correlation with the pseudotime, run the commands below:

```
cors1<-NULL
data<-assays(sce4)$expression
cors1<-cor(t(data),colData(sce4)$Pseudotime,method="spearman")
abscors1<-abs(cors1)
geneord1<-order(abscors1,decreasing=TRUE)
pseudogene1<-rownames(data)[geneord1]
head(pseudogene1)
```
The following commands plots first two genes with highest correlation with the estimated pseudotime:

```
par(mar = c(5, 5, 4, 2))
plot(rank(colData(sce4)$Pseudotime),data[geneord1[1],],col="red",main="Expression of Car2",xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)

par(mar = c(5, 5, 4, 2))
plot(rank(colData(sce4)$Pseudotime),data[geneord1[2],],col="red",main="Expression of Rtn1",xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)
```

![](https://github.com/pronoymondal/pseudogadata/blob/main/marques_gene1.png)
![](https://github.com/pronoymondal/pseudogadata/blob/main/marques_gene2.png)

