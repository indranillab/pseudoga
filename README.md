- [pseudoga](#pseudoga)
- [Usage:](#usage-)
- [Output](#output)
- [Recommendations](#recommendations)
- [Analysis of Human midbrain cell development data](#analysis-of-human-midbrain-cell-development-data)
  * [Gene Selection](#gene-selection)
  * [Pseudotime estimation by PseudoGA](#pseudotime-estimation-by-pseudoga)
  * [Visualization of estimated pseudotime with respect to actual pseudotime](#visualization-of-estimated-pseudotime-with-respect-to-actual-pseudotime)
- [Genes with highest correlation with the pseudotime](#genes-with-highest-correlation-with-the-pseudotime)

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

The data can be downloaded from <https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/manno_human.rds> . The dataset contains whole transcriptome profiles of human midbrain cells at different developmental stages with age ranging from 0 day to 77 days. More details about the dataset can be found here: http://dx.doi.org/10.1016/j.cell.2016.09.027 . For a collection of similar datasets one can visit <https://hemberg-lab.github.io/scRNA.seq.datasets/> . The dataset is available in the form of a cell by gene expression matrix. The read counts have been estimated from raw reads in an single cell RNA-seq experiment. 

## Gene Selection

All genes may not change with pseudotime. So, gene filtering before pseudotime estimation should improve the accuracy the algorithm. We divide the cells into two clusters and find top 2000 differentally expressed genes between these two clusters. However, any other feature selection method can also be used to at this step. 

```
a<-readRDS("manno_human.rds")
library(pseudoga)
sce<-SingleCellExperiment(list(expression=assays(a)$logcounts))
sce1<-select_genes(sce,numgenes=2000,type="expression")
```

## Pseudotime estimation by PseudoGA

Next, we perform pseudotime estimation by PseudoGA. Since this can be considered a large daatset, PseudoGA based on subsampling should be used here.

```
sce2<-pseudoga_parallel(sce1,type="expression",normalization="TMM",repl=20,subsample=300)

```

## Visualization of estimated pseudotime with respect to actual pseudotime

The actual age of the cells are available. Next we compare them with the estimated pseudotime:

```
time<-colData(a)$age
time1<-rep(NA,length(time))

time1[(time=="day_0")]<-1
time1[(time=="day_12")]<-2
time1[(time=="day_17")]<-3
time1[(time=="day_35")]<-4
time1[(time=="day_42")]<-5
time1[(time=="day_63")]<-8
time1[(time=="week_6")]<-5
time1[(time=="week_7")]<-6
time1[(time=="week_8")]<-7
time1[(time=="week_9")]<-8
time1[(time=="week_10")]<-9
time1[(time=="week_11")]<-10
time<-time1

library(viridis)
cols<-viridis(12)[time]
ndays<-c(0,12,17,35,42,49,56,63,70,77)
days<-c("DAY 0","DAY 12","DAY 17","DAY 35","DAY 42","DAY 49","DAY 56","DAY 63","DAY 70","DAY 77")

oripath1<-order(colData(sce2)$Pseudotime)
par(mar=c(5,6,4,2)+0.1)
plot(ndays[time][oripath1],col=cols[oripath1],ylab="Days of development",xlab="Pseudotime",main="Human brain development data",yaxt="n",cex.lab=2,cex.main=2)
axis(side=2, at=ndays[c(1,3,5,7,9)], labels =days[c(1,3,5,7,9)])

```
![](https://github.com/pronoymondal/pseudogadata/blob/main/manno_path1_png.png)

# Genes with highest correlation with the pseudotime
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
plot(rank(colData(sce2)$Pseudotime),data[geneord1[1],],col="red",main="Expression of LIN28A",xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)

par(mar = c(5, 5, 4, 2))
plot(rank(colData(sce2)$Pseudotime),data[geneord1[2],],col="red",main="Expression of RPL23A",xlab="Pseudotime",ylab="Expression",cex.lab=1.8,font.lab=2,cex.main=2,cex.axis=2)

```
A | B
- | - 
![alt](https://github.com/pronoymondal/pseudogadata/blob/main/manno_cluster1_gene1%20.png) | ![alt](https://github.com/pronoymondal/pseudogadata/blob/main/manno_cluster1_gene2.png)
