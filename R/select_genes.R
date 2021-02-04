#' Function to perform gene selection
#'
#' This function performs gene selection on a single cell expression dataset. Given a dataset containing cell cluster ids, this function performs differential expression analysis between all pairs of clusters. The final gene set is chosen as the union of top differentially expressed genes obtained from all comparisons between pairs of clusters.
#' @param sce A SingleCellExpression object. The slot "assays" must contain one of the followings: "count" (count data), "normalized" (normalized data), "expression" (neither count data nor normalized data). The slot "colData" must contain a variable "cluster" with cluster ids of the cells. If no cluster id is provided, the program divides the dataset into two clusters and performs the analysis.
#' @param type Type of expression data to be analyzed. It should assume one of these values: "counts", "normalized", "expression". Default is "counts".
#' @param numgenes Number of top genes included in the final dataset from each of the differential expression analyses between pairs of clusters. Default value is 500.
#' @keywords
#' @export
#' @return A SingleCellExperiment object with expression matrix containing selected genes.
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' counts <- matrix(rpois(3000000, lambda = 10), ncol=300, nrow=10000)
#' rownames(counts)<-1:(dim(counts)[1])
#' colnames(counts)<-1:(dim(counts)[2])
#' sce <- SingleCellExperiment(list(counts=counts))
#' km<-kmeans(t(counts),3)
#' colData(sce)$cluster<-km$cluster
#' sce1<-select_genes(sce)
#' }

select_genes<-function(sce,type=c("counts","nomralized","expression"),numgenes=500,seed=12345)
{
set.seed(seed)
type<-type[1]
if(type=="counts")
{
  mat<-assays(sce)$counts
}else if(type=="normalized")
{
  mat<-assays(sce)$normalized
}else if(type=="expression")
{
  mat<-assays(sce)$expression
}else
{
  stop("The \"normalization\" must be one of the followwings: \"TMM\" / \"quant\" / \"cpm\"")
}

if(is.null(colData(sce)$cluster))
{
  km<-kmeans(t(mat),2)
  colData(sce)$cluster<-km$cluster
}


clust<-colData(sce)$cluster
names(clust)<-colnames(sce)
tab<-table(clust)
mat<-mat[(!is.na(rowSums(mat))),]

indd<-length(tab)

library(exactRankTests)
vec<-NULL
ind1<-NULL
ind2<-NULL
k1<-1
for(i in 1:(indd-1))
{
  for(j in (i+1):indd)
  {
    print(c(i,j))
    mat1<-mat[,(clust==i)]
    mat2<-mat[,(clust==j)]
    pval<-NULL
    for(k in 1:(dim(mat1)[1]))
    {
      fit<-wilcox.exact(mat1[k,],mat2[k,])
      pval[k]<-fit$p.value
    }
    vec<-c(vec,order(pval)[1:numgenes])
    ind1[k1]<-i
    ind2[k1]<-j
    k1<-k1+1
  }
}
vec<-na.omit(vec)
mat_modi<-mat[(unique(vec)),]
  
sce1<-sce  

if(type=="counts")
{
  assays(sce1)$counts<-mat_modi
  colData(sce1)$cluster<-clust
}else if(type=="normalized")
{
  assays(sce1)$normalized<-mat_modi
  colData(sce1)$cluster<-clust
}else if(type=="expression")
{
  assays(sce1)$expression<-mat_modi
  colData(sce1)$cluster<-clust
}

return(sce1)
}


