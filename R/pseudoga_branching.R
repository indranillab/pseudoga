#' Function to perform branching
#'
#' This function produces a tree structure with cells as nodes.This graph can be used to detect branching.
#' @param sce A SingleCellExpression object. The slot "assays" must contain one of the followings: "count" (count data), "normalized" (normalized data), "expression" (neither count data nor normalized data). The slot "colData" must contain a variable named "cluster" that contains cluster ids of individual cells.
#' @param type Type of expression data to be analyzed. It should assume one of these values: "counts", "normalized", "expression". Default is "counts".
#' @param ntest Size of each generation considered in genetic algorithm is (8*ntest) . Default value for ntest is 50.
#' @param repl Number of independent pseudotime estimates used to calculate final trajectory. Default value is 3.
#' @param subset Number of cells used in each of the independent pseudotime estimation. Default value is 100. For larger datasets, larger values are recommended.
#' @param minit Minimum number of generations considered in genetic algorithm. Default value is 30.
#' @param epsilon The tolerance value used in convergence. Default value is \eqn{10^{-4}}.
#' @param normalization Normalization method applied on the input dataset. It must be one of the following "default" (DESeq Normalization), "TMM" (TMM Normalization), "quant" (Quantile Normalization) or "cpm" (Counts Per Million).
#' @param nnprop Proportion of sample size used as value of k for kNN smoothing. If subsample size is N, (N*nnprop) is considered as value of k. Default value is 0.3.
#' @keywords
#' @export
#' @return An igraph object with cells as nodes. Consecutive cells connected with edges can be considered as a path.
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' counts <- matrix(rpois(30000, lambda = 10), ncol=300, nrow=100)
#' sce <- SingleCellExperiment(list(counts=counts))
#' km<-kmeans(t(counts),3)
#' colData(sce)$cluster<-km$cluster
#' gr<-pseudoga_branching(sce)
#' }


pseudoga_branching<-function(sce,type=c("counts","nomralized","expression"),ntest=50,repl=3,subsample=100,minit=30,epsilon=0.0001,normalization=c("default","TMM","quant","cpm"),nnprop=0.3)
{
  library(SingleCellExperiment)
  type<-type[1]
  normalization<-normalization[1]
  if(type=="counts")
  {
    if(normalization=="default")
    {
      library(DESeq)
      library(Matrix)
      data<-assays(sce)$counts
      data<-data[(rowSums(abs(data))>0),]
      cds<-newCountDataSet(data,conditions=rep(1,dim(data)[2]))
      cds<-estimateSizeFactors(cds)
      data1<-normalise(cds)
    }else if(normalization=="TMM")
    {
      library(edgeR)
      data<-assays(sce)$counts
      y<-DGEList(counts=data)
      y<-calcNormFactors(y)
      data1<-t(apply(data,1,"/",y$samples[,3]))
    }else if(normalization=="quant")
    {
      library(HEM)
      data<-assays(sce)$counts
      data1<-quant.normal(data)
    }else if(normalization=="cpm")
    {
      data<-assays(sce)$counts
      y<-colSums(data)
      data1<-t(apply(data,1,"/",y))
    }else
    {
      stop("The \"normalization\" must be one of the followwings: \"TMM\" / \"quant\" / \"cpm\"")
    }

  }else if(type=="normalized")
  {
    data<-assays(sce)$normalized
    data1<-as.matrix(data)
  }else if(type=="expression")
  {
    if(normalization=="default")
    {
      stop("The \"normalization\" must be one of the followwings: \"TMM\" / \"quant\" / \"cpm\"")
    }else if(normalization=="TMM")
    {
      library(edgeR)
      data<-assays(sce)$expression
      y<-DGEList(counts=data)
      y<-calcNormFactors(y)
      data1<-t(apply(data,1,"/",y$samples[,3]))
    }else if(normalization=="quant")
    {
      library(HEM)
      data<-assays(sce)$expression
      data1<-quant.normal(data)
    }else if(normalization=="cpm")
    {
      data<-assays(sce)$expression
      y<-colSums(data)
      data1<-t(apply(data,1,"/",y))
    }else
    {
      stop("The \"normalization\" must be one of the followwings: \"TMM\" / \"quant\" / \"cpm\"")
    }
  }else
  {
    stop("The \"type\" must be one of the followwings: \"counts\" / \"normalized\" / \"expression\"")
  }


  clust<-colData(sce)$cluster
  names(clust)<-colnames(sce)
  tab<-table(clust)
  if(length(tab)>2)
  {
    for(k in 1:(length(tab)))
    {
      data2<-data1[,(clust==(names(tab)[k]))]
      ord<-findorders_normalized_parallel(data2,ntest,repl,subsample,minit,epsilon,nnprop)
      ord<-order(ord)
      assign(paste0("ord",k),ord)
    }
  }
  if(length(tab)>2)
  {
    for(k in 1:(length(tab)))
    {
      ord<-get(paste0("ord",k))
      path<-order(ord)
      assign(paste0("path",k),path)
      y<-which(clust==k)
      assign(paste0("oripath",k),y[path])
    }
  }


  cluster_id<-rep(NA,(dim(data2)[2]))
  pseudotime<-rep(NA,(dim(data2)[2]))

  for(i in 1:(length(tab)))
  {
    pst<-get(paste0("path",i))
    pos<-get(paste0("oripath",i))
    cluster_id[pos]<-i
    pseudotime[pos]<-(1:length(pos))
  }

  #source("join_cluster.R")
  graph<-join_cluster(cluster_id,pseudotime,t(data1))

  M<-length(tab)
  A<-array(0,dim<-c(2*M,2*M))
  colnames(A)<-1:(2*M)
  rownames(A)<-1:(2*M)

  library(igraph)
  for(i in 1:(M-1))
  {
    for(j in (i+1):M)
    {
      if(graph[i,j]==1)
      {
        A[i,j]<-1
      }
      if(graph[i,j]==2)
      {
        A[(M+i),j]<-1
      }
      if(graph[i,j]==3)
      {
        A[(i),(M+j)]<-1
      }
      if(graph[i,j]==4)
      {
        A[(M+i),(M+j)]<-1
      }
    }
  }
  gr<-graph_from_adjacency_matrix(A)
  gr<-as.undirected(gr)
  V(gr)$name<-paste0("V",V(gr)$name)
  y<-colnames(data1)
  library(RColorBrewer)
  cols<-brewer.pal(M,"Accent")
  gr<-gr+vertices(y,col=cols[clust])

  for(i in 1:M)
  {
    if(!is.null(names(get(paste0("oripath",i)))))
    {
    seq<-c(paste0("V",i),names(get(paste0("oripath",i))),paste0("V",(M+i)))
    }else
    {
      seq<-c(paste0("V",i),get(paste0("oripath",i)),paste0("V",(M+i)))
    }
    gr<-gr+path(seq,col=cols[i])
  }
  V(gr)$color<-c(rep(1,2*M),cols[clust])
  return(gr)

}






