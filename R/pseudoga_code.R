set.seed(100)
#' Function to estimate pseudotime
#'
#' This function estimates pseudotime using PseudoGA algorithm. Ordering of cells is performed using genetic algorithm. Subsequently, pseudotimes are updated using nearest neighbor smoothing.
#' @param sce A SingleCellExpression object. The slot "assays" must contain one of the followings: "count" (count data), "normalized" (normalized data), "expression" (neither count data nor normalized data).
#' @param type Type of expression data to be analyzed. It should assume one of these values: "counts", "normalized", "expression". Default is "counts".
#' @param ntest Size of each generation considered in genetic algorithm is (8*ntest) . Default value for ntest is 50.
#' @param minit Minimum number of generations considered in genetic algorithm. Default value is 30.
#' @param epsilon The tolerance value used in convergence. Default value is \eqn{10^{-4}}.
#' @param normalization Normalization method applied on the input dataset. It must be one of the following "default" (DESeq2 Normalization), "TMM" (TMM Normalization), "quant" (Quantile Normalization) or "cpm" (Counts Per Million).
#' @param nnprop Proportion of sample size used as value of k for kNN smoothing. If sample size is N, (N*nnprop) is considered as value of k. Default value is 0.3.
#' @keywords
#' @export
#' @return An updated SingleCellExperiment object in which colData contains values for pseudotime for each cell with slot name "Pseudotime" .
#' @examples
#' \dontrun{
#' library(SingleCellExperiment)
#' counts <- matrix(rpois(10000, lambda = 10), ncol=100, nrow=100)
#' sce <- SingleCellExperiment(list(counts=counts))
#' sce<-pseudoga(sce)
#' }

pseudoga<-function(sce,type=c("counts","nomralized","expression"),ntest=50,minit=30,epsilon=0.0001,normalization=c("default","TMM","quant","cpm"),nnprop=0.3)
{
  library(SingleCellExperiment)
  type<-type[1]
  normalization<-normalization[1]
  if(type=="counts")
  {
    if(normalization=="default")
    {
      data<-assays(sce)$counts
      ord<-findorders_counts(data,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else if(normalization=="TMM")
    {
      library(edgeR)
      data<-assays(sce)$counts
      y<-DGEList(counts=data)
      y<-calcNormFactors(y)
      data1<-t(apply(data,1,"/",y$samples[,3]))
      ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else if(normalization=="quant")
    {
      library(HEM)
      data<-assays(sce)$counts
      data1<-quant.normal(data)
      ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else if(normalization=="cpm")
    {
      data<-assays(sce)$counts
      y<-colSums(data)
      data1<-t(apply(data,1,"/",y))
      ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else
    {
      stop("The \"normalization\" must be one of the followwings: \"TMM\" / \"quant\" / \"cpm\"")
    }

  }else if(type=="normalized")
  {
    data<-assays(sce)$normalized
    data1<-as.matrix(data)
    ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
    colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
    return(sce)
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
      ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else if(normalization=="quant")
    {
      library(HEM)
      data<-assays(sce)$expression
      data1<-quant.normal(data)
      ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else if(normalization=="cpm")
    {
      data<-assays(sce)$expression
      y<-colSums(data)
      data1<-t(apply(data,1,"/",y))
      ord<-findorders_normalized(data1,ntest,minit,epsilon,nnprop)
      colData(sce)$Pseudotime<-(ord-min(ord))/(max(ord)-min(ord))
      return(sce)
    }else
    {
      stop("The \"normalization\" must be one of the followwings: \"TMM\" / \"quant\" / \"cpm\"")
    }
  }else
  {
    stop("The \"type\" must be one of the followwings: \"counts\" / \"normalized\" / \"expression\"")
  }
}


#Normalizes expression data based on. The data can be optionally log-transformed, pseudoc-counts can optionally added and one can choose to use relative expression.or actual expression.
normalize_expr_data <- function(cds,
                                norm_method = c("log",  "none"),
                                pseudo_expr = 1,
                                relative_expr = TRUE){
  FM <- counts(cds)
  norm_method <- match.arg(norm_method)
  #checkSizeFactors(cds)

  if (norm_method == "log") {
    # If we are using log, normalize by size factor before log-transforming

    if (relative_expr)
      FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))

    FM <- FM + pseudo_expr
    FM <- log2(FM)
  }else if (norm_method == "none"){
    # If we are using log, normalize by size factor before log-transforming
    if (relative_expr)
    {
      FM <- Matrix::t(Matrix::t(FM)/sizeFactors(cds))
    }
    FM <- FM + pseudo_expr
  }

  return (FM)
}

#Performs normalization on expression matrix.
normalise<-function(cds)
{
  data<-counts(cds)
  data3<-normalize_expr_data(cds,norm_method="log",pseudo_expr=0.5,relative_expr=TRUE)
  return(data3)
}

#The main function to compute pseudotime along with normalization on count data
findorders_counts<-function(data,ntest,minit,epsilon,nnprop)
{
  library(DESeq2)
  library(Matrix)
  data<-data[(rowSums(abs(data))>0),]
  conditions<-rep(1,dim(data)[2])
  cds<-DESeqDataSetFromMatrix(data,DataFrame(conditions),~1)
  cds<-estimateSizeFactors(cds)
  data<-normalise(cds)
  data<-t(data)
	data<-data[,!is.nan(colSums(data))]
	data<-data[,!is.na(colSums(data))]
	data1<-data
	data1<-apply(data,2,"rank",tie="average")
	data<-data1
	n<-dim(data)[1]
	Y<<-data
	yy<-colSums(Y^2)
	sumsqr<<-replicate((8*ntest),yy)
	yy<-(colMeans(Y))^2
	sqrsum<<-n*replicate((8*ntest),yy)
	matt<-polyfit(data,ntest,minit,epsilon)
	pseudotime<-match(1:(dim(data)[1]),matt)
	dd<-as.matrix(dist(data))
	dd1<-t(apply(dd,1,"order"))
	dd1<-dd1[,(1:(floor(nnprop*n)+1))]
	pst<-apply(dd1,1,"perm",data=t(t(pseudotime)))
	return(colMeans(pst))
}

#The main function to compute pseudotime on data which has been already normalized.
findorders_normalized<-function(data,ntest,minit,epsilon,nnprop)
{
  library(Matrix)
  data<-data[(rowSums(abs(data))>0),]
  data<-t(data)
  data<-data[,!is.nan(colSums(data))]
  data<-data[,!is.na(colSums(data))]
  data1<-data
  data1<-apply(data,2,"rank",tie="average")
  data<-data1
  n<-dim(data)[1]
  Y<<-data
  yy<-colSums(Y^2)
  sumsqr<<-replicate((8*ntest),yy)
  yy<-(colMeans(Y))^2
  sqrsum<<-n*replicate((8*ntest),yy)
  matt<-polyfit(data,ntest,minit,epsilon)
  pseudotime<-match(1:(dim(data)[1]),matt)
  dd<-as.matrix(dist(data))
  dd1<-t(apply(dd,1,"order"))
  dd1<-dd1[,(1:(floor(nnprop*n)+1))]
  pst<-apply(dd1,1,"perm",data=t(t(pseudotime)))
  return(colMeans(pst))
}

#Function to apply genetic algorithm
polyfit<-function(data,ntest,minit,epsilon)
{
  data<-data[,(apply(data,2,"sd")>0)]
  n<-dim(data)[1]
  Y<<-data
  yy<-colSums(Y^2)
  sumsqr<<-replicate((8*ntest),yy)
  yy<-(colMeans(Y))^2
  sqrsum<<-n*replicate((8*ntest),yy)
  X<-initialize(data,ntest)
  X4<-X
  iter<-1
  costmat<-array(0,dim<-c(minit,(8*ntest)))
  while(iter<=minit)
  {
    #print(iter)
    X<-X4
    X1<-X
    X21<-recombination(X1)
    X2<-rbind(X21,X1)
    X3<-mutation(X2)
    if(iter==1)
    {
      X4<-selection1(rbind(X2,X3),data)
    }else
    {
      X4<-selection2(rbind(X2,X3),data)
    }
    costmat[iter,]<-cost
    iter<-iter+1
  }
  eps<-abs(min(costmat[(minit-1),])-min(costmat[minit,]))/prod(dim(data))
  cost1<-min(costmat[(minit-1),])
  cost2<-min(costmat[minit,])
  while(eps>epsilon)
  {
    X<-X4
    X1<-X
    X21<-recombination(X1)
    X2<-rbind(X21,X1)
    X3<-mutation(X2)
    if(iter==1)
    {
      X4<-selection1(rbind(X2,X3),data)
    }else
    {
      X4<-selection2(rbind(X2,X3),data)
    }
    cost1<-cost2
    cost2<-min(cost)
    eps<-abs(cost1-cost2)/prod(dim(data))
    iter<-iter+1
  }

  return(X4[1,])
}

#Initializing random population
initialize<-function(data,ntest)
{
	n<-dim(data)[1]
	X<-array(0,dim<-c((2*ntest),n))
	for(i in 1:(dim(X)[1]))
	{
		X[i,]<-sample(1:n)
	}
	return(X)
}

#Helper function to permute rows of a matrix
perm<-function(y,data)
{
  return(data[y,])
}

#Calculates BIC from residuals
findbicnew<-function(res,k=1,n=300)
{
  w<-rep(1,n)
  ll<-0.5 * (sum(log(w)) - n * (log(2 * pi) + 1 - log(n) + log(res)))
  val<-(k+2)*log(n)-2*ll
  return(val)
}

#Selection performed on a population based on fitness on first generation
selection1<-function(X,data)
{
N<-dim(X)[1]
cost<-matrix(nrow=N,ncol=1)
n<-dim(data)[1]
k<-dim(data)[2]
XX<-apply(X,1,"invPerm")
dat<-as.numeric(data)
x1<-1:n
basis1<-(x1-mean(x1))/sqrt(sum((x1-mean(x1))^2))
expls1<-apply(t(XX),1,"perm",data=t(t(basis1)))
x2<-(1:n)^2
resid2<-lm(x2~x1)$residuals
basis2<-resid2/sqrt(sum(resid2^2))
expls2<-apply(t(XX),1,"perm",data=t(t(basis2)))
x3<-(1:n)^3
resid3<-lm(x3~x1+x2)$residuals
basis3<-resid3/sqrt(sum(resid3^2))
expls3<-apply(t(XX),1,"perm",data=t(t(basis3)))

beta<-crossprod(Y,cbind(expls1,expls2,expls3))

resid1<-sumsqr-sqrsum-(beta[,(1:N)])^2
resid2<-resid1-(beta[,(N+(1:N))])^2
resid3<-resid2-(beta[,(2*N+(1:N))])^2

bic1<-findbicnew(resid1,1,n)
bic2<-findbicnew(resid2,2,n)
bic3<-findbicnew(resid3,3,n)

bic<-pmin(bic1,bic2,bic3)

cost<<-colSums(bic)

X1<-X[(order(cost))[(1:(N/4))],]
return(X1)
}

#Selection performed on a population based on fitness on second generation onwards
selection2<-function(X,data)
{
  N<-dim(X)[1]
  X1<-X[(1:(N/4)),]
  X<-X[((N/4+1):N),]
  cost1<-cost[1:(length(cost)/4)]
  n<-dim(data)[1]
  k<-dim(data)[2]
  XX<-apply(X,1,"invPerm")
  dat<-as.numeric(data)
  x1<-1:n
  basis1<-(x1-mean(x1))/sqrt(sum((x1-mean(x1))^2))
  expls1<-apply(t(XX),1,"perm",data=t(t(basis1)))
  x2<-(1:n)^2
  resid2<-lm(x2~x1)$residuals
  basis2<-resid2/sqrt(sum(resid2^2))
  expls2<-apply(t(XX),1,"perm",data=t(t(basis2)))
  x3<-(1:n)^3
  resid3<-lm(x3~x1+x2)$residuals
  basis3<-resid3/sqrt(sum(resid3^2))
  expls3<-apply(t(XX),1,"perm",data=t(t(basis3)))

  beta<-crossprod(Y,cbind(expls1,expls2,expls3))

  resid1<-sumsqr[,(1:(3*N/4))]-sqrsum[,(1:(3*N/4))]-(beta[,(1:(3*N/4))])^2
  resid2<-resid1[,(1:(3*N/4))]-(beta[,((3*N/4)+(1:(3*N/4)))])^2
  resid3<-resid2[,(1:(3*N/4))]-(beta[,(2*(3*N/4)+(1:(3*N/4)))])^2

  bic1<-findbicnew(resid1,1,n)
  bic2<-findbicnew(resid2,2,n)
  bic3<-findbicnew(resid3,3,n)

  bic<-pmin(bic1,bic2,bic3)

  cost<<-c(cost1,colSums(bic))

  X<-rbind(X1,X)
  X1<-X[(order(cost))[(1:(N/4))],]
  return(X1)
}

#Recombination applied on parents
recomb<-function(y1)
{
  y<-y1[(1:(length(y1)/2))]
  z<-y1[(length(y1)/2+1):(length(y1))]
  N<-length(y)
  if(cor(y,z)<0)
  {
    z<-rev(z)
  }
  u<-rpois(1,1)
  while(u==0)
  {
    u<-rpois(1,3)
  }
  stepsize<-floor(N/(u+1))
  y1<-rep(NA,length(y))
  z1<-rep(NA,length(z))
  for(j in 1:(u+1))
  {
    if(j%%2==1)
    {
      if(j<(u+1))
      {
        y1[((j-1)*stepsize+1):(j*stepsize)]<-y[((j-1)*stepsize+1):(j*stepsize)]
      }else
      {
        y1[((j-1)*stepsize+1):(length(y))]<-y[((j-1)*stepsize+1):(length(y))]
      }
    }else
    {
      if(j<(u+1))
      {
        z1[((j-1)*stepsize+1):(j*stepsize)]<-z[((j-1)*stepsize+1):(j*stepsize)]
      }else
      {
        z1[((j-1)*stepsize+1):(length(z))]<-z[((j-1)*stepsize+1):(length(z))]
      }
    }
  }
  sety1<-setdiff(sample(length(y)),na.omit(y1))
  setz1<-setdiff(sample(length(z)),na.omit(z1))
  posz<-match(sety1,z)
  posy<-match(setz1,y)
  oriposy<-match(sety1,y)
  oriposz<-match(setz1,z)

  y1[sort(oriposy)]<-sety1[order(posz)]
  z1[sort(oriposz)]<-setz1[order(posy)]
  return(c(y1,z1))
}


#Recombination applied on a matrix whose rows are permutations
recombination<-function(X1)
{
	N<-dim(X1)[2]
	n<-dim(X1)[1]
  samp<-sample(1:n,(n/2))
	X2<-X1[samp,]
	X3<-X1[setdiff(1:n,samp),]
  X<-t(apply(cbind(X2,X3),1,"recomb"))
	return(rbind(X[,(1:N)],X[,((N+1):(2*N))]))
}

#Mutation operation applied on a population
mutation<-function(X2)
{
	N<-dim(X2)[1]
	n<-dim(X2)[2]
	for(i in 1:N)
	{
		samp<-sample(1:n,2)
		u1<-min(samp)
		u2<-max(samp)
		X2[i,(u1:u2)]<-rev(X2[i,(u1:u2)])
	}
	return(X2)
}







