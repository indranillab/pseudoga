# Function to determine which termination points of paths from different clusters should be joined.
#
# @param clster_id Cluster ids of cells.
# @param pseudotime Pseudotime information of individuals within their own clusters.
# @param data Matrix of gene expressions.
#
# @return Returns the adjacent matrix indicating the terminationation points that should be joined in the tree.
#
join_cluster<-function(cluster_id,pseudotime,data)
{
	M<-max(cluster_id)
	A<-array(0,dim<-c(M,M))
	B<-array(0,dim<-c((2*M),(2*M)))
	join_index<-array(0,dim<-c(M,M))
	for(i in 1:M)
	{
		for(j in 1:M)
		{
			if(i==j)
			{
				next
			}
			ps1<-which(cluster_id==i)
			ps2<-which(cluster_id==j)
			pst1<-pseudotime[ps1]
			pst2<-pseudotime[ps2]
			l1<-length(pst1)
			l2<-length(pst2)
			x11<-match((1:(floor(l1/20))),pst1)
			x12<-match(((floor(19*l1/20)):(floor(l1))),pst1)
			x21<-match((1:(floor(l2/20))),pst2)
			x22<-match(((floor(19*l2/20)):(floor(l2))),pst2)
			dat1<-data[ps1,]
			dat2<-data[ps2,]
			y11<-dat1[x11,]
			y12<-dat1[x12,]
			y21<-dat2[x21,]
			y22<-dat2[x22,]
			s1<-0
			count1<-0
			for(i1 in 1:(dim(y11)[1]))
			{
				for(j1 in 1:(dim(y21)[1]))
				{
					s1<-s1+sqrt(sum((y11[i1,]-y21[j1,])^2))
					count1<-count1+1
				}
			}
			B[i,j]<-s1/count1
			s2<-0
			count2<-0
			for(i1 in 1:(dim(y11)[1]))
			{
				for(j1 in 1:(dim(y22)[1]))
				{
					s2<-s2+sqrt(sum((y11[i1,]-y22[j1,])^2))
					count2<-count2+1
				}
			}
			B[i,(M+j)]<-s2/count2
			s3<-0
			count3<-0
			for(i1 in 1:(dim(y12)[1]))
			{
				for(j1 in 1:(dim(y21)[1]))
				{
					s3<-s3+sqrt(sum((y12[i1,]-y21[j1,])^2))
					count3<-count3+1
				}
			}
			B[(M+i),j]<-s3/count3
			s4<-0
			count4<-0
			for(i1 in 1:(dim(y12)[1]))
			{
				for(j1 in 1:(dim(y22)[1]))
				{
					s4<-s4+sqrt(sum((y12[i1,]-y22[j1,])^2))
					count4<-count4+1
				}
			}
			B[(M+i),(M+j)]<-s4/count4
		}
	}

	MM<-max(B)
	for(i in 1:M)
	{
		for(j in 1:M)
		{
			if(i==j)
			{
				B[i,j]<-2*MM
				B[i,(M+j)]<-2*MM
				B[(M+i),j]<-2*MM
				B[(M+i),(M+j)]<-2*MM
			}
		}
	}

	for(i in 1:M)
	{
		for(j in 1:M)
		{
			if(i!=j)
			{
			A[i,j]<-min(B[i,j],B[(M+i),j],B[i,(M+j)],B[(M+i),(M+j)])
			vec<-c(B[i,j],B[(M+i),j],B[i,(M+j)],B[(M+i),(M+j)])
			join_index[i,j]<-which(vec==min(vec))
			}
		}
	}


	MM<-max(A)
	diag(A)<-2*MM
	graph<-array(0,dim<-dim(A))
	index<-cbind(rep(1:(dim(A)[1]),(dim(A)[2])),rep(1:(dim(A)[1]),each=(dim(A)[2])),as.numeric(A))
	index1<-index[(order(index[,3])),]
	kk<-1
	while(1)
	{
		xxx<-index1[kk,]
		graph[(xxx[1]),(xxx[2])]<-1
		graph[(xxx[2]),(xxx[1])]<-1
		if(form_cycle(graph))
		{
			graph[(xxx[1]),(xxx[2])]<-0
		    graph[(xxx[2]),(xxx[1])]<-0
		}
		if(sum(graph)==(2*(dim(A)[1]-1)))
		{
			break
		}
		kk<-kk+1
	}


	return(graph*join_index)
}

# Function to check whether there is any cycle in a given tree structure of pseudotime paths.
#
# @param graph The adjacent matrix representing the tree structure of pseudotime paths.
#
# @return Returns a flag denoting whether there is any cycle in the given tree structure.
#

form_cycle<-function(graph)
{
	dd<-dim(graph)
	flag<-0
	touch<-rep(-1,(dim(graph)[1]))
	for(i in 1:(dim(graph)[1]-1))
	{
		for(j in (i+1):(dim(graph)[2]))
		{
			if(graph[i,j]>0)
			{
				if((touch[i]==1)&(touch[j]==1))
				{
					flag<-1
				}
				touch[i]<-1
				touch[j]<-1
			}
		}
	}
	return(flag)

}



