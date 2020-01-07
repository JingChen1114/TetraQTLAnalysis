amvce <-
function(y,X,A,tolerance=1.0e-3,maxIter=100,verbose=FALSE){
#Animal model variance component estimation
#Author: Minghui Wang <m.h.wang@live.com>
	y=as.vector(y);X=as.matrix(X);A=as.matrix(A)
	n=dim(X)[1]
	NP=dim(X)[2]
	if(n!=length(y) || n!=nrow(A) || n!=ncol(A)) stop('Invalid input\n')
	R=diag(n)		#200 x 200 with 1's on the diagonal
	ireduce<-function(X){
		n=nrow(X)
		NP=ncol(X)
		ij=do.call("order", split(X, col(X)))
		iRemv=rep(FALSE,n)
		nUr=1
		lstEt=NA
		for(i in 2:n){
			if(all(X[ij[i],]==X[ij[i-1],])){
				lstEt=ij[i]
				if(i==n) iRemv[lstEt]=T
			}else{
				if(!is.na(lstEt)) iRemv[lstEt]=T
				lstEt=NA
				nUr=nUr+1
			}
		}
		if(n-nUr<NP){
			iRemv[(1:n)[!iRemv][1:(NP+nUr-n)]]=T
		}
		return(iRemv)
	}
	ginv <- function(X, tol = sqrt(.Machine$double.eps)){
	#borrowed from library MASS
		if(length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
			stop("'X' must be a numeric or complex matrix")
		if(!is.matrix(X)) X <- as.matrix(X)
		Xsvd <- svd(X)
		if(is.complex(X)) Xsvd$u <- Conj(Xsvd$u)
		Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
		if (all(Positive)) Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
		else if(!any(Positive)) array(0, dim(X)[2L:1L])
		else Xsvd$v[, Positive, drop=FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop=FALSE]))
	}
	ptr<-function(A,B){
		s=0.0
		for(i in 1:nrow(A)){
			s=s+sum(A[i,]*B[,i])
		}
		return(s)
	}
	logdet<-function(a) determinant(a,logarithm=TRUE)[[1]][1]
	M=diag(n) - X %*% tcrossprod(ginv(crossprod(X,X)),X)
	iRemoval=ireduce(X)
	K1=M[!iRemoval,]
	K=t(K1)
	Xr=X[!iRemoval,]
	nIter=0
	if(verbose==T) cat('Running REML algorithm...\nIter.\tlogL\tV(u)\tV(e)\n')
	loglik0=NA
	loglik=NA
	V=NA
	K3=K %*% solve(K1 %*% K,K1)
	VC=var(y)*(1:2/3)
	while(nIter<=maxIter){
		V=VC[1] * A + VC[2] * R
		V=K1 %*% V %*% K
		Vi=solve(V)
		P=crossprod(K1,Vi %*% K1)
		P1=P %*% y
		tum=VC[1]*n-VC[1]^2 * ptr(P,A)
		tum=tum+VC[1]^2.0 * (crossprod(P1,A %*% P1))
		tem=VC[2]*sum(diag(K3))-VC[2]^2*ptr(K3,P)
		tem=tem+VC[2]^2*(crossprod(P1,K3 %*% P1))
		VC[1]=as.vector(tum/n)
		VC[2]=as.vector(tem/(n-NP))
		nIter=nIter+1
		X3=crossprod(Xr,Vi) %*% Xr
		loglik=-0.5*(logdet(V)+logdet(X3)+crossprod(y,P1)[,1])
		if(verbose==T) cat(nIter,'\t',sprintf("%.4f",loglik),'\t',paste(sprintf("%.3f",VC),collapse='\t'),'\n',sep='')
		if(nIter>1 && abs(loglik-loglik0)<tolerance) break
		loglik0=loglik
	}
	if(nIter>maxIter){
		warning(paste('Failed to converge after',nIter,'iterations'))
	}
	Vi=solve(VC[1] * A + VC[2] * R)
	b=tcrossprod(solve(crossprod(X,Vi) %*% X),X) %*% (Vi %*% y)[,1]
	res=list(coefficients=b,VC=VC,loglik=loglik)
	return(res)
}
