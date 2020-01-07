Heritability <-
function(ir,trait,alpha){
  n <- rep(0,ir)
  y <- as.matrix(as.numeric(trait[,2]))
  y <- na.omit(y)
  n[1]=nrow(y)
  for(i in 2:ir){
    y1 <- as.matrix(as.numeric(trait[,i*2]))
    y1 <- na.omit(y1)
    n[i]=nrow(y1)
    y <- rbind(y,y1,deparse.level=1) 
  }
  n1=0
  for(i in 1:ir){
    n1=n1+n[i]
  }
  X=matrix(rep(0,ir*n1),ncol=ir)
  for(i in 1:n[1]){
    X[i,1]=1
  }
  for(i in (n[1]+1):(n[1]+n[2])){
    X[i,2]=1
  }
  for(i in (n[1]+n[2]+1):(n[1]+n[2]+n[3])){
    X[i,3]=1
  }	
  L <- as.matrix(as.numeric(trait[,1]))
  L <- na.omit(L)
  for(i in 2:ir){
    L1 <- as.matrix(as.numeric(trait[,2*i-1]))
    L1 <- na.omit(L1)
    L <- rbind(L,L1,deparse.level=1) 
  }
  r=0.5-0.5*alpha+0.25*alpha*alpha
  A=matrix(r,nrow=n1,ncol=n1)
  for(i in 1:n1){
    for(j in 1:n1){
      if(L[i]==L[j]){
        A[i,j]=1
      }
    }
  }
  fit=amvce(y,X,A)
  h2n=(fit$VC[1]/(fit$VC[1] + fit$VC[2]))
  cat("Genetic variance:\n")
  print(fit$VC[1])
  cat("Residual variance:\n")
  print(fit$VC[2])
  cat("Estimated heritability:\n")
  print(h2n)
  return(list(fit$VC[1],fit$VC[2],h2n))
}
