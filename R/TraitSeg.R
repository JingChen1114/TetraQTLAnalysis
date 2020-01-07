TraitSeg <-
function(n,trait){
  kmeans1 <- matrix(seq(1,20),nrow=4,ncol=5)
  y <- seq(1,n)
  for(i in 1:n){
    y[i] <- trait[i]
  }
  for(l in 1:4){
    k <- l+1
    a <- kmeans(y,k,iter.max=100)
    
    for(m in k:2){                    #bubble sort
      i1 <- 1
      while(i1<m){
        if(a$centers[i1]>a$centers[i1+1]){
          temp<-a$centers[i1+1]
          a$centers[i1+1] <- a$centers[i1]
          a$centers[i1] <- temp
        }
        i1 <- i1+1
      }
    }
    for(i in 1:k){
      kmeans1[l,i] <- a$centers[i]
    }
  }
  pa1 <- rep(0,4)
  pa2 <- rep(0,4)
  fre <- rep(0.0,5)
  g <- rep(0.0,5)
  trait1 <- rep(0.0,1000)
  for(i in 1:n){
    trait1[i]=trait[i]
  }
  out <- .Fortran("traitseg",n=as.integer(n),trait=as.double(trait1),kmean=as.double(kmeans1),pa1=as.integer(pa1),pa2=as.integer(pa2),fre=as.double(fre),g=as.double(g),e=as.double(0.0),alpha=as.double(0.0),df=as.integer(0),sLRT=as.double(0.0),PACKAGE="TetraQTLAnalysis")
  hist(trait,xlim=range(trait),xlab="",nclass=30,freq=F)
  max1=max(trait)
  min1=min(trait)
  x1=seq(min1,max1,0.01)
  f0=out$fre[1]
  f1=out$fre[2]
  f2=out$fre[3]
  f3=out$fre[4]
  f4=out$fre[5]
  g0=out$g[1]
  g1=out$g[2]
  g2=out$g[3]
  g3=out$g[4]
  g4=out$g[5]
  e=out$e
  y1=(1/(((2*pi)^(1/2))*e))*(f0*exp(-((x1-g0)^2)/(2*(e^2)))+
      f1*exp(-((x1-g1)^2)/(2*(e^2)))+f2*exp(-((x1-g2)^2)/(2*
      (e^2)))+f3*exp(-((x1-g3)^2)/(2*(e^2)))+f4*exp(-((x1-g4)^2)
      /(2*(e^2))))
  lines(x1,y1,type="l",col="red",lwd=2)

  y10=(1/(((2*pi)^(1/2))*e))*(f0*exp(-((x1-g0)^2)/(2*(e^2))))
  lines(x1,y10,type="p",col="blue",cex=0.1)

  y11=(1/(((2*pi)^(1/2))*e))*(f1*exp(-((x1-g1)^2)/(2*(e^2))))
  lines(x1,y11,type="p",col="blue",cex=0.1)

  y12=(1/(((2*pi)^(1/2))*e))*(f2*exp(-((x1-g2)^2)/(2*(e^2))))
  lines(x1,y12,type="p",col="blue",cex=0.1)

  y13=(1/(((2*pi)^(1/2))*e))*(f3*exp(-((x1-g3)^2)/(2*(e^2))))
  lines(x1,y13,type="p",col="blue",cex=0.1)

  y14=(1/(((2*pi)^(1/2))*e))*(f4*exp(-((x1-g4)^2)/(2*(e^2))))
  lines(x1,y14,type="p",col="blue",cex=0.1)
  p <- pchisq(out$sLRT,df=out$df)
  if(p>0.5){
    p_value=2*(1-p)
  }else if(p<=0.5){
    p_value=2*p
  }
  sum=0.0
  for(i in 1:n){
    sum=sum+trait[i]
  }
  tr_mean=sum/n
  sum=0.0
  for(i in 1:n){
    sum=sum+(trait[i]-tr_mean)^2
  }
  variance=sum/(n-1)
  g_var=variance-e^2
  cat("Estimated parental configuration:\n")
  print(out$pa1)
  print(out$pa2)
  cat("Estimated coefficient of double reduction:\n")
  print(out$alpha)
  cat("Significance of major QTL segregation (P value):\n")
  print(p_value)
  cat("Total variance:\n")
  print(variance)
  cat("Variance contributed by this major QTL:\n")
  print(g_var)
  return(list(out$p1,out$p2,out$alpha,p_value,variance,g_var))
}
