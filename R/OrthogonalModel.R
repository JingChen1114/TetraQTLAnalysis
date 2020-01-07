OrthogonalModel <-
function(P1,P2,G,alpha){
  fre_qtl <- rep(0,5)
  out <- .Fortran("qtldistri",fqtl=as.integer(P1),mqtl=as.integer(P2),a=as.double(alpha),fre_qtl=as.double(fre_qtl),PACKAGE="TetraQTLAnalysis")
  igenotype=0
  for(i in 1:5){
    if(out$fre_qtl[i]>0.0){
      igenotype=igenotype+1
    }
  }
  fre <- out$fre_qtl
  w <- matrix(rep(0,20),nrow=4,ncol=5)
  x <- matrix(rep(0,25),nrow=5,ncol=5)
  if(igenotype==2){
    orthogonal <- matrix(rep(1,4),nrow=2,ncol=2)
    y <- rep(0,2)   
    u=-4.0*fre[5]-3.0*fre[4]-2.0*fre[3]-fre[2]
    for(i in 1:5){
      w[1,i]=u+i-1
    }
    ik=1
    for(i in 1:5){
      if(fre[i]>0){
        orthogonal[ik,2]=w[1,i]
        y[ik]=G[i]
        ik=ik+1
      }
    }
    Geffects <- solve(orthogonal,y)
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[1,i]
      s2=s2+fre[i]*w[1,i]*w[1,i]
    }
    Gvariance=s1*s1/s2
  }else if(igenotype==3){
    orthogonal <- matrix(rep(1,9),nrow=3,ncol=3)
    y <- rep(0,3)
    u=-4.0*fre[5]-3.0*fre[4]-2.0*fre[3]-fre[2]
    for(i in 1:5){
      w[1,i]=u+i-1
    }
    x[1,1]=1.0
    x[1,2]=-2.0
    x[1,3]=1.0
    x[2,2]=1.0
    x[2,3]=-2.0
    x[2,4]=1.0
    x[3,3]=1.0
    x[3,4]=-2.0
    x[3,5]=1.0
    for(i in 1:5){
      x[4,i]=fre[i]
      x[5,i]=w[1,i]*fre[i]
    }
    yi <- c(1,1,1,0,0)
    wi <- solve(x,yi)
    for(i in 1:5){
      w[2,i]=wi[i]
    }
    ik=1
    for(i in 1:5){
      if(fre[i]>0){
        orthogonal[ik,2]=w[1,i]
        orthogonal[ik,3]=w[2,i]
        y[ik]=G[i]
        ik=ik+1
      }
    }
    Geffects <- solve(orthogonal,y)
    Gvariance <- c(0,0)
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[1,i]
      s2=s2+fre[i]*w[1,i]*w[1,i]
    }
    Gvariance[1]=s1*s1/s2
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[2,i]
      s2=s2+fre[i]*w[2,i]*w[2,i]
    }
    Gvariance[2]=s1*s1/s2   
  }else if(igenotype==4){
    orthogonal <- matrix(rep(1,16),nrow=4,ncol=4)
    y <- rep(0,4)
    u=-4.0*fre[5]-3.0*fre[4]-2.0*fre[3]-fre[2]
    for(i in 1:5){
      w[1,i]=u+i-1
    }
    x[1,1]=1.0
    x[1,2]=-2.0
    x[1,3]=1.0
    x[2,2]=1.0
    x[2,3]=-2.0
    x[2,4]=1.0
    x[3,3]=1.0
    x[3,4]=-2.0
    x[3,5]=1.0
    for(i in 1:5){
      x[4,i]=fre[i]
      x[5,i]=w[1,i]*fre[i]
    }
    yi <- c(1,1,1,0,0)
    wi <- solve(x,yi)
    for(i in 1:5){
      w[2,i]=wi[i]
    }
    x <- matrix(rep(0,25),nrow=5,ncol=5)
    x[1,1]=-1.0
    x[1,2]=3.0
    x[1,3]=-3.0
    x[1,4]=1.0
    x[2,2]=-1.0
    x[2,3]=3.0
    x[2,4]=-3.0
    x[2,5]=1.0
    for(i in 1:5){
      x[3,i]=fre[i]
      x[4,i]=w[1,i]*fre[i]
      x[5,i]=w[2,i]*fre[i]
    }
    yi <- c(1,1,0,0,0)
    wi <- solve(x,yi)
    for(i in 1:5){
      w[3,i]=wi[i]
    }
    ik=1
    for(i in 1:5){
      if(fre[i]>0){
        orthogonal[ik,2]=w[1,i]
        orthogonal[ik,3]=w[2,i]
        orthogonal[ik,4]=w[3,i]
        y[ik]=G[i]
        ik=ik+1
      }
    }
    Geffects <- solve(orthogonal,y)
    Gvariance <- c(0,0,0)
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[1,i]
      s2=s2+fre[i]*w[1,i]*w[1,i]
    }
    Gvariance[1]=s1*s1/s2
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[2,i]
      s2=s2+fre[i]*w[2,i]*w[2,i]
    }
    Gvariance[2]=s1*s1/s2 
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[3,i]
      s2=s2+fre[i]*w[3,i]*w[3,i]
    }
    Gvariance[3]=s1*s1/s2
  }else if(igenotype==5){
    orthogonal <- matrix(rep(1,25),nrow=5,ncol=5)
    y <- rep(0,5)
    u=-4.0*fre[5]-3.0*fre[4]-2.0*fre[3]-fre[2]
    for(i in 1:5){
      w[1,i]=u+i-1
    }
    x[1,1]=1.0
    x[1,2]=-2.0
    x[1,3]=1.0
    x[2,2]=1.0
    x[2,3]=-2.0
    x[2,4]=1.0
    x[3,3]=1.0
    x[3,4]=-2.0
    x[3,5]=1.0
    for(i in 1:5){
      x[4,i]=fre[i]
      x[5,i]=w[1,i]*fre[i]
    }
    yi <- c(1,1,1,0,0)
    wi <- solve(x,yi)
    for(i in 1:5){
      w[2,i]=wi[i]
    }
    x <- matrix(rep(0,25),nrow=5,ncol=5)
    x[1,1]=-1.0
    x[1,2]=3.0
    x[1,3]=-3.0
    x[1,4]=1.0
    x[2,2]=-1.0
    x[2,3]=3.0
    x[2,4]=-3.0
    x[2,5]=1.0
    for(i in 1:5){
      x[3,i]=fre[i]
      x[4,i]=w[1,i]*fre[i]
      x[5,i]=w[2,i]*fre[i]
    }
    yi <- c(1,1,0,0,0)
    wi <- solve(x,yi)
    for(i in 1:5){
      w[3,i]=wi[i]
    }
    x <- matrix(rep(0,25),nrow=5,ncol=5)
    x[1,1]=1.0
 	  x[1,2]=-4.0
	  x[1,3]=6.0
	  x[1,4]=-4.0
	  x[1,5]=1.0
    for(i in 1:5){
      x[2,i]=fre[i]
      x[3,i]=w[1,i]*fre[i]
      x[4,i]=w[2,i]*fre[i]
      x[5,i]=w[3,i]*fre[i]
    }
    yi <- c(1,0,0,0,0)
    wi <- solve(x,yi)
    for(i in 1:5){
      w[4,i]=wi[i]
    }
    ik=1
    for(i in 1:5){
      if(fre[i]>0){
        orthogonal[ik,2]=w[1,i]
        orthogonal[ik,3]=w[2,i]
        orthogonal[ik,4]=w[3,i]
        orthogonal[ik,5]=w[4,i]
        y[ik]=G[i]
        ik=ik+1
      }
    }
    Geffects <- solve(orthogonal,y)
    Gvariance <- c(0,0,0,0)
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[1,i]
      s2=s2+fre[i]*w[1,i]*w[1,i]
    }
    Gvariance[1]=s1*s1/s2
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[2,i]
      s2=s2+fre[i]*w[2,i]*w[2,i]
    }
    Gvariance[2]=s1*s1/s2 
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[3,i]
      s2=s2+fre[i]*w[3,i]*w[3,i]
    }
    Gvariance[3]=s1*s1/s2
    s1=0
    s2=0
    for(i in 1:5){
      s1=s1+fre[i]*G[i]*w[4,i]
      s2=s2+fre[i]*w[4,i]*w[4,i]
    }
    Gvariance[4]=s1*s1/s2
  }
  cat("Genetic effects (mean,monogniec,digenic,trigenic and quadrigenic effects):\n")
  print(Geffects)
  cat("The corresponding genetic component variance:\n")
  print(Gvariance)
  return(list(Geffects,Gvariance))
}
