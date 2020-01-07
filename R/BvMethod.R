BvMethod <-
function(n,mloci,ftype,mtype,rec,o,trait){
  mxind=500
  mxloci=400
  temp1 <- matrix(rep(0,mxloci*4),nrow=mxloci,ncol=4)
  temp2 <- matrix(rep(0,mxloci*4),nrow=mxloci,ncol=4)
  for(i in 1:mloci){
    for(j in 1:4){
      temp1[i,j]=ftype[i,j]
      temp2[i,j]=mtype[i,j]
    }
  }
  ftype <- temp1
  mtype <- temp2
  temp4 <- as.double(rep(0.0,mxloci))
  for(i in 1:mloci){
    temp4[i] <- rec[i]
  }
  rec <- temp4
  temp5 <- array(rep(0,mxind*mxloci*8),dim=c(mxind,mxloci,8))
  temp6 <- as.double(rep(0.0,mxind))
  for(i in 1:n){
    temp6[i] <- trait[i]
    for(j in 1:mloci){
      for(k in 1:8){
        temp5[i,j,k] <- o[i,j,k]
      }
    }
  }
  trait <- temp6
  o <- temp5
  fqtl_max <- rep(0,4)
  mqtl_max <- rep(0,4)
  ProfileLODs <- as.double(rep(0.0,mxloci))
  GeneticDistance <- as.double(rep(0.0,mxloci))
  u_max <- rep(0.0,5)
  y <- seq(1,n)
  for(i in 1:n){
    y[i] <- trait[i]
  }
  kmeans1 <- matrix(seq(1,20),nrow=4,ncol=5)
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
   out <- .Fortran("bvmethod",n=as.integer(n),mloci=as.integer(mloci),ftype=as.integer(ftype),mtype=as.integer(mtype),rec=as.double(rec),o=as.integer(o),trait=as.double(trait),kmeans1=as.double(kmeans1),fqtl_max=as.integer(fqtl_max),mqtl_max=as.integer(mqtl_max),BIC_min=as.double(0.0),ProfileLODs=as.double(ProfileLODs),GeneticDistance=as.double(GeneticDistance),gDistance_max=as.double(0.0),u_max=as.double(u_max),sigma_max=as.double(0.0),PACKAGE="TetraQTLAnalysis")
  index <- 0
  for(i in 1:mxloci){
    if(!is.na(out$ProfileLODs[i])){
      if(out$ProfileLODs[i]>0.0){
        index <- index+1
      }
    }
  }
  GD <- rep(0.0,index)
  PL <- rep(0.0,index)
  for(i in 1:index){
    GD[i] <- out$GeneticDistance[i]
    PL[i] <- out$ProfileLODs[i]
  }
  plot(GD,PL,type="l",main="Profile of LOD scores",ylab="LOD Scores",xlab="Position (cM)",lwd=2)
  cat("The most likely predicted parental QTL genotypes:\n")
  print(out$fqtl_max)
  print(out$mqtl_max)
  cat("BIC:\n")
  print(out$BIC_min)
  cat("Predicted location (cM):\n")
  print(out$gDistance_max)
  cat("G1:\n")
  print(out$u_max[1])
  cat("G2:\n")
  print(out$u_max[2])
  cat("G3:\n")
  print(out$u_max[3])
  cat("G4:\n")
  print(out$u_max[4])
  cat("G5:\n")
  print(out$u_max[5])
  cat("Environmental residual:\n")
  print(out$sigma_max)                         
  return(list(out$fqtl_max,out$mqtl_max,out$BIC_min,out$gDistance_max,out$u_max,out$sigma_max))
}
