BvThreshold <-
function(n,mloci,ftype,mtype,rec,o,trait,fqtl,mqtl,ir){
  mxind=500
  mxloci=400
  iRS=1000
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
  y <- matrix(rep(0,iRS*mxind),nrow=iRS,ncol=mxind)
  for(i in 1:ir){
    for(j in 1:n){
     y[i,j]=trait[j]
    }
    for(j in seq(n,2,-1)){
      k <- sample(1:j,1)
      temp7 <- y[i,j]
      y[i,j] <- y[i,k]
      y[i,k] <- temp7 
    }
  }
  sLOD_max <- rep(0,iRS)
  out <- .Fortran("bvthreshold",n=as.integer(n),mloci=as.integer(mloci),ftype=as.integer(ftype),mtype=as.integer(mtype),rec=as.double(rec),o=as.integer(o),y=as.double(y),fqtl=as.integer(fqtl),mqtl=as.integer(mqtl),ir=as.integer(ir),sLOD_max=as.double(sLOD_max),ilocus=as.integer(0),PACKAGE="TetraQTLAnalysis")  
  sLOD <- rep(0,out$ilocus)
  for(i in 1:out$ilocus){
    sLOD[i] <- out$sLOD_max[i]
  }
  threshold <- quantile(sLOD,0.95) 
  return(threshold)
}
