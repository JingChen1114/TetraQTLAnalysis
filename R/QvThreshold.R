QvThreshold <-
function(n,mloci,ftype,mtype,alpha,rec,o,trait,fqtl,mqtl,ir){
  mxind=500
  mxloci=400
  iRS=1000
  max_g=136*136
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
  temp3 <- as.double(rep(0.0,mxloci))
  temp4 <- as.double(rep(0.0,mxloci))
  for(i in 1:mloci){
    temp4[i] <- rec[i]
    temp3[i] <- alpha[i]
  }
  alpha <- temp3
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
  g_o <- matrix(rep(0,mxind*max_g),nrow=mxind,ncol=max_g)
  q_g <- array(rep(0,mxind*max_g*5),dim=c(mxind,max_g,5))
  n_ig <- rep(0,mxind)
  ip1 <- array(rep(0,mxind*max_g*2),dim=c(mxind,max_g,2))
  ip2 <- ip1
  sLOD_t <- rep(0,iRS)
  sLOD_tt <- matrix(rep(0,iRS*500),nrow=iRS,ncol=500)
  sLOD_max <- rep(0,iRS)
  f_type <- matrix(rep(0,12),nrow=3,ncol=4)
  m_type <- matrix(rep(0,12),nrow=3,ncol=4)
  qtl_ge <- rep(0,5)
  ilocus=0
  for(i in 2:mloci){
    out1 <- .Fortran("qvthreshold1",i=as.integer(i),n=as.integer(n),mloci=as.integer(mloci),ftype=as.integer(ftype),mtype=as.integer(mtype),alpha=as.double(alpha),rec=as.double(rec),o=as.integer(o),ip1=as.integer(ip1),ip2=as.integer(ip2),g_o=as.double(g_o),n_ig=as.integer(n_ig),PACKAGE="TetraQTLAnalysis") 
    rec1=0.005
    r=rec[i]
    ip1 <- out1$ip1
    ip2 <- out1$ip2
    g_o <- out1$g_o
    n_ig <- out1$n_ig
    while(rec1<r){
      ilocus=ilocus+1
      ii=i-1
      for(k in 1:4){
        f_type[1,k] <- ftype[ii,k]
        f_type[3,k] <- ftype[i,k]
        m_type[1,k] <- mtype[ii,k]
        m_type[3,k] <- mtype[i,k]
        f_type[2,k] <- fqtl[k]
        m_type[2,k] <- mqtl[k]
      }
      for(j in 1:n){
        i_g <- n_ig[j]
        for(iC in 1:i_g){
          if1 <- ip1[j+mxind*(iC-1)]
          if2 <- ip1[j+mxind*(max_g+iC-1)]
          im1 <- ip2[j+mxind*(iC-1)]
          im2 <- ip2[j+mxind*(max_g+iC-1)]
          out2 <- .Fortran("qvthreshold2",f_type=as.integer(f_type),m_type=as.integer(m_type),r=as.double(r),rec1=as.double(rec1),if1=as.integer(if1),if2=as.integer(if2),im1=as.integer(im1),im2=as.integer(im2),iC=as.integer(iC),qtl_ge=as.double(qtl_ge),PACKAGE="TetraQTLAnalysis")
          for(k in 1:5){
            q_g[j,iC,k] <- out2$qtl_ge[k]
          }
        }
      }
      out3 <- .Fortran("qvthreshold3",ir=as.integer(ir),fqtl=as.integer(fqtl),mqtl=as.integer(mqtl),n=as.integer(n),y=as.double(y),g_o=as.double(g_o),q_g=as.double(q_g),n_ig=as.integer(n_ig),sLOD_t=as.double(sLOD_t),PACKAGE="TetraQTLAnalysis")
      for(i1 in 1:ir){
        sLOD_tt[i1,ilocus]=out3$sLOD_t[i1]
      }
      rec1=rec1+0.01
    }
  }
  for(i in 1:ir){
    sLOD_max[i]=max(sLOD_tt[i,])
  }
  threshold <- quantile(sLOD_max,0.95)
  return(threshold)
}
