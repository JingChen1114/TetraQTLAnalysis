C     Purpose:
C             working out the threshold
C             quadrivalent pairing during meiosis 
C     Records of revision:
C     Date                Programmer
C     03/01/2019          Jing Chen
      subroutine qvthreshold3(ir,fqtl,mqtl,n,y,g_o,q_g,n_ig,
     /           sLOD_t)
      implicit double precision(a-h,o-z)
      
      parameter(iRS=1000,mxind=500,max_g=136*136)
      
      integer ir,n,mloci,kqtl
      integer n_ig(mxind)
      integer fqtl(4),mqtl(4),nqtl(5)
      
      dimension y(iRS,mxind),trait(mxind)
      dimension g_o(mxind,max_g)
      dimension q_g(mxind,max_g,5)
      dimension kmeans(5),u(5),sLOD_t(iRS)
      dimension fre_qtl(5)
      
      a=0.01d0
      call qtl_distribution1(fqtl,mqtl,a,fre_qtl,nqtl)
      kqtl=0
      do k=1,5
        if(nqtl(k)>0)then
          kqtl=kqtl+1
        end if
      end do
      
C     EM ALGORITHM FOR INTERVAL MAPPING 
      do k=1,ir
        do l=1,n
          trait(l)=y(k,l)
        end do  
        call emalgorithm1(fqtl,mqtl,n,trait,g_o,q_g,n_ig,kqtl,
     /       kmeans,u,sigma,slikelihood,sLOD)
        sLOD_t(k)=sLOD
      end do
      
	return
      end  
C     ***************************************************************
      subroutine emalgorithm1(fqtl,mqtl,n,trait,g_o,q_g,n_ig,kqtl,
     /           kmeans,u,sigma,slikelihood,sLOD)
      implicit double precision(a-h,o-z)
      
      parameter(mxind=500,max_g=136*136,pi=3.1415926d0)
      
      integer n_ig(mxind),kqtl
      integer fqtl(4),mqtl(4),nqtl(5)
      
      dimension g_o(mxind,max_g),q_g(mxind,max_g,5)
      dimension u1(5),u2(5),trait(mxind),g_theo(5)
      dimension ww(mxind,5),q(mxind,5)
      dimension fre_qtl(5)
      dimension u(5),kmeans(5)
            
    
      do i=1,5               
        u1(i)=kmeans(i)
      end do
      sum=0.0d0
      do i=1,n
        sum=sum+trait(i)
      end do
      gmean=sum/n
      u1=gmean
      sum=0.0d0
      do i=1,n                     
        sum=sum+(trait(i)-gmean)**2
      end do
      sigma1=sqrt(sum/(n-1))
      q=0.0d0
      do k=1,5
        do i=1,n
          ig=n_ig(i)
          do m=1,ig
            q(i,k)=q(i,k)+q_g(i,m,k)*g_o(i,m)
          end do
        end do
      end do        
      slikelihood1=0.0d0
      do i=1,n
        s2=0.0d0
        do j=1,n_ig(i)
          s1=0.0d0
          do k=1,5
            s1=s1+(1.0d0/sqrt(2.0d0*pi*(sigma1**2)))*exp(-((trait(i)
     /         -u1(k))**2)/(2.0d0*(sigma1**2)))*q_g(i,j,k)*g_o(i,j)
	    end do
	    s2=s2+s1
	  end do
	  slikelihood1=slikelihood1+log(s2)
	end do  	
	difference=1.0d0
	do while(difference>0.0001d0)  
	  ww=0.0d0
	  do i=1,n
	    do j=1,5
	      s2=0.0d0
	      do k=1,n_ig(i)
	        s1=0.0d0
	        do l=1,5
	          s1=s1+(1.0d0/sqrt(2.0d0*pi*(sigma1**2)))*exp(-((
     /             trait(i)-u1(l))**2)/(2.0d0*(sigma1**2)))
     /             *q_g(i,k,l)
	        end do
	        s2=s2+(((1.0d0/sqrt(2.0d0*pi*(sigma1**2)))*exp(-((
     /           trait(i)-u1(j))**2)/(2.0d0*(sigma1**2))))*
     /           q_g(i,k,j)*g_o(i,k))/s1
	      end do
	      ww(i,j)=s2        
	    end do
	  end do	  
	  s3=0.0d0
        u2=0.0d0
	  do j=1,5
	    s1=0.0d0
          do i=1,n
	      s1=s1+ww(i,j)
	    end do
	    s2=0.0d0
	    do i=1,n
	      s2=s2+ww(i,j)*trait(i)
	    end do
	    if(s1.ne.0.0d0)then
	      u2(j)=s2/s1      
	      do i=1,n
	        s3=s3+((trait(i)-u2(j))**2)*ww(i,j)
	      end do
	    end if
	  end do
        sigma2=sqrt(s3/(n*1.0d0))
        
        slikelihood2=0.0d0
	  do i=1,n
	    s2=0.0d0
	    do j=1,n_ig(i)
	      s1=0.0d0
	      do k=1,5
	        s1=s1+(1.0d0/sqrt(2.0d0*pi*(sigma2**2)))*exp(-((
     /           trait(i)-u2(k))**2)/(2.0d0*(sigma2**2)))*
     /           q_g(i,j,k)*g_o(i,j)
	      end do
	      s2=s2+s1
	    end do
	    slikelihood2=slikelihood2+log(s2)
	  end do
	  
	  difference=slikelihood2-slikelihood1
	  
	  slikelihood1=slikelihood2
	  do k=1,5
	    u1(k)=u2(k)
	  end do
	  sigma1=sigma2
      end do         
      do i=1,5        
        u(i)=u1(i)    
      end do
      sigma=sigma1
      slikelihood=slikelihood1  
      s1=0.0d0
	do i=1,n
	  s1=s1+trait(i)
	end do
	u_mean=s1/(n*1.0d0)
	s2=0.0d0
      do i=1,n
	  s2=s2+(trait(i)-u_mean)**2
	end do
	var=sqrt(s2/(n*1.0d0)) 
      s_likelihood=0.0d0
      do i=1,n
	  s2=0.0d0
	  do j=1,n_ig(i)
	    s1=0.0d0
	    do k=1,5
 	      s1=s1+(1.0d0/sqrt(2.0d0*pi*(var**2)))*exp(-((trait(i)-
     /         u_mean)**2)/(2.0d0*(var**2)))*q_g(i,j,k)*g_o(i,j)
	    end do
	    s2=s2+s1
	  end do
	  s_likelihood=s_likelihood+(log(s2))
      end do
	sLOD=slikelihood-s_likelihood   

      return
      end
C     ****************************************************************
      subroutine qtl_distribution1(fqtl,mqtl,a,fre_qtl,nqtl)
      implicit double precision(a-h,o-z)
      
      integer fqtl(4),mqtl(4)      
      integer g1(10,2),g2(10,2)    
      integer zygote(100,4)        
      integer nqtl(5)             
                                   
      dimension fre_qtl(5)         
      dimension fre1(10),fre2(10)  
      dimension fre(100)           
      
      call gametes1(g1,a,fqtl,fre1) 
      call gametes1(g2,a,mqtl,fre2)     
      do i=1,10
        do j=1,10
          k=(i-1)*10+j
          zygote(k,1)=g1(i,1)
          zygote(k,2)=g1(i,2)
          zygote(k,3)=g2(j,1)
          zygote(k,4)=g2(j,2)
          fre(k)=fre1(i)*fre2(j)
        end do
      end do      
      fre_qtl=0.0d0
      do k=1,100
        iq=0
        do i=1,4
          iq=iq+zygote(k,i)
        end do
        if(iq.eq.4)then       
          fre_qtl(1)=fre_qtl(1)+fre(k)
        else if(iq.eq.5)then  
          fre_qtl(2)=fre_qtl(2)+fre(k)
        else if(iq.eq.6)then   
          fre_qtl(3)=fre_qtl(3)+fre(k)
        else if(iq.eq.7)then   
          fre_qtl(4)=fre_qtl(4)+fre(k)
        else if(iq.eq.8)then  
          fre_qtl(5)=fre_qtl(5)+fre(k)
        end if
      end do    
      nqtl=0
      do i=1,5
        if(fre_qtl(i)>0.0d0)then
          nqtl(i)=1
        end if
      end do
      
      return
      end     
C     ***********************************************************
      subroutine gametes1(g,a,qtl,fre)
      implicit double precision(a-h,o-z)
      
      integer g(10,2),qtl(4)
      
      dimension fre(10)
      
      g(1,1)=qtl(1)
      g(1,2)=qtl(1)
      fre(1)=a/4.0d0          
      
      g(2,1)=qtl(2)
      g(2,2)=qtl(2)
      fre(2)=a/4.0d0         
      
      g(3,1)=qtl(3)
      g(3,2)=qtl(3)
      fre(3)=a/4.0d0          
      
      g(4,1)=qtl(4)
      g(4,2)=qtl(4)
      fre(4)=a/4.0d0         
      
      g(5,1)=qtl(1)
      g(5,2)=qtl(2)
      fre(5)=(1.0d0-a)/6.0d0 
      
      g(6,1)=qtl(1)
      g(6,2)=qtl(3)
      fre(6)=(1.0d0-a)/6.0d0   
      
      g(7,1)=qtl(1)
      g(7,2)=qtl(4)
      fre(7)=(1.0d0-a)/6.0d0   
      
      g(8,1)=qtl(2)
      g(8,2)=qtl(3)
      fre(8)=(1.0d0-a)/6.0d0   
      
      g(9,1)=qtl(2)
      g(9,2)=qtl(4)
      fre(9)=(1.0d0-a)/6.0d0  
      
      g(10,1)=qtl(3)
      g(10,2)=qtl(4)
      fre(10)=(1.0d0-a)/6.0d0  
      
      return
      end  
