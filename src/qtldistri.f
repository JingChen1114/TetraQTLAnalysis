
C     ****************************************************************
      subroutine qtldistri(fqtl,mqtl,a,fre_qtl)
      implicit double precision(a-h,o-z)
      
      integer fqtl(4),mqtl(4)      
      integer g1(10,2),g2(10,2)    
      integer zygote(100,4)        
      integer nqtl(5)             
                                   
      dimension fre_qtl(5)         
      dimension fre1(10),fre2(10)  
      dimension fre(100)           
      
      call ga(g1,a,fqtl,fre1) 
      call ga(g2,a,mqtl,fre2)     
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
      
      return
      end     
C     ***********************************************************
      subroutine ga(g,a,qtl,fre)
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
