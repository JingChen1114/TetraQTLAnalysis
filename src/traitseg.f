      subroutine traitseg(n,trait,kmean,pa1,pa2,fre,g,e,alpha,df,sLRT)
      implicit double precision(a-h,o-z)
      parameter(maxind=1000,pi=3.1415926d0)
      
      integer df,p1(4),p2(4),ioff(5)
      integer pa1(4),pa2(4)
      
      dimension trait(maxind),fre_off(5)
      dimension g_esti1(5),g_esti2(5),posterior(maxind,5)
      dimension kmean(4,5)
      dimension fre(5),g(5)
      
      sLRT=0.0d0
      do i1=1,12
        do ia=1,50
          call offspring(ia,i1,p1,p2,ioff,fre_off)
          igenotype=0
          do i=1,5
            igenotype=igenotype+ioff(i)
          end do  
          sum=0.0d0        
          do i=1,n
            sum=sum+trait(i)
          end do
          g_esti1(3)=sum/(n*1.0d0)
          g_esti1(1)=g_esti1(3)*0.6d0
          g_esti1(2)=g_esti1(3)*0.8d0
          g_esti1(4)=g_esti1(3)*1.2d0
          g_esti1(5)=g_esti1(3)*1.4d0
          sum1=0.0d0
          do i=1,n
            sum1=sum1+trait(i)/n
          end do
          sum2=0.0d0
          do i=1,n
            sum2=sum2+(trait(i)-sum1)**2
          end do
          var_esti1=sqrt(sum2/(n-1))
          
          slikelihood1=0.0d0
          do i=1,n
            s=0.0d0
            do k=1,5
              s=s+fre_off(k)*exp(-((trait(i)-g_esti1(k))**2)/(2.0d0
     &          *(var_esti1**2)))/(var_esti1*sqrt(2.0d0*pi))
            end do
            slikelihood1=slikelihood1+log(s)
          end do
          
          s=1.0d0
          do while(s>0.000001d0)
            do i=1,5
              do j=1,n
                p=0.0d0
                do k=1,5
                  p=p+((1.0d0/(sqrt(2.0d0*(var_esti1**2)*pi))*(exp(
     &              -((trait(j)-g_esti1(k))**2)/(2.0d0*(var_esti1**2
     &              ))))))*fre_off(k)
                end do
            
                posterior(j,i)=(((1.0d0/(sqrt(2.0d0*pi*(var_esti1**2
     &                          )))*(exp(-((trait(j)-g_esti1(i))**2)/
     &                         (2.0d0 *(var_esti1**2))))))*
     &                         fre_off(i))/p
              end do
            end do
            do i=1,5
              ss=0.0d0
              do j=1,n
                ss=ss+posterior(j,i)*trait(j)
              end do
              tt=0.0d0
              do j=1,n
                tt=tt+posterior(j,i)
              end do
              if(tt==0.0d0)then
                g_esti2(i)=0.0d0
              else
                g_esti2(i)=ss/tt
              end if
            end do
        
            p=0.0d0
            do i=1,5
              do j=1,n
                p=p+posterior(j,i)*((trait(j)-g_esti2(i))**2)
              end do
            end do
            var_esti2=sqrt(p/(n*1.0d0))
        
            slikelihood2=0.0d0
            do i=1,n
              ss=0.0d0
              do k=1,5
                ss=ss+fre_off(k)*exp(-((trait(i)-g_esti2(k))**2)/
     &             (2.0d0*(var_esti2**2)))/(var_esti2*sqrt(2.0d0*pi))
              end do
              slikelihood2=slikelihood2+log(ss)
            end do
        
            s=slikelihood2-slikelihood1
        
            slikelihood1=slikelihood2
        
            do i=1,5
              g_esti1(i)=g_esti2(i)
            end do
            var_esti1=var_esti2  
          end do
          sum=0.0d0        !mean
          do i=1,n
            sum=sum+trait(i)
          end do
          smean=sum/(n*1.0d0)
            
          sum=0.0d0        !variance
          do i=1,n
            sum=sum+(trait(i)-smean)**2
          end do
          svar=sqrt(sum/(n-1))
              
          slikelihood=0.0d0
          do i=1,n
            s=0.0d0
            do k=1,5
              s=s+fre_off(k)*exp(-((trait(i)-smean)**2)/(2.0d0
     &          *(svar**2)))/(svar*sqrt(2.0d0*pi))
            end do
            slikelihood=slikelihood+log(s)
          end do
              
          s=2.0d0*(slikelihood1-slikelihood)
          if(s>sLRT)then
            sLRT=s
            do i=1,4
              pa1(i)=p1(i)
              pa2(i)=p2(i)
            end do
            alpha=(ia-1)*0.005d0
            do i=1,5
              fre(i)=fre_off(i)
              g(i)=g_esti1(i)
              e=var_esti1
            end do
            df=igenotype-1
          end if
        end do
      end do
      
      return
      end
          
          
C     *********************************************************    
      subroutine offspring(ia,i1,p1,p2,ioff,fre_off)
      
      implicit double precision(a-h,o-z)
      
      integer p1(4),p2(4),g1(10,2),g2(10,2)
      integer zygote(100,4),ioff(5)
      
      dimension fre_off(5),freg1(10),freg2(10)
      dimension frezygote(100)
      
      alpha=(ia-1)*0.005d0
      
      call genotype(i1,p1,p2)
      call gametogesis(alpha,p1,g1,freg1)
      call gametogesis(alpha,p2,g2,freg2)
      
      fre_off=0.0d0
      do i=1,10
        do j=1,10
          k=(i-1)*10+j
          zygote(k,1)=g1(i,1)
          zygote(k,2)=g1(i,2)
          zygote(k,3)=g2(j,1)
          zygote(k,4)=g2(j,2)
          frezygote(k)=freg1(i)*freg2(j)
        end do
      end do
      do i=1,100
        isum=0
        do k=1,4
          isum=isum+zygote(i,k)
        end do
        if(isum==4)then
          fre_off(1)=fre_off(1)+frezygote(i)
        else if(isum==5)then
          fre_off(2)=fre_off(2)+frezygote(i)
        else if(isum==6)then
          fre_off(3)=fre_off(3)+frezygote(i)
        else if(isum==7)then
          fre_off(4)=fre_off(4)+frezygote(i)
        else if(isum==8)then
          fre_off(5)=fre_off(5)+frezygote(i)
        end if
      end do
      
      ioff=0
      do i=1,5
        if(fre_off(i)>0.0d0)then
          ioff(i)=1
        end if
      end do
          
      return
      end
      
c     *******************************************************************
      subroutine genotype(i,p1,p2)
      
      implicit double precision(a-h,o-z)
      
      integer p1(4),p2(4)
      
      if(i==1)then
        p1(1)=1
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=1
        p2(3)=1
        p2(4)=1
      else if(i==2)then
        p1(1)=1
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=1
        p2(4)=1
      else if(i==3)then
        p1(1)=1
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      else if(i==4)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=1
        p2(4)=1
      else if(i==5)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      else if(i==6)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=2
      else if(i==7)then
        p1(1)=2
        p1(2)=2
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      else if(i==8)then
        p1(1)=2
        p1(2)=2
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=2
      else if(i==9)then
        p1(1)=2
        p1(2)=2
        p1(3)=2 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=2
      else if(i==10)then
        p1(1)=2
        p1(2)=1
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=1
        p2(3)=1
        p2(4)=1  
      else if(i==11)then
        p1(1)=2
        p1(2)=2
        p1(3)=1 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=1
        p2(4)=1
      else if(i==12)then
        p1(1)=2
        p1(2)=2
        p1(3)=2 
        p1(4)=1
        
        p2(1)=2
        p2(2)=2
        p2(3)=2
        p2(4)=1
      
      end if
      
      return
      end
C     ******************************************************
      subroutine gametogesis(alpha,p,g,freg)
      
      implicit double precision(a-h,o-z)
      
      integer p(4),g(10,2)
      dimension freg(10)
      
      do i=1,4
        g(i,1)=p(i)
        g(i,2)=p(i)
        freg(i)=alpha/4.0d0
      end do
      
      g(5,1)=p(1)
      g(5,2)=p(2)
      freg(5)=(1.0d0-alpha)/6.0d0
      
      g(6,1)=p(1)
      g(6,2)=p(3)
      freg(6)=(1.0d0-alpha)/6.0d0
      
      g(7,1)=p(1)
      g(7,2)=p(4)
      freg(7)=(1.0d0-alpha)/6.0d0
      
      g(8,1)=p(2)
      g(8,2)=p(3)
      freg(8)=(1.0d0-alpha)/6.0d0
      
      g(9,1)=p(2)
      g(9,2)=p(4)
      freg(9)=(1.0d0-alpha)/6.0d0
      
      g(10,1)=p(3)
      g(10,2)=p(4)
      freg(10)=(1.0d0-alpha)/6.0d0
      
      return
      end 

          
         
