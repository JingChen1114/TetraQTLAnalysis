C     Purpose:
C             working out the threshold
C             quadrivalent pairing during meiosis 
C     Records of revision:
C     Date                Programmer
C     03/01/2019          Jing Chen
      subroutine qvthreshold2(f_type,m_type,r,rec1,if1,if2,im1,
     /           im2,iC,qtl_ge)
      parameter(max_q=20)
      
      implicit double precision(a-h,o-z)
      
      integer iC,if1,if2,im1,im2
      integer f_type(3,4),m_type(3,4) 
      integer mq1(max_q,2),mq2(max_q,2)
      
      dimension q_g1(max_q),q_g2(max_q),qtl_ge(5)
       
      call proba_q_g1(if1,if2,f_type,r,rec1,mq1,iq1,q_g1)
      call proba_q_g1(im1,im2,m_type,r,rec1,mq2,iq2,q_g2)
      do k1=1,iq1
        do k2=1,iq2
          n_allele=mq1(k1,1)+mq1(k1,2)+mq2(k2,1)+mq2(k2,2)
          if(n_allele==4)then
            qtl_ge(1)=qtl_ge(1)+q_g1(k1)*q_g2(k2)
          else if(n_allele==5)then
	      qtl_ge(2)=qtl_ge(2)+q_g1(k1)*q_g2(k2)
	    else if(n_allele==6)then
	      qtl_ge(3)=qtl_ge(3)+q_g1(k1)*q_g2(k2)
	    else if(n_allele==7)then
	      qtl_ge(4)=qtl_ge(4)+q_g1(k1)*q_g2(k2)
	    else if(n_allele==8)then
	      qtl_ge(5)=qtl_ge(5)+q_g1(k1)*q_g2(k2)
	    end if
	  end do
	end do
	
	return
	end
C     ***********************************************************
      subroutine proba_q_g1(ii,mm,mp,r,r1,mq,q,q_g)
      implicit double precision(a-h,o-z)
      
      parameter(max_q=20)
      
      integer q,mp(3,4),mq(max_q,2)
      
      dimension q_g(max_q)
      
      q=0
      mq=0
      m=1
      r2=(r-r1)/(1.0d0-4.0d0*r1/3.0d0)
      q_g=0.0d0

      if(ii==1)then
        do i=1,4
          if(i==mm)then
	      q=q+1
	      mq(q,1)=mp(2,i)
	      mq(q,2)=mp(2,i)
		    q_g(q)=((1.0d0-r1)*(1.0d0-r2)/(1.0d0-r))**2          
		  
	      if(q_g(q)==0.0d0)then
	        q=q-1
	      end if
	  
		    do j=1,4
	        if(j.ne.i)then
	          q=q+1
	          mq(q,1)=mp(2,i)
	          mq(q,2)=mp(2,j)
	          q_g(q)=2.0d0*r1*(1.0d0-r1)*r2*(1.0d0-r2)/(9.0d0*
     /          ((1.0d0-r)**2))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	        end if	
	      end do
                                                
		    do j=1,4
		      if(j.ne.i)then
			    q=q+1
			    mq(q,1)=mp(2,j)
	          mq(q,2)=mp(2,j)
			    q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
			  end if
	      end do
	      do j=1,4
	        do k=1,4
	          if((j.ne.i).and.(j.ne.k).and.(i.ne.k))then
	            q=q+1
				  mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	          end if
	        end do
	      end do                                                 
		  end if
	    m=m+1
	  end do

	else if(ii==2)then
	  do i=1,4
	    do j=1,4
	      if(i.ne.j)then
	        if(m==mm)then
	          q=q+1
	          mq(q,1)=mp(2,j)
	          mq(q,2)=mp(2,j)
			    q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	          q=q+1
	          mq(q,1)=mp(2,i)
	          mq(q,2)=mp(2,j)
	          q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	          q=q+1
	          mq(q,1)=mp(2,j)
	          mq(q,2)=mp(2,i)
	          q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
		  	    q=q+1
			    mq(q,1)=mp(2,i)
			    mq(q,2)=mp(2,i)
			    q_g(q)=((1.0d0-r1)**2)*(r2**2)/(r**2)
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if                                      
			    do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,i)
				    mq(q,2)=mp(2,k)
	              q_g(q)=r1*(1.0d0-r1)*(r2**2)/(3.0d0*(r**2))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
	          end do
			    do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,k)
				    mq(q,2)=mp(2,i)
	              q_g(q)=r1*(1.0d0-r1)*(r2**2)/(3.0d0*(r**2))
  		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
	          end do                                                  
			    do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,k)
				    mq(q,2)=mp(2,j)
	              q_g(q)=r2*(1.0d0-r2)*(r1**2)/(3.0d0*(r**2))
     		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
	          end do
			    do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,j)
				    mq(q,2)=mp(2,k)
	              q_g(q)=r2*(1.0d0-r2)*(r1**2)/(3.0d0*(r**2))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
	          end do			  			                          
			    do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
				    mq(q,1)=mp(2,k)
				    mq(q,2)=mp(2,k)
				    q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
	          end do
			    do k=1,4
	            do l=1,4
	              if((k.ne.l).and.(k.ne.i).and.(k.ne.j).and.
     /                (l.ne.i).and.(l.ne.j))then
					  q=q+1
				   	  mq(q,1)=mp(2,k)
					  mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
				    end if
				  end do
	          end do
			  end if                                                    
              m=m+1
	      end if
	    end do
	  end do

	else if(ii==3)then
        do i=1,4
	    do j=1,4
		    if(i.ne.j)then
			  if(m==mm)then
	          q=q+1
	          mq(q,1)=mp(2,i)
			    mq(q,2)=mp(2,j)
	          q_g(q)=r1*(1.0d0-r1)*((1.0d0-r2)**2)/(r*(1.0d0-r))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	          q=q+1
	          mq(q,1)=mp(2,i)
	          mq(q,2)=mp(2,i)
			    q_g(q)=((1.0d0-r1)**2)*r2*(1.0d0-r2)/(r*(1.0d0-r))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if                                                   
	          do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
				    mq(q,1)=mp(2,i) 
				    mq(q,2)=mp(2,k)
				    q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(3.0d0*r*
     /                     (1.0d0-r))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
                  end if
	          end do                                                  
	          q=q+1
			    mq(q,1)=mp(2,j)
			    mq(q,2)=mp(2,j)
			    q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
			    q=q+1
	          mq(q,1)=mp(2,j)
	          mq(q,2)=mp(2,i)
	          q_g(q)=r1*(1.0d0-r1)*(r2**2)/(9.0d0*r*(1.0d0-r))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	          do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
				    q=q+1
	              mq(q,1)=mp(2,k)
                    mq(q,2)=mp(2,j)
	              q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
                    q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,i)
	              q_g(q)=r1*(1.0d0-r1)*(r2**2)/(9.0d0*r*(1.0d0-r))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	            end if
	          end do                                                  
                do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,k)
	              q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	              q=q+1
	              mq(q,1)=mp(2,j)
	              mq(q,2)=mp(2,k)
	              q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	            end if
	          end do
	          do k=1,4
	            do l=1,4
	              if((k.ne.i).and.(k.ne.j).and.(k.ne.l).and.
     /                (l.ne.i).and.(l.ne.j))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,k)
	                q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do
	          end do                                                  
	        end if
	        m=m+1
		    end if
		  end do
	  end do

      else if(ii==4)then
	  do i=1,4
	    do j=1,4
		    do k=j+1,4
		      if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
			    if(m==mm)then
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,i)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,k)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,i)
	            q_g(q)=((1.0d0-r1)**2)*(r2**2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if                                             
	            q=q+1
                  mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  q=q+1
                  mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
	            if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
                      mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,k)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
                      mq(q,1)=mp(2,j)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                                
	            q=q+1
                  mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
                  if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
                  mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
                      mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,i)
	                q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
		              q=q+1
                      mq(q,1)=mp(2,i)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                             
		          q=q+1
                  mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
                      mq(q,1)=mp(2,k)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
                      mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,j)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
                      mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                              
			    end if 
	          m=m+1
	        end if
	      end do
	    end do
	  end do
	  
	else if(ii==5)then
	  n_mode=12
	  do i=1,4
	    do j=1,4
	      if(i.ne.j)then
	        if(m==mm)then
	          q=q+1
	          mq(q,1)=mp(2,i)
	          mq(q,2)=mp(2,i)
	          q_g(q)=r1*(1.0d0-r1)*((1.0d0-r2)**2)/(r*(1.0d0-r))
			    if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	          q=q+1
	          mq(q,1)=mp(2,i)
	          mq(q,2)=mp(2,j)
	          q_g(q)=r2*(1.0d0-r2)*((1.0d0-r1)**2)/(r*(1.0d0-r))
			    if(q_g(q)==0.0d0)then
	            q=q-1
	          end if       
			    do k=1,4
	            if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,i)
	              mq(q,2)=mp(2,k)
	              q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(3.0d0*r*
     /                     (1.0d0-r))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	            end if
	          end do                                                 
	          q=q+1
	          mq(q,1)=mp(2,j)
	          mq(q,2)=mp(2,i)
	          q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
		        if(q_g(q)==0.0d0)then
	            q=q-1
	          end if
	          q=q+1
	          mq(q,1)=mp(2,j)
	          mq(q,2)=mp(2,j)
	          q_g(q)=(r2**2)*r1*(1.0d0-r1)/(9.0d0*r*(1.0d0-r))
			    if(q_g(q)==0.0d0)then
	            q=q-1
	          end if			  
			    do k=1,4
				  if((k.ne.i).and.(k.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,i)
	              q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	              q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,j)
	              q_g(q)=(r2**2)*r1*(1.0d0-r1)/(9.0d0*r*(1.0d0-r))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	    		  end if
			    end do                                                  
			    do k=1,4
			      if((k.ne.i).and.(k.ne.j))then
				    q=q+1
	              mq(q,1)=mp(2,j)
	              mq(q,2)=mp(2,k)
	              q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
				    if(q_g(q)==0.0d0)then
	                q=q-1
	              end if            
				    q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,k)
	              q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
				    if(q_g(q)==0.0d0)then
	                q=q-1
	              end if 
	            end if
	          end do
	          do k=1,4
	            do l=1,4
	              if((k.ne.i).and.(k.ne.j).and.(k.ne.l).and.
     /                (l.ne.i).and.(l.ne.j))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,k)
	                q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
				    end if
				  end do
			    end do                                             
	        end if
			  m=m+1
	      end if
	    end do
	  end do

	else if(ii==6)then
	  do j=1,4
	    do i=1,4
	      do k=i+1,4
              if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          if(m==mm)then
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,j)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
		          q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,k)
	            q_g(q)=((1.0d0-r1)**2)*(r2**2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if				                   
		          q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
		          q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	  		      q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,j)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
			          q=q+1
	                mq(q,1)=mp(2,j)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do
			      q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,k)
	            q_g(q)=r1*(1.0d0-r1)*(r2**2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
			      q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,i)
	            q_g(q)=r1*(1.0d0-r1)*(r2**2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
				  do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,k)
	                q_g(q)=r1*(1.0d0-r1)*(r2**2)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
			          q=q+1
	                mq(q,1)=mp(2,i)
	                mq(q,2)=mp(2,l)
	                q_g(q)=r1*(1.0d0-r1)*(r2**2)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                                  
                  q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,k)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,i)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                                  
	          end if
	          m=m+1
	        end if
	      end do
	    end do
	  end do

	else if(ii==7)then
	  do i=1,4
	    do j=i+1,4
	      if(m==mm)then
	        q=q+1
	        mq(q,1)=mp(2,i)
	        mq(q,2)=mp(2,j)
	        q_g(q)=((1.0d0-r1)**2)*((1.0d0-r2)**2)/((1.0d0-r)**2)         
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if		                                           	 
	        q=q+1
	        mq(q,1)=mp(2,j)
	        mq(q,2)=mp(2,j)
	        q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/
     /               (9.0d0*((1.0d0-r)**2))
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
	        q=q+1
	        mq(q,1)=mp(2,i)
	        mq(q,2)=mp(2,i)
	        q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/
     /               (9.0d0*((1.0d0-r)**2))
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
              do k=1,4
	          if((k.ne.i).and.(k.ne.j))then
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,j)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/
     /                   (9.0d0*((1.0d0-r)**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,k)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/
     /                   (9.0d0*((1.0d0-r)**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                end if
			  end do                                                        
              q=q+1
	        mq(q,1)=mp(2,j)
	        mq(q,2)=mp(2,i)
	        q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
	        do k=1,4
	          if((k.ne.i).and.(k.ne.j))then
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if				 
			    end if
			  end do
			  do k=1,4
			    do l=1,4
			      if((k.ne.i).and.(k.ne.j).and.(k.ne.l).and.(l.ne.i).
     /              and.(l.ne.j))then
                    q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,l)
	              q_g(q)=(r1**2)*(r2**2)/(81.0d0*((1.0d0-r)**2))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
			    end do
			  end do                                                            
	      end if
	      m=m+1
	    end do
	  end do

	else if(ii==8)then
	  do i=1,4
	    do j=1,4
	      do k=1,4
	        if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          if(m==mm)then
                  q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,k)
	            q_g(q)=r1*(1.0d0-r1)*((1.0d0-r2)**2)/(r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if					            
                  q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,j)
	            q_g(q)=r2*(1.0d0-r2)*((1.0d0-r1)**2)/(r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if                          					  
                  q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,i)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(3.0d0*r*
     /                   (1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
                      q=q+1
	                mq(q,1)=mp(2,i)
	                mq(q,2)=mp(2,l)
	                q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(3.0d0*r*
     /                       (1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                                 
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,k)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(9.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do
		          q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(9.0d0*r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
		          q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(9.0d0*r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
		              q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,j)
	                q_g(q)=(r2**2)*r1*(1.0d0-r1)/(9.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                                 
			      q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
			      q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,i)
	                q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
			          q=q+1
	                mq(q,1)=mp(2,j)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
				      q=q+1
	                mq(q,1)=mp(2,k)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
				      q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(27.0d0*r*(1.0d0-r))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                               
	          end if
	          m=m+1
	        end if
	      end do
	    end do
	  end do

	else if(ii==9)then
	  do i=1,4
	    do j=i+1,4
	      if(m==mm)then
	        q=q+1
	        mq(q,1)=mp(2,j)
	        mq(q,2)=mp(2,i)
	        q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
	        q=q+1
	        mq(q,1)=mp(2,j)
	        mq(q,2)=mp(2,j)
	        q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
		      q=q+1
	        mq(q,1)=mp(2,i)
	        mq(q,2)=mp(2,i)
	        q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
		      q=q+1
	        mq(q,1)=mp(2,i)
	        mq(q,2)=mp(2,j)
	        q_g(q)=((1.0d0-r1)**2)*(r2**2)/(r**2)
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if			                        
			  do k=1,4
	          if((k.ne.i).and.(k.ne.j))then
			      q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
		  	      q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
			      q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
			      q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	          end if
	        end do                                                        
              do k=1,4
	          if((k.ne.i).and.(k.ne.j))then
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
	            if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	          end if
	        end do
	        do k=1,4
	          do l=1,4
	            if((k.ne.i).and.(k.ne.j).and.(k.ne.l).and.(l.ne.i)
     /              .and.(l.ne.j))then
	              q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,l)
	              q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	            end if
	          end do
	        end do                                                       
		    end if
	      m=m+1
	    end do
	  end do

	else if(ii==10)then
	  do i=1,4
	    do j=1,4
	      do k=1,4
	        if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          if(m==mm)then
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,j)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,k)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,j)
	            q_g(q)=((1.0d0-r1)**2)*(r2**2)/(r**2)
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if                                           
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,k)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,j)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,k)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
                      q=q+1
	                mq(q,1)=mp(2,j)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
		              if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,j)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1
	            mq(q,1)=mp(2,i)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
		          if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  do l=1,4
	              if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,j)
	                q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
  	                q=q+1
	                mq(q,1)=mp(2,i)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do												  
	            q=q+1
	            mq(q,1)=mp(2,k)
	            mq(q,2)=mp(2,i)
	            q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			      if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
				do l=1,4
				  if((l.ne.i).and.(l.ne.j).and.(l.ne.k))then
	                q=q+1
	                mq(q,1)=mp(2,k)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,i)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1
	                mq(q,1)=mp(2,l)
	                mq(q,2)=mp(2,l)
	                q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			          if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do                                                
	          end if
	          m=m+1
	        end if
		  end do
	    end do
	  end do

	else if(ii==11)then
	  do i=1,4
	    do j=1,4
	      do k=i+1,4
	        do l=1,4
	          if((i.ne.j).and.(i.ne.k).and.(i.ne.l).and.(j.ne.k)
     /            .and.(j.ne.l).and.(k.ne.l))then
	            if(m==mm)then
                    q=q+1
	              mq(q,1)=mp(2,j)
	              mq(q,2)=mp(2,l)
	              q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
                    q=q+1
	              mq(q,1)=mp(2,j)
	              mq(q,2)=mp(2,k)
	              q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
                    q=q+1
	              mq(q,1)=mp(2,i)
	              mq(q,2)=mp(2,l)
	              q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(r**2)
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
                    q=q+1
	              mq(q,1)=mp(2,i)
	              mq(q,2)=mp(2,k)
	              q_g(q)=((1.0d0-r1)**2)*(r2**2)/(r**2)
				    if(q_g(q)==0.0d0)then
	                q=q-1
	              end if                                                
				    q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,l)
	              q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,l)
	              mq(q,2)=mp(2,l)
	              q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,j)
	              mq(q,2)=mp(2,i)
	              q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,j)
	              mq(q,2)=mp(2,j)
	              q_g(q)=(r1**2)*r2*(1.0d0-r2)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,l)
	              mq(q,2)=mp(2,k)
	              q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,k)
	              q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,i)
	              mq(q,2)=mp(2,j)
	              q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,i)
	              mq(q,2)=mp(2,i)
	              q_g(q)=(r2**2)*r1*(1.0d0-r1)/(3.0d0*(r**2))
				    if(q_g(q)==0.0d0)then
	                q=q-1
	              end if                                                  
				    q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,i)
	              q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
			        if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
	              q=q+1
	              mq(q,1)=mp(2,k)
	              mq(q,2)=mp(2,j)
	              q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
	              if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,l)
	              mq(q,2)=mp(2,i)
	              q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))
		            if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				    q=q+1
	              mq(q,1)=mp(2,l)
	              mq(q,2)=mp(2,j)
	              q_g(q)=(r1**2)*(r2**2)/(9.0d0*(r**2))              
	      	    if(q_g(q)==0.0d0)then
	                q=q-1
	              end if
				  end if
	            m=m+1
	          end if
	        end do
	      end do
	    end do
	  end do
      end if

      return
	end
