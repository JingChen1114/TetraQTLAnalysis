C     Purpose:
C             working out the threshold
C             quadrivalent pairing during meiosis 
C     Records of revision:
C     Date                Programmer
C     03/01/2019          Jing Chen
      subroutine qvthreshold1(i,n,mloci,ftype,mtype,alpha,rec,o
     /           ,ip1,ip2,g_o,n_ig)
      implicit double precision(a-h,o-z)
      
      parameter(mxind=500,mxloci=400,max_g=136*136)
      
      integer i,n,mloci
      integer ftype(mxloci,4),mtype(mxloci,4)  
      integer f_type(3,4),m_type(3,4) 
      integer o(mxind,mxloci,8)
      integer flankmark(2,8),o_ind(mxloci,8)
      integer n_ig(mxind)
      integer ip1(mxind,max_g,2),ip2(mxind,max_g,2)
      
      dimension alpha(mxloci)
      dimension rec(mxloci)
      dimension zygote_g_o(max_g),g_o(mxind,max_g)
      dimension proba_r(max_g),proba_l(max_g)
      
c      i_s1=1
c      ilocus=0
c      rr=0.0d0
c      do i=2,mloci
        do j=1,n
          ii=i-1    
          r=rec(i)  
          alpha1=alpha(ii)   
          imark=ii
          do k=1,4
            f_type(1,k)=ftype(ii,k)  
            f_type(3,k)=ftype(i,k)  
            m_type(1,k)=mtype(ii,k)  
            m_type(3,k)=mtype(i,k)  
          end do 
          do k=1,8  
            flankmark(1,k)=o(j,ii,k)
            flankmark(2,k)=o(j,i,k)
          end do
          do l=1,mloci   
            do k=1,8
              o_ind(l,k)=o(j,l,k)
            end do
          end do
          call pr_chrconfig_markers1(j,mloci,imark,f_type,m_type,rec,
     /         alpha1,r,flankmark,o_ind,i_g,zygote_g_o,proba_r,
     /         proba_l,ftype,mtype,ip1,ip2)
          s=0.0d0
          do m=1,i_g
            s=s+zygote_g_o(m)*proba_l(m)*proba_r(m)
          end do
          do m=1,i_g
            g_o(j,m)=zygote_g_o(m)*proba_l(m)*proba_r(m)/s
	    end do
          n_ig(j)=i_g
	  end do                   
c        i_s=1
          
c        rec1=0.005d0   
c        do k=1,4
c          f_type(2,k)=fqtl(k)
c          m_type(2,k)=mqtl(k)
c        end do
c        it=0
c        do while(rec1<r)
c          it=it+1
c          ilocus=ilocus+1
c          re(ilocus)=rec1+rr-4.0d0*rec1*rr/3.0d0
c          ge_distance(ilocus)=-75.0d0*log(1.0d0-4.0d0*re(ilocus
c     /      )/3.0d0)
c          alpha_t(ilocus)=alpha(ii)
c          do j=1,n
c            i_g=n_ig(j)
c            qtl_g=0.0d0
c            do iC=1,i_g
c              if1=ip1(j,iC,1)
c              if2=ip1(j,iC,2)
c              im1=ip2(j,iC,1)
c              im2=ip2(j,iC,2)
c              call pr_qtl_chrconfig1(f_type,m_type,r,rec1,if1,if2,im1
c     /             ,im2,iC,qtl_g)
c              do k=1,5
c                q_g(j,iC,k)=qtl_g(iC,k)
c              end do 
c            end do
c          end do   
c          s=0.0d0
c          do j=1,n
c            do iC=1,n_ig(j)
c              do k=1,5
c                s=s+q_g(j,iC,k)
c              end do
c            end do
c          end do
c          q=0.0d0
c          do k1=1,5
c            do k2=1,n
c              ig=n_ig(k2)
c              do m=1,ig
c                q(k2,k1)=q(k2,k1)+q_g(k2,m,k1)*g_o(k2,m)
c              end do
c            end do
c          end do 
c          sum=0.0d0
c          do k2=1,n
c            do k1=1,5
c              sum=sum+q(k2,k1)
c            end do
c          end do      
c          rec1=rec1+0.01d0  
c        end do  
c        rr=rec(i)+rr-4.0d0*rec(i)*rr/3.0d0
c	end do     
	
	return
      end  
C     *********************************************************************
      subroutine pr_chrconfig_markers1(individual,mloci,imark,f_type
     /           ,m_type,rec,alpha1,r,flankmark,o_ind,i_g,zyg_g_o,
     /           proba_r,proba_l,ftype,mtype,ip1,ip2)
     
      implicit double precision(a-h,o-z)
      
      parameter(mxind=500,mxloci=400,max_g=136*136)
      
      integer individual,mloci,imark,i_g
      integer f_type(3,4),m_type(3,4) 
      integer flankmark(2,8)   
      integer o_ind(mxloci,8)  
      integer ftype(mxloci,4),mtype(mxloci,4)  
      integer mg1(24,2,2),mg2(24,2,2)
      integer zygote(3,4),zyg_p(2,8)
      integer ip1(mxind,max_g,2),ip2(mxind,max_g,2)
      
      dimension rec(mxloci)    
      dimension zyg_g_o(max_g)
      dimension proba_r(max_g),proba_l(max_g)
      
      i_g=0
      zyg_g_o=0.0d0
      proba_r=0.0d0
      proba_l=0.0d0
      do i=1,11              
        do j=1,11
          call possgamete1(f_type,i,mg1,n_mode1)
          call possgamete1(m_type,j,mg2,n_mode2)
          do m=1,n_mode1
            do n=1,n_mode2
              zygote(1,1)=mg1(m,1,1)
	          zygote(1,2)=mg1(m,1,2)
	          zygote(1,3)=mg2(n,1,1)
	          zygote(1,4)=mg2(n,1,2)
	          zygote(3,1)=mg1(m,2,1)
	          zygote(3,2)=mg1(m,2,2)
	          zygote(3,3)=mg2(n,2,1)
	          zygote(3,4)=mg2(n,2,2)
	          zyg_p=0
	          do k=1,4
	            if(zygote(1,k).gt.0) zyg_p(1,zygote(1,k))=zyg_p
     /              (1,zygote(1,k))+1
	            if(zygote(3,k).gt.0) zyg_p(2,zygote(3,k))=zyg_p
     /              (2,zygote(3,k))+1
	          end do
	          sum=0.0d0
	          do l=1,2
	            do k=1,8
      	          sum=sum+abs(zyg_p(l,k)-flankmark(l,k))
      	        end do
      	      end do
      	      if(sum.eq.0.0d0)then
      	        call proba_g_o1(i,alpha1,r,g_o1)  
      	        call proba_g_o1(j,alpha1,r,g_o2)
      	        call proba_g_o_multilocus1(rec,o_ind,i,j,m,n,imark,
     /               mloci,ftype,mtype,pro_l,pro_r)
                if((pro_r.ne.0.0d0).and.(pro_l.ne.0.0d0).and.
     /            (g_o1.ne.0.0d0).and.(g_o2.ne.0.0d0))then
                  i_g=i_g+1  
                  ip1(individual,i_g,1)=i
                  ip1(individual,i_g,2)=m
                  ip2(individual,i_g,1)=j
                  ip2(individual,i_g,2)=n
                  zyg_g_o(i_g)=g_o1*g_o2
                  proba_r(i_g)=pro_r
                  proba_l(i_g)=pro_l 
                end if
              end if 
            end do
          end do
        end do
      end do
      
      return  
      end
C     *******************************************************
      subroutine possgamete1(mp,mode_i,mg,n_mode)
      implicit double precision(a-h,o-z)
      
      integer mode_i,n_mode,mp(3,4),mg(24,2,2)
      
      mg=0  
	m=1
      if(mode_i==1)then
	    n_mode=4
        do i=1,4
	    mg(i,1,1)=mp(1,i)
          mg(i,1,2)=mp(1,i)
	    mg(i,2,1)=mp(3,i)
	    mg(i,2,2)=mp(3,i)
	  end do

	else if(mode_i==2)then
	  n_mode=12
	  do i=1,4
          do j=1,4
	      if(i.ne.j)then
              mg(m,1,1)=mp(1,i)
	        mg(m,1,2)=mp(1,i)
			  mg(m,2,1)=mp(3,j)
	        mg(m,2,2)=mp(3,j) 
	        m=m+1
	      end if
	    end do
	  end do

	else if(mode_i==3)then
	  n_mode=12
        do i=1,4
	    do j=1,4
		    if(i.ne.j)then
		      mg(m,1,1)=mp(1,i)
  	 		  mg(m,1,2)=mp(1,i)
			  mg(m,2,1)=mp(3,i)
			  mg(m,2,2)=mp(3,j)
	        m=m+1
		    end if
		  end do
	  end do
	   
	else if(mode_i==4)then
	  n_mode=12
	  do i=1,4
	    do j=1,4
		    do k=j+1,4
		      if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
			    mg(m,1,1)=mp(1,i)
			    mg(m,1,2)=mp(1,i)
			    mg(m,2,1)=mp(3,j)
	          mg(m,2,2)=mp(3,k)
	          m=m+1
	        end if
	      end do
	    end do
	  end do

	else if(mode_i==5)then
	  n_mode=12
	  do i=1,4
	    do j=1,4
	      if(i.ne.j)then
	        mg(m,1,1)=mp(1,i)
	        mg(m,1,2)=mp(1,j)
			  mg(m,2,1)=mp(3,i)
	        mg(m,2,2)=mp(3,i)
			  m=m+1
	      end if
	    end do
	  end do

	else if(mode_i==6)then
	  n_mode=12
	  do j=1,4
	    do i=1,4
	      do k=i+1,4
              if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          mg(m,1,1)=mp(1,i)
			    mg(m,1,2)=mp(1,k)
	          mg(m,2,1)=mp(3,j)
	          mg(m,2,2)=mp(3,j)
	          m=m+1
	        end if
	      end do
	    end do
	  end do

	else if(mode_i==7)then
	  n_mode=6
	  do i=1,4
	    do j=i+1,4
	      mg(m,1,1)=mp(1,i)
	      mg(m,1,2)=mp(1,j)
	      mg(m,2,1)=mp(3,i)
	      mg(m,2,2)=mp(3,j)
	      m=m+1
	    end do
	  end do

      else if(mode_i==8)then
	  n_mode=24
	  do i=1,4
	    do j=1,4
	      do k=1,4
	        if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          mg(m,1,1)=mp(1,i)
	          mg(m,1,2)=mp(1,j)
			    mg(m,2,1)=mp(3,i)
	          mg(m,2,2)=mp(3,k)
	          m=m+1
	        end if
	      end do
	    end do
	  end do

	else if(mode_i==9)then
	  n_mode=6
	  do i=1,4
	    do j=i+1,4
	      mg(m,1,1)=mp(1,i)
	      mg(m,1,2)=mp(1,j)
	      mg(m,2,1)=mp(3,j)
		    mg(m,2,2)=mp(3,i)
	      m=m+1
	    end do
	  end do

	else if(mode_i==10)then
	  n_mode=24
	  do i=1,4
	    do j=1,4
	      do k=1,4
	        if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          mg(m,1,1)=mp(1,i)
	          mg(m,1,2)=mp(1,j)
	          mg(m,2,1)=mp(3,j)
	          mg(m,2,2)=mp(3,k)
	          m=m+1
	        end if
		    end do
	    end do
	  end do

      else if(mode_i==11)then
	  n_mode=12
	  do i=1,4
	    do j=1,4
	      do k=i+1,4
	        do l=1,4
	          if((i.ne.j).and.(i.ne.k).and.(i.ne.l).and.(j.ne.k)
     /            .and.(j.ne.l).and.(k.ne.l))then
	            mg(m,1,1)=mp(1,i)
	            mg(m,1,2)=mp(1,k)
	            mg(m,2,1)=mp(3,j)
	            mg(m,2,2)=mp(3,l)
	            m=m+1
	          end if
	        end do
	      end do
	    end do
	  end do
      end if
	
	return     
	end
C     ******************************************************
      subroutine proba_g_o1(ii,a,r,g_o)
      implicit double precision(a-h,o-z)
      
      integer ii
      
      if(ii==1)then
		g_o=a*((1.0d0-r)**2)/4.0d0
		 
	  else if(ii==2)then
	    g_o=a*(r**2)/36.0d0

	  else if(ii==3)then
	    g_o=a*r*(1.0d0-r)/6.0d0
	  
	  else if(ii==4)then
		g_o=a*(r**2)/18.0d0

	  else if(ii==5)then
	    g_o=(1.0d0-a)*r*(1.0d0-r)/18.0d0

	  else if(ii==6)then
	    g_o=(1.0d0-a)*(r**2)/54.0d0

	  else if(ii==7)then
	    g_o=(1.0d0-a)*((1.0d0-r)**2)/6.0d0

	  else if(ii==8)then
	    g_o=(1.0d0-a)*r*(1.0d0-r)/18.0d0

	  else if(ii==9)then
	    g_o=(1.0d0-a)*(r**2)/54.0d0

	  else if(ii==10)then
	    g_o=(1.0d0-a)*(r**2)/54.0d0

	  else if(ii==11)then
	    g_o=(1.0d0-a)*(r**2)/54.0d0
      end if

	return
	end
C     ***************************************************************
      subroutine proba_g_o_multilocus1(rec,o_ind,i1,j1,m1,n1,imark,
     /           mloci,ftype,mtype,pro_l,pro_r)
      implicit double precision(a-h,o-z)
      
      parameter(mxloci=400)
      
      integer o_ind(mxloci,8),ftype(mxloci,4),mtype(mxloci,4)
      
      dimension trans(10,10),trans_t(10,10)
      dimension rec(mxloci)
      dimension zy_matrix_l1(100),zy_matrix_l2(100)
      dimension zy_matrix_r1(100),zy_matrix_r2(100)
      dimension trans_zy(100,100),trans_t_zy(100,100)
      
      call ge_interval1(i1,m1,i_ga_f1,i_ga_f2)      
      call ge_interval1(j1,n1,i_ga_m1,i_ga_m2)      
      zy_matrix_l1=0.0d0
      zy_matrix_r1=0.0d0
      i=(i_ga_f1-1)*10+i_ga_m1
	j=(i_ga_f2-1)*10+i_ga_m2
	zy_matrix_l1(i)=1.0d0
	zy_matrix_r1(j)=1.0d0 

      if(imark.ne.(mloci-1))then 
	  do i=imark+1,mloci-1       
	    ii=i+1
	    r=rec(ii)
	    call tra_matrix1(r,trans) 
	    do k1=1,10
	      do l1=1,10
	        do k2=1,10
	          do l2=1,10
	            jj=(k1-1)*10+l1
	            ll=(k2-1)*10+l2
	            trans_zy(jj,ll)=trans(k1,k2)*trans(l1,l2)
	          end do
              end do
	      end do
          end do
          s=0.0d0
	    do k1=1,100
	      do l1=1,100
	        s=s+trans_zy(k1,l1)
	      end do
          end do
          zy_matrix_r2=0.0d0
          do l=1,100
	      do j=1,100
	        zy_matrix_r2(l)=zy_matrix_r2(l)+trans_zy(j,l)*
     /                        zy_matrix_r1(j)
	      end do
	    end do
	    s=0.0d0
	    do l=1,100
	      s=s+zy_matrix_r2(l)
	    end do
	    call check_ph1(ii,zy_matrix_r2,ftype,mtype,o_ind,
     /         zy_matrix_r1)
	  end do
	  pro_r=0.0d0
	  do i=1,100
	    pro_r=pro_r+zy_matrix_r1(i)
	  end do
	else if(imark.eq.(mloci-1))then
	  pro_r=1.0d0
	end if
      if(imark.ne.1)then
	  do i=imark,2,-1        
	    ii=i-1
	    r=rec(i)
	    call tra_matrix1(r,trans)
	    do l=1,10
	      do j=1,10
	        trans_t(l,j)=trans(j,l)
	      end do
	    end do
	    do k1=1,10
	      do l1=1,10
	        do k2=1,10
                do l2=1,10
	            jj=(k1-1)*10+l1
	            ll=(k2-1)*10+l2
	            trans_t_zy(jj,ll)=trans_t(k1,k2)*trans_t(l1,l2)
	          end do
	        end do
	      end do
	    end do
	    zy_matrix_l2=0.0d0
          do l=1,100
	      do j=1,100
	        zy_matrix_l2(l)=zy_matrix_l2(l)+trans_t_zy(l,j)*
     /        zy_matrix_l1(j)
	      end do
	    end do
          call check_ph1(ii,zy_matrix_l2,ftype,mtype,o_ind,
     /         zy_matrix_l1)
        end do
        pro_l=0.0d0
	  do i=1,100
	    pro_l=pro_l+zy_matrix_l1(i)
	  end do
      else if(imark.eq.1)then
	  pro_l=1.0d0
	end if
	
	return
	end 
C     ***************************************************************
      subroutine ge_interval1(ii,im,i_gamete1,i_gamete2)
      implicit double precision(a-h,o-z)
      integer ii,im,i_gamete1,i_gamete2
      m=0
      if(ii==1)then
        do i=1,4
	    if(i==im)then
	      i_gamete1=i
	      i_gamete2=i
	    end if
	  end do
	else if(ii==2)then
	  do i=1,4
	    do j=1,4
	      if(i.ne.j)then
	        m=m+1
	        if(m==im)then
	          i_gamete1=i
	          i_gamete2=j
	        end if
	      end if
	    end do
	  end do
	else if(ii==3)then
        do i=1,4
	    do j=1,4
		    if(i.ne.j)then
	        m=m+1
	        if(m==im)then
	          i_gamete1=i
	          call num_gamete1(i,j,i_gamete2)
	        end if
		    end if
		  end do
	  end do	  
	else if(ii==4)then
	  do i=1,4
	    do j=1,4
		    do k=j+1,4
		      if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          m=m+1
	          if(m==im)then
	            i_gamete1=i
	            call num_gamete1(j,k,i_gamete2)
	          end if
	        end if
	      end do
	    end do
	  end do
      else if(ii==5)then
	  do i=1,4
	    do j=1,4
	      if(i.ne.j)then
			  m=m+1
	        if(m==im)then
	          call num_gamete1(i,j,i_gamete1)
	          i_gamete2=i
	        end if
	      end if
	    end do
	  end do
	else if(ii==6)then
	  do j=1,4
	    do i=1,4
	      do k=i+1,4
              if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          m=m+1
	          if(m==im)then
	            call num_gamete1(i,k,i_gamete1)
	            i_gamete2=j
	          end if
	        end if
	      end do
	    end do
	  end do
	else if(ii==7)then
	  do i=1,4
	    do j=i+1,4
	      m=m+1
	      if(m==im)then
	        call num_gamete1(i,j,i_gamete1)
	        call num_gamete1(i,j,i_gamete2)
	      end if
	    end do
	  end do
	else if(ii==8)then
	  do i=1,4
	    do j=1,4
	      do k=1,4
	        if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          m=m+1
	          if(m==im)then
	            call num_gamete1(i,j,i_gamete1)
	            call num_gamete1(i,k,i_gamete2)
	          end if
	        end if
	      end do
	    end do
	  end do
	else if(ii==9)then
	  do i=1,4
	    do j=i+1,4
	      m=m+1
	      if(m==im)then
	        call num_gamete1(i,j,i_gamete1)
	        call num_gamete1(j,i,i_gamete2)
	      end if
	    end do
	  end do
	else if(ii==10)then
	  do i=1,4
	    do j=1,4
	      do k=1,4
	        if((i.ne.j).and.(i.ne.k).and.(j.ne.k))then
	          m=m+1
	          if(m==im)then
	            call num_gamete1(i,j,i_gamete1)
	            call num_gamete1(j,k,i_gamete2)
	          end if
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
	            m=m+1
	            if(m==im)then
	              call num_gamete1(i,k,i_gamete1)
	              call num_gamete1(j,l,i_gamete2)
	            end if
	          end if
	        end do
	      end do
	    end do
	  end do
      end if

      return
	end
C     ****************************************************************
      subroutine num_gamete1(i1,i2,num_g)
	implicit double precision(a-h,o-z)

	integer i1,i2,num_g
	integer j1,j2
	
	if(i1<i2)then
	  j1=i1
	  j2=i2
	else if(i1>i2)then
	  j1=i2
	  j2=i1
	end if

	if((j1==1).and.(j2==2))then
	  num_g=5
	else if((j1==1).and.(j2==3))then
	  num_g=6
	else if((j1==1).and.(j2==4))then
	  num_g=7
	else if((j1==2).and.(j2==3))then
	  num_g=8
	else if((j1==2).and.(j2==4))then
	  num_g=9
	else if((j1==3).and.(j2==4))then
        num_g=10
	end if
         
      return
	end 
C     ****************************************************************
      subroutine tra_matrix1(r,trans)
	implicit double precision(a-h,o-z)

	dimension trans(10,10)

	a=(1.0d0-r)**2
	b=r*(1.0d0-r)
	c=r**2
      do i=1,4
	  do j=1,4
	    trans(i,j)=c/9.0d0
	  end do
	end do
	do i=1,4
	  trans(i,i)=a
	end do                       
	do i=5,10
	  do j=5,10
	    trans(i,j)=b/3.0d0+c/9.0d0
	  end do
	end do
	do i=5,10
	  trans(i,i)=a+c/9.0d0
	end do
	trans(10,5)=2.0d0*c/9.0d0
      trans(9,6)=2.0d0*c/9.0d0
	trans(8,7)=2.0d0*c/9.0d0
	trans(7,8)=2.0d0*c/9.0d0
	trans(6,9)=2.0d0*c/9.0d0
	trans(5,10)=2.0d0*c/9.0d0   
      do i=1,4
	  do j=5,10
	    trans(i,j)=2.0d0*c/9.0d0
	  end do
	end do
	trans(1,5)=2.0d0*b/3.0d0
      trans(1,6)=2.0d0*b/3.0d0
	trans(1,7)=2.0d0*b/3.0d0
	trans(2,5)=2.0d0*b/3.0d0
	trans(2,8)=2.0d0*b/3.0d0
	trans(2,9)=2.0d0*b/3.0d0
	trans(3,6)=2.0d0*b/3.0d0
	trans(3,8)=2.0d0*b/3.0d0
	trans(3,10)=2.0d0*b/3.0d0
	trans(4,7)=2.0d0*b/3.0d0
	trans(4,9)=2.0d0*b/3.0d0
	trans(4,10)=2.0d0*b/3.0d0         
	do i=5,10
	  do j=1,4
	    trans(i,j)=trans(j,i)/2.0d0
	  end do
	end do                  
	s=0.0d0
	do i=1,10
	  do j=1,10
	    s=s+trans(i,j)
	  end do
	end do

      return
	end
C     ****************************************************************
      subroutine check_ph1(m,zy_matri,fp,mp,o_ind,genotype)
	implicit double precision(a-h,o-z)
	parameter(mxloci=400)

	integer m,mp(mxloci,4),fp(mxloci,4),o_ind(mxloci,8)
	integer nz(8),n1
	integer gamete1(2),gamete2(2),zygote(4)

	dimension zy_matri(100),genotype(100)

	genotype=0.0d0
	do i=1,10
	  do j=1,10
	    call off1(m,i,fp,gamete1)
	    call off1(m,j,mp,gamete2)
	    do k=1,2
	      k1=k+2
	      zygote(k)=gamete1(k)
	      zygote(k1)=gamete2(k)
	    end do
	    nz=0
	    do k=1,4
	      if(zygote(k)>0)then
	        nz(zygote(k))=nz(zygote(k))+1
	      end if
	    end do
          n1=0
	    do k=1,8
	      n1=n1+abs(nz(k)-o_ind(m,k))
	    end do
	    if(n1==0)then
	      ii=(i-1)*10+j
	      genotype(ii)=zy_matri(ii)
	    end if
	  end do
	end do
	
	return
	end
C     ************************************************************
      subroutine off1(m,i,parent,gamete)
	implicit double precision(a-h,o-z)
	parameter(mxloci=400)

	integer m,i
	integer parent(mxloci,4)
	integer gamete(2)

      if(i==1)then
	  gamete(1)=parent(m,1)
	  gamete(2)=parent(m,1)
	else if(i==2)then
	  gamete(1)=parent(m,2)
	  gamete(2)=parent(m,2)
	else if(i==3)then
	  gamete(1)=parent(m,3)
	  gamete(2)=parent(m,3)
	else if(i==4)then
	  gamete(1)=parent(m,4)
	  gamete(2)=parent(m,4)
	else if(i==5)then
	  gamete(1)=parent(m,1)
	  gamete(2)=parent(m,2)
	else if(i==6)then
	  gamete(1)=parent(m,1)
	  gamete(2)=parent(m,3)
	else if(i==7)then
	  gamete(1)=parent(m,1)
	  gamete(2)=parent(m,4)
	else if(i==8)then
        gamete(1)=parent(m,2)
	  gamete(2)=parent(m,3)
  	else if(i==9)then
	  gamete(1)=parent(m,2)
	  gamete(2)=parent(m,4)
	else if(i==10)then
	  gamete(1)=parent(m,3)
	  gamete(2)=parent(m,4)
	end if
      
      return
	end 
      
