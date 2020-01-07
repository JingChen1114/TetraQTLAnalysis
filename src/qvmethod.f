C     Purpose:
C             Interval mapping for QTL in autotetraploids by assuming 
C             quadrivalent pairing during meiosis 
C     Records of revision:
C     Date                Programmer
C     10/10/2018          Jing Chen
      subroutine qvmethod(n,mloci,ftype,mtype,alpha,rec,o,trait,kmeans1
     /           ,fqtl_max,mqtl_max,BIC_min,ProfileLODs,GeneticDistance
     /           ,gDistance_max,alpha_max,u_max,sigma_max)
      implicit double precision(a-h,o-z)
      
      parameter(mxind=500,mxloci=400,max_g=136*136,nqtype=92)
      
      integer n,mloci,kqtl
      integer ftype(mxloci,4),mtype(mxloci,4)   
      integer fqtype(nqtype,4),mqtype(nqtype,4)
      integer f_type(3,4),m_type(3,4) 
      integer o(mxind,mxloci,8)
      integer flankmark(2,8),o_ind(mxloci,8)
      integer n_ig(mxind)
      integer fqtl(4),mqtl(4)
      integer fqtl_max(4),mqtl_max(4),fqtl_max1(4),mqtl_max1(4)
      integer ip1(mxind,max_g,2),ip2(mxind,max_g,2)
      integer ilocus(nqtype),nqtl(5)
      
      dimension alpha(mxloci)
      dimension rec(mxloci),trait(mxind)
      dimension zygote_g_o(max_g),g_o(mxind,max_g)
      dimension q_g(mxind,max_g,5),qtl_g(max_g,5)
      dimension kmeans1(4,5),kmeans(5),u(5)
      dimension proba_r(max_g),proba_l(max_g)
      dimension re(500),ge_distance(500),alpha_t(500)
      dimension u_max(5),rr(92)
      dimension fre_qtl(5)
      dimension BIC_t(nqtype,1000),sLOD_t(nqtype,1000)
      dimension u_t(nqtype,1000,5),residual_t(nqtype,1000)
      dimension ProfileLODs(mxloci),GeneticDistance(mxloci)
      dimension q(mxind,5),u1(5)
      
      call parent_qtl_genotype(fqtype,mqtype)       
      i_s1=1
      ilocus=0
      rr=0.0d0
      do i=2,mloci
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
          call pr_chrconfig_markers(j,mloci,imark,f_type,m_type,rec,
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
        i_s=1
        do im=1,nqtype
          do k=1,4
            fqtl(k)=fqtype(im,k)
            mqtl(k)=mqtype(im,k)
          end do
          a=alpha(ii)
          call qtl_distribution(fqtl,mqtl,a,fre_qtl,nqtl)
          kqtl=0
          do k=1,5
            if(nqtl(k)>0)then
              kqtl=kqtl+1
            end if
          end do
          ik=-1
          do k=1,5
            ik=ik+nqtl(k)
          end do
          i1=1
          kmeans=0.0d0
          do i2=1,5
            if(fre_qtl(i2)>0.0d0)then
              kmeans(i2)=kmeans1(ik,i1)
              i1=i1+1
            end if
          end do
          rec1=0.005d0   
          do k=1,4
            f_type(2,k)=fqtype(im,k)
            m_type(2,k)=mqtype(im,k)
          end do
          it=0
          do while(rec1<r)
            it=it+1
            ilocus(im)=ilocus(im)+1
            re(ilocus(im))=rec1+rr(im)-4.0d0*rec1*rr(im)/3.0d0
            ge_distance(ilocus(im))=-75.0d0*log(1.0d0-4.0d0*re(ilocus
     /      (im))/3.0d0)
            alpha_t(ilocus(im))=alpha(ii)
            do j=1,n
              i_g=n_ig(j)
              qtl_g=0.0d0
              do iC=1,i_g
                if1=ip1(j,iC,1)
                if2=ip1(j,iC,2)
                im1=ip2(j,iC,1)
                im2=ip2(j,iC,2)
                call pr_qtl_chrconfig(f_type,m_type,r,rec1,if1,if2,im1
     /               ,im2,iC,qtl_g)
                do k=1,5
                  q_g(j,iC,k)=qtl_g(iC,k)
                end do 
              end do
            end do   
            s=0.0d0
            do j=1,n
              do iC=1,n_ig(j)
                do k=1,5
                  s=s+q_g(j,iC,k)
               end do
              end do
            end do
            q=0.0d0
            do k1=1,5
              do k2=1,n
                ig=n_ig(k2)
                do m=1,ig
                  q(k2,k1)=q(k2,k1)+q_g(k2,m,k1)*g_o(k2,m)
                end do
              end do
            end do 
            sum=0.0d0
            do k2=1,n
              do k1=1,5
                sum=sum+q(k2,k1)
              end do
            end do      
C           EM ALGORITHM FOR INTERVAL MAPPING 
            call emalgorithm(fqtl,mqtl,n,trait,g_o,q_g,n_ig,kqtl,
     /           kmeans,u,sigma,slikelihood,sLOD)
            BIC=log(n*1.0d0)*(kqtl+1)-2.0d0*slikelihood
            BIC_t(im,ilocus(im))=BIC
            sLOD_t(im,ilocus(im))=sLOD
            do k=1,5
              u_t(im,ilocus(im),k)=u(k)
            end do
            residual_t(im,ilocus(im))=sigma
            rec1=rec1+0.01d0  
          end do  
          rr(im)=rec(i)+rr(im)-4.0d0*rec(i)*rr(im)/3.0d0
	  end do  
	end do    
	i_s=1
      do im=1,nqtype
        do k=1,4
          fqtl(k)=fqtype(im,k)
          mqtl(k)=mqtype(im,k)
        end do
        do i=1,ilocus(im)
          do k=1,5
            u(k)=u_t(im,i,k)
          end do
          icheck=0
	    call check(fqtl,mqtl,u,icheck)
	    if(icheck==1)then
	      if(i_s==1)then
	        BIC_min=BIC_t(im,i)
	        im_max=im
	      else if(i_s>1)then
	        if(BIC_t(im,i)<BIC_min)then
	          BIC_min=BIC_t(im,i)
	          im_max=im
	        end if
	      end if
	      i_s=i_s+1
	    end if
	  end do
	end do
	do k=1,4
	  fqtl_max(k)=fqtype(im_max,k)
	  mqtl_max(k)=mqtype(im_max,k)
	end do
	
	i_s=1
	do i=1,ilocus(im_max)
	  if(i_s==1)then
	    sLOD_max=sLOD_t(im_max,i)
	    gDistance_max=ge_distance(i)
	    alpha_max=alpha_t(i)
	    do k=1,5
	      u_max(k)=u_t(im_max,i,k)
	    end do
	    sigma_max=residual_t(im_max,i)
	    i_s=i_s+1
        else
          if(sLOD_t(im_max,i)>sLOD_max)then
            sLOD_max=sLOD_t(im_max,i)
	      gDistance_max=ge_distance(i)
	      alpha_max=alpha_t(i)
	      do k=1,5
	        u_max(k)=u_t(im_max,i,k)
	      end do
	      sigma_max=residual_t(im_max,i)
	      i_s=i_s+1
	    end if
	  end if
	end do
	GeneticDistance=0.0d0
	ProfileLODs=0.0d0
	do i=1,ilocus(im_max)
	  GeneticDistance(i)=ge_distance(i)
	  ProfileLODs(i)=sLOD_t(im_max,i)
	end do
      
	return
      end  
      
C     **********************************************************
      subroutine parent_qtl_genotype(fqtype,mqtype)
      
      implicit double precision(a-h,o-z)
      
      parameter(nqtype=92)
     
      integer fqtype(nqtype,4),mqtype(nqtype,4)
      integer p1(6,4),p2(6,4)
      
      iqtype=1
      do i=1,4
        do j=i+1,5
          if(i.eq.1)then
            if(j.ne.5)then
              call iqgenotype(i,p1,iqtl1)
              call iqgenotype(j,p2,iqtl2)
              do k=1,iqtl1
                do l=1,iqtl2
                  do m=1,4
                    fqtype(iqtype,m)=p1(k,m)
                    mqtype(iqtype,m)=p2(l,m)
                  end do
                  iqtype=iqtype+1
                end do
              end do
            end if
          else
            call iqgenotype(i,p1,iqtl1)
            call iqgenotype(j,p2,iqtl2)
            do k=1,iqtl1
              do l=1,iqtl2
                do m=1,4
                  fqtype(iqtype,m)=p1(k,m)
                  mqtype(iqtype,m)=p2(l,m)
                end do
                iqtype=iqtype+1
              end do
            end do
          end if
        end do
      end do
      
      return
      end
C     **********************************************************
      subroutine iqgenotype(i,p,iqtl)
      
      implicit double precision(a-h,o-z)
      
      integer i,iqtl
      integer p(6,4)
      
      if(i.eq.1)then      
        iqtl=1
        do k=1,4
          p(1,k)=2
        end do
      else if(i.eq.2)then 
        iqtl=4
        p=2
        p(1,4)=1         
        p(2,3)=1       
        p(3,2)=1         
        p(4,1)=1         
      else if(i.eq.3)then
        iqtl=6
        p=2
        p(1,3)=1
        p(1,4)=1         
        p(2,2)=1
        p(2,4)=1        
        p(3,2)=1
        p(3,3)=1         
        p(4,1)=1
        p(4,4)=1        
        p(5,1)=1
        p(5,3)=1         
        p(6,1)=1
        p(6,2)=1          
      else if(i.eq.4)then 
        iqtl=4
        p=1
        p(1,1)=2         
        p(2,2)=2         
        p(3,3)=2         
        p(4,4)=2         
      else if(i.eq.5)then 
        iqtl=1
        p=1
      end if
      
      return
      end
C     ***************************************************************
      subroutine pr_chrconfig_markers(individual,mloci,imark,f_type,
     /           m_type,rec,alpha1,r,flankmark,o_ind,i_g,zyg_g_o,
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
          call possgamete(f_type,i,mg1,n_mode1)
          call possgamete(m_type,j,mg2,n_mode2)
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
      	        call proba_g_o(i,alpha1,r,g_o1)  
      	        call proba_g_o(j,alpha1,r,g_o2)
      	        call proba_g_o_multilocus(rec,o_ind,i,j,m,n,imark,
     /               mloci,ftype,mtype,pro_l,pro_r)
                if((pro_r.ne.0.0d0).and.(pro_l.ne.0.0d0))then
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
      subroutine possgamete(mp,mode_i,mg,n_mode)
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
      subroutine proba_g_o(ii,a,r,g_o)
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
      subroutine proba_g_o_multilocus(rec,o_ind,i1,j1,m1,n1,imark,
     /           mloci,ftype,mtype,pro_l,pro_r)
      implicit double precision(a-h,o-z)
      
      parameter(mxloci=400)
      
      integer o_ind(mxloci,8),ftype(mxloci,4),mtype(mxloci,4)
      
      dimension trans(10,10),trans_t(10,10)
      dimension rec(mxloci)
      dimension zy_matrix_l1(100),zy_matrix_l2(100)
      dimension zy_matrix_r1(100),zy_matrix_r2(100)
      dimension trans_zy(100,100),trans_t_zy(100,100)
      
      call ge_interval(i1,m1,i_ga_f1,i_ga_f2)      
      call ge_interval(j1,n1,i_ga_m1,i_ga_m2)      
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
	    call tra_matrix(r,trans) 
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
	    call check_ph(ii,zy_matrix_r2,ftype,mtype,o_ind,zy_matrix_r1)
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
	    call tra_matrix(r,trans)
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
          call check_ph(ii,zy_matrix_l2,ftype,mtype,o_ind,zy_matrix_l1)
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
      subroutine ge_interval(ii,im,i_gamete1,i_gamete2)
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
	          call num_gamete(i,j,i_gamete2)
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
	            call num_gamete(j,k,i_gamete2)
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
	          call num_gamete(i,j,i_gamete1)
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
	            call num_gamete(i,k,i_gamete1)
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
	        call num_gamete(i,j,i_gamete1)
	        call num_gamete(i,j,i_gamete2)
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
	            call num_gamete(i,j,i_gamete1)
	            call num_gamete(i,k,i_gamete2)
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
	        call num_gamete(i,j,i_gamete1)
	        call num_gamete(j,i,i_gamete2)
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
	            call num_gamete(i,j,i_gamete1)
	            call num_gamete(j,k,i_gamete2)
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
	          if((i.ne.j).and.(i.ne.k).and.(i.ne.l).and.(j.ne.k).and.
     /            (j.ne.l).and.(k.ne.l))then
	            m=m+1
	            if(m==im)then
	              call num_gamete(i,k,i_gamete1)
	              call num_gamete(j,l,i_gamete2)
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
      subroutine num_gamete(i1,i2,num_g)
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
      subroutine tra_matrix(r,trans)
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
      subroutine check_ph(m,zy_matri,fp,mp,o_ind,genotype)
	implicit double precision(a-h,o-z)
	parameter(mxloci=400)

	integer m,mp(mxloci,4),fp(mxloci,4),o_ind(mxloci,8)
	integer nz(8),n1
	integer gamete1(2),gamete2(2),zygote(4)

	dimension zy_matri(100),genotype(100)

	genotype=0.0d0
	do i=1,10
	  do j=1,10
	    call off(m,i,fp,gamete1)
	    call off(m,j,mp,gamete2)
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
      subroutine off(m,i,parent,gamete)
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
C     ****************************************************************
      subroutine qtl_distribution(fqtl,mqtl,a,fre_qtl,nqtl)
      implicit double precision(a-h,o-z)
      
      integer fqtl(4),mqtl(4)      
      integer g1(10,2),g2(10,2)    
      integer zygote(100,4)        
      integer nqtl(5)             
                                   
      dimension fre_qtl(5)         
      dimension fre1(10),fre2(10)  
      dimension fre(100)           
      
      call gametes(g1,a,fqtl,fre1) 
      call gametes(g2,a,mqtl,fre2)     
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
      subroutine gametes(g,a,qtl,fre)
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
C     ****************************************************************
      subroutine pr_qtl_chrconfig(f_type,m_type,r,rec1,if1,if2,im1,im2,
     /           iC,qtl_ge)
      parameter(mxind=500,max_g=136*136,max_q=20)
      
      implicit double precision(a-h,o-z)
      
      integer f_type(3,4),m_type(3,4) 
      integer mq1(max_q,2),mq2(max_q,2)
      
      dimension q_g1(max_q),q_g2(max_q),qtl_ge(max_g,5)
       
      call proba_q_g(if1,if2,f_type,r,rec1,mq1,iq1,q_g1)
      call proba_q_g(im1,im2,m_type,r,rec1,mq2,iq2,q_g2)
      do k1=1,iq1
        do k2=1,iq2
          n_allele=mq1(k1,1)+mq1(k1,2)+mq2(k2,1)+mq2(k2,2)
          if(n_allele==4)then
            qtl_ge(iC,1)=qtl_ge(iC,1)+q_g1(k1)*q_g2(k2)
          else if(n_allele==5)then
	      qtl_ge(iC,2)=qtl_ge(iC,2)+q_g1(k1)*q_g2(k2)
	    else if(n_allele==6)then
	      qtl_ge(iC,3)=qtl_ge(iC,3)+q_g1(k1)*q_g2(k2)
	    else if(n_allele==7)then
	      qtl_ge(iC,4)=qtl_ge(iC,4)+q_g1(k1)*q_g2(k2)
	    else if(n_allele==8)then
	      qtl_ge(iC,5)=qtl_ge(iC,5)+q_g1(k1)*q_g2(k2)
	    end if
	  end do
	end do
	
	return
	end
C     ***********************************************************
      subroutine proba_q_g(ii,mm,mp,r,r1,mq,q,q_g)
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
	              if((k.ne.l).and.(k.ne.i).and.(k.ne.j).and.(l.ne.i)
     /                .and.(l.ne.j))then
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
	              if((k.ne.i).and.(k.ne.j).and.(k.ne.l).and.(l.ne.i)
     /                .and.(l.ne.j))then
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
	              if((k.ne.i).and.(k.ne.j).and.(k.ne.l).and.(l.ne.i)
     /                .and.(l.ne.j))then
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
	        q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(9.0d0*((1.0d0-r)**2))
		      if(q_g(q)==0.0d0)then
	          q=q-1
	        end if
	        q=q+1
	        mq(q,1)=mp(2,i)
	        mq(q,2)=mp(2,i)
	        q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/(9.0d0*((1.0d0-r)**2))
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
	          if((i.ne.j).and.(i.ne.k).and.(i.ne.l).and.(j.ne.k).and.
     /            (j.ne.l).and.(k.ne.l))then
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
C     ***************************************************************
      subroutine emalgorithm(fqtl,mqtl,n,trait,g_o,q_g,n_ig,kqtl,
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
      subroutine check(fqtl,mqtl,u,icheck)
      implicit double precision(a-h,o-z)
      
      integer fqtl(4),mqtl(4)
      
      dimension u(5)
      
      call number_parent(fqtl,ip1)
      call number_parent(mqtl,ip2)
      
      if((u(ip1).eq.0.0d0).or.(u(ip2).eq.0.0d0))then
        icheck=1
      else 
        if(u(ip1)>u(ip2))then
          icheck=1
        else
          icheck=0
        end if
      end if
      
      return
      end
      
C     ****************************************************************
      subroutine number_parent(pqtl,ip)
      
      implicit double precision(a-h,o-z)
      integer ip,pqtl(4)
      
      isum=0
      do i=1,4
        isum=isum+pqtl(i)
      end do
      if(isum.eq.4)then
        ip=1
      else if(isum.eq.5)then
        ip=2
      else if(isum.eq.6)then
        ip=3
      else if(isum.eq.7)then
        ip=4
      else if(isum.eq.8)then
        ip=5
      end if
      
      return
      end  
