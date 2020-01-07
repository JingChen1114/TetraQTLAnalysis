C     Purpose:
C             working out the threshold
C             bivalent pairing during meiosis 
C 
C     Records of revision:
C     Date                Programmer
C     03/01/2019          Jing Chen
      subroutine bvthreshold(n,mloci,ftype,mtype,rec,o,y,
     /           fqtl,mqtl,ir,sLOD_max,ilocus)
      implicit double precision(a-h,o-z)
     
      parameter(iRS=1000,mxind=500,mxloci=400,max_g=48*48)
      
      integer n,mloci,inum
      integer ftype(mxloci,4),mtype(mxloci,4)   
      integer f_type(3,4),m_type(3,4) 
      integer o(mxind,mxloci,8)
      integer flankmark(2,8),o_ind(mxloci,8)
      integer n_ig(mxind,3,3),i_g(3,3)
      integer fqtl(4),mqtl(4)
      integer fqtl_max(4),mqtl_max(4),fqtl_max1(4),mqtl_max1(4)
      integer nqtl(5)
      
      dimension rec(mxloci),y(iRS,mxind),trait(mxind)
      dimension zygote_g_o(3,3,max_g),g_o(mxind,3,3,max_g)
      dimension q_g(mxind,3,3,max_g,5),qtl_g(3,3,max_g,5)
      dimension u(5),kmeans(5)
      dimension proba_r(3,3,max_g),proba_l(3,3,max_g)
      dimension re(500),ge_distance(500),alpha_t(500)
      dimension fre_qtl(5)
      dimension sLOD_t(iRS,1000),sLOD_max(iRS)
            
      i_s1=1
      ilocus=0
      rr=0.0d0
      call qtl_distribution_b1(fqtl,mqtl,fre_qtl,nqtl)
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
      do i=2,mloci
        rec1=0.005d0
        r=rec(i)     
        ii=i-1
        imark=ii 
        do k=1,4
          f_type(1,k)=ftype(ii,k)
          f_type(2,k)=fqtl(k)
          f_type(3,k)=ftype(i,k)
          m_type(1,k)=mtype(ii,k)
          m_type(2,k)=mqtl(k)
          m_type(3,k)=mtype(i,k)
        end do   
        it=0
        do while(rec1<r)
          it=it+1
          ilocus=ilocus+1
          re(ilocus)=rec1+rr-2.0d0*rec1*rr
          ge_distance(ilocus)=-50.0d0*log(1.0d0-2.0d0*
     /                            re(ilocus))    
          do j=1,n
            do k=1,8  
              flankmark(1,k)=o(j,ii,k)
              flankmark(2,k)=o(j,i,k)
            end do
            do l=1,mloci   
              do k=1,8
                o_ind(l,k)=o(j,l,k)
              end do
            end do
            call  posszygote_b1(mloci,imark,f_type,m_type,rec,r,rec1
     /            ,flankmark,o_ind,i_g,zygote_g_o,proba_r,proba_l
     /            ,qtl_g,ftype,mtype)
            s=0.0d0
            do ip1=1,3
              do ip2=1,3 
                ig=i_g(ip1,ip2)
                do m=1,ig
                  s=s+zygote_g_o(ip1,ip2,m)*proba_l
     /              (ip1,ip2,m)*proba_r(ip1,ip2,m)
                end do
              end do
            end do
            do ip1=1,3
              do ip2=1,3
                ig=i_g(ip1,ip2)
                do m=1,ig
                  g_o(j,ip1,ip2,m)=zygote_g_o(ip1,ip2,m)*proba_l
     /                          (ip1,ip2,m)*proba_r(ip1,ip2,m)/s
                  do k=1,5
                    q_g(j,ip1,ip2,m,k)=qtl_g(ip1,ip2,m,k)
	            end do
	          end do
	          n_ig(j,ip1,ip2)=ig
	        end do
	      end do     
	    end do 
         
C         EM ALGORITHM FOR INTERVAL MAPPING 
          do k=1,ir
            do l=1,n
              trait(l)=y(k,l)
            end do   
            call emalgorithm_b1(fqtl,mqtl,n,trait,g_o,q_g,n_ig,kqtl,
     /      kmeans,u,sigma,slikelihood,sLOD)
            
            sLOD_t(k,ilocus)=sLOD
          end do
          rec1=rec1+0.01d0  
        end do  
        rr=rec(i)+rr-2.0d0*rec(i)*rr
      end do    
	
	sLOD_max=0.0d0
	do k=1,ir
	  smax=0.0d0
	  do i=1,ilocus
	    if(sLOD_t(k,i)>smax)then
	      smax=sLOD_t(k,i)
	    end if
	  end do
	  sLOD_max(k)=smax
	end do
	
      return
      end      
C     ***************************************************************
      subroutine posszygote_b1(mloci,imark,mp1,mp2,rec,r,rec1,flankmark,
     /                      o_ind,i_g,zyg_g_o,proba_r,proba_l,qtl_ge,
     /                      ftype,mtype)
     
      implicit double precision(a-h,o-z)
      
      parameter(mxloci=400,max_g=48*48)
      
      integer mp1(3,4),mp2(3,4)  
      integer flankmark(2,8)     
      integer o_ind(mxloci,8)    
      integer ftype(mxloci,4),mtype(mxloci,4)  
      integer pair11(3,2),pair12(3,2),pair21(3,2),pair22(3,2)
      integer mg1(4,2,2),mg2(4,2,2),zygote(3,4),zyg_p(2,8)
      integer i_g(3,3),mq1(4,2),mq2(4,2)
      
      dimension rec(mxloci)     
      dimension zyg_g_o(3,3,max_g)
      dimension proba_r(3,3,max_g),proba_l(3,3,max_g)
      dimension q_g1(4),q_g2(4),qtl_ge(3,3,max_g,5)
      
      zyg_g_o=0.0d0
      proba_r=0.0d0
      proba_l=0.0d0
      qtl_ge=0.0d0
      i_g=0
      do i=1,3 
        do j=1,3
          ig=0
          call pairing1(i,mp1,pair11,pair12)   
          call pairing1(j,mp2,pair21,pair22)   
          do i1=1,4 
            do j1=1,4
              call possgamete_b1(pair11,pair12,i1,mg1,n_mode1)
              call possgamete_b1(pair21,pair22,j1,mg2,n_mode2)
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
                    if(zygote(1,k).gt.0) zyg_p(1,zygote(1,k))=
     /                                   zyg_p(1,zygote(1,k))+1
                    if(zygote(3,k).gt.0) zyg_p(2,zygote(3,k))=
     /                                   zyg_p(2,zygote(3,k))+1
                  end do                 
                  sum=0.0d0
                  do l=1,2  
                    do k=1,8 
                      sum=sum+abs(zyg_p(l,k)-flankmark(l,k))
                    end do
                  end do
                  if(sum.eq.0.0d0)then    
                    call proba_g_o_b1(i1,r,g_o1)  
                    call proba_g_o_b1(j1,r,g_o2)
                    call proba_g_o_multilocus_b1(rec,o_ind,i,j,i1,j1,
     /                   m,n,imark,mloci,ftype,mtype,pro_l,pro_r)
                    if((pro_r.ne.0.0d0).and.(pro_l.ne.0.0d0))then
                      ig=ig+1  
                      zyg_g_o(i,j,ig)=g_o1*g_o2
                      proba_r(i,j,ig)=pro_r
                      proba_l(i,j,ig)=pro_l
                      call proba_q_g_b1(i1,m,pair11,pair12,r,rec1,mq1,
     /                     iq1,q_g1)  
                      call proba_q_g_b1(j1,n,pair21,pair22,r,rec1,mq2,
     /                     iq2,q_g2) 
                      do k1=1,iq1
                        do k2=1,iq2
                          n_allele=mq1(k1,1)+mq1(k1,2)+mq2(k2,1)+mq2
     /                             (k2,2)
                          if(n_allele==4)then
                            qtl_ge(i,j,ig,1)=qtl_ge(i,j,ig,1)+q_g1(k1)*
     /                                       q_g2(k2)
                          else if(n_allele==5)then
	                      qtl_ge(i,j,ig,2)=qtl_ge(i,j,ig,2)+q_g1(k1)*
     /                                       q_g2(k2)
	                    else if(n_allele==6)then
	                      qtl_ge(i,j,ig,3)=qtl_ge(i,j,ig,3)+q_g1(k1)*
     /                                       q_g2(k2)
	                    else if(n_allele==7)then
	                      qtl_ge(i,j,ig,4)=qtl_ge(i,j,ig,4)+q_g1(k1)*
     /                                       q_g2(k2)
	                    else if(n_allele==8)then
	                      qtl_ge(i,j,ig,5)=qtl_ge(i,j,ig,5)+q_g1(k1)*
     /                                       q_g2(k2)
	                    end if
	                  end do
	                end do 
	                s=0.0d0
	                do k1=1,5
	                  s=s+qtl_ge(i,j,ig,k1)
	                end do
                      ss=ss+s
	              end if 
	            end if 
	          end do
	        end do
	      end do
	    end do
	    i_g(i,j)=ig
	  end do
	end do

      return
      end
C     ****************************************************************
      subroutine pairing1(ip,mp,pair1,pair2)
      
      implicit double precision(a-h,o-z)
      
      integer mp(3,4),pair1(3,2),pair2(3,2)
      
      if(ip==1)then        
        do i=1,3
          pair1(i,1)=mp(i,1)
          pair1(i,2)=mp(i,2)
          pair2(i,1)=mp(i,3)
          pair2(i,2)=mp(i,4)
        end do
      else if(ip==2)then    
        do i=1,3
          pair1(i,1)=mp(i,1)
          pair1(i,2)=mp(i,3)
          pair2(i,1)=mp(i,2)
          pair2(i,2)=mp(i,4)
        end do
      else if(ip==3)then    
        do i=1,3
          pair1(i,1)=mp(i,1)
          pair1(i,2)=mp(i,4)
          pair2(i,1)=mp(i,2)
          pair2(i,2)=mp(i,3)
        end do
      end if
      
      return
      end
C     ****************************************************************
      subroutine possgamete_b1(pair1,pair2,mode_i,mg,n_mode)
      
      implicit double precision(a-h,o-z)
      
      integer pair1(3,2),pair2(3,2),mg(4,2,2)
      
      mg=0
      m=0
      if(mode_i==1)then      
        n_mode=4
        do i=1,2
          do j=1,2
            m=m+1
            mg(m,1,1)=pair1(1,i)
            mg(m,1,2)=pair2(1,j)
            mg(m,2,1)=pair1(3,i)
            mg(m,2,2)=pair2(3,j)
          end do
        end do
        
      else if(mode_i==2)then  
        n_mode=4
        do i=1,2
          do j=1,2
            do k=1,2
              if(j.ne.k)then
                m=m+1
                mg(m,1,1)=pair1(1,i)
                mg(m,1,2)=pair2(1,j)
                mg(m,2,1)=pair1(3,i)
                mg(m,2,2)=pair2(3,k)
              end if
            end do
          end do
        end do
        
      else if(mode_i==3)then  
        n_mode=4
        do i=1,2
          do j=1,2
            do k=1,2
              if(i.ne.k)then
                m=m+1
                mg(m,1,1)=pair1(1,i)
                mg(m,1,2)=pair2(1,j)
                mg(m,2,1)=pair1(3,k)
                mg(m,2,2)=pair2(3,j)
              end if
            end do
          end do
        end do
      
      else if(mode_i==4)then  
        n_mode=4
        do i=1,2
          do j=1,2
            do k=1,2
              do l=1,2
                if(i.ne.k)then
                  if(j.ne.l)then
                    m=m+1
                    mg(m,1,1)=pair1(1,i)
                    mg(m,1,2)=pair2(1,j)
                    mg(m,2,1)=pair1(3,k)
                    mg(m,2,2)=pair2(3,l)
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
      subroutine proba_g_o_b1(ii,r,g_o)
      
      implicit double precision(a-h,o-z)
      
      if(ii==1)then
        g_o=((1.0d0-r)**2)/12.0d0
      
      else if(ii==2)then
        g_o=r*(1.0d0-r)/12.0d0
        
      else if(ii==3)then
        g_o=r*(1.0d0-r)/12.0d0
        
      else if(ii==4)then
        g_o=(r**2)/12.0d0
      end if
      
      return
      end
C     ****************************************************************      
      subroutine proba_g_o_multilocus_b1(rec,o_ind,i,j,i1,j1,m,n,
     /           imark,mloci,ftype,mtype,pro_l,pro_r)
     
      implicit double precision(a-h,o-z)
      
      parameter(mxloci=400)
      
      integer o_ind(mxloci,8),ftype(mxloci,4),mtype(mxloci,4)
      
      dimension trans(4,4),trans_t(4,4)
      dimension rec(mxloci)
      dimension zy_matrix_l1(16),zy_matrix_l2(16)
      dimension zy_matrix_r1(16),zy_matrix_r2(16)
      dimension trans_zy(16,16),trans_t_zy(16,16)
      
      
      call ge_interval_b1(i1,m,i_ga_f1,i_ga_f2)   
      call ge_interval_b1(j1,n,i_ga_m1,i_ga_m2)   
      zy_matrix_l1=0.0d0
      zy_matrix_r1=0.0d0
      i2=(i_ga_f1-1)*4+i_ga_m1
      j2=(i_ga_f2-1)*4+i_ga_m2
      zy_matrix_l1(i2)=1.0d0
      zy_matrix_r1(j2)=1.0d0
      if(imark.ne.(mloci-1))then
        do is=imark+1,mloci-1
          ii=is+1
          r=rec(ii)
          call tra_matrix_b1(r,trans)   
          do k1=1,4
            do l1=1,4
              do k2=1,4
                do l2=1,4
                  jj=(k1-1)*4+l1
                  ll=(k2-1)*4+l2
                  trans_zy(jj,ll)=trans(k1,k2)*trans(l1,l2)
                end do              
              end do
            end do
          end do
          s=0.0d0
          do k1=1,16
	      do l1=1,16
	        s=s+trans_zy(k1,l1)
	      end do
          end do 
          zy_matrix_r2=0.0d0
          do l=1,16
	      do k=1,16
	        zy_matrix_r2(l)=zy_matrix_r2(l)+trans_zy(k,l)*
     /                        zy_matrix_r1(k)
	      end do
	    end do
	    call check_ph_b1(i,j,ii,zy_matrix_r2,ftype,mtype,o_ind,
     /         zy_matrix_r1)
	  end do
	  pro_r=0.0d0
	  do k=1,16
	    pro_r=pro_r+zy_matrix_r1(k)
	  end do
	else if(imark.eq.(mloci-1))then
	  pro_r=1.0d0
	end if 
C     ***************************THE LEFT PART************************
      if(imark.ne.1)then
	  do is=imark,2,-1        
	    ii=is-1
	    r=rec(is)
	    call tra_matrix_b1(r,trans)
	    do l=1,4
	      do k=1,4
	        trans_t(l,k)=trans(k,l)
	      end do
	    end do
	    do k1=1,4
	      do l1=1,4
	        do k2=1,4
	          do l2=1,4
	            jj=(k1-1)*4+l1
	            ll=(k2-1)*4+l2
	            trans_t_zy(jj,ll)=trans_t(k1,k2)*trans_t(l1,l2)
	          end do
	        end do
	      end do
	    end do
	    zy_matrix_l2=0.0d0
	    do l=1,16
	      do k=1,16
	        zy_matrix_l2(l)=zy_matrix_l2(l)+trans_t_zy(l,k)*
     /                        zy_matrix_l1(k)
	      end do
	    end do
          call check_ph_b1(i,j,ii,zy_matrix_l2,ftype,mtype,o_ind,
     /                   zy_matrix_l1)
        end do
        pro_l=0.0d0
	  do k=1,16
	    pro_l=pro_l+zy_matrix_l1(k)
	  end do
      else if(imark.eq.1)then
	  pro_l=1.0d0
	end if
	
	return
	end 
C     ****************************************************************
      subroutine ge_interval_b1(ii,im,i_gamete1,i_gamete2)
      
      implicit double precision(a-h,o-z)
      
      m=0
      if(ii==1)then  
        do i=1,2
          do j=1,2
            m=m+1
            if(m==im)then
              call num_gamete_b1(i,j,i_gamete1) 
              call num_gamete_b1(i,j,i_gamete2)
            end if
          end do
        end do
        
      else if(ii==2)then
        do i=1,2
          do j=1,2
            do k=1,2
              if(j.ne.k)then
                m=m+1
                if(m==im)then
                  call num_gamete_b1(i,j,i_gamete1)
                  call num_gamete_b1(i,k,i_gamete2)
                end if
              end if
            end do
          end do
        end do
        
      else if(ii==3)then  
        do i=1,2
          do j=1,2
            do k=1,2
              if(i.ne.k)then
                m=m+1
                if(m==im)then
                  call num_gamete_b1(i,j,i_gamete1)
                  call num_gamete_b1(k,j,i_gamete2)
                end if
              end if
            end do
          end do
        end do
        
      else if(ii==4)then  
        do i=1,2
          do j=1,2
            do k=1,2
              do l=1,2
                if(i.ne.k)then
                  if(j.ne.l)then
                    m=m+1
                    if(m==im)then
                      call num_gamete_b1(i,j,i_gamete1)
                      call num_gamete_b1(k,l,i_gamete2)
                    end if
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
      subroutine num_gamete_b1(i1,i2,num_g)
      
      implicit double precision(a-h,o-z)
      
      integer i1,i2,num_g,j1,j2
      
      if((i1.eq.1).and.(i2.eq.1))then
        num_g=1
      else if((i1.eq.1).and.(i2.eq.2))then
        num_g=2
      else if((i1.eq.2).and.(i2.eq.1))then
        num_g=3
      else if((i1.eq.2).and.(i2.eq.2))then
        num_g=4
      end if
      
      return
      end	       
C     **************************************************************
      subroutine proba_q_g_b1(ii,mm,pair1,pair2,r,r1,mq,q,q_g)
      
      implicit double precision(a-h,o-z)
      
      integer q,pair1(3,2),pair2(3,2),mq(4,2)
      
      dimension q_g(4)
      
      q=0
      mq=0
      m=0
      r2=(r-r1)/(1-2.0d0*r1)
      q_g=0.0d0
      if(ii==1)then
        do i=1,2
          do j=1,2
            m=m+1
            if(m==mm)then
              q=q+1
              mq(q,1)=pair1(2,i)
              mq(q,2)=pair2(2,j)
              q_g(q)=((1.0d0-r1)*(1.0d0-r2)/(1.0d0-r))**2          
              if(q_g(q)==0.0d0)then                                     
	          q=q-1
	        end if
	        
	        do k=1,2
	          if(k.ne.i)then
	            q=q+1
	            mq(q,1)=pair1(2,k)
	            mq(q,2)=pair2(2,j)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/((1.0d0-r)**2) 
	            if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                end if
              end do
              
              do k=1,2
	          if(k.ne.j)then
	            q=q+1
	            mq(q,1)=pair1(2,i)
	            mq(q,2)=pair2(2,k)
	            q_g(q)=r1*(1.0d0-r1)*r2*(1.0d0-r2)/((1.0d0-r)**2) 
	            if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                end if
              end do
              
              do k=1,2
                do l=1,2
 	            if(k.ne.i)then
 	              if(l.ne.j)then
	                q=q+1
	                mq(q,1)=pair1(2,k)
	                mq(q,2)=pair2(2,l)
	                q_g(q)=(r1**2)*(r2**2)/((1.0d0-r)**2)         
	                if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
                    end if
                  end if
                end do
              end do
            end if
              
          end do
        end do
      else if(ii==2)then
        do i=1,2
          do j=1,2
            do k=1,2
              if(j.ne.k)then
                m=m+1
                if(m==mm)then
                  q=q+1                                           
                  mq(q,1)=pair1(2,i)
                  mq(q,2)=pair2(2,k)
                  q_g(q)=r1*(1.0d0-r1)*((1.0d0-r2)**2)/(r*(1.0d0-r))
                  if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1                                            
                  mq(q,1)=pair1(2,i)
                  mq(q,2)=pair2(2,j)
                  q_g(q)=r2*(1.0d0-r2)*((1.0d0-r1)**2)/(r*(1.0d0-r))
                  if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  do l=1,2
                    if(i.ne.l)then
                      q=q+1                                         
                      mq(q,1)=pair1(2,l)
                      mq(q,2)=pair2(2,k)
                      q_g(q)=(r1**2)*(r2*(1.0d0-r2))/(r*(1.0d0-r))
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1                                        
                      mq(q,1)=pair1(2,l)
                      mq(q,2)=pair2(2,j)
                      q_g(q)=(r2**2)*(r1*(1.0d0-r1))/(r*(1.0d0-r))
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do
	          end if
	        end if
	      end do
	    end do
	  end do
      else if(ii==3)then
        do i=1,2
          do j=1,2
            do k=1,2
              if(i.ne.k)then
                m=m+1
                if(m==mm)then
                  q=q+1                                            
                  mq(q,1)=pair1(2,k)
                  mq(q,2)=pair2(2,j)
                  q_g(q)=r1*(1.0d0-r1)*((1.0d0-r2)**2)/(r*(1.0d0-r))
                  if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
	            q=q+1                                            
                  mq(q,1)=pair1(2,i)
                  mq(q,2)=pair2(2,j)
                  q_g(q)=r2*(1.0d0-r2)*((1.0d0-r1)**2)/(r*(1.0d0-r))
                  if(q_g(q)==0.0d0)then
	              q=q-1
	            end if
                  do l=1,2
                    if(j.ne.l)then
                      q=q+1                                         
                      mq(q,1)=pair1(2,k)
                      mq(q,2)=pair2(2,l)
                      q_g(q)=(r1**2)*(r2*(1.0d0-r2))/(r*(1.0d0-r))
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1                                         
                      mq(q,1)=pair1(2,i)
                      mq(q,2)=pair2(2,l)
                      q_g(q)=(r2**2)*(r1*(1.0d0-r1))/(r*(1.0d0-r))
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end do
	          end if
	        end if
	      end do
	    end do
	  end do
      else if(ii==4)then
        do i=1,2
          do j=1,2
            do k=1,2
              do l=1,2
                if(i.ne.k)then
                  if(j.ne.l)then
                    m=m+1
                    if(m==mm)then
                      q=q+1                                        
                      mq(q,1)=pair1(2,k)
                      mq(q,2)=pair2(2,l)
                      q_g(q)=(r1**2)*((1.0d0-r2)**2)/(r**2)
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1                                        
                      mq(q,1)=pair1(2,k)
                      mq(q,2)=pair2(2,j)
                      q_g(q)=r2*(1.0d0-r2)*r1*(1.0d0-r1)/(r**2)
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
                      q=q+1                                         
                      mq(q,1)=pair1(2,i)
                      mq(q,2)=pair2(2,l)
                      q_g(q)=(r1*(1.0d0-r1))*(r2*(1.0d0-r2))/(r**2)
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	                q=q+1                                         
                      mq(q,1)=pair1(2,i)
                      mq(q,2)=pair2(2,j)
                      q_g(q)=(r2**2)*((1.0d0-r1)**2)/(r**2)
                      if(q_g(q)==0.0d0)then
	                  q=q-1
	                end if
	              end if
	            end if
	          end if
	        end do
	      end do
	    end do
	  end do
      end if
	return
	end
C     ************************************************************
      subroutine qtl_distribution_b1(fqtl,mqtl,fre_qtl,nqtl)
      
      implicit double precision(a-h,o-z)
      
      integer fqtl(4),mqtl(4)    
      integer g1(6,2),g2(6,2)  
      integer zygote(36,4)       
      integer nqtl(5)           
      
      dimension fre_qtl(5)        
      dimension fre1(6),fre2(6)  
      dimension fre(36)          
      
      call gametes_b1(g1,fqtl,fre1) 
      call gametes_b1(g2,mqtl,fre2)
      
      do i=1,6
        do j=1,6
          k=(i-1)*6+j
          zygote(k,1)=g1(i,1)
          zygote(k,2)=g1(i,2)
          zygote(k,3)=g2(j,1)
          zygote(k,4)=g2(j,2)
          fre(k)=fre1(i)*fre2(j)
        end do
      end do
      
      fre_qtl=0.0d0
      do k=1,36
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
      
C     ***************************************************************
      subroutine gametes_b1(g,qtl,fre)
      
      implicit double precision(a-h,o-z)
      
      integer g(6,2),qtl(4)
      
      dimension fre(6)
      
      g(1,1)=qtl(1)
      g(1,2)=qtl(2)
      fre(1)=1.0d0/6.0d0   
      
      g(2,1)=qtl(1)
      g(2,2)=qtl(3)
      fre(2)=1.0d0/6.0d0   
      
      g(3,1)=qtl(1)
      g(3,2)=qtl(4) 
      fre(3)=1.0d0/6.0d0   
      
      g(4,1)=qtl(2)
      g(4,2)=qtl(3)
      fre(4)=1.0d0/6.0d0   
      
      g(5,1)=qtl(2)
      g(5,2)=qtl(4)
      fre(5)=1.0d0/6.0d0   
      
      g(6,1)=qtl(3)
      g(6,2)=qtl(4)
      fre(6)=1.0d0/6.0d0   
      
      return
      end	
C     *********************************************/****************
      subroutine emalgorithm_b1(fqtl,mqtl,n,trait,g_o,q_g,n_ig,kqtl,
     /                         kmeans,u,sigma,slikelihood,sLOD)
      
      implicit double precision(a-h,o-z)
      
      parameter(mxind=500,max_g=48*48,pi=3.1415926d0)
      
      integer n_ig(mxind,3,3)
      integer fqtl(4),mqtl(4),nqtl(5)
      
      dimension g_o(mxind,3,3,max_g),q_g(mxind,3,3,max_g,5)
      dimension u1(5),u2(5),trait(mxind),g_theo(5)
      dimension w(mxind,5),q(mxind,5)
      dimension fre_qtl(5)
      dimension u(5),kmeans(5)
      
      call qtl_distribution_b1(fqtl,mqtl,fre_qtl,nqtl)
      u1=0.0d0
      do i=1,5                     
        if(nqtl(i).eq.1)then
          u1(i)=kmeans(i)
        end if
      end do
      u1=5.0d0
      sum=0.0d0
      do i=1,n
        sum=sum+trait(i)
      end do
      gmean=sum/n
      sum=0.0d0
      do i=1,n                     
        sum=sum+(trait(i)-gmean)**2
      end do
      sigma1=sqrt(sum/(n-1)) 
      q=0.0d0
      do k=1,5
        do i=1,n
          do ip1=1,3
            do ip2=1,3
              ig=n_ig(i,ip1,ip2)
              do m=1,ig
                q(i,k)=q(i,k)+q_g(i,ip1,ip2,m,k)*g_o(i,ip1,ip2,m)
              end do
            end do
          end do
        end do
      end do   
      slikelihood1=0.0d0
      do i=1,n
        s2=0.0d0
        do ip1=1,3
          do ip2=1,3
            do j=1,n_ig(i,ip1,ip2)
              s1=0.0d0
              do k=1,5
                s1=s1+(1.0d0/sqrt(2.0d0*pi*(sigma1**2)))*
     /             exp(-((trait(i)-u1(k))**2)/(2.0d0*(sigma1**2)))
     /             *q_g(i,ip1,ip2,j,k)*g_o(i,ip1,ip2,j)
	        end do
	        s2=s2+s1
	      end do
	    end do
	  end do
	  slikelihood1=slikelihood1+log(s2)
	end do  
	ss=0.0d0
	do i=1,n
	  do ip1=1,3
	    do ip2=1,3
	      ss1=0.0d0
	      do j=1,n_ig(i,ip1,ip2)
              ss2=0.0d0
	        do k=1,5
	          ss2=ss2+q_g(i,ip1,ip2,j,k)*g_o(i,ip1,ip2,j)
	        end do
	        ss1=ss1+ss2	    
	      end do
            ss=ss+ss1
          end do
        end do
	end do	
	difference=1.0d0
	do while(difference>0.0001d0) 
	  w=0.0d0
	  do i=1,n
	    do j=1,5
	      s2=0.0d0
	      do ip1=1,3
	        do ip2=1,3
	          do k=1,n_ig(i,ip1,ip2)
	            s1=0.0d0
	            do l=1,5
	              s1=s1+(1.0d0/sqrt(2.0d0*pi*(sigma1**2)))*
     /                 exp(-((trait(i)-u1(l))**2)/(2.0d0*(sigma1**2)))
     /                 *q_g(i,ip1,ip2,k,l)
	            end do
	            s2=s2+(((1.0d0/sqrt(2.0d0*pi*(sigma1**2)))*
     /               exp(-((trait(i)-u1(j))**2)/(2.0d0*(sigma1**2))))
     /               *q_g(i,ip1,ip2,k,j)*g_o(i,ip1,ip2,k))/s1
	          end do
	        end do
	      end do
	      w(i,j)=s2         
	    end do
	  end do
	  s=0.0d0
        do j=1,5
	    do i=1,n
	      s=s+w(i,j)
	    end do
	  end do 
	  
	  s3=0.0d0
        u2=0.0d0
	  do j=1,5
	    s1=0.0d0
          do i=1,n
	      s1=s1+w(i,j)
	    end do
	    s2=0.0d0
	    do i=1,n
	      s2=s2+w(i,j)*trait(i)
	    end do
	    if(s1.ne.0.0d0)then
	      u2(j)=s2/s1       
	      do i=1,n
	        s3=s3+((trait(i)-u2(j))**2)*w(i,j)
	      end do
	    end if
	  end do
        sigma2=sqrt(s3/(n*1.0d0))
        
        slikelihood2=0.0d0
	  do i=1,n
	    s2=0.0d0
	    do ip1=1,3
	      do ip2=1,3
	        do j=1,n_ig(i,ip1,ip2)
	          s1=0.0d0
	          do k=1,5
	            s1=s1+(1.0d0/sqrt(2.0d0*pi*(sigma2**2)))*
     /              exp(-((trait(i)-u2(k))**2)/(2.0d0*(sigma2**2)))*
     /              q_g(i,ip1,ip2,j,k)*g_o(i,ip1,ip2,j)
	          end do
	          s2=s2+s1
	        end do
	      end do
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
	  do ip1=1,3
	    do ip2=1,3
	      do j=1,n_ig(i,ip1,ip2)
	        s1=0.0d0
	        do k=1,5
 	          s1=s1+(1.0d0/sqrt(2.0d0*pi*(var**2)))*exp
     /             (-((trait(i)-u_mean)**2)/(2.0d0*(var**2)))*
     /              q_g(i,ip1,ip2,j,k)*g_o(i,ip1,ip2,j)
	        end do
	        s2=s2+s1
	      end do
	    end do
	  end do
	  s_likelihood=s_likelihood+(log(s2))
      end do
	sLOD=slikelihood-s_likelihood   
      return
      end     
C     ****************************************************************
      subroutine check_ph_b1(i,j,m,zy_matri,fp,mp,o_ind,genotype)
      
      implicit double precision(a-h,o-z)
      parameter(mxloci=400)
      
      integer m,mp(mxloci,4),fp(mxloci,4),o_ind(mxloci,8)
      integer gamete1(2),gamete2(2),zygote(4),nz(8)
      
      dimension zy_matri(16),genotype(16)
      
      genotype=0.0d0
      
      do i1=1,4
        do j1=1,4
          call off_b1(m,i,i1,fp,gamete1)
          call off_b1(m,j,j1,mp,gamete2)
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
            ii=(i1-1)*4+j1
            genotype(ii)=zy_matri(ii)
          end if
        end do
      end do
      
      return
      end     
          
C     ******************************************************
      subroutine off_b1(m,i,i1,parent,gamete)
      
      implicit double precision(a-h,o-z)
      
      parameter(mxloci=400)
      
      integer parent(mxloci,4),gamete(2),pair1(2),pair2(2)
      
      if(i==1)then          
        pair1(1)=parent(m,1)
        pair1(2)=parent(m,2)
        pair2(1)=parent(m,3)
        pair2(2)=parent(m,4)
      else if(i==2)then     
        pair1(1)=parent(m,1)
        pair1(2)=parent(m,3)
        pair2(1)=parent(m,2)
        pair2(2)=parent(m,4)
      else if(i==3)then     
        pair1(1)=parent(m,1)
        pair1(2)=parent(m,4)
        pair2(1)=parent(m,2)
        pair2(2)=parent(m,3)
      end if
      
      if(i1==1)then
        gamete(1)=pair1(1)
        gamete(2)=pair2(1)
      else if(i1==2)then
        gamete(1)=pair1(1)
        gamete(2)=pair2(2)
      else if(i1==3)then
        gamete(1)=pair1(2)
        gamete(2)=pair2(1)
      else if(i1==4)then
        gamete(1)=pair1(2)
        gamete(2)=pair2(2)
      end if
      
      return
      end     
C     ***********************************************************
      subroutine tra_matrix_b1(r,trans)
      
      implicit double precision(a-h,o-z)
      
      dimension trans(4,4)
      
      a=(1.0d0-r)**2
      b=r*(1.0d0-r)
      c=r**2
      
      trans=b
      do i=1,4
        trans(i,i)=a
      end do 
      trans(1,4)=c
      trans(2,3)=c
      trans(3,2)=c
      trans(4,1)=c
      
      return 
      end   
C     ***********************************************************
      subroutine parent_qtl_genotype_b1(fqtype,mqtype)
      
      implicit double precision(a-h,o-z)
      
      parameter(nqtype=92)
     
      integer fqtype(nqtype,4),mqtype(nqtype,4)
      integer p1(6,4),p2(6,4)
      
      iqtype=1
      do i=1,4
        do j=i+1,5
          if(i.eq.1)then
            if(j.ne.5)then
              call iqgenotype_b1(i,p1,iqtl1)
              call iqgenotype_b1(j,p2,iqtl2)
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
            call iqgenotype_b1(i,p1,iqtl1)
            call iqgenotype_b1(j,p2,iqtl2)
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
C     *****************************************************
      subroutine iqgenotype_b1(i,p,iqtl)
      
      implicit double precision(a-h,o-z)
      
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
C     ********************************************************
      subroutine check_b1(fqtl,mqtl,u,icheck)
      
      implicit double precision(a-h,o-z)
      
      integer fqtl(4),mqtl(4)
      
      dimension u(5)
      
      call number_parent_b1(fqtl,ip1)
      call number_parent_b1(mqtl,ip2)
      
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
      
C     ***********************************************************
      subroutine number_parent_b1(pqtl,ip)
      
      integer pqtl(4)
      
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
                  	
 
