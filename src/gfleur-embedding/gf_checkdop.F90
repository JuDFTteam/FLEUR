!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_checkdop 
      use m_juDFT 
      IMPLICIT NONE
      !-----------------------------------------------                  
      ! Testing the potential...                                        
      !-----------------------------------------------                  
      PRIVATE 
                                                 !maximal allowed charge
      REAL,PARAMETER :: maxerr_cdn = 200000.0 
                                            !error                      
                                                  !maximal allowed poten
      REAL,PARAMETER :: maxerr_pot = 1000000.0 
      PUBLIC gf_checkdop 
      CONTAINS 
      !<--S: gf_checkdop                                                
      SUBROUTINE gf_checkdop(                                           &
     &     atoms,cell,stars,sphhar,sym,                                 &
     &     cdn,fpw,fr)                                                  
      !-----------------------------------------------                  
      !  The gf-version of checkdop                                     
      !  The key part of of the code is taken from FLEUR checkdop       
      !              (last modified: 04-07-03) D. Wortmann              
      !-----------------------------------------------                  
      USE m_constants, ONLY: pi_const, oUnit 
      USE m_ylm
      USE m_starf 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      TYPE(t_atoms),INTENT(IN)   :: atoms 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_stars),INTENT(IN)   :: stars 
      TYPE(t_sphhar),INTENT(IN)  :: sphhar 
      TYPE(t_sym),INTENT(IN)     :: sym 
                                                                        
      LOGICAL,INTENT(IN)         :: cdn 
      COMPLEX,INTENT(IN)         :: fpw(:,:) 
      REAL   ,INTENT(IN)         :: fr(:,0:,:,:) 
                                                                        
      !>                                                                
      !<--Locals                                                        
                                                                        
                                                                    !loo
      INTEGER             :: na,n,nn,j,i,jspin,jspins,lh,mem,ll1,lm 
                                                                    !sta
      REAL                :: av,dms(2),rms(2),s,maxerr 
      COMPLEX             :: sf3(stars%ng3) 
      COMPLEX             :: ylm((MAXVAL(atoms%lmax)+1)**2 ) 
      REAL                :: rcc(3),x(3) 
                                     !No of testing points              
      INTEGER,PARAMETER   :: np = 500
      REAL                :: p(3,np) 
      REAL                :: v1(np,SIZE(fpw,2)),v2(np,SIZE(fpw,2)) 
      maxerr=0.0 
      jspins = MIN(SIZE(fr,4),SIZE(fpw,2)) 
                                                                        
      !>                                                                
      WRITE(oUnit,*) 
      IF (cdn) THEN 
         WRITE(oUnit,*) "Checking continuity of the charge density" 
      ELSE 
         WRITE(oUnit,*) "Checking continuity of the potential" 
      ENDIF 
      !<-- Loop over all atoms                                          
                                                                        
                                                                        
      na=0 
      DO n = 1,atoms%ntype 
         DO nn = 1,atoms%neq(n) 
            na = na+1 
            !<--Write info header                                       
            IF (jspins == 1) THEN 
               WRITE(oUnit,8000) n,nn 
            ELSE 
               WRITE(oUnit,8010) n,nn 
            ENDIF 
 8000       FORMAT (/,'    int.-m.t. boundary: atom type =',i2,' atom=' &
     &           ,i2,/,t2,'      int-coord',t20,'      cart-coord',t42  &
     &           ,' inter.   m. t. ')                                   
 8010       FORMAT (/,'    int.-m.t. boundary: atom type =',i2,' atom=' &
     &           ,i2,/,t2,'      int-coord',t20,'      cart-coord',t42  &
     &           ,' inter. (up) m. t. ',t62,' inter. (down) m. t. ')    
            !>                                                          
            !<-- Generate testing points...                             
                                                                        
            CALL priv_sphpts(p,atoms%rmt(n),atoms%pos(:,na)) 
                                                                        
            !>                                                          
            !<-- Loop over all points..                                 
                                                                        
            DO j = 1,np 
               !<-- Calculate values in Interstitial                    
                                                                        
               rcc=matmul(cell%bmat,p(:,j))/2./pi_const 
               !CALL cotra1(p(:,j),rcc,cell%bmat)                       
               CALL starf3(&
     &              sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot      &
     &              ,sym%tau,rcc,sym%invtab,sf3)                        
               DO jspin = 1,jspins 
                  v1(j,jspin) = SUM(REAL(fpw(:,jspin)*sf3(:))           &
     &                 *stars%nstr(:))                                  
               ENDDO 
                                                                        
               !>                                                       
               !<-- Values in MT                                        
                                                                        
               x(:) = p(:,j) - atoms%pos(:,na) 
               IF (sym%ngopr(na) /= 1) THEN 
                  !switch to internal units, rotate and go back to      
                  !cartesian                                            
                  rcc = MATMUL(cell%bmat,x)/2./pi_const 
                  !CALL cotra1(x,rcc,cell%bmat)                         
                  rcc = MATMUL(sym%mrot(:,:,sym%ngopr(na)),rcc) 
                  x=matmul(cell%amat,rcc) 
                  !CALL cotra0(rcc,x,cell%amat)                         
               ENDIF 
               CALL ylm4(atoms%lmax(n),x,ylm)                                                
               v2(j,:) = 0.0 
               DO lh   = 0,sphhar%nlh(sym%ntypsy(na)) 
                  s    = 0.0 
                  ll1  = sphhar%llh(lh,sym%ntypsy(na))                &
     &                 * (sphhar%llh(lh,sym%ntypsy(na)) + 1 ) + 1     
                  DO mem = 1,sphhar%nmem(lh,sym%ntypsy(na)) 
                     lm  = ll1 + sphhar%mlh(mem,lh,sym%ntypsy(na)) 
                     s   = s + REAL(sphhar%clnu(mem,lh                  &
     &                    ,sym%ntypsy(na))* ylm(lm))                  
                  ENDDO 
                  IF (cdn) s = s / (atoms%rmt(n)*atoms%rmt(n)) 
                  DO jspin=1,jspins 
                     v2(j,jspin) = v2(j,jspin) + fr(atoms%jri(n),lh,n   &
     &                    ,jspin)*s                                     
                  ENDDO 
               ENDDO 
               !<--Write info                                           
                                                                        
               !IF (j <= 8) THEN
                  rcc=matmul(cell%bmat,p(:,j))/2./pi_const 
                  !CALL cotra1(p(1,j),rcc,cell%bmat)                    
                  IF (jspins == 1) THEN 
                     WRITE (oUnit,8020) rcc,(p(i,j),i = 1,3),v1(j,1),v2(j,1)
                     write(66,8020) rcc,(p(i,j),i = 1,3),v1(j,1),v2(j,1)
                  ELSE 
                     WRITE (oUnit,8020) rcc,(p(i,j),i = 1,3),v1(j,1),v2(j,1)&
     &                    ,v1(j,2),v2(j,2)                              
                  ENDIF 

 8020             FORMAT (1x,2(3f8.4,1x),4f10.6)
               !ENDIF
                                                                        
               !>                                                       
                                                                        
               !>                                                       
            ENDDO 
                                                                        
            !>                                                          
            !<-- Calculate statistics                                   
                                                                        
            DO jspin = 1,jspins 
               CALL priv_fitchk(v1(:,jspin),v2(:,jspin),av              &
     &              ,rms(jspin),dms(jspin))                             
               maxerr=max(maxerr,abs(rms(jspin))*av) 
            ENDDO 
            IF (jspins == 1) THEN 
               WRITE(oUnit,8030) SUM(v1(:,1))/np,SUM(v2(:,1))/np 
            ELSE 
               WRITE(oUnit,8030) SUM(v1(:,1))/np,SUM(v2(:,1))/np,SUM(v1(:,2)&
     &              )/np,SUM(v2(:,2))/np                                
            ENDIF 
            WRITE(oUnit,8040) (rms(jspin),jspin = 1,jspins) 
            WRITE(oUnit,8050) (dms(jspin),jspin = 1,jspins) 
 8030       FORMAT (29x,'averages: ',4f10.6) 
 8040       FORMAT (5x,'rms (percent): ',2f15.6) 
 8050       FORMAT (5x,'dms (percent): ',2f15.6) 
                                                                        
            !>                                                          
                                                                        
         ENDDO 
      ENDDO 
                                                                        
      !>                                                                
      !<--check is missmatch is too big                                 
                                                                        
      IF(cdn.AND.maxerr>maxerr_cdn) THEN 
         WRITE(oUnit,*) "Too large missmatch in Charge" 
         WRITE(oUnit,*) maxerr,">",maxerr_cdn 
         WRITE(oUnit,*) "If this is acceptable, you might change the ",     &
     &        "compile-time value in gf_checkdop"                       
         CALL juDFT_error("gf_checkdop:INT-MT charge-test failed") 
      ENDIF 
      IF((.NOT.cdn).AND.maxerr>maxerr_pot) THEN 
         WRITE(oUnit,*) "Too large missmatch in Potential" 
         WRITE(oUnit,*) maxerr,">",maxerr_pot 
         WRITE(oUnit,*) "If this is acceptable, you might change the ",     &
     &        "compile-time value in gf_checkdop"                       
         CALL juDFT_error("gf_checkdop:INT-MT Potential-test failed") 
      ENDIF 
                                                                        
      !>                                                                
      RETURN 
      END SUBROUTINE gf_checkdop 
      !>                                                                
                                                                        
      !<-- S:priv_fitchk(f1,f2,n,av,rms,dmx)                            
      SUBROUTINE priv_fitchk(f1,f2,av,rms,dmx) 
!-----------------------------------------------                        
!     compare functions f1 and f2                                       
!         code taken from FLEUR                                         
!           (last modified: 05-03-22) D. Wortmann                       
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(OUT)    ::  av,dmx,rms 
      REAL   ,INTENT(IN)     ::  f1(:),f2(:) 
      !>                                                                
      !<--Locals                                                        
                                                                        
      REAL    :: d 
      INTEGER :: i,n 
      INTRINSIC abs,max,sqrt 
                                                                        
      !>                                                                
      n = SIZE(f1) 
      av = 0. 
      rms = 0. 
      dmx = 0. 
      DO 10 i = 1,n 
         av = av + f1(i) 
         d = (f1(i)-f2(i))**2 
         dmx = max(d,dmx) 
         rms = rms + d 
   10 END DO 
      av = av/n 
      IF (abs(av)<1.E-30) THEN 
         rms = 0. 
         dmx = 0. 
         RETURN 
      END IF 
      rms = sqrt(rms/n)/av*100. 
      dmx = sqrt(dmx)/av*100. 
      RETURN 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:priv_sphpts(p,r,pos)                                       
                                                                        
      SUBROUTINE priv_sphpts(p,r,pos) 
!-----------------------------------------------                        
!     generates points on a sphere                                      
!     e. wimmer     feb. 1980,m. weinert    jan. 1982                   
!         code taken from FLEUR                                         
!           (last modified: 05-03-22) D. Wortmann                       
!-----------------------------------------------                        
      USE m_constants, ONLY: pi_const, oUnit 
      IMPLICIT NONE 
      REAL,INTENT(IN)     :: r 
      REAL   ,INTENT(IN)  :: pos(3) 
      REAL   ,INTENT(OUT) :: p(:,:) 
                                                                        
                                                                        
      REAL                :: tc,phi,t,x(3) 
      REAL                :: sqrt13,sqrt7 
      INTEGER             :: i 
      sqrt13= SQRT(13.);sqrt7 = SQRT(7.0) 
!     ..                                                                
      DO i = 1,SIZE(p,2) 
         tc = 2.0*(real(i)*sqrt13-INT(i*sqrt13)) - 1.0 
         phi = 2.0*pi_const*(real(i)*sqrt7-AINT(i*sqrt7)) 
         t = sqrt(1.0-tc*tc) 
         x(1) = t*cos(phi) 
         x(2) = t*sin(phi) 
         x(3) = tc 
         p(:,i) = r*x + pos(:) 
      ENDDO 

      DO i = 1,SIZE(p,2)
         call random_number(tc)
         call random_number(phi)
         if (phi>0.5) tc=-1.*tc
         call random_number(phi)
         phi=phi*2.0*pi_const
         t = sqrt(1.0-tc*tc)
         x(1) = t*cos(phi)
         x(2) = t*sin(phi)
         x(3) = tc
         p(:,i) = r*x + pos(:)
      ENDDO

      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
