!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_plot 
      use m_juDFT
      IMPLICIT NONE
      PRIVATE 
      INTEGER,PARAMETER ::GF_PLOT_HARTREE=1 
      INTEGER,PARAMETER ::GF_PLOT_TOTPOT=2 
      INTEGER,PARAMETER ::GF_PLOT_CHARGE=3 
      INTEGER,PARAMETER ::GF_PLOT_DIFFRHO=4 
      INTEGER,PARAMETER ::GF_PLOT_DIFFPOT=5 
                                                                        
      PUBLIC gf_plot,GF_PLOT_HARTREE,GF_PLOT_TOTPOT,GF_PLOT_CHARGE      &
     &     ,gf_plot_cdn_INT,GF_PLOT_DIFFRHO,GF_PLOT_DIFFPOT             
      CONTAINS 
      !<--S: gf_plot                                                    
      SUBROUTINE gf_plot(layer,stars,cell,atoms,sym,jspind,             &
     &     fpw,ptype,sphhar,fr)                                         
!*****************************************************************      
! DESC:Provides an interface between the gf_subroutines and plot        
!                          Daniel Wortmann, Tue Mar  5 13:42:12 2002    
!*****************************************************************      
      USE m_gf_types 
      USE m_outcdn 
      USE m_xsf_io 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,INTENT(IN)        :: layer 
      TYPE(t_cell),INTENT(IN)   ::cell 
      TYPE(t_atoms),INTENT(IN)  ::atoms 
      TYPE(t_sym),INTENT(IN)    ::sym 
      TYPE(t_stars),INTENT(IN)  ::stars 
      INTEGER,INTENT(IN)        ::jspind 
      COMPLEX, INTENT (IN)      :: fpw(:,:) 
      INTEGER, INTENT(IN)       :: ptype 
                                                                        
      TYPE(t_sphhar),INTENT(IN),OPTIONAL :: sphhar 
      REAL,    INTENT (IN),OPTIONAL      :: fr(:,0:,:,:) 
                                                                        
                                                                        
      !<-- locals                                                       

      INTEGER                 :: nx,ny,nz,i,ii,iii,n,l 
      INTEGER                 :: jspins,jspin,na,itype 
      REAL                    :: pt(3),pt2(3) 
      INTEGER                 :: npoints(3) 
      REAL                    :: box(3,3),pos(3),s 
      CHARACTER               :: ctype 
      CHARACTER(LEN = 2)        :: app(2) 
      CHARACTER(LEN=30)       :: filename 
      LOGICAL                 :: l_plot,plpot 
      REAL,ALLOCATABLE        :: plot(:,:,:,:) 

      !>                                                                
      INQUIRE(FILE="gf_plot",EXIST=l_plot) 
      IF (.NOT.l_plot) RETURN 
      !OPEN the gf_plot-file                                            

      OPEN(99,FILE="gf_plot",FORM="formatted",STATUS="old") 
      DO 
         DO 
                                                !End jumps to end of sub
            READ(99,'(i1,a1)',END = 99) l,ctype 
            READ(99,'(a)') filename 
            READ(99,*) (pos(i),i=1,3) 
            READ(99,*) (box(i,1),i=1,3) 
            READ(99,*) (box(i,2),i=1,3) 
            READ(99,*) (box(i,3),i=1,3) 
            READ(99,*) (npoints(i),i=1,3) 
            !Now check if these settings are for the current type of plo
            IF (l == layer.AND.((ctype =='h'.OR.ctype=='H').AND.(ptype== GF_PLOT_HARTREE)).OR.((ctype =='v'.OR.ctype =='V')   &
     &           .AND.(ptype == GF_PLOT_Totpot)).OR.((ctype =='c'.OR.    &
     &           ctype =='C').AND.(ptype == GF_PLOT_CHARGE))) EXIT
            IF (l == layer.AND.(ctype =='p'.OR.ctype =='P').AND.(ptype==&
     &           GF_PLOT_DIFFPOT)) EXIT                                 
            IF (l == layer.AND.(ctype =='d'.OR.ctype =='D').AND.(ptype==&
     &           GF_PLOT_DIFFrho)) EXIT                                 
         ENDDO 
      !change to cart. coords                                           
         pos(1:2) = MATMUL(cell%amat(1:2,1:2),pos(1:2)) 
         box(1:2,:) = MATMUL(cell%amat(1:2,1:2),box(1:2,:)) 
      !for some plots a spin-loop might be needed                       
      IF ((ptype==GF_PLOT_Totpot).OR.(ptype==GF_PLOT_CHARGE)) THEN 
         jspins=MIN(jspind,2) 
         app(1)='.1' 
         app(2)='.2' 
      ELSE 
         jspins=1 
         app(1)='  ' 
      ENDIF 
      !Sometimes you plot a charge                                      
      IF (ptype==GF_PLOT_CHARGE) THEN 
         plpot=.FALSE. 
      ELSE 
         plpot=.TRUE. 
      ENDIF 

      ALLOCATE(plot(npoints(1),npoints(2),npoints(3),jspins)) 
      DO jspin=1,jspins 
         DO nz=1,npoints(3) 
            DO ny=1,npoints(2) 
               pointloop:DO nx=1,npoints(1) 
                  pt=pos+box(:,1)/(npoints(1)+1.)*nx                    &
     &                 +box(:,2)/(npoints(2)+1.)*ny                     &
     &                 +box(:,3)/(npoints(3)+1.)*nz                     
                  !Check if pt is in some atom                          
                  IF (PRESENT(fr).AND.PRESENT(sphhar)) THEN 
                     na = 0 
                     DO itype = 1,atoms%ntype 
                        DO n = 1,atoms%neq(itype) 
                           na = na+1 
                           DO i =-2,2 
                              DO ii =-2,2 
                                 pt2 = pt+i*cell%amat(:,1)+ii           &
     &                                *cell%amat(:,2)                   
                                 s = SQRT(dot_PRODUCT(pt2-atoms%pos(:,na&
     &                                ),pt2-atoms%pos(:,na)))           
                                 IF (s<atoms%rmt(itype)) THEN 
                                    plot(nx,ny,nz,jspin) =              &
     &                                   priv_cdn_MT(pt2,itype,na,plpot &
     &                                   ,atoms,cell,sphhar,sym,fr(:,0: &
     &                                   ,:,jspin))                     
                                    CYCLE pointloop 
                                 ENDIF 
                              ENDDO 
                           ENDDO 
                        ENDDO 
                     ENDDO 
                  ENDIF 
                  !In interstitial                                      
                  plot(nx,ny,nz,jspin) = gf_plot_cdn_INT(pt,cell%bmat   &
     &                 ,stars,sym,fpw(:,jspin))                         
               ENDDO pointloop 
            ENDDO 
         ENDDO 
      ENDDO 

                                                                        
      !<-- Output the plot                                              

      OPEN(98,FILE = TRIM(filename),FORM ='formatted'                   &
     &     ,STATUS ='replace')                                          
      CALL xsf_WRITE_atoms(98,atoms,.TRUE.,cell%amat)                                       
      CALL xsf_WRITE_header(98,.FALSE.,filename,box(:,1),box(:,2)       &
     &     ,box(:,3),pos,npoints)

      DO jspin = 1,jspins 
         WRITE(98,100) (((plot(i,ii,iii,jspin),i = 1,npoints(1)),ii =   &
     &        1,npoints(2)),iii = 1,npoints(3))                         
         IF (jspin /= jspins) THEN
             CALL xsf_WRITE_newblock(98,.FALSE.,box(:,1),box(:,2),box(:,3),pos,npoints)
         ENDIF
      ENDDO 
      CALL xsf_WRITE_endblock(98,.FALSE.) 
      CLOSE(98) 

      !>                                                                
      DEALLOCATE(plot) 
!       CALL juDFT_error("gf_plot",calledby="gf_plot.F90")
      ENDDO 
   99 CLOSE(99) 
  100 FORMAT (5E14.6) 
  101 FORMAT (3(i4,1x),3(f13.6,1x)) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- F: priv_cdn_MT()                                             
                                                                        
      FUNCTION  priv_cdn_MT(p,n,na,plpot,atoms,cell,sphhar,sym,rho)     &
     &     RESULT(xdnout)                                               
!-----------------------------------------------                        
!     calculte charge density (plpot=F) or potential in                 
!     muffin tin of atom type n, atom nr na at position p               
!     Code taken from FLEUR                                             
!             (last modified: 05-03-22) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_constants 
      USE m_ylm 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER, INTENT (IN)       :: n,na 
      LOGICAL, INTENT (IN)       :: plpot 
      TYPE(t_atoms),INTENT(IN)   :: atoms 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_sphhar),INTENT(IN)  :: sphhar 
      TYPE(t_sym),INTENT(IN)     :: sym 
      REAL,    INTENT (IN)       :: rho(:,0:,:) 
      REAL,    INTENT (INOUT)    :: p(3) 
      REAL ::xdnout 
      !>                                                                
      !<-- Locals                                                       
                                                                        
      REAL    :: s,sx,xd1,xd2 
      INTEGER :: j,jr,lh,mem,nd,nopa,ll1,lm 
      COMPLEX :: ylm((MAXVAL(atoms%lmax)+1)**2) 
      REAL    :: rcc(3),x(3) 
                                                                        
      !>                                                                
                                                                        
      nd = sym%ntypsy(na) 
      nopa = sym%ngopr(na) 
      x(:) = p(:) - atoms%pos(:,na) 
      sx = sqrt(dot_product(x,x)) 
                                                                        
      IF (nopa/=1) THEN 
         !switch to internal units                                      
         rcc = MATMUL(cell%bmat,x)/2.0/pi_const 
         !rotate into representative                                    
         p = MATMUL(sym%mrot(:,:,nopa),rcc) 
         !switch back to cartesian units                                
         x = MATMUL(cell%amat,p) 
      END IF 
      DO  j = atoms%jri(n),2,-1 
         IF (sx>=atoms%rmsh(j,n)) EXIT 
      ENDDO 
      jr = j 
      CALL ylm4(atoms%lmax(n),x,ylm)                                                    
      xd1 = 0.0 
      xd2 = 0.0 
      DO lh = 0, sphhar%nlh(nd) 
         ll1 = sphhar%llh(lh,nd) * ( sphhar%llh(lh,nd) + 1 ) + 1 
         s = 0.0 
         DO mem = 1,sphhar%nmem(lh,nd) 
           lm = ll1 + sphhar%mlh(mem,lh,nd) 
           s = s + real( sphhar%clnu(mem,lh,nd)*ylm(lm) ) 
         ENDDO 
         IF (plpot) THEN 
            xd1 = xd1 + rho(jr,lh,n)*s 
         ELSE 
            xd1 = xd1 + rho(jr,lh,n)*s/ (atoms%rmsh(jr,n)               &
     &           *atoms%rmsh(jr,n))                                     
         END IF 
         IF (jr==1) CYCLE 
         IF (plpot) THEN 
            xd2 = xd2 + rho(jr+1,lh,n)*s 
         ELSE 
            xd2 = xd2 + rho(jr+1,lh,n)*s/                               &
     &            (atoms%rmsh(jr+1,n)*atoms%rmsh(jr+1,n))               
         END IF 
      ENDDO 
      IF (jr==1) THEN 
         xdnout = xd1 
      ELSE 
         xdnout = xd1 + (xd2-xd1) *                                     &
     &        (sx-atoms%rmsh(jr,n)) / (atoms%rmsh(jr+1,n)-atoms%rmsh(jr &
     &        ,n))                                                      
      END IF 
      RETURN 
      END FUNCTION 
                                                                        
      !>                                                                
                                                                        
      !<-- F: gf_plot_cdn_INT()                                         
                                                                        
      FUNCTION gf_plot_cdn_INT(p,bmat,stars,sym,qpw) RESULT(xdnout) 
!-----------------------------------------------                        
!  calculate charge density/potential at point p in interstitial        
!   Code taken from FLEUR                                               
!             (last modified: 05-03-22) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_constants 
      USE m_starf 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)         :: bmat(3,3) 
      TYPE(t_stars),INTENT(IN)   :: stars 
      TYPE(t_sym),INTENT(IN)     :: sym 
      COMPLEX,INTENT(IN)         :: qpw(:) 
      REAL   ,INTENT(IN)         :: p(3) 
      REAL                       :: xdnout 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: rcc(3) 
      COMPLEX             :: sf3(SIZE(qpw)) 
      INTEGER             :: k 
      !>                                                                
      rcc=matmul(bmat,p)/2./pi_const 
      CALL starf3(                                                &
     &     sym%nop,stars%ng3,sym%symor,stars%kv3,sym%mrot,sym%tau,rcc   &
     &     ,sym%invtab,sf3)                                             
!                                                                       
      xdnout = 0.0 
      DO k = 1,stars%ng3 
         xdnout = xdnout + real(qpw(k)*sf3(k))*stars%nstr(k) 
      ENDDO 
      END FUNCTION 
                                                                        
      !>                                                                
      END                                           
