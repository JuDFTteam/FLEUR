!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE gf_hsint 
      IMPLICIT NONE
      PUBLIC :: hsintstep
      PRIVATE ::dsjn, gammln, ISHIDA_SBSL, DFFCTR 
      CONTAINS 
      ! The subroutines in this module generate the interstitial part of
      ! Hamiltonian and the Overlap-Matrix. The Values are stored in a  
      ! 1-D array in the following arrangement:                         
      !      -g'                                                        
      !   (->    )                                                      
      ! +g(--->  )    Values <g'|H|g>                                   
      !   (----->)                                                      
                                                                        
      SUBROUTINE hsintstep(                                             &
     &     matsize,lapw,stars,                                          &
     &     napw,                                                        &
     &     stepf,jspin,mpi,                                             &
     &     bkpt,bbmat,nlotot,vpw,l_noco,                                &
     &     aa,bb)                                                       
!*********************************************************************  
!     Version of hsint which uses the potential setup of gf_potconv     
!     in a GF-calculation, in this CASE the potential is set            
!     to zero in the auxillary volumes                                  
!                                                                       
!     D. Wortmann, Juelich, 2001                                        
!*********************************************************************  
      USE m_gf_types 
      IMPLICIT NONE 
!     ..                                                                
!     .. Scalar Arguments ..                                            
      INTEGER, INTENT (IN)          :: jspin 
      INTEGER, INTENT (IN)          :: matsize,napw 
      INTEGER, INTENT (IN)          :: nlotot 
      TYPE(t_LAPW),INTENT(INOUT)    :: lapw 
      TYPE(t_stars),INTENT(IN)      :: stars 
      TYPE(t_mpi), INTENT(IN)       :: mpi 
      LOGICAL, INTENT(IN)           :: l_noco 
      COMPLEX, INTENT (IN)          :: stepf(-stars%mx1:,-stars%mx2:,-2 &
     &     *napw:)                                                      
      REAL,    INTENT (IN)          :: bkpt(3),bbmat(3,3) 
      COMPLEX, INTENT (IN)          :: vpw(:,:) 
      COMPLEX, INTENT (OUT)   :: aa(matsize),bb(matsize) 
!     ..                                                                
!     .. Local Scalars ..                                               
      COMPLEX th,ts 
      REAL b1(3),b2(3),r2 
      COMPLEX ust1,vp1 
      INTEGER i,i1,i2,i3,ii,in,j,ispin 
                                                                        
      INTEGER :: nc,istart 
                                                                        
!     ..                                                                
!     .. INTRINSIC Functions ..                                         
      INTRINSIC REAL 
!     ..                                                                
                                                                        
                                                                        
                                                                        
                                                                        
      aa = CMPLX(0.0,0.0) 
      bb = CMPLX(0.0,0.0) 
                                                                        
                                                                        
      ust1 = stepf(0,0,0) 
      ispin=jspin 
      lapw%nmat = lapw%nv(ispin) 
      IF (l_noco) THEN 
!         READ (npotmatfile) (vpw(ig3),ig3=1,ng3)                       
         ispin = 1 
         lapw%nmat=lapw%nv(1)+lapw%nv(2) 
      ENDIF 
                                                                        
      vp1 = vpw(1,ispin) 
                                                                        
!---  > loop over (k+g')                                                
      ii = 0 
      DO 30 i = mpi%n_rank+1, lapw%nv(ispin),mpi%n_size 
!---  >    loop over (k+g)                                              
         DO 20 j = 1,i - 1 
            ii = ii + 1 
!--   >     determine index and phase factor                            
            i1 = lapw%k%k1(i,ispin) - lapw%k%k1(j,ispin) 
            i2 = lapw%k%k2(i,ispin) - lapw%k%k2(j,ispin) 
            i3 = lapw%k%k3(i,ispin) - lapw%k%k3(j,ispin) 
!     +APW_LO                                                           
            b1(1) = bkpt(1)+lapw%k%k1(i,ispin) ;b2(1) = bkpt(1)         &
     &           +lapw%k%k1(j,ispin)                                    
            b1(2) = bkpt(2)+lapw%k%k2(i,ispin) ;b2(2) = bkpt(2)         &
     &           +lapw%k%k2(j,ispin)                                    
            b1(3) = bkpt(3)+lapw%k%k3(i,ispin) ;b2(3) = bkpt(3)         &
     &           +lapw%k%k3(j,ispin)                                    
            r2 = dot_product(b1,matmul(b2,bbmat)) 
!            r2 = dotirp(b1,b2,bbmat)                                   
            th = (0.5*r2*stepf(i1,i2,i3)) 
!                                                                       
!         add the potential if there is a corresponding matrix element  
!                                                                       
            IF (abs(i3)<=ubound(stars%ig,3)) THEN 
               in=stars%ig(i1,i2,i3) 
               IF (in>0) THEN 
!                  IF (stars%sk3(in)<=3*lapw%rkmax) THEN                
                     th = th+vpw(in,ispin) 
!                  ENDIF                                                
               ENDIF 
            ENDIF 
            ts = stepf(i1,i2,i3) 
            aa(ii) = th 
            bb(ii) = ts 
   20    CONTINUE 
!---  >    diagonal term (g-g'=0 always first star)                     
         ii = ii + 1 
         aa(ii) = 0.5*lapw%k%rk(i,ispin)*lapw%k%rk(i,ispin)*ust1 + vp1 
!     aa(ii) = vp1                                                      
         bb(ii) = ust1 
                                                                        
   30 END DO 
                                                                        
      IF (l_noco) THEN 
!+gb99                                                                  
      nc = int( 1. + (lapw%nv(1)+nlotot - mpi%n_rank - 1)/mpi%n_size ) 
      istart = mpi%n_rank + nc*mpi%n_size - (lapw%nv(1)+nlotot) 
      ii = (lapw%nv(1)+nlotot+1)*(lapw%nv(1)+nlotot+2)/2 - 1 
      ispin = 2 
!---> determine spin-down spin-down part of Hamiltonian- and ovlp-matrix
!---> reload V_22                                                       
!      READ (npotmatfile) (vpw(ig3),ig3=1,ng3)                          
      vp1 = real(vpw(1,ispin)) 
                                                                        
         DO i = istart+1, lapw%nv(ispin),mpi%n_size 
            nc=nc+1 
!---  >    loop over (k+g)                                              
            DO  j = 1,i - 1 
            ii = (nc-1)*( mpi%n_rank - mpi%n_size + 1 ) +               &
     &          mpi%n_size*(nc-1)*nc/2 +  lapw%nv(1)+nlotot + j         
!--   >     determine index and phase factor                            
            i1 = lapw%k%k1(i,ispin) - lapw%k%k1(j,ispin) 
            i2 = lapw%k%k2(i,ispin) - lapw%k%k2(j,ispin) 
            i3 = lapw%k%k3(i,ispin) - lapw%k%k3(j,ispin) 
!     +APW_LO                                                           
            b1(1) = bkpt(1)+lapw%k%k1(i,ispin) ;b2(1) = bkpt(1)         &
     &           +lapw%k%k1(j,ispin)                                    
            b1(2) = bkpt(2)+lapw%k%k2(i,ispin) ;b2(2) = bkpt(2)         &
     &           +lapw%k%k2(j,ispin)                                    
            b1(3) = bkpt(3)+lapw%k%k3(i,ispin) ;b2(3) = bkpt(3)         &
     &           +lapw%k%k3(j,ispin)                                    
            r2 = dot_PRODUCT(b1,MATMUL(b2,bbmat)) 
!            r2 = dotirp(b1,b2,bbmat)                                   
            th = (0.5*r2*stepf(i1,i2,i3)) 
!                                                                       
!         add the potential if there is a corresponding matrix element  
!                                                                       
            IF (ABS(i3) <= UBOUND(stars%ig,3)) THEN 
               in = stars%ig(i1,i2,i3) 
               IF (in>0) THEN 
!                  IF (stars%sk3(in) <= 2*lapw%rkmax) THEN              
                     th = th+vpw(in,ispin) 
!                  ENDIF                                                
               ENDIF 
            ENDIF 
            ts = stepf(i1,i2,i3) 
            aa(ii) = th 
            bb(ii) = ts 
         ENDDO 
         ii = ii + 1 
         aa(ii) = 0.5*lapw%k%rk(i,ispin)*lapw%k%rk(i,ispin)*ust1 + vp1 
         bb(ii) = ust1 
      ENDDO 
                                                                        
!---> determine spin-down spin-up part of Hamiltonian- and ovlp-matrix  
!---> reload real part of V_21                                          
!      READ (npotmatfile) (vpw(ig3),ig3=1,ng3)                          
      nc = int( 1. + (lapw%nv(1)+nlotot - mpi%n_rank - 1)/mpi%n_size ) 
!                                                                       
!---> loop over (k+g')                                                  
!                                                                       
      DO i = istart+1, lapw%nv(2), mpi%n_size 
         nc = nc + 1 
!--->    loop over (k+g)                                                
         DO j = 1,lapw%nv(1) 
!-gb99      ii = (lapw%nv(1)+i-1)*(lapw%nv(1)+i)/2 + j                  
            ii = (nc-1)*( mpi%n_rank - mpi%n_size + 1 ) +               &
     &            mpi%n_size*(nc-1)*nc/2 + j                            
!--->       determine index and phase factor                            
            i1 = lapw%k%k1(i,2) - lapw%k%k1(j,1) 
            i2 = lapw%k%k2(i,2) - lapw%k%k2(j,1) 
            i3 = lapw%k%k3(i,2) - lapw%k%k3(j,1) 
!                                                                       
!         add the potential if there is a corresponding matrix element  
!                                                                       
                                                                        
            IF (abs(i3) <= ubound(stars%ig,3)) THEN 
               in = stars%ig(i1,i2,i3) 
               IF (in>0) THEN 
!                  if (stars%sk3(in)<=2*lapw%rkmax) THEN                
                     aa(ii)  = vpw(in,3) 
!                  ENDIF                                                
               ENDIF 
            ENDIF 
         ENDDO 
      ENDDO 
      ENDIF 
                                                                        
                                                                        
      RETURN 
      END SUBROUTINE hsintstep 
                                                                        
#ifdef CPP_NO_LONGER_SUPPORTED                                          
      SUBROUTINE hsintsolwil(lapw,matsize,atoms,                        &
     &    ispin,cell,bkpt,                                              &
     &    k1d,k2d,vpw,                                                  &
     &    gfinp,                                                        &
     &    H,S             )                                             
!*********************************************************************  
!     modified version of the interstital setup                         
!     The setup is not done WITH the step-FUNCTION but by               
!     subtracting the MT-contribution afterwards                        
!     The Contribution of the Overlap+Kinetic Energy in the             
!     MT are subtracted ONLY up to lmax                                 
!                                                                       
!                                                                       
!     WARNING: THIS SUBROUTINE does not work correctly WITH non-symmorph
!     symmetry                                                          
!                                                                       
!                       Daniel Wortmann, Tokio 2001                     
!*********************************************************************  
                                                                        
      USE m_gf_types 
      USE m_constants, ONLY: pimach 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)::matsize 
      INTEGER,INTENT(IN)::k1d,k2d 
      INTEGER,INTENT(IN)::ispin 
                                                                        
      REAL,INTENT(IN)::bkpt(3) 
      TYPE(t_gfinp),INTENT(IN) ::gfinp 
      TYPE(t_cell),INTENT(IN)  ::cell 
      TYPE(t_atoms),INTENT(IN) ::atoms 
      TYPE(t_lapw),INTENT(IN)  ::lapw 
      COMPLEX,INTENT(IN)::                                              &
     &     vpw(-2*k1d:2*k1d,-2*k2d:2*k2d,-2*gfinp%napw:2*gfinp%napw)    
      COMPLEX,INTENT(OUT)::H(matsize),S(matsize) 
                                                                        
!locals                                                                 
      INTEGER::n1,n2,in,ik,nt,na,l 
      REAL   ::eps,bz,g,rk1,rk2,norm,fpi,jl,djl,cosg 
      REAL   ::b1(3),b2(3),gk(3),gk1(3),pos(3) 
      COMPLEX::form 
      COMPLEX,ALLOCATABLE::FF(:,:),O(:,:) 
      REAL,ALLOCATABLE::p(:) 
! for potential                                                         
      INTEGER::i1,i2,i3 
                                                                        
                                                                        
                                                                        
      INTRINSIC cmplx,epsilon,exp 
                                                                        
                                                                        
      ALLOCATE(FF(atoms%ntype,lapw%nvd)) 
      ALLOCATE(O(-gfinp%napw:gfinp%napw,-gfinp%napw:gfinp%napw)) 
      ALLOCATE(p(0:maxval(atoms%lmax)+2)) 
!                                                                       
      fpi=4*pimach() 
      eps=EPSILON(1.0) 
      bz=cell%bmat(3,3) 
                                                                        
                       !init arrays                                     
      H=CMPLX(0.0,0.0) 
      S=CMPLX(0.0,0.0) 
                                                                        
                                                                        
!                                                                       
!     first the kinetic energy and the overlapp                         
!                                                                       
                                                                        
! calculate z-part of overlapp matrix first                             
      DO n1=-gfinp%napw,gfinp%napw 
         DO n2=-gfinp%napw,gfinp%napw 
            IF (n1==n2) THEN 
               O(n1,n2)=2*cell%z1/cell%amat(3,3) 
            ELSE 
               g=(n2-n1)*bz 
               O(n1,n2)=1/cell%amat(3,3)*CMPLX(0,-1.0/g)*               &
     &      (EXP(CMPLX(0,g*cell%z1))-EXP(CMPLX(0.0,g*(-1.0)*cell%z1)))  
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
! calculate Form-Faktors                                                
      DO ik=1,lapw%nvd 
                                             ! + qss1                   
         gk1(1)= bkpt(1)+lapw%k%k1(ik,ispin) 
         gk1(2)= bkpt(2)+lapw%k%k2(ik,ispin) 
         gk1(3)= bkpt(3)+lapw%k%k3(ik,ispin) 
         gk=MATMUL(cell%bmat,gk1) 
         na=1 
         DO nt=1,atoms%ntype 
                                                    !cartesian coordinat
            pos=MATMUL(cell%amat,atoms%taual(:,na)) 
            na=na+atoms%neq(nt) 
            FF(nt,ik)=fpi*EXP(CMPLX(0.0,dot_PRODUCT(gk,pos))) 
                                ! atoms                                 
         ENDDO 
                                ! g-vectors                             
      ENDDO 
                                                                        
! Now the 3d-matrix                                                     
      in=0 
      DO n1=1,lapw%nv(ispin) 
         DO n2=1,n1 
            in=in+1 
            b1(1) = bkpt(1)+lapw%k%k1(n1,ispin) ; b2(1) =bkpt(1)        &
     &           +lapw%k%k1(n2,ispin)                                   
            b1(2) = bkpt(2)+lapw%k%k2(n1,ispin) ; b2(2) =bkpt(2)        &
     &           +lapw%k%k2(n2,ispin)                                   
            b1(3) = bkpt(3)+lapw%k%k3(n1,ispin) ; b2(3) =bkpt(3)        &
     &           +lapw%k%k3(n2,ispin)                                   
            rk1 = dot_product(b1,matmul(b1,cell%bbmat)) 
            rk2=dot_product(b2,matmul(b2,cell%bbmat)) 
            cosg = dot_product(b1,matmul(b2,cell%bbmat))/MAX(SQRT(rk1   &
     &           *rk2),eps)                                             
                                                                        
            !rk1 = dotirp(b1,b1,cell%bbmat)                             
            !rk2 = dotirp(b2,b2,cell%bbmat)                             
            !cosg = dotirp(b1,b2,cell%bbmat)/MAX(SQRT(rk1*rk2),eps)     
            IF (lapw%k%k1(n1,ispin)==lapw%k%k1(n2,ispin)              &
     &           .AND.lapw%k%k2(n1,ispin)==lapw%k%k2(n2,ispin)) THEN  
               S(in)=O(lapw%k%k3(n1,ispin),lapw%k%k3(n2,ispin)) 
               H(in)=0.5*dot_product(b1,matmul(b2,cell%bbmat))*S(in) 
               !H(in)=0.5*dotirp(b1,b2,cell%bbmat)*S(in)                
            ENDIF 
!            GOTO 333                                                   
            ! Adding the potential                                      
             i1 = lapw%k%k1(n1,ispin) - lapw%k%k1(n2,ispin) 
            i2 = lapw%k%k2(n1,ispin) - lapw%k%k2(n2,ispin) 
            i3 = lapw%k%k3(n1,ispin) - lapw%k%k3(n2,ispin) 
            H(in)=H(in)+vpw(i1,i2,i3) 
                                                                        
            ! subtract the contributions of the MT                      
            ! calculate legendre ploynoms                               
                                                                        
            p(0)=1.0 
            p(1)=cosg 
            DO l=1,maxval(atoms%lmax)-1 
               p(l+1)=REAL(l+l+1)/REAL(l+1)*cosg*                       &
     &              p(l) -REAL(l)/REAL(l+1)*p(l-1)                      
            ENDDO 
                                                                        
            DO nt=1,atoms%ntype 
               form=atoms%neq(nt)*CONJG(FF(nt,n1))*FF(nt,n2)/cell%area  &
     &              /cell%amat(3,3)                                     
               DO l=0,atoms%lmax(nt) 
                  norm=(2.0*l+1)/fpi*p(l) 
                  CALL ishida_sbsl(atoms%rmt(nt),SQRT(rk1),SQRT(rk2),jl &
     &                 ,l)                                              
                  CALL ishida_tbsl(atoms%rmt(nt),SQRT(rk1),SQRT(rk2),djl&
     &                 ,l)                                              
                  !update H,S                                           
                  S(in)=S(in)-norm*form*jl 
                  H(in)=H(in)-norm*form*djl 
                                !l-loop                                 
               ENDDO 
                                ! Atoms                                 
            ENDDO 
! 333        continue                                                   
                               !n2                                      
        ENDDO 
                                !n1                                     
      ENDDO 
                                                                        
      DEALLOCATE(ff,O,p) 
                                                                        
      RETURN 
      END SUBROUTINE hsintsolwil 
#endif                                                                  
                                                                        
                                                                        
                                                                        
! Now additional subroutines written by H. Ishida                       
! These are only needed in hsintsolwil                                  
!---*----1----*----2----*----3----*----4----*----5----*----6----*----7--
!     def. sbsl=\int_0^R (r^2 dr) j_l(ak1*r) j_l(ak2*r)                 
!                                                                       
      SUBROUTINE ISHIDA_SBSL(R,Ak1,Ak2,S,L) 
      IMPLICIT NONE 
!     /* input & output variables */                                    
      INTEGER L 
      REAL R, Ak1, Ak2, S 
!     /* working variables */                                           
      REAL f1, f2, ep1, ep2 
      f1 = Ak1*R 
      f2 = Ak2*R 
      ep1 = 1.D-10 
      ep2 = 1.D-6 
      IF ( ABS(Ak1-Ak2)>ep1 ) THEN 
         S = R*R*(Ak1*DSJN(L,f2)*DSJN(L+1,f1)-Ak2*DSJN(L,f1)            &
     &       *DSJN(L+1,f2))/(Ak1*Ak1-Ak2*Ak2)                           
      ELSEIF ( L>=1 ) THEN 
         S = 0.5D0*(DSJN(L,f1)**2-DSJN(L-1,f1)*DSJN(L+1,f1))*R**3 
      ELSEIF ( ABS(Ak1)>ep2 ) THEN 
         S = (0.5D0*R-SIN(2.D0*f1)*0.25D0/Ak1)/Ak1**2 
      ELSE 
         S = R**3/3.D0 
      ENDIF 
      RETURN 
      END SUBROUTINE ISHIDA_SBSL 
!*==emb_ishida_tbsl.f    processed by SPAG 5.11R  at 16:43 on 16 May 200
!---*----1----*----2----*----3----*----4----*----5----*----6----*----7--
!     def. tbsl=\int_0^R (r^2 dr)                                       
!       j_l(ak1*r)[-\Delta/2 + (1/2)*\delta(r-R)d/dr ] j_l(ak2*r)       
!                                                                       
      SUBROUTINE ISHIDA_TBSL(R,Ak1,Ak2,T,L) 
      IMPLICIT NONE 
!     /* input & output variables */                                    
      INTEGER L 
      REAL R, Ak1, Ak2, T 
!     /* working variables */                                           
      REAL f1, f2, ep1, ep2, c1, c2 
      ep1 = 1.D-10 
      ep2 = 1.D-6 
      f1 = Ak1*R 
      f2 = Ak2*R 
      IF ( ABS(Ak1-Ak2)>ep1 ) THEN 
         c1 = DSJN(L,f1) 
         c2 = DSJN(L,f2) 
         T = (Ak1*Ak2*Ak2*c2*DSJN(L+1,f1)-Ak2*Ak1*Ak1*c1*DSJN(L+1,f2))  &
     &       *0.5D0*R*R/(Ak1*Ak1-Ak2*Ak2) + 0.5D0*R*DBLE(L)*c1*c2       
      ELSEIF ( L>=1 ) THEN 
         c1 = DSJN(L,f1) 
         c2 = DSJN(L+1,f1) 
         T = 0.25D0*(Ak1**2)*(R**3)*(c1*c1-DSJN(L-1,f1)*c2)             &
     &       + 0.5D0*R*DBLE(L)*c1*c1 - 0.5D0*R*R*Ak1*c1*c2              
      ELSEIF ( ABS(Ak1)>ep2 ) THEN 
         c1 = DSJN(L,f1) 
         c2 = DSJN(L+1,f1) 
         T = (R-0.5D0*SIN(2.D0*f1)/Ak1)*0.25D0 - 0.5D0*R*R*Ak1*c1*c2 
      ELSE 
         T = 0.D0 
      ENDIF 
      RETURN 
      END SUBROUTINE ISHIDA_TBSL 
                                                                        
                                                                        
!*==dsjn.f    processed by SPAG 5.11R  at 16:43 on 16 May 2001          
                                                                        
!---*----1----*----2----*----3----*----4----*----5----*----6----*----7  
!#numpac#dsjn                revised on 1984-11-30                      
      FUNCTION DSJN(N,X) 
      IMPLICIT NONE 
!*** Start of declarations inserted by SPAG                             
      REAL allt, alt, am, c, cp, cq, d, DSJN, eps, EPSILON, f,          &
     &     fn, g, p, q, r, rpb2, s                                      
      REAL sg, t, u, X, xh, xmax, xmin, xx, y 
      INTEGER l, N 
!*** END of declarations inserted by SPAG                               
      REAL v 
      DATA eps, xmax, xmin/1.D-16, 3.53D15, -180.218D0/ 
      DATA rpb2/0.886226925452758014D+00/ 
      IF ( N<0 .OR. X<0.D0 .OR. X>=xmax ) THEN 
!     argument error                                                    
         DSJN = 0.D0 
!     CALL mgidd(5hdsjn ,n,x,dsjn                                       
!    *,36harg1 or arg2 lt 0 or arg2 gt 3.53d15,9)                       
         RETURN 
      ELSE 
         IF ( N>0 .AND. X==0.D0 ) GOTO 200 
         IF ( N==0 .AND. X==0.D0 ) THEN 
            DSJN = 1.D0 
            RETURN 
         ELSE 
            fn = N 
            v = fn + 0.5D0 
            xh = X*0.5D0 
            xx = -xh*xh 
            IF ( v+1.D0+xx>=0.D0 ) THEN 
!     power series                                                      
               IF ( v>48.D0 ) THEN 
!  10 allt=LOG(xh)*fn-dlgama(v+1.d0)                                    
                  allt = LOG(xh)*fn - GAMMLN(v+1.0) 
                  IF ( allt<xmin ) GOTO 200 
                  alt = EXP(allt)*rpb2 
               ELSE 
                  alt = X**N/DFFCTR(N+N+1) 
                  IF ( alt==0.D0 ) GOTO 200 
               ENDIF 
               d = 0.D0 
               s = 1.D0 
               t = 1.D0 
   10          d = d + 1.D0 
               t = t*xx/((v+d)*d) 
               s = t + s 
               IF ( ABS(t)>=eps ) GOTO 10 
               DSJN = s*alt 
               RETURN 
            ELSEIF ( X<18.D0 .OR. X<v*v*0.55D0 ) THEN 
!     backward recurrence                                               
               IF ( X<=10.D0 ) THEN 
                  u = 4.7D0*X + 43.D0 
                  am = 2.1D0*X + 17.D0 
               ELSEIF ( X>100.D0 ) THEN 
                  IF ( X>200.D0 ) GOTO 300 
                  u = 1.5D0*X + 116.D0 
                  am = 1.14*X + 42.D0 
               ELSE 
                  u = 1.5D0*X + 116.D0 
                  am = 1.26D0*X + 30.D0 
               ENDIF 
               IF ( v>=u ) GOTO 200 
               d = DINT((DMAX1(v-(am+X)*0.5D0,0.D0)+am+1.D0)*0.5D0)*2.D0 
               p = 0.D0 
               t = 0.D0 
               s = EPSILON(1.0) 
               c = 1.D0 
               l = 1 
            ELSE 
!     asymptotic expansion                                              
               r = X*8.D0 
               t = 1.D0 
               p = t 
               q = 0.D0 
               g = v + v 
               f = -1.D0 
   20          f = f + 2.D0 
               t = (g+f+f-1.D0)*(g-f-f+1.D0)*t/(f*r) 
               q = q + t 
               t = (f+f+1.D0+g)*(f+f+1.D0-g)*t/((f+1.D0)*r) 
               p = p + t 
               IF ( ABS(t)>=eps ) GOTO 20 
               s = SIN(X) 
               c = COS(X) 
               sg = 1 - MOD(N/2,2)*2 
               IF ( MOD(N,2)==1 ) THEN 
                  cp = -c 
                  cq = s 
               ELSE 
                  cp = s 
                  cq = c 
               ENDIF 
               DSJN = (p*cp+q*cq)*sg/X 
               RETURN 
            ENDIF 
         ENDIF 
      ENDIF 
  100 IF ( l>=0 ) THEN 
         p = c*s + p 
         c = (d-1.5D0)*d*c/((d+0.5D0)*(d-1.D0)) 
      ENDIF 
      r = (d+0.5D0)*s/xh - t 
      d = d - 1.D0 
      IF ( d==fn ) y = r 
      t = s 
      s = r 
      l = -l 
      IF ( d>0.D0 ) GOTO 100 
      DSJN = c*y/(c*s+p) 
      RETURN 
!     underflow                                                         
  200 DSJN = 0.D0 
      RETURN 
  300 DSJN = 0.D0 
!     CALL mgidd(5hdsjn ,n,x,dsjn                                       
!    *,40harg1 gt SQRT(arg2)*1.384 and arg2 gt 200,10)                  
      RETURN 
      END FUNCTION DSJN 
!*==dffctr.f    processed by SPAG 5.11R  at 16:44 on 16 May 2001        
                                                                        
!---*----1----*----2----*----3----*----4----*----5----*----6----*----7  
!#numpac#dffctr              revised on 1984-11-30                      
      FUNCTION DFFCTR(N) 
      IMPLICIT NONE 
!*** Start of declarations inserted by SPAG                             
      REAL DFFCTR, fact 
      INTEGER j, N 
!*** END of declarations inserted by SPAG                               
      DIMENSION fact(97) 
      DATA (fact(j),j=1,38)/0.10000000000000000D+01,                    &
     &      0.10000000000000000D+01, 0.20000000000000000D+01,           &
     &      0.30000000000000000D+01, 0.80000000000000000D+01,           &
     &      0.15000000000000000D+02, 0.48000000000000000D+02,           &
     &      0.10500000000000000D+03, 0.38400000000000000D+03,           &
     &      0.94500000000000000D+03, 0.38400000000000000D+04,           &
     &      0.10395000000000000D+05, 0.46080000000000000D+05,           &
     &      0.13513500000000000D+06, 0.64512000000000000D+06,           &
     &      0.20270250000000000D+07, 0.10321920000000000D+08,           &
     &      0.34459425000000000D+08, 0.18579456000000000D+09,           &
     &      0.65472907500000000D+09, 0.37158912000000000D+10,           &
     &      0.13749310575000000D+11, 0.81749606400000000D+11,           &
     &      0.31623414322500000D+12, 0.19619905536000000D+13,           &
     &      0.79058535806250000D+13, 0.51011754393600000D+14,           &
     &      0.21345804667687500D+15, 0.14283291230208000D+16,           &
     &      0.61902833536293750D+16, 0.42849873690624000D+17,           &
     &      0.19189878396251062D+18, 0.13711959580999680D+19,           &
     &      0.63326598707628506D+19, 0.46620662575398912D+20,           &
     &      0.22164309547669977D+21, 0.16783438527143608D+22,           &
     &      0.82007945326378916D+22/                                    
      DATA (fact(j),j=39,76)/0.63777066403145712D+23,                   &
     &      0.31983098677287777D+24, 0.25510826561258285D+25,           &
     &      0.13113070457687989D+26, 0.10714547155728480D+27,           &
     &      0.56386202968058351D+27, 0.47144007485205310D+28,           &
     &      0.25373791335626258D+29, 0.21686243443194443D+30,           &
     &      0.11925681927744341D+31, 0.10409396852733332D+32,           &
     &      0.58435841445947272D+32, 0.52046984263666662D+33,           &
     &      0.29802279137433109D+34, 0.27064431817106664D+35,           &
     &      0.15795207942839548D+36, 0.14614793181237599D+37,           &
     &      0.86873643685617512D+37, 0.81842841814930553D+38,           &
     &      0.49517976900801982D+39, 0.47468848252659721D+40,           &
     &      0.29215606371473169D+41, 0.28481308951595832D+42,           &
     &      0.17821519886598633D+43, 0.17658411549989416D+44,           &
     &      0.11227557528557139D+45, 0.11301383391993226D+46,           &
     &      0.72979123935621403D+46, 0.74589130387155294D+47,           &
     &      0.48896013036866340D+48, 0.50720608663265600D+49,           &
     &      0.33738248995437775D+50, 0.35504426064285920D+51,           &
     &      0.23954156786760820D+52, 0.25563186766285862D+53,           &
     &      0.17486534454335399D+54, 0.18916758207051538D+55,           &
     &      0.13114900840751549D+56/                                    
      DATA (fact(j),j=77,97)/0.14376736237359169D+57,                   &
     &      0.10098473647378693D+58, 0.11213854265140152D+59,           &
     &      0.79777941814291672D+59, 0.89710834121121214D+60,           &
     &      0.64620132869576255D+61, 0.73562883979319396D+62,           &
     &      0.53634710281748291D+63, 0.61792822542628292D+64,           &
     &      0.45589503739486048D+65, 0.53141827386660331D+66,           &
     &      0.39662868253352861D+67, 0.46764808100261092D+68,           &
     &      0.35299952745484047D+69, 0.42088327290234982D+70,           &
     &      0.32122956998390482D+71, 0.38721261107016184D+72,           &
     &      0.29874350008503149D+73, 0.36397985440595213D+74,           &
     &      0.28380632508077991D+75, 0.34942066022971404D+76/           
      IF ( N>=0 .AND. N<=96 ) THEN 
         DFFCTR = fact(N+1) 
         RETURN 
      ENDIF 
      DFFCTR = 0. 
!     CALL mgid(6hdffctr,n,dffctr,17harg lt 0 or gt 96,5)               
      RETURN 
      END FUNCTION DFFCTR 
!*==gammln.f    processed by SPAG 5.11R  at 16:44 on 16 May 2001        
!---*----1----*----2----*----3----*----4----*----5----*----6----*----7--
!  (c) copr. 1986-92 numerical recipes software s21o)ws$!.              
      FUNCTION GAMMLN(Xx) 
      IMPLICIT NONE 
      REAL GAMMLN, Xx 
      INTEGER j 
      REAL ser, stp, tmp, x, y, cof(6) 
      SAVE cof, stp 
      DATA cof, stp/76.18009172947146D0, -86.50532032941677D0,          &
     &     24.01409824083091D0, -1.231739572450155D0,                   &
     &     .1208650973866179D-2, -.5395239384953D-5,                    &
     &     2.5066282746310005D0/                                        
      x = Xx 
      y = x 
      tmp = x + 5.5D0 
      tmp = (x+0.5D0)*LOG(tmp) - tmp 
      ser = 1.000000000190015D0 
      DO j = 1, 6 
         y = y + 1.D0 
         ser = ser + cof(j)/y 
      ENDDO 
      GAMMLN = tmp + LOG(stp*ser/x) 
      RETURN 
      END FUNCTION GAMMLN 
                                                                        
                                                                        
                                                                        
      END                                           
