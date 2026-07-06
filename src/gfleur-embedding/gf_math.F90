!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_math 
      use m_juDFT
      USE m_constants, ONLY: pi_const, tpi_const, oUnit
      IMPLICIT NONE
!-------------------------------------------------------------------    
!     This module creates easy to use interfaces to useful math and     
!     some lapack routines, to be used with care.                       
!     It simplifies the code but might decrease performance!            
!-------------------------------------------------------------------    
      PRIVATE 
      REAL,PARAMETER :: sqrt_pi_const = 1.77245385090551602729816748334 
      INTERFACE mat_inverse 
         MODULE PROCEDURE cmat_inverse 
         MODULE PROCEDURE rmat_inverse 
      END INTERFACE 
      INTERFACE lin_equation 
         MODULE PROCEDURE cmat_lineq 
         MODULE PROCEDURE cmat_lineq_single 
         MODULE PROCEDURE rmat_lineq 
         MODULE PROCEDURE rmat_lineq_single 
      END INTERFACE 
      INTERFACE mat_sqrt 
         MODULE PROCEDURE cmat_sqrt 
         MODULE PROCEDURE rmat_sqrt 
      END INTERFACE 
      INTERFACE trace 
         MODULE PROCEDURE trace_real 
         MODULE PROCEDURE trace_complex 
      END INTERFACE 
      INTERFACE eigenvalues 
         MODULE PROCEDURE cmat_eigenvalues 
         MODULE PROCEDURE rmat_eigenvalues 
      END INTERFACE 
      PUBLIC  mat_inverse,lin_equation,mat_sqrt,eigenvalues,trace 
      PUBLIC sort,imag2d,interpolate_analytic,triang,sort2 
      PUBLIC pi_const,tpi_const,sqrt_pi_const 
      PUBLIC cblas_matmul_sq2,cblas_matmul_sq 
      PUBLIC zmat_product 
      PUBLIC simple_broyden !,make_charge2D(mat,lapw,stars,n)
      CONTAINS 
      !<-- inversion of a real/complex matrix                           
                                                                        
      !<-- F:cmat_inverse(A)                                            

      FUNCTION cmat_inverse(A) RESULT(Ainv) 
!     *****************************************                         
!     Use the LAPACK to invert the matrix A                             
!     *****************************************                         
      IMPLICIT NONE 
      COMPLEX, INTENT(IN) :: A(:,:) 
      COMPLEX             :: Ainv(SIZE(A,1),SIZE(A,1)) 
!     locals                                                            
      INTEGER :: nv2,info 
      COMPLEX :: work(4*SIZE(A,1)) 
      INTEGER :: ipiv(SIZE(A,1)) 
                                                                        
      nv2 = SIZE(A,1) 
      Ainv = A 
      IF (nv2 /= SIZE(A,2)) CALL juDFT_error                              &
     &     ('NO SQUARE MATRIX in matrix inverse')                       
      CALL zgetrf(Nv2,Nv2,Ainv,Nv2,ipiv,info) 
      IF ( info /= 0 )  THEN
            write(999,*) a
            CALL juDFT_error("cgetrf (inversion of A)")
      ENDIF
      CALL zgetri(Nv2,Ainv,Nv2,ipiv,work,4*nv2,info) 
      IF ( info /= 0 )  CALL juDFT_error("cgetri (inversion of A)") 
      END FUNCTION 

      !>                                                                
      !<-- F:rmat_inverse(A)                                            
      FUNCTION rmat_inverse(A) RESULT(Ainv) 
!     *****************************************                         
!     Use the LAPACK to invert the matrix A                             
!     *****************************************                         
      IMPLICIT NONE 
      REAL, INTENT(IN) :: A(:,:) 
      REAL             :: Ainv(SIZE(A,1),SIZE(A,1)) 
!     locals                                                            
      INTEGER :: nv2,info 
      REAL    :: work(4*SIZE(A,1)) 
      INTEGER :: ipiv(SIZE(A,1)) 
                                                                        
      nv2 = SIZE(A,1) 
      Ainv = A 
      IF (nv2 /= SIZE(A,2)) CALL juDFT_error                              &
     &     ('NO SQUARE MATRIX in matrix inverse')                       
      CALL dgetrf(Nv2,Nv2,Ainv,Nv2,ipiv,info) 
      IF ( info /= 0 )  CALL juDFT_error("sgetrf (inversion of A)") 
      CALL dgetri(Nv2,Ainv,Nv2,ipiv,work,4*nv2,info) 
      IF ( info /= 0 )  CALL juDFT_error("sgetri (inversion of A)") 
      END FUNCTION 
      !>                                                                
                                                                        
      !>                                                                
                                                                        
      !<-- functions for solving linear equations                       
      !<-- F:cmat_lineq(A)                                              
      FUNCTION cmat_lineq(A,b)RESULT(x) 
!     *****************************************                         
!     Use the LAPACK to get the solution of                             
!     Ax = b where A and b are complex matrices                         
!     x and b are m vectors of dimension n                              
!     *****************************************                         
      IMPLICIT NONE 
      COMPLEX,INTENT(IN) :: A(:,:),b(:,:) 
      COMPLEX            :: x(SIZE(b,1),SIZE(b,2)) 
                                                                        
      INTEGER           :: m,n,info,ipiv(SIZE(b,1)) 
      COMPLEX           :: AA(SIZE(A,1),SIZE(A,2)) 
      IF (SIZE(b,1) /= SIZE(a,2).OR.SIZE(a,1) /= SIZE(a,2)) CALL        &
     &     juDFT_error('Wrong dimensions in mat_lineq')                   
                   !a,b are intent in                                   
      x = b 
      AA = A 
      CALL zgesv(SIZE(aa,1),SIZE(b,2),AA,SIZE(AA,1),ipiv,x   &
     &     ,SIZE(x),info)                                               
                                                                        
      IF (info /= 0) CALL juDFT_error('mat_lineq failed') 
      END FUNCTION 
      !>                                                                
      !<-- F:cmat_lineq(A)                                              
      FUNCTION cmat_lineq_single(A,b)RESULT(x) 
!     *****************************************                         
!     Use the LAPACK to get the solution of                             
!     Ax = b where A and b are complex  matrices                        
!     x and b are m vectors of dimension n                              
!     *****************************************                         
      IMPLICIT NONE 
      COMPLEX,INTENT(IN) :: A(:,:),b(:) 
      COMPLEX            :: x(SIZE(b,1)) 
                                                                        
      INTEGER           :: m,n,info,ipiv(SIZE(b,1)) 
      COMPLEX           :: AA(SIZE(A,1),SIZE(A,2)) 
      IF (SIZE(b,1) /= SIZE(a,2).OR.SIZE(a,1) /= SIZE(a,2)) CALL        &
     &     juDFT_error('Wrong dimensions in mat_lineq')                   
                   !a,b are intent in                                   
      x = b 
      AA = A 
      CALL zgesv(SIZE(aa,1),1,AA,SIZE(AA,1),ipiv,x           &
     &     ,SIZE(x),info)                                               
                                                                        
      IF (info /= 0) CALL juDFT_error('mat_lineq failed') 
      END FUNCTION 
      !>                                                                
      !<-- F:rmat_lineq(A)                                              
      FUNCTION rmat_lineq(A,b)RESULT(x) 
!     *****************************************                         
!     Use the LAPACK to get the solution of                             
!     Ax = b where A and b are real matrices                            
!     x and b are m vectors of dimension n                              
!     *****************************************                         
      IMPLICIT NONE 
      REAL,INTENT(IN) :: A(:,:),b(:,:) 
      REAL            :: x(SIZE(b,1),SIZE(b,2)) 
                                                                        
      INTEGER         :: m,n,info,ipiv(SIZE(b,1)) 
      REAL            :: AA(SIZE(A,1),SIZE(A,2)) 
      IF (SIZE(b,1) /= SIZE(a,2).OR.SIZE(a,1) /= SIZE(a,2)) CALL        &
     &     juDFT_error('Wrong dimensions in mat_lineq')                   
                   !a,b are intent in                                   
      x = b 
      AA = A 
      CALL dgesv(SIZE(aa,1),SIZE(b,2),AA,SIZE(AA,1),ipiv,x   &
     &     ,SIZE(x),info)                                               
                                                                        
      IF (info /= 0) CALL juDFT_error('mat_lineq failed') 
      END FUNCTION 
      !>                                                                
      !<-- F:rmat_lineq(A)                                              
      FUNCTION rmat_lineq_single(A,b)RESULT(x) 
!     *****************************************                         
!     Use the LAPACK to get the solution of                             
!     Ax = b where A and b are real matrices                            
!     x and b are m vectors of dimension n                              
!     *****************************************                         
      IMPLICIT NONE 
      REAL,INTENT(IN) :: A(:,:),b(:) 
      REAL            :: x(SIZE(b,1)) 
                                                                        
      INTEGER         :: m,n,info,ipiv(SIZE(b,1)) 
      REAL            :: AA(SIZE(A,1),SIZE(A,2)) 
      IF (SIZE(b,1) /= SIZE(a,2).OR.SIZE(a,1) /= SIZE(a,2)) CALL        &
     &     juDFT_error('Wrong dimensions in mat_lineq')                   
                   !a,b are intent in                                   
      x = b 
      AA = A 
      CALL dgesv(SIZE(aa,1),1,AA,SIZE(AA,1),ipiv,x           &
     &     ,SIZE(x),info)                                               
                                                                        
      IF (info /= 0) CALL juDFT_error('mat_lineq failed') 
      END FUNCTION 
      !>                                                                
      !>                                                                
                                                                        
      !<-- Eigenproblem solvers                                         
      !<-- S:cmat_eigenvalues(A,ew,ev)                                  

      SUBROUTINE cmat_eigenvalues(A,ew,ev) 
!     *****************************************                         
!     Use the LAPACK to get the eigenvalues ew                          
!     and the right eigenvectors ev of the complex                      
!     square matrix A                                                   
!     *****************************************                         
      IMPLICIT NONE 
                                    !requires A_local                   
      COMPLEX,INTENT(IN)  :: A(:,:) 
      COMPLEX,INTENT(OUT) :: ew(:),ev(:,:) 
                                                                        
      INTEGER :: n,info 
      COMPLEX :: A_local(SIZE(A,1),SIZE(A,2)) 
      COMPLEX :: work(2*SIZE(a,1)) 
      REAL    :: rwork(2*SIZE(a,1)) 
      LOGICAL :: l_condnum 
      REAL b 
      INTEGER ilo,ihi 
      REAL scale(size(a,1)) 
      REAL abnrm 
      REAL rconde(size(a,1)) 
      REAL rcondv(size(a,1)) 
      COMPLEX,ALLOCATABLE:: work2(:) 
                                                                        
      INQUIRE(FILE='condnum',EXIST=l_condnum) 
                                                                        
      A_local = A 
      n = SIZE(A,1) 
      IF (SIZE(A,2) /= n) CALL juDFT_error                                &
     &     ('No square matrix in cmat_eigenvalues')                     
      IF (n /= SIZE(ew).OR.n /= SIZE(ev,1).OR.n /= SIZE(ev,2))          &
     &     CALL juDFT_error('Wrong dimensions in cmat_eigenvalues')       
                                                                        
                                                                        
      IF(.NOT.l_condnum)THEN 
         CALL zgeev('N','V',n,A_local,n,ew,                  &
     &        Ev,n,Ev,n,work,2*n,rwork,info)                            
         IF (info /= 0) CALL juDFT_error('cmat_eigenvalues failed') 
      ELSE 
         ALLOCATE( work2(n*n+2*n)) 
         CALL zgeevx('B','N','V','N',n,a_local,n,ew,b,1,ev,n,           &
     &               ilo,ihi,scale,abnrm,rconde,rcondv,work2,           &
     &               n*n+2*n,rwork,info)                                
         WRITE(*,*) "Cond.No:",abnrm 
         IF(info/=0)  CALL juDFT_error("cgeevx",calledby="gf_math.F90")
      ENDIF 
                                                                        
      END SUBROUTINE 

      !>                                                                
      !<-- S:rmat_eigenvalues(A,ew,ev)                                  
                                                                        
      SUBROUTINE rmat_eigenvalues(A,ew,ev) 
!     *****************************************                         
!     Use the LAPACK to get the eigenvalues ew                          
!     and the right eigenvectors ev of the real                         
!     square matrix A                                                   
!     *****************************************                         
      IMPLICIT NONE 
                                 !requires A_local                      
      REAL,INTENT(IN)  :: A(:,:) 
      REAL,INTENT(OUT) :: ew(:),ev(:,:) 
                                                                        
      INTEGER :: n,info 
      REAL    :: A_local(SIZE(A,1),SIZE(A,2)) 
      REAL    :: work(2*SIZE(a,1)) 
      REAL    :: rwork(2*SIZE(a,1)) 
                                                                        
      A_local = A 
      n = SIZE(A,1) 
      IF (SIZE(A,2) /= n) CALL juDFT_error                                &
     &     ('No square matrix in cmat_eigenvalues')                     
      IF (n /= SIZE(ew).OR.n /= SIZE(ev,1).OR.n /= SIZE(ev,2))          &
     &     CALL juDFT_error('Wrong dimensions in rmat_eigenvalues')       
                                                                        
      CALL dgeev('N','V',n,A_local,n,ew,                     &
     &     Ev,n,Ev,n,work,2*n,rwork,info)                               
      IF (info /= 0) CALL juDFT_error('rmat_eigenvalues failed') 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !>                                                                
                                                                        
      !<-- square root of a matrix                                      
      !<-- F: cmat_sqrt(A)                                              
                                                                        
      FUNCTION cmat_SQRT(A)RESULT(sqrtA) 
!-----------------------------------------------                        
!   Calculates the square-root of a hermitian matrix                    
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)     :: A(:,:) 
      COMPLEX                :: sqrtA(SIZE(A,1),SIZE(A,2)) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      COMPLEX             :: ew(SIZE(A,1)) 
      COMPLEX             :: U(SIZE(A,1),SIZE(A,2)) 
      !>                                                                
      IF (SIZE(A,1) /= SIZE(A,2)) THEN
       CALL juDFT_error("Non-square matrix in cmat_sqrt")
      ENDIF
      CALL cmat_eigenvalues(A,ew,U) 
      IF (ANY(AIMAG(ew)>1E-7)) THEN
       CALL juDFT_error("Non-Hermitian Matrix in cmat_sqrt")
      ENDIF
                                                                        
      sqrtA = 0.0 
      DO n  = 1,SIZE(A,1) 
         sqrtA(n,n) = SQRT(ew(n)) 
      ENDDO 
      sqrtA = MATMUL(U,MATMUL(sqrtA,mat_inverse(U))) 

      END FUNCTION 
                                                                        
      !>                                                                
      !<-- F: rmat_sqrt(A)                                              
                                                                        
      FUNCTION rmat_SQRT(A)RESULT(sqrtA) 
!-----------------------------------------------                        
!   Calculates the square-root of a hermitian matrix                    
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL,INTENT(IN)     :: A(:,:) 
      REAL                :: sqrtA(SIZE(A,1),SIZE(A,2)) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      REAL             :: ew(SIZE(A,1)) 
      REAL             :: U(SIZE(A,1),SIZE(A,2)) 
      !>                                                                
      IF (SIZE(A,1) /= SIZE(A,2)) THEN
               CALL juDFT_error("Non-square matrix in cmat_sqrt")
      ENDIF
      CALL rmat_eigenvalues(A,ew,U) 
                                                                        
      sqrtA = 0.0 
      DO n  = 1,SIZE(A,1) 
         sqrtA(n,n) = SQRT(ew(n)) 
      ENDDO 
      sqrtA = MATMUL(U,MATMUL(sqrtA,mat_inverse(U))) 
                                                                        
                                                                        
      END FUNCTION 
                                                                        
      !>                                                                
      !>                                                                
                                                                        
      !<-- Trace of a matrix                                            
      !<-- F:trace_complex(A)                                           
      FUNCTION trace_COMPLEX(A)RESULT(trace) 
!-----------------------------------------------                        
!  the trace of A                                                       
!             (last modified: 04-09-03) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)     :: A(:,:) 
      COMPLEX                :: trace 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: i 
      !>                                                                
      trace=0.0 
      DO i=1,size(A,1) 
         trace=trace+A(i,i) 
      ENDDO 
      END FUNCTION 
      !>                                                                
      !<-- F:trace_real(A)                                              
      FUNCTION trace_real(A)RESULT(trace) 
!-----------------------------------------------                        
!  the trace of A                                                       
!             (last modified: 04-09-03) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL,INTENT(IN)     :: A(:,:) 
      REAL                :: trace 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: i 
      !>                                                                
      trace=0.0 
      DO i=1,size(A,1) 
         trace=trace+A(i,i) 
      ENDDO 
      END FUNCTION 
      !>                                                                
      !>                                                                
                                                                        
                                                                        
      !misc stuff                                                       
      !<-- S:sort2(s1,s2,idx)                                           
      SUBROUTINE sort2(s1,s2,idx) 
!-----------------------------------------------                        
!   sort subroutine which uses second array for sorting                 
!   elements of same size                                               
!           (last modified:08-06-24) D. Wortmann                        
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     :: s1(:),s2(:) 
      INTEGER,INTENT(OUT)    :: idx(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n_same,n,nn 
      REAL                :: s(SIZE(idx)) 
      INTEGER             :: i(SIZE(idx)) 
      !>                                                                
                                                                        
      !do a standard sort                                               
      idx(1)=1 
      IF (SIZE(s1)<2) RETURN 
      CALL sort(s1,idx) 
      IF (SIZE(s1) == 1) RETURN 
      n = 1 
      nn = 1 
      loop:DO 
         IF (n>SIZE(idx)) EXIT loop 
         DO WHILE (s1(idx(n)) == s1(idx(nn))) 
            nn = nn+1 
            IF (nn>SIZE(idx)) EXIT loop 
         ENDDO 
         IF (nn-n == 1) THEN 
            n = n+1 
            nn = n 
            CYCLE loop 
         ENDIF 
         !ok there are elements of same length                          
         n_same = nn-n+1 
         s(:n_same) = s2(idx(n:nn)) 
         CALL sort(s(:n_same),i(:n_same)) 
         idx(n:nn) = idx(n-1+i(:n_same)) 
         n = nn 
      ENDDO loop 
                                                                        
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
      !<-- S:sort(s,idx)                                                
                                                                        
      SUBROUTINE sort(s,idx) 
!-----------------------------------------------                        
!     Generates an array idx containing sorted indices of s             
!     Compare sort from FLEUR                                           
!           (last modified: 05-03-22) D. Wortmann                       
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     :: s(:) 
      INTEGER,INTENT(OUT)    :: idx(:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,i,ind,ir,j,l 
      REAL                :: q 
      REAL,PARAMETER      :: eps         = 1.E-10 
      !>                                                                
      n = size(s) 
                                                                        
      IF (n /= SIZE(idx)) CALL juDFT_error                                &
     &     ("Wrong dimensions in gf_tools:sort")                        
      IF (n == 0) CALL juDFT_error                                        &
     &     ("Null array in gf_math:sort")                               
      idx = (/(i,i = 1,n)/) 
      l = n/2+1 
      ir = n 
      DO 
         IF (l>1) THEN 
            l = l - 1 
            ind = idx(l) 
            q   = s(ind) 
         ELSE 
            ind = idx(ir) 
            q   = s(ind) 
            idx(ir) = idx(1) 
            ir      = ir - 1 
            IF (ir  == 1) THEN 
               idx(1) = ind 
               RETURN 
            END IF 
         END IF 
         i = l 
         j = l + l 
         DO WHILE (j <= ir) 
            IF (j<ir) THEN 
               IF ((s(idx(j+1))-s(idx(j)))>eps) j = j + 1 
            END IF 
            IF ((s(idx(j))-q)>eps) THEN 
               idx(i) = idx(j) 
               i      = j 
               j      = j + j 
            ELSE 
               j = ir + 1 
            END IF 
         ENDDO 
         idx(i) = ind 
      ENDDO 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- F:imag2d(A,B)                                                
      FUNCTION imag2d(A,B) 
!-----------------------------------------------                        
!     Calculates the imaginary part of a matrix A(g,g')                 
!     which is the 2D  transform of a matrix A(r,r') = A(r',r)          
!     Optionally uses a second matrix B with A(r,r') =B(r',r)           
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)          :: A(:,:) 
      COMPLEX,INTENT(IN),OPTIONAL :: B(:,:) 
      COMPLEX                     :: imag2d(SIZE(A,1),SIZE(A,2)) 
      !>                                                                
      IF (PRESENT(B)) THEN 
         imag2d = CMPLX(0.0,-0.5)*(A-TRANSPOSE(CONJG(B))) 
      ELSE 
         imag2d = CMPLX(0.0,-0.5)*(A-TRANSPOSE(CONJG(A))) 
      ENDIF 
      END FUNCTION 
      !>                                                                
      !<-- F: interpolate_analytic(f,pos,x)                             
      FUNCTION interpolate_analytic(f,pos,x) 
!-----------------------------------------------                        
!    Interpolates analytic function given at positions pos with value f 
!    to point x                                                         
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(IN)     :: f(:),pos(:) 
      COMPLEX,INTENT(IN)     :: x 
      COMPLEX                :: interpolate_analytic 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      COMPLEX             :: fx,dfx 
      !>                                                                
      n = SIZE(f) 
      IF (n /= SIZE(pos)) CALL juDFT_error                                &
     &     ("Wrong sizes in interpolate_analytic")                      
      CALL priv_ratint(pos,f,n,x,fx) 
      interpolate_analytic = fx 
                                                                        
      CONTAINS 
      !<-- S: priv_ratint                                               
                                                                        
                                                                        
!     Given values of an anlytic function on N mesh points,             
!     interpolate the function to the entire z plane (Numerical recipes)
      SUBROUTINE priv_RATINT(XA,YA,N,X,Y)
      IMPLICIT NONE 
      INTEGER,INTENT(IN)  :: n 
      COMPLEX,INTENT(IN)  :: XA(n),YA(n),X 
      COMPLEX,INTENT(OUT) :: y 
      COMPLEX             :: C(N),D(N),w,ph,T,dd,dy 
      INTEGER             :: i,ns,m 
      REAL                :: HH,H 
      REAL,PARAMETER      :: TINY = 1.D-25 
      !<--Find position closest to X                                    
                             !nearest pos                               
      NS = 1 
                             !distance                                  
      HH = ABS(X-XA(1)) 
      DO  I = 1,N 
        H = ABS(X-XA(I)) 
        IF (H == 0.D0)THEN 
          Y = YA(I) 
          DY = (0.D0,0.D0) 
          RETURN 
        ELSE IF (H<HH) THEN 
          NS = I 
          HH = H 
        ENDIF 
        C(I) = YA(I) 
        D(I) = YA(I)+TINY 
      ENDDO 
      !>                                                                
      Y = YA(NS) 
      NS = NS-1 
      DO  M = 1,N-1 
         DO  I = 1,N-M 
            W = C(I+1)-D(I) 
            PH = XA(I+M)-X 
            T = (XA(I)-X)*D(I)/PH 
            DD = T-C(I+1) 
            IF(ABS(DD) <1000.*epsilon(1.0)) THEN 
               dd = HUGE(1.0)/1000.0 
            ELSE 
               DD = W/DD 
            ENDIF 
            D(I) = C(I+1)*DD 
            C(I) = T*DD 
         ENDDO 
         IF (2*NS<N-M)THEN 
            DY = C(NS+1) 
         ELSE 
            DY = D(NS) 
            NS = NS-1 
         ENDIF 
         Y = Y+DY 
      ENDDO 
      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      END FUNCTION 
      !>                                                                
      !<-- S: TRIANG(V,It,Ntria,At,Att,angle)                           

      SUBROUTINE TRIANG(V,It,Ntria,At,Att,angle) 
!-------------------------------------------------------------------    
!     find a triangular decomposition of the irreducible wedge of       
!     the first brillouin zone for a given k-mesh. k-points at all      
!     vertices are assumed.                                             
!     erich wimmer     july 1981                                        
!                                                                       
!     adapted 2D version for gfleur                                     
!     Be aware, this code only works for sufficiently regular k-grids   
!-------------------------------------------------------------------    
      IMPLICIT NONE 
!     Arguments                                                         
      REAL , INTENT(IN)        :: V(:,:) 
      REAL , INTENT(OUT)       :: Att,At(:),angle(:) 
      INTEGER , INTENT(OUT)    :: Ntria , It(:,:) 
!     locals                                                            
      REAL         :: a , a1 , d , dm , s0 , s1 
      INTEGER      :: j , j1 , j2 , jc , jj , k , k1 , k2 , k3 , kk 
      INTEGER      :: l1, l2 , l3 , n , nt0 , nkpts 
      LOGICAL      :: new 
!     constants                                                         
      REAL , PARAMETER :: BIG = 1.E8 , TOL = 1.E-5 
!                                                                       
      nkpts = SIZE(V,2) 
      IF ( Nkpts<3 ) RETURN 
!                                                                       
!     create first triangle                                             
!                                                                       
      Att = 0.0 
      k1  = 1 
      dm  = BIG 
      DO k = 2 , nkpts 
         d = VD2(V(:,k1),V(:,k)) 
         IF ( d<dm ) THEN 
            dm = d 
            k2 = k 
         ENDIF 
      ENDDO 
!                                                                       
      dm = BIG 
      new = .FALSE. 
      DO k = 1 , nkpts 
         IF ( k /= k1 .AND. k /= k2 ) THEN 
            d    = VD2(V(:,k1),V(:,k)) + VD2(V(:,k2),V(:,k)) 
            a    = AREA(V(:,k1),V(:,k2),V(:,k)) 
            IF ( ABS(a) >= TOL ) THEN 
               IF ( d<dm ) THEN 
                  new = .TRUE. 
                  dm  = d 
                  a1  = a 
                  k3 = k 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDDO 
      IF ( .NOT.new ) CALL juDFT_error("triang: Collinear kpts") 
      Ntria = 1 
      IF ( k3 <= k2 ) THEN 
         It(:,Ntria) = (/k1,k3,k2/) 
      ELSE 
         It(:,Ntria) = (/k1,k2,k3/) 
      ENDIF 
      angle(ntria) = min_angle(v(:,k1),v(:,k2),v(:,k3)) 
      At(Ntria)   = ABS(a1)/2 
      Att         = Att + At(Ntria) 
!                                                                       
!     create, if possible, a new triangle from each side of the mother  
!     triangle nt0 with minimal sum of sides                            
!                                                                       
      nt0 = 0 
      DO 
         nt0 = nt0 + 1 
         IF ( nt0>Ntria ) RETURN 
         DO l1 = 1 , 3 
            l2 = MOD(l1,3) + 1 
            l3 = MOD(l2,3) + 1 
            k1 = It(l1,nt0) 
            k2 = It(l2,nt0) 
            IF ( k2 <= k1 ) THEN 
               kk    = k1 
               k1    = k2 
               k2    = kk 
            ENDIF 
!---  >>     check if side k1-k2 belongs already to another triangle    
            new = .TRUE. 
            DO n = 1 , Ntria 
               IF ( n /= nt0 ) THEN 
                  IF ( (k1 == It(1,n) .AND. k2 == It(2,n)) .OR.         &
     &                 (k1 == It(1,n) .AND. k2 == It(3,n)) .OR.         &
     &                 (k1 == It(2,n) .AND. k2 == It(3,n)) )            &
     &                 new = .FALSE.                                    
               ENDIF 
            ENDDO 
            IF ( new ) THEN 
               k3 = It(l3,nt0) 
               a  = AREA(V(:,k1),V(:,k2),V(:,k3)) 
               s0 = SIGN(1.,a) 
!---  >>     a new triangle sharing the side k1-k2 with the mother      
!---  >>     triangle nt0 has to ly on the other side, i.e. the cross   
!---  >>     products (k2-k1)x(k3(old)-k1) and (k2-k1)x(k3(new)-k1)     
!---  >>     have to have opposite signs                                
               dm = BIG 
               new = .FALSE. 
               kloop:DO k = 1 , nkpts 
                  IF ( k /= k1 .AND. k/=k2 ) THEN 
!---  >>     check if a new side, (k1,k) or (k2,k), belongs             
!---  >>     already to an older triangle                               
                     j1 = k1 
                     j2 = k 
                     DO j = 1 , 2 
                        IF ( j2 <= j1 ) THEN 
                           jj  = j1 
                           j1  = j2 
                           j2  = jj 
                        ENDIF 
                        jc = 0 
                        DO n = 1 , Ntria 
                           IF ( j1==It(1,n) .AND. j2==It(2,n)         &
     &                          .OR. j1==It(1,n) .AND.                  &
     &                          j2==It(3,n) .OR. j1==It(2,n)          &
     &                          .AND. j2==It(3,n) ) jc = jc + 1         
                        ENDDO 
                        IF ( jc==2 ) CYCLE kloop 
                        j1  = k2 
                        j2  = k 
                     ENDDO 
                     a = AREA(V(:,k1),V(:,k2),V(:,k)) 
                     IF ( ABS(a)>=TOL ) THEN 
                        s1 = SIGN(1.,a) 
                        IF ( ABS(s1-s0)>=TOL ) THEN 
                           d = VD2(V(:,k1),V(:,k))+ VD2(V(:,k2),V(:,k)) 
                           IF ( d<dm ) THEN 
                              new = .TRUE. 
                              dm = d 
                              a1 = a 
                              k3 = k 
                           ENDIF 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDDO kloop 
               IF ( new ) THEN 
                  Ntria = Ntria + 1 
                  IF ( Ntria>SIZE(it,2)) CALL                           &
     &                 juDFT_error("Too many triangle in triang")         
                  IF ( k2<=k1 ) THEN 
                     kk   = k1 
                     k1 = k2 
                     k2 = kk 
                  ENDIF 
                  IF ( k3<=k1 ) THEN 
                     kk = k1 
                     k1 = k3 
                     k3 = kk 
                  ENDIF 
                  IF ( k3<=k2 ) THEN 
                     kk = k2 
                     k2 = k3 
                     k3 = kk 
                  ENDIF 
                  It(:,Ntria) = (/k1,k2,k3/) 
                  At(Ntria) = ABS(a1)/2 
                  Att    = Att + At(Ntria) 
                  angle(ntria) = min_angle(v(:,k1),v(:,k2),v(:,k3)) 
               ENDIF 
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
      CONTAINS 
!-----------------------------------------------------                  
                                   ! distance between vec1 = (x1,y1) and
      REAL FUNCTION VD2(vec1,vec2) 
      IMPLICIT NONE 
      REAL,INTENT(IN) :: vec1(2),vec2(2) 
      VD2 = (vec2(1)-vec1(1))**2 + (vec2(2)-vec1(2))**2 
      END FUNCTION 
                                          ! area of triangle (x1,y1),(x2
      REAL FUNCTION AREA(vec1,vec2,vec3) 
      IMPLICIT NONE 
      REAL,INTENT(IN) :: vec1(2),vec2(2),vec3(2) 
      AREA = (vec2(1)-vec1(1))*(vec3(2)-vec1(2)) - (vec3(1)-vec1(1))    &
     &     *(vec2(2)-vec1(2))                                           
      END FUNCTION 
                                               ! minimal angle of the tr
      REAL FUNCTION min_angle(vec1,vec2,vec3) 
      IMPLICIT NONE 
      REAL,INTENT(IN) :: vec1(2),vec2(2),vec3(2) 
      min_angle = ABS(ACOS(dot_PRODUCT(vec1,vec2)/SQRT(dot_PRODUCT(vec1 &
     &     ,vec1)*dot_PRODUCT(vec2,vec2))))                             
      min_angle = MIN(min_angle,ABS(ACOS(dot_PRODUCT(vec1,vec3)         &
     &     /SQRT(dot_PRODUCT(vec1,vec1)*dot_PRODUCT(vec3,vec3)))))      
      min_angle = MIN(min_angle,ABS(ACOS(dot_PRODUCT(vec2,vec3)         &
     &     /SQRT(dot_PRODUCT(vec2,vec2)*dot_PRODUCT(vec3,vec3)))))      
      END FUNCTION 
!-----------------------------------------------------                  
      END SUBROUTINE 
      SUBROUTINE cblas_matmul(a,b,c,transa,transb) 
!*************************************                                  
!     Call the blas routine zgemm in                                    
!     order to perform the product of                                   
!     matrices a and b.                                                 
!     Frank Freimuth, November 2007                                     
!*************************************                                  
      IMPLICIT NONE 
      COMPLEX,INTENT(IN) :: a(:,:) 
      COMPLEX,INTENT(IN)::b(:,:) 
      COMPLEX,INTENT(OUT)::c(:,:) 
                                                                        
      CHARACTER(LEN=1),INTENT(IN),OPTIONAL::transa 
      CHARACTER(LEN=1),INTENT(IN),OPTIONAL::transb 
      CHARACTER(LEN=1) transa2 
      CHARACTER(LEN=1) transb2 
      INTEGER m,n,k,lda,ldb,ldc,k2 
      lda=size(a,1) 
      ldb=size(b,1) 
      ldc=size(c,1) 
      IF(present(transa).AND.(transa=='C'.OR.transa=='c'.OR.        &
     &             transa=='T'.OR.transa=='t'))THEN                 
         m=size(a,2) 
         k=size(a,1) 
         transa2=transa 
      ELSE 
         transa2='N' 
         m=size(a,1) 
         k=size(a,2) 
      ENDIF 
      IF(present(transb).AND.(transb=='C'.OR.transb=='c'.OR.        &
     &             transb=='T'.OR.transb=='t'))THEN                 
         k2=size(b,2) 
         n=size(b,1) 
         transb2=transb 
      ELSE 
         transb2='N' 
         k2=size(b,1) 
         n=size(b,2) 
      ENDIF 
      IF(k2/=k)  CALL juDFT_error("cblas_matmul1",calledby="gf_math.F90")
      IF(size(c,1)/=m)  CALL juDFT_error("cblas_matmul2",calledby="gf_math.F90")
      IF(size(c,2)/=n)  CALL juDFT_error("cblas matmul3",calledby="gf_math.F90")
      CALL zgemm(transa,transb,m,n,k,                                   &
     &           cmplx(1.0,0.0),a,lda,b,ldb,                            &
     &           cmplx(0.0,0.0),c,ldc)                                  
      END SUBROUTINE cblas_matmul 

      !>                                                                
      !<-- S:  cblas_matmul_sq(a,b)                                     
      SUBROUTINE cblas_matmul_sq(a,b) 
!*********************************************                          
!     Multiply the square matrices a and b.                             
!     The matrix product a*b is written to a.                           
!     Frank Freimuth, November 2007                                     
!*********************************************                          
      IMPLICIT NONE 
      COMPLEX,INTENT(INOUT)::a(:,:) 
      COMPLEX,INTENT(IN)::b(:,:) 
      INTEGER m 
      COMPLEX,ALLOCATABLE::c(:,:) 
      m=size(a,1) 
      IF(size(a,2)/=m)  CALL juDFT_error("matmul_sq1",calledby="gf_math.F90")
      IF(size(b,1)/=m)  CALL juDFT_error("matmul_sq2",calledby="gf_math.F90")
      IF(size(b,2)/=m)  CALL juDFT_error("matmul_sq3",calledby="gf_math.F90")
      ALLOCATE( c(m,m) ) 
      CALL zgemm('N','N',m,m,m,cmplx(1.0,0.0),                          &
     &            a,m,b,m,cmplx(0.0,0.0),c,m)                           
      a=c 
      DEALLOCATE( c ) 
      END SUBROUTINE cblas_matmul_sq 
      !>                                                                
      !<-- S: cblas_matmul_sq2(a,b,c)                                   
      SUBROUTINE cblas_matmul_sq2(a,b,c) 
!*********************************************                          
!     Multiply the square matrices a and b.                             
!     The matrix product a*b is written to c.                           
!     Frank Freimuth, November 2007                                     
!*********************************************                          
      IMPLICIT NONE 
      COMPLEX,INTENT(IN)::a(:,:) 
      COMPLEX,INTENT(IN)::b(:,:) 
      COMPLEX,INTENT(OUT)::c(:,:) 
      INTEGER m 
      m=size(a,1) 
      CALL zgemm('N','N',m,m,m,cmplx(1.0,0.0),                          &
     &            a,m,b,m,cmplx(0.0,0.0),c,m)                           
      END SUBROUTINE cblas_matmul_sq2 
      !>                                                                
      !<--F: zmat_product(a,b)                                          
      FUNCTION zmat_product(a,b) 
!*********************************************                          
!     Multiply the matrices a and b.                                    
!     Frank Freimuth, November 2007                                     
!*********************************************                          
      IMPLICIT NONE 
      COMPLEX,INTENT(IN)::a(:,:) 
      COMPLEX,INTENT(IN)::b(:,:) 
      COMPLEX zmat_product(size(a,1),size(a,1)) 
      INTEGER m 
      m=size(a,1) 
      CALL zgemm('N','N',m,m,m,cmplx(1.0,0.0),                          &
     &            a,m,b,m,cmplx(0.0,0.0),zmat_product,m)                
      END FUNCTION 
      !>                                                                
                                                                        
      !<-- F: simple_broyden(x,y)                                       
      FUNCTION simple_broyden(x,y) 
!-----------------------------------------------                        
!   very simple (and slow) implementation of broydens method            
!   Input: x (:,:n_iter)                                                
!          y (:,:n_iter)                                                
!                                                                       
!    Solves for y=f(x)=0                                                
!                                                                       
!           (last modified:09-11-30) D. Wortmann                        
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<-- Arguments                                                    
      REAL   ,INTENT(IN)       :: x(:,:) 
      REAL   ,INTENT(IN)       :: y(:,:) 
      REAL                     :: simple_broyden(SIZE(x,1)) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      REAL                :: IJ(SIZE(y,1),SIZE(y,1)) 
      REAL                :: norm 
      REAL                :: dx(SIZE(y,1)),df(SIZE(y,1)) 
      REAL                :: left(SIZE(y,1),1) 
      REAL                :: right(1,SIZE(y,1)) 
      !>                                                                
      IF ((SIZE(y,1)/=SIZE(x,1)).OR.(SIZE(x,2)/=SIZE(y,2))          &
     &     ) CALL juDFT_error                                             &
     &     ("Bug in gf_math:simple_broyden")                            
                                                                        
      !<-- Initial inverse Jacobian = 1                                 
      IJ = 0.0 
      DO n = 1,SIZE(y,1) 
         IJ(n,n) = 1.0 
      ENDDO 
      !>                                                                
      DO n = 2,SIZE(y,1) 
         dx = x(:,n)-x(:,n-1) 
         df = y(:,n)-y(:,n-1) 
         norm = dot_PRODUCT(dx,MATMUL(IJ,df)) 
         left(:,1) = dx-MATMUL(IJ,df) 
         right(1,:) = MATMUL(dx,IJ) 
         IJ = IJ+MATMUL(left,right)/norm 
      ENDDO 
                                                                        
      simple_broyden = x(:,SIZE(x,2))-MATMUL(IJ,y(:,SIZE(x,2))) 
                                                                        
      END FUNCTION 
      !>                                                                
#if (1==2)
      subroutine make_charge2D(mat,lapw,stars,n)
      !-----------------------------------------------
      !   Take a 2d matrix in lapw-basis and return the charge in stars
      !-----------------------------------------------
      use m_types
      implicit none
      type(t_lapw),intent(in) :: lapw
      type(t_stars),intent(in):: stars
      complex,intent(in)      :: mat(:,:)
      complex,intent(out)     :: n(:)

      integer :: n1,n2,ng
      integer :: n3_star,n2_star

      n=0.0
      DO ng=1,stars%nq2
         DO n1=1,lapw%nv2(1)
            DO n2=1,lapw%nv2(1)
               !determine if n1-n2 are in star ng
               n3_star=stars%ig(lapw%k%kp1(n1,1)-lapw%k%kp1(n2,1),lapw%k%kp2(n1,1)-lapw%k%kp2(n2,1),0)
               if (n3_star/=0) then
                 if (stars%ig2(n3_star)==ng) n(ng)=n(ng)+mat(n1,n2)
               endif
               !determine if n2-n1 are in star ng
               n3_star=stars%ig(lapw%k%kp1(n2,1)-lapw%k%kp1(n1,1),lapw%k%kp2(n2,1)-lapw%k%kp2(n1,1),0)
               if (n3_star/=0) then
                 if (stars%ig2(n3_star)==ng) n(ng)=n(ng)-conjg(mat(n1,n2))
               endif
           enddo
        enddo
      enddo
      n=n*cmplx(0,-0.5)
      end subroutine
#endif

                                                                        
      END                                           
