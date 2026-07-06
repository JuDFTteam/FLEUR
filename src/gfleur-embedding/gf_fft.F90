!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_fft 
      use m_juDFT
      use m_juDFT 
      IMPLICIT NONE
!*****************************************************************      
!     Easy to use interfaces for fft subroutines                        
!     The FLEUR cfft subroutines are not used anymore but instead       
!     the newer f90 fft routines from the module included in this file  
!     Daniel Wortmann, Thu Oct 16 10:33:01 2003                         
!*****************************************************************      
      PRIVATE 
      INTEGER,PARAMETER ::GF_FFT_TO_G_SPACE=-1                          &
     &     ,GF_FFT_TO_REAL_SPACE=1                                      
      INTERFACE gf_fft3d 
      MODULE PROCEDURE gf_fft3dc,gf_fft3dr 
      END INTERFACE 
      INTERFACE gf_fft2d 
      MODULE PROCEDURE gf_fft2dr,gf_fft2dc 
      END INTERFACE 
      PUBLIC  gf_fft2d,gf_fft3d,GF_FFT_TO_G_SPACE,GF_FFT_TO_REAL_SPACE  &
     &     ,gf_fft23d,gf_fft2d_small,gf_convol,cfft                     
      CONTAINS 
      !<--S:gf_fft3dr(r,g,stars,dir)                                    
                                                                        
      SUBROUTINE gf_fft3dr(r,g,stars,dir) 
!******************************************                             
!     real to complex fft                                               
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL                     :: r(:) 
      COMPLEX,INTENT(INOUT)    :: g(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      INTEGER,INTENT(IN)       :: dir 
      !>                                                                
      !<--Locals                                                        
                                                                        
      REAL ::ri(0:SIZE(r)-1,2) 
      INTEGER :: ifftd,i 
                                                                        
      !>                                                                
      !<-- Tests                                                        

      IF(SIZE(r)/=27*stars%mx1*stars%mx2*stars%mx3) CALL juDFT_error      &
     &     ('Dimension of real-mesh wrong in gf_fft3dr')                
      IF(SIZE(g)/=stars%nq3) CALL juDFT_error                           &
     &     ('Dimension of star-array wrong in gf_fft3dr')               
      IF(abs(dir)/=1)  CALL juDFT_error                                 &
     &     ('Invalid direction specified in gf_fft3dr')                 

      !>                                                                
      ifftd    = 27*stars%mx1*stars%mx2*stars%mx3 
                                                                        
                                                                        
      IF(dir==GF_FFT_TO_G_SPACE) THEN 
         ri(:,1)=r 
         ri(:,2)=0 
      ELSE 
         !<-- put stars onto the fft-grid                               
         ri=0.0 
         DO i = 0,stars%kimax-1 
          ri(stars%igfft(i,2),1)=real(g(stars%igfft(i,1)))*             &
     &                                   stars%pgfft(i)                 
          ri(stars%igfft(i,2),2)=aimag(g(stars%igfft(i,1)))*            &
     &                                   stars%pgfft(i)                 
        ENDDO 
        !>                                                              
      ENDIF 
                                                                        
      !<--now do the fft (dir=+1 : G -> r ; dir=-1 : r -> G)            
                                                                        
      CALL cfft(ri(:,1),ri(:,2),ifftd,3*stars%mx1,3*stars%mx1,dir) 
      CALL cfft(ri(:,1),ri(:,2),ifftd,3*stars%mx2,9*stars%mx1*stars%mx2 &
     &     ,dir)                                                        
      CALL cfft(ri(:,1),ri(:,2),ifftd,3*stars%mx3,ifftd,dir) 
                                                                        
      !>                                                                
      IF (dir==GF_FFT_TO_G_SPACE) THEN 
                                                                        
         !<-- collect stars from the fft-grid                           
                                                                        
         g = CMPLX(0.0,0.0) 
         DO i = 0,stars%kimax-1 
            g(stars%igfft(i,1)) = g(stars%igfft(i,1))+stars%pgfft(i)*   &
     &           CMPLX(ri(stars%igfft(i,2),1),ri(stars%igfft(i,2),2))   
        ENDDO 
        g = g/stars%nstr/ifftd 
                                                                        
        !>                                                              
      ELSE 
         r = ri(:,1) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<--S:gf_fft3dc(r,g,stars,dir)                                    
                                                                        
      SUBROUTINE gf_fft3dc(r,g,stars,dir) 
!******************************************                             
!     COMPLEX to COMPLEX fft-code                                       
!     code addaped from fft3d of fleur                                  
!     D. Wortmann                                                       
!******************************************                             
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(INOUT)    :: r(:) 
      COMPLEX,INTENT(INOUT)    :: g(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      INTEGER,INTENT(IN)       :: dir 
      !>                                                                
      !<--Locals                                                        
                                                                        
      REAL    :: ri(0:SIZE(r)-1,2) 
      INTEGER :: ifftd,i 
                                                                        
      !>                                                                
      !<-- Tests                                                        

      IF(SIZE(r) /= 27*stars%mx1*stars%mx2*stars%mx3) CALL juDFT_error(   &
     &     'Dimension of real-mesh wrong in gf_fft3dc')                 
      IF(SIZE(g) /= stars%nq3) CALL juDFT_error(                          &
     &     'Dimension of star-array wrong in gf_fft3dc')                
      IF(abs(dir) /= 1) CALL juDFT_error(                                 &
     &     'Invalid direction specified in gf_fft3dc')                  

      !>                                                                
      ifftd = 27*stars%mx1*stars%mx2*stars%mx3 
                                                                        
                                                                        
      IF(dir==GF_FFT_TO_G_SPACE) THEN 
         ri(:,1)=REAL(r) 
         ri(:,2)=AIMAG(r) 
      ELSE 
         !<-- put stars onto the fft-grid                               
         ri=0.0 
         DO i = 0,stars%kimax-1 
          ri(stars%igfft(i,2),1)=real(g(stars%igfft(i,1)))*             &
     &                                   stars%pgfft(i)                 
          ri(stars%igfft(i,2),2)=aimag(g(stars%igfft(i,1)))*            &
     &                                   stars%pgfft(i)                 
        ENDDO 
        !>                                                              
      ENDIF 
                                                                        
      !<--now do the fft (dir=+1 : G -> r ; dir=-1 : r -> G)            
                                                                        
      CALL cfft(ri(:,1),ri(:,2),ifftd,3*stars%mx1,3*stars%mx1,dir) 
      CALL cfft(ri(:,1),ri(:,2),ifftd,3*stars%mx2,9*stars%mx1*stars%mx2 &
     &     ,dir)                                                        
      CALL cfft(ri(:,1),ri(:,2),ifftd,3*stars%mx3,ifftd,dir) 
                                                                        
      !>                                                                
      IF (dir==GF_FFT_TO_G_SPACE) THEN 
                                                                        
         !<-- collect stars from the fft-grid                           
         g = CMPLX(0.0,0.0) 
         DO i = 0,stars%kimax-1 
            g(stars%igfft(i,1)) = g(stars%igfft(i,1))+stars%pgfft(i)*   &
     &           CMPLX(ri(stars%igfft(i,2),1),ri(stars%igfft(i,2),2))   
        ENDDO 
        g = g/stars%nstr/ifftd 
        !>                                                              
      ELSE 
         r = CMPLX(ri(:,1),ri(:,2)) 
      ENDIF 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S:gf_fft23d(r,g,stars)                                        
                                                                        
      SUBROUTINE gf_fft23d(r,g,stars) 
!******************************************                             
!     COMPLEX to COMPLEX fft,                                           
!     only the x und y components are                                   
!     transformed to real space!                                        
!     D. Wortmann                                                       
!******************************************                             
      USE m_gf_fft_singleton,ONLY:fft 
      USE m_gf_types 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      COMPLEX,INTENT(INOUT) :: r(:) 
      COMPLEX,INTENT(INOUT)    :: g(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX   :: ri(0:SIZE(r)-1) 
      INTEGER ::i 
      !>                                                                
      !<-- Tests                                                        

      IF(SIZE(r) /= 27*stars%mx1*stars%mx2*stars%mx3) CALL juDFT_error(   &
     &     'Dimension of real-mesh wrong in gf_fft3dc')                 
      IF(SIZE(g) /= stars%nq3) CALL juDFT_error(                          &
     &     'Dimension of star-array wrong in gf_fft3dc')                

      !>                                                                
      ri=0.0 
      DO i=0,stars%kimax-1 
         ri(stars%igfft(i,2))=g(stars%igfft(i,1))*stars%pgfft(i) 
      ENDDO 
                                                                        
      r = reshape(fft(reshape(ri,(/3*stars%mx1,3*stars%mx2,3*stars%mx3/)&
     &     ),(/1,2/)),(/size(r)/))                                      
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<--S:gf_fft2dr(r,g,stars,dir)                                    
                                                                        
      SUBROUTINE gf_fft2dr(r,g,stars,dir) 
!******************************************                             
!     real to complex fft                                               
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(INOUT)    :: r(:) 
      COMPLEX,INTENT(INOUT)    :: g(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      INTEGER,INTENT(IN)       :: dir 
      !>                                                                
      !<--Locals                                                        
      REAL    ::ri(0:SIZE(r)-1,2) 
      INTEGER :: i 
      !>                                                                
      !<-- Tests                                                        

      IF(SIZE(r) /= 9*stars%mx1*stars%mx2) CALL juDFT_error(              &
     &     'Dimension of real-mesh wrong in gf_fft2dr')                 
      IF(SIZE(g)/=stars%nq2) CALL juDFT_error(                            &
     &     'Dimension of star-array wrong in gf_fft2dr')                
      IF(abs(dir)/=1) CALL juDFT_error(                                 &
     &     'Invalid direction specified in gf_fft2dr')                  

      !>                                                                
!                                                                       
      IF (dir==GF_FFT_TO_REAL_SPACE) THEN 
         !<-- put stars onto the fft-grid                               
         ri = 0.0 
         DO i = 0,stars%kimax2-1 
            ri(stars%igfft2(i,2),1) = REAL(g(stars%igfft2(i,1)))*       &
     &           stars%pgfft2(i)                                        
            ri(stars%igfft2(i,2),2) = AIMAG(g(stars%igfft2(i,1)))*      &
     &           stars%pgfft2(i)                                        
         ENDDO 
         !>                                                             
      ELSE 
         ri(:,1) = r 
      ENDIF 
                                                                        
      !<--now do the fft (dir =+1 : G -> r ; dir =-1 : r -> G)          
                                                                        
      CALL cfft(ri(:,1),ri(:,2),9*stars%mx1*stars%mx2,3*stars%mx1,3     &
     &     *stars%mx1,dir)                                              
      CALL cfft(ri(:,1),ri(:,2),9*stars%mx1*stars%mx2,3*stars%mx2,9     &
     &     *stars%mx1*stars%mx2,dir)                                    
                                                                        
      !>                                                                
      IF (dir == GF_FFT_TO_G_SPACE) THEN 
         !<-- collect stars from the fft-grid                           
         g =  CMPLX(0.0,0.0) 
         DO i = 0,stars%kimax2-1 
            g(stars%igfft2(i,1)) = g(stars%igfft2(i,1))+stars%pgfft2(i)*&
     &           CMPLX(ri(stars%igfft2(i,2),1),ri(stars%igfft2(i,2),2)) 
         ENDDO 
         g = g/stars%nstr2/9/stars%mx1/stars%mx2 
         !>                                                             
      ELSE 
         r = ri(:,1) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_fft2dc(r,g,stars,dir)                                   
                                                                        
      SUBROUTINE gf_fft2dc(r,g,stars,dir) 
!******************************************                             
!     real to complex fft                                               
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(INOUT)    :: r(:) 
      COMPLEX,INTENT(INOUT)    :: g(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      INTEGER,INTENT(IN)       :: dir 
      !>                                                                
      !<--Locals                                                        
      REAL ::ri(0:SIZE(r)-1,2) 
      INTEGER :: i 
      !>                                                                
                                                                        
      !<--Tests                                                         

      IF(SIZE(r)/=9*stars%mx1*stars%mx2) CALL juDFT_error(              &
     &     'Dimension of real-mesh wrong in gf_fft2dc')                 
      IF(SIZE(g)/=stars%nq2) CALL juDFT_error(                          &
     &     'Dimension of star-array wrong in gf_fft2dc')                
      IF(abs(dir)/=1) CALL juDFT_error(                                 &
     &     'Invalid direction specified in gf_fft2dc')                  

      !>                                                                
                                                                        
      !<-- Tests                                                        

      IF(SIZE(r) /= 9*stars%mx1*stars%mx2) CALL juDFT_error(              &
     &     'Dimension of real-mesh wrong in gf_fft2dr')                 
      IF(SIZE(g)/=stars%nq2) CALL juDFT_error(                            &
     &     'Dimension of star-array wrong in gf_fft2dr')                
      IF(abs(dir)/=1) CALL juDFT_error(                                 &
     &     'Invalid direction specified in gf_fft2dr')                  

      !>                                                                
!                                                                       
      IF (dir==GF_FFT_TO_REAL_SPACE) THEN 
         !<-- put stars onto the fft-grid                               
         ri = 0.0 
         DO i = 0,stars%kimax2-1 
            ri(stars%igfft2(i,2),1) = REAL(g(stars%igfft2(i,1)))*       &
     &           stars%pgfft2(i)                                        
            ri(stars%igfft2(i,2),2) = AIMAG(g(stars%igfft2(i,1)))*      &
     &           stars%pgfft2(i)                                        
         ENDDO 
         !>                                                             
      ELSE 
         ri(:,1) = REAL(r) 
         ri(:,2) = AIMAG(r) 
      ENDIF 
                                                                        
      !<--now do the fft (dir =+1 : G -> r ; dir =-1 : r -> G)          
                                                                        
      CALL cfft(ri(:,1),ri(:,2),9*stars%mx1*stars%mx2,3*stars%mx1,3     &
     &     *stars%mx1,dir)                                              
      CALL cfft(ri(:,1),ri(:,2),9*stars%mx1*stars%mx2,3*stars%mx2,9     &
     &     *stars%mx1*stars%mx2,dir)                                    
                                                                        
      !>                                                                
      IF (dir == GF_FFT_TO_G_SPACE) THEN 
         !<-- collect stars from the fft-grid                           
         g =  CMPLX(0.0,0.0) 
         DO i = 0,stars%kimax2-1 
            g(stars%igfft2(i,1)) = g(stars%igfft2(i,1))+stars%pgfft2(i)*&
     &           CMPLX(ri(stars%igfft2(i,2),1),ri(stars%igfft2(i,2),2)) 
         ENDDO 
         g = g/stars%nstr2/9/stars%mx1/stars%mx2 
         !>                                                             
      ELSE 
         r = CMPLX(ri(:,1),ri(:,2)) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_fft2d_small(r,g,stars,dir)                              

      SUBROUTINE gf_fft2d_small(r,g,stars,dir) 
!-----------------------------------------------                        
!     complex to complex FFT with small real-space grid                 
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      COMPLEX,INTENT(INOUT)    :: r(:) 
      COMPLEX,INTENT(INOUT)    :: g(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      INTEGER,INTENT(IN)       :: dir 
      !>                                                                
      !<--Locals                                                        
      REAL ::grid(0:SIZE(r)-1,2) 
      INTEGER :: griddim,x,y,index 
      !>                                                                
      griddim = SIZE(r) 
      IF(dir==GF_FFT_TO_REAL_SPACE) THEN 
         !Put stars on fft-grid                                         
         IF (ANY(stars%pgfft<0.0)) CALL juDFT_error("Symmetry problem") 
         grid=0.0 
         DO x=0,2*stars%mx1-1 
            DO y = 0,2*stars%mx2-1 
               IF (x>stars%mx1) THEN 
                  IF (y>stars%mx2) THEN 
                     index = stars%ig(x-2*stars%mx1,y-2*stars%mx2,0) 
                  ELSE 
                     index = stars%ig(x-2*stars%mx1,y,0) 
                  ENDIF 
               ELSE 
                  IF (y>stars%mx2) THEN 
                     index = stars%ig(x,y-2*stars%mx2,0) 
                  ELSE 
                     index = stars%ig(x,y,0) 
                  ENDIF 
               ENDIF 
               IF (index<1) CYCLE 
               index = stars%ig2(index) 
               grid(x+2*stars%mx1*y,1) = REAL(g(index)) 
               grid(x+2*stars%mx1*y,2) = AIMAG(g(index)) 
            ENDDO 
         ENDDO 
      ELSE 
         grid(0:,1) = REAL(r) 
         grid(0:,2) = AIMAG(r) 
      ENDIF 
                                                                        
                                                                        
      CALL cfft(grid(:,1),grid(:,2),griddim,2*stars%mx1,2*stars%mx1,dir) 
      CALL cfft(grid(:,1),grid(:,2),griddim,2*stars%mx2,griddim,dir) 
                                                                        
      IF(dir==GF_FFT_TO_G_SPACE) THEN 
         !Collect stars                                                 
         g = 0.0 
         DO x = 0,2*stars%mx1-1 
            DO y = 0,2*stars%mx2-1 
               IF (x>stars%mx1) THEN 
                  IF (y>stars%mx2) THEN 
                     index = stars%ig(x-2*stars%mx1,y-2*stars%mx2,0) 
                  ELSE 
                     index = stars%ig(x-2*stars%mx1,y,0) 
                  ENDIF 
               ELSE 
                  IF (y>stars%mx2) THEN 
                     index = stars%ig(x,y-2*stars%mx2,0) 
                  ELSE 
                     index = stars%ig(x,y,0) 
                  ENDIF 
               ENDIF 
               IF (index<1) CYCLE 
               index = stars%ig2(index) 
                                                                        
               g(index) = cmplx(grid(x+2*stars%mx1*y,1),grid(x+2        &
     &              *stars%mx1*y,2))+g(index)                           
            ENDDO 
         ENDDO 
         g=g/real(griddim)/stars%nstr2 
      ELSE 
         r=CMPLX(grid(:,1),grid(:,2)) 
      ENDIF 
      END SUBROUTINE 

      !>                                                                
                                                                        
      !<-- S: gf_convol(f,ufft,stars,f_out)                             
      SUBROUTINE gf_convol(f_in,ufft,stars,f_out) 
!-----------------------------------------------                        
!  convolutes the g-space function f with the real space function ufft  
!  see convol from FLEUR                                                
!           (last modified: 05-03-23) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                          !This is intent(in) only!     
      COMPLEX,INTENT(INOUT)    :: f_in(:) 
      COMPLEX,INTENT(OUT)      :: f_out(:) 
      REAL   ,INTENT(IN)       :: ufft(:) 
      TYPE(t_stars),INTENT(IN) :: stars 
      !>                                                                
      !<-- Locals                                                       
      REAL :: r_fft(size(ufft)) 
      !>                                                                
                                                                        
      CALL gf_fft3d(r_fft,f_in,stars,GF_FFT_TO_REAL_SPACE) 
                                                                        
      r_fft = r_fft*ufft 
                                                                        
      CALL gf_fft3d(r_fft,f_out,stars,GF_FFT_TO_G_SPACE) 
                               !needed for compatibility with FLEURs def
      f_out = f_out*stars%nstr 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: cfft(a,b,ntot,n,nspan,isn)                                
                                                                        
      SUBROUTINE cfft(a,b,ntot,n,nspan,isn) 
!-----------------------------------------------                        
!   fft-subroutine from FLEUR                                           
!           (last modified: 05-03-23) D. Wortmann                       
!-----------------------------------------------                        
!     ***************************************************************   
!     multivariate complex fourier transform, computed in place         
!     using mixed-radix fast fourier transform algorithm.               
!     by r. c. singleton, stanford research institute, oct. 1968        
!     arrays a and b originally hold the real and imaginary             
!     components of the data, and return the real and                   
!     imaginary components of the resulting fourier coefficients.       
!     multivariate data is indexed according to the fortran             
!     array element successor function, without limit                   
!     on the number of implied multiple subscripts.                     
!     the subroutine is called once for each variate.                   
!     the calls for a multivariate transform may be in any order.       
!     ntot is the total number of complex data values.                  
!     n is the dimension of the current variable.                       
!     nspan/n is the spacing of consucutive data values                 
!     while indexing the current variable.                              
!     the sign of isn determines the sign of the complex                
!     exponential, and the magnitude of isn is normally one.            
!     for a single-variate transform,                                   
!     ntot = n = nspan = (number of complex data values), f.g.          
!     call cft(a,b,n,n,n,1)                                             
!     a tri-variate transform with a(n1,n2,n3), b(n1,n2,n3)             
!     is computed by                                                    
!     call cft(a,b,n1*n2*n3,n1,n1,1)                                    
!     call cft(a,b,n1*n2*n3,n2,n1*n2,1)                                 
!     call cft(a,b,n1*n2*n3,n3,n1*n2*n3,1)                              
!     the data may alternatively be stored in a single complex          
!     array a, then the magnitude of isn changed to two to              
!     give the correct indexing increment and the second parameter      
!     used to pass the initial address for the sequence of              
!     imaginary values, e.g.                                            
!        real s(2)                                                      
!        equivalence (a,s)                                              
!        ....                                                           
!        ....                                                           
!        call cft(a,s(2),ntot,n,nspan,2)                                
!     arrays at(maxf), ck(maxf), bt(maxf), sk(maxf), and np(maxp)       
!     are used for temporary storage. if the available storage          
!     is insufficient, the program is terminated by a stop.             
!     maxf must be .ge. the maximum prime factor of n.                  
!     maxp must be .gt. the number of prime factors of n.               
!     in addition, if the square-free portion k cf n has two or         
!     more prime factors, then maxp must be .ge. k-1.                   
!     array storage in nfac for a maximum of 11 factors of n.           
!     if n has more than one square-free factor, the product of the     
!     square-free factors must be .le. 210                              
!     ******************************************************************
!     array storage for maximum prime factor of 199                     
!     the following two constants should agree with the array dimensions
                                                                        
      IMPLICIT NONE 
!     .. Scalar Arguments ..                                            
      INTEGER isn,n,nspan,ntot 
!     ..                                                                
!     .. Array Arguments ..                                             
      REAL a(*),b(*) 
!     ..                                                                
!     .. Local Scalars ..                                               
      REAL aa,aj,ajm,ajp,ak,akm,akp,bb,bj,bjm,bjp,bk,bkm,bkp,c1,c2,c3,  &
     &     c72,cd,rad,radf,s1,s120,s2,s3,s72,sd                         
      INTEGER i,ii,inc,j,jc,jf,jj,k,k1,k2,k3,k4,kk,ks,kspan,kspnn,kt,m, &
     &        maxf,maxp,nn,nt                                           
!     ..                                                                
!     .. Local Arrays ..                                                
      REAL at(199),bt(199),ck(199),sk(199) 
      INTEGER nfac(11),np(503) 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
      INTRINSIC cos,real,mod,sin 
!     ..                                                                
!     .. Equivalences ..                                                
      EQUIVALENCE (i,ii) 
!     ..                                                                
      maxf = 199 
      maxp = 503 
      IF (n<2) RETURN 
      inc = isn 
!     the following constants are rad = 2.*pi , s72 = sin(0.4*pi) ,     
!     c72 = cos(0.4*pi) and s120 = sqrt(0.75)                           
      rad = 6.2831853071796 
      s72 = 0.95105651629515 
      c72 = 0.30901699437495 
      s120 = 0.86602540378444 
      IF (isn>=0) GO TO 10 
      s72 = -s72 
      s120 = -s120 
      rad = -rad 
      inc = -inc 
   10 nt = inc*ntot 
      ks = inc*nspan 
      kspan = ks 
      nn = nt - inc 
      jc = ks/n 
      radf = rad*real(jc)*0.5 
      i = 0 
      jf = 0 
!     determine the factors of n                                        
      m = 0 
      k = n 
      GO TO 30 
   20 m = m + 1 
      nfac(m) = 4 
      k = k/16 
   30 IF (k- (k/16)*16==0) GO TO 20 
      j = 3 
      jj = 9 
      GO TO 50 
   40 m = m + 1 
      nfac(m) = j 
      k = k/jj 
   50 IF (mod(k,jj)==0) GO TO 40 
      j = j + 2 
      jj = j**2 
      IF (jj<=k) GO TO 50 
      IF (k>4) GO TO 60 
      kt = m 
      nfac(m+1) = k 
      IF (k/=1) m = m + 1 
      GO TO 100 
   60 IF (k- (k/4)*4/=0) GO TO 70 
      m = m + 1 
      nfac(m) = 2 
      k = k/4 
   70 kt = m 
      j = 2 
   80 IF (mod(k,j)/=0) GO TO 90 
      m = m + 1 
      nfac(m) = j 
      k = k/j 
   90 j = ((j+1)/2)*2 + 1 
      IF (j<=k) GO TO 80 
  100 IF (kt==0) GO TO 120 
      j = kt 
  110 m = m + 1 
      nfac(m) = nfac(j) 
      j = j - 1 
      IF (j/=0) GO TO 110 
!     compute fourier transform                                         
  120 sd = radf/real(kspan) 
      cd = 2.0*sin(sd)**2 
      sd = sin(sd+sd) 
      kk = 1 
      i = i + 1 
      IF (nfac(i)/=2) GO TO 170 
!     transform for factor of 2 (including rotation factor)             
      kspan = kspan/2 
      k1 = kspan + 2 
  130 k2 = kk + kspan 
      ak = a(k2) 
      bk = b(k2) 
      a(k2) = a(kk) - ak 
      b(k2) = b(kk) - bk 
      a(kk) = a(kk) + ak 
      b(kk) = b(kk) + bk 
      kk = k2 + kspan 
      IF (kk<=nn) GO TO 130 
      kk = kk - nn 
      IF (kk<=jc) GO TO 130 
      IF (kk>kspan) GO TO 360 
  140 c1 = 1.0 - cd 
      s1 = sd 
  150 k2 = kk + kspan 
      ak = a(kk) - a(k2) 
      bk = b(kk) - b(k2) 
      a(kk) = a(kk) + a(k2) 
      b(kk) = b(kk) + b(k2) 
      a(k2) = c1*ak - s1*bk 
      b(k2) = s1*ak + c1*bk 
      kk = k2 + kspan 
      IF (kk<nt) GO TO 150 
      k2 = kk - nt 
      c1 = -c1 
      kk = k1 - k2 
      IF (kk>k2) GO TO 150 
      ak = c1 - (cd*c1+sd*s1) 
      s1 = (sd*c1-cd*s1) + s1 
!     the following three statements compensate for truncation          
!     error. if rounded arithmetic is used, they may be deleted.        
!     c1=0.5/(ak**2+s1**2)+0.5                                          
!     s1=c1*s1                                                          
!     c1=c1*ak                                                          
!     next statement should be deleted if non-rounded arithmetic is used
      c1 = ak 
      kk = kk + jc 
      IF (kk<k2) GO TO 150 
      k1 = k1 + inc + inc 
      kk = (k1-kspan)/2 + jc 
      IF (kk<=jc+jc) GO TO 140 
      GO TO 120 
!     transform for factor of 3 (optional code)                         
  160 k1 = kk + kspan 
      k2 = k1 + kspan 
      ak = a(kk) 
      bk = b(kk) 
      aj = a(k1) + a(k2) 
      bj = b(k1) + b(k2) 
      a(kk) = ak + aj 
      b(kk) = bk + bj 
      ak = -0.5*aj + ak 
      bk = -0.5*bj + bk 
      aj = (a(k1)-a(k2))*s120 
      bj = (b(k1)-b(k2))*s120 
      a(k1) = ak - bj 
      b(k1) = bk + aj 
      a(k2) = ak + bj 
      b(k2) = bk - aj 
      kk = k2 + kspan 
      IF (kk<nn) GO TO 160 
      kk = kk - nn 
      IF (kk<=kspan) GO TO 160 
      GO TO 320 
!     transform for factor of 4                                         
  170 IF (nfac(i)/=4) GO TO 260 
      kspnn = kspan 
      kspan = kspan/4 
  180 c1 = 1.0 
      s1 = 0 
  190 k1 = kk + kspan 
      k2 = k1 + kspan 
      k3 = k2 + kspan 
      akp = a(kk) + a(k2) 
      akm = a(kk) - a(k2) 
      ajp = a(k1) + a(k3) 
      ajm = a(k1) - a(k3) 
      a(kk) = akp + ajp 
      ajp = akp - ajp 
      bkp = b(kk) + b(k2) 
      bkm = b(kk) - b(k2) 
      bjp = b(k1) + b(k3) 
      bjm = b(k1) - b(k3) 
      b(kk) = bkp + bjp 
      bjp = bkp - bjp 
      IF (isn<0) GO TO 220 
      akp = akm - bjm 
      akm = akm + bjm 
      bkp = bkm + ajm 
      bkm = bkm - ajm 
      IF (s1==0.0) GO TO 230 
  200 a(k1) = akp*c1 - bkp*s1 
      b(k1) = akp*s1 + bkp*c1 
      a(k2) = ajp*c2 - bjp*s2 
      b(k2) = ajp*s2 + bjp*c2 
      a(k3) = akm*c3 - bkm*s3 
      b(k3) = akm*s3 + bkm*c3 
      kk = k3 + kspan 
      IF (kk<=nt) GO TO 190 
  210 c2 = c1 - (cd*c1+sd*s1) 
      s1 = (sd*c1-cd*s1) + s1 
!     the following three statements compensate for truncation          
!     error. if rounded arithmetic is used, they may be deleted.        
!     c1=0.5/(c2**2+s1**2)+0.5                                          
!     s1=c1*s1                                                          
!     c1=c1*c2                                                          
!     next statement should be deleted if non-rounded arithmetic is used
      c1 = c2 
      c2 = c1**2 - s1**2 
      s2 = 2.0*c1*s1 
      c3 = c2*c1 - s2*s1 
      s3 = c2*s1 + s2*c1 
      kk = kk - nt + jc 
      IF (kk<=kspan) GO TO 190 
      kk = kk - kspan + inc 
      IF (kk<=jc) GO TO 180 
      IF (kspan==jc) GO TO 360 
      GO TO 120 
  220 akp = akm + bjm 
      akm = akm - bjm 
      bkp = bkm - ajm 
      bkm = bkm + ajm 
      IF (s1/=0.0) GO TO 200 
  230 a(k1) = akp 
      b(k1) = bkp 
      a(k2) = ajp 
      b(k2) = bjp 
      a(k3) = akm 
      b(k3) = bkm 
      kk = k3 + kspan 
      IF (kk<=nt) GO TO 190 
      GO TO 210 
!     transform for factor of 5 (optional code)                         
  240 c2 = c72**2 - s72**2 
      s2 = 2.0*c72*s72 
  250 k1 = kk + kspan 
      k2 = k1 + kspan 
      k3 = k2 + kspan 
      k4 = k3 + kspan 
      akp = a(k1) + a(k4) 
      akm = a(k1) - a(k4) 
      bkp = b(k1) + b(k4) 
      bkm = b(k1) - b(k4) 
      ajp = a(k2) + a(k3) 
      ajm = a(k2) - a(k3) 
      bjp = b(k2) + b(k3) 
      bjm = b(k2) - b(k3) 
      aa = a(kk) 
      bb = b(kk) 
      a(kk) = aa + akp + ajp 
      b(kk) = bb + bkp + bjp 
      ak = akp*c72 + ajp*c2 + aa 
      bk = bkp*c72 + bjp*c2 + bb 
      aj = akm*s72 + ajm*s2 
      bj = bkm*s72 + bjm*s2 
      a(k1) = ak - bj 
      a(k4) = ak + bj 
      b(k1) = bk + aj 
      b(k4) = bk - aj 
      ak = akp*c2 + ajp*c72 + aa 
      bk = bkp*c2 + bjp*c72 + bb 
      aj = akm*s2 - ajm*s72 
      bj = bkm*s2 - bjm*s72 
      a(k2) = ak - bj 
      a(k3) = ak + bj 
      b(k2) = bk + aj 
      b(k3) = bk - aj 
      kk = k4 + kspan 
      IF (kk<nn) GO TO 250 
      kk = kk - nn 
      IF (kk<=kspan) GO TO 250 
      GO TO 320 
!     transform for odd factors                                         
  260 k = nfac(i) 
      kspnn = kspan 
      kspan = kspan/k 
      IF (k==3) GO TO 160 
      IF (k==5) GO TO 240 
      IF (k==jf) GO TO 280 
      jf = k 
      s1 = rad/real(k) 
      c1 = cos(s1) 
      s1 = sin(s1) 
      IF (jf>maxf) GO TO 590 
      ck(jf) = 1.0 
      sk(jf) = 0.0 
      j = 1 
  270 ck(j) = ck(k)*c1 + sk(k)*s1 
      sk(j) = ck(k)*s1 - sk(k)*c1 
      k = k - 1 
      ck(k) = ck(j) 
      sk(k) = -sk(j) 
      j = j + 1 
      IF (j<k) GO TO 270 
  280 k1 = kk 
      k2 = kk + kspnn 
      aa = a(kk) 
      bb = b(kk) 
      ak = aa 
      bk = bb 
      j = 1 
      k1 = k1 + kspan 
  290 k2 = k2 - kspan 
      j = j + 1 
      at(j) = a(k1) + a(k2) 
      ak = at(j) + ak 
      bt(j) = b(k1) + b(k2) 
      bk = bt(j) + bk 
      j = j + 1 
      at(j) = a(k1) - a(k2) 
      bt(j) = b(k1) - b(k2) 
      k1 = k1 + kspan 
      IF (k1<k2) GO TO 290 
      a(kk) = ak 
      b(kk) = bk 
      k1 = kk 
      k2 = kk + kspnn 
      j = 1 
  300 k1 = k1 + kspan 
      k2 = k2 - kspan 
      jj = j 
      ak = aa 
      bk = bb 
      aj = 0.0 
      bj = 0.0 
      k = 1 
  310 k = k + 1 
      ak = at(k)*ck(jj) + ak 
      bk = bt(k)*ck(jj) + bk 
      k = k + 1 
      aj = at(k)*sk(jj) + aj 
      bj = bt(k)*sk(jj) + bj 
      jj = jj + j 
      IF (jj>jf) jj = jj - jf 
      IF (k<jf) GO TO 310 
      k = jf - j 
      a(k1) = ak - bj 
      b(k1) = bk + aj 
      a(k2) = ak + bj 
      b(k2) = bk - aj 
      j = j + 1 
      IF (j<k) GO TO 300 
      kk = kk + kspnn 
      IF (kk<=nn) GO TO 280 
      kk = kk - nn 
      IF (kk<=kspan) GO TO 280 
!     multiply by rotation factor (except for factors of 2 and 4)       
  320 IF (i==m) GO TO 360 
      kk = jc + 1 
  330 c2 = 1.0 - cd 
      s1 = sd 
  340 c1 = c2 
      s2 = s1 
      kk = kk + kspan 
  350 ak = a(kk) 
      a(kk) = c2*ak - s2*b(kk) 
      b(kk) = s2*ak + c2*b(kk) 
      kk = kk + kspnn 
      IF (kk<=nt) GO TO 350 
      ak = s1*s2 
      s2 = s1*c2 + c1*s2 
      c2 = c1*c2 - ak 
      kk = kk - nt + kspan 
      IF (kk<=kspnn) GO TO 350 
      c2 = c1 - (cd*c1+sd*s1) 
      s1 = s1 + (sd*c1-cd*s1) 
!     the following three statements compensate for truncation          
!     error. if rounded arithmetic is used, they may                    
!     be deleted.                                                       
!     c1=0.5/(c2**2+s1**2)+0.5                                          
!     s1=c1*s1                                                          
!     c2=c1*c2                                                          
      kk = kk - kspnn + jc 
      IF (kk<=kspan) GO TO 340 
      kk = kk - kspan + jc + inc 
      IF (kk<=jc+jc) GO TO 330 
      GO TO 120 
!     permute the results to normal order---done in two stages          
!     permutation for square factors of n                               
  360 np(1) = ks 
      IF (kt==0) GO TO 450 
      k = kt + kt + 1 
      IF (m<k) k = k - 1 
      j = 1 
      np(k+1) = jc 
  370 np(j+1) = np(j)/nfac(j) 
      np(k) = np(k+1)*nfac(j) 
      j = j + 1 
      k = k - 1 
      IF (j<k) GO TO 370 
      k3 = np(k+1) 
      kspan = np(2) 
      kk = jc + 1 
      k2 = kspan + 1 
      j = 1 
      IF (n/=ntot) GO TO 410 
!     permutation for single-variate transform (optional code)          
  380 ak = a(kk) 
      a(kk) = a(k2) 
      a(k2) = ak 
      bk = b(kk) 
      b(kk) = b(k2) 
      b(k2) = bk 
      kk = kk + inc 
      k2 = kspan + k2 
      IF (k2<ks) GO TO 380 
  390 k2 = k2 - np(j) 
      j = j + 1 
      k2 = np(j+1) + k2 
      IF (k2>np(j)) GO TO 390 
      j = 1 
  400 IF (kk<k2) GO TO 380 
      kk = kk + inc 
      k2 = kspan + k2 
      IF (k2<ks) GO TO 400 
      IF (kk<ks) GO TO 390 
      jc = k3 
      GO TO 450 
!     permutation for multivariate transform                            
  410 k = kk + jc 
  420 ak = a(kk) 
      a(kk) = a(k2) 
      a(k2) = ak 
      bk = b(kk) 
      b(kk) = b(k2) 
      b(k2) = bk 
      kk = kk + inc 
      k2 = k2 + inc 
      IF (kk<k) GO TO 420 
      kk = kk + ks - jc 
      k2 = k2 + ks - jc 
      IF (kk<nt) GO TO 410 
      k2 = k2 - nt + kspan 
      kk = kk - nt + jc 
      IF (k2<ks) GO TO 410 
  430 k2 = k2 - np(j) 
      j = j + 1 
      k2 = np(j+1) + k2 
      IF (k2>np(j)) GO TO 430 
      j = 1 
  440 IF (kk<k2) GO TO 410 
      kk = kk + jc 
      k2 = kspan + k2 
      IF (k2<ks) GO TO 440 
      IF (kk<ks) GO TO 430 
      jc = k3 
  450 IF (2*kt+1>=m) RETURN 
      kspnn = np(kt+1) 
!     permutation for square-free factors of n                          
      j = m - kt 
      nfac(j+1) = 1 
  460 nfac(j) = nfac(j)*nfac(j+1) 
      j = j - 1 
      IF (j/=kt) GO TO 460 
      kt = kt + 1 
      nn = nfac(kt) - 1 
      IF (nn>maxp) GO TO 590 
      jj = 0 
      j = 0 
      GO TO 490 
  470 jj = jj - k2 
      k2 = kk 
      k = k + 1 
      kk = nfac(k) 
  480 jj = kk + jj 
      IF (jj>=k2) GO TO 470 
      np(j) = jj 
  490 k2 = nfac(kt) 
      k = kt + 1 
      kk = nfac(k) 
      j = j + 1 
      IF (j<=nn) GO TO 480 
!     determine the permutation cycles of length greater than 1         
      j = 0 
      GO TO 510 
  500 k = kk 
      kk = np(k) 
      np(k) = -kk 
      IF (kk/=j) GO TO 500 
      k3 = kk 
  510 j = j + 1 
      kk = np(j) 
      IF (kk<0) GO TO 510 
      IF (kk/=j) GO TO 500 
      np(j) = -j 
      IF (j/=nn) GO TO 510 
      maxf = inc*maxf 
!     reorder a and b, following the permutation cycles                 
      GO TO 580 
  520 j = j - 1 
      IF (np(j)<0) GO TO 520 
      jj = jc 
  530 kspan = jj 
      IF (jj>maxf) kspan = maxf 
      jj = jj - kspan 
      k = np(j) 
      kk = jc*k + ii + jj 
      k1 = kk + kspan 
      k2 = 0 
  540 k2 = k2 + 1 
      at(k2) = a(k1) 
      bt(k2) = b(k1) 
      k1 = k1 - inc 
      IF (k1/=kk) GO TO 540 
  550 k1 = kk + kspan 
      k2 = k1 - jc* (k+np(k)) 
      k = -np(k) 
  560 a(k1) = a(k2) 
      b(k1) = b(k2) 
      k1 = k1 - inc 
      k2 = k2 - inc 
      IF (k1/=kk) GO TO 560 
      kk = k2 
      IF (k/=j) GO TO 550 
      k1 = kk + kspan 
      k2 = 0 
  570 k2 = k2 + 1 
      a(k1) = at(k2) 
      b(k1) = bt(k2) 
      k1 = k1 - inc 
      IF (k1/=kk) GO TO 570 
      IF (jj/=0) GO TO 530 
      IF (j/=1) GO TO 520 
  580 j = k3 + 1 
      nt = nt - kspnn 
      ii = nt - inc + 1 
      IF (nt>=0) GO TO 520 
      RETURN 
!     error finish, insufficient array storage                          
  590 isn = 0 
      WRITE (6,FMT=8000) 
      CALL juDFT_error("array bounds exceeded within subroutine cft",calledby="gf_fft.F90")
 8000 FORMAT ('array bounds exceeded within subroutine cft') 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
