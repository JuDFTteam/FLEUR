!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_energies
      use m_juDFT
      USE m_constants, ONLY: oUnit
!***************************************************                    
! This module supplies all energies for which the GF                    
! should be calculated                                                  
!                  Daniel Wortman, June 2001                            
!***************************************************                    
!DOKU                                                                   
!The energies are read in from gf_en                                    
!Three different ways of specifying energies are possible in this       
!file                                                                   
!1. Listing all values                                                  
!Format-Beispiel:                                                       
!list    3                                                              
!1.0000000 0.00000000                                                   
!1.1000000 0.00000000                                                   
!1.2000000 0.00000000                                                   
!                                                                       
!2. Specifing an intervall                                              
!inter   3                                                              
!1.0000000 0.00000000                                                   
!1.2000000 0.00000000                                                   
!                                                                       
!3. Specifing an intervall for a complex semi-circle-loop               
!circle  3                                                              
!1.0000000 0.00000000                                                   
!1.2000000 0.00000000!                                                  
!                                                                       
! Formats can be mixed in one gf_en file!                               
!                                                                       
! Additionally in column 30+ the following can be specified:            
! 0.000000 0.00000000  +L-R                                             
! indicating that in self-consistency mode the charge due to incomming  
! bloch waves from the left should be added and the charge dur to       
! states incomming from the right should be subtracted for this         
! energy. This is possible for each energy in a list, or can be         
! specified for the first energy for the whole interval. All            
! combinations like +R-L,+R+L are possible                              
!     in the first line a shift=... can be specified to shift all       
!     following energies by some offset                                 
                                                                        
                                                                        
      IMPLICIT NONE 
      PRIVATE 
      COMPLEX,ALLOCATABLE,SAVE::z(:) 
      COMPLEX,ALLOCATABLE,SAVE::wz(:) 
      INTEGER,ALLOCATABLE,SAVE::direction(:) 
      INTEGER,PARAMETER::FROM_LEFT=1 
      INTEGER,PARAMETER::FROM_RIGHT=10 
      INTEGER,SAVE     :: energies,bias_layer 
      REAL,SAVE        :: bias 
      PUBLIC      ::init,gf_NoEn,gf_Z,gf_weightZ,gf_allz,FROM_LEFT      &
     &     ,FROM_RIGHT,direction,gf_bias_layer                          
      CONTAINS 
      !<-- F:gf_bias_layer(bias)                                        
      FUNCTION gf_bias_layer(en_bias) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:09-07-22) D. Wortmann                        
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(OUT),OPTIONAL    :: en_bias 
      INTEGER                         :: gf_bias_layer 
      !>                                                                
                                                                        
      IF (PRESENT(en_bias)) en_bias = bias 
                                                                        
      gf_bias_layer = bias_layer 
                                                                        
      END FUNCTION 
      !>                                                                
      !<--                                                              
      FUNCTION gf_NoEn() 
!******************************************                             
!     return the total No of energies                                   
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      INTEGER :: gf_NoEn 
      gf_NoEn = energies 
      END FUNCTION 
      !>                                                                
      !<--                                                              

      FUNCTION gf_weightZ(en) 
!******************************************                             
!     returns the weight no en                                          
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en 
      COMPLEX :: gf_weightZ 
      !>                                                                
      IF (en>energies) CALL juDFT_error('ERROR in m_gf_energies%gf_Z') 
      gf_weightZ= wz(en) 
      END FUNCTION 

      !>                                                                
      !<--                                                              

      FUNCTION gf_Z(en,layer) 
!******************************************                             
!     returns the energy no en                                          
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: en 
      INTEGER,INTENT(IN) :: layer 
      COMPLEX            :: gf_Z 
      !>                                                                
      IF (en>energies) CALL juDFT_error('ERROR in m_gf_energies%gf_Z') 
      IF (layer<bias_layer) THEN 
         gf_Z = Z(en) 
      ELSE 
         gf_Z = Z(en)+bias 
      ENDIF 
      END FUNCTION 

      !>                                                                
      !<--                                                              
                                                                        
      FUNCTION gf_ALLZ(layer) 
!******************************************                             
!     returns all the energies                                          
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: layer 
      COMPLEX                :: gf_ALLZ(energies) 
      !>                                                                
      IF (layer<bias_layer) THEN 
         gf_ALLZ = z 
      ELSE 
         gf_allz = z+bias 
      ENDIF 
      END FUNCTION 
                                                                        
      !>                                                                
  ! this subroutine reads the gf_en file and counts the total           
  ! no of energies used. It sets energies                               
      SUBROUTINE init(pe0) 
      LOGICAL,INTENT(IN)::pe0 
      INTEGER:: i,en,no 
      REAL   :: shift 
      CHARACTER::type 
      en=0 
      OPEN(37,FILE='gf_en',FORM='formatted',STATUS='old',ERR=200) 
                                          !check for shift              
      READ(37,'(a)') TYPE 
      REWIND 37 
      IF (TYPE=='S'.OR.TYPE=='s') THEN 
         READ(37,'(7x,f20.0)') shift 
      ELSE 
         shift=0.0 
      ENDIF 
      !set bias layer is present                                        
                          !check for shift                              
      READ(37,'(a)') TYPE 
      REWIND 37 
      IF (TYPE=='b'.OR.TYPE=='b') THEN 
         READ(37,'(7x,i20)') bias_layer 
         READ(37,'(7x,f20.0)') bias 
      ELSE 
         bias_layer = 9999 
         bias = 0.0 
      ENDIF 
      DO 
         READ(37,'(a,7x,i10)',END=200,ERR=200) type,no 
                                !skip the next lines                    
         IF (TYPE=='l'.OR.TYPE=='L') THEN 
            DO i=1,no 
               READ(37,*,END=200,ERR=200) 
            ENDDO 
         ELSEIF(type=='k'.OR.type=='K') THEN 
            READ(37,*,END=200,ERR=200) 
            READ(37,*,END=200,ERR=200) 
            READ(37,*,END=200,ERR=200) 
            READ(37,*,END=200,ERR=200) 
            READ(37,*,END=200,ERR=200) 
         ELSE 
            READ(37,*,END=200,ERR=200) 
            READ(37,*,END=200,ERR=200) 
         ENDIF 
         en=en+no 
      ENDDO 
      !here the file has been read completely                           
  200 CLOSE(37) 
      IF (en==0) CALL juDFT_error('GF_energies: No energies specified!') 
      energies=en 
      ! now read the energies!                                          
      CALL getenergies() 
                  !apply shift                                          
      z=z+shift 
      IF (pe0) THEN 
         WRITE(oUnit,*) 'ENERGY MESH' 
         WRITE(oUnit,*) 'No of energies:',energies 
         WRITE(oUnit,*) 'Energy-point             Weight' 
         DO i=1,energies 
            WRITE(oUnit,'(2(2(e15.8,1x),2x))') z(i),wz(i) 
         ENDDO 
      ENDIF 
      IF (ANY(direction>0.AND.abs(aimag(z))>epsilon(0.))) THEN 
         CALL juDFT_error('gf_energies:direction>0&Im(z)>0') 
      ENDIF 
      END SUBROUTINE init 
                                                                        
      SUBROUTINE getenergies() 
      USE m_gf_kkrpoints 
      IMPLICIT NONE 
      INTEGER::i,no,en 
      COMPLEX :: en_slide,en_start,en_STOP 
      CHARACTER::type 
      CHARACTER(LEN=10)::a 
      REAL             :: ef,temp 
      INTRINSIC INDEX 
      REAL simpsonw 
      IF (allocated(z)) THEN 
          ! the energies are read only once!!                           
         RETURN 
      ENDIF 
      ALLOCATE(z(energies),wz(energies),direction(energies)) 
      direction=0 
      en=0 
      OPEN(37,FILE='gf_en',FORM='formatted',STATUS='old',ERR=200) 
      !Check if there is a shift in the first line                      
      READ(37,'(a,5x)',END=200,ERR=200) TYPE 
      REWIND 37 
      IF (TYPE=='s'.OR.TYPE=='S') THEN 
         READ(37,*) 
      ENDIF 
      IF (TYPE=='b'.OR.TYPE=='B') THEN 
         READ(37,*) 
         READ(37,*) 
      ENDIF 
      DO 
         READ(37,'(a,7x,i10)',END = 200,ERR=200) TYPE,no 
                                !read the energies                      
         IF (TYPE=='l'.OR.TYPE=='L') THEN 
            DO i=1,no 
               en=en+1 
               a='          ' 
               READ(37,'(4f10.5,a10)',END = 200,ERR=200) z(en),wz(en),a 
               IF (wz(en) == 0.0)  wz(en)=1.0 
               IF (index(a,"-R")>0.OR.index(a,"-r")>0) THEN 
                  direction(en)=direction(en)-FROM_RIGHT 
               ENDIF 
               IF (index(a,'+r')>0.OR.index(a,'+R')>0) THEN 
                  direction(en)=direction(en)+FROM_RIGHT 
               ENDIF 
               IF (index(a,'-L')>0.OR.index(a,'-l')>0) THEN 
                  direction(en)=direction(en)-FROM_LEFT 
               ENDIF 
               IF (index(a,'+L')>0.OR.index(a,'+l')>0) THEN 
                  direction(en)=direction(en)+FROM_LEFT 
               ENDIF 
            ENDDO 
         ELSE IF(type=='i'.OR.type=='I') THEN 
                                ! construct the energies in the interval
            a='          ' 
            READ(37,'(2f10.8,a10)',END=200,ERR=200) en_start,a 
            READ(37,'(2f10.8)',END=200,ERR=200) en_stop

            IF (no<2)                                                &
     &           CALL juDFT_error                                         &
     &           ('Energy-Interval-Mode with less than two energies' )  
            en_slide=(en_stop-en_start)/(no-1.0) 
            DO i=1,no 
               en=en+1 
               z(en)=en_start+(i-1)*en_slide 
               wz(en)=(en_stop-en_start)/no 
               IF (index(a,'-R')>0.OR.index(a,'-r')>0) THEN 
                  direction(en)=direction(en)-FROM_RIGHT 
               ENDIF 
               IF (index(a,'+r')>0.OR.index(a,'+R')>0) THEN 
                  direction(en)=direction(en)+FROM_RIGHT 
               ENDIF 
               IF (index(a,'-L')>0.OR.index(a,'-l')>0) THEN 
                  direction(en)=direction(en)-FROM_LEFT 
               ENDIF 
               IF (index(a,'+L')>0.OR.index(a,'+l')>0) THEN 
                  direction(en)=direction(en)+FROM_LEFT 
               ENDIF 
            ENDDO 
                                             !use simpson               
         ELSE IF(type=='n'.OR.type=='N')THEN 
                                   ! construct the energies in the inter
            IF(MOD(no,2) == 0) WRITE(*,*)                               &
     &           "Number of energies should be odd"                     
            a='          ' 
            READ(37,'(2f10.8,a10)',END=200,ERR=200) en_start,a 
            READ(37,'(2f10.8)',END=200,ERR=200) en_stop
 
            IF (no<2)                                                &
     &           CALL juDFT_error                                         &
     &           ('Energy-Interval-Mode with less than two energies' )  
            en_slide=(en_stop-en_start)/(no-1.0) 
            simpsonw=en_slide/3.0 
            DO i=1,no 
               en=en+1 
               z(en)=en_start+(i-1)*en_slide 
               IF(i==1)THEN 
                  wz(en)=simpsonw 
               ELSEIF(i==no)THEN 
                  wz(en)=simpsonw 
               ELSEIF(mod(i,2)==0)THEN 
                  wz(en)=4.0*simpsonw 
               ELSE 
                  wz(en)=2.0*simpsonw 
               ENDIF 
               IF (index(a,'-R')>0.OR.index(a,'-r')>0) THEN 
                  direction(en)=direction(en)-FROM_RIGHT 
               ENDIF 
               IF (index(a,'+r')>0.OR.index(a,'+R')>0) THEN 
                  direction(en)=direction(en)+FROM_RIGHT 
               ENDIF 
               IF (index(a,'-L')>0.OR.index(a,'-l')>0) THEN 
                  direction(en)=direction(en)-FROM_LEFT 
               ENDIF 
               IF (index(a,'+L')>0.OR.index(a,'+l')>0) THEN 
                  direction(en)=direction(en)+FROM_LEFT 
               ENDIF 
            ENDDO 
         ELSEIF(type=='k'.OR.type=='K') THEN 
            CALL gf_kkrpoints(37,z(en+1:en+no),wz(en+1:en+no)) 
            en=en+no 
         ELSEIF(TYPE =='o'.OR.TYPE =='O') THEN 
            READ(37,*) ef 
            READ(37,*) temp 
            IF (no<2) CALL juDFT_error("Ozaki's method must use n>1") 
            CALL priv_ozaki(temp,ef,(no)*2,z(en+1:en+no),wz(en+1:en     &
     &           +no))                                                  
!            z(en+no) = CMPLX(ef,1e10)                                  
!            wz(en+no) = 1e10                                           
            en=en+no 
         ELSE 
            ! construct the energies on a semicircle in the complex plan
            READ(37,'(2f10.8)',END=200,ERR=200) en_start 
            READ(37,'(2f10.8)',END = 200,ERR=200) en_STOP 
 
            ! check if no=2**n-1                                        
            i=2 
            test:DO 
               IF (i-1==no) EXIT test 
               IF (i>no) THEN 
                  WRITE(*,*) 'WARNING, no of energy points in',         &
     &           ' semi-circle must be modified to be of form 2**n-1'   
                  WRITE(*,*) 'Old:',no,'New:',i-1 
                  WRITE(oUnit,*) 'WARNING, no of energy points in',         &
     &           'semi-circle must be modified to be of form 2**n-1'    
                  WRITE(oUnit,*) 'Old:',no,'New:',i-1 
                  no=i-1 
                  CALL juDFT_error('gf_energies') 
               ENDIF 
               i=2*i 
            ENDDO test 
            IF (abs(aimag(en_start))>epsilon(1.0)                    &
     &           .OR.abs(aimag(en_stop))>epsilon(1.0)) THEN          
               WRITE(*,*) 'Upper and lower energy for semi-circle',     &
     &              ' must be real energies!'                           
               CALL juDFT_error('gf_energies') 
            ENDIF 
            CALL spherical_points(real(en_start),real(en_stop),         &
     &           no,z(en+1:en+no),                                      &
     &           wz(en+1:en+no))                                        
            en=en+no 
         ENDIF 
      ENDDO 
                                !here the file has been read completely 
  200 CLOSE(37) 
                                                                        
      END SUBROUTINE getenergies 
                                                                        
!*********************************************************************  
!     Two subroutines BASED on subroutines of S. Crampin                
!*********************************************************************  
                                                                        
      SUBROUTINE spherical_points(emin,ef,npts,ze,we) 
!                                                                       
! determine coordinates and weights for integration around              
! semicircular contour in upper half plane starting at (emin,0)         
! ending at (ef,0) using modified gauss-chebyshev with npts             
! points                                                                
!                                                                       
! on entry - emin,ef,npts specify range and required number             
!            of points (npts=2**n-1 for positive integer n)             
!          - x,w and real*8 workspace arrays of dimension               
!            >= npts. on exit contain coordinates and weights           
!            for integral along unit interval.                          
! on exit  - complex*8 arrays ze,we contain coordinates and             
!            weights for contour integral                               
! s crampin                                                             
!                                                                       
      IMPLICIT NONE 
      INTEGER npts,k,l 
      COMPLEX ze(npts),we(npts) 
      REAL emin,ef,a,r,pi,xx,ww 
      REAL,ALLOCATABLE :: x(:),w(:) 
      ALLOCATE (x(npts),w(npts)) 
      CALL spherical_weights(npts,x,w) 
      IF(1<0) THEN 
      DO k = 1,npts 
      DO l = k+1,npts 
      IF(x(l) < x(k)) THEN 
        xx = x(k) 
        ww = w(k) 
        x(k) = x(l) 
        w(k) = w(l) 
        x(l) = xx 
        w(l) = ww 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDIF 
      a=0.5d0*(ef+emin) 
      r=0.5d0*(ef-emin) 
      pi=4.0d0*atan(1.0d0) 
      DO k=1,npts 
       ze(k)=a+r*exp(cmplx(0.0d0,pi*(1-x(k)))) 
       we(k)=w(k)*cmplx(0.0d0,pi)*(a-ze(k)) 
      ENDDO 
      DEALLOCATE (x,w) 
      RETURN 
      END SUBROUTINE 
!---*----1----*----2----*----3----*----4----*----5----*----6----*----7--
      SUBROUTINE spherical_weights(npts,x,w) 
!                                                                       
! determine coordinates and weights for integration between 0 and 1 by  
! modified Gauss-Chebyshev technique.                                   
!                                                                       
! On input - npts (=2**n-1 for some positive integer n) is number       
!                  sampling points required. Unchanged upon exit        
!            x(),w() arrays of dimension >= npts                        
!                                                                       
! On exit - x and w contains coordinates and weights                    
!                                                                       
! s crampin                                                             
!                                                                       
      IMPLICIT NONE 
      INTEGER npts,k,mm,m,n,i 
      REAL twopi,a,b,c,theta 
      REAL x(npts),w(npts) 
      twopi=8.0d0*atan(1.0d0) 
      b=1.0d0/(npts+1.0d0) 
      x(1)=0.5d0 
      w(1)=b+b 
      IF(npts>1)THEN 
       k=2 
       a=0.25D0 
       mm=nint(dlog(npts+1.0d0)/dlog(2.0d0)) 
       DO m=2,mm 
        n=2**m-1 
        DO i=1,n,2 
         c=a*i 
         theta=twopi*c 
         x(k)=c-sin(theta)/twopi 
         w(k)=b*(1.0d0-cos(theta)) 
         k=k+1 
        ENDDO 
        a=0.5*a 
       ENDDO 
      ENDIF 
      RETURN 
      END SUBROUTINE 
                                                                        
                                                                        
      !<-- S:priv_ozaki(temp,ef,n)                                      
      SUBROUTINE priv_ozaki(temp,ef,num,z,wz) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:09-06-26) D. Wortmann                        
!-----------------------------------------------                        
      USE m_constants 
      IMPLICIT NONE 
      !<--Arguments                                                     
      REAL   ,INTENT(IN)     :: temp,ef 
      INTEGER,INTENT(IN)     :: num 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: a(num,num),b(num,num) 
      REAL                :: ew(num),r(num),work(10*num+16) 
      INTEGER             :: n,info 
      REAL,PARAMETER      :: k_b=0.316679E-5 
      COMPLEX,INTENT(OUT) :: z(:),wz(:) 
      !>                                                                
                                                                        
      !<-- construct the A and B matrices                               
      A = 0 
      B = 0 
      DO n = 1,num 
         a(n,n) = (2*n-1) 
         IF (n>1)   b(n-1,n) =-0.5 
         IF (n<num) b(n+1,n) =-0.5 
      ENDDO 
      !>                                                                
                                                                        
      !<-- solve generalized Eigenvalue problem                         
      CALL dsygv(1,'V','U',num,b,num,a,num,ew,work,10*n+16,info) 
      !>                                                                
      !<-- eigenvalues of inverse                                       
      ew = 1./ew 
      !>                                                                
      !<-- calculate residues                                           
      DO n = 1,num 
         r(n) = -b(1,n)*b(1,n)*ew(n)*ew(n)/4.0 
      ENDDO 
      !>                                                                
      info = 1 
      DO n = 1,num 
         IF (ew(n)>0) THEN 
            IF (info>SIZE(z)) CALL juDFT_error("Error in priv_ozaki") 
            z(info) = CMPLX(ef,ew(n)*k_b*temp) 
            wz(info) = -2.*r(n)*CMPLX(0.0,k_b*temp*pi_const) 
            info = info+1 
         ENDIF 
      ENDDO 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
