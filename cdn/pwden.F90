!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_pwden
CONTAINS
  SUBROUTINE pwden(stars,kpts,banddos,oneD, input,mpi,noco,cell,atoms,sym, &
       ikpt,jspin,lapw,ne, igq_fft,we,z, eig,bkpt, qpw,cdom, qis,forces,f_b8)
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !     In this subroutine the star function expansion coefficients of
    !     the plane wave charge density is determined.
    !
    !     This subroutine is called for each k-point and each spin.
    !
    !
    !     Two methods are implemented to calculate the charge density
    !     1) which uses the FFT. The effort in calculating the charge
    !        density is proportional to M * N * log(N) , M being number of
    !        states and N being number of plane waves. This is the method
    !        which we use for production runs
    !     2) the traditional method for calculating the charge density
    !        using the double summation. In this case the effort scales as
    !        M * N * N. The method is only used for test purposes or for
    !        special cases.
    !
    !
    !     INPUT:    eigen vectors
    !               reciprocal lattice information
    !               Brillouine zone sampling
    !               FFT information
    !
    !     OUTPUT:   qpw(s)
    !               1) using FFT
    !
    !                2) traditional method
    !
    !                             -1             ef
    !                qpw  (g) = vol * sum{ sum{ sum{ sum{ w(k) * f(nu) *
    !                                  sp   k    nu   g'
    !                                     *
    !                                    c(g'-g,nu,k) * c(g',nu,k) } } } }
    !                or :
    !                             -1             ef
    !                qpw  (g) = vol * sum{ sum{ sum{ sum{ w(k) * f(nu) *
    !                                  sp   k    nu   g'
    !                                     *
    !                                    c(g',nu,k) * c(g'+g,nu,k) } } } }
    !
    !                qpw(g) are actuall 
    ! 
    !                the weights w(k) are normalized: sum{w(k)} = 1
    !                                                  k                -6
    !                         a) 1                           for kT < 10
    !                f(nu) = {                           -1             -6
    !                         b){ 1+exp(e(k,nu) -ef)/kt) }   for kt >=10
    !
    !
    !                                      Stefan Bl"ugel, JRCAT, Feb. 1997
    !                                      Gustav Bihlmayer, UniWien       
    !
    !     In non-collinear calculations the density becomes a hermitian 2x2
    !     matrix. This subroutine generates this density matrix in the 
    !     interstitial region. The diagonal elements of this matrix 
    !     (n_11 & n_22) are stored in qpw, while the real and imaginary part
    !     of the off-diagonal element are store in cdom. 
    !
    !     Philipp Kurz 99/07
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !
    !
#include"cpp_double.h"
    USE m_forceb8
    USE m_pwint
    USE m_juDFT
    USE m_rfft
    USE m_cfft
    USE m_types
    IMPLICIT NONE
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_mpi),INTENT(IN)      :: mpi
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_banddos),INTENT(IN)  :: banddos
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_stars),INTENT(IN)    :: stars
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_kpts),INTENT(IN)     :: kpts
    TYPE(t_atoms),INTENT(IN)    :: atoms
    INTEGER, INTENT (IN)        :: igq_fft(0:stars%kq1d*stars%kq2d*stars%kq3d-1)
    REAL,INTENT(IN)   :: we(:) !(nobd) 
    REAL,INTENT(IN)   :: eig(:)!(dimension%neigd)
    REAL,INTENT(IN)   :: bkpt(3)
    !----->  BASIS FUNCTION INFORMATION
    INTEGER,INTENT(IN):: ne
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
    COMPLEX, INTENT (IN) :: z(:,:) !(dimension%nbasfcn,dimension%neigd)
#else
    REAL,    INTENT (IN) :: z(:,:) !(dimension%nbasfcn,dimension%neigd)
#endif
    !----->  CHARGE DENSITY INFORMATION
    INTEGER,INTENT(IN)    :: ikpt,jspin 
    COMPLEX,INTENT(INOUT) :: qpw(:,:) !(stars%n3d,dimension%jspd)
    COMPLEX,INTENT(INOUT) :: cdom(:)!(stars%n3d)
    REAL,INTENT(OUT)      :: qis(:,:,:) !(dimension%neigd,kpts%nkptd,dimension%jspd)
    COMPLEX, INTENT (INOUT) ::  f_b8(3,atoms%ntypd)
    REAL,    INTENT (INOUT) :: forces(:,:,:) !(3,atoms%ntypd,dimension%jspd)
    !
    !-----> LOCAL VARIABLES
    !
    !----->  FFT  INFORMATION
    INTEGER :: ifftq2d,ifftq3d

    INTEGER  isn,nu,iv,ir,ik,il,im,in,istr,nw1,nw2,nw3,i,j
    INTEGER  ifftq1,ifftq2,ifftq3
    INTEGER  idens,ndens,ispin,jkpt,jsp_start,jsp_end
    REAL     q0,q0_11,q0_22,scale,xk(3)
    REAL     s
    COMPLEX  x
    INTEGER,PARAMETER::  ist(-1:1)=(/1,0,0/)
    REAL,PARAMETER:: zero   = 0.00,  tol_3=1.0e-3 
    !
    INTEGER  iv1d(SIZE(lapw%k1,1),input%jspins)
    REAL wtf(ne),wsave(stars%kq3d+15)
    REAL,    ALLOCATABLE :: psir(:),psii(:),rhon(:)
    REAL,    ALLOCATABLE :: psi1r(:),psi1i(:),psi2r(:),psi2i(:)
    REAL,    ALLOCATABLE :: rhomat(:,:)
    REAL,    ALLOCATABLE :: kpsir(:),kpsii(:)
    REAL,    ALLOCATABLE :: ekin(:)
    COMPLEX, ALLOCATABLE :: cwk(:),ecwk(:)
    !
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     REAL     CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot
#else
    COMPLEX  CPP_BLAS_cdotc
    EXTERNAL CPP_BLAS_cdotc
#endif

    !------->          ABBREVIATIONS
    !
    !     rhon  : charge density in real space
    !     ne    : number of occupied states
    !     nv    : number of g-components in eigenstate
    !     cv=z  : wavefunction in g-space (reciprocal space)
    !     psir   : wavefunction in r-space (real-space)
    !     cwk   : complex work array: charge density in g-space (as stars)
    !     qpw   : charge density stored as stars
    !     trdchg: logical key, determines the mode of charge density
    !             calculation: false (default) : fft
    !                          true            : double sum over stars
    !     we    : weights for the BZ-integration for a particular k-point
    !     omtil : volume (slab) unit cell, between -.5*D_tilde and +.5*D_tilde
    !     k1   : reciprocal lattice vectors G=G(k1,k2,k3) for wavefunction
    !     k2   :                             =k1*a_1 + k2*a_2 + k3*a_3
    !     k3   : where a_i= Bravais lattice vectors in reciprocal space
    !             kwi, depend on k-point.                            
    !     kq1d  : dimension of the charge density FFT box in the pos. domain
    !     kq2d  : defined in dimens.f program (subroutine apws).1,2,3 indicate
    !     kq3d  ; a_1, a_2, a_3 directions.
    !     kq(i) : i=1,2,3 actual length of the fft-box for which FFT is done.
    !     nstr  : number of members (arms) of reciprocal lattice (g) vector
    !             of each star.
    !     ng3_fft: number of stars in the  charge density  FFT-box
    !     ng3   : number of 3 dim. stars in the charge density sphere defined
    !             by gmax
    !     kmxq_fft: number of g-vectors forming the ng3_fft stars in the
    !               charge density sphere 
    !     kimax : number of g-vectors forming the ng3 stars in the gmax-sphere
    !     iv1d  : maps vector (k1,k2,k3) of wave function into one
    !             dimensional vector of cdn-fft box in positive domain.
    !     ifftq3d: elements (g-vectors) in the charge density  FFT-box
    !     igfft : pointer from the g-sphere (stored as stars) to fft-grid 
    !             and     from fft-grid to g-sphere (stored as stars)
    !     pgfft : contains the phases of the g-vectors of sph.     
    !     isn   : isn = +1, FFT transform for g-space to r-space
    !             isn = -1, vice versa
    !


    ALLOCATE(cwk(stars%n3d),ecwk(stars%n3d))

    IF (noco%l_noco) THEN
       ALLOCATE ( psi1r(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            psi1i(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            psi2r(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            psi2i(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            rhomat(0:stars%kq1d*stars%kq2d*stars%kq3d-1,4) )
    ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     ALLOCATE ( psir(-stars%kq1d*stars%kq2d:2*stars%kq1d*stars%kq2d*(stars%kq3d+1)-1),&
     psii(1),&
     rhon(-stars%kq1d*stars%kq2d:stars%kq1d*stars%kq2d*(stars%kq3d+1)-1) )
       IF (input%l_f) ALLOCATE ( kpsii(1),&
            kpsir(-stars%kq1d*stars%kq2d:2*stars%kq1d*stars%kq2d*(stars%kq3d+1)-1),&
            ekin(-stars%kq1d*stars%kq2d:2*stars%kq1d*stars%kq2d*(stars%kq3d+1)-1))
#else
       ALLOCATE ( psir(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            psii(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            rhon(0:stars%kq1d*stars%kq2d*stars%kq3d-1) )
       IF (input%l_f) ALLOCATE ( kpsir(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            kpsii(0:stars%kq1d*stars%kq2d*stars%kq3d-1),&
            ekin(0:stars%kq1d*stars%kq2d*stars%kq3d-1) )
#endif
    ENDIF
    !
    !=======>  CALCULATE CHARGE DENSITY USING FFT
    ! 
    !
    !------> setup FFT
    !
    ifftq1  = stars%kq1_fft
    ifftq2  = stars%kq1_fft*stars%kq2_fft
    ifftq3  = stars%kq1_fft*stars%kq2_fft*stars%kq3_fft
    ifftq3d = stars%kq1d*stars%kq2d*stars%kq3d
    ifftq2d = stars%kq1d*stars%kq2d
    !
    nw1=NINT(stars%kq1_fft/4.+0.3)
    nw2=NINT(stars%kq2_fft/4.+0.3)
    nw3=NINT(stars%kq3_fft/4.+0.3)
    !
    !------> g=0 star: calculate the charge for this k-point and spin
    !                  analytically to test the quality of FFT
    !
    q0 = zero
    q0_11 = zero
    q0_22 = zero
    IF (noco%l_noco) THEN
       q0_11 = zero
       q0_22 = zero
       DO nu = 1 , ne
#if ( !defined(CPP_INVERSION) )
          q0_11 = q0_11 + we(nu) *&
               CPP_BLAS_cdotc(lapw%nv(1),z(1,nu),1,z(1,nu),1)
          q0_22 = q0_22 + we(nu) *&
               CPP_BLAS_cdotc(lapw%nv(2),z(lapw%nv(1)+atoms%nlotot+1,nu),1,&
               z(lapw%nv(1)+atoms%nlotot+1,nu),1)
#endif
       ENDDO
       q0_11 = q0_11/cell%omtil
       q0_22 = q0_22/cell%omtil
    ELSE
       DO nu = 1 , ne
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     q0=q0+we(nu)*CPP_BLAS_sdot(lapw%nv(jspin),z(1,nu),1,z(1,nu),1)
#else
          q0=q0+we(nu)&
               *REAL(CPP_BLAS_cdotc(lapw%nv(jspin),z(1,nu),1,z(1,nu),1))
#endif
       ENDDO
       q0 = q0/cell%omtil
    ENDIF
    !
    !--------> initialize charge density with zero
    !
    IF (noco%l_noco) THEN
       rhomat = 0.0
       IF (ikpt.LE.mpi%isize) THEN
          qis=0.0
       ENDIF
    ELSE
       rhon=0.0
       IF (input%l_f) ekin=0.0
    ENDIF
    !
    !------> calculate:  wtf(nu,k) =  w(k)*f(nu,k)/vol
    !
    wtf(:ne) = we(:ne)/cell%omtil
    !
    !------> prepare mapping from wave function box to cdn FFT box
    !
    IF (noco%l_ss) THEN
       jsp_start = 1
       jsp_end   = 2
    ELSE
       jsp_start = jspin
       jsp_end   = jspin
    ENDIF
    DO ispin = jsp_start,jsp_end
       DO iv = 1 , lapw%nv(ispin)
          !                                              -k1d <= L <= k1d
          !                                              -k2d <= M <= k2d
          !                                              -k3d <= N <= k3d
          il = lapw%k1(iv,ispin)
          im = lapw%k2(iv,ispin)
          in = lapw%k3(iv,ispin)
          !
          !------>  L,M,N LATTICE POINTS OF G-VECTOR IN POSITIVE DOMAIN
          !         (since charge density box = two times charge density box
          !          wrap arround error should not occur )
          !                                           0<= L <=2*k1-1 = kq1_fft-1
          !                                           0<= M <=2*k2-1 = kq2_fft-1
          !                                           0<= N <=2*k3-1 = kq3_fft-1
          !
          il = il  +  stars%kq1_fft * ist( isign(1,il) )
          im = im  +  stars%kq2_fft * ist( isign(1,im) )
          in = in  +  stars%kq3_fft * ist( isign(1,in) )
          !
          iv1d(iv,ispin) =  in*ifftq2 + im*ifftq1 + il
       ENDDO
    ENDDO

    !
    !------------> LOOP OVER OCCUPIED STATES
    !
    DO  nu = 1 , ne
       !
       !---> FFT transform c_nu,k(g) --> psi_nu,k(r), for each k-point
       !                                              and each nu-state
       IF (noco%l_noco) THEN
          psi1r=0.0
          psi1i=0.0
          psi2r=0.0
          psi2i=0.0
          !------> map WF into FFTbox
#if ( !defined(CPP_INVERSION) )
          ! the compiler complains about the aimag if z is real, though these
          ! lines are never reached in a collinear calculation
          IF (noco%l_ss) THEN
             DO iv = 1 , lapw%nv(1)
                psi1r( iv1d(iv,1) )   = REAL( z(iv,nu) )
                psi1i( iv1d(iv,1) )   = AIMAG( z(iv,nu) )
             ENDDO
             DO iv = 1 , lapw%nv(2)
                psi2r( iv1d(iv,2) ) =  REAL(z(lapw%nv(1)+atoms%nlotot+iv,nu))
                psi2i( iv1d(iv,2) ) = AIMAG(z(lapw%nv(1)+atoms%nlotot+iv,nu))
             ENDDO
          ELSE
             DO iv = 1 , lapw%nv(jspin)
                psi1r( iv1d(iv,jspin) ) = REAL( z(iv,nu) )
                psi1i( iv1d(iv,jspin) ) = AIMAG( z(iv,nu) )
                psi2r(iv1d(iv,jspin))=REAL( z(lapw%nv(1)+atoms%nlotot+iv,nu))
                psi2i(iv1d(iv,jspin))=AIMAG(z(lapw%nv(1)+atoms%nlotot+iv,nu))
             ENDDO
          ENDIF
#endif
       ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     psir=0.0
#else
          psir=0.0
          psii=0.0
#endif
          !------> map WF into FFTbox
          DO iv = 1 , lapw%nv(jspin)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     psir( iv1d(iv,jspin) ) = z(iv,nu)
#else
             psir( iv1d(iv,jspin) ) =  REAL(z(iv,nu))
             psii( iv1d(iv,jspin) ) = AIMAG(z(iv,nu))
#endif
          ENDDO
       ENDIF
       !
       !------> do (real) inverse FFT; notice that the array psir is filled from
       !        0 to ifftq3-1, but starts at -ifftq2 to give work space for rfft
       !
       IF (noco%l_noco) THEN
          isn = 1

          CALL cfft(psi1r,psi1i,ifftq3,stars%kq1_fft,ifftq1,isn)
          CALL cfft(psi1r,psi1i,ifftq3,stars%kq2_fft,ifftq2,isn)
          CALL cfft(psi1r,psi1i,ifftq3,stars%kq3_fft,ifftq3,isn)

          CALL cfft(psi2r,psi2i,ifftq3,stars%kq1_fft,ifftq1,isn)
          CALL cfft(psi2r,psi2i,ifftq3,stars%kq2_fft,ifftq2,isn)
          CALL cfft(psi2r,psi2i,ifftq3,stars%kq3_fft,ifftq3,isn)
       ELSE
          isn = 1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     CALL rfft(isn,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft+1,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,&
     nw1,nw2,nw3,wsave,psir(ifftq3d), psir(-ifftq2))

     ! GM forces part
          IF (input%l_f) THEN
             DO in=-1,stars%kq3_fft,2
                DO im=0,ifftq2-1
                   ir = ifftq2 * in + im
                   ekin(ir) = ekin(ir) - wtf(nu) * eig(nu) *&
                        (psir(ir)**2 + psir(ir+ifftq2)**2)
                ENDDO
             ENDDO

             DO j = 1,3
                kpsir(ifftq3d:)=0.0
                kpsir(-ifftq2d:ifftq3d)=0.0
                DO iv = 1 , lapw%nv(jspin)
                   xk(1)=lapw%k1(iv,jspin)+bkpt(1)
                   xk(2)=lapw%k2(iv,jspin)+bkpt(2)
                   xk(3)=lapw%k3(iv,jspin)+bkpt(3)
                   s = 0.0
                   DO i = 1,3
                      s = s + xk(i)*cell%bmat(i,j)
                   ENDDO
                   kpsir( iv1d(iv,jspin) ) = s * z(iv,nu)
                ENDDO
                CALL rfft(isn,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft+1,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,&
                     nw1,nw2,nw3,wsave,kpsir(ifftq3d), kpsir(-ifftq2))
                DO in=-1,stars%kq3_fft,2
                   DO im=0,ifftq2-1
                      ir = ifftq2 * in + im
                      ekin(ir) = ekin(ir) + wtf(nu) * 0.5 * &
                           (kpsir(ir)**2 + kpsir(ir+ifftq2)**2)
                   ENDDO
                ENDDO
             ENDDO
          ENDIF
#else
          CALL cfft(psir,psii,ifftq3,stars%kq1_fft,ifftq1,isn)
          CALL cfft(psir,psii,ifftq3,stars%kq2_fft,ifftq2,isn)
          CALL cfft(psir,psii,ifftq3,stars%kq3_fft,ifftq3,isn)
          ! GM forces part
          IF (input%l_f) THEN
             DO ir = 0,ifftq3d-1
                ekin(ir) = ekin(ir) - wtf(nu)*eig(nu)*&
                     (psir(ir)**2+psii(ir)**2)
             ENDDO

             DO j = 1,3
                kpsir=0.0
                kpsii=0.0
                DO iv = 1 , lapw%nv(jspin)
                   xk(1)=lapw%k1(iv,jspin)+bkpt(1)
                   xk(2)=lapw%k2(iv,jspin)+bkpt(2)
                   xk(3)=lapw%k3(iv,jspin)+bkpt(3)
                   s = 0.0
                   DO i = 1,3
                      s = s + xk(i)*cell%bmat(i,j)
                   ENDDO
                   kpsir( iv1d(iv,jspin) ) = s *  REAL(z(iv,nu))
                   kpsii( iv1d(iv,jspin) ) = s * AIMAG(z(iv,nu))
                ENDDO

                CALL cfft(kpsir,kpsii,ifftq3,stars%kq1_fft,ifftq1,isn)
                CALL cfft(kpsir,kpsii,ifftq3,stars%kq2_fft,ifftq2,isn)
                CALL cfft(kpsir,kpsii,ifftq3,stars%kq3_fft,ifftq3,isn)

                DO ir = 0,ifftq3d-1
                   ekin(ir) = ekin(ir) + wtf(nu) * 0.5 *&
                        (kpsir(ir)**2+kpsii(ir)**2)
                ENDDO
             ENDDO
          ENDIF
#endif
       ENDIF
       !----> calculate rho(r) = sum w(k)*f(nu)*conjg(psi_nu,k(r))*psi_nu,k(r)
       !                         k,nu
       !      again, we fill rhon() from -ifftq2 to ifftq3-1 for the rfft
       !
       IF (noco%l_noco) THEN
          !--->             in the non-collinear case rho becomes a hermitian 2x2
          !--->             matrix (rhomat).
          DO ir = 0,ifftq3d-1
             rhomat(ir,1) = rhomat(ir,1) &
                  + wtf(nu)*( psi1r(ir)**2 + psi1i(ir)**2 )
             rhomat(ir,2) = rhomat(ir,2) &
                  + wtf(nu)*( psi2r(ir)**2 + psi2i(ir)**2 )
             rhomat(ir,3) = rhomat(ir,3) + wtf(nu)*&
                  (psi2r(ir)*psi1r(ir)+psi2i(ir)*psi1i(ir))
             rhomat(ir,4) = rhomat(ir,4) + wtf(nu)*&
                  (psi2r(ir)*psi1i(ir)-psi2i(ir)*psi1r(ir))
          ENDDO
          !--->             in a non-collinear calculation the interstitial charge
          !--->             cannot be calculated by a simple substraction if the
          !--->             muffin-tin (and vacuum) charge is know, because the
          !--->             total charge does not need to be one in each spin-
          !--->             channel. Thus it has to be calculated explicitly, if
          !--->             it is needed.
          IF (banddos%dos .OR. banddos%vacdos .OR. input%cdinf) THEN
             DO ir = 0,ifftq3d-1
                psi1r(ir) = (psi1r(ir)**2 + psi1i(ir)**2)
                psi2r(ir) = (psi2r(ir)**2 + psi2i(ir)**2)
             ENDDO
             isn = -1
             psi1i=0.0
             CALL cfft(psi1r,psi1i,ifftq3,stars%kq1_fft,ifftq1,isn)
             CALL cfft(psi1r,psi1i,ifftq3,stars%kq2_fft,ifftq2,isn)
             CALL cfft(psi1r,psi1i,ifftq3,stars%kq3_fft,ifftq3,isn)
             psi2i=0.0
             CALL cfft(psi2r,psi2i,ifftq3,stars%kq1_fft,ifftq1,isn)
             CALL cfft(psi2r,psi2i,ifftq3,stars%kq2_fft,ifftq2,isn)
             CALL cfft(psi2r,psi2i,ifftq3,stars%kq3_fft,ifftq3,isn)
             cwk=0.0
             DO ik = 0 , stars%kmxq_fft - 1
                cwk(stars%igfft(ik,1))=cwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
                     CMPLX(psi1r(igq_fft(ik)),psi1i(igq_fft(ik)))
             ENDDO
             DO istr = 1,stars%ng3_fft
                CALL pwint(&
                     stars,atoms,sym,&
                     oneD,cell,stars%kv3(1,istr),x)
                qis(nu,ikpt,1) = qis(nu,ikpt,1)&
                     + REAL(cwk(istr)*x)/cell%omtil/REAL(ifftq3)
             ENDDO

             cwk=0.0
             DO ik = 0 , stars%kmxq_fft - 1
                cwk(stars%igfft(ik,1))=cwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
                     CMPLX(psi2r(igq_fft(ik)),psi2i(igq_fft(ik)))
             ENDDO
             DO istr = 1,stars%ng3_fft
                CALL pwint(&
                     stars,atoms,sym,&
                     oneD,cell,&
                     stars%kv3(1,istr),&
                     x)
                qis(nu,ikpt,input%jspins) = qis(nu,ikpt,input%jspins)&
                     + REAL(cwk(istr)*x)/cell%omtil/REAL(ifftq3)
             ENDDO
          ENDIF
       ELSE
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     DO in=-1,stars%kq3_fft,2
          DO im=0,ifftq2-1
             ir = ifftq2 * in + im
             rhon(ir) = rhon(ir) + wtf(nu) * ( psir(ir)**2 +&
                  psir(ir+ifftq2)**2 )
          ENDDO
       ENDDO
#else
       DO ir = 0,ifftq3d-1
          rhon(ir)=rhon(ir)+wtf(nu)*(psir(ir)**2+psii(ir)**2)
       ENDDO
#endif
    ENDIF
    !              DO ir = -ifftq2 , ifftq3-1
    !     +                      + wtf(nu)*(psi(ir+ifftq3d) * psi(ir+ifftq3d)
    !     +                               + psi(ir  ) * psi(ir  )
    !     +                                 )
    !              ENDDO

 ENDDO
 !
 !<<<<<<<<<<<<<< END OUTER LOOP OVER STATES NU  >>>>>>>>>>>>>>>>>>
 !
 !
 !----> perform back  FFT transform: rho(r) --> chgn(star)
 !        ( do direct FFT)                    = cwk(star)

 !--->  In a collinear calculation pwden is calles once per spin.
 !--->  However in a non-collinear calculation pwden is only called once
 !--->  and all four components of the density matrix (n_11 n_22 n_12
 !--->  n_21) have to be calculated at once.
 ndens = 1
 IF (noco%l_noco) ndens = 4
 DO idens = 1,ndens
    IF (noco%l_noco) THEN
       psi1r=0.0
       isn = -1
       CALL cfft(rhomat(0,idens),psi1r,ifftq3,stars%kq1_fft,ifftq1,isn)
       CALL cfft(rhomat(0,idens),psi1r,ifftq3,stars%kq2_fft,ifftq2,isn)
       CALL cfft(rhomat(0,idens),psi1r,ifftq3,stars%kq3_fft,ifftq3,isn)
    ELSE
       !--->  psir is used here as work array, charge is real ,but fft complex
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     psir(ifftq3d:)=0.0
       IF (input%l_f) kpsir(ifftq3d:)=0.0
#else
       psir=0.0
       psii=0.0
       IF (input%l_f) kpsir=0.0
       IF (input%l_f) kpsii=0.0
#endif

       isn = -1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     CALL rfft(isn,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft+1,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,&
     stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,wsave,psir(ifftq3d), rhon(-ifftq2))
       IF (input%l_f) CALL rfft(isn,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft+1,stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,&
            stars%kq1_fft,stars%kq2_fft,stars%kq3_fft,wsave,kpsir(ifftq3d), ekin(-ifftq2))
#else
       CALL cfft(rhon,psir,ifftq3,stars%kq1_fft,ifftq1,isn)
       CALL cfft(rhon,psir,ifftq3,stars%kq2_fft,ifftq2,isn)
       CALL cfft(rhon,psir,ifftq3,stars%kq3_fft,ifftq3,isn)
       !+apw
       IF (input%l_f) THEN 
          CALL cfft(ekin,psii,ifftq3,stars%kq1_fft,ifftq1,isn)
          CALL cfft(ekin,psii,ifftq3,stars%kq2_fft,ifftq2,isn)
          CALL cfft(ekin,psii,ifftq3,stars%kq3_fft,ifftq3,isn)
       ENDIF
#endif
    ENDIF
    !  ---> collect stars from the fft-grid
    !
    cwk=0.0
    ecwk=0.0
    IF (noco%l_noco) THEN
       DO ik = 0 , stars%kmxq_fft - 1
          cwk(stars%igfft(ik,1))=cwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
               CMPLX(rhomat(igq_fft(ik),idens),psi1r(igq_fft(ik)))
       ENDDO
    ELSE
       DO ik = 0 , stars%kmxq_fft - 1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     cwk(stars%igfft(ik,1))=cwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
     CMPLX(rhon(igq_fft(ik)),zero)
#else
          cwk(stars%igfft(ik,1))=cwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
               CMPLX(rhon(igq_fft(ik)),psir(igq_fft(ik)))
#endif
       ENDDO
       !+apw
       IF (input%l_f) THEN 
          DO ik = 0 , stars%kmxq_fft - 1
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
     ecwk(stars%igfft(ik,1))=ecwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
     CMPLX(ekin(igq_fft(ik)),zero)
#else
             ecwk(stars%igfft(ik,1))=ecwk(stars%igfft(ik,1))+CONJG(stars%pgfft(ik))*&
                  CMPLX(ekin(igq_fft(ik)),psii(igq_fft(ik)))
#endif
          ENDDO
       ENDIF
       !-apw
    ENDIF
    !
    scale=1.0/ifftq3
    DO istr = 1 , stars%ng3_fft
       cwk(istr) = scale * cwk(istr) / REAL( stars%nstr(istr) )
    ENDDO
    IF (input%l_useapw) THEN

       IF (input%l_f) THEN
          DO istr = 1 , stars%ng3_fft
             ecwk(istr) = scale * ecwk(istr) / REAL( stars%nstr(istr) )
          ENDDO
          CALL force_b8(&
               atoms,ecwk,stars,&
               sym,cell,&
               jspin,&
               forces,f_b8)
       ENDIF
    ENDIF
    !
    !---> check charge neutralilty
    !
    IF ((idens.EQ.1).OR.(idens.EQ.2)) THEN
       IF (noco%l_noco) THEN
          IF (idens.EQ.1) THEN
             q0 = q0_11
          ELSE
             q0 = q0_22
          ENDIF
       ENDIF
       IF ( ABS( q0 ) .GT. 1.0e-9) THEN
          IF ( ABS( q0 - REAL(cwk(1)) )/q0 .GT. tol_3 ) THEN
             WRITE(99,*) "XX:",ne,lapw%nv
             DO istr=1,SIZE(z,2)
                WRITE(99,*) "X:",istr,z(:,istr)
             ENDDO
             WRITE ( 6,'(''bad quality of charge density'',&
                  2f13.8)')q0, REAL( cwk(1) )
             WRITE (16,'(''bad quality of charge density'',&
                  2f13.8)')q0, REAL( cwk(1) )
             CALL juDFT_warn('pwden: bad quality of charge')
          ENDIF
       ENDIF
    ENDIF
    !
    !---> add charge density to existing one
    !
    IF(idens.LE.2) THEN
       !--->       add to spin-up or -down density (collinear & non-collinear)
       ispin = jspin
       IF (noco%l_noco) ispin = idens
       DO istr = 1 , stars%ng3_fft
          qpw(istr,ispin) = qpw(istr,ispin) + cwk(istr)
       ENDDO
    ELSE IF (idens.EQ.3) THEN
       !--->       add to off-diag. part of density matrix (only non-collinear)
       DO istr = 1 , stars%ng3_fft
          cdom(istr) = cdom(istr) + cwk(istr)
       ENDDO
    ELSE
       !--->       add to off-diag. part of density matrix (only non-collinear)
       DO istr = 1 , stars%ng3_fft
          cdom(istr) = cdom(istr) + CMPLX(0.0,1.0)*cwk(istr)
       ENDDO
    ENDIF

 ENDDO

 DEALLOCATE(cwk,ecwk)

 IF (noco%l_noco) THEN
    DEALLOCATE ( psi1r,psi1i,psi2r,psi2i,rhomat )
 ELSE
    DEALLOCATE ( psir,psii,rhon )
    IF (input%l_f) DEALLOCATE ( kpsir,kpsii,ekin)
 ENDIF

END SUBROUTINE pwden
END MODULE m_pwden
