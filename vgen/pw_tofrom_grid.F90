!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_pw_tofrom_grid
  USE m_types
  PRIVATE
  REAL,PARAMETER:: d_15=1.e-15
  COMPLEX,PARAMETER:: ci=CMPLX(0.,1.)

  INTEGER :: ifftd,ifftxc3d,ifftxc3
  !----->  fft  information  for xc potential + energy
  INTEGER, ALLOCATABLE :: igxc_fft(:)
  REAL,    ALLOCATABLE :: gxc_fft(:,:)
  
  PUBLIC :: init_pw_grid,pw_to_grid,pw_from_grid,finish_pw_grid
CONTAINS
  SUBROUTINE init_pw_grid(stars,sym,cell)
    USE m_prpxcfftmap
    USE m_types
     IMPLICIT NONE
      TYPE(t_stars),INTENT(IN)      :: stars
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_cell),INTENT(IN)       :: cell
      
    !---> set up pointer for backtransformation of from g-vector in
    !     positive domain of xc density fftbox into stars.
    !     also the x,y,z components of the g-vectors are set up to calculate
    !     derivatives.
    !     in principle this can also be done in main program once.
    !     it is done here to save memory.
    !
    ifftd=27*stars%mx1*stars%mx2*stars%mx3
    ifftxc3d = stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft
    ALLOCATE ( igxc_fft(0:ifftxc3d-1),gxc_fft(0:ifftxc3d-1,3) )
    CALL prp_xcfft_map(stars,sym, cell, igxc_fft,gxc_fft)
    !
    ifftxc3=stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft
       
  END SUBROUTINE init_pw_grid
  
  SUBROUTINE pw_to_grid(xcpot,input,noco,stars,cell,den,rho,grad)
    !.....------------------------------------------------------------------
    !------->          abbreviations
    !
    !     ph_wrk: work array containing phase * g_x,gy...... 
    !     den%pw: charge density stored as stars
    !     rho   : charge density stored in real space
    !     v_xc   : exchange-correlation potential in real space
    !     exc   : exchange-correlation energy density in real space
    !     kxc1d  : dimension of the charge density fft box in the pos. domain
    !     kxc2d  : defined in dimens.f program (subroutine apws).1,2,3 indic
    !     kxc3d  ; a_1, a_2, a_3 directions.
    !     kq(i) : i=1,2,3 actual length of the fft-box for which fft is done
    !     nstr  : number of members (arms) of reciprocal lattice (g) vector
    !             of each star.
    !     nxc3_fft: number of stars in the  charge density  fft-box
    !     ng3   : number of 3 dim. stars in the charge density sphere define
    !             by gmax
    !     kmxxc_fft: number of g-vectors forming the nxc3_fft stars in the
    !               charge density or xc-density sphere
    !     kimax : number of g-vectors forming the ng3 stars in the gmax-sphe
    !     ifftxc3d: elements (g-vectors) in the charge density  fft-box
    !     igfft : pointer from the g-sphere (stored as stars) to fft-grid
    !             and     from fft-grid to g-sphere (stored as stars)
    !     pgfft : contains the phases of the g-vectors of sph.
    !     isn   : isn = +1, fft transform for g-space to r-space
    !             isn = -1, vice versa
    !
    !-------------------------------------------------------------------
    USE m_grdrsis
    USE m_mkgxyz3
    USE m_fft3dxc
    USE m_types
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_potden),INTENT(IN)     :: den
    REAL,ALLOCATABLE,INTENT(out)  :: rho(:,:)
    TYPE(t_gradients),INTENT(OUT)  :: grad


    INTEGER      :: js,i,idm,ig,ndm,jdm
    REAL         :: rhotot
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: cqpw(:,:),ph_wrk(:)
    REAL,    ALLOCATABLE :: bf3(:)
    REAL,    ALLOCATABLE :: rhd1(:,:,:),rhd2(:,:,:)
    REAL,    ALLOCATABLE :: mx(:),my(:)
    REAL,    ALLOCATABLE :: magmom(:),dmagmom(:,:),ddmagmom(:,:,:) 
    
    ! Allocate arrays
    ALLOCATE( bf3(0:ifftd-1), rho(0:ifftxc3d-1,input%jspins))
    IF (xcpot%is_gga()) ALLOCATE( ph_wrk(0:ifftxc3d-1),rhd1(0:ifftxc3d-1,input%jspins,3), &
         rhd2(0:ifftxc3d-1,input%jspins,6) )
    IF (noco%l_noco)  THEN
       ALLOCATE( mx(0:ifftxc3-1),my(0:ifftxc3-1),magmom(0:ifftxc3-1))
       IF (xcpot%is_gga()) ALLOCATE(dmagmom(0:ifftxc3-1,3),ddmagmom(0:ifftxc3-1,3,3) )
    END IF
    !Put den%pw on grid and store into rho(:,1:2)
    DO js=1,input%jspins
       CALL fft3dxc(rho(0:,js),bf3, den%pw(:,js), stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
    END DO

    IF (noco%l_noco) THEN  
       !  Get mx,my on real space grid and recalculate rho and magmom
       CALL fft3dxc(mx,my, den%pw(:,3), stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
       DO i=0,ifftxc3-1 
          rhotot= 0.5*( rho(i,1) + rho(i,2) )
          magmom(i)= SQRT(  (0.5*(rho(i,1)-rho(i,2)))**2 + mx(i)**2 + my(i)**2 )
          rho(i,1)= rhotot+magmom(i)
          rho(i,2)= rhotot-magmom(i)
       END DO
    ENDIF

    IF (xcpot%is_gga()) THEN  

    ! In collinear calculations all derivatives are calculated in g-spce,
    ! in non-collinear calculations the derivatives of |m| are calculated in real space. 

    !-->   for d(rho)/d(x,y,z) = rhd1(:,:,idm) (idm=1,2,3).
    !
    !         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z).

       ALLOCATE(cqpw(stars%ng3,input%jspins))

       cqpw(:,:)= ci*den%pw(:,:)
   
       DO idm=1,3
          DO ig = 0 , stars%kmxxc_fft - 1
             ph_wrk(ig) = stars%pgfft(ig) * gxc_fft(ig,idm)
          END DO

          DO js=1,input%jspins
             CALL fft3dxc(rhd1(0:,js,idm),bf3, cqpw(:,js), stars%kxc1_fft,stars%kxc2_fft,&
                  stars%kxc3_fft,stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,ph_wrk,stars%nstr)
          END DO
       END DO

       IF (noco%l_noco) THEN

          CALL grdrsis(magmom,cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,dmagmom )

          DO i=0,ifftxc3-1
             DO idm=1,3
                rhotot= rhd1(i,1,idm)/2.+rhd1(i,2,idm)/2.
                rhd1(i,1,idm)= rhotot+dmagmom(i,idm) 
                rhd1(i,2,idm)= rhotot-dmagmom(i,idm) 
             END DO
          END DO
       END IF

       !-->   for dd(rho)/d(xx,xy,yy,zx,yz,zz) = rhd2(:,:,idm) (idm=1,2,3,4,5,6)
       !
       !         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z) * g_(x,y,z)

       cqpw(:,:)= -den%pw(:,:)
   
       ndm = 0
       DO idm = 1,3
          DO jdm = 1,idm
             ndm = ndm + 1
             DO ig = 0 , stars%kmxxc_fft-1
                ph_wrk(ig) = stars%pgfft(ig)*gxc_fft(ig,idm)*gxc_fft(ig,jdm)
             ENDDO
             
             DO js=1,input%jspins
                CALL fft3dxc(rhd2(0:,js,ndm),bf3, cqpw(:,js), stars%kxc1_fft,stars%kxc2_fft,&
                     stars%kxc3_fft,stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,ph_wrk,stars%nstr)
             END DO
          END DO ! jdm 
       END DO   ! idm 

       DEALLOCATE(cqpw)

       IF (noco%l_noco) THEN
          DO idm = 1,3
             CALL grdrsis(dmagmom(0,idm),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,ddmagmom(0,1,idm) )
          END DO
          ndm= 0
          DO idm = 1,3
             DO jdm = 1,idm
                ndm = ndm + 1  
                DO i=0,ifftxc3-1
                   rhotot= rhd2(i,1,ndm)/2.+rhd2(i,2,ndm)/2.
                   rhd2(i,1,ndm)= rhotot + ( ddmagmom(i,jdm,idm) + ddmagmom(i,idm,jdm) )/2. 
                   rhd2(i,2,ndm)= rhotot - ( ddmagmom(i,jdm,idm) + ddmagmom(i,idm,jdm) )/2. 
                END DO
             ENDDO !jdm
          ENDDO   !idm 
       END IF
       CALL xcpot%alloc_gradients(ifftxc3d,input%jspins,grad)
 
       !
       !     calculate the quantities such as abs(grad(rho)),.. used in
       !     evaluating the gradient contributions to potential and energy.
       !
       CALL mkgxyz3 (ifftxc3d,input%jspins,ifftxc3,input%jspins,rho, rhd1(0,1,1),rhd1(0,1,2),rhd1(0,1,3),&
            rhd2(0,1,1),rhd2(0,1,3),rhd2(0,1,6), rhd2(0,1,5),rhd2(0,1,4),rhd2(0,1,2), grad)
       
    ENDIF
    rho(i,js)=MAX(rho(i,js),d_15)
   
  END SUBROUTINE pw_to_grid


  SUBROUTINE pw_from_grid(stars,l_pw_w,v_in,v_out)
    USE m_fft3d
    USE m_fft3dxc
    USE m_types
    IMPLICIT NONE
    TYPE(t_stars),INTENT(IN)      :: stars
    REAL,INTENT(IN)               :: v_in(:,:)
    LOGICAL,INTENT(in)            :: l_pw_w
    TYPE(t_potden),INTENT(INOUT)  :: v_out
    
    
    INTEGER              :: js,k,i
    REAL,ALLOCATABLE     :: bf3(:),vcon(:) 
    COMPLEX, ALLOCATABLE :: fg3(:)
    ALLOCATE( bf3(0:ifftd-1),fg3(stars%ng3))
    DO js = 1,SIZE(v_in,2)
       bf3=0.0
       CALL fft3dxc(v_in(0:,js),bf3, fg3, stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
            stars%nxc3_fft,stars%kmxxc_fft,-1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
       
       DO k = 1,stars%nxc3_fft
          v_out%pw(k,js) = v_out%pw(k,js) + fg3(k)
       ENDDO

       IF (l_pw_w) THEN
          !----> Perform fft transform: v_xc(star) --> vxc(r) 
          !     !Use large fft mesh for convolution
          fg3(stars%nxc3_fft+1:)=0.0
          ALLOCATE ( vcon(0:ifftd-1) )
          CALL fft3d(vcon(0),bf3, fg3, stars,+1)
          !
          !----> Convolute with step function
          !
          DO i=0,ifftd-1
             vcon(i)=stars%ufft(i)*vcon(i)
          ENDDO
          bf3=0.0
          CALL fft3d(vcon(0),bf3, fg3, stars,-1,.FALSE.)
          DEALLOCATE ( vcon )
          !
          !----> add to warped coulomb potential
          !
          DO k = 1,stars%ng3
             v_out%pw_w(k,js) = v_out%pw_w(k,js) + fg3(k)
          ENDDO
       ENDIF
    END DO
  END SUBROUTINE pw_from_grid
    
  SUBROUTINE finish_pw_grid()
    IMPLICIT NONE
    DEALLOCATE(igxc_fft,gxc_fft)
  END SUBROUTINE finish_pw_grid

END MODULE m_pw_tofrom_grid
