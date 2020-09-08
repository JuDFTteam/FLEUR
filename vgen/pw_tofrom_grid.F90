!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_pw_tofrom_grid
   USE m_types
   PRIVATE
   REAL,PARAMETER:: d_15=1.e-15

   INTEGER :: ifftd,ifftxc3
   !----->  fft  information  for xc potential + energy
   INTEGER, ALLOCATABLE :: igxc_fft(:)
   REAL,    ALLOCATABLE :: gxc_fft(:,:) !gxc_fft(ig,idm)

   PUBLIC :: init_pw_grid, pw_to_grid, pw_from_grid, finish_pw_grid
CONTAINS
  SUBROUTINE init_pw_grid(dograds,stars,sym,cell)
    USE m_prpxcfftmap
    USE m_types
    IMPLICIT NONE
    LOGICAL,INTENT(IN)            :: dograds
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell

      !---> set up pointer for backtransformation of from g-vector in
      !     positive domain of xc density fftbox into stars.
      !     also the x,y,z components of the g-vectors are set up to calculate
      !     derivatives.
      !     in principle this can also be done in main program once.
      !     it is done here to save memory.

    ifftd=27*stars%mx1*stars%mx2*stars%mx3
    ifftxc3  = stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft
    IF (dograds) THEN
       CALL prp_xcfft_map(stars,sym, cell, igxc_fft,gxc_fft)
    ENDIF

  END SUBROUTINE init_pw_grid

  SUBROUTINE pw_to_grid(dograds,jspins,l_noco,stars,cell,den_pw,grad,xcpot,rho)
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
    !     igfft : pointer from the g-sphere (stored as stars) to fft-grid
    !             and     from fft-grid to g-sphere (stored as stars)
    !     pgfft : contains the phases of the g-vectors of sph.
    !     isn   : isn = +1, fft transform for g-space to r-space
    !             isn = -1, vice versa
    !     l_rdm : if true, calculate noco gradients in local frame of density matrix
    !             improves convergence in case of GGA + noco (best with large G_max_xc)
    !
    !-------------------------------------------------------------------
    USE m_grdrsis
    USE m_polangle
    USE m_mkgxyz3
    USE m_fft3dxc
    USE m_fft3d
    USE m_types
    USE m_constants
    IMPLICIT NONE

    LOGICAL,INTENT(IN)                    :: dograds
    INTEGER,INTENT(IN)                    :: jspins
    LOGICAL,INTENT(IN)                    :: l_noco
    TYPE(t_stars),INTENT(IN)              :: stars
    TYPE(t_cell),INTENT(IN)               :: cell
    COMPLEX,INTENT(IN)                    :: den_pw(:,:)
    TYPE(t_gradients),INTENT(OUT)         :: grad
    CLASS(t_xcpot), INTENT(IN),OPTIONAL   :: xcpot
    REAL,ALLOCATABLE,INTENT(OUT),OPTIONAL :: rho(:,:)


    INTEGER      :: js,i,idm,ig,ndm,jdm,j
    REAL         :: rhotot,mmx,mmy,mmz,theta,phi
    COMPLEX      :: ci,rho21
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: cqpw(:,:),ph_wrk(:)
    REAL,    ALLOCATABLE :: bf3(:)
    REAL,    ALLOCATABLE :: rhd1(:,:,:),rhd2(:,:,:)
    REAL,    ALLOCATABLE :: mx(:),my(:)
    REAL,    ALLOCATABLE :: magmom(:),dmagmom(:,:),ddmagmom(:,:,:)
    REAL,    ALLOCATABLE :: rhodiag(:,:),der(:,:,:),dder(:,:,:,:)
    REAL,    ALLOCATABLE :: sinsqu(:),cossqu(:),sincos(:),rhdd(:,:,:,:)
    COMPLEX, ALLOCATABLE :: exi(:)

    LOGICAL, PARAMETER :: l_rdm=.true.
    ci=cmplx(0.,1.)

    ! Allocate arrays
    ALLOCATE( bf3(0:ifftd-1))
    IF (dograds) THEN
       IF (PRESENT(rho)) ALLOCATE(rho(0:ifftxc3-1,jspins))
       ALLOCATE( ph_wrk(0:ifftxc3-1),rhd1(0:ifftxc3-1,jspins,3))
       ALLOCATE( rhd2(0:ifftxc3-1,jspins,6) )
     ELSE
        IF (PRESENT(rho)) ALLOCATE(rho(0:ifftd-1,jspins))
     ENDIF
    IF (l_noco)  THEN
       IF (dograds) THEN
          ALLOCATE( mx(0:ifftxc3-1),my(0:ifftxc3-1),magmom(0:ifftxc3-1))
          IF (l_rdm) THEN
           ALLOCATE( rhodiag(0:ifftxc3-1,jspins),der(0:ifftxc3-1,3,4),dder(0:ifftxc3-1,3,3,4),rhdd(0:ifftxc3-1,2,3,3) )
           ALLOCATE( sinsqu(0:ifftxc3-1),cossqu(0:ifftxc3-1),sincos(0:ifftxc3-1),exi(0:ifftxc3-1) )
          ELSE
            ALLOCATE( dmagmom(0:ifftxc3-1,3),ddmagmom(0:ifftxc3-1,3,3) )
          ENDIF
       ELSE
          ALLOCATE( mx(0:ifftd-1),my(0:ifftd-1),magmom(0:ifftd-1))
       ENDIF
    END IF

    IF (PRESENT(rho)) THEN
    !Put den_pw on grid and store into rho(:,1:2)
       DO js=1,jspins
          IF (dograds) THEN
             CALL fft3dxc(rho(0:,js),bf3, den_pw(:,js), stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
                  stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
          ELSE
             CALL fft3d(rho(0,js),bf3, den_pw(:,js), stars,+1)
          ENDIF
       END DO

       IF (l_noco) THEN
          !  Get mx,my on real space grid and recalculate rho and magmom
          IF (dograds) THEN
             CALL fft3dxc(mx,my, den_pw(:,3), stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
                  stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
          ELSE
             CALL fft3d(mx,my, den_pw(:,3), stars,+1)
          ENDIF
          DO i=0,MIN(SIZE(rho,1),size(mx))-1
             rhotot= 0.5*( rho(i,1) + rho(i,2) )
             magmom(i)= SQRT(  (0.5*(rho(i,1)-rho(i,2)))**2 + mx(i)**2 + my(i)**2 )
             IF (l_rdm.AND.dograds) THEN
                rhodiag(i,1) = rho(i,1)
                rhodiag(i,2) = rho(i,2)
             ENDIF
             rho(i,1)= rhotot+magmom(i)
             rho(i,2)= rhotot-magmom(i)
             IF (l_rdm.AND.dograds) THEN              ! prepare rotation matrix
                mmx = 2.0 * mx(i)
                mmy = 2.0 * my(i)
                mmz = rhodiag(i,1) - rhodiag(i,2)
                CALL pol_angle(mmx,mmy,mmz,theta,phi)
                cossqu(i) = cos(0.5*theta)**2
                sinsqu(i) = sin(0.5*theta)**2
                sincos(i) = 2.0*sin(0.5*theta)*cos(0.5*theta)
                exi(i)    = exp(-ci*phi)
             ENDIF
          END DO
       ENDIF
    ENDIF
    IF (dograds) THEN

    ! In collinear calculations all derivatives are calculated in g-spce,
    ! in non-collinear calculations the derivatives of |m| are calculated in real space.

    !-->   for d(rho)/d(x,y,z) = rhd1(:,:,idm) (idm=1,2,3).
    !
    !         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z).

       ALLOCATE(cqpw(stars%ng3,jspins))

       cqpw(:,:)= ImagUnit*den_pw(:,:jspins)

       DO idm=1,3
          DO ig = 0 , stars%kmxxc_fft - 1
             ph_wrk(ig) = stars%pgfft(ig) * gxc_fft(ig,idm)
          END DO

          DO js=1,jspins
             CALL fft3dxc(rhd1(0:,js,idm),bf3, cqpw(:,js), stars%kxc1_fft,stars%kxc2_fft,&
                  stars%kxc3_fft,stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,ph_wrk,stars%nstr)
          END DO
       END DO

       IF (l_noco) THEN

          IF (l_rdm) THEN

             CALL grdrsis( rhodiag(0,1),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,der(0,1,1) )
             CALL grdrsis( rhodiag(0,2),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,der(0,1,2) )
             CALL grdrsis( mx(0),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,der(0,1,3) )
             CALL grdrsis( my(0),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,der(0,1,4) )

             DO i=0,ifftxc3-1 ! project on magnetization axis
                DO idm=1,3
                   rho21 = der(i,idm,3) + ci * der(i,idm,4)
                   rhd1(i,1,idm) = cossqu(i) * der(i,idm,1) + &
                                   sincos(i) * real( exi(i)*rho21 ) + &
                                   sinsqu(i) * der(i,idm,2)
                   rhd1(i,2,idm) = sinsqu(i) * der(i,idm,1) - &
                                   sincos(i) * real( exi(i)*rho21 ) + &
                                   cossqu(i) * der(i,idm,2)
                ENDDO
             ENDDO

          ELSE

             CALL grdrsis(magmom,cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,dmagmom )

             DO i=0,ifftxc3-1
                DO idm=1,3
                   rhotot= rhd1(i,1,idm)/2.+rhd1(i,2,idm)/2.
                   rhd1(i,1,idm)= rhotot+dmagmom(i,idm)
                   rhd1(i,2,idm)= rhotot-dmagmom(i,idm)
                END DO
             END DO

          END IF

       END IF

       !-->   for dd(rho)/d(xx,xy,yy,zx,yz,zz) = rhd2(:,:,idm) (idm=1,2,3,4,5,6)
       !
       !         ph_wrk: exp(i*(g_x,g_y,g_z)*tau) * g_(x,y,z) * g_(x,y,z)

       cqpw(:,:)= -den_pw(:,:jspins)

       ndm = 0
       DO idm = 1,3
          DO jdm = 1,idm
             ndm = ndm + 1
             DO ig = 0 , stars%kmxxc_fft-1
                ph_wrk(ig) = stars%pgfft(ig)*gxc_fft(ig,idm)*gxc_fft(ig,jdm)
             ENDDO

             DO js=1,jspins
                CALL fft3dxc(rhd2(0:,js,ndm),bf3, cqpw(:,js), stars%kxc1_fft,stars%kxc2_fft,&
                     stars%kxc3_fft,stars%nxc3_fft,stars%kmxxc_fft,+1, stars%igfft(0:,1),igxc_fft,ph_wrk,stars%nstr)
             END DO
          END DO ! jdm
       END DO   ! idm

       DEALLOCATE(cqpw)

       IF (l_noco) THEN

          IF (l_rdm) THEN
             DO idm = 1,3
               CALL grdrsis( der(0,idm,1),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft, dder(0,1,idm,1) )
               CALL grdrsis( der(0,idm,2),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft, dder(0,1,idm,2) )
               CALL grdrsis( der(0,idm,3),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft, dder(0,1,idm,3) )
               CALL grdrsis( der(0,idm,4),cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft, dder(0,1,idm,4) )
               DO i=0,ifftxc3-1 ! project on magnetization axis
                 DO jdm=1,3
                   rho21 = dder(i,jdm,idm,3) + ci * dder(i,jdm,idm,4)
                   rhdd(i,1,jdm,idm) = cossqu(i) * dder(i,jdm,idm,1) + &
                                   sincos(i) * real( exi(i)*rho21 ) + &
                                   sinsqu(i) * dder(i,jdm,idm,2)
                   rhdd(i,2,jdm,idm) = sinsqu(i) * dder(i,jdm,idm,1) - &
                                   sincos(i) * real( exi(i)*rho21 ) + &
                                   cossqu(i) * dder(i,jdm,idm,2)
                 ENDDO
               ENDDO
             ENDDO
             DO j=1,2
               DO i=0,ifftxc3-1
                 rhd2(i,j,1) = rhdd(i,j,1,1)
                 rhd2(i,j,2) = 0.5*(rhdd(i,j,1,2)+rhdd(i,j,2,1)) ! xy - averaging should be unneccessary
                 rhd2(i,j,3) = rhdd(i,j,2,2)
                 rhd2(i,j,4) = 0.5*(rhdd(i,j,1,3)+rhdd(i,j,3,1)) ! zx
                 rhd2(i,j,5) = 0.5*(rhdd(i,j,2,3)+rhdd(i,j,3,2)) ! yz
                 rhd2(i,j,6) = rhdd(i,j,3,3)
               ENDDO
             ENDDO
             DEALLOCATE (rhodiag,der,rhdd,dder,sinsqu,cossqu,sincos,exi)
          ELSE

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
             DEALLOCATE(dmagmom,ddmagmom)
          END IF
       END IF

       IF (PRESENT(xcpot)) THEN
          CALL xcpot%alloc_gradients(ifftxc3,jspins,grad)
       END IF

       !
       !     calculate the quantities such as abs(grad(rho)),.. used in
       !     evaluating the gradient contributions to potential and energy.
       !
       IF (PRESENT(rho)) THEN
          CALL mkgxyz3 (rho,rhd1(0:,:,1),rhd1(0:,:,2),rhd1(0:,:,3),&
               rhd2(0:,:,1),rhd2(0:,:,3),rhd2(0:,:,6), rhd2(0:,:,5),rhd2(0:,:,4),rhd2(0:,:,2),0,grad)
       ELSE
          !Dummy rho (only possible if grad is used for libxc mode)
          !CALL mkgxyz3 (RESHAPE((/0.0/),(/1,1/)),rhd1(0:,:,1),rhd1(0:,:,2),rhd1(0:,:,3),&
          !     rhd2(0:,:,1),rhd2(0:,:,3),rhd2(0:,:,6), rhd2(0:,:,5),rhd2(0:,:,4),rhd2(0:,:,2),grad)

          IF (dograds.and.(.not.PRESENT(xcpot))) THEN
             ALLOCATE(grad%gr(3,ifftxc3,1))
          END IF

          CALL mkgxyz3 (0*rhd1(0:,:,1),rhd1(0:,:,1),rhd1(0:,:,2),rhd1(0:,:,3),&
               rhd2(0:,:,1),rhd2(0:,:,3),rhd2(0:,:,6), rhd2(0:,:,5),rhd2(0:,:,4),rhd2(0:,:,2),0,grad)
       END IF

    ENDIF
    IF (PRESENT(rho)) THEN
       WHERE(ABS(rho) < d_15) rho = d_15
    ENDIF

  END SUBROUTINE pw_to_grid


  SUBROUTINE pw_from_grid(dograds,stars,l_pw_w,v_in,v_out_pw,v_out_pw_w)
    USE m_fft3d
    USE m_fft3dxc
    USE m_types
    IMPLICIT NONE
    LOGICAL,INTENT(IN)            :: dograds
    TYPE(t_stars),INTENT(IN)      :: stars
    REAL,INTENT(INOUT)            :: v_in(0:,:)
    LOGICAL,INTENT(in)            :: l_pw_w
    COMPLEX,INTENT(INOUT)         :: v_out_pw(:,:)
    COMPLEX,INTENT(INOUT),OPTIONAL:: v_out_pw_w(:,:)


    INTEGER              :: js,k,i
    REAL,ALLOCATABLE     :: bf3(:),vcon(:)
    COMPLEX, ALLOCATABLE :: fg3(:)
    ALLOCATE( bf3(0:ifftd-1),fg3(stars%ng3))
    ALLOCATE ( vcon(0:ifftd-1) )
    DO js = 1,SIZE(v_in,2)
       bf3=0.0
       IF (dograds) THEN
          CALL fft3dxc(v_in(0:,js),bf3, fg3, stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,&
               stars%nxc3_fft,stars%kmxxc_fft,-1, stars%igfft(0:,1),igxc_fft,stars%pgfft,stars%nstr)
       ELSE
          vcon(0:)=v_in(0:,js)
          CALL fft3d(v_in(0:,js),bf3, fg3, stars,-1)
       ENDIF
       DO k = 1,MERGE(stars%nxc3_fft,stars%ng3,dograds)
          v_out_pw(k,js) = v_out_pw(k,js) + fg3(k)
       ENDDO

       IF (l_pw_w) THEN
          IF (dograds) THEN
             !----> Perform fft transform: v_xc(star) --> vxc(r)
             !     !Use large fft mesh for convolution
             fg3(stars%nxc3_fft+1:)=0.0
             CALL fft3d(vcon(0),bf3, fg3, stars,+1)
          ENDIF
          !
          !----> Convolute with step function
          !
          DO i=0,ifftd-1
             vcon(i)=stars%ufft(i)*vcon(i)
          ENDDO
          bf3=0.0
          CALL fft3d(vcon(0),bf3, fg3, stars,-1)
          fg3=fg3*stars%nstr
          !
          !----> add to warped coulomb potential
          !
          IF (PRESENT(v_out_pw_w)) THEN
             DO k = 1,stars%ng3
                v_out_pw_w(k,js) = v_out_pw_w(k,js) + fg3(k)
             ENDDO
          END IF
       ENDIF
    END DO
  END SUBROUTINE pw_from_grid

  SUBROUTINE finish_pw_grid()
    IMPLICIT NONE
    IF (ALLOCATED(igxc_fft)) DEALLOCATE(igxc_fft,gxc_fft)
  END SUBROUTINE finish_pw_grid

END MODULE m_pw_tofrom_grid
