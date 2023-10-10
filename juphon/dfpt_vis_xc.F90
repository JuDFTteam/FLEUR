!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vis_xc
   USE m_juDFT
   use m_convol
   !     ******************************************************
   !     subroutine generates the exchange-correlation potential
   !     in the interstitial region    c.l.fu
   !     including gradient corrections. t.a. 1996.
   !     ******************************************************
CONTAINS
   SUBROUTINE dfpt_vis_xc(stars,starsq,sym,cell,den,den1,xcpot,input,vTot)

      !     ******************************************************
      !     instead of visxcor.f: the different exchange-correlation
      !     potentials defined through the key icorr are called through
      !     the driver subroutine vxcallg.f,for the energy density - excallg
      !     subroutines vectorized
      !     ** r.pentcheva 22.01.96
      !     *********************************************************
      !     in case of total = .true. calculates the ex-corr. energy
      !     density
      !     ** r.pentcheva 08.05.96
      !     ******************************************************************
      USE m_pw_tofrom_grid
      USE m_types
      USE m_types_xcpot_libxc
      USE m_libxc_postprocess_gga
      USE m_metagga
      IMPLICIT NONE

      CLASS(t_xcpot),INTENT(IN)     :: xcpot
      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_stars),INTENT(IN)      :: stars, starsq
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_potden),INTENT(IN)     :: den, den1
      TYPE(t_potden),INTENT(INOUT)  :: vTot

      TYPE(t_gradients) :: grad
      TYPE(t_potden) :: vTotim

      REAL, ALLOCATABLE :: rho(:,:), rho1re(:,:), rho1im(:,:), ED_rs(:,:), vTot_rs(:,:)
      REAL, ALLOCATABLE :: rho_conv(:,:), ED_conv(:,:), vTot_conv(:,:)
      REAL, ALLOCATABLE :: v_xc1re(:,:),v_xc1im(:,:),f_xc(:,:)
      INTEGER           :: iSpin, jSpin, fxcSpin, i, js, nfxc
      LOGICAL           :: perform_MetaGGA, l_libxc

      nfxc = 2 * input%jspins - 1

      l_libxc=.FALSE.

      IF (ALLOCATED(vTotim%pw)) DEALLOCATE(vTotim%pw)
      ALLOCATE(vTotim%pw,mold=vTot%pw)
      vTotim%pw = CMPLX(0.0,0.0)

      call timestart("init_pw_grid")
      CALL init_pw_grid(stars,sym,cell,xcpot)
      call timestop("init_pw_grid")

      !Put the charge on the grid, in GGA case also calculate gradients
      call timestart("pw_to_grid")
      CALL pw_to_grid(.FALSE.,input%jspins,.FALSE.,stars,cell,den%pw,grad,xcpot,rho)
      CALL pw_to_grid(.FALSE.,input%jspins,.FALSE.,starsq,cell,den1%pw,grad,xcpot,rho1re,rho1im)
      call timestop("pw_to_grid")

      ALLOCATE(f_xc(SIZE(rho,1),nfxc))
      ALLOCATE(v_xc1re,mold=rho)
      ALLOCATE(v_xc1im,mold=rho)

      !call timestart("apply_cutoffs")
      !CALL xcpot%apply_cutoffs(1.E-6,rho,grad)
      !call timestop("apply_cutoffs")
#ifdef CPP_LIBXC
      CALL xcpot%get_fxc(input%jspins, rho, f_xc)
#else
      CALL judft_error("You compiled Fleur without libxc but want to use DFPT. Please fix that.")
      !CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)
      !TODO: Maybe place the old way with x-Alpha here for fun.
#endif

      v_xc1re = 0.0
      v_xc1im = 0.0
      DO iSpin = 1, input%jspins
          DO jSpin = 1, input%jspins
              fxcSpin = iSpin + jSpin - 1
              v_xc1re(:, iSpin) = v_xc1re(:, iSpin) + f_xc(:, fxcSpin) * rho1re(:, jSpin)
              v_xc1im(:, iSpin) = v_xc1im(:, iSpin) + f_xc(:, fxcSpin) * rho1im(:, jSpin)
          END DO
      END DO

      !Put the potentials in rez. space.
      call timestart("pw_from_grid")
      CALL  pw_from_grid(starsq,v_xc1re,vTot%pw)
      CALL  pw_from_grid(starsq,v_xc1im,vTotim%pw)
      vTot%pw = vTot%pw + ImagUnit * vTotim%pw
      call timestop("pw_from_grid")

!      call timestart("finish_pw_grid")
!      CALL finish_pw_grid()
!      call timestop("finish_pw_grid")
   END SUBROUTINE dfpt_vis_xc
END MODULE m_dfpt_vis_xc
