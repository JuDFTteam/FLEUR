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
MODULE m_vis_xc
   USE m_juDFT
   use m_convol
   !     ******************************************************
   !     subroutine generates the exchange-correlation potential
   !     in the interstitial region    c.l.fu
   !     including gradient corrections. t.a. 1996.
   !     ******************************************************
CONTAINS
   SUBROUTINE vis_xc(stars,sym,cell,den,xcpot,input,noco,EnergyDen,vTot,vx,exc)

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

      CLASS(t_xcpot),INTENT(INOUT)  :: xcpot
      TYPE(t_input),INTENT(IN)      :: input
      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_stars),INTENT(IN)      :: stars
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_potden),INTENT(IN)  :: den, EnergyDen
      TYPE(t_potden),INTENT(INOUT)  :: vTot,vx,exc

      TYPE(t_gradients) :: grad, tmp_grad
      REAL, ALLOCATABLE :: rho(:,:), ED_rs(:,:), vTot_rs(:,:)
      REAL, ALLOCATABLE :: rho_conv(:,:), ED_conv(:,:), vTot_conv(:,:)
      REAL, ALLOCATABLE :: v_x(:,:),v_xc(:,:),e_xc(:,:)
      INTEGER           :: jspin, i, js
      LOGICAL           :: perform_MetaGGA

      perform_MetaGGA = ALLOCATED(EnergyDen%mt) &
                      .AND. (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA())
      CALL init_pw_grid(xcpot,stars,sym,cell)

      !Put the charge on the grid, in GGA case also calculate gradients
      CALL pw_to_grid(xcpot,input%jspins,noco%l_noco,stars,cell,den%pw,grad,rho)

      ALLOCATE(v_xc,mold=rho)
      ALLOCATE(v_x,mold=rho)

      if(perform_MetaGGA .and. xcpot%kinED%set) then
         CALL xcpot%get_vxc(input%jspins,rho,v_xc, v_x,grad, kinED_KS=xcpot%kinED%is)
      else
         CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)
      endif

      IF (xcpot%needs_grad()) THEN
         SELECT TYPE(xcpot)
         TYPE IS (t_xcpot_libxc)
            CALL libxc_postprocess_gga_pw(xcpot,stars,cell,v_xc,grad)
         END SELECT
      ENDIF
      !Put the potentials in rez. space.
      CALL  pw_from_grid(xcpot,stars,input%total,v_xc,vTot%pw,vTot%pw_w)
      CALL  pw_from_grid(xcpot,stars,input%total,v_x,vx%pw,vx%pw_w)

      !calculate the ex.-cor energy density
      IF (ALLOCATED(exc%pw_w)) THEN
         ALLOCATE ( e_xc(SIZE(rho,1),1) ); e_xc=0.0

         IF(xcpot%kinED%set) THEN
            CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, xcpot%kinED%is, mt_call=.False.)
            xcpot%lapl%is = grad%laplace
         ELSE
            CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, mt_call=.False.)
         ENDIF
         CALL pw_from_grid(xcpot,stars,.TRUE.,e_xc,exc%pw,exc%pw_w)
      ENDIF

      CALL finish_pw_grid()
   END SUBROUTINE vis_xc
END MODULE m_vis_xc
