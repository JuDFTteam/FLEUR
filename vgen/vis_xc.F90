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
   SUBROUTINE vis_xc(stars,sym,cell,den,xcpot,input,noco,EnergyDen,kinED,vTot,vx,exc,vxc)

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
      TYPE(t_noco),INTENT(IN)       :: noco
      TYPE(t_sym),INTENT(IN)        :: sym
      TYPE(t_stars),INTENT(IN)      :: stars
      TYPE(t_cell),INTENT(IN)       :: cell
      TYPE(t_potden),INTENT(IN)  :: den, EnergyDen
      TYPE(t_potden),INTENT(INOUT)  :: vTot,vx,exc,vxc
      TYPE(t_kinED),INTENT(IN)      ::kinED

      TYPE(t_gradients) :: grad
      REAL, ALLOCATABLE :: rho(:,:), ED_rs(:,:), vTot_rs(:,:)
      REAL, ALLOCATABLE :: rho_conv(:,:), ED_conv(:,:), vTot_conv(:,:)
      REAL, ALLOCATABLE :: v_x(:,:),v_xc(:,:),v_xc2(:,:),e_xc(:,:)
      INTEGER           :: jspin, i, js
      LOGICAL           :: perform_MetaGGA, l_libxc

      l_libxc=.FALSE.

      perform_MetaGGA = ALLOCATED(EnergyDen%mt) &
                      .AND. (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA())

      call timestart("init_pw_grid")
      CALL init_pw_grid(xcpot%needs_grad(),stars,sym,cell)
      call timestop("init_pw_grid")

      !Put the charge on the grid, in GGA case also calculate gradients
      call timestart("pw_to_grid")
      CALL pw_to_grid(xcpot%needs_grad(),input%jspins,noco%l_noco,stars,cell,den%pw,grad,xcpot,rho)
      call timestop("pw_to_grid")

      ALLOCATE(v_xc,mold=rho)
      ALLOCATE(v_xc2,mold=rho)
      ALLOCATE(v_x,mold=rho)

      call timestart("apply_cutoffs")
      CALL xcpot%apply_cutoffs(1.E-6,rho,grad)
      call timestop("apply_cutoffs")
#ifdef CPP_LIBXC
      if(perform_MetaGGA .and. kinED%set) then
         CALL xcpot%get_vxc(input%jspins,rho,v_xc, v_x,grad, kinEnergyDen_KS=kinED%is)
      else
         CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)
      endif
#else
      call timestart("get_vxc")   
      CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)
      call timestop("get_vxc")
#endif
      
      SELECT TYPE(xcpot)
      TYPE IS (t_xcpot_libxc)
         l_libxc=.TRUE.
         IF (xcpot%needs_grad()) THEN
            CALL libxc_postprocess_gga_pw(xcpot,stars,cell,v_xc,grad)
            CALL libxc_postprocess_gga_pw(xcpot,stars,cell,v_x,grad)
         END IF
      END SELECT

      !IF (l_libxc.AND.xcpot%needs_grad()) THEN
      !   CALL save_npy('vxc_gga_ir_libxc.npy',v_xc)
      !ELSE IF (l_libxc.AND.(.NOT.xcpot%needs_grad())) THEN
      !  CALL save_npy('vxc_lda_ir_libxc.npy',v_xc)
      !ELSE IF ((.NOT.l_libxc).AND.xcpot%needs_grad()) THEN
      !   CALL save_npy('vxc_gga_ir_inbuild.npy',v_xc)
      !ELSE
      !  CALL save_npy('vxc_lda_ir_inbuild.npy',v_xc)
      !END IF

      v_xc2=v_xc
      !Put the potentials in rez. space.
      call timestart("pw_from_grid")
      CALL  pw_from_grid(stars,v_xc,vTot%pw,vTot%pw_w)
      CALL  pw_from_grid(stars,v_xc2,vxc%pw)
      CALL  pw_from_grid(stars,v_x,vx%pw,vx%pw_w)
      call timestop("pw_from_grid")

      !calculate the ex.-cor energy density
      IF (ALLOCATED(exc%pw_w)) THEN
         ALLOCATE ( e_xc(SIZE(rho,1),1) ); e_xc=0.0
#ifdef CPP_LIBXC
         IF(kinED%set) THEN
            CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, kinED%is, mt_call=.False.)
         ELSE
            CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, mt_call=.False.)
         ENDIF

#else
         call timestart("get_exc")  
         CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad, mt_call=.False.)
         call timestop("get_exc")
#endif
         call timestart("pw_from_grid")
         CALL pw_from_grid(stars,e_xc,exc%pw,exc%pw_w)
         call timestop("pw_from_grid")
      ENDIF

      call timestart("finish_pw_grid")
      CALL finish_pw_grid()
      call timestop("finish_pw_grid")
   END SUBROUTINE vis_xc
END MODULE m_vis_xc
