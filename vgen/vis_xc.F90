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
  !     ******************************************************
  !     subroutine generates the exchange-correlation potential
  !     in the interstitial region    c.l.fu
  !     including gradient corrections. t.a. 1996.
  !     ******************************************************
CONTAINS
  SUBROUTINE vis_xc(stars,sym,cell,den,xcpot,input,noco, vxc,vx,exc)

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
    IMPLICIT NONE

    CLASS(t_xcpot),INTENT(IN)     :: xcpot
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_stars),INTENT(IN)      :: stars
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_potden),INTENT(IN)     :: den
    TYPE(t_potden),INTENT(INOUT)  :: vxc,vx,exc
    
    TYPE(t_gradients) :: grad
    REAL, ALLOCATABLE :: rho(:,:)
    REAL, ALLOCATABLE :: v_x(:,:),v_xc(:,:),e_xc(:,:) 
  
    CALL init_pw_grid(xcpot,stars,sym,cell)

    !Put the charge on the grid, in GGA case also calculate gradients
    CALL pw_to_grid(xcpot,input,noco,stars,cell,den,rho,grad)
    ALLOCATE(v_xc,v_x,mold=rho)

    CALL xcpot%get_vxc(input%jspins,rho,v_xc,v_x,grad)

    !Put the potentials in rez. space. 
    CALL  pw_from_grid(xcpot,stars,input%total,v_xc,vxc)
    CALL  pw_from_grid(xcpot,stars,input%total,v_x,vx)

    !calculate the ex.-cor energy density 
    IF (ALLOCATED(exc%pw_w)) THEN
       ALLOCATE ( e_xc(SIZE(rho,1),1) );e_xc=0.0
       CALL xcpot%get_exc(input%jspins,rho,e_xc(:,1),grad)
       CALL pw_from_grid(xcpot,stars,.TRUE.,e_xc,exc)
    ENDIF

    CALL finish_pw_grid()
    
  END SUBROUTINE vis_xc
END MODULE m_vis_xc
