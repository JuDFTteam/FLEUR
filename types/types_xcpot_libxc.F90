!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>This module contains the xcpot-type providing an interface to libxc
MODULE m_types_xcpot_libxc
#ifdef CPP_LIBXC  
  USE xc_f03_lib_m
#endif
  USE m_types_xcpot
  USE m_judft
  IMPLICIT NONE
  PRIVATE

  TYPE,EXTENDS(t_xcpot):: t_xcpot_libxc 
#ifdef CPP_LIBXC    
     TYPE(xc_f03_func_t) :: xc_func_x,xc_func_c
     TYPE(xc_f03_func_info_t) :: xc_info_x,xc_info_c
#endif
     INTEGER          :: func_id_c,func_id_x,jspins
   CONTAINS
     PROCEDURE        :: is_gga=>xcpot_is_gga
     PROCEDURE        :: is_hybrid=>xcpot_is_hybrid 
     PROCEDURE        :: get_exchange_weight=>xcpot_get_exchange_weight
     PROCEDURE        :: get_vxc=>xcpot_get_vxc
     PROCEDURE        :: get_exc=>xcpot_get_exc
     PROCEDURE,NOPASS :: alloc_gradients=>xcpot_alloc_gradients
     !Not overloeaded...
     PROCEDURE        :: init=>xcpot_init
  END TYPE t_xcpot_libxc
  PUBLIC t_xcpot_libxc
CONTAINS

  SUBROUTINE xcpot_init(xcpot,jspins,id_x,id_c)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(OUT)    :: xcpot
    INTEGER,INTENT(IN)                  :: jspins,id_x,id_c
    
#ifdef CPP_LIBXC
    INTEGER :: err
    xcpot%jspins=jspins
    xcpot%func_id_x=id_x
    xcpot%func_id_c=id_c
 
    IF (jspins==1) THEN
       CALL xc_f03_func_init(xcpot%xc_func_x, xcpot%func_id_x, XC_UNPOLARIZED)
       IF (xcpot%func_id_c>0) CALL xc_f03_func_init(xcpot%xc_func_c, xcpot%func_id_c, XC_UNPOLARIZED)
    ELSE
       CALL xc_f03_func_init(xcpot%xc_func_x, xcpot%func_id_x, XC_POLARIZED)
       IF (xcpot%func_id_c>0) CALL xc_f03_func_init(xcpot%xc_func_c, xcpot%func_id_c, XC_POLARIZED)
    END IF
    xcpot%xc_info_x=xc_f03_func_get_info(xcpot%xc_func_x)
    IF (xcpot%func_id_c>0) xcpot%xc_info_c=xc_f03_func_get_info(xcpot%xc_func_c)
    
    PRINT *, "TODO: some info and output on libxc functionals"

#else
    CALL judft_error("You specified a libxc-exchange correlation potential but FLEUR is not linked against libxc",hint="Please recompile FLEUR with libxc support")
#endif
  END SUBROUTINE xcpot_init


  LOGICAL FUNCTION xcpot_is_gga(xcpot)
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN):: xcpot
#ifdef CPP_LIBXC    
    xcpot_is_gga=ANY((/XC_FAMILY_GGA, XC_FAMILY_HYB_GGA/)==xc_f03_func_info_get_family(xcpot%xc_info_x))
#endif
  END FUNCTION xcpot_is_gga

  LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN):: xcpot
#ifdef CPP_LIBXC
    xcpot_is_hybrid=ANY((/XC_FAMILY_HYB_MGGA, XC_FAMILY_HYB_GGA/)==xc_f03_func_info_get_family(xcpot%xc_info_x))
#endif
  END FUNCTION xcpot_is_hybrid

  FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN):: xcpot

    REAL:: a_ex
#ifdef CPP_LIBXC    
    a_ex=xc_f03_hyb_exx_coef(xcpot%xc_func_x)
#endif
  END FUNCTION xcpot_get_exchange_weight

  !***********************************************************************
  SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh, vxc,vx, grad)
    !***********************************************************************
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN) :: xcpot
    INTEGER, INTENT (IN)     :: jspins
    REAL,INTENT (IN) :: rh(:,:)
    REAL, INTENT (OUT) :: vx (:,:)
    REAL, INTENT (OUT) :: vxc(:,:)
    ! optional arguments for GGA
    TYPE(t_gradients),OPTIONAL,INTENT(IN)::grad
#ifdef CPP_LIBXC    
    REAL,ALLOCATABLE::vxc_tmp(:,:),vx_tmp(:,:),vsigma(:)
    vx (:,:) = 0.0  !exchange potential not calculated....
    !libxc uses the spin as a first index, hence we have to transpose....
    ALLOCATE(vxc_tmp(SIZE(vxc,2),SIZE(vxc,1)));vxc_tmp=0.0
    ALLOCATE(vx_tmp(SIZE(vx,2),SIZE(vx,1)));vx_tmp=0.0
    IF (xcpot%is_gga()) THEN
       IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_vxc for a GGA potential without providing derivatives")
       ALLOCATE(vsigma(SIZE(grad%sigma)))
       CALL xc_f03_gga_vxc(xcpot%xc_func_x, SIZE(rh,1), TRANSPOSE(rh),grad%sigma,vx_tmp,vsigma)
       IF (xcpot%func_id_c>0) CALL xc_f03_gga_vxc(xcpot%xc_func_c, SIZE(rh,1), TRANSPOSE(rh),grad%sigma,vxc_tmp,vsigma)
    ELSE  !LDA potentials
       CALL xc_f03_lda_vxc(xcpot%xc_func_x, SIZE(rh,1), TRANSPOSE(rh), vxc_tmp)
       IF (xcpot%func_id_c>0) CALL xc_f03_lda_vxc(xcpot%xc_func_c, SIZE(rh,1), TRANSPOSE(rh), vxc_tmp)
    ENDIF
    vx=TRANSPOSE(vx_tmp) 
    vxc_tmp=vxc_tmp+vx_tmp
    vxc=TRANSPOSE(vxc_tmp) 
    
#endif
  END SUBROUTINE xcpot_get_vxc

  !***********************************************************************
  SUBROUTINE xcpot_get_exc(xcpot,jspins,rh,exc,grad)
    !***********************************************************************
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN) :: xcpot
    INTEGER, INTENT (IN)     :: jspins
    REAL,INTENT (IN) :: rh(:,:)
    REAL, INTENT (OUT) :: exc(:)
    ! optional arguments for GGA      
    TYPE(t_gradients),OPTIONAL,INTENT(IN)::grad
#ifdef CPP_LIBXC    
    IF (xcpot%is_gga()) THEN
       IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_vxc for a GGA potential without providing derivatives")
       CALL xc_f03_gga_exc(xcpot%xc_func_x, SIZE(rh,1), TRANSPOSE(rh),grad%sigma,exc)
       IF (xcpot%func_id_c>0) CALL xc_f03_gga_exc(xcpot%xc_func_c, SIZE(rh,1), TRANSPOSE(rh),grad%sigma,exc)
    ELSE  !LDA potentials
       CALL xc_f03_lda_exc(xcpot%xc_func_x, SIZE(rh,1), TRANSPOSE(rh), exc)
       IF (xcpot%func_id_c>0) CALL xc_f03_lda_exc(xcpot%xc_func_c, SIZE(rh,1), TRANSPOSE(rh), exc)
    ENDIF
#endif
  END SUBROUTINE xcpot_get_exc

  SUBROUTINE xcpot_alloc_gradients(ngrid,jspins,grad)
    INTEGER, INTENT (IN)         :: jspins,ngrid
    TYPE(t_gradients),INTENT(OUT):: grad
    !For libxc we only need the sigma array...
    ALLOCATE(grad%sigma(MERGE(1,3,jspins==1),ngrid))
  END SUBROUTINE xcpot_alloc_gradients


END MODULE m_types_xcpot_libxc
