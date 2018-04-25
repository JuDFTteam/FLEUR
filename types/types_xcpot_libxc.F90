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
     TYPE(xc_f03_func_t) :: xc_func
     TYPE(xc_f03_func_info_t) :: xc_info
#endif
     INTEGER          :: func_id
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

  SUBROUTINE xcpot_init(xcpot,name)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(OUT)    :: xcpot
    CHARACTER(len=*),INTENT(IN)         :: name

#ifdef CPP_LIBXC
    INTEGER :: err

    READ(name(7:),*,stat=err) xcpot%func_id
    IF (err.NE.0) CALL judft_error("Error in specifying libxc xc-functional: '"//name//"' is not a valid specification")
    IF (name(1:6).NE."libxc:") CALL judft_error("BUG:inconsistency in xcpot%init")

    IF (jspins==1) THEN
       CALL xc_f03_func_init(xcpot%xc_func, xcpot%func_id, XC_UNPOLARIZED)
    ELSE
       CALL xc_f03_func_init(xcpot%xc_func, xcpot%func_id, XC_POLARIZED)
    END IF
    xcpot%xc_info=xc_f03_func_get_info(xcpot%xc_func)

    PRINT * "TODO: some info and output on libxc functionals"

#else
    CALL judft_error("You specified a libxc-exchange correlation potential but FLEUR is not linked against libxc",hint="Please recompile FLEUR with libxc support")
#endif
  END SUBROUTINE xcpot_init


  LOGICAL FUNCTION xcpot_is_gga(xcpot)
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN):: xcpot
#ifdef CPP_LIBXC    
    xcpot_is_gga=ANY((/XC_FAMILY_GGA, XC_FAMILY_HYB_GGA/)==xc_f03_func_info_get_family(xcpot%xc_info))
#endif
  END FUNCTION xcpot_is_gga

  LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN):: xcpot
#ifdef CPP_LIBXC
    xcpot_is_hybrid=ANY((/XC_FAMILY_HYB_MGGA, XC_FAMILY_HYB_GGA/)==xc_f03_func_info_get_family(xcpot%xc_info))
#endif
  END FUNCTION xcpot_is_hybrid

  FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot_libxc),INTENT(IN):: xcpot

    REAL:: a_ex
#ifdef CPP_LIBXC    
    a_ex=xc_f03_hyb_exx_coef(xcpot%xc_func)
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
    REAL,ALLOCATABLE::vxc_temp(:,:)
    vx (:,:) = 0.0  !exchange potential not calculated....
    !libxc uses the spin as a first index, hence we have to transpose....
    ALLOCATE(vxc_tmp(SIZE(vxc,2),SIZE(vxc,1)))

    IF (xcpot%is_gga()) THEN
       IF (.NOT.PRESENT(grad)) CALL judft_error("Bug: You called get_vxc for a GGA potential without providing derivatives")
       CALL xc_f03_gga_vxc(xc_func, SIZE(rh,1), TRANSPOSE(rh),grad%sigma,vxc_tmp)
    ELSE  !LDA potentials
       CALL xc_f03_lda_vxc(xc_func, SIZE(rh,1), TRANSPOSE(rh), vxc_tmp)
    ENDIF
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
       CALL xc_f03_gga_exc(xc_func, SIZE(rh,1), TRANSPOSE(rh),grad%sigma,exc)
    ELSE  !LDA potentials
       CALL xc_f03_lda_exc(xc_func, SIZE(rh,1), TRANSPOSE(rh), exc)
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
