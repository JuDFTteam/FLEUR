!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
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
   use m_types_misc
   IMPLICIT NONE

#ifdef CPP_LIBXC
   PRIVATE :: write_xc_info
#endif

   TYPE,EXTENDS(t_xcpot):: t_xcpot_libxc
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_t)      :: vxc_func_x, vxc_func_c
      TYPE(xc_f03_func_t)      :: exc_func_x, exc_func_c
      TYPE(xc_f03_func_t)      :: aux_func_x, aux_func_c
      LOGICAL                  :: l_has_aux = .FALSE.
#endif
      INTEGER                  :: jspins

   CONTAINS
      !PROCEDURE        :: vxc_is_LDA          => xcpot_vxc_is_LDA
      !PROCEDURE        :: vxc_is_gga          => xcpot_vxc_is_gga

      PROCEDURE        :: vx_is_LDA => xcpot_vx_is_LDA
      PROCEDURE        :: vx_is_GGA => xcpot_vx_is_GGA
      PROCEDURE        :: vx_is_MetaGGA => xcpot_vx_is_MetaGGA

      PROCEDURE        :: vc_is_LDA => xcpot_vc_is_LDA
      PROCEDURE        :: vc_is_GGA => xcpot_vc_is_GGA

      PROCEDURE        :: exc_is_LDA => xcpot_exc_is_LDA
      PROCEDURE        :: exc_is_gga => xcpot_exc_is_gga
      PROCEDURE        :: exc_is_MetaGGA => xcpot_exc_is_MetaGGA

      PROCEDURE        :: is_hybrid => xcpot_is_hybrid
      PROCEDURE        :: get_exchange_weight => xcpot_get_exchange_weight
      PROCEDURE        :: get_vxc => xcpot_get_vxc
      PROCEDURE        :: get_exc => xcpot_get_exc
      PROCEDURE        :: get_fxc => xcpot_get_fxc
      PROCEDURE, NOPASS :: alloc_gradients => xcpot_alloc_gradients
      !Not             overloeaded...
      PROCEDURE        :: init => xcpot_init
      PROCEDURE        :: create_from_aux => xcpot_create_from_aux
      PROCEDURE,NOPASS :: apply_cutoffs
   END TYPE t_xcpot_libxc
   PUBLIC t_xcpot_libxc
CONTAINS
  subroutine apply_cutoffs(density_cutoff,rh,grad)
    real,intent(INOUT) :: rh(:,:)
    real,INTENT(IN)    :: density_cutoff
    type(t_gradients),INTENT(INOUT),OPTIONAL :: grad



    integer:: i,j
    DO j=1,size(rh,2)
      DO i=1,size(rh,1)
        if (abs(rh(i,j))<density_cutoff) THEN
          rh(i,j)=density_cutoff
          if (present(grad)) Then
            if (allocated(grad%sigma)) grad%sigma(:,i)=0.0 !if one spin is small, apply cutoff to all gradients!
            if (allocated(grad%gr)) grad%gr(:,i,j)=0.0
            if (allocated(grad%laplace)) grad%laplace(i,j)=0.0
          endif
        endif
        if (allocated(grad%sigma)) then
         where (abs(grad%sigma(:,i))<1E-5) grad%sigma(:,i)=0.0
        endif 
      ENDDO
    ENDDO

  end subroutine

   SUBROUTINE xcpot_init(xcpot, func_vxc_id_x, func_vxc_id_c, func_exc_id_x, func_exc_id_c, jspins, l_bj)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(INOUT)    :: xcpot
      INTEGER, INTENT(IN)                 :: jspins, func_vxc_id_x, func_vxc_id_c, func_exc_id_x, func_exc_id_c
      LOGICAL, INTENT(IN)                 :: l_bj
      LOGICAL                             :: same_functionals   ! are vxc and exc equal
      INTEGER                             :: errors(4)

#ifdef CPP_LIBXC
      errors = -1
      xcpot%jspins = jspins
      xcpot%l_bj = l_bj
      xcpot%func_vxc_id_x = func_vxc_id_x
      xcpot%func_exc_id_x = func_exc_id_x
      xcpot%func_vxc_id_c = func_vxc_id_c
      xcpot%func_exc_id_c = func_exc_id_c

      IF (xcpot%func_vxc_id_x == 0 .OR. xcpot%func_exc_id_x == 0) THEN
         CALL judft_error("LibXC exchange- and correlation-function indicies need to be set" &
                          , hint='Try this: '//ACHAR(10)// &
                          '<xcFunctional name="libxc" relativisticCorrections="F">'//ACHAR(10)// &
                          '  <libXC  exchange="1" correlation="1" /> '//ACHAR(10)// &
                          '</xcFunctional> ')
      ENDIF

      IF (jspins==1) THEN
         ! potential functionals
         CALL xc_f03_func_init(xcpot%vxc_func_x, xcpot%func_vxc_id_x, XC_UNPOLARIZED, err=errors(1))
         IF (xcpot%func_vxc_id_c>0) CALL xc_f03_func_init(xcpot%vxc_func_c, xcpot%func_vxc_id_c, &
                                                                 XC_UNPOLARIZED, err=errors(2))

         ! energy functionals
         CALL xc_f03_func_init(xcpot%exc_func_x, xcpot%func_exc_id_x, XC_UNPOLARIZED, err=errors(3))
         IF (xcpot%func_exc_id_c>0) CALL xc_f03_func_init(xcpot%exc_func_c, xcpot%func_exc_id_c, &
                                                                  XC_UNPOLARIZED, err=errors(4))

      ELSE
         ! potential functionals
         CALL xc_f03_func_init(xcpot%vxc_func_x, xcpot%func_vxc_id_x, XC_POLARIZED, err=errors(1))
         IF (xcpot%func_vxc_id_c>0) CALL xc_f03_func_init(xcpot%vxc_func_c, xcpot%func_vxc_id_c, &
                                                                  XC_POLARIZED, err=errors(2))

         !energy functionals
         CALL xc_f03_func_init(xcpot%exc_func_x, xcpot%func_exc_id_x, XC_POLARIZED, err=errors(3))
         IF (xcpot%func_exc_id_c>0) CALL xc_f03_func_init(xcpot%exc_func_c, xcpot%func_exc_id_c, &
                                                                  XC_POLARIZED, err=errors(4))
      END IF

      !IF(errors(1) /= 0) call juDFT_error("Exchange potential functional not in LibXC")
      !IF(errors(2) /= 0) call juDFT_error("Correlation potential functional not in LibXC")
      !IF(errors(3) /= 0) call juDFT_error("Exchange energy functional not in LibXC")
      !IF(errors(4) /= 0) call juDFT_error("Correlation energy functional not in LibXC")

      !check if any potental is a MetaGGA -- now allowed for self-consistent MetaGGA
      IF (ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA] == xc_get_family(xcpot%vxc_func_x))) THEN
         WRITE(*,*) "MetaGGA potential functional detected - V_tau will be computed"
      ELSEIF (xcpot%func_vxc_id_c > 0) THEN
         IF (ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA] == xc_get_family(xcpot%vxc_func_c))) THEN
            WRITE(*,*) "MetaGGA correlation potential functional detected - V_tau will be computed"
         ENDIF
      ENDIF

      ! Initialize auxiliary GGA functional for radial basis generation in MetaGGA
      xcpot%l_has_aux = (xcpot%func_aux_id_x > 0)
      IF (xcpot%l_has_aux) THEN
         WRITE(*,*) "Auxiliary Potential for radial basis: exchange=", xcpot%func_aux_id_x, &
                    " correlation=", xcpot%func_aux_id_c
         IF (jspins == 1) THEN
            CALL xc_f03_func_init(xcpot%aux_func_x, xcpot%func_aux_id_x, XC_UNPOLARIZED)
            IF (xcpot%func_aux_id_c > 0) &
               CALL xc_f03_func_init(xcpot%aux_func_c, xcpot%func_aux_id_c, XC_UNPOLARIZED)
         ELSE
            CALL xc_f03_func_init(xcpot%aux_func_x, xcpot%func_aux_id_x, XC_POLARIZED)
            IF (xcpot%func_aux_id_c > 0) &
               CALL xc_f03_func_init(xcpot%aux_func_c, xcpot%func_aux_id_c, XC_POLARIZED)
         ENDIF
      ENDIF

      CALL write_xc_info(xcpot%vxc_func_x)

      IF (xcpot%func_vxc_id_c > 0) THEN
         CALL write_xc_info(xcpot%vxc_func_c)
      ELSE
         WRITE (*, *) "No Correlation functional"
      END IF

      same_functionals = (xcpot%func_vxc_id_x == xcpot%func_exc_id_x) &
                         .AND. (xcpot%func_vxc_id_c == xcpot%func_exc_id_c)
      IF (.NOT. same_functionals) THEN
         CALL write_xc_info(xcpot%exc_func_x)
         IF (xcpot%func_exc_id_c > 0) THEN
            CALL write_xc_info(xcpot%exc_func_c)
         ELSE
            WRITE (*, *) "No Correlation functional for TotalE"
         ENDIF
      ELSE
         WRITE (*, *) "Using same functional for VXC and EXC"
      END IF
#else
      CALL judft_error("You specified a libxc-exchange correlation potential but FLEUR is not linked against libxc", &
                       hint="Please recompile FLEUR with libxc support")
#endif
   END SUBROUTINE xcpot_init

   ! LDA
   LOGICAL FUNCTION xcpot_vx_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info
      if (xcpot%func_vxc_id_x == 0) then
         xcpot_vx_is_LDA = .false.
         return
      endif   
      xc_info = xc_f03_func_get_info(xcpot%vxc_func_x)
      xcpot_vx_is_LDA =  XC_FAMILY_LDA == xc_f03_func_info_get_family(xc_info)
#else
      xcpot_vx_is_LDA = .false.
#endif
   END FUNCTION xcpot_vx_is_LDA

   LOGICAL FUNCTION xcpot_vc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info

      if (xcpot%func_vxc_id_c == 0) then
         xcpot_vc_is_LDA = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%vxc_func_c)
      xcpot_vc_is_LDA =  XC_FAMILY_LDA == xc_f03_func_info_get_family(xc_info)
#else
      xcpot_vc_is_LDA = .false.
#endif
   END FUNCTION xcpot_vc_is_LDA

   LOGICAL FUNCTION xcpot_exc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info

      if (xcpot%func_exc_id_x == 0) then
         xcpot_exc_is_LDA = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%exc_func_x)
      xcpot_exc_is_LDA = (XC_FAMILY_LDA == xc_f03_func_info_get_family(xc_info))
#else
      xcpot_exc_is_LDA = .false.
#endif
   END FUNCTION xcpot_exc_is_LDA

   ! GGA
   LOGICAL FUNCTION xcpot_vc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info

      if (xcpot%func_vxc_id_c == 0) then
         xcpot_vc_is_gga = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%vxc_func_c)
      xcpot_vc_is_gga =  ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA]==xc_f03_func_info_get_family(xc_info))
#else
      xcpot_vc_is_gga = .false.
#endif
   END FUNCTION xcpot_vc_is_gga

   LOGICAL FUNCTION xcpot_vx_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info

      if (xcpot%func_vxc_id_x == 0) then
         xcpot_vx_is_gga = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%vxc_func_x)
      xcpot_vx_is_gga =  ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA]==xc_f03_func_info_get_family(xc_info))
#else
      xcpot_vx_is_gga = .false.
#endif
   END FUNCTION xcpot_vx_is_gga

   LOGICAL FUNCTION xcpot_vx_is_MetaGGA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info
      if (xcpot%func_vxc_id_x == 0) then
         xcpot_vx_is_MetaGGA = .false.
         return
      endif
      if (xcpot%l_bj) then
         xcpot_vx_is_MetaGGA = .TRUE.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%vxc_func_x)
      xcpot_vx_is_MetaGGA =  ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA]==xc_f03_func_info_get_family(xc_info))
#else
      xcpot_vx_is_MetaGGA = .false.
#endif
   END FUNCTION xcpot_vx_is_MetaGGA

   LOGICAL FUNCTION xcpot_exc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info
      if (xcpot%func_exc_id_x == 0) then
         xcpot_exc_is_gga = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%exc_func_x)
      xcpot_exc_is_gga =  ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA]==xc_f03_func_info_get_family(xc_info))
#else
      xcpot_exc_is_gga = .false.
#endif
   END FUNCTION xcpot_exc_is_gga

   LOGICAL FUNCTION xcpot_exc_is_MetaGGA(xcpot)
      IMPLICIT NONE
   CLASS(t_xcpot_libxc),INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info
      if (xcpot%func_exc_id_x == 0) then
         xcpot_exc_is_MetaGGA = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%exc_func_x)
      xcpot_exc_is_MetaGGA=ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA]==xc_f03_func_info_get_family(xc_info))
#else
      xcpot_exc_is_MetaGGA = .False.
#endif
   END FUNCTION xcpot_exc_is_MetaGGA

   LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)        :: xc_info
      if (xcpot%func_vxc_id_x == 0) then
         xcpot_is_hybrid = .false.
         return
      endif
      xc_info = xc_f03_func_get_info(xcpot%vxc_func_x)
      xcpot_is_hybrid=ANY([XC_FAMILY_HYB_MGGA, XC_FAMILY_HYB_GGA]==xc_f03_func_info_get_family(xc_info))
#else
      xcpot_is_hybrid = .False.
#endif
   END FUNCTION xcpot_is_hybrid

   FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot

      REAL:: a_ex
#ifdef CPP_LIBXC
      a_ex=xc_f03_hyb_exx_coef(xcpot%vxc_func_x)
#endif
   END FUNCTION xcpot_get_exchange_weight

   !***********************************************************************
   SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh, vxc,vx, grad, kinenergyden_ks, vtau, l_aux)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      use iso_c_binding
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN) :: xcpot
      INTEGER, INTENT(IN)     :: jspins
      REAL, INTENT(IN)         :: rh(:, :)   !Dimensions here
      REAL, INTENT(OUT)       :: vx(:, :)  !points,spin
      REAL, INTENT(OUT)     :: vxc(:, :)  !
      ! optional arguments for GGA
      TYPE(t_gradients),OPTIONAL,INTENT(INOUT)::grad
      REAL, INTENT(IN), OPTIONAL     :: kinenergyden_ks(:, :)
      ! optional output for MetaGGA: V_tau = dE_xc/d(tau)
      REAL, INTENT(OUT), OPTIONAL    :: vtau(:, :)
      ! optional: use auxiliary GGA functional instead of primary functional
      LOGICAL, INTENT(IN), OPTIONAL  :: l_aux
#ifdef CPP_LIBXC
   REAL, ALLOCATABLE :: vxc_tmp(:,:), vx_tmp(:,:)
   REAL, ALLOCATABLE :: kinED_libxc(:,:), vtau_tmp(:,:)
   REAL, ALLOCATABLE :: sigma(:,:), vsigma(:,:), laplace(:,:)
   TYPE(xc_f03_func_t) :: func_x, func_c
   LOGICAL :: use_aux, has_c

      use_aux = .FALSE.
      IF (PRESENT(l_aux)) use_aux = l_aux
      IF (.NOT. use_aux) use_aux = xcpot%vx_is_MetaGGA() .AND. (.NOT. PRESENT(kinenergyden_ks)) .AND. xcpot%l_has_aux

      ! libxc uses spin as first index, hence transpose on output
      ALLOCATE(vxc_tmp(SIZE(vxc, 2), SIZE(vxc, 1))); vxc_tmp = 0.0
      ALLOCATE(vx_tmp(SIZE(vx, 2), SIZE(vx, 1))); vx_tmp = 0.0
      ALLOCATE(vtau_tmp(SIZE(vx, 2), SIZE(vx, 1))); vtau_tmp = 0.0

      IF (PRESENT(vtau)) vtau = 0.0

      ! Select exchange and correlation functionals
      IF (use_aux) THEN
         IF (.NOT. xcpot%l_has_aux) &
            CALL judft_error("get_vxc with l_aux=.TRUE. but no auxiliary GGA configured")
         IF (.NOT. PRESENT(grad)) &
            CALL judft_error("get_vxc with l_aux=.TRUE. requires gradients")
         func_x = xcpot%aux_func_x
         has_c = xcpot%func_aux_id_c > 0
         IF (has_c) func_c = xcpot%aux_func_c
      ELSE
         func_x = xcpot%vxc_func_x
         has_c = xcpot%func_vxc_id_c > 0
         IF (has_c) func_c = xcpot%vxc_func_c
      ENDIF

      ! Prepare gradient arrays (leave unallocated for LDA)
      IF (PRESENT(grad)) THEN
         sigma = grad%sigma
         ALLOCATE(vsigma, mold=grad%vsigma); vsigma = 0.0
         IF (ALLOCATED(grad%laplace)) laplace = grad%laplace
      ENDIF

      ! Prepare kinetic energy density for MetaGGA
      IF (PRESENT(kinenergyden_ks)) kinED_libxc = transpose(kinenergyden_ks)

      ! Evaluate exchange functional (auto-detects LDA/GGA/MetaGGA)
      CALL eval_vxc(func_x, vx_tmp, sigma, vsigma, laplace, kinED_libxc, vtau_tmp)

      ! Accumulate correlation on top of exchange
      vxc_tmp = vx_tmp
      IF (has_c) CALL eval_vxc(func_c, vxc_tmp, sigma, vsigma, laplace, kinED_libxc, vtau_tmp)
      ! Apply BJ correction if needed
      IF (xcpot%l_bj .AND. PRESENT(kinenergyden_ks)) CALL eval_BJ_correction(vx_tmp, vxc_tmp)
      ! Copy back gradient results
      IF (PRESENT(grad) .AND. ALLOCATED(vsigma)) grad%vsigma = vsigma
      vx = TRANSPOSE(vx_tmp)
      vxc = TRANSPOSE(vxc_tmp)
      IF (PRESENT(vtau)) vtau = TRANSPOSE(vtau_tmp)

   CONTAINS

      !> Unified evaluation of a single xc functional (exchange OR correlation).
      !! Auto-detects LDA/GGA/MetaGGA family via libxc and calls the appropriate
      !! xc_f03_*_vxc routine. Results are accumulated (added) into the output arrays.
      !! Uses rh from host association (parent scope).
      SUBROUTINE eval_vxc(func, vxc_out, sigma, vsigma_out, laplace, kinED, vtau_out)
         TYPE(xc_f03_func_t), INTENT(IN)      :: func
         REAL, INTENT(INOUT)                   :: vxc_out(:, :)    ! (spin, points) - accumulated
         REAL, INTENT(IN), ALLOCATABLE         :: sigma(:, :)      ! (nsigma, npoints)
         REAL, INTENT(INOUT), ALLOCATABLE      :: vsigma_out(:, :) ! (nsigma, npoints) - accumulated
         REAL, INTENT(IN), ALLOCATABLE         :: laplace(:, :)    ! (npoints, jspin)
         REAL, INTENT(IN), ALLOCATABLE         :: kinED(:, :)      ! kinetic energy density
         REAL, INTENT(INOUT), ALLOCATABLE      :: vtau_out(:, :)   ! (spin, npoints) - accumulated

         REAL, ALLOCATABLE               :: vxc(:, :), vsigma(:, :), vtau(:, :), vlapl(:, :)
         TYPE(xc_f03_func_info_t)        :: info
         INTEGER                         :: family

         ALLOCATE(vxc, mold=vxc_out); vxc = 0.0
         info = xc_f03_func_get_info(func)
         family = xc_f03_func_info_get_family(info)

         IF (ANY([XC_FAMILY_LDA, XC_FAMILY_HYB_LDA] == family)) THEN
            CALL xc_f03_lda_vxc(func, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), vxc)

         ELSEIF (ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA] == family)) THEN
            IF (.NOT. ALLOCATED(sigma)) CALL judft_error("eval_vxc: GGA functional requires sigma (gradients)")
            ALLOCATE(vsigma, mold=vsigma_out); vsigma = 0.0
            CALL xc_f03_gga_vxc(func, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), sigma, vxc, vsigma)
            vsigma_out = vsigma_out + vsigma

         ELSEIF (ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA] == family)) THEN
            IF (.NOT. ALLOCATED(sigma)) CALL judft_error("eval_vxc: MetaGGA functional requires sigma (gradients)")
            IF (.NOT. ALLOCATED(kinED)) CALL judft_error("eval_vxc: MetaGGA functional requires kinetic energy density")
            ALLOCATE(vsigma, mold=vsigma_out); vsigma = 0.0
            ALLOCATE(vtau, mold=vtau_out); vtau = 0.0
            ALLOCATE(vlapl(SIZE(rh, 2), SIZE(rh, 1))); vlapl = 0.0
            CALL xc_f03_mgga_vxc(func, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), sigma, &
                                 TRANSPOSE(laplace), kinED, vxc, vsigma, vlapl, vtau)
            vsigma_out = vsigma_out + vsigma
            vtau_out = vtau_out + vtau

         ELSE
            CALL judft_error("Functional family not supported in eval_vxc")
         ENDIF

         vxc_out = vxc_out + vxc
      END SUBROUTINE eval_vxc

      !> Becke-Johnson correction term, Eq. (2) of
      !! Tran, Blaha, Schwarz, J. Phys.: Condens. Matter 19, 196208 (2007)
      !!
      !! Adds (1/π)√(5/12) √(2τ_σ(r)/ρ_σ(r)) to already-evaluated vx and vxc.
      !! Uses rh, jspins, kinenergyden_ks from host association (parent scope).
      SUBROUTINE eval_BJ_correction(vx_out, vxc_out)
         USE m_constants
         REAL, INTENT(INOUT)    :: vx_out(:, :), vxc_out(:, :) ! (spin, points)

         REAL, PARAMETER :: BJ_prefactor = (1.0/pi_const) * SQRT(5.0/12.0)
         REAL    :: rho_sigma, tau_sigma, bj_corr
         INTEGER :: i, jspin

         DO jspin = 1, jspins
            DO i = 1, SIZE(rh, 1)
               rho_sigma = rh(i, jspin)
               tau_sigma = kinenergyden_ks(i, jspin)
               IF (rho_sigma > 1.0e-10 .AND. tau_sigma > 0.0) THEN
                  bj_corr = BJ_prefactor * SQRT(2.0 * tau_sigma / rho_sigma)
                  vx_out(jspin, i) = vx_out(jspin, i) + bj_corr
                  vxc_out(jspin, i) = vxc_out(jspin, i) + bj_corr
               ENDIF
            ENDDO
         ENDDO
      END SUBROUTINE eval_BJ_correction

#endif
   END SUBROUTINE xcpot_get_vxc

   SUBROUTINE xcpot_get_exc(xcpot, jspins, rh, exc, grad, kinEnergyDen_KS, mt_call)
      use m_constants
      use ISO_C_BINDING
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN)          :: xcpot
      INTEGER, INTENT(IN)                  :: jspins
      REAL, INTENT(IN)                      :: rh(:, :)  !points,spin
      REAL, INTENT(OUT)                    :: exc(:) !points
      ! optional arguments for GGA
      TYPE(t_gradients), OPTIONAL, INTENT(IN) :: grad
      LOGICAL, OPTIONAL, INTENT(IN)         :: mt_call

      ! kinED from Kohn-Sham equations:
      ! tau = sum[phi_i(r)^dag nabla phi_i(r)]
      ! see eq (2) in https://doi.org/10.1063/1.1565316
      ! (-0.5 is applied below)
      REAL, INTENT(IN), OPTIONAL     :: kinEnergyDen_KS(:, :)

#ifdef CPP_LIBXC
      TYPE(xc_f03_func_info_t)       :: xc_info
      REAL  :: excc(SIZE(exc))
      REAL  :: cut_ratio = 0.1
      INTEGER :: cut_idx
      LOGICAL :: is_mt

      ! tau = 0.5 * sum[|grad phi_i(r)|²]
      ! see eq (3) in https://doi.org/10.1063/1.1565316
      REAL, ALLOCATABLE              :: kinEnergyDen_libXC(:, :), pkzb_ratio(:, :), pkzb_zaehler(:, :), pkzb_nenner(:, :)

      is_mt = merge(mt_call, .False., present(mt_call))
      IF (xcpot%exc_is_gga()) THEN
         IF (.NOT. PRESENT(grad)) CALL judft_error("Bug: You called get_exc for a GGA potential without providing derivatives")
         CALL xc_f03_gga_exc(xcpot%exc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, exc)
         IF (xcpot%func_exc_id_c > 0) THEN
            CALL xc_f03_gga_exc(xcpot%exc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, excc)
            exc = exc + excc
         END IF
      ELSEIF (xcpot%exc_is_LDA()) THEN  !LDA potentials
         CALL xc_f03_lda_exc(xcpot%exc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), exc)
         IF (xcpot%func_exc_id_c > 0) THEN
            CALL xc_f03_lda_exc(xcpot%exc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), excc)
            exc = exc + excc
         END IF
      ELSEIF (xcpot%exc_is_MetaGGA()) THEN
         IF (PRESENT(kinEnergyDen_KS)) THEN
            
            kinEnergyDen_libXC = transpose(kinEnergyDen_KS )

            !only cut core of muffin tin
            cut_idx = MERGE(NINT(size(rh, 1)*cut_ratio), 0, is_mt)

            exc = 0.0
            excc = 0.0
            call xc_f03_mgga_exc(xcpot%exc_func_x, SIZE(rh(cut_idx + 1:, :), 1, kind=c_size_t), &
                                 TRANSPOSE(rh(cut_idx + 1:, :)), &
                                 grad%sigma(:, cut_idx + 1:), &
                                 transpose(grad%laplace(cut_idx + 1:, :)), &
                                 kinEnergyDen_libXC(:, cut_idx + 1:), &
                                 exc(cut_idx + 1:))

            call xc_f03_gga_exc(xcpot%vxc_func_x, SIZE(rh(:cut_idx, :), 1, kind=c_size_t), &
                                TRANSPOSE(rh(:cut_idx, :)), &
                                grad%sigma(:, :cut_idx), &
                                exc(:cut_idx))

            IF (xcpot%func_exc_id_c > 0) THEN
               call xc_f03_mgga_exc(xcpot%exc_func_c, SIZE(rh(cut_idx + 1:, :), 1, kind=c_size_t), &
                                    TRANSPOSE(rh(cut_idx + 1:, :)), &
                                    grad%sigma(:, cut_idx + 1:), &
                                    transpose(grad%laplace(cut_idx + 1:, :)), &
                                    kinEnergyDen_libXC(:, cut_idx + 1:), &
                                    excc(cut_idx + 1:))

               call xc_f03_gga_exc(xcpot%vxc_func_c, SIZE(rh(:cut_idx, :), 1, kind=c_size_t), &
                                   TRANSPOSE(rh(:cut_idx, :)), &
                                   grad%sigma(:, :cut_idx), &
                                   excc(:cut_idx))
               exc = exc + excc
            END IF

         ELSE ! first iteration is GGA
            IF (.NOT. PRESENT(grad)) CALL judft_error("Bug: You called get_exc for a MetaGGA potential without providing derivatives")
            CALL xc_f03_gga_exc(xcpot%vxc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, exc)
            IF (xcpot%func_exc_id_c > 0) THEN
               CALL xc_f03_gga_exc(xcpot%vxc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, excc)
               exc = exc + excc
            END IF
         ENDIF

      ELSE
         call juDFT_error("exc is part of a known Family", calledby="xcpot_get_exc@libxc")
      ENDIF

#endif
   END SUBROUTINE xcpot_get_exc

   SUBROUTINE xcpot_get_fxc(xcpot, jspins, rh, fxc)
      USE, INTRINSIC :: IEEE_ARITHMETIC
      use iso_c_binding

      IMPLICIT NONE

      CLASS(t_xcpot_libxc), INTENT(IN)  :: xcpot
      INTEGER,              INTENT(IN)  :: jspins
      REAL,                 INTENT(IN)  :: rh(:, :)
      REAL,                 INTENT(OUT) :: fxc(:, :)

#ifdef CPP_LIBXC
      REAL,ALLOCATABLE  :: fxc_tmp(:,:), fx_tmp(:,:)

      integer(kind=c_size_t)           :: idx

      !libxc uses the spin as a first index, hence we have to transpose....
      ALLOCATE (fxc_tmp(SIZE(fxc, 2), SIZE(fxc, 1))); fxc_tmp = 0.0
      ALLOCATE (fx_tmp(SIZE(fxc, 2), SIZE(fxc, 1))); fx_tmp = 0.0

      IF (xcpot%needs_grad().OR.xcpot%exc_is_MetaGGA()) THEN
         CALL judft_error("Bug: You called get_fxc for a (meta)GGA potential. This is not implemented (yet?).")
      ELSE  !LDA potentials
         CALL xc_f03_lda_fxc(xcpot%vxc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), fx_tmp)
         IF (xcpot%func_vxc_id_c > 0) THEN
            CALL xc_f03_lda_fxc(xcpot%vxc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), fxc_tmp)
            fxc_tmp = fxc_tmp + fx_tmp
         ENDIF
      ENDIF
      fxc = TRANSPOSE(fxc_tmp)

#endif
   END SUBROUTINE xcpot_get_fxc

   SUBROUTINE xcpot_alloc_gradients(ngrid, jspins, grad)
      INTEGER, INTENT(IN)         :: jspins, ngrid
      TYPE(t_gradients), INTENT(INOUT):: grad
      !For libxc we only need the sigma array...
      if (allocated(grad%gr).and..not.allocated(grad%sigma)) return !externally allocated grad%gr
      IF (ALLOCATED(grad%sigma)) DEALLOCATE (grad%sigma, grad%gr, grad%laplace, grad%vsigma)
      ALLOCATE (grad%sigma(MERGE(1, 3, jspins == 1), ngrid))
      ALLOCATE (grad%gr(3, ngrid, jspins))
      ALLOCATE (grad%laplace(ngrid, jspins))
      ALLOCATE (grad%vsigma(MERGE(1, 3, jspins == 1), ngrid))

   END SUBROUTINE xcpot_alloc_gradients

   !> Create a new t_xcpot_libxc in which the auxiliary GGA functionals of a MetaGGA xcpot
   !! are used to initialize the func_vxc_id_* and func_exc_id_* fields.
   !! This is useful to evaluate the auxiliary GGA potential/energy as a standalone xcpot.
   !! Overrides the error-stub in the base class t_xcpot.
   FUNCTION xcpot_create_from_aux(xcpot) RESULT(aux_libxc)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN) :: xcpot
      CLASS(t_xcpot), ALLOCATABLE      :: aux_libxc

#ifdef CPP_LIBXC
    
      IF (.NOT. xcpot%l_has_aux) &
         CALL judft_error("create_from_aux: no auxiliary GGA functional configured in this MetaGGA xcpot")

      ALLOCATE(t_xcpot_libxc :: aux_libxc)
      select type(aux_libxc)
      type is (t_xcpot_libxc)
         aux_libxc%l_libxc        = xcpot%l_libxc
         aux_libxc%l_relativistic = xcpot%l_relativistic
         CALL aux_libxc%init(xcpot%func_aux_id_x, xcpot%func_aux_id_c, &
                          xcpot%func_aux_id_x, xcpot%func_aux_id_c, &
                          xcpot%jspins, .false.)
      end select !No other type possible due to the ALLOCATE statement above
      
#else
      CALL judft_error("create_from_aux requires FLEUR compiled with libxc support")
#endif
   END FUNCTION xcpot_create_from_aux

   subroutine mpi_bc_xcpot_libxc(This, Mpi_comm, Irank)
      Use M_mpi_bc_tool
      Class(t_xcpot_libxc), Intent(Inout)::This
      Integer, Intent(In):: Mpi_comm
      Integer, Intent(In), Optional::Irank
      Integer ::Rank
      If (Present(Irank)) Then
         Rank = Irank
      Else
         Rank = 0
      End If

      ! Bcasts for abstract base class t_xcpot
      CALL mpi_bc(this%l_libxc, rank, mpi_comm)
      CALL mpi_bc(this%func_vxc_id_c, rank, mpi_comm)
      CALL mpi_bc(this%func_vxc_id_x, rank, mpi_comm)
      CALL mpi_bc(this%func_exc_id_c, rank, mpi_comm)
      CALL mpi_bc(this%func_exc_id_x, rank, mpi_comm)
      CALL mpi_bc(this%l_inbuild, rank, mpi_comm)
      CALL mpi_bc(rank, mpi_comm, this%inbuild_name)
      CALL mpi_bc(this%l_relativistic, rank, mpi_comm)
      CALL mpi_bc(this%l_bj, rank, mpi_comm)
      CALL mpi_bc(this%func_aux_id_x, rank, mpi_comm)
      CALL mpi_bc(this%func_aux_id_c, rank, mpi_comm)

   END SUBROUTINE mpi_bc_xcpot_libxc
#ifdef CPP_LIBXC
   SUBROUTINE write_xc_info(xc_func, is_E_func)
      IMPLICIT NONE
      LOGICAL, INTENT(IN), OPTIONAL         :: is_E_func
      INTEGER                             :: i
      CHARACTER(len=120)                  :: kind, family
      LOGICAL                             :: is_energy_func

      TYPE(xc_f03_func_t),INTENT(IN)      :: xc_func
      TYPE(xc_f03_func_info_t)            :: xc_info

      xc_info = xc_f03_func_get_info(xc_func)

      is_energy_func = .FALSE.
      IF (PRESENT(is_E_func)) is_energy_func = is_E_func

      SELECT CASE(xc_f03_func_info_get_kind(xc_info))
      CASE (XC_EXCHANGE)
         WRITE (kind, '(a)') 'an exchange functional'
      CASE (XC_CORRELATION)
         WRITE (kind, '(a)') 'a correlation functional'
      CASE (XC_EXCHANGE_CORRELATION)
         WRITE (kind, '(a)') 'an exchange-correlation functional'
      CASE (XC_KINETIC)
         WRITE (kind, '(a)') 'a kinetic energy functional'
      CASE default
         WRITE (kind, '(a)') 'of unknown kind'
      END SELECT
      SELECT CASE (xc_f03_func_info_get_family(xc_info))
      CASE (XC_FAMILY_LDA);
         WRITE (family, '(a)') "LDA"
      CASE (XC_FAMILY_GGA);
         WRITE (family, '(a)') "GGA"
      CASE (XC_FAMILY_HYB_GGA);
         WRITE (family, '(a)') "hybrid GGA"
      CASE (XC_FAMILY_MGGA);
         WRITE (family, '(a)') "MGGA"
      CASE (XC_FAMILY_HYB_MGGA);
         WRITE (family, '(a)') "hybrid MGGA"
      CASE default;
         WRITE (family, '(a)') "unknown"
      END SELECT

      IF(.not. is_energy_func) THEN
         WRITE(*,'("The functional ''", a, "'' is ", a, ", it belongs to the ''", a, "'' family and is defined in the reference(s):")') &
            TRIM(xc_f03_func_info_get_name(xc_info)), TRIM(kind), TRIM(family)
      ELSE
         WRITE(*,'("The functional used for TotalE ''", a, "'' is ", a, ", it belongs to the ''", a, "'' family and is defined in the reference(s):")') &
            TRIM(xc_f03_func_info_get_name(xc_info)), TRIM(kind), TRIM(family)
      ENDIF

      i = 0
      DO WHILE(i >= 0)
         WRITE(*, '(a,i1,2a)') '[', i+1, '] ', TRIM(xc_f03_func_reference_get_ref(xc_f03_func_info_get_references(xc_info, i)))
      END DO
   END SUBROUTINE write_xc_info

   FUNCTION xc_get_family(xc_func) result(family)
      IMPLICIT NONE
      TYPE(xc_f03_func_t)  :: xc_func
      integer              :: family
      family = xc_f03_func_info_get_family(xc_f03_func_get_info(xc_func))
   END FUNCTION xc_get_family
#endif

END MODULE m_types_xcpot_libxc
