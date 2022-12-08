!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!>This module contains the xcpot-type providing an interface to libxc
MODULE m_types_xcpot_libxc
#ifdef CPP_LIBXC
   USE xc_f90_lib_m
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
      TYPE(xc_f90_func_t)      :: vxc_func_x, vxc_func_c
      TYPE(xc_f90_func_t)      :: exc_func_x, exc_func_c
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
      ENDDO
    ENDDO

  end subroutine

   SUBROUTINE xcpot_init(xcpot, func_vxc_id_x, func_vxc_id_c, func_exc_id_x, func_exc_id_c, jspins)
      USE m_judft
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(INOUT)    :: xcpot
      INTEGER, INTENT(IN)                 :: jspins, func_vxc_id_x, func_vxc_id_c, func_exc_id_x, func_exc_id_c
      LOGICAL                             :: same_functionals   ! are vxc and exc equal
      INTEGER                             :: errors(4)

#ifdef CPP_LIBXC
      errors = -1
      xcpot%jspins = jspins
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
         CALL xc_f90_func_init(xcpot%vxc_func_x, xcpot%func_vxc_id_x, XC_UNPOLARIZED, err=errors(1))
         IF (xcpot%func_vxc_id_c>0) CALL xc_f90_func_init(xcpot%vxc_func_c, xcpot%func_vxc_id_c, &
                                                                 XC_UNPOLARIZED, err=errors(2))

         ! energy functionals
         CALL xc_f90_func_init(xcpot%exc_func_x, xcpot%func_exc_id_x, XC_UNPOLARIZED, err=errors(3))
         IF (xcpot%func_exc_id_c>0) CALL xc_f90_func_init(xcpot%exc_func_c, xcpot%func_exc_id_c, &
                                                                  XC_UNPOLARIZED, err=errors(4))

      ELSE
         ! potential functionals
         CALL xc_f90_func_init(xcpot%vxc_func_x, xcpot%func_vxc_id_x, XC_POLARIZED, err=errors(1))
         IF (xcpot%func_vxc_id_c>0) CALL xc_f90_func_init(xcpot%vxc_func_c, xcpot%func_vxc_id_c, &
                                                                  XC_POLARIZED, err=errors(2))

         !energy functionals
         CALL xc_f90_func_init(xcpot%exc_func_x, xcpot%func_exc_id_x, XC_POLARIZED, err=errors(3))
         IF (xcpot%func_exc_id_c>0) CALL xc_f90_func_init(xcpot%exc_func_c, xcpot%func_exc_id_c, &
                                                                  XC_POLARIZED, err=errors(4))
      END IF

      !IF(errors(1) /= 0) call juDFT_error("Exchange potential functional not in LibXC")
      !IF(errors(2) /= 0) call juDFT_error("Correlation potential functional not in LibXC")
      !IF(errors(3) /= 0) call juDFT_error("Exchange energy functional not in LibXC")
      !IF(errors(4) /= 0) call juDFT_error("Correlation energy functional not in LibXC")

      !check if any potental is a MetaGGA
      IF (ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA] == xc_get_family(xcpot%vxc_func_x))) THEN
         CALL juDFT_error("vxc_x: MetaGGA is not implemented for potentials")
      ELSEIF (xcpot%func_vxc_id_c > 0) THEN
         IF (ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA] == xc_get_family(xcpot%vxc_func_c))) THEN
            CALL juDFT_error("vxc_x: MetaGGA is not implemented for potentials")
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
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%vxc_func_x)
      xcpot_vx_is_LDA =  XC_FAMILY_LDA == xc_f90_func_info_get_family(xc_info)
#else
      xcpot_vx_is_LDA = .false.
#endif
   END FUNCTION xcpot_vx_is_LDA

   LOGICAL FUNCTION xcpot_vc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%vxc_func_c)
      xcpot_vc_is_LDA =  XC_FAMILY_LDA == xc_f90_func_info_get_family(xc_info)
#else
      xcpot_vc_is_LDA = .false.
#endif
   END FUNCTION xcpot_vc_is_LDA

   LOGICAL FUNCTION xcpot_exc_is_LDA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%exc_func_x)
      xcpot_exc_is_LDA = (XC_FAMILY_LDA == xc_f90_func_info_get_family(xc_info))
#else
      xcpot_exc_is_LDA = .false.
#endif
   END FUNCTION xcpot_exc_is_LDA

   ! GGA
   LOGICAL FUNCTION xcpot_vc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%vxc_func_c)
      xcpot_vc_is_gga =  ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA]==xc_f90_func_info_get_family(xc_info))
#else
      xcpot_vc_is_gga = .false.
#endif
   END FUNCTION xcpot_vc_is_gga

   LOGICAL FUNCTION xcpot_vx_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%vxc_func_x)
      xcpot_vx_is_gga =  ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA]==xc_f90_func_info_get_family(xc_info))
#else
      xcpot_vx_is_gga = .false.
#endif
   END FUNCTION xcpot_vx_is_gga

   LOGICAL FUNCTION xcpot_vx_is_MetaGGA(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%vxc_func_x)
      xcpot_vx_is_MetaGGA =  ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA]==xc_f90_func_info_get_family(xc_info))
#else
      xcpot_vx_is_MetaGGA = .false.
#endif
   END FUNCTION xcpot_vx_is_MetaGGA

   LOGICAL FUNCTION xcpot_exc_is_gga(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%exc_func_x)
      xcpot_exc_is_gga =  ANY([XC_FAMILY_GGA, XC_FAMILY_HYB_GGA]==xc_f90_func_info_get_family(xc_info))
#else
      xcpot_exc_is_gga = .false.
#endif
   END FUNCTION xcpot_exc_is_gga

   LOGICAL FUNCTION xcpot_exc_is_MetaGGA(xcpot)
      IMPLICIT NONE
   CLASS(t_xcpot_libxc),INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%exc_func_x)
      xcpot_exc_is_MetaGGA=ANY([XC_FAMILY_MGGA, XC_FAMILY_HYB_MGGA]==xc_f90_func_info_get_family(xc_info))
#else
      xcpot_exc_is_MetaGGA = .False.
#endif
   END FUNCTION xcpot_exc_is_MetaGGA

   LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
      IMPLICIT NONE
      CLASS(t_xcpot_libxc), INTENT(IN):: xcpot
#ifdef CPP_LIBXC
      TYPE(xc_f90_func_info_t)        :: xc_info

      xc_info = xc_f90_func_get_info(xcpot%vxc_func_x)
      xcpot_is_hybrid=ANY([XC_FAMILY_HYB_MGGA, XC_FAMILY_HYB_GGA]==xc_f90_func_info_get_family(xc_info))
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
      a_ex=xc_f90_hyb_exx_coef(xcpot%vxc_func_x)
#endif
   END FUNCTION xcpot_get_exchange_weight

   !***********************************************************************
   SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh, vxc,vx, grad, kinenergyden_ks)
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
#ifdef CPP_LIBXC
      REAL,ALLOCATABLE  :: vxc_tmp(:,:),vx_tmp(:,:),vsigma(:,:), &
                           tmp_vsig(:,:), tmp_vlapl(:,:), tmp_vtau(:,:), &
                           kinED_libxc(:,:)
      integer(kind=c_size_t)           :: idx
      !libxc uses the spin as a first index, hence we have to transpose....
      ALLOCATE (vxc_tmp(SIZE(vxc, 2), SIZE(vxc, 1))); vxc_tmp = 0.0
      ALLOCATE (vx_tmp(SIZE(vx, 2), SIZE(vx, 1))); vx_tmp = 0.0

      IF (xcpot%needs_grad()) THEN
         IF (.NOT. PRESENT(grad)) CALL judft_error("Bug: You called get_vxc for a GGA potential without providing derivatives")
         ALLOCATE (vsigma, mold=grad%vsigma)
         !where(abs(grad%sigma)<1E-9) grad%sigma=1E-9
         CALL xc_f90_gga_vxc(xcpot%vxc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, vx_tmp, vsigma)
         IF (xcpot%func_vxc_id_c > 0) THEN
            CALL xc_f90_gga_vxc(xcpot%vxc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, vxc_tmp, grad%vsigma)
            grad%vsigma = grad%vsigma + vsigma
            vxc_tmp = vxc_tmp + vx_tmp
         ELSE
            vxc_tmp = vx_tmp
         ENDIF
      ELSE  !LDA potentials
         CALL xc_f90_lda_vxc(xcpot%vxc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), vx_tmp)
         IF (xcpot%func_vxc_id_c > 0) THEN
            CALL xc_f90_lda_vxc(xcpot%vxc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), vxc_tmp)
            vxc_tmp = vxc_tmp + vx_tmp
         ENDIF
      ENDIF
      vx = TRANSPOSE(vx_tmp)
      vxc = TRANSPOSE(vxc_tmp)

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
      TYPE(xc_f90_func_info_t)       :: xc_info
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
         CALL xc_f90_gga_exc(xcpot%exc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, exc)
         IF (xcpot%func_exc_id_c > 0) THEN
            CALL xc_f90_gga_exc(xcpot%exc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, excc)
            exc = exc + excc
         END IF
      ELSEIF (xcpot%exc_is_LDA()) THEN  !LDA potentials
         CALL xc_f90_lda_exc(xcpot%exc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), exc)
         IF (xcpot%func_exc_id_c > 0) THEN
            CALL xc_f90_lda_exc(xcpot%exc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), excc)
            exc = exc + excc
         END IF
      ELSEIF (xcpot%exc_is_MetaGGA()) THEN
         IF (PRESENT(kinEnergyDen_KS)) THEN
            ! apply correction in  eq (4) in https://doi.org/10.1063/1.1565316
            kinEnergyDen_libXC = transpose(kinEnergyDen_KS + 0.25*grad%laplace)

            !only cut core of muffin tin
            cut_idx = MERGE(NINT(size(rh, 1)*cut_ratio), 0, is_mt)

            exc = 0.0
            excc = 0.0
            call xc_f90_mgga_exc(xcpot%exc_func_x, SIZE(rh(cut_idx + 1:, :), 1, kind=c_size_t), &
                                 TRANSPOSE(rh(cut_idx + 1:, :)), &
                                 grad%sigma(:, cut_idx + 1:), &
                                 transpose(grad%laplace(cut_idx + 1:, :)), &
                                 kinEnergyDen_libXC(:, cut_idx + 1:), &
                                 exc(cut_idx + 1:))

            call xc_f90_gga_exc(xcpot%vxc_func_x, SIZE(rh(:cut_idx, :), 1, kind=c_size_t), &
                                TRANSPOSE(rh(:cut_idx, :)), &
                                grad%sigma(:, :cut_idx), &
                                exc(:cut_idx))

            IF (xcpot%func_exc_id_c > 0) THEN
               call xc_f90_mgga_exc(xcpot%exc_func_c, SIZE(rh(cut_idx + 1:, :), 1, kind=c_size_t), &
                                    TRANSPOSE(rh(cut_idx + 1:, :)), &
                                    grad%sigma(:, cut_idx + 1:), &
                                    transpose(grad%laplace(cut_idx + 1:, :)), &
                                    kinEnergyDen_libXC(:, cut_idx + 1:), &
                                    excc(cut_idx + 1:))

               call xc_f90_gga_exc(xcpot%vxc_func_c, SIZE(rh(:cut_idx, :), 1, kind=c_size_t), &
                                   TRANSPOSE(rh(:cut_idx, :)), &
                                   grad%sigma(:, :cut_idx), &
                                   excc(:cut_idx))
               exc = exc + excc
            END IF

         ELSE ! first iteration is GGA
            IF (.NOT. PRESENT(grad)) CALL judft_error("Bug: You called get_exc for a MetaGGA potential without providing derivatives")
            CALL xc_f90_gga_exc(xcpot%vxc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, exc)
            IF (xcpot%func_exc_id_c > 0) THEN
               CALL xc_f90_gga_exc(xcpot%vxc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), grad%sigma, excc)
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
         CALL xc_f90_lda_fxc(xcpot%vxc_func_x, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), fx_tmp)
         IF (xcpot%func_vxc_id_c > 0) THEN
            CALL xc_f90_lda_fxc(xcpot%vxc_func_c, SIZE(rh, 1, kind=c_size_t), TRANSPOSE(rh), fxc_tmp)
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
      IF (ALLOCATED(grad%sigma)) DEALLOCATE (grad%sigma, grad%gr, grad%laplace, grad%vsigma)
      ALLOCATE (grad%sigma(MERGE(1, 3, jspins == 1), ngrid))
      ALLOCATE (grad%gr(3, ngrid, jspins))
      ALLOCATE (grad%laplace(ngrid, jspins))
      ALLOCATE (grad%vsigma(MERGE(1, 3, jspins == 1), ngrid))

   END SUBROUTINE xcpot_alloc_gradients

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

   END SUBROUTINE mpi_bc_xcpot_libxc
#ifdef CPP_LIBXC
   SUBROUTINE write_xc_info(xc_func, is_E_func)
      IMPLICIT NONE
      LOGICAL, INTENT(IN), OPTIONAL         :: is_E_func
      INTEGER                             :: i
      CHARACTER(len=120)                  :: kind, family
      LOGICAL                             :: is_energy_func

      TYPE(xc_f90_func_t),INTENT(IN)      :: xc_func
      TYPE(xc_f90_func_info_t)            :: xc_info

      xc_info = xc_f90_func_get_info(xc_func)
      is_energy_func = merge(is_E_func, .False., PRESENT(is_E_func))

      SELECT CASE(xc_f90_func_info_get_kind(xc_info))
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
      SELECT CASE (xc_f90_func_info_get_family(xc_info))
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
            TRIM(xc_f90_func_info_get_name(xc_info)), TRIM(kind), TRIM(family)
      ELSE
         WRITE(*,'("The functional used for TotalE ''", a, "'' is ", a, ", it belongs to the ''", a, "'' family and is defined in the reference(s):")') &
            TRIM(xc_f90_func_info_get_name(xc_info)), TRIM(kind), TRIM(family)
      ENDIF

      i = 0
      DO WHILE(i >= 0)
         WRITE(*, '(a,i1,2a)') '[', i+1, '] ', TRIM(xc_f90_func_reference_get_ref(xc_f90_func_info_get_references(xc_info, i)))
      END DO
   END SUBROUTINE write_xc_info

   FUNCTION xc_get_family(xc_func) result(family)
      IMPLICIT NONE
      TYPE(xc_f90_func_t)  :: xc_func
      integer              :: family
      family = xc_f90_func_info_get_family(xc_f90_func_get_info(xc_func))
   END FUNCTION xc_get_family
#endif

END MODULE m_types_xcpot_libxc
