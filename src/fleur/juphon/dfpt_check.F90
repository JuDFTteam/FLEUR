!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_check
   implicit none

IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_check(fi, xcpot)

        USE m_types_fleurinput
        USE m_types_xcpot_libxc
        USE m_juDFT_stop, only : juDFT_error

        TYPE(t_fleurinput), INTENT(IN) :: fi
        CLASS(t_xcpot),     INTENT(IN) :: xcpot

        LOGICAL :: l_libxc
        LOGICAL :: l_exist

        l_libxc = .FALSE.
        l_exist = .FALSE.
        !Symmetry
        IF (fi%sym%nop.GT.1) CALL judft_error("juPhon uses only unit symmetry. Please redo the calculation without symmetry.",calledby="dfpt_check.F90")

        !Coretails
        IF (fi%input%ctail) CALL judft_error("Coretails are not supported. Please consider moving the problematic states to LOs",calledby="dfpt_check.F90")

        !Noco
        IF (fi%noco%l_noco) CALL judft_error("juPhon doesn't do non-collinear systems yet.",calledby="dfpt_check.F90")

        !libxc
        SELECT TYPE(xcpot)
        TYPE IS (t_xcpot_libxc)
            l_libxc=.TRUE.
        END SELECT

        IF (.NOT.l_libxc) CALL judft_error("juPhon needs libxc functionals.",calledby="dfpt_check.F90")

        !GGA
        IF (xcpot%needs_grad()) CALL judft_error("GGA functionals are not supported yet.",calledby="dfpt_check.F90")

        !MetaGGA
        IF (xcpot%is_MetaGGA()) CALL judft_error("juPhon doesn't do MetaGGA functionals.",calledby="dfpt_check.F90")

        !DFTU etc.
        IF ((fi%atoms%n_u.GT.0).OR.(fi%atoms%n_hia.GT.0).OR.(fi%atoms%n_opc.GT.0)) CALL judft_error("Currently juPhon doesn't support DFT+X.",calledby="dfpt_check.F90")

        !SOC:
        IF (fi%noco%l_soc) CALL judft_error("juPhon doesn't support spin-orbit coupling.",calledby="dfpt_check.F90")

        !Spin spirals:
        IF (fi%noco%l_ss) CALL judft_error("juPhon doesn't support spin-spiral systems.",calledby="dfpt_check.F90")

        !vdW
        IF (fi%input%vdw.GT.0) CALL judft_error("juPhon doesn't contain van-der-Waals corrections.",calledby="dfpt_check.F90")

        !Hybrid/RDMFT
        IF (fi%hybinp%l_hybrid .OR. fi%input%l_rdmft) CALL judft_error("juPhon doesn't support hybrid or RDMFT.",calledby="dfpt_check.F90")

        IF (fi%juPhon%l_polar) THEN
            inquire(file="born_eff_charge", exist=l_exist)
            IF (.NOT. l_exist) CALL judft_error("born_eff_charge file not present, required for LO-TO splitting",calledby="dfpt_check.F90")
            inquire(file="diel_tensor", exist=l_exist)
            IF (.NOT. l_exist) CALL judft_error("diel_tensor file not present, required for LO-TO splitting",calledby="dfpt_check.F90")
        END IF
    END SUBROUTINE dfpt_check
END MODULE m_dfpt_check
