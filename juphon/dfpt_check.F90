!--------------------------------------------------------------------------------
! Copyright (c) 2021 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_check

IMPLICIT NONE

CONTAINS
    SUBROUTINE dfpt_check(fi, xcpot)

        USE m_types_fleurinput
        USE m_juDFT_stop, only : juDFT_error

        TYPE(t_fleurinput), INTENT(IN) :: fi
        CLASS(t_xcpot),     INTENT(IN) :: xcpot

        !Symmetry
        IF (fi%sym%nop.GT.1) THEN
            CALL judft_error("juPhon uses only unit symmetry.")
        END IF

        !Coretails
        IF (fi%input%ctail) THEN
            CALL judft_error("juPhon coretails are problematic at the moment.")
        END IF

        !LOs
        IF (ANY(fi%atoms%nlo.GT.0)) THEN
            CALL judft_error("juPhon doesn't do local orbitals yet.")
        END IF

        !Magnetic
        IF (fi%input%jspins.GT.1) THEN
            CALL judft_error("juPhon doesn't do spin polarized systems yet.")
        END IF

        !Polyatomic
        IF (fi%atoms%nat.GT.1) THEN
            CALL judft_error("juPhon doesn't do polyatomic systems yet.")
        END IF

        !Noco
        IF (fi%noco%l_noco) THEN
            CALL judft_error("juPhon doesn't do non-collinear systems yet.")
        END IF

        !GGA
        IF (xcpot%needs_grad()) THEN
            CALL judft_error("juPhon doesn't do GGA functionals [yet].")
        END IF

        !MetaGGA
        IF (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA()) THEN
            CALL judft_error("juPhon doesn't do MetaGGA functionals.")
        END IF

        !DFTU etc.
        IF ((fi%atoms%n_u.GT.0).OR.(fi%atoms%n_hia.GT.0).OR.(fi%atoms%n_opc.GT.0)) THEN
            CALL judft_error("juPhon doesn't do DFT+X [yet].")
        END IF

        !SOC:
        IF (fi%noco%l_soc) THEN
            CALL judft_error("juPhon doesn't do spin-orbit coupling [yet].")
        END IF

        !Spin spirals:
        IF (fi%noco%l_ss) THEN
            CALL judft_error("juPhon doesn't do spin-spiral systems [yet].")
        END IF

        !vdW
        IF (fi%input%vdw.GT.0) THEN
            CALL judft_error("juPhon doesn't do van-der-Waals systems.")
        END IF
        !Film
        IF (fi%input%film) THEN
            CALL judft_error("juPhon doesn't do film systems.")
        END IF

        !oneD
        IF (fi%oneD%odi%d1) THEN
            CALL judft_error("juPhon doesn't do 1D systems.")
        END IF

        !Hybrid/RDMFT
        IF (fi%hybinp%l_hybrid .OR. fi%input%l_rdmft) THEN
            CALL judft_error("juPhon doesn't do hybrid or RDMFT.")
        END IF
    END SUBROUTINE dfpt_check
END MODULE m_dfpt_check
