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
        USE m_types_xcpot_libxc
        USE m_juDFT_stop, only : juDFT_error

        TYPE(t_fleurinput), INTENT(IN) :: fi
        CLASS(t_xcpot),     INTENT(IN) :: xcpot

        LOGICAL :: l_libxc

        l_libxc = .FALSE.

        !Symmetry
        IF (fi%sym%nop.GT.1) CALL judft_error("juPhon uses only unit symmetry.")

        !Coretails
        IF (fi%input%ctail) CALL judft_error("juPhon coretails are problematic at the moment.")

        !!LOs
        !IF (ANY(fi%atoms%nlo.GT.0)) THEN
        !    CALL judft_error("juPhon doesn't do local orbitals yet.")
        !END IF

        !Noco
        IF (fi%noco%l_noco) CALL judft_error("juPhon doesn't do non-collinear systems yet.")

        !libxc
        SELECT TYPE(xcpot)
        TYPE IS (t_xcpot_libxc)
            l_libxc=.TRUE.
        END SELECT

        IF (.NOT.l_libxc) CALL judft_error("juPhon needs libxc functionals.")

        !GGA
        IF (xcpot%needs_grad()) CALL judft_error("juPhon doesn't do GGA functionals yet.")

        !MetaGGA
        IF (xcpot%exc_is_MetaGGA() .or. xcpot%vx_is_MetaGGA()) CALL judft_error("juPhon doesn't do MetaGGA functionals.")

        !DFTU etc.
        IF ((fi%atoms%n_u.GT.0).OR.(fi%atoms%n_hia.GT.0).OR.(fi%atoms%n_opc.GT.0)) CALL judft_error("juPhon doesn't do DFT+X [yet].")

        !SOC:
        IF (fi%noco%l_soc) CALL judft_error("juPhon doesn't do spin-orbit coupling yet.")

        !Spin spirals:
        IF (fi%noco%l_ss) CALL judft_error("juPhon doesn't do spin-spiral systems [yet].")

        !vdW
        IF (fi%input%vdw.GT.0) CALL judft_error("juPhon doesn't do van-der-Waals systems.")

        !Film
        IF (fi%input%film) CALL judft_error("juPhon doesn't do film systems.")

        !Hybrid/RDMFT
        IF (fi%hybinp%l_hybrid .OR. fi%input%l_rdmft) CALL judft_error("juPhon doesn't do hybrid or RDMFT.")
    END SUBROUTINE dfpt_check
END MODULE m_dfpt_check
