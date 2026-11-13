!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_distspins
  USE m_juDFT
  USE m_types_mat
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: hsmt_distspins
CONTAINS
  !> Adds chi(igSpinPr,igSpin)*mat_tmp to the spin blocks mat(igSpinPr,igSpin)
  !!
  !! Blocks with a vanishing spinor factor are skipped.
  !! l_upper: only the upper triangle (including the diagonal) is added. Only valid for
  !! undistributed matrices whose lower triangle is not used later.
  SUBROUTINE hsmt_distspins(chi,mat_tmp,mat,l_upper)
    COMPLEX,INTENT(in)          :: chi(2,2)
    CLASS(t_mat),INTENT(IN)     :: mat_tmp
    CLASS(t_mat),INTENT(INOUT)  :: mat(:,:)
    LOGICAL,INTENT(IN),OPTIONAL :: l_upper

    REAL,PARAMETER :: tol_chi = 1.0e-15
    INTEGER :: igSpinPr,igSpin,j,nrow
    LOGICAL :: l_up,l_use(2,2)

#ifndef _OPENACC
    IF (.NOT.mat_tmp%l_real) THEN
       l_up = .FALSE.
       IF (PRESENT(l_upper)) l_up = l_upper
       l_use = ABS(chi) > tol_chi
       DO igSpin = 1,2
          DO igSpinPr = 1,2
             IF (.NOT.l_use(igSpinPr,igSpin)) CYCLE
             IF (mat(igSpinPr,igSpin)%l_real) CALL judft_error("hsmt_distspins: matrix type mismatch")
             IF (mat(igSpinPr,igSpin)%matsize1 /= mat_tmp%matsize1 .OR. mat(igSpinPr,igSpin)%matsize2 /= mat_tmp%matsize2) &
                CALL judft_error("hsmt_distspins: matrix size mismatch")
          ENDDO
       ENDDO
       !$OMP PARALLEL DO DEFAULT(none) SHARED(mat,mat_tmp,chi,l_use,l_up) PRIVATE(j,nrow,igSpinPr,igSpin) SCHEDULE(dynamic,16)
       DO j = 1,mat_tmp%matsize2
          nrow = MERGE(MIN(j,mat_tmp%matsize1),mat_tmp%matsize1,l_up)
          DO igSpin = 1,2
             DO igSpinPr = 1,2
                IF (.NOT.l_use(igSpinPr,igSpin)) CYCLE
                mat(igSpinPr,igSpin)%data_c(:nrow,j) = mat(igSpinPr,igSpin)%data_c(:nrow,j) &
                                                     + chi(igSpinPr,igSpin)*mat_tmp%data_c(:nrow,j)
             ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
       RETURN
    ENDIF
#endif
    DO igSpinPr=1,2
      DO igSpin=1,2
        call mat(igSpinPr,igSpin)%add(mat_tmp,chi(igSpinPr,igSpin))
      enddo
    ENDDO

  END SUBROUTINE hsmt_distspins
END MODULE m_hsmt_distspins
