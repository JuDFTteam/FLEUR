!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_cell
  USE m_judft
  IMPLICIT NONE
  PRIVATE
  TYPE t_cell
      !name of 2D-lattice type
      CHARACTER*3::latnam
      !vol of dtilde box
      REAL::omtil
      !2D area
      REAL::area
      !bravais matrix
      REAL::amat(3, 3)
      !rez. bravais matrx
      REAL::bmat(3, 3)
      !square of bbmat
      REAL::bbmat(3, 3),aamat(3,3)
      !d-value
      REAL::z1
      !volume of cell
      REAL::vol
      !volume of interstitial
      REAL::volint
    CONTAINS
      PROCEDURE :: init
   END TYPE t_cell
   PUBLIC t_cell
 CONTAINS
   SUBROUTINE init(cell,dvac)
     !initialize cell, only input it cell%amat and dvac in case of a film
     USE m_inv3
     USE m_constants,ONLY:tpi_const
     CLASS (t_cell),INTENT(INOUT):: cell
     REAL,INTENT(IN)             :: dvac


     CALL inv3(cell%amat,cell%bmat,cell%omtil)
     IF (cell%omtil<0) CALL judft_warn("Negative volume! You are using a left-handed coordinate system")
     cell%omtil=ABS(cell%omtil)
     
     cell%bmat=tpi_const*cell%bmat
     IF (dvac>0) THEN
        cell%vol = (cell%omtil/cell%amat(3,3))*dvac
        cell%area = cell%omtil/cell%amat(3,3)
     ELSE
        cell%vol = cell%omtil
        cell%area = cell%amat(1,1)*cell%amat(2,2)-cell%amat(1,2)*cell%amat(2,1)
        IF (cell%area < 1.0e-7) THEN
           cell%area = 1.
           CALL juDFT_warn("area = 0",calledby ="types_cell")
        END IF
     END IF
     cell%z1=dvac/2

     cell%bbmat=matmul(cell%bmat,transpose(cell%bmat))
     cell%aamat=matmul(transpose(cell%amat),cell%amat)
   END SUBROUTINE init
 END MODULE m_types_cell
