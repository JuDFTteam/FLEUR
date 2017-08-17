!--------------------------------------------------------------------------------
! Copyright (c) 2017 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_bs_comfort

      CONTAINS

      SUBROUTINE bs_comfort(eig_id,DIMENSION,input,noco,nkpt,param)

      USE m_types
      USE m_eig66_io, ONLY : read_eig

      IMPLICIT NONE

      TYPE(t_dimension), INTENT(IN)   :: DIMENSION
      TYPE(t_input),     INTENT(IN)   :: input
      TYPE(t_noco),      INTENT(IN)   :: noco

C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: eig_id
      INTEGER, INTENT (IN) :: nkpt
      INTEGER, INTENT (IN) :: param

C     .. Local Scalars ..
      INTEGER              :: i, jsp, k, nspins

C     .. Local Arrays ..
      REAL,    ALLOCATABLE :: eig(:,:,:)
      INTEGER, ALLOCATABLE :: ne(:,:)

c---> pk non-collinear
      IF (noco%l_noco) THEN
         nspins = 1
      ELSE
         nspins = input%jspins
      ENDIF
c---> pk non-collinear

      ALLOCATE (eig(DIMENSION%neigd,nkpt,nspins))
      ALLOCATE (ne(nkpt,nspins))

      DO jsp = 1, nspins
         DO k = 1, nkpt
            CALL read_eig(eig_id,k,jsp,neig=ne(k,jsp),eig=eig(:,k,jsp))
         END DO ! k = 1,nkpt

         DO i = 1, DIMENSION%neigd
            DO k = 1, nkpt
               WRITE(776+jsp,*) param,k,eig(i,k,jsp)
            END DO
            WRITE(776+jsp,*)
         END DO

         DO i = 1, DIMENSION%neigd
            DO k = 1, nkpt
               IF (k.eq.param) WRITE(778+jsp,*) param,k,eig(i,k,jsp)
            END DO
            WRITE(778+jsp,*)
         END DO
      END DO ! jsp = 1,nspins

      DEALLOCATE (ne,eig)

      END SUBROUTINE bs_comfort
      END MODULE m_bs_comfort
