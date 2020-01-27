!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_od_kptsgen
      CONTAINS 
      SUBROUTINE od_kptsgen (nkpt)

!          generates a kpts file in a half of the 
!       one-dimensional Brillouin zone, uniform distribution
!                      Y. Mokrousov Dec 2005 

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nkpt

!  local
      INTEGER i,j
      REAL   , ALLOCATABLE :: kpts(:,:)
      REAL    scale

      scale = 1.0   

      ALLOCATE ( kpts(3,nkpt) )

      OPEN (41,file='kpts',form='formatted',status='new')
      WRITE (41,FMT=8110) nkpt,scale,.false.

      DO i = 1,nkpt

         kpts(1,i) = 0.
         kpts(2,i) = 0.
         IF (nkpt.EQ.1) THEN
            kpts(3,i) = 0.0
         ELSE
            kpts(3,i) = scale*(i-1)/(2.*(nkpt-1))
         ENDIF

         WRITE (41,FMT=8040) (kpts(j,i),j=1,3)

      END DO

 8110 FORMAT (i5,f20.10,3x,l1)
 8040 FORMAT (4f10.5)

      CLOSE (41)

      DEALLOCATE ( kpts )

      END SUBROUTINE od_kptsgen
      END MODULE m_od_kptsgen
