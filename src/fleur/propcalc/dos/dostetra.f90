!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dostetra

   USE m_types_kpts
   USE m_types_input
   USE m_tetrahedronInit

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: dostetra

   CONTAINS

   SUBROUTINE dostetra(kpts,input,eMesh,eig,qal,g,energyShift)

      TYPE(t_kpts),        INTENT(IN)     :: kpts
      TYPE(t_input),       INTENT(IN)     :: input
      REAL,                INTENT(IN)     :: eMesh(:)
      REAL,                INTENT(IN)     :: eig(:,:,:)
      REAL,                INTENT(IN)     :: qal(:,:,:)
      REAL,                INTENT(INOUT)  :: g(:,:)
      REAL, OPTIONAL,      INTENT(IN)     :: energyShift

      INTEGER :: ikpt,iBand,ie,ispin,neig(kpts%nkpt)
      REAL    :: w(size(eMesh),size(qal,1))
      REAL    :: shift
      REAL, ALLOCATABLE :: eig_shifted(:,:,:)

      shift = 0.0
      IF(PRESENT(energyShift)) shift = energyShift

      eig_shifted = eig - shift

      g = 0.0
      DO ispin = 1, size(qal,3)
         DO ikpt = 1, kpts%nkpt
            neig=count(eig(:,ikpt,ispin)<1E99)
         ENDDO
         DO ikpt = 1, kpts%nkpt
            !------------------------------------------------------
            ! Calculate the weights for the DOS on the energy Grid
            !------------------------------------------------------
            CALL tetrahedronInit(kpts,input,ikpt,eig_shifted(:,:,ispin),MINVAL(neig),eMesh,&
                                 w,dos=.TRUE.)
            DO iBand = 1, neig(ikpt)
               DO ie = 1, SIZE(eMesh)
                  g(ie,ispin) = g(ie,ispin) + w(ie,iBand) * 2.0/input%jspins * qal(iBand,ikpt,ispin)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE dostetra
END MODULE m_dostetra
