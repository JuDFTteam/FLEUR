MODULE m_dostetra

   USE m_types_kpts
   USE m_types_input
   USE m_tetrahedronInit

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE dostetra(kpts,input,eMesh,eig,qal,g)

      TYPE(t_kpts),        INTENT(IN)     :: kpts
      TYPE(t_input),       INTENT(IN)     :: input
      REAL,                INTENT(IN)     :: eMesh(:)
      REAL,                INTENT(IN)     :: eig(:,:,:)
      REAL,                INTENT(IN)     :: qal(:,:,:)
      REAL,                INTENT(INOUT)  :: g(:,:)

      INTEGER :: ikpt,iBand,ie,ispin,neig(kpts%nkpt)
      REAL    :: w(size(eMesh),size(qal,1))

      g = 0.0
      DO ispin = 1, size(qal,3)
         DO ikpt = 1, kpts%nkpt
            neig=count(eig(:,ikpt,ispin)<1E99)
         ENDDO
         DO ikpt = 1, kpts%nkpt
            !------------------------------------------------------
            ! Calculate the weights for the DOS on the energy Grid
            !------------------------------------------------------
            CALL tetrahedronInit(kpts,input,ikpt,eig(:,:,ispin),MINVAL(neig),eMesh,&
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
