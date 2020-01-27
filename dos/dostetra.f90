MODULE m_dostetra

   USE m_types
   USE m_tetrahedronInit

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE dostetra(kpts,input,ndos,ne,eMesh,neig,eig,qal,g)

      TYPE(t_kpts),        INTENT(IN)     :: kpts
      TYPE(t_input),       INTENT(IN)     :: input
      INTEGER,             INTENT(IN)     :: ndos
      INTEGER,             INTENT(IN)     :: ne
      REAL,                INTENT(IN)     :: eMesh(:)
      INTEGER,             INTENT(IN)     :: neig(:)
      REAL,                INTENT(IN)     :: eig(:,:)
      REAL,                INTENT(IN)     :: qal(:,:,:)
      REAL,                INTENT(INOUT)  :: g(:,:)

      INTEGER :: ikpt,iBand,idos,ie
      REAL    :: w(ne,MAXVAL(neig))

      g = 0.0
      DO ikpt = 1, kpts%nkpt
         !------------------------------------------------------
         ! Calculate the weights for the DOS on the energy Grid
         !------------------------------------------------------
         CALL tetrahedronInit(kpts,ikpt,eig,neig,eMesh,ne,&
                              input%film,w,dos=.TRUE.)
         DO iBand = 1, neig(ikpt)
            DO idos = 1, ndos
               DO ie = 1, ne
                  g(ie,idos) = g(ie,idos) + w(ie,iBand) * 2.0/input%jspins * qal(idos,iBand,ikpt)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

   END SUBROUTINE dostetra
END MODULE m_dostetra