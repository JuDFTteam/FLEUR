!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_tetrahedron_regular

   !--------------------------------------------------------------------
   !Decomposes an equidistant k-point grid into tetrahedra
   !This happens by looking at a cube of 8 neighbouring k-points
   !and dividing this cube into 6 tetrahedra along the shortest diagonal
   ! Largely equivalent to spex routine
   !--------------------------------------------------------------------

   USE m_types_kpts
   USE m_types_cell
   USE m_juDFT
   USE m_constants

   IMPLICIT NONE

   CONTAINS

   SUBROUTINE tetrahedron_regular(kpts,cell,grid,ntetra,voltet)

      TYPE(t_kpts),           INTENT(INOUT)  :: kpts
      TYPE(t_cell),           INTENT(IN)     :: cell
      INTEGER,                INTENT(IN)     :: grid(:)
      INTEGER, ALLOCATABLE,   INTENT(INOUT)  :: ntetra(:,:)
      REAL,    ALLOCATABLE,   INTENT(INOUT)  :: voltet(:)


      INTEGER :: ntetraCube,k1,k2,k3,ikpt,itetra
      REAL    :: vol,sumvol,volbz
      INTEGER :: tetra(4,24),iarr(3),kcorn(8)
      INTEGER, ALLOCATABLE :: p(:,:,:)

      volbz = ABS(det(cell%bmat))
      !Choose the tetrahedra decomposition along the shortest diagonal
      CALL get_tetra(kpts,cell,grid,ntetraCube,vol,tetra)
      !Set up pointer array for the kpts
      ALLOCATE(p(0:grid(1),0:grid(2),0:grid(3)),source=0)
      p = 0
      DO ikpt = 1, kpts%nkptf
         iarr = nint(kpts%bkf(:,ikpt)*grid)
         p(iarr(1),iarr(2),iarr(3)) = ikpt
      ENDDO
      p(grid(1),:,:) = p(0,:,:)
      p(:,grid(2),:) = p(:,0,:)
      p(:,:,grid(3)) = p(:,:,0)


      !Check for invalid indices
      IF(ANY(p<=0).OR.ANY(p>kpts%nkptf)) THEN
         CALL juDFT_error("Invalid kpoint index in pointer array",calledby="tetrahedron_regular")
      ENDIF

      !Temporary Size
      ALLOCATE(ntetra(4,kpts%nkptf*6),source=0)
      ALLOCATE(voltet(kpts%nkptf*6),source=0.0)

      kpts%ntet = 0
      sumvol = 0.0
      !Set up the tetrahedrons
      DO k3 = 0, grid(3)-1
         DO k2 = 0, grid(2)-1
            DO k1 = 0, grid(1)-1
               !Corners of the current cube
               kcorn(1) = p(k1  ,k2  ,k3  )
               kcorn(2) = p(k1+1,k2  ,k3  )
               kcorn(3) = p(k1  ,k2+1,k3  )
               kcorn(4) = p(k1+1,k2+1,k3  )
               kcorn(5) = p(k1  ,k2  ,k3+1)
               kcorn(6) = p(k1+1,k2  ,k3+1)
               kcorn(7) = p(k1  ,k2+1,k3+1)
               kcorn(8) = p(k1+1,k2+1,k3+1)

               !Now divide the cube into tetrahedra
               DO itetra = 1, ntetraCube
                  !Drop all tetrahedra without kpoints inside the IBZ
                  sumvol = sumvol + vol
                  IF(ALL(kcorn(tetra(1:4,itetra)).GT.kpts%nkpt)) CYCLE
                  kpts%ntet = kpts%ntet+1
                  ntetra(1:4,kpts%ntet) = kcorn(tetra(1:4,itetra))
                  voltet(kpts%ntet) = vol
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      !Has the whole brillouin zone been covered?
      IF(ABS(sumvol-volbz).GT.1E-10) THEN
            CALL juDFT_error("tetrahedron_regular failed", calledby="tetrahedron_regular")
      ENDIF
      voltet = voltet/volbz

      !Rescale volumes for IO to inp.xml
      voltet = voltet * kpts%ntet

   END SUBROUTINE tetrahedron_regular


   SUBROUTINE get_tetra(kpts,cell,grid,ntetra,vol,tetra)

      TYPE(t_kpts),  INTENT(IN)     :: kpts
      TYPE(t_cell),  INTENT(IN)     :: cell
      INTEGER,       INTENT(IN)     :: grid(:)
      INTEGER,       INTENT(INOUT)  :: ntetra
      REAL,          INTENT(INOUT)  :: vol
      INTEGER,       INTENT(INOUT)  :: tetra(:,:)

      REAL rlv(3,3),diag(4),d(3)
      INTEGER idmin

      !Calculate the lengths of the three diagonals
      rlv(:,1) = cell%bmat(:,1) / grid
      rlv(:,2) = cell%bmat(:,2) / grid
      rlv(:,3) = cell%bmat(:,3) / grid

      vol = 1/6.0*ABS(det(rlv))
      d = rlv(:,1) + rlv(:,3) - rlv(:,2)
      diag(1) = sum(d*d)
      d = rlv(:,2) + rlv(:,3) - rlv(:,1)
      diag(2) = sum(d*d)
      d = rlv(:,1) + rlv(:,2) + rlv(:,3)
      diag(3) = sum(d*d)
      d = rlv(:,1) + rlv(:,2) - rlv(:,3)
      diag(4) = sum(d*d)
      idmin = minloc(diag,1)

      ntetra = 0
      !From spex tetrahedron.f (For now we only choose one decomposition)
      if(idmin==1) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 1,2,3,6, 5,7,3,6, 1,5,3,6, 2,4,3,6, 4,8,3,6, 7,8,3,6 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(idmin==2) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 5,6,2,7, 1,5,2,7, 1,3,2,7, 8,6,2,7, 4,3,2,7, 8,4,2,7 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(idmin==3) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 2,6,1,8, 2,4,1,8, 3,4,1,8, 3,7,1,8, 5,7,1,8, 5,6,1,8 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif
      if(idmin==4) then
        tetra(:,ntetra+1:ntetra+6) = reshape ( [ 2,6,4,5, 1,2,4,5, 1,3,4,5, 3,7,4,5, 7,8,4,5, 6,8,4,5 ], [ 4,6 ] )
        ntetra                     = ntetra + 6
      endif

   END SUBROUTINE get_tetra

   REAL FUNCTION det(m)
      REAL m(:,:)
      det = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) + &
            m(2,1)*m(3,2)*m(1,3) - m(1,3)*m(2,2)*m(3,1) - &
            m(2,3)*m(3,2)*m(1,1) - m(2,1)*m(1,2)*m(3,3)
   END FUNCTION det

END MODULE m_tetrahedron_regular