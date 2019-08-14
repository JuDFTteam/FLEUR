MODULE m_calc_tetra

   !Calculate tetraeder from equidistant grid ref Phys. Rev. B 49, (16223).
   !Largely equivalent to spex routines 

   CONTAINS 

   SUBROUTINE calc_tetra(kpts,cell,input,sym)

      USE m_types
      USE m_juDFT
      USE m_constants

      IMPLICIT NONE

      TYPE(t_kpts),  INTENT(INOUT) :: kpts 
      TYPE(t_input), INTENT(INOUT)    :: input
      TYPE(t_sym),   INTENT(IN)    :: sym
      TYPE(t_cell),  INTENT(IN)    :: cell

      INTEGER p(0:kpts%nkpt3(1),0:kpts%nkpt3(2),0:kpts%nkpt3(3))
      INTEGER tetra(4,24), ntetra,k1,k2,k3,ikpt,itetra,iarr(3)
      REAL kcorn(8),vol
      REAL sumvol,volbz

      CALL timestart("Calculation of Tetrahedra")

      volbz = sqrt(sum(cell%bmat(:,1)*cell%bmat(:,1))*sum(cell%bmat(:,2)*cell%bmat(:,2))*sum(cell%bmat(:,3)*cell%bmat(:,3)))
      !Choose the tetrahedra decomposition along the shortest diagonal
      CALL get_tetra(tetra,ntetra,kpts,cell,vol)
      !MISSING: Generate all kpts 
      !Set up pointer array for the kpts 
      p = 0
      DO ikpt = 1, kpts%nkptf 
         iarr = nint(kpts%bkf(:,ikpt)*kpts%nkpt3)
         p(iarr(1),iarr(2),iarr(3)) = ikpt
      ENDDO
      !Wrap around at the end
      IF(ALL(p(kpts%nkpt3(1),:,:).EQ.0)) THEN
         p(kpts%nkpt3(1),:,:) = p(0,:,:)
         p(:,kpts%nkpt3(2),:) = p(:,0,:)
         p(:,:,kpts%nkpt3(3)) = p(:,:,0)
      ELSE
         p(0,:,:) = p(kpts%nkpt3(1),:,:)
         p(:,0,:) = p(:,kpts%nkpt3(2),:)
         p(:,:,0) = p(:,:,kpts%nkpt3(3))
      ENDIF
      !Check for invalid indices
      IF(ANY(p<=0).OR.ANY(p>kpts%nkptf)) CALL juDFT_error("Invalid kpoint index in pointer array",calledby="calc_tetra")

      IF (ALLOCATED(kpts%ntetra)) THEN
         DEALLOCATE(kpts%ntetra)
      END IF
      IF (ALLOCATED(kpts%voltet)) THEN
         DEALLOCATE(kpts%voltet)
      END IF
      ALLOCATE(kpts%ntetra(4,kpts%nkptf*6))
      ALLOCATE(kpts%voltet(kpts%nkptf*6))

      kpts%ntet = 0
      sumvol = 0.0
      !Set up the tetrahedrons
      DO k3 = 0, kpts%nkpt3(3)-1
         DO k2 = 0, kpts%nkpt3(2)-1
            DO k1 = 0, kpts%nkpt3(1)-1
               !Corners of the current cube
               kcorn(1) = p(k1  ,k2  ,k3  )
               kcorn(2) = p(k1+1,k2  ,k3  )
               kcorn(3) = p(k1  ,k2+1,k3  )
               kcorn(4) = p(k1+1,k2+1,k3  )
               kcorn(5) = p(k1  ,k2  ,k3+1)
               kcorn(6) = p(k1+1,k2  ,k3+1)
               kcorn(7) = p(k1  ,k2+1,k3+1)
               kcorn(8) = p(k1+1,k2+1,k3+1)

               !Now write the information about the tetrahedron into the corresponding arrays in kpts
               DO itetra = 1, 6
                  kpts%ntet = kpts%ntet+1
                  kpts%ntetra(1:4,kpts%ntet) = kcorn(tetra(1:4,itetra))
                  kpts%voltet(kpts%ntet) = vol
                  sumvol = sumvol + vol
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      IF(ABS(sumvol-volbz).GT.1E-12) CALL juDFT_error("calc_tetra failed", calledby="calc_tetra")
      input%gfTet = .TRUE.
      kpts%voltet = kpts%voltet*kpts%ntet/volbz
      CALL timestop("Calculation of Tetrahedra")

   END SUBROUTINE calc_tetra


   SUBROUTINE get_tetra(tetra,ntetra,kpts,cell,vol)

      USE m_types

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: tetra(4,*)
      INTEGER, INTENT(OUT) :: ntetra
      TYPE(t_kpts), INTENT(IN) :: kpts 
      TYPE(t_cell), INTENT(IN) :: cell
      REAL,         INTENT(OUT):: vol

      REAL rlv(3,3),diag(4),d(3)
      INTEGER idmin

      !Calculate the lengths of the three diagonals
      rlv(:,1) = cell%bmat(:,1) / kpts%nkpt3(1)
      rlv(:,2) = cell%bmat(:,2) / kpts%nkpt3(2)
      rlv(:,3) = cell%bmat(:,3) / kpts%nkpt3(3)
      vol = 1/6.0*sqrt(sum(rlv(:,1)*rlv(:,1))*sum(rlv(:,2)*rlv(:,2))*sum(rlv(:,3)*rlv(:,3))) 
      d = rlv(:,1) + rlv(:,3) - rlv(:,2)
      diag(1) = sum(d*d)
      d = rlv(:,2) + rlv(:,3) - rlv(:,1)
      diag(1) = sum(d*d)
      d = rlv(:,1) + rlv(:,2) + rlv(:,3)
      diag(1) = sum(d*d)
      d = rlv(:,1) + rlv(:,2) - rlv(:,3)
      diag(1) = sum(d*d)
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


END MODULE m_calc_tetra