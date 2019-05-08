MODULE m_sel_sites

   CONTAINS 

   SUBROUTINE n_neighbours(n,atoms,cell,i_nn,nn,vec_nn,at_nn)

      USE m_types

      IMPLICIT NONE

      INTEGER,       INTENT(IN)  :: n
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell),  INTENT(IN)  :: cell 
      INTEGER,       INTENT(IN)  :: i_nn  !number of nearest neighbour layers
      INTEGER,       INTENT(OUT) :: nn
      INTEGER,       INTENT(OUT) :: at_nn(27*atoms%nat)
      REAL,          INTENT(OUT) :: vec_nn(3,27*atoms%nat)
      !WORK ARRAYS
      INTEGER ind(27*atoms%nat)
      REAL    dist(27*atoms%nat)
      REAL    vec(3,27*atoms%nat)

      INTEGER index, iat, ix, iy, iz, i
      REAL tol, min_dist
      REAL shift(3),pos1(3)
      !We need to "stitch" multiple unit cells together it should be a 3 cells x 3 cells x 3 cells cube (ARE 2x2x2 enough )
      !Calculate all the distances
      index = 0 
      DO iat = 1, atoms%nat
         !cycle through the unit cells 
         DO ix = -1, 1
            DO iy = -1, 1
               DO iz = -1, 1
                  !construct the shift vector
                  shift = matmul(cell%amat(:,:),(/ix,iy,iz/))
                  index = index + 1
                  ind(index) = iat
                  !positions are in 
                  pos1 = atoms%pos(:,iat) - shift
                  dist(index) = NORM2(pos1-atoms%pos(:,n))
                  vec(:,index)= pos1-atoms%pos(:,n)
                  !Exclude the atom itself
                  IF(dist(index).EQ.0.0) dist(index) = 9.9E99
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      tol = 1E-8
      !Find the smallest distances
      nn = 0
      at_nn = 0
      vec_nn = 0.0
      DO i = 1, i_nn
         min_dist = MINVAL(dist)
         DO index = 1, 27*atoms%nat
            IF(ABS(dist(index)-min_dist).LT.tol) THEN
               nn = nn + 1
               at_nn(nn) = ind(index)
               !Exclude for the next layer
               dist(index) = 9.9E99
               vec_nn(:,nn) = vec(:,index)
            ENDIF 
         ENDDO
      ENDDO
      WRITE(*,9000) i_nn, nn 

9000  FORMAT("Coordination number with ", I3, " nearest-neighbour layers: ", I3)
   END SUBROUTINE n_neighbours

END MODULE m_sel_sites