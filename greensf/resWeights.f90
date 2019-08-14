MODULE m_resWeights

   CONTAINS

   SUBROUTINE resWeights(ikpt,kpts,neig,eig,ez,nz,weights,boundInd)

      !Weights for analytical tetrahedron method for spectral functions 

      USE m_types
      USE m_juDFT

      IMPLICIT NONE

      INTEGER,          INTENT(IN)     :: ikpt
      TYPE(t_kpts),     INTENT(IN)     :: kpts
      INTEGER,          INTENT(IN)     :: neig
      REAL,             INTENT(IN)     :: eig(neig,kpts%nkpt)
      INTEGER,          INTENT(IN)     :: nz
      COMPLEX,          INTENT(IN)     :: ez(nz)
      COMPLEX,          INTENT(INOUT)  :: weights(nz,neig)
      INTEGER,          INTENT(INOUT)  :: boundInd(neig,2)


      INTEGER itet,ib,i,j,iz,icorn,ind(4),k(4)
      REAL e(4),vol,tmp
      COMPLEX z(4),rk(4),sk(4)

      CALL juDFT_error("Not yet implemented",calledby="resWeights")


      !Here we do no truncation for now
      boundInd(:,1) = 1
      boundInd(:,1) = nz
      weights = 0.0

      DO itet = 1, kpts%ntet

         IF(ALL(kpts%ntetra(1:4,itet).NE.ikpt)) CYCLE

         IF(kpts%nkptf.NE.0) THEN
            DO i = 1, 4
               IF(kpts%ntetra(i,itet).GT.kpts%nkpt) THEN
                  k(i) = kpts%bkp(kpts%ntetra(i,itet))
               ELSE
                  k(i) = kpts%ntetra(i,itet)
               ENDIF
            ENDDO
         ENDIF

         !$OMP PARALLEL DEFAULT(none) &
         !$OMP SHARED(ikpt,itet,neig,nz,k) &
         !$OMP SHARED(kpts,eig,ez,weights) &
         !$OMP PRIVATE(ib,iz,i,j,icorn,tmp,vol) &
         !$OMP PRIVATE(ind,e,z,rk,sk)

         !$OMP DO
         DO ib = 1, neig

            e(1:4) = eig(ib,k(1:4)) 
            ind=(/1,2,3,4/)
            !Sort the energies in the tetrahedron in ascending order
            DO i = 1, 3
               DO j = i+1, 4
                  IF (e(ind(i)).GT.e(ind(j))) THEN
                     tmp = ind(i)
                     ind(i) = ind(j)
                     ind(j) = tmp
                  ENDIF
               ENDDO
            ENDDO
            !search for the corner ikpt in the sorted array
            DO i = 1, 4
               IF(kpts%ntetra(ind(i),itet).EQ.ikpt) icorn = i
            ENDDO

            vol = kpts%voltet(itet)/kpts%ntet
            DO iz = 1, nz
               z(1:4) = ez(iz)-e(ind(1:4))
               !CALL resWeightTetra(rk,sk,z)
               weights(iz,ib) = weights(iz,ib) + rk(icorn)*vol
            ENDDO

         ENDDO
         !$OMP END DO
         !$OMP END PARALLEL
      ENDDO

   END SUBROUTINE resWeights


END MODULE m_resWeights