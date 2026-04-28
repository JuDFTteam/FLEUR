!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_types_hyperfine
   implicit none

   PRIVATE

   PUBLIC :: t_hyperfine

   TYPE t_hyperfine

      LOGICAL              :: l_hyperfine

      REAL,    ALLOCATABLE :: dipolContribs(:,:,:)

      COMPLEX, ALLOCATABLE :: dipolSumACoeffs(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dipolSumACoeffsLO(:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dipolSumACoeffsLOLO(:,:,:,:,:)

      COMPLEX, ALLOCATABLE :: dipolSumBCoeffs(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dipolSumBCoeffsLO(:,:,:,:,:,:)
      COMPLEX, ALLOCATABLE :: dipolSumBCoeffsLOLO(:,:,:,:,:,:)

      CONTAINS

      PROCEDURE, PASS :: init => hyperfineInit
      PROCEDURE, PASS :: collect => hyperfineCollect
      PROCEDURE, PASS :: calcDipolTermCoeffs => hyperfineDipolarTermCoeffs
      PROCEDURE, PASS :: calcDipolRadFunIntegrals => hyperfineDipolarRadFunIntegrals
      PROCEDURE, PASS :: printValenceHyperfine => hyperfinePrintValence

      PROCEDURE, PASS :: calcPrintIsomerShifts => hyperfineCalcPrintIsomerShifts

   END TYPE t_hyperfine

   CONTAINS

   SUBROUTINE hyperfineInit(thisHF, input, atoms)

      USE m_types

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms

      IF (ALLOCATED(thisHF%dipolContribs)) DEALLOCATE(thisHF%dipolContribs)
      IF (ALLOCATED(thisHF%dipolSumACoeffs)) DEALLOCATE(thisHF%dipolSumACoeffs)
      IF (ALLOCATED(thisHF%dipolSumACoeffsLO)) DEALLOCATE(thisHF%dipolSumACoeffsLO)
      IF (ALLOCATED(thisHF%dipolSumACoeffsLOLO)) DEALLOCATE(thisHF%dipolSumACoeffsLOLO)
      IF (ALLOCATED(thisHF%dipolSumBCoeffs)) DEALLOCATE(thisHF%dipolSumBCoeffs)
      IF (ALLOCATED(thisHF%dipolSumBCoeffsLO)) DEALLOCATE(thisHF%dipolSumBCoeffsLO)
      IF (ALLOCATED(thisHF%dipolSumBCoeffsLOLO)) DEALLOCATE(thisHF%dipolSumBCoeffsLOLO)

      thisHF%l_hyperfine = (input%kcrel.EQ.1).AND.(input%jspins.EQ.2)

      IF (.NOT.thisHF%l_hyperfine) RETURN

      ALLOCATE(thisHF%dipolContribs(-1:3,atoms%ntype,2)) ! Last index is the spin index
      ALLOCATE(thisHF%dipolSumACoeffs(0:1, 0:1, 0:(atoms%lmaxd+1)**2-1, atoms%ntype, input%jspins))
      ALLOCATE(thisHF%dipolSumACoeffsLO(0:1, atoms%nlod, -atoms%lmaxd:atoms%lmaxd, atoms%ntype, input%jspins))
      ALLOCATE(thisHF%dipolSumACoeffsLOLO(atoms%nlod, atoms%nlod, -atoms%lmaxd:atoms%lmaxd, atoms%ntype, input%jspins))
      ALLOCATE(thisHF%dipolSumBCoeffs(0:1, 0:1, 0:(atoms%lmaxd+1)**2-1, 0:(atoms%lmaxd+1)**2-1, atoms%ntype, input%jspins))
      ALLOCATE(thisHF%dipolSumBCoeffsLO(0:1, atoms%nlod, -atoms%lmaxd:atoms%lmaxd, 0:(atoms%lmaxd+1)**2-1, atoms%ntype, input%jspins))
      ALLOCATE(thisHF%dipolSumBCoeffsLOLO(atoms%nlod, atoms%nlod, -atoms%lmaxd:atoms%lmaxd, -atoms%lmaxd:atoms%lmaxd, atoms%ntype, input%jspins))


      thisHF%dipolContribs = 0.0
      thisHF%dipolSumACoeffs = 0.0
      thisHF%dipolSumACoeffsLO = 0.0
      thisHF%dipolSumACoeffsLOLO = 0.0
      thisHF%dipolSumBCoeffs = 0.0
      thisHF%dipolSumBCoeffsLO = 0.0
      thisHF%dipolSumBCoeffsLOLO = 0.0
   END SUBROUTINE hyperfineInit

   SUBROUTINE hyperfineCollect(thisHF, fmpi, atoms, iSpin)

#ifdef CPP_MPI
      USE mpi
#endif
      USE m_juDFT
      USE m_types

      CLASS(t_hyperfine), INTENT(INOUT) :: thisHF
      TYPE(t_mpi),        INTENT(IN)    :: fmpi
      TYPE(t_atoms),      INTENT(IN)    :: atoms

      INTEGER,            INTENT(IN)    :: iSpin


      INTEGER :: length, ierr
      COMPLEX, ALLOCATABLE :: cBuffer(:)

#ifdef CPP_MPI
      IF (.NOT.thisHF%l_hyperfine) RETURN

      length = 2 * 2 * (atoms%lmaxd+1)**2 * atoms%ntype
      ALLOCATE(cBuffer(length))
      CALL MPI_ALLREDUCE(thisHF%dipolSumACoeffs(0:,0:,0:,:,iSpin),cBuffer,length,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL zcopy(length, cBuffer, 1, thisHF%dipolSumACoeffs(0:,0:,0:,:,iSpin), 1)
      DEALLOCATE(cBuffer)

      length = 2 * atoms%nlod * (2*atoms%lmaxd+1) * atoms%ntype
      ALLOCATE(cBuffer(length))
      CALL MPI_ALLREDUCE(thisHF%dipolSumACoeffsLO(0:1,:,-atoms%lmaxd:atoms%lmaxd,:,iSpin),cBuffer,length,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL zcopy(length, cBuffer, 1, thisHF%dipolSumACoeffsLO(0:1,:,-atoms%lmaxd:atoms%lmaxd,:,iSpin), 1)
      DEALLOCATE(cBuffer)

      length = atoms%nlod * atoms%nlod * (2*atoms%lmaxd+1) * atoms%ntype
      ALLOCATE(cBuffer(length))
      CALL MPI_ALLREDUCE(thisHF%dipolSumACoeffsLOLO(:,:,-atoms%lmaxd:atoms%lmaxd,:,iSpin),cBuffer,length,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL zcopy(length, cBuffer, 1, thisHF%dipolSumACoeffsLOLO(:,:,-atoms%lmaxd:atoms%lmaxd,:,iSpin), 1)
      DEALLOCATE(cBuffer)

      length = 2 * 2 * (atoms%lmaxd+1)**2 * (atoms%lmaxd+1)**2 * atoms%ntype
      ALLOCATE(cBuffer(length))
      CALL MPI_ALLREDUCE(thisHF%dipolSumBCoeffs(0:,0:,0:,0:,:,iSpin),cBuffer,length,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL zcopy(length, cBuffer, 1, thisHF%dipolSumBCoeffs(0:,0:,0:,0:,:,iSpin), 1)
      DEALLOCATE(cBuffer)

      length = 2 * atoms%nlod * (2*atoms%lmaxd+1) * (atoms%lmaxd+1)**2 * atoms%ntype
      ALLOCATE(cBuffer(length))
      CALL MPI_ALLREDUCE(thisHF%dipolSumBCoeffsLO(0:1,:,-atoms%lmaxd:atoms%lmaxd,0:,:,iSpin),cBuffer,length,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL zcopy(length, cBuffer, 1, thisHF%dipolSumBCoeffsLO(0:1,:,-atoms%lmaxd:atoms%lmaxd,0:,:,iSpin), 1)
      DEALLOCATE(cBuffer)

      length = atoms%nlod * atoms%nlod * (2*atoms%lmaxd+1) * (2*atoms%lmaxd+1) * atoms%ntype
      ALLOCATE(cBuffer(length))
      CALL MPI_ALLREDUCE(thisHF%dipolSumBCoeffsLOLO(:,:,-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,:,iSpin),cBuffer,length,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
      CALL zcopy(length, cBuffer, 1, thisHF%dipolSumBCoeffsLOLO(:,:,-atoms%lmaxd:atoms%lmaxd,-atoms%lmaxd:atoms%lmaxd,:,iSpin), 1)
      DEALLOCATE(cBuffer)
#endif

   END SUBROUTINE hyperfineCollect

   SUBROUTINE hyperfineDipolarTermCoeffs(thisHF, atoms, eigVecCoeffs, numOccBands, weights, iSpin)

      USE m_constants
      USE m_types
      USE m_gaunt, ONLY: gaunt1

      CLASS(t_hyperfine),   INTENT(INOUT) :: thisHF
      TYPE(t_atoms),        INTENT(IN)    :: atoms
      TYPE(t_eigVecCoeffs), INTENT(IN)    :: eigVecCoeffs

      INTEGER,              INTENT(IN)    :: numOccBands, iSpin

      REAL,                 INTENT(IN)    :: weights(:) ! weights for the states (occupation times k-weight)

      INTEGER :: l, m, lm, lp, mp, lmp, iAtom, iType, iBand, iRadFun, jRadFun, iLO, jLO, matSize
      REAL    :: factorA

      REAL    :: dirOpZMat(0:(atoms%lmaxd+1)**2-1,0:(atoms%lmaxd+1)**2-1)
      REAL    :: dirOpZMatSquare(0:(atoms%lmaxd+1)**2-1,0:(atoms%lmaxd+1)**2-1)

      !
      ! Delta E_{hf}^{dipolar} = mu_{B} * Sum_{iBand} < Psi_{iBand} | 1/(r^3) [ sigma*mu_{N} - 3 * (sigma*hat{r}) * (mu_{N}*hat{r}) ] | Psi_{iBand} >
      !
      ! sigma and mu_{N} are assumed to point along z direction. The magnitude of mu_{N} is not considered. The radial integration is not part of this subroutine.
      ! The spin loop is outside this subroutine. Only the prefactors to the two summands in the dipolar term for the respective spin channel are calculated here.
      !

      ! Plan for this subroutine:
      !
      ! 1. Calculate the prefactors to the combinations of the radial functions for the 1st summand of the dipolar term
      ! 2. Calculate a matrix of size lm_max x lm_max that covers the application of the z component of the direction operator (\hat{r}) on spherical harmonics.
      ! 3. Square the matrix and multiply the product by -3.0.
      ! 4. Calculate the prefactors to the combinations of the radial functions for the 2nd summand of the dipolar term



      ! 0. General initializations

      ! 1. Calculate the prefactors to the combinations of the radial functions for the 1st summand of the dipolar term

      IF (.NOT.thisHF%l_hyperfine) RETURN

      DO iAtom = 1,atoms%nat
         iType = atoms%itype(iAtom)

         ! 1.1 Radial functions u and u-dot
         DO l = 0,atoms%lmax(iType)
            DO m = -l, l
               lm = l*(l+1) + m
               DO iBand = 1, numOccBands
                  DO jRadFun = 0, 1
                     DO iRadFun = 0, 1
                        thisHF%dipolSumACoeffs(iRadFun, jRadFun, lm, iType, iSpin) = thisHF%dipolSumACoeffs(iRadFun, jRadFun, lm, iType, iSpin) + &
                           weights(iBand) * eigVecCoeffs%abcof(iBand,lm,iRadFun,iAtom,iSpin) * CONJG(eigVecCoeffs%abcof(iBand,lm,jRadFun,iAtom,iSpin))
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! 1.2 Local orbial radial functions
         DO iLO = 1, atoms%nlo(iType)
            l = atoms%llo(iLO,iType)

            ! 1.2.1 LO - LAPW
            DO m = -l, l
               lm = l* (l+1) + m
               DO iBand = 1, numOccBands
                  DO iRadFun = 0, 1
                     ! Note: LO-LAPW and LAPW-LO contributions are combined in the same indices
                     thisHF%dipolSumACoeffsLO(iRadFun, iLO, m, iType, iSpin) = thisHF%dipolSumACoeffsLO(iRadFun, iLO, m, iType, iSpin) + &
                        weights(iBand) * (eigVecCoeffs%abcof(iBand,lm,iRadFun,iAtom,iSpin) * CONJG(eigVecCoeffs%ccof(m,iBand,iLO,iAtom,iSpin)) + &
                                         eigVecCoeffs%ccof(m,iBand,iLO,iAtom,iSpin) * CONJG(eigVecCoeffs%abcof(iBand,lm,iRadFun,iAtom,iSpin)))
                  END DO
               END DO
            END DO

            ! 1.2.2 LO - LO
            DO jLO = 1, atoms%nlo(iType)
               IF (atoms%llo(jLO,iType).EQ.l) THEN
                  DO m = -l, l
                     DO iBand = 1, numOccBands
                        thisHF%dipolSumACoeffsLOLO(jLO, iLO, m, iType, iSpin) = thisHF%dipolSumACoeffsLOLO(jLO, iLO, m, iType, iSpin) + &
                           weights(iBand) * eigVecCoeffs%ccof(m,iBand,jLO,iAtom,iSpin) * CONJG(eigVecCoeffs%ccof(m,iBand,iLO,iAtom,iSpin))
                     END DO
                  END DO
               END IF
            END DO
         END DO
      END DO

      ! 2. Calculate a matrix of size lm_max x lm_max that covers the application of the z component of the direction operator (\hat{r}) on spherical harmonics.
      matSize = (atoms%lmaxd+1)**2
      dirOpZMat = 0.0
      dirOpZMatSquare = 0.0
      ! The z component of the direction operator is 2 * sqrt(Pi/3) * Y_{1,0}
      factorA = 2.0 * SQRT(pi_const/3.0)
      DO l = 0, atoms%lmaxd
         DO m = -l, l
            lm = l*(l+1) + m
            DO lp = 0, atoms%lmaxd
               DO mp = -lp, lp
                  lmp = lp*(lp+1) + mp
                  dirOpZMat(lmp,lm) = dirOpZMat(lmp,lm) + factorA * gaunt1(lp,1,l,mp,0,m,atoms%lmaxd)
               END DO
            END DO
         END DO
      END DO

      ! 3. Square the matrix and multiply the product by -3.0.
      ! Note the factor -3.0 for alpha in dgemm. This is the factor -3.0 also present in the dipolar term. It is not from the direction operator.
      CALL dgemm("N","N",matSize,matSize,matSize,-3.0,dirOpZMat(0:,0:),matSize,dirOpZMat(0:,0:),matSize,0.0,dirOpZMatSquare(0:,0:),matSize)

      ! 4. Calculate the prefactors to the combinations of the radial functions for the 2nd summand of the dipolar term
      DO iAtom = 1,atoms%nat
         iType = atoms%itype(iAtom)
         ! 4.1 Radial functions u and u-dot
         DO  l = 0, atoms%lmax(iType)
            DO m = -l, l
               lm = l*(l+1) + m
               DO lp = 0, atoms%lmax(iType)
                  DO mp = -lp, lp
                     lmp = lp*(lp+1) + mp
                     DO iBand = 1, numOccBands
                        DO jRadFun = 0, 1
                           DO iRadFun = 0, 1
                              thisHF%dipolSumBCoeffs(iRadFun, jRadFun, lm, lmp, iType, iSpin) = thisHF%dipolSumBCoeffs(iRadFun, jRadFun, lm, lmp, iType, iSpin) + &
                                 weights(iBand) * eigVecCoeffs%abcof(iBand,lm,iRadFun,iAtom,iSpin) * CONJG(eigVecCoeffs%abcof(iBand,lmp,jRadFun,iAtom,iSpin)) * dirOpZMat(lmp,lm)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! 4.2 Local orbial radial functions

         DO iLO = 1, atoms%nlo(iType)
            l = atoms%llo(iLO,iType)
            DO m = -l, l
               lm = l* (l+1) + m

               ! 1.2.1 LO - LAPW
               DO lp = 0, atoms%lmax(iType)
                  DO mp = -lp, lp
                     lmp = lp*(lp+1) + mp
                     DO iBand = 1, numOccBands
                        DO iRadFun = 0, 1
                           thisHF%dipolSumBCoeffsLO(iRadFun, iLO, m, lmp, iType, iSpin) = thisHF%dipolSumBCoeffsLO(iRadFun, iLO, m, lmp, iType, iSpin) + &
                              weights(iBand) * (eigVecCoeffs%abcof(iBand,lmp,iRadFun,iAtom,iSpin) * CONJG(eigVecCoeffs%ccof(m,iBand,iLO,iAtom,iSpin) * dirOpZMat(lm,lmp) ) + &
                                               eigVecCoeffs%ccof(m,iBand,iLO,iAtom,iSpin) * CONJG(eigVecCoeffs%abcof(iBand,lmp,iRadFun,iAtom,iSpin) * dirOpZMat(lmp,lm)))
                        END DO
                     END DO
                  END DO
               END DO

               ! 1.2.2 LO - LO
               DO jLO = 1, atoms%nlo(iType)
                  IF (atoms%llo(jLO,iType).EQ.l) THEN
                     DO mp = -l, l
                        DO iBand = 1, numOccBands
                           thisHF%dipolSumBCoeffsLOLO(jLO, iLO, m, mp, iType, iSpin) = thisHF%dipolSumBCoeffsLOLO(jLO, iLO, m, mp, iType, iSpin) + &
                              weights(iBand) * eigVecCoeffs%ccof(mp,iBand,jLO,iAtom,iSpin) * CONJG(eigVecCoeffs%ccof(m,iBand,iLO,iAtom,iSpin)) * dirOpZMat(lm,lmp)
                        END DO
                     END DO
                  END IF
               END DO

            END DO
         END DO
      END DO

   END SUBROUTINE hyperfineDipolarTermCoeffs

   SUBROUTINE hyperfineDipolarRadFunIntegrals(thisHF, atoms, enpara, vTot, fmpi, usdus, jsp_start, jsp_end)

      USE m_types
      USE m_intgr
      USE m_genMTBasis

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_enpara),      INTENT(IN)    :: enpara
      TYPE(t_potden),      INTENT(IN)    :: vTot
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      TYPE(t_usdus),       INTENT(INOUT) :: usdus

      INTEGER,             INTENT(IN)    :: jsp_start, jsp_end

      INTEGER :: iType, iAtom, iSpin, l, lp, m, mp, lm, lmp, iRadFun, jRadFun, iRad, iLO, jLO
      REAL    :: integral

      REAL    :: radFun(atoms%jmtd,2,0:atoms%lmaxd,0:1)
      REAL    :: flo(atoms%jmtd,2,atoms%nlod)
      REAL    :: productRadFun(atoms%jmtd,-1:3)
      REAL    :: contribution(-1:3)

      IF (.NOT.thisHF%l_hyperfine) RETURN

      DO iSpin = jsp_start, jsp_end
         DO iType = 1,atoms%ntype

            radFun = 0.0
            flo = 0.0
            contribution = 0.0

            CALL genMTBasis(atoms,enpara,vTot,fmpi,iType,iSpin,usdus,radFun(:,:,:,0),radFun(:,:,:,1),flo,l_writeArg=.FALSE.)

            iAtom = atoms%firstAtom(itype)

            ! A terms:

            ! Radial functions u and u-dot

            DO l = 0,atoms%lmax(itype)
               DO jRadFun = 0, 1
                  DO iRadFun = 0, 1
                     productRadFun = 0.0
                     DO m = -l, l
                        lm = l*(l+1) + m
                        DO iRad = 1, atoms%jri(iType)
                           productRadFun(iRad,-1) = productRadFun(iRad,-1) + thisHF%dipolSumACoeffs(iRadFun, jRadFun, lm, iType, iSpin) * (radFun(iRad,1,l,iRadFun) * radFun(iRad,1,l,jRadFun) + radFun(iRad,2,l,iRadFun) * radFun(iRad,2,l,jRadFun)) / (atoms%rmsh(iRad,iType)**3.0)
                           IF(l.LE.3) THEN
                              productRadFun(iRad,l) = productRadFun(iRad,l) + thisHF%dipolSumACoeffs(iRadFun, jRadFun, lm, iType, iSpin) * (radFun(iRad,1,l,iRadFun) * radFun(iRad,1,l,jRadFun) + radFun(iRad,2,l,iRadFun) * radFun(iRad,2,l,jRadFun)) / (atoms%rmsh(iRad,iType)**3.0)
                           END IF
                        END DO
                     END DO
                     DO lp = -1, 3
                        CALL intgr3(productRadFun(:,lp),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),integral)
                        contribution(lp) = contribution(lp) + integral
                     END DO
                  END DO
               END DO
            END DO


            ! Local orbial radial functions

            DO iLO = 1, atoms%nlo(iType)
               l = atoms%llo(iLO,iType)

               ! LO - LAPW

               DO iRadFun = 0, 1
                  DO m = -l, l
                     lm = l* (l+1) + m
                     productRadFun = 0.0
                     DO iRad = 1, atoms%jri(iType)
                        productRadFun(iRad,-1) = productRadFun(iRad,-1) + thisHF%dipolSumACoeffsLO(iRadFun, iLO, m, iType, iSpin) * (radFun(iRad,1,l,iRadFun) * flo(iRad,1,iLO) + radFun(iRad,2,l,iRadFun) * flo(iRad,2,iLO)) / (atoms%rmsh(iRad,iType)**3.0)
                        IF(l.LE.3) THEN
                           productRadFun(iRad,l) = productRadFun(iRad,l) + thisHF%dipolSumACoeffsLO(iRadFun, iLO, m, iType, iSpin) * (radFun(iRad,1,l,iRadFun) * flo(iRad,1,iLO) + radFun(iRad,2,l,iRadFun) * flo(iRad,2,iLO)) / (atoms%rmsh(iRad,iType)**3.0)
                        END IF
                     END DO
                     DO lp = -1, 3
                        CALL intgr3(productRadFun(:,lp),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),integral)
                        contribution(lp) = contribution(lp) + integral
                     END DO
                  END DO
               END DO

               ! LO - LO

               DO jLO = 1,atoms%nlo(itype)
                  IF (atoms%llo(jLO,itype)==l) THEN
                     DO m = -l, l
                        productRadFun = 0.0
                        DO iRad = 1, atoms%jri(iType)
                           productRadFun(iRad,-1) = productRadFun(iRad,-1) + thisHF%dipolSumACoeffsLOLO(jLO, iLO, m, iType, iSpin) * (flo(iRad,1,jLO) * flo(iRad,1,iLO) + flo(iRad,2,jLO) * flo(iRad,2,iLO)) / (atoms%rmsh(iRad,iType)**3.0)
                           IF(l.LE.3) THEN
                              productRadFun(iRad,l) = productRadFun(iRad,l) + thisHF%dipolSumACoeffsLOLO(jLO, iLO, m, iType, iSpin) * (flo(iRad,1,jLO) * flo(iRad,1,iLO) + flo(iRad,2,jLO) * flo(iRad,2,iLO)) / (atoms%rmsh(iRad,iType)**3.0)
                           END IF
                        END DO
                        DO lp = -1, 3
                           CALL intgr3(productRadFun(:,lp),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),integral)
                           contribution(lp) = contribution(lp) + integral
                        END DO
                     END DO
                  END IF
               END DO

            END DO

            ! B terms:

            ! Radial functions u and u-dot

            DO lp = 0,atoms%lmax(itype)
               DO jRadFun = 0, 1
                  DO mp = -lp, lp
                     lmp = lp*(lp+1) + mp
                     DO l = 0,atoms%lmax(itype)
                        DO iRadFun = 0, 1
                           DO m = -l, l
                              lm = l*(l+1) + m
                              productRadFun = 0.0
                              DO iRad = 1, atoms%jri(iType)
                                 productRadFun(iRad,-1) = productRadFun(iRad,-1) + thisHF%dipolSumBCoeffs(iRadFun, jRadFun, lm, lmp, iType, iSpin) * (radFun(iRad,1,l,iRadFun) * radFun(iRad,1,lp,jRadFun) + radFun(iRad,2,l,iRadFun) * radFun(iRad,2,lp,jRadFun)) / (atoms%rmsh(iRad,iType)**3.0)
                                 ! For the time being we don't l-resolve the contributions here. It is not clear what to do there (l and lp are involved).
                              END DO
                              CALL intgr3(productRadFun(:,-1),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),integral)
                              contribution(-1) = contribution(-1) + integral
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            ! Local orbial radial functions

            DO iLO = 1, atoms%nlo(iType)
               l = atoms%llo(iLO,iType)
               DO m = -l, l

                  ! LO - LAPW

                  DO lp = 0,atoms%lmax(itype)
                     DO mp = -lp, lp
                        lmp = lp*(lp+1) + mp
                        DO jRadFun = 0, 1
                           productRadFun = 0.0
                           DO iRad = 1, atoms%jri(iType)
                              productRadFun(iRad,-1) = productRadFun(iRad,-1) + thisHF%dipolSumBCoeffsLO(iRadFun, iLO, m, lmp, iType, iSpin) * (radFun(iRad,1,lp,iRadFun) * flo(iRad,1,iLO) + radFun(iRad,2,lp,iRadFun) * flo(iRad,2,iLO)) / (atoms%rmsh(iRad,iType)**3.0)
                              ! For the time being we don't l-resolve the contributions here. It is not clear what to do there (l and lp are involved).
                           END DO
                           CALL intgr3(productRadFun(:,-1),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),integral)
                           contribution(-1) = contribution(-1) + integral
                        END DO
                     END DO
                  END DO

                  ! LO - LO

                  DO jLO = 1, atoms%nlo(iType)
                     lp = atoms%llo(jLO,iType)
                     DO mp = -lp, lp
                        productRadFun = 0.0
                        DO iRad = 1, atoms%jri(iType)
                           productRadFun(iRad,-1) = productRadFun(iRad,-1) + thisHF%dipolSumBCoeffsLOLO(jLO, iLO, m, mp, iType, iSpin) * (flo(iRad,1,jLO) * flo(iRad,1,iLO) + flo(iRad,2,jLO) * flo(iRad,2,iLO)) / (atoms%rmsh(iRad,iType)**3.0)
                           ! For the time being we don't l-resolve the contributions here. It is not clear what to do there (l and lp are involved).
                        END DO
                        CALL intgr3(productRadFun(:,-1),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),integral)
                        contribution(-1) = contribution(-1) + integral
                     END DO
                  END DO
               END DO
            END DO

            thisHF%dipolContribs(-1:3,iType,iSpin) = contribution(-1:3)

         END DO

      END DO

   END SUBROUTINE hyperfineDipolarRadFunIntegrals

   SUBROUTINE hyperfinePrintValence(thisHF, input, atoms, fmpi, moments)
      USE m_constants
      USE m_types

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      TYPE(t_moments),     INTENT(INOUT) :: moments

      INTEGER :: iType

      REAL    :: a0, e0, cautog, bohrMagInCGS
      REAL    :: hyperfineResults(-1:3), hyperfineResultsTotal(-1:3)

      IF ((fmpi%irank.EQ.0).AND.thisHF%l_hyperfine) THEN
         ! Print out valence contributions to the hyperfine field
         a0 = bohr_to_angstrom_const * 1.0e-8
         e0 = 1.6021892e-19 * 2.997930e+09
         cautog = e0 / (a0*a0)
         bohrMagInCGS = 1.0/(2.0*c_light(1.0))
         WRITE(oUnit,*) ''
         WRITE(ounit,*) ' Hyperfine field valence contributions in kG '
         WRITE(ounit,*) ' ========================================================== '
         WRITE(ounit,*) ' atom type                          contribution'
         WRITE(ounit,*) '                total         s           p           d           f'
         DO iType = 1, atoms%ntype
            moments%hypFineContribs(:,iType,1,1) = moments%hypFineContribs(:,iType,1,1) - moments%hypFineContribs(:,iType,2,1)
            hyperfineResults(:) = moments%hypFineContribs(:,iType,1,1) * cautog * 0.001 * sfp_const * bohrMagInCGS * 8.0 * pi_const / 3.0
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResults(-1:3), 'contact term'
            hyperfineResultsTotal(:) = hyperfineResults(:)
            ! dipolar term is not yet correct
!            moments%hypFineContribs(:,iType,1,2) = thisHF%dipolContribs(:,iType,1) - thisHF%dipolContribs(:,iType,2)
!            hyperfineResults(:) = moments%hypFineContribs(:,iType,1,2) * cautog * 0.001 * sfp_const * bohrMagInCGS
!            WRITE(oUnit,'(i7,5x,f12.5,53x,a)') iType, hyperfineResults(-1), 'dipolar term'
!            hyperfineResultsTotal(-1) = hyperfineResultsTotal(-1) + hyperfineResults(-1)
            moments%hypFineContribs(:,iType,1,3) = moments%hypFineContribs(:,iType,1,3) + moments%hypFineContribs(:,iType,2,3)
            hyperfineResults(:) = moments%hypFineContribs(:,iType,1,3) * cautog * 0.001 / c_light(1.0)
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResults(-1:3), 'orbital term'
            hyperfineResultsTotal(:) = hyperfineResultsTotal(:) + hyperfineResults(:)
            WRITE(oUnit,'(i7,5x,5f12.5,5x,a)') iType, hyperfineResultsTotal(-1:3), 'all terms'
         END DO
         WRITE(ounit,*) ' ========================================================== '
      END IF

   END SUBROUTINE hyperfinePrintValence

   SUBROUTINE hyperfineCalcPrintIsomerShifts(thisHF, input, atoms, fmpi, den)

      USE m_constants
      USE m_types
      USE m_intgr, ONLY : intgr2

      CLASS(t_hyperfine),  INTENT(INOUT) :: thisHF
      TYPE(t_input),       INTENT(IN)    :: input
      TYPE(t_atoms),       INTENT(IN)    :: atoms
      TYPE(t_mpi),         INTENT(IN)    :: fmpi
      TYPE(t_potden),      INTENT(IN)    :: den

      INTEGER :: iType, iRad

      REAL    :: nucRad, alpha, contactDenVal, smallRadDenVal, largeRadDenVal, averageDen
      REAL    :: indefInteg(atoms%jmtd), sphrDen(atoms%jmtd)

      LOGICAL :: isSmaller

      IF ((fmpi%irank.EQ.0) ) THEN !.AND.thisHF%l_hyperfine) THEN
         WRITE(oUnit,*) ''
         WRITE(ounit,*) ' Isomer shift output (contact charge density) '
         WRITE(ounit,*) ' ============================================================== '
         WRITE(ounit,*) ' atom type     radius of nucleus (Bohr)   smallest radial mesh point (Bohr)   charge density value (e/Bohr^3) '
         DO iType = 1, atoms%ntype
            indefInteg = 0.0
            nucRad = r0_const*(atomicMasses_const(atoms%nz(iType))**(1.0/3.0))
            IF (atoms%rmsh(1,iType).LT.nucRad) THEN
               iRad = 0
               isSmaller = .TRUE.
               DO WHILE (isSmaller)
                  iRad = iRad + 1
                  IF (atoms%rmsh(iRad,iType).GE.nucRad) isSmaller = .FALSE.
               END DO
               ! Poor man's solution: linear interpolation. :)
               alpha = (nucRad - atoms%rmsh(iRad-1,iType)) / (atoms%rmsh(iRad,iType) - atoms%rmsh(iRad-1,iType))
               IF (input%jspins.EQ.1) THEN
                  smallRadDenVal = den%mt(iRad-1,0,iType,1) / ((atoms%rmsh(iRad-1,iType)**2.0)*sfp_const)
                  largeRadDenVal = den%mt(iRad,0,iType,1) / ((atoms%rmsh(iRad,iType)**2.0)*sfp_const)
               ELSE
                  smallRadDenVal = (den%mt(iRad-1,0,iType,1) + den%mt(iRad-1,0,iType,2)) / ((atoms%rmsh(iRad-1,iType)**2.0)*sfp_const)
                  largeRadDenVal = (den%mt(iRad,0,iType,1) + den%mt(iRad,0,iType,2)) / ((atoms%rmsh(iRad,iType)**2.0)*sfp_const)
               END IF
               contactDenVal = (1.0-alpha)*smallRadDenVal + alpha*largeRadDenVal

!               WRITE(oUnit,*) alpha, smallRadDenVal, largeRadDenVal, atoms%rmsh(iRad-1,iType), atoms%rmsh(iRad,iType)

               sphrDen = 0.0
               IF (input%jspins.EQ.1) THEN
                  sphrDen(:) = den%mt(:,0,iType,1)
               ELSE
                  sphrDen(:) = den%mt(:,0,iType,1) + den%mt(:,0,iType,2)
               END IF
               indefInteg = 0.0
               CALL intgr2(sphrDen(:),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),indefInteg)

               ! Again linear interpolation. :)
               smallRadDenVal = indefInteg(iRad-1)*sfp_const / ((4.0/3.0)*pi_const*(atoms%rmsh(iRad-1,iType)**3.0))
               largeRadDenVal = indefInteg(iRad)*sfp_const / ((4.0/3.0)*pi_const*(atoms%rmsh(iRad,iType)**3.0))
               averageDen = (1.0-alpha)*smallRadDenVal + alpha*largeRadDenVal
               averageDen = averageDen

               WRITE(oUnit,'(i7,10x,f15.8,15x,f15.8,20x,f17.8,5x,a)') iType, nucRad, atoms%rmsh(1,iType), contactDenVal, 'density at radius of nuncleus'
               WRITE(oUnit,'(i7,10x,15x,50x,f17.8,5x,a)') iType, averageDen, 'average density over nuncleus'
            ELSE
               WRITE(oUnit,'(i7,10x,f15.8,20x,a)') iType, nucRad, 'smallest radial mesh point outside radius of nucleus'
            END IF
         END DO
         WRITE(ounit,*) ' ============================================================== '
      END IF


   END SUBROUTINE hyperfineCalcPrintIsomerShifts

END MODULE m_types_hyperfine
