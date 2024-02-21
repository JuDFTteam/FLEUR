!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_genMTBasis

CONTAINS

  SUBROUTINE genMTBasis(atoms,enpara,vTot,fmpi,iType,jspin,usdus,f,g,flo,hub1data,l_writeArg)
    USE m_types
    USE m_constants
    USE m_radfun
    USE m_radflo
    USE m_find_enpara
    USE m_radsra
    USE m_intgr, ONLY : intgr0
    !$  use omp_lib

    IMPLICIT NONE

    TYPE(t_atoms),  INTENT(IN)    :: atoms
    TYPE(t_enpara), INTENT(IN)    :: enpara
    TYPE(t_potden), INTENT(IN)    :: vTot
    TYPE(t_mpi),    INTENT(IN)    :: fmpi
    TYPE(t_usdus),  INTENT(INOUT) :: usdus

    INTEGER,        INTENT(IN)    :: iType
    INTEGER,        INTENT(IN)    :: jspin

    TYPE(t_hub1data), OPTIONAL, INTENT(IN)  :: hub1data

    REAL,           INTENT(INOUT) :: f(atoms%jmtd,2,0:atoms%lmaxd)
    REAL,           INTENT(INOUT) :: g(atoms%jmtd,2,0:atoms%lmaxd)
    REAL,           INTENT(INOUT) :: flo(atoms%jmtd,2,atoms%nlod)
    LOGICAL, OPTIONAL, INTENT(IN) :: l_writeArg


    INTEGER                       :: l, nodeu, noded, iLO, nLOsForL, iLOForL, nMTBasisFcts, iFunA, iFunB, iFunC
    REAL                          :: wronk, overlap, energy, e_lo, e_up, us, dus, coreNorm


    LOGICAL    :: l_write,l_hia,l_performSpinavg
    REAL       :: vrTmp(atoms%jmtd)
    REAL       :: productVals(atoms%jmtd)
    INTEGER    :: i, info, addFunA, addFunB, iState, nQuantumNumber
    
    REAL, ALLOCATABLE :: coreFun(:,:)
    REAL, ALLOCATABLE :: funValues(:,:,:)
    REAL, ALLOCATABLE :: bestRadFun(:,:)
    REAL, ALLOCATABLE :: mtOverlapMat(:,:), invMTOverlapMat(:,:)
    REAL, ALLOCATABLE :: rWorkArray(:), projVector(:)
    INTEGER, ALLOCATABLE :: ipiv(:)
    
    INTEGER :: largestCoreMainQuantumNumbers(0:3)

    l_write = .TRUE.
    IF(PRESENT(l_writeArg)) l_write = l_writeArg

    l_performSpinavg = .false.
    IF(PRESENT(hub1data)) l_performSpinavg = hub1data%l_performSpinavg

    l_write=l_write .and. fmpi%irank==0
    !$ l_write = l_write .and. omp_get_num_threads()==1
    

    IF (l_write) WRITE (oUnit,FMT=8000) iType

    DO l = 0,atoms%lmax(iType)
       !Check if the orbital is to be treated with Hubbard 1
       l_hia=.FALSE.
       DO i = atoms%n_u+1, atoms%n_u+atoms%n_hia
          IF(atoms%lda_u(i)%atomType.EQ.itype.AND.atoms%lda_u(i)%l.EQ.l) THEN
             l_hia=.TRUE.
          ENDIF
       ENDDO
       !In the case of a spin-polarized calculation with Hubbard 1 we want to treat
       !the correlated orbitals with a non-spin-polarized basis
       IF((l_hia.AND.SIZE(vTot%mt,4).GT.1 .AND. l_performSpinavg).or.atoms%l_nonpolbas(itype)) THEN
         vrTmp = (vTot%mt(:,0,iType,1) + vTot%mt(:,0,iType,2))/2.0
       ELSE
         vrTmp = vTot%mt(:,0,iType,jspin)
       ENDIF
   
       CALL radfun(l,iType,jspin,enpara%el0(l,iType,jspin),vrTmp,atoms,&
            f(1,1,l),g(1,1,l),usdus,nodeu,noded,wronk)
       IF (l_write) THEN
          WRITE (oUnit,FMT=8010) l,enpara%el0(l,iType,jspin),usdus%us(l,iType,jspin),usdus%dus(l,iType,jspin),&
               nodeu,usdus%uds(l,iType,jspin),usdus%duds(l,iType,jspin),noded,usdus%ddn(l,iType,jspin),wronk
       END IF
    END DO

    ! Generate the extra wavefunctions for the local orbitals, if there are any.
    IF (atoms%nlo(iType).GE.1) THEN
       CALL radflo(atoms,iType,jspin,enpara%ello0(1,1,jspin),vrtmp,f,g,fmpi,&
            usdus,usdus%uuilon(1,1,jspin),usdus%duilon(1,1,jspin),usdus%ulouilopn(1,1,1,jspin),flo)
    END IF
    
    largestCoreMainQuantumNumbers(:) = -1
    DO iState = 1, atoms%econf(iType)%num_core_states
       nQuantumNumber = atoms%econf(iType)%nprnc(iState)
       l = atoms%econf(iType)%get_state_l(iState)
       IF (nQuantumNumber.GT.largestCoreMainQuantumNumbers(l)) largestCoreMainQuantumNumbers(l) = nQuantumNumber
    END DO
    
    ! As quality measure: Calculate and write out overlap in MT spheres for each l for which LOs exist
    IF (l_write) THEN
    
       WRITE(oUnit,*) ''
       WRITE(oUnit,'(a,i5,a)') 'Representation of radial basis functions, core electron states in terms of the other radial basis functions for atom type ',iType,':'
       WRITE(oUnit,*) ''

       DO l = 0, atoms%lmax(iType)
          nLOsForL = 0
          DO iLO = 1, atoms%nlo(iType)
             If (atoms%llo(iLO, iType).EQ.l) nLOsForL = nLOsForL + 1
          END DO
          nMTBasisFcts = nLOsForL + 2
          ! Fill radial function array for each function
          ALLOCATE (funValues(atoms%jri(iType),2,nMTBasisFcts))
          ALLOCATE (mtOverlapMat(nMTBasisFcts,nMTBasisFcts))
          ALLOCATE (invMTOverlapMat(nMTBasisFcts-1,nMTBasisFcts-1))
          ALLOCATE (rWorkArray(4*nMTBasisFcts))
          ALLOCATE (ipiv(nMTBasisFcts))
          ALLOCATE (bestRadFun(atoms%jri(iType),2))
          ALLOCATE (coreFun(atoms%jri(iType),2))
          ALLOCATE (projVector(nMTBasisFcts-1))
          funValues(:,:,1) = f(1:atoms%jri(iType),:,l)
          funValues(:,:,2) = g(1:atoms%jri(iType),:,l)
          iLOForL = 0
          DO iLO = 1, atoms%nlo(iType)
             If (atoms%llo(iLO, iType).EQ.l) THEN
                iLOForL = iLOForL + 1
                funValues(:,:,iLOForL+2) = flo(1:atoms%jri(iType),:,iLO)
             END IF
          END DO
          ! calculate MT overlap
          DO iFunA = 1, nMTBasisFcts
             DO iFunB = 1, nMTBasisFcts
                productVals = 0.0
                DO i = 1, atoms%jri(iType)
                   productVals(i) = funValues(i,1,iFunA)*funValues(i,1,iFunB)+funValues(i,2,iFunA)*funValues(i,2,iFunB)
                END DO
                CALL intgr0(productVals,atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),overlap)
                mtOverlapMat(iFunA, iFunB) = overlap
             END DO
          END DO
          
          IF(nLOsForL.GT.0) THEN
             DO iFunC = 1, nMTBasisFcts
                DO iFunA = 1, nMTBasisFcts
                   addFunA = 0
                   IF(iFunA.EQ.iFunC) CYCLE
                   IF(iFunA.GT.iFunC) addFunA = -1
                   projVector(iFunA+addFunA) = mtOverlapMat(iFunA,iFunC)
                   DO iFunB = 1, nMTBasisFcts
                      addFunB = 0
                      IF(iFunB.EQ.iFunC) CYCLE
                      IF(iFunB.GT.iFunC) addFunB = -1
                      invMTOverlapMat(iFunA+addFunA,iFunB+addFunB) = mtOverlapMat(iFunA,iFunB)
                   END DO
                END DO
                ipiv = 0
                CALL DGETRF(nMTBasisFcts-1,nMTBasisFcts-1,invMTOverlapMat,nMTBasisFcts-1,ipiv(1:nMTBasisFcts-1),info)
                CALL DGETRI(nMTBasisFcts-1,invMTOverlapMat,nMTBasisFcts-1,ipiv(1:nMTBasisFcts-1),rWorkArray,SIZE(rWorkArray),info)
                projVector(:) = MATMUL(invMTOverlapMat(:,:),projVector(:))
                bestRadFun = 0.0
                DO iFunA = 1, nMTBasisFcts
                   addFunA = 0
                   IF(iFunA.EQ.iFunC) CYCLE
                   IF(iFunA.GT.iFunC) addFunA = -1
                   ! Calculate best representation of iFunC with all other basis functions.
                   bestRadFun(:,:) = bestRadFun(:,:) + projVector(iFunA+addFunA)*funValues(:,:,iFunA)
                END DO
                ! Subtract iFunC from its best representation
                bestRadFun(:,:) = bestRadFun(:,:) - funValues(:,:,iFunC)
                ! Calculate norm of resulting function as representation error
                productVals(:) = bestRadFun(:,1)*bestRadFun(:,1) + bestRadFun(:,2)*bestRadFun(:,2)
                overlap = 0.0
                CALL intgr0(productVals,atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),overlap)
                overlap = SQRT(overlap / mtOverlapMat(iFunC, iFunC))
                WRITE(oUnit,'(a,i2,a,i2,a,f12.8)') 'Representation error for l = ', l, ', radial basis function ', iFunC, ': ', overlap
             END DO
          END IF
          
          DEALLOCATE (projVector)
          ALLOCATE(projVector(nMTBasisFcts))

          IF (l.LE.3) THEN
             IF (largestCoreMainQuantumNumbers(l).GT.0) THEN
                coreFun = 0.0
                energy = find_enpara(.TRUE.,l,iType,jspin,largestCoreMainQuantumNumbers(l),atoms,vrTmp,e_lo,e_up,.TRUE.)
                CALL radsra(energy,l,vrTmp(:),atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),c_light(1.0), us,dus,nodeu,coreFun(:,1),coreFun(:,2))
                DO iFunA = 1, nMTBasisFcts
                   productVals = 0.0
                   DO i = 1, atoms%jri(iType)
                      productVals(i) = funValues(i,1,iFunA)*coreFun(i,1)+funValues(i,2,iFunA)*coreFun(i,2)
                   END DO
                   CALL intgr0(productVals,atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),overlap)
                   projVector(iFunA) = overlap
                END DO
                ipiv = 0
                CALL DGETRF(nMTBasisFcts,nMTBasisFcts,mtOverlapMat,nMTBasisFcts,ipiv,info)
                CALL DGETRI(nMTBasisFcts,mtOverlapMat,nMTBasisFcts,ipiv,rWorkArray,SIZE(rWorkArray),info)
                projVector(:) = MATMUL(mtOverlapMat(:,:),projVector(:))
                bestRadFun = 0.0
                DO iFunA = 1, nMTBasisFcts
                   bestRadFun(:,:) = bestRadFun(:,:) + projVector(iFunA)*funValues(:,:,iFunA)
                END DO
                bestRadFun(:,:) = bestRadFun(:,:) - coreFun(:,:)
                productVals(:) = bestRadFun(:,1)*bestRadFun(:,1) + bestRadFun(:,2)*bestRadFun(:,2)
                overlap = 0.0
                CALL intgr0(productVals,atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),overlap)
                productVals(:) = coreFun(:,1)*coreFun(:,1) + coreFun(:,2)*coreFun(:,2)
                CALL intgr0(productVals,atoms%rmsh(1,iType),atoms%dx(iType),atoms%jri(iType),coreNorm)
                overlap = SQRT(overlap / coreNorm)
                WRITE(oUnit,'(a,i2,a,i2,a,f12.8)') 'Valence basis representation error for core state with n = ', largestCoreMainQuantumNumbers(l), ', l = ', l, ': ', overlap
             END IF
          END IF
          
          DEALLOCATE (projVector)
          DEALLOCATE (coreFun)
          DEALLOCATE (bestRadFun)
          DEALLOCATE (ipiv)
          DEALLOCATE (rWorkArray)
          DEALLOCATE (invMTOverlapMat)
          DEALLOCATE (mtOverlapMat)
          DEALLOCATE (funValues)
       END DO
    END IF

8000 FORMAT (1x,/,/,' wavefunction parameters for atom type',i5,':',&
         /,t32,'radial function',t79,'energy derivative',/,t3,&
         'l',t8,'energy',t26,'value',t39,'derivative',t53,&
         'nodes',t68,'value',t81,'derivative',t95,'nodes',t107,&
         'norm',t119,'wronskian')
8010 FORMAT (i3,f10.5,2 (5x,1p,2e16.7,i5),1p,2e16.7)

  END SUBROUTINE genMTBasis

END MODULE m_genMTBasis
