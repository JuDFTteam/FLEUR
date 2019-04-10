!--------------------------------------------------------------------------------
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_bfgs_b2

CONTAINS

! The subroutine bfgs_b2 is an implementation of a BFGS algorithm with
! box constraints. It is motivated by the paper:
! Byrd et al., "A Limited Memory Algorithm for Bound Constrained Optimization", 
! SIAM Journal on Scientific Computing 16, 1190 (1995)
! However, it is substantially different from that paper.
!
!                            GM'2019

! Note: lastGradient and lastParameters have to be 0.0 if this is the first call to the subroutine and these
!       vectors don't exist yet. Also historyLength has to be 0 in that case.

SUBROUTINE bfgs_b2(vecLength,gradient,lastGradient,minConstraints,maxConstraints,enabledConstraints,parameters,&
                   lastParameters,equalityLinCombi,equalityCriterion,maxHistoryLength,paramCorrections,&
                   gradientCorrections,iStep,mixParam,l_converged,convCrit)

   USE m_juDFT

   IMPLICIT NONE

   INTEGER,        INTENT(IN)    :: vecLength
   INTEGER,        INTENT(IN)    :: maxHistoryLength
   INTEGER,        INTENT(INOUT) :: iStep ! The index of the bfgs_b iteration
   REAL,           INTENT(IN)    :: mixParam
   REAL,           INTENT(IN)    :: equalityCriterion
   REAL,           INTENT(IN)    :: gradient(vecLength)
   REAL,           INTENT(INOUT) :: lastGradient(vecLength)
   REAL,           INTENT(IN)    :: minConstraints(vecLength)
   REAL,           INTENT(IN)    :: maxConstraints(vecLength)
   LOGICAL,        INTENT(IN)    :: enabledConstraints(vecLength)
   REAL,           INTENT(INOUT) :: parameters(vecLength)
   REAL,           INTENT(INOUT) :: lastParameters(vecLength)
   REAL,           INTENT(IN)    :: equalityLinCombi(vecLength)
   REAL,           INTENT(INOUT) :: paramCorrections(vecLength,maxHistoryLength)
   REAL,           INTENT(INOUT) :: gradientCorrections(vecLength,maxHistoryLength)
   LOGICAL,        INTENT(OUT)   :: l_converged
   REAL,           INTENT(IN)    :: convCrit

   INTEGER                :: i, j, stepIndex, numFreeParams, numFixedParams
   INTEGER                :: historyLength
   REAL, PARAMETER        :: eps = 1.0e-10
   REAL, PARAMETER        :: largeEps = 1.0e-6
   REAL                   :: temp, norm, subspaceEqualityCrit
   REAL                   :: ddot

   LOGICAL                :: converged

   REAL                   :: lastGradCorrection(vecLength)
   REAL                   :: lastParamCorrection(vecLength)
   REAL                   :: prelimNextPoint(vecLength)
   REAL                   :: prelimGradient(vecLength)
   REAL                   :: tempVecA(vecLength)
   REAL                   :: subspaceGrad(vecLength)
   REAL                   :: subspaceNextPoint(vecLength)
   REAL                   :: subspaceLinCombi(vecLength)
   REAL                   :: diffPoint(vecLength)
   REAL                   :: diffGrad(vecLength)

   LOGICAL                :: freeParams(vecLength)

   REAL, ALLOCATABLE      :: hessMat(:,:), invHessMat(:,:)
   REAL, ALLOCATABLE      :: tempMatB(:,:), tempMatC(:,:)
   REAL, ALLOCATABLE      :: subspaceBasis(:,:)
   REAL, ALLOCATABLE      :: subspaceInvHessMat(:,:)

   WRITE(*,*) 'bfgs_b2 - Point A'

   l_converged = .FALSE.

   ! 0. Add last data pair to history iff it is not obviously destructive for the positive definiteness of the Hessian matrix

   WRITE(2000,*) '=============================================='
   WRITE(2000,'(a,8f15.10)') 'parameters(:):     ', parameters(:)
   WRITE(2000,'(a,8f15.10)') 'lastParameters(:): ', lastParameters(:)
   WRITE(2000,'(a,8f15.10)') 'gradient(:):       ', gradient(:)
   WRITE(2000,'(a,8f15.10)') 'lastGradient(:):   ', lastGradient(:)
   IF(ANY(lastGradient(:).NE.0.0).OR.ANY(lastParameters(:).NE.0.0)) THEN
   WRITE(*,*) 'bfgs_b2 - Point A - 1'
      lastGradCorrection(:) = gradient(:) - lastGradient(:)
      lastParamCorrection(:) = parameters(:) - lastParameters(:)
      temp = ddot(vecLength,lastParamCorrection(:),1,lastGradCorrection(:),1)
      IF (temp.GT.eps) THEN
         ! Add data pair to history
         istep = istep + 1
         stepIndex = MOD(iStep - 1,maxHistoryLength) + 1
         paramCorrections(:,stepIndex) = lastParamCorrection(:)
         gradientCorrections(:,stepIndex) = lastGradCorrection(:)
      END IF

      ! Check for convergence and return if converged
      IF(ALL(ABS(lastGradCorrection(:)).LT.convCrit).AND.ALL(ABS(lastParamCorrection(:)).LT.convCrit)) THEN
         l_converged = .TRUE.
         RETURN
      END IF
   END IF

   ! 1. Construct Hessian matrix and its inverse.

   historyLength = MIN(istep,maxHistoryLength)

   ALLOCATE(hessMat(vecLength,vecLength))
   ALLOCATE(invHessMat(vecLength,vecLength))
   ALLOCATE(tempMatB(vecLength,vecLength))

   hessMat(:,:) = 0.0
   invHessMat(:,:) = 0.0
   tempMatB(:,:) = 0.0

   DO i = 1, vecLength
      hessMat(i,i) = 1.0
      invHessMat(i,i) = 1.0
   END DO

   IF(historyLength.EQ.0) THEN
      DO i = 1, vecLength
         hessMat(i,i) = 1.0 / mixParam
         invHessMat(i,i) = mixParam
      END DO
   END IF

   WRITE(*,*) 'bfgs_b2 - Point B'

   DO i = 1, historyLength
      stepIndex = MOD(iStep - historyLength + i - 1,maxHistoryLength) + 1

      ! First the Hessian matrix
      ! second correction term
      tempMatB(:,:) = 0.0
      tempVecA(:) = 0.0
      CALL dgemv('N',vecLength,vecLength,1.0,hessMat(:,:),vecLength,paramCorrections(:,stepIndex),1,0.0,tempVecA(:),1)
      norm = ddot(vecLength,paramCorrections(:,stepIndex),1,tempVecA(:),1)
      CALL dgemm('N','N',vecLength,vecLength,1,-1.0/norm,tempVecA(:),vecLength,tempVecA(:),1,0.0,tempMatB,vecLength)

      WRITE(2002,*) '======================================='
      WRITE(2002,'(a,2i7)') 'stepIndex, iStep: ', stepIndex, iStep
      WRITE(2002,'(a,8f15.8)') 'paramCorrections(:,stepIndex): ', paramCorrections(:,stepIndex)
      WRITE(2002,*) 'HessMat:'
      DO j = 1, vecLength
         WRITE(2002,'(8f15.10)') hessMat(:,j)
      END DO
      WRITE(2002,'(a,f15.8)') 'norm: ', norm
      WRITE(2002,'(a,8f15.8)') 'tempVecA(:): ', tempVecA(:)
      WRITE(2002,*) 'tempMatB:'
      DO j = 1, vecLength
         WRITE(2002,'(8f15.10)') tempMatB(:,j)
      END DO

      ! first correction term
      norm = ddot(vecLength,gradientCorrections(:,stepIndex),1,paramCorrections(:,stepIndex),1)
      CALL dgemm('N','N',vecLength,vecLength,1,1.0/norm,gradientCorrections(:,stepIndex),vecLength,gradientCorrections(:,stepIndex),1,1.0,tempMatB,vecLength)
      hessMat(:,:) = hessMat(:,:) + tempMatB(:,:)

      WRITE(2002,*) 'HessMat:'
      DO j = 1, vecLength
         WRITE(2002,'(8f15.10)') hessMat(:,j)
      END DO

      ! Now the inverse of the Hessian matrix
      tempMatB(:,:) = 0.0
      tempVecA(:) = 0.0
      CALL dgemv('N',vecLength,vecLength,1.0,invHessMat(:,:),vecLength,gradientCorrections(:,stepIndex),1,0.0,tempVecA(:),1)
      norm = ddot(vecLength,paramCorrections(:,stepIndex),1,gradientCorrections(:,stepIndex),1)

      WRITE(2003,*) '======================================='
      WRITE(2003,'(a,2i7)') 'stepIndex, iStep: ', stepIndex, iStep
      WRITE(2003,'(a,8f15.8)') 'paramCorrections(:,stepIndex): ', paramCorrections(:,stepIndex)
      WRITE(2003,'(a,8f15.8)') 'gradientCorrections(:,stepIndex): ', gradientCorrections(:,stepIndex)
      WRITE(2003,*) 'invHessMat:'
      DO j = 1, vecLength
         WRITE(2003,'(8f15.10)') invHessMat(:,j)
      END DO
      WRITE(2003,'(a,f15.8)') 'norm: ', norm
      WRITE(2003,'(a,8f15.8)') 'tempVecA(:): ', tempVecA(:)

      ! second correction term
      CALL dgemm('N','N',vecLength,vecLength,1,1.0,tempVecA(:),vecLength,paramCorrections(:,stepIndex),1,0.0,tempMatB,vecLength)
      WRITE(2003,*) 'tempMatB:'

      DO j = 1, vecLength
         WRITE(2003,'(8f15.10)') tempMatB(:,j)
      END DO

      tempMatB(:,:) = tempMatB(:,:) + TRANSPOSE(tempMatB(:,:))
      tempMatB(:,:) = -tempMatB(:,:) / norm

      WRITE(2003,*) 'tempMatB:'
      DO j = 1, vecLength
         WRITE(2003,'(8f15.10)') tempMatB(:,j)
      END DO

      ! first correction term
      temp = norm + ddot(vecLength,gradientCorrections(:,stepIndex),1,tempVecA(:),1)
      CALL dgemm('N','N',vecLength,vecLength,1,temp/(norm*norm),paramCorrections(:,stepIndex),vecLength,paramCorrections(:,stepIndex),1,1.0,tempMatB,vecLength)
      invHessMat(:,:) = invHessMat(:,:) + tempMatB(:,:)

      WRITE(2003,*) 'invHessMat:'
      DO j = 1, vecLength
         WRITE(2003,'(8f15.10)') invHessMat(:,j)
      END DO

      !Check whether invHessMat is the inverse of hessMat (for debugging only)
      tempMatB(:,:) = 0.0
      CALL dgemm('N','N',vecLength,vecLength,vecLength,1.0,hessMat(:,:),vecLength,invHessMat(:,:),vecLength,0.0,tempMatB,vecLength)      
      WRITE(2003,*) 'identityMatrix (hessMat*invHessMat):'
      DO j = 1, vecLength
         WRITE(2003,'(8f15.10)') tempMatB(:,j)
      END DO

   END DO

      WRITE(2001,*) '========================================================='
      WRITE(2001,*) 'HessMat:'
      DO i = 1, vecLength
         WRITE(2001,'(8f15.10)') hessMat(:,i)
      END DO
      WRITE(2001,*) 'invHessMat:'
      DO i = 1, vecLength
         WRITE(2001,'(8f15.10)') invHessMat(:,i)
      END DO

   WRITE(*,*) 'bfgs_b2 - Point C'

   ! 2. Determine next point based on the gradient, the inverse of the Hessian, and the box constraints

   freeParams(:) = .TRUE.
   prelimNextPoint(:) = parameters(:)
   prelimGradient(:) = gradient(:)

   converged = .FALSE.
   DO WHILE (.NOT.converged)

   WRITE(*,*) 'bfgs_b2 - Point D'
   WRITE(1999,'(a,8f15.10)') 'prelimNextPoint(:vecLength): ', prelimNextPoint(:vecLength)
   WRITE(1999,'(a,8f15.10)') 'prelimGradient(:vecLength):  ', prelimGradient(:vecLength)

      ! Determine free parameters (parameters not on boundaries with a gradient pointing outwards)
      numFixedParams = 0
      DO i = 1, vecLength
         IF(enabledConstraints(i)) THEN
            IF((prelimNextPoint(i).LE.minConstraints(i)+eps).AND.(prelimGradient(i).GT.0.0)) THEN
               freeParams(i) = .FALSE.
               numFixedParams = numFixedParams + 1
            END IF
            IF((prelimNextPoint(i).GE.maxConstraints(i)-eps).AND.(prelimGradient(i).LT.0.0)) THEN
               freeParams(i) = .FALSE.
               numFixedParams = numFixedParams + 1
            END IF
         END IF
      END DO
      numFreeParams = vecLength - numFixedParams

   WRITE(*,*) 'freeParams(:vecLength): ', freeParams(:vecLength)

      ! Perform basis transformation to obtain inverse Hessian matrix for the free parameter subspace

      ALLOCATE (subspaceBasis(vecLength,numFreeParams))
      subspaceBasis = 0.0
      subspaceGrad = 0.0
      j = 1
      DO i = 1, vecLength
         IF(freeParams(i)) THEN
            subspaceBasis(i,j) = 1.0
            subspaceGrad(j) = prelimGradient(i)
            j = j + 1
         END IF
      END DO

      WRITE(*,*) 'subspaceBasis:'
      DO i = 1, numFreeParams
         WRITE(*,'(8f10.5)') subspaceBasis(:,i)
      END DO

   WRITE(*,*) 'bfgs_b2 - Point E'

      ALLOCATE (tempMatC(vecLength,numFreeParams))
      ALLOCATE (subspaceInvHessMat(numFreeParams,numFreeParams))
      tempMatC(:,:) = 0.0
      subspaceInvHessMat(:,:) = 0.0
      CALL dgemm('N','N',vecLength,numFreeParams,vecLength,1.0,invHessMat(:,:),vecLength,subspaceBasis(:,:),vecLength,0.0,tempMatC(:,:),vecLength)

      WRITE(*,*) 'tempMatC:'
      DO i = 1, numFreeParams
         WRITE(*,'(8f10.5)') tempMatC(:,i)
      END DO

      CALL dgemm('N','N',numFreeParams,numFreeParams,vecLength,1.0,TRANSPOSE(subspaceBasis(:,:)),numFreeParams,tempMatC(:,:),vecLength,0.0,subspaceInvHessMat(:,:),numFreeParams)
      
      DEALLOCATE(tempMatC)
      DEALLOCATE(subspaceBasis)

   WRITE(*,*) 'bfgs_b2 - Point F'

      ! Calculate next point (ignoring constraints)

      WRITE(*,*) 'vecLength: ', vecLength
      WRITE(*,*) 'numFreeParams: ', numFreeParams
      WRITE(*,*) 'SHAPE(subspaceInvHessMat): ', SHAPE(subspaceInvHessMat)
      WRITE(*,*) 'SHAPE(subspaceGrad): ', SHAPE(subspaceGrad)
      WRITE(*,*) 'SHAPE(tempVecA): ', SHAPE(tempVecA)
      WRITE(*,*) 'subspaceGrad(:numFreeParams): ', subspaceGrad(:numFreeParams)
      WRITE(*,*) 'subspaceInvHessMat:'
      DO i = 1, numFreeParams
         WRITE(*,'(8f10.5)') subspaceInvHessMat(:,i)
      END DO

      tempVecA(:) = 0.0
      CALL dgemv('N',numFreeParams,numFreeParams,1.0,subspaceInvHessMat(:,:),numFreeParams,-subspaceGrad(:numFreeParams),1,0.0,tempVecA(:numFreeParams),1)

      WRITE(*,*) 'tempVecA(:numFreeParams): ', tempVecA(:numFreeParams)

   WRITE(*,*) 'bfgs_b2 - Point F - 0'

      j = 1
      DO i = 1, vecLength
         IF(freeParams(i)) THEN
            prelimNextPoint(i) = prelimNextPoint(i) + tempVecA(j)
            j = j + 1
         END IF
      END DO

   WRITE(*,*) 'bfgs_b2 - Point F - 1'

      DEALLOCATE(subspaceInvHessMat)

   WRITE(*,*) 'bfgs_b2 - Point F - 2'

      ! Constrain next point to box and check convergence
      converged = .TRUE.
      DO i = 1, vecLength
         IF (enabledConstraints(i).AND.freeParams(i)) THEN
            IF(prelimNextPoint(i).LT.minConstraints(i)-eps) THEN
               converged = .FALSE.
            END IF
            IF(prelimNextPoint(i).GT.maxConstraints(i)+eps) THEN
               converged = .FALSE.
            END IF
            IF(prelimNextPoint(i).LT.minConstraints(i)) THEN
               prelimNextPoint(i) = minConstraints(i)
            END IF
            IF(prelimNextPoint(i).GT.maxConstraints(i)) THEN
               prelimNextPoint(i) = maxConstraints(i)
            END IF
         END IF
      END DO

   WRITE(*,*) 'bfgs_b2 - Point F - 3'

      ! If not converged calculate prelimGradient for next point according to quadratic model
      IF(.NOT.converged) THEN
         diffPoint(:) = prelimNextPoint(:) - parameters(:)
         diffGrad(:) = 0.0
         CALL dgemv('N',vecLength,vecLength,1.0,hessMat(:,:),vecLength,diffPoint(:),1,0.0,diffGrad(:),1)
         prelimGradient(:) = prelimGradient(:) + diffGrad(:)
      END IF

   WRITE(*,*) 'prelimNextPoint(:vecLength): ', prelimNextPoint(:vecLength)
   WRITE(*,*) 'prelimGradient(:vecLength): ', prelimGradient(:vecLength)

   WRITE(*,*) 'bfgs_b2 - Point G'

   END DO

   ! 3. Project prelimNextPoint onto the allowed subspace given by the equality constraints
   ! Note: I only consider a single equality constraint so this part becomes simpler

   ! Determine free parameters (parameters not on boundaries with a gradient pointing outwards)

   WRITE(*,*) 'bfgs_b2 - Point H'

   converged = .FALSE.
   DO WHILE (.NOT.converged)

      freeParams(:) = .TRUE.
      numFixedParams = 0
      DO i = 1, vecLength
         IF(enabledConstraints(i)) THEN
            IF(prelimNextPoint(i).LT.minConstraints(i) - eps) THEN
               freeParams(i) = .FALSE.
               numFixedParams = numFixedParams + 1
            END IF
            IF(prelimNextPoint(i).GT.maxConstraints(i) + eps) THEN
               freeParams(i) = .FALSE.
               numFixedParams = numFixedParams + 1
            END IF
            IF(prelimNextPoint(i).LT.minConstraints(i)) THEN
               prelimNextPoint(i) = minConstraints(i)
            END IF
            IF(prelimNextPoint(i).GT.maxConstraints(i)) THEN
               prelimNextPoint(i) = maxConstraints(i)
            END IF
         END IF
      END DO
      numFreeParams = vecLength - numFixedParams

      subspaceNextPoint(:) = 0.0
      subspaceLinCombi(:) = 0.0
      subspaceEqualityCrit = equalityCriterion
      j = 1
      DO i = 1, vecLength
         IF(freeParams(i)) THEN
            subspaceNextPoint(j) = prelimNextPoint(i)
            subspaceLinCombi(j) = equalityLinCombi(i)
            j = j + 1
         ELSE
            subspaceEqualityCrit = subspaceEqualityCrit - prelimNextPoint(i) * equalityLinCombi(i)
         END IF
      END DO
      temp = ddot(numFreeParams,subspaceLinCombi(:),1,subspaceNextPoint(:),1)
      converged = .TRUE.
      IF(ABS(temp - subspaceEqualityCrit).GT.largeEps) THEN
         converged = .FALSE.
         subspaceNextPoint(:) = subspaceNextPoint(:) * subspaceEqualityCrit / temp
      END IF
      j = 1
      DO i = 1, vecLength
         IF(freeParams(i)) THEN
            prelimNextPoint(i) = subspaceNextPoint(j)
            j = j + 1
         END IF
      END DO
   END DO

   WRITE(*,*) 'bfgs_b2 - Point I'

   ! 4. Update lastParameters, lastGradient, parameters

   lastParameters(:) = parameters(:)
   lastGradient(:) = gradient(:)
   parameters(:) = prelimNextPoint(:)

   DEALLOCATE(tempMatB)
   DEALLOCATE(invHessMat)
   DEALLOCATE(hessMat)

   WRITE(*,*) 'parameters(:vecLength): ', parameters(:vecLength)

END SUBROUTINE bfgs_b2

END MODULE m_bfgs_b2
