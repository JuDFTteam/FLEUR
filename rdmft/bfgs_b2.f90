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

   USE m_constants

   IMPLICIT NONE

   INTEGER,        INTENT(IN)    :: vecLength
   INTEGER,        INTENT(IN)    :: maxHistoryLength
   INTEGER,        INTENT(INOUT) :: iStep ! The index of the bfgs_b iteration
   REAL,           INTENT(IN)    :: mixParam
   REAL,           INTENT(IN)    :: equalityCriterion
   REAL,           INTENT(INOUT)    :: gradient(vecLength)
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
   REAL                   :: temp, norm, subspaceEqualityCrit, maxDisp, dampingFactor
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

   WRITE(*,*) 'bfgs_b2 (input): parameters(:vecLength): ', parameters(:vecLength)
   WRITE(*,*) 'bfgs_b2 (input): gradient(:vecLength): ', gradient(:vecLength)

   WRITE(5000,*) '========================================================================='
   WRITE(5000,*) 'ENTERED BFGS_B2: '

   !!! Test start
!   DO i = 1, vecLength - 1
!      IF(ABS(gradient(i)).LT.eps) CYCLE
!      gradient(i) = sqrt(ABS(gradient(i)))*2.0*ATAN(gradient(i))/pi_const
!   END DO
   !!! Test end

   l_converged = .FALSE.

   ! 0. Add last data pair to history iff it is not obviously destructive for the positive definiteness of the Hessian matrix

   IF(ANY(lastGradient(:).NE.0.0).OR.ANY(lastParameters(:).NE.0.0)) THEN

      WRITE(2501,'(a,12f10.6)'), '    gradient: ', gradient(:)
      WRITE(2501,'(a,12f10.6)'), 'lastGradient: ', lastGradient(:)
      WRITE(2502,'(a,12f10.6)'), '    parameters: ', parameters(:)
      WRITE(2502,'(a,12f10.6)'), 'lastParameters: ', lastParameters(:)

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
      WRITE(*,*) 'MAXVAL(ABS(lastGradCorrection(:))): ', MAXVAL(ABS(lastGradCorrection(:)))
      WRITE(*,*) 'MAXVAL(ABS(lastParamCorrection(:))): ', MAXVAL(ABS(lastParamCorrection(:)))
      WRITE(*,*) 'MAXVAL(ABS(lastGradCorrection(:vecLength-1))): ', MAXVAL(ABS(lastGradCorrection(:vecLength-1)))
      WRITE(*,*) 'MAXVAL(ABS(lastParamCorrection(:vecLength-1))): ', MAXVAL(ABS(lastParamCorrection(:vecLength-1)))
!      IF(ALL(ABS(lastGradCorrection(:)).LT.convCrit).AND.ALL(ABS(lastParamCorrection(:)).LT.convCrit)) THEN
      IF(ALL(ABS(lastParamCorrection(:vecLength-1)).LT.convCrit)) THEN
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

   DO i = 1, historyLength
      stepIndex = MOD(iStep - historyLength + i - 1,maxHistoryLength) + 1

      ! First the Hessian matrix
      ! second correction term
      tempMatB(:,:) = 0.0
      tempVecA(:) = 0.0
      CALL dgemv('N',vecLength,vecLength,1.0,hessMat(:,:),vecLength,paramCorrections(:,stepIndex),1,0.0,tempVecA(:),1)
      norm = ddot(vecLength,paramCorrections(:,stepIndex),1,tempVecA(:),1)
      CALL dgemm('N','N',vecLength,vecLength,1,-1.0/norm,tempVecA(:),vecLength,tempVecA(:),1,0.0,tempMatB,vecLength)

      ! first correction term
      norm = ddot(vecLength,gradientCorrections(:,stepIndex),1,paramCorrections(:,stepIndex),1)
      CALL dgemm('N','N',vecLength,vecLength,1,1.0/norm,gradientCorrections(:,stepIndex),vecLength,gradientCorrections(:,stepIndex),1,1.0,tempMatB,vecLength)
      hessMat(:,:) = hessMat(:,:) + tempMatB(:,:)

      ! Now the inverse of the Hessian matrix
      tempMatB(:,:) = 0.0
      tempVecA(:) = 0.0
      CALL dgemv('N',vecLength,vecLength,1.0,invHessMat(:,:),vecLength,gradientCorrections(:,stepIndex),1,0.0,tempVecA(:),1)
      norm = ddot(vecLength,paramCorrections(:,stepIndex),1,gradientCorrections(:,stepIndex),1)

      ! second correction term
      CALL dgemm('N','N',vecLength,vecLength,1,1.0,tempVecA(:),vecLength,paramCorrections(:,stepIndex),1,0.0,tempMatB,vecLength)

      tempMatB(:,:) = tempMatB(:,:) + TRANSPOSE(tempMatB(:,:))
      tempMatB(:,:) = -tempMatB(:,:) / norm

      ! first correction term
      temp = norm + ddot(vecLength,gradientCorrections(:,stepIndex),1,tempVecA(:),1)
      CALL dgemm('N','N',vecLength,vecLength,1,temp/(norm*norm),paramCorrections(:,stepIndex),vecLength,paramCorrections(:,stepIndex),1,1.0,tempMatB,vecLength)
      invHessMat(:,:) = invHessMat(:,:) + tempMatB(:,:)

      !Check whether invHessMat is the inverse of hessMat (for debugging only)
      tempMatB(:,:) = 0.0
      CALL dgemm('N','N',vecLength,vecLength,vecLength,1.0,hessMat(:,:),vecLength,invHessMat(:,:),vecLength,0.0,tempMatB,vecLength)      

   END DO

   ! 2. Determine next point based on the gradient, the inverse of the Hessian, and the box constraints

   freeParams(:) = .TRUE.
   prelimNextPoint(:) = parameters(:)
   prelimGradient(:) = gradient(:)

   WRITE(5300,*) '=========================================================='
   WRITE(5300,*) 'bfgs_b2: Starting next point search: '
   WRITE(5300,*) '=========================================================='

   converged = .FALSE.
   DO WHILE (.NOT.converged)

      WRITE(5300,*) 'Point A: '
      WRITE(5300,*) 'prelimNextPoint: '
      WRITE(5300,'(5f20.10)') prelimNextPoint(:)
      WRITE(5300,*) 'prelimGradient: '
      WRITE(5300,'(5f20.10)') prelimGradient(:)


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

      ALLOCATE (tempMatC(vecLength,numFreeParams))
      ALLOCATE (subspaceInvHessMat(numFreeParams,numFreeParams))
      tempMatC(:,:) = 0.0
      subspaceInvHessMat(:,:) = 0.0

      WRITE(5000,*) '========================================================================='
      WRITE(5000,*) 'prelimNextPoint: '
      WRITE(5000,'(5f20.10)') prelimNextPoint(:)
      WRITE(5000,*) 'prelimGradient: '
      WRITE(5000,'(5f20.10)') prelimGradient(:)
      WRITE(5000,*) 'subspaceBasis: '
      DO i = 1, SIZE(subspaceBasis,2)
         WRITE(5000,'(5f20.10)') subspaceBasis(:,i)
      END DO

      CALL dgemm('N','N',vecLength,numFreeParams,vecLength,1.0,invHessMat(:,:),vecLength,subspaceBasis(:,:),vecLength,0.0,tempMatC(:,:),vecLength)
      CALL dgemm('N','N',numFreeParams,numFreeParams,vecLength,1.0,TRANSPOSE(subspaceBasis(:,:)),numFreeParams,tempMatC(:,:),vecLength,0.0,subspaceInvHessMat(:,:),numFreeParams)
      
      DEALLOCATE(tempMatC)
      DEALLOCATE(subspaceBasis)

      ! Calculate next point (ignoring constraints)

      tempVecA(:) = 0.0

      WRITE(5000,*) 'numFreeParams: ', numFreeParams
      WRITE(5000,*) 'subspaceInvHessMat: '
      DO i = 1, SIZE(subspaceInvHessMat,2)
         WRITE(5000,'(5f20.10)') subspaceInvHessMat(:,i)
      END DO
      WRITE(5000,*) '-subspaceGrad: '
      WRITE(5000,'(5f20.10)') -subspaceGrad(:numFreeParams)
      WRITE(5000,*) 'tempVecA: '
      WRITE(5000,'(5f20.10)') tempVecA(:numFreeParams)

      CALL dgemv('N',numFreeParams,numFreeParams,1.0,subspaceInvHessMat(:,:),numFreeParams,-subspaceGrad(:numFreeParams),1,0.0,tempVecA(:numFreeParams),1)

      WRITE(5000,*) '----------------------------'
      WRITE(5000,*) 'tempVecA: '
      WRITE(5000,'(5f20.10)') tempVecA(:numFreeParams)
      WRITE(5000,*) 'prelimNextPoint: '
      WRITE(5000,'(5f20.10)') prelimNextPoint(:)

      !!! TEST START
!      temp = SQRT(ddot(numFreeParams,tempVecA(:numFreeParams),1,tempVecA(:numFreeParams),1))
      temp = MAXVAL(ABS(tempVecA(:numFreeParams-1)))
      maxDisp = 0.05
      IF (temp.GT.maxDisp) tempVecA(:numFreeParams) = tempVecA(:numFreeParams) / (temp/maxDisp)
      dampingFactor = 0.8
      !!! TEST END

      j = 1
      DO i = 1, vecLength
         IF(freeParams(i)) THEN
            prelimNextPoint(i) = prelimNextPoint(i) + dampingFactor * tempVecA(j)
!            prelimNextPoint(i) = prelimNextPoint(i) + tempVecA(j)
            j = j + 1
         END IF
      END DO

      WRITE(5000,*) 'prelimNextPoint: '
      WRITE(5000,'(5f20.10)') prelimNextPoint(:)

      DEALLOCATE(subspaceInvHessMat)

      ! Constrain next point to box and check convergence
      converged = .TRUE.
      DO i = 1, vecLength
         IF (enabledConstraints(i).AND.freeParams(i)) THEN
            IF(prelimNextPoint(i).LT.minConstraints(i)-eps) THEN
!               converged = .FALSE.
            END IF
            IF(prelimNextPoint(i).GT.maxConstraints(i)+eps) THEN
!               converged = .FALSE.
            END IF
            IF(prelimNextPoint(i).LT.minConstraints(i)) THEN
               prelimNextPoint(i) = minConstraints(i)
            END IF
            IF(prelimNextPoint(i).GT.maxConstraints(i)) THEN
               prelimNextPoint(i) = maxConstraints(i)
            END IF
         END IF
      END DO

      WRITE(5000,*) 'prelimNextPoint: '
      WRITE(5000,'(5f20.10)') prelimNextPoint(:)
      WRITE(5000,*) 'prelimGradient: '
      WRITE(5000,'(5f20.10)') prelimGradient(:)

      ! If not converged calculate prelimGradient for next point according to quadratic model
      IF(.NOT.converged) THEN
         diffPoint(:) = prelimNextPoint(:) - parameters(:)
         diffGrad(:) = 0.0
         CALL dgemv('N',vecLength,vecLength,1.0,hessMat(:,:),vecLength,diffPoint(:),1,0.0,diffGrad(:),1)
         prelimGradient(:) = prelimGradient(:) + diffGrad(:)
      END IF

      ! Test: Set the gradient entry for the Lagrange multiplier explicitly
      ! (I hope this makes the approach more stable. The here introduced data is not automatically obtained)
      temp = ddot(vecLength,equalityLinCombi(:vecLength),1,prelimNextPoint(:vecLength),1)
      prelimGradient(vecLength) = temp - equalityCriterion   ! derivative with respect to the Lagrange parameter

      WRITE(5000,*) 'prelimGradient: '
      WRITE(5000,'(5f20.10)') prelimGradient(:)

   END DO

   ! 3. Project prelimNextPoint onto the allowed subspace given by the equality constraints
   ! Note: I only consider a single equality constraint so this part becomes simpler

   ! Determine free parameters (parameters not on boundaries with a gradient pointing outwards)

!   WRITE(5200,*) '=========================================================='
!   WRITE(5200,*) 'bfgs_b2: Starting charge constraint ensurance: '
!   WRITE(5200,*) '=========================================================='

   converged = .FALSE.
   DO WHILE (.NOT.converged)

!      WRITE(5200,*) 'Point A: '
!      WRITE(5200,*) 'prelimNextPoint: '
!      WRITE(5200,'(5f20.10)') prelimNextPoint(:)

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

!         WRITE(5100,*) 'subspaceNextPoint(:): '
!         WRITE(5100,*) subspaceNextPoint(:)
!         WRITE(5100,*) 'subspaceEqualityCrit: ', subspaceEqualityCrit
!         WRITE(5100,*) 'temp: ', temp

         subspaceNextPoint(:) = subspaceNextPoint(:) * subspaceEqualityCrit / temp
      END IF
      j = 1
      DO i = 1, vecLength
         IF(freeParams(i)) THEN
            prelimNextPoint(i) = subspaceNextPoint(j)
            j = j + 1
         END IF
      END DO

      ! Also ensure that the equality criterion is fulfilled for the whole space.
      ! (This is only for pathological cases)

      temp = ddot(vecLength,equalityLinCombi(:vecLength),1,prelimNextPoint(:vecLength),1)
      IF(ABS(temp - equalityCriterion).GT.largeEps) THEN
         converged = .FALSE.
         prelimNextPoint(:vecLength) = prelimNextPoint(:vecLength) * equalityCriterion / temp
      END IF

   END DO

   ! 4. Update lastParameters, lastGradient, parameters

   lastParameters(:) = parameters(:)
   lastGradient(:) = gradient(:)
   parameters(:) = prelimNextPoint(:)
!   parameters(:) = 0.25*(prelimNextPoint(:) - lastParameters(:)) + lastParameters(:)

   DEALLOCATE(tempMatB)
   DEALLOCATE(invHessMat)
   DEALLOCATE(hessMat)

   WRITE(*,*) 'bfgs_b2 (output): parameters(:vecLength): ', parameters(:vecLength)

END SUBROUTINE bfgs_b2

END MODULE m_bfgs_b2
