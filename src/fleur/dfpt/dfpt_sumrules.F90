!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_sumrules
   use m_juDFT
   use m_types
   use m_constants

   implicit none
   private

   public :: dfpt_born_huang

contains

   subroutine dfpt_born_huang(atoms, cell, ft_lim, nNZ, Rvecs, indStored, weightNZ, fcm)
      !! Projects the real-space force constants onto the translational, rotational
      !! and Born-Huang invariance conditions. fcm holds the de-normalized force
      !! constants on the coarse supercell grid and is overwritten in place.
      type(t_atoms), intent(in)    :: atoms
      type(t_cell),  intent(in)    :: cell
      integer,       intent(in)    :: ft_lim(2,3)
      integer,       intent(in)    :: nNZ
      integer,       intent(in)    :: Rvecs(:,:), indStored(:,:)
      real,          intent(in)    :: weightNZ(:)
      complex,       intent(inout) :: fcm(:,:,0:,0:,0:)

      integer, parameter :: maxiter = 200
      real,    parameter :: tol     = 1e-8
      ! Voigt order of the second moments: xx, yy, zz, yz, xz, xy
      integer, parameter :: voigt(2,6) = reshape([1,1, 2,2, 3,3, 2,3, 1,3, 1,2],[2,6])
      ! Cartesian pairs of the rotational conditions: (y,z), (z,x), (x,y)
      integer, parameter :: rpair(2,3) = reshape([2,3, 3,1, 1,2],[2,3])

      integer :: dynDim, nGrid, nR(3), nBH
      integer :: iter, iBond, iRot, iBH, jBH, iGrid, iAtom, jAtom, iRow, iVoigt, jVoigt, iC, jC, ix, iy, iz
      integer :: alpha, beta, gamma
      real    :: rCart(3), bondVec(3), maxBondPairLength, fcScale, violationAfter(3), violationBefore(3)
      logical :: l_conv

      real,    allocatable :: atomPos(:,:)
      real,    allocatable :: mom0(:,:,:), mom1(:,:,:,:), mom2(:,:,:,:)
      real,    allocatable :: consTR(:,:,:), consBH(:,:,:), consBH2(:,:,:), gramMat(:,:)
      complex, allocatable :: phi(:,:,:), phiOld(:,:,:), consPhiOverlap(:), lambda(:)
      integer, allocatable :: negRIdx(:), bhPairs(:,:)

      dynDim = size(fcm,1)
      nR     = [size(fcm,3), size(fcm,4), size(fcm,5)]
      nGrid  = nR(1)*nR(2)*nR(3)
      nBH    = 15 ! number of Born-Huang conditions

      allocate(mom0(0:nGrid-1,atoms%nat,atoms%nat), mom1(3,0:nGrid-1,atoms%nat,atoms%nat), mom2(6,0:nGrid-1,atoms%nat,atoms%nat))
      mom0 = 0.0
      mom1 = 0.0
      mom2 = 0.0
 
      ! go cartesian, constrained moments are defined in that subspace 
      ! we reproduce atomPos so it matches with the Rvecs here that we compute
      ! this is only necessary as we use 2 inputs fullsym and normal one
      ! remove once we do dfpt with only one set of input files. 
      allocate(atomPos(3,atoms%nat))
      do iAtom = 1, atoms%nat
         atomPos(:,iAtom) = matmul(cell%amat, atoms%taual(:,iAtom))
      end do

      ! d and -d give the same conditions, but R and the offset have to carry the
      ! same sign: mixing them pairs a force constant with the wrong bond.
      do iBond = 1, nNZ
         iGrid = indStored(1,iBond) + nR(1)*(indStored(2,iBond) + nR(2)*indStored(3,iBond))
         rCart = matmul(cell%amat, real(Rvecs(:,iBond)))
         do iAtom = 1, atoms%nat
            do jAtom = 1, atoms%nat
               bondVec = rCart + atomPos(:,jAtom) - atomPos(:,iAtom)
               mom0(iGrid,jAtom,iAtom)   = mom0(iGrid,jAtom,iAtom)   + weightNZ(iBond)
               mom1(:,iGrid,jAtom,iAtom) = mom1(:,iGrid,jAtom,iAtom) + weightNZ(iBond)*bondVec
               do iVoigt = 1, 6
                  mom2(iVoigt,iGrid,jAtom,iAtom) = mom2(iVoigt,iGrid,jAtom,iAtom) &
                                                  + weightNZ(iBond)*bondVec(voigt(1,iVoigt))*bondVec(voigt(2,iVoigt))
               end do
            end do
         end do
      end do

      ! get maximum Bond pair length 
      ! used for comparing the magniutde of each change introduced by the moments
      maxBondPairLength = max(sqrt(maxval(mom2(1,:,:,:)+mom2(2,:,:,:)+mom2(3,:,:,:))), 1e-12) 

      ! -R for hermicity 
      allocate(negRIdx(0:nGrid-1))
      do iz = 0, nR(3)-1
         do iy = 0, nR(2)-1
            do ix = 0, nR(1)-1
               iGrid = ix + nR(1)*(iy + nR(2)*iz)
               negRIdx(iGrid) = modulo(-ix-2*ft_lim(1,1),nR(1)) + nR(1)*(modulo(-iy-2*ft_lim(1,2),nR(2)) &
                              + nR(2)*modulo(-iz-2*ft_lim(1,3),nR(3)))
            end do
         end do
      end do

      allocate(bhPairs(2,nBH))
      iBH = 0
      do iVoigt = 1, 5
         do jVoigt = iVoigt+1, 6
            iBH = iBH + 1
            bhPairs(:,iBH) = [iVoigt,jVoigt]
         end do
      end do

      allocate(phi(dynDim,dynDim,0:nGrid-1), phiOld(dynDim,dynDim,0:nGrid-1))
      phi = reshape(fcm, [dynDim,dynDim,nGrid])
      phiOld = phi
      fcScale = max(maxval(abs(phi)), 1e-30)

      call check_conditions(dynDim, atoms%nat, nGrid, nBH, bhPairs, voigt, rpair, mom0, mom1, mom2, phi, violationBefore)
      violationBefore = violationBefore/[fcScale, fcScale*maxBondPairLength, fcScale*maxBondPairLength**2]

      allocate(consTR(dynDim,0:nGrid-1,6), gramMat(max(6,nBH),max(6,nBH)), consPhiOverlap(max(6,nBH)), lambda(max(6,nBH)))
      allocate(consBH(dynDim,dynDim,0:nGrid-1), consBH2(dynDim,dynDim,0:nGrid-1))

      ! The hermicity at the end of this loop can not be enforced by one shot, it is 
      ! not within lagrange multipliers. 
      ! We opt to do 0th and 1st moment and afterwards the 2nd moment iteratively 
      ! the matrix scales linear, instead of N**2 if we do all at once. 
      ! Also a reason for iter, as the space of the moments dont commute. 
      l_conv = .false.
      do iter = 1, maxiter
         ! Translational and rotational conditions; both act inside a single row of phi
         ! and share the same constraint set for all three Cartesian components.
         do iAtom = 1, atoms%nat
            consTR = 0.0
            do beta = 1, 3
               do jAtom = 1, atoms%nat
                  consTR(3*(jAtom-1)+beta,:,beta) = mom0(:,jAtom,iAtom)
               end do
            end do
            do iRot = 1, 3
               beta  = rpair(1,iRot)
               gamma = rpair(2,iRot)
               do jAtom = 1, atoms%nat
                  consTR(3*(jAtom-1)+beta ,:,3+iRot) = consTR(3*(jAtom-1)+beta ,:,3+iRot) + mom1(gamma,:,jAtom,iAtom)
                  consTR(3*(jAtom-1)+gamma,:,3+iRot) = consTR(3*(jAtom-1)+gamma,:,3+iRot) - mom1(beta ,:,jAtom,iAtom)
               end do
            end do

            do iC = 1, 6
               do jC = 1, 6
                  gramMat(iC,jC) = sum(consTR(:,:,iC)*consTR(:,:,jC))
               end do
            end do

            do alpha = 1, 3
               iRow = 3*(iAtom-1) + alpha
               do iC = 1, 6
                  consPhiOverlap(iC) = sum(consTR(:,:,iC)*phi(iRow,:,:))
               end do
               call lsq_solve(6, gramMat, consPhiOverlap, lambda)
               do iC = 1, 6
                  phi(iRow,:,:) = phi(iRow,:,:) - lambda(iC)*consTR(:,:,iC)
               end do
            end do
         end do

         ! Born-Huang conditions; these couple all rows compute big matrix
         do iBH = 1, nBH
            call build_bh_constraint(dynDim, atoms%nat, nGrid, bhPairs(:,iBH), voigt, mom2, consBH)
            consPhiOverlap(iBH) = sum(consBH*phi)
            gramMat(iBH,iBH) = sum(consBH*consBH)
            do jBH = iBH+1, nBH
               call build_bh_constraint(dynDim, atoms%nat, nGrid, bhPairs(:,jBH), voigt, mom2, consBH2)
               gramMat(iBH,jBH) = sum(consBH*consBH2)
               gramMat(jBH,iBH) = gramMat(iBH,jBH)
            end do
         end do
         call lsq_solve(nBH, gramMat, consPhiOverlap, lambda)
         do iBH = 1, nBH
            call build_bh_constraint(dynDim, atoms%nat, nGrid, bhPairs(:,iBH), voigt, mom2, consBH)
            phi = phi - lambda(iBH)*consBH
         end do

         ! Hermicity
         phiOld = phi
         do iGrid = 0, nGrid-1
            phi(:,:,iGrid) = 0.5*(phiOld(:,:,iGrid) + conjg(transpose(phiOld(:,:,negRIdx(iGrid)))))
         end do

         call check_conditions(dynDim, atoms%nat, nGrid, nBH, bhPairs, voigt, rpair, mom0, mom1, mom2, phi, violationAfter)
         violationAfter = violationAfter/[fcScale, fcScale*maxBondPairLength, fcScale*maxBondPairLength**2]
         if (all(violationAfter < tol)) then
            l_conv = .true.
            exit
         end if
      end do

      phiOld = reshape(fcm, [dynDim,dynDim,nGrid])
      write (oUnit,'(a)')       ' Born-Huang projection of the force constants'
      write (oUnit,'(a,i0,a)')  '   iterations                   : ', min(iter,maxiter)
      write (oUnit,'(a)')       '                                   before        after'
      write (oUnit,'(a,2es13.3)') '   translational violation      : ', violationBefore(1), violationAfter(1)
      write (oUnit,'(a,2es13.3)') '   rotational violation         : ', violationBefore(2), violationAfter(2)
      write (oUnit,'(a,2es13.3)') '   Born-Huang violation         : ', violationBefore(3), violationAfter(3)
      write (oUnit,'(a,es13.3)')  '   relative change of the FCM   : ', sqrt(sum(abs(phi-phiOld)**2)/sum(abs(phiOld)**2))

      fcm = reshape(phi, [dynDim,dynDim,nR(1),nR(2),nR(3)])

      if (.not.l_conv) call juDFT_warn("Born-Huang projection did not fully converge.", calledby="dfpt_sumrules.F90")

   end subroutine dfpt_born_huang

   subroutine build_bh_constraint(dynDim, nAtoms, nGrid, pair, voigt, mom2, consBH)
      !! Constraint array of one Born-Huang condition [ab,cd] = [cd,ab].
      integer, intent(in)  :: dynDim, nAtoms, nGrid, pair(2), voigt(2,6)
      real,    intent(in)  :: mom2(:,0:,:,:)
      real,    intent(out) :: consBH(:,:,0:)

      integer :: iAtom, jAtom, alpha, beta, gamma, delta

      alpha = voigt(1,pair(1))
      beta  = voigt(2,pair(1))
      gamma = voigt(1,pair(2))
      delta = voigt(2,pair(2))

      consBH = 0.0
      do iAtom = 1, nAtoms
         do jAtom = 1, nAtoms
            consBH(3*(iAtom-1)+alpha,3*(jAtom-1)+beta ,:) = consBH(3*(iAtom-1)+alpha,3*(jAtom-1)+beta ,:) + mom2(pair(2),:,jAtom,iAtom)
            consBH(3*(iAtom-1)+gamma,3*(jAtom-1)+delta,:) = consBH(3*(iAtom-1)+gamma,3*(jAtom-1)+delta,:) - mom2(pair(1),:,jAtom,iAtom)
         end do
      end do
   end subroutine build_bh_constraint

   subroutine check_conditions(dynDim, nAtoms, nGrid, nBH, bhPairs, voigt, rpair, mom0, mom1, mom2, phi, violation)
      !! Compute the conditions and evaluate how strong each is fulfilled 
      integer, intent(in)  :: dynDim, nAtoms, nGrid, nBH, bhPairs(:,:), voigt(2,6), rpair(2,3)
      real,    intent(in)  :: mom0(0:,:,:), mom1(:,0:,:,:), mom2(:,0:,:,:)
      complex, intent(in)  :: phi(:,:,0:)
      real,    intent(out) :: violation(3)

      integer :: iAtom, jAtom, alpha, beta, gamma, delta, iRow, iBH, iRot
      complex :: temp

      violation = 0.0

      do iAtom = 1, nAtoms
         do alpha = 1, 3
            iRow = 3*(iAtom-1) + alpha
            do beta = 1, 3
               temp = cmplx(0.0,0.0)
               do jAtom = 1, nAtoms
                  temp = temp + sum(phi(iRow,3*(jAtom-1)+beta,:)*mom0(:,jAtom,iAtom))
               end do
               violation(1) = max(violation(1), abs(temp))
            end do
            do iRot = 1, 3
               beta  = rpair(1,iRot)
               gamma = rpair(2,iRot)
               temp = cmplx(0.0,0.0)
               do jAtom = 1, nAtoms
                  temp = temp + sum(phi(iRow,3*(jAtom-1)+beta ,:)*mom1(gamma,:,jAtom,iAtom)) &
                            - sum(phi(iRow,3*(jAtom-1)+gamma,:)*mom1(beta ,:,jAtom,iAtom))
               end do
               violation(2) = max(violation(2), abs(temp))
            end do
         end do
      end do

      do iBH = 1, nBH
         alpha = voigt(1,bhPairs(1,iBH))
         beta  = voigt(2,bhPairs(1,iBH))
         gamma = voigt(1,bhPairs(2,iBH))
         delta = voigt(2,bhPairs(2,iBH))
         temp = cmplx(0.0,0.0)
         do iAtom = 1, nAtoms
            do jAtom = 1, nAtoms
               temp = temp + sum(phi(3*(iAtom-1)+alpha,3*(jAtom-1)+beta ,:)*mom2(bhPairs(2,iBH),:,jAtom,iAtom)) &
                         - sum(phi(3*(iAtom-1)+gamma,3*(jAtom-1)+delta,:)*mom2(bhPairs(1,iBH),:,jAtom,iAtom))
            end do
         end do
         violation(3) = max(violation(3), abs(temp))
      end do
   end subroutine check_conditions

   subroutine lsq_solve(nCon, gramMat, rhs, lambda)
      !! Least-norm solution of gramMat*lambda = rhs via a pseudo-inverse; the constraint
      !! set is linearly dependent whenever the geometry degenerates, e.g. for films.
      integer, intent(in)  :: nCon
      real,    intent(in)  :: gramMat(:,:)
      complex, intent(in)  :: rhs(:)
      complex, intent(out) :: lambda(:)

      integer :: iC, info, lwork
      real    :: emax
      real,    allocatable :: evec(:,:), eigval(:), work(:)
      complex, allocatable :: coef(:)

      allocate(evec(nCon,nCon), eigval(nCon), coef(nCon))
      evec = gramMat(1:nCon,1:nCon)
      lwork = 8*nCon
      allocate(work(lwork))
      call dsyev('V','U',nCon,evec,nCon,eigval,work,lwork,info)
      if (info /= 0) call juDFT_error("Constraint Gram matrix diagonalization failed.", calledby="dfpt_sumrules.F90")

      do iC = 1, nCon
         coef(iC) = sum(evec(:,iC)*rhs(1:nCon))
      end do
      emax = maxval(abs(eigval))
      do iC = 1, nCon
         if (abs(eigval(iC)) > 1e-10*emax) then
            coef(iC) = coef(iC)/eigval(iC)
         else
            coef(iC) = cmplx(0.0,0.0)
         end if
      end do
      do iC = 1, nCon
         lambda(iC) = sum(evec(iC,:)*coef)
      end do
   end subroutine lsq_solve

end module m_dfpt_sumrules
