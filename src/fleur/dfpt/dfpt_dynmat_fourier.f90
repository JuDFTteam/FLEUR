!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_dynmat_fourier
   USE m_juDFT
   USE m_constants
   USE m_inv3
   USE m_dfpt_sumrules
   USE m_types_atoms
   USE m_types_cell
   USE m_types_kpts
   USE m_types_sym
   !USE m_npy

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: ft_dyn, ft_ws_forward, build_ws_ft, ft_fcm_weight_packed, rotate_dynmat, ift_dyn, &
      interpolate_dynmat, remap_fullsym_atoms, build_ws_shifted, build_ws_pairs, ft_fcm_weight_pairs
CONTAINS
   SUBROUTINE ft_dyn(atoms, qpts, sym, nNZ, Rvecs, amat, dyn_mat_q, dyn_mat_r, dyn_mat_q_full)
      !! Transforms the dynamical matrices for a set of q vectors in the
      !! irreducible Brillouin zone onto the full set of q vector in the BZ
      !! and subsequently onto the Wigner-Seitz lattice vectors Rvecs, to
      !! calculate the mass-normalized Force Constant Matrix.
      type(t_atoms), intent(in)       :: atoms
      type(t_kpts),  intent(in)       :: qpts
      type(t_sym),   intent(in)       :: sym
      integer,       intent(in)       :: nNZ
      integer,       intent(in)       :: Rvecs(:,:)
      real,          intent(in)       :: amat(3,3)
      complex,       intent(in)       :: dyn_mat_q(:,:,:) ! (dyn_dim,dyn_dim,nqpt)
      complex,intent(out)             :: dyn_mat_r(:,:,0:) ! (dyn_dim,dyn_dim,nNZ)
      complex, allocatable, intent(out) :: dyn_mat_q_full(:,:,:)

      integer :: mrot(3,3),invmrot(3,3)
      logical :: l_inv

      integer :: iq, dyn_dim, iqfull
      integer :: isym
      real    :: q_full(3), trans(3)
      complex, allocatable :: dyn_mat_qsym(:,:)

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_qsym(dyn_dim,dyn_dim))
      allocate(dyn_mat_q_full(dyn_dim,dyn_dim,qpts%nkptf))
      dyn_mat_r(:,:,:) = cmplx(0.0,0.0)

      do iqfull = 1, qpts%nkptf
         ! Get q vector index and that of its representative in the irreducible wedge
         iq = qpts%bkp(iqfull)
         ! Fold vector back to 1st BZ if necessary
         q_full = qpts%bkf(:,iqfull)
         isym = qpts%bksym(iqfull)
         call sym%get_sym_operation_int_coord(isym,mrot,invmrot,trans,l_inv)
         if (l_inv) isym = isym - sym%nop ! the corresponding symmetry operation 
         dyn_mat_qsym(:,:) = cmplx(0.0,0.0)
         call rotate_dynmat(atoms,sym,isym,invmrot,l_inv,amat,trans,q_full,dyn_mat_q(:,:,iq),dyn_mat_qsym)
         dyn_mat_qsym(:,:) = dyn_mat_qsym(:,:)
         dyn_mat_q_full(:,:,iqfull) = dyn_mat_qsym

         call ft_ws_forward(nNZ,Rvecs,q_full,dyn_mat_qsym,dyn_mat_r)
      end do

      dyn_mat_r(:,:,:)=dyn_mat_r(:,:,:)/qpts%nkptf

   END SUBROUTINE

   subroutine ft_ws_forward(nNZ,Rvecs,bqpt,mat_q,mat_r)
      !! Adds e^{+i 2pi q.R} mat_q to mat_r at every lattice vector R = Rvecs(:,i),
      !! stored at index i-1.
      integer, intent(in)    :: nNZ
      integer, intent(in)    :: Rvecs(:,:)
      real,    intent(in)    :: bqpt(3)
      complex, intent(in)    :: mat_q(:,:)
      complex, intent(inout) :: mat_r(:,:,0:)

      integer :: i
      real    :: phas

      do i = 1, nNZ
         phas = tpi_const*(bqpt(1)*Rvecs(1,i)+bqpt(2)*Rvecs(2,i)+bqpt(3)*Rvecs(3,i))
         mat_r(:,:,i-1) = mat_r(:,:,i-1) + cmplx(cos(phas),sin(phas))*mat_q(:,:)
      end do
   end subroutine

   subroutine build_ws_ft(bigBoxLim,weights,nNZ,Rvecs,weightNZ)
      !! Compact list of the lattice vectors Rvecs(:,i) with nonzero Wigner-Seitz
      !! weight weightNZ(i); the real-space data of Rvecs(:,i) is stored at index i-1.
      integer, intent(in)  :: bigBoxLim(2,3)
      real,    intent(in)  :: weights(:)
      integer, intent(out) :: nNZ
      integer, allocatable, intent(out) :: Rvecs(:,:)
      real,    allocatable, intent(out) :: weightNZ(:)

      integer :: iGrid, ix, iy, iz, i

      nNZ = count(weights /= 0.0)
      allocate(Rvecs(3,nNZ), weightNZ(nNZ))

      iGrid = 0
      i = 0
      do iz=bigBoxLim(1,3),bigBoxLim(2,3)
         do iy=bigBoxLim(1,2),bigBoxLim(2,2)
            do ix=bigBoxLim(1,1),bigBoxLim(2,1)
               iGrid = iGrid+1
               if (weights(iGrid) == 0.0) cycle
               i = i+1
               Rvecs(:,i) = (/ix,iy,iz/)
               weightNZ(i) = weights(iGrid)
            end do
         end do
      end do
   end subroutine

   subroutine ft_fcm_weight_packed(isn,nNZ,Rvecs,weightNZ,bqpt,dyn_mat_q,dyn_mat_r)
      !! Weighted backward FT r -> k (isn = -1) from the lattice vectors of build_ws_ft.
      integer, intent(in) :: isn
      integer, intent(in) :: nNZ
      integer, intent(in) :: Rvecs(:,:)
      real,    intent(in) :: weightNZ(:)
      real,    intent(in) :: bqpt(3)

      complex, intent(inout) :: dyn_mat_q(:,:)
      complex, intent(in)    :: dyn_mat_r(:,:,0:)

      integer :: i
      real    :: phas
      complex :: phase_fac

      do i = 1, nNZ
         phas = isn*tpi_const*(bqpt(1)*Rvecs(1,i)+bqpt(2)*Rvecs(2,i)+bqpt(3)*Rvecs(3,i))
         phase_fac = cmplx(cos(phas),sin(phas))
         dyn_mat_q(:,:) = dyn_mat_q(:,:) + phase_fac*weightNZ(i)*dyn_mat_r(:,:,i-1)
      end do
   end subroutine

   subroutine build_ws_shifted(amat, scaleSupercell, l_film, shifts, nNZ, Rvecs, weightNZ)
      !! Compact list of the lattice vectors Rvecs(:,i) for which R + shifts(:,iShift) (cartesian) lies in the
      !! Wigner-Seitz cell of the supercell, with weight weightNZ(iShift,i); the data of Rvecs(:,i) is stored at index i-1.
      !! For films only in-plane lattice vectors and supercell images are used.
      real,    intent(in)  :: amat(3,3)
      integer, intent(in)  :: scaleSupercell(3)
      logical, intent(in)  :: l_film
      real,    intent(in)  :: shifts(:,:)
      integer, intent(out) :: nNZ
      integer, allocatable, intent(out) :: Rvecs(:,:)
      real,    allocatable, intent(out) :: weightNZ(:,:)

      integer, parameter :: imageLim = 3 ! supercell images searched in each direction
      integer :: bigBoxLim(2,3), boxSize, nShift, nImage, imageLimZ, iShift, iGrid, iImage, iMem, i, ix, iy, iz
      real    :: tol, vecRCart(3), vecShift(3), distanceR, distance
      logical :: l_inside
      integer, allocatable :: supercellR(:,:)
      real,    allocatable :: FTweight(:,:), imageCart(:,:)

      nShift = size(shifts,2)
      tol    = 1e-6*dot_product(amat(:,1),amat(:,1))

      imageLimZ = imageLim
      if (l_film) imageLimZ = 0
      nImage = (2*imageLim+1)**2*(2*imageLimZ+1)
      allocate(imageCart(3,nImage))
      iImage = 0
      do iz = -imageLimZ, imageLimZ
         do iy = -imageLim, imageLim
            do ix = -imageLim, imageLim
               iImage = iImage + 1
               imageCart(:,iImage) = matmul(amat, real([ix,iy,iz]*scaleSupercell))
            end do
         end do
      end do

      bigBoxLim(2,:) =  2*scaleSupercell(:)
      bigBoxLim(1,:) = -2*scaleSupercell(:)
      if (l_film) bigBoxLim(:,3) = 0
      boxSize = product(bigBoxLim(2,:) - bigBoxLim(1,:) + 1)
      allocate(supercellR(3,boxSize), FTweight(nShift,boxSize))
      FTweight = 0.0

      iGrid = 0
      do iz = bigBoxLim(1,3), bigBoxLim(2,3)
         do iy = bigBoxLim(1,2), bigBoxLim(2,2)
            do ix = bigBoxLim(1,1), bigBoxLim(2,1)
               iGrid = iGrid + 1
               supercellR(:,iGrid) = [ix,iy,iz]
               vecRCart = matmul(amat, real([ix,iy,iz]))
               do iShift = 1, nShift
                  vecShift = vecRCart + shifts(:,iShift)
                  distanceR = dot_product(vecShift,vecShift)
                  l_inside = .true.
                  iMem = 0
                  do iImage = 1, nImage
                     distance = dot_product(vecShift+imageCart(:,iImage), vecShift+imageCart(:,iImage))
                     if (distance < distanceR - tol) then
                        l_inside = .false.
                        exit
                     end if
                     if (distance < distanceR + tol) iMem = iMem + 1
                  end do
                  if (l_inside) FTweight(iShift,iGrid) = 1.0/iMem
               end do
            end do
         end do
      end do

      nNZ = count(any(FTweight /= 0.0, dim=1))
      allocate(Rvecs(3,nNZ), weightNZ(nShift,nNZ))
      i = 0
      do iGrid = 1, boxSize
         if (all(FTweight(:,iGrid) == 0.0)) cycle
         i = i + 1
         Rvecs(:,i) = supercellR(:,iGrid)
         weightNZ(:,i) = FTweight(:,iGrid)
      end do

      do iShift = 1, nShift
         if (abs(sum(weightNZ(iShift,:)) - product(scaleSupercell)) > 1e-6) then
            call juDFT_error("Wigner-Seitz weights do not sum to the supercell size.", calledby="build_ws_shifted")
         end if
      end do
   end subroutine

   subroutine build_ws_pairs(atoms, amat, scaleSupercell, l_film, nNZ, Rvecs, weightPair)
      !! Wigner-Seitz list of the force constants with weightPair(i,j,iR) for the bond R + tau_j - tau_i
      !! between atom i and atom j, on the union of all atom pairs.
      type(t_atoms), intent(in)  :: atoms
      real,          intent(in)  :: amat(3,3)
      integer,       intent(in)  :: scaleSupercell(3)
      logical,       intent(in)  :: l_film
      integer,       intent(out) :: nNZ
      integer, allocatable, intent(out) :: Rvecs(:,:)
      real,    allocatable, intent(out) :: weightPair(:,:,:)

      logical, parameter :: l_bondWS = .true. ! .false.: plain lattice Wigner-Seitz weights for every pair
      integer :: iAtom, jAtom
      real, allocatable :: shifts(:,:), weightNZ(:,:)

      allocate(shifts(3,atoms%nat**2))
      shifts = 0.0
      if (l_bondWS) then
         do jAtom = 1, atoms%nat
            do iAtom = 1, atoms%nat
               shifts(:,iAtom+atoms%nat*(jAtom-1)) = matmul(amat, atoms%taual(:,jAtom) - atoms%taual(:,iAtom))
            end do
         end do
      end if
      call build_ws_shifted(amat, scaleSupercell, l_film, shifts, nNZ, Rvecs, weightNZ)
      weightPair = reshape(weightNZ, [atoms%nat, atoms%nat, nNZ])
   end subroutine

   subroutine ft_fcm_weight_pairs(isn,nNZ,Rvecs,weightPair,bqpt,dyn_mat_q,dyn_mat_r)
      !! Backward FT r -> q (isn = -1) with a separate Wigner-Seitz weight for every 3x3 atom block.
      integer, intent(in) :: isn
      integer, intent(in) :: nNZ
      integer, intent(in) :: Rvecs(:,:)
      real,    intent(in) :: weightPair(:,:,:)
      real,    intent(in) :: bqpt(3)

      complex, intent(inout) :: dyn_mat_q(:,:)
      complex, intent(in)    :: dyn_mat_r(:,:,0:)

      integer :: i, iAtom, jAtom, iAlpha, iBeta
      real    :: phas
      complex :: phase_fac

      do i = 1, nNZ
         phas = isn*tpi_const*(bqpt(1)*Rvecs(1,i)+bqpt(2)*Rvecs(2,i)+bqpt(3)*Rvecs(3,i))
         phase_fac = cmplx(cos(phas),sin(phas))
         do jAtom = 1, size(weightPair,2)
            iBeta = 3*(jAtom-1)
            do iAtom = 1, size(weightPair,1)
               if (weightPair(iAtom,jAtom,i) == 0.0) cycle
               iAlpha = 3*(iAtom-1)
               dyn_mat_q(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3) = dyn_mat_q(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3) + phase_fac*weightPair(iAtom,jAtom,i)*dyn_mat_r(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3,i-1)
            end do
         end do
      end do
   end subroutine

   SUBROUTINE rotate_dynmat(atoms,sym,isym,invmrot,l_inv,amat,trans,q_full,dyn,dyn_mat_qsym)
      !! Applies a symmetry operation to the dynamical matrix of an IBZ q vector
      !! to find the matrix of its mapped q vector in the full BZ:
      !! $$\underline{D}_{p(\alpha),p(\beta)}(\boldsymbol{q}_{\mathrm{full}})=f(\alpha)f^{*}(\beta)\,\underline{B}\,\underline{D}_{\alpha,\beta}(\boldsymbol{q}_{\mathrm{rep}})\,\underline{B}^{-1}$$
      !! with \(\underline{B}=\underline{A}S\underline{A}^{-1}\) and the lattice-vector phase
      !! $$f(\alpha)=e^{+i2\pi\boldsymbol{q}_{\mathrm{full}}\cdot\boldsymbol{R}_{S\alpha}},\quad\boldsymbol{R}_{S\alpha}=S\boldsymbol{\tau}_{\alpha}+\boldsymbol{v}-\boldsymbol{\tau}_{p(\alpha)}$$
      !! read at the source atom \(\alpha\). For `isym > nop` the spatial \(S\) is
      !! `-invmrot`, the source block is conjugated, and `q_full` already carries
      !! the \(\boldsymbol{q}\to-\boldsymbol{q}\) sign.
      !!
      !! `dfpt_lambda::dfpt_unfold_gmat` carries \(f^{*}\) on the same
      !! \(\boldsymbol{R}_{S\alpha}\). That is not an inconsistency: the two meet
      !! on opposite sides of \(g^{\nu}=\sum_{\kappa\alpha}e^{\nu}_{\kappa\alpha}
      !! g^{\kappa\alpha}\), so the phases have to cancel there. See
      !! `changed_phase_convention.md`.
      type(t_atoms), INTENT(IN)    :: atoms
      type(t_sym),   INTENT(IN)    :: sym
      INTEGER,       INTENT(IN)    :: isym
      INTEGER,       INTENT(IN)    :: invmrot(3,3)
      LOGICAL,       INTENT(IN)    :: l_inv
      REAL,          INTENT(IN)    :: amat(3,3)
      REAL,          INTENT(IN)    :: trans(3)
      REAL,          INTENT(IN)    :: q_full(3)
      COMPLEX,       INTENT(IN)    :: dyn(:,:)
      COMPLEX,       INTENT(INOUT) :: dyn_mat_qsym(:,:)

      INTEGER :: iAtom, jAtom
      INTEGER :: iAlpha, iBeta, iAlpha_map, iBeta_map
      INTEGER :: rotReal(3,3)
      REAL    :: invamat(3,3), phas, det, rlat(3)
      COMPLEX :: brot(3,3), temp_mat_1(3,3), temp_mat_2(3,3)
      COMPLEX :: phase_fac

      INTEGER :: map(atoms%nat)
      COMPLEX :: phase_map(atoms%nat)

      CALL inv3(amat,invamat,det)

      rotReal = invmrot
      IF (l_inv) rotReal = -invmrot ! get_sym_operation_int_coord negates the spatial part for isym > nop

      temp_mat_1 = MATMUL(REAL(rotReal),invamat)
      brot = MATMUL(amat,temp_mat_1)

      DO iAtom = 1, atoms%nat
         jAtom = sym%mapped_atom(isym,iAtom)
         map(iAtom) = jAtom
         rlat = MATMUL(REAL(rotReal),atoms%taual(:,iAtom)) + trans - atoms%taual(:,jAtom)
         phas = tpi_const*dot_product(q_full,rlat)
         phase_map(iAtom) = cmplx(cos(phas),sin(phas))
      END DO

      ! Transform the dynamical matrix from the representative atom and q vector to the unfolded ones
      DO iAtom=1, atoms%nat
         iAlpha = 3*(iAtom-1)
         iAlpha_map = 3*(map(iAtom)-1)
         DO jAtom=1, atoms%nat
            iBeta = 3*(jAtom-1)
            iBeta_map = 3*(map(jAtom)-1)
            temp_mat_1 = dyn(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3)
            if (l_inv) temp_mat_1 = conjg(temp_mat_1) ! inversion maps q -> -q  , which results into a conjugate
            temp_mat_2 = MATMUL(brot,temp_mat_1)
            temp_mat_1 = MATMUL(temp_mat_2,TRANSPOSE(brot))
            phase_fac = phase_map(iAtom)*conjg(phase_map(jAtom))

            dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) &
          = dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) + phase_fac*temp_mat_1
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE ift_dyn(atoms,qpts,nNZ,Rvecs,weightPair,bqpt,dyn_mat_r,dyn_mat_q)
      !! Transforms the dynamical matrix on a real space lattice vector grid
      !! (--> mass-normalized FCM) back onto a specific q vector provided as
      !! input (bqpt) by the inverse Fourier Transformation as compared to
      type(t_atoms), intent(in) :: atoms
      type(t_kpts),  intent(in) :: qpts
      integer,       intent(in) :: nNZ
      integer,       intent(in) :: Rvecs(:,:)
      real,          intent(in) :: weightPair(:,:,:)
      real,          intent(in) :: bqpt(3)
      complex,       intent(in) :: dyn_mat_r(:,:,0:)
      complex, allocatable, intent(out) :: dyn_mat_q(:,:)

      integer :: dyn_dim

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_q(dyn_dim,dyn_dim))
      dyn_mat_q(:,:) = cmplx(0.0,0.0)

      call ft_fcm_weight_pairs(-1,nNZ,Rvecs,weightPair,bqpt,dyn_mat_q,dyn_mat_r)

   END SUBROUTINE

   SUBROUTINE interpolate_dynmat(atoms, sym, cell, qpts_coarse, dyn_mat_coarse,l_WSinterpol, q_target, dyn_mat_interp, l_bornhuang)
      !! Fourier interpolation of the dynamical matrix from a coarse q-mesh onto an
      !! arbitrary set of target q-points. The coarse (IBZ, "raw") dynamical matrices
      !! are unfolded and transformed to the real-space (mass-normalized) FCM once via
      !! ft_dyn, then de-normalized and inverse-transformed (ift_dyn) to every target q.
      !! The result is NOT diagonalized: it is returned unscaled, so the caller applies
      !! hermitize_scale_dynmat before the diagonalization.

      type(t_atoms), intent(in) :: atoms
      type(t_sym),   intent(in) :: sym
      type(t_cell),  intent(in) :: cell
      type(t_kpts),  intent(in) :: qpts_coarse                  ! coarse q-mesh
      complex,       intent(in) :: dyn_mat_coarse(:,:,:)        ! (3nat,3nat,nq_coarse)
      logical,       intent(in) :: l_WSinterpol
      real,          intent(in) :: q_target(:,:)               ! (3,nTarget)
      complex, allocatable, intent(out) :: dyn_mat_interp(:,:,:) ! (3nat,3nat,nTarget), NOT diagonalized
      logical, optional, intent(in) :: l_bornhuang

      integer :: dyn_dim, iq, iDir, iDir2, iR
      integer :: ft_lim(2,3)
      real    :: mass_mat(3*atoms%nat, 3*atoms%nat)
      real,    allocatable :: FTweight(:)
      integer :: nNZ                              ! packed nonzero-weight list (built once)
      integer, allocatable :: Rvecs(:,:)
      real,    allocatable :: weightNZ(:), weightPair(:,:,:)
      complex, allocatable :: dyn_mat_r(:,:,:), dyn_mat_q_full(:,:,:), dyn_mat_q(:,:)
      logical :: l_bh

      dyn_dim = 3*atoms%nat

      l_bh = .false.
      if (present(l_bornhuang)) l_bh = l_bornhuang

      if (l_WSinterpol) then
         call build_ws_pairs(atoms, cell%amat, qpts_coarse%nkpt3, cell%z1 > 0.0, nNZ, Rvecs, weightPair)
      else
         ! simple back-transformation on the plain Fourier box
         ft_lim(2,:) =  qpts_coarse%nkpt3(:)/2
         ft_lim(1,:) = ft_lim(2,:) - qpts_coarse%nkpt3(:) + 1
         allocate(FTweight(product(qpts_coarse%nkpt3)))
         FTweight = 1.0
         call build_ws_ft(ft_lim, FTweight, nNZ, Rvecs, weightNZ)
         allocate(weightPair(atoms%nat,atoms%nat,nNZ))
         do iR = 1, nNZ
            weightPair(:,:,iR) = weightNZ(iR)
         end do
      end if

      ! coarse q dynamical matrices --> real-space (mass-normalized) FCM
      allocate(dyn_mat_r(dyn_dim,dyn_dim,0:nNZ-1))
      call ft_dyn(atoms, qpts_coarse, sym, nNZ, Rvecs, cell%amat, dyn_mat_coarse, dyn_mat_r, dyn_mat_q_full)

      ! De-normalize the FCM: hermitize_scale_dynmat re-applies the mass scaling,
      ! so the FCM must be non-normalized here.
      do iDir = 1, dyn_dim
         do iDir2 = 1, dyn_dim
            mass_mat(iDir,iDir2) = massInElectronMasses * &
               SQRT(atomicMasses_const(atoms%nz(atoms%itype(CEILING(iDir/3.0))))*atomicMasses_const(atoms%nz(atoms%itype(CEILING(iDir2/3.0)))))
         end do
      end do
      do iR = 0, nNZ-1
         dyn_mat_r(:,:,iR) = dyn_mat_r(:,:,iR) * mass_mat(:,:)
      end do

      if (l_bh) then
         call timestart("Born-Huang projection")
         call dfpt_born_huang(atoms, cell, qpts_coarse%nkpt3, nNZ, Rvecs, weightPair, dyn_mat_r)
         call timestop("Born-Huang projection")
      end if

      ! TODO: NAC / polar long-range term (currently inactive in dfpt_interpolation)

      ! inverse Fourier transform onto every target q-point
      allocate(dyn_mat_interp(dyn_dim,dyn_dim,size(q_target,2)))
      do iq = 1, size(q_target,2)
         call ift_dyn(atoms, qpts_coarse, nNZ, Rvecs, weightPair, q_target(:,iq), dyn_mat_r, dyn_mat_q)
         dyn_mat_interp(:,:,iq) = dyn_mat_q
         deallocate(dyn_mat_q)
      end do

   END SUBROUTINE interpolate_dynmat

   subroutine remap_fullsym_atoms(atoms, atoms_full, sym_full)
      !! Rewrites sym_full%mapped_atom, built on the atom order of atoms_full, into the
      !! atom order of atoms. Atoms are matched by their position modulo lattice vectors.
      type(t_atoms), intent(in)    :: atoms
      type(t_atoms), intent(in)    :: atoms_full
      type(t_sym),   intent(inout) :: sym_full

      real, parameter :: tol = 1e-7

      integer :: iAtom, jAtom
      real    :: diff(3)
      integer, allocatable :: perm(:), permInv(:), mappedFull(:,:)

      if (atoms_full%nat /= atoms%nat) then
         call juDFT_error("fullsym_ cell and DFPT cell have a different number of atoms.", calledby="remap_fullsym_atoms")
      end if

      ! perm(iAtom): index in atoms_full of atom iAtom of atoms
      allocate(perm(atoms%nat), permInv(atoms%nat))
      perm = 0
      permInv = 0
      do iAtom = 1, atoms%nat
         do jAtom = 1, atoms%nat
            diff = atoms_full%taual(:,jAtom) - atoms%taual(:,iAtom)
            if (all(abs(diff - nint(diff)) < tol) .and. atoms_full%nz(atoms_full%itype(jAtom)) == atoms%nz(atoms%itype(iAtom))) then
               perm(iAtom) = jAtom
               exit
            end if
         end do
         if (perm(iAtom) == 0 .or. permInv(max(perm(iAtom),1)) /= 0) then
            call juDFT_error("fullsym_ cell and DFPT cell contain different atoms.", calledby="remap_fullsym_atoms")
         end if
         permInv(perm(iAtom)) = iAtom
      end do

      mappedFull = sym_full%mapped_atom
      do iAtom = 1, atoms%nat
         sym_full%mapped_atom(:,iAtom) = permInv(mappedFull(:,perm(iAtom)))
      end do
   end subroutine remap_fullsym_atoms

END MODULE m_dfpt_dynmat_fourier
