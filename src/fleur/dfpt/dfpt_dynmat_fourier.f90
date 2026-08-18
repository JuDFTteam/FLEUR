MODULE m_dfpt_dynmat_fourier
   USE m_juDFT
   USE m_types
   USE m_constants
   !USE m_npy

   IMPLICIT NONE
CONTAINS
   SUBROUTINE ft_dyn(atoms, qpts, sym, ft_lim, amat, dyn_mat_q, dyn_mat_r, dyn_mat_q_full)
      !! Transforms the dynamical matrices for a set of q vectors in the
      !! irreducible Brillouin zone onto the full set of q vector in the BZ
      !! and subsequently transforms it to real space (lattice vector grid), to
      !! calculate the mass-normalized Force Constant Matrix.
      type(t_atoms), intent(in)       :: atoms
      type(t_kpts),  intent(in)       :: qpts
      type(t_sym),   intent(in)       :: sym
      integer,       intent(in)       :: ft_lim(2,3)
      real,          intent(in)       :: amat(3,3)
      complex,       intent(in)       :: dyn_mat_q(:,:,:) ! (dyn_dim,dyn_dim,nqpt)
      complex,intent(out)             :: dyn_mat_r(:,:,0:,0:,0:) ! (dyn_dim,dyn_dim,n1,n2,n3)
      complex, allocatable, intent(out) :: dyn_mat_q_full(:,:,:)

      integer :: mrot(3,3),invmrot(3,3)
      logical :: l_inv

      integer :: iq, dyn_dim, iqfull
      integer :: isym
      integer :: iz, iy, ix, nx, ny, nz 
      real    :: q_full(3), trans(3)
      complex, allocatable :: fft_grid(:,:,:) ! (dyn_dim,dyn_dim,nqptf)
      complex, allocatable :: dyn_mat_qsym(:,:)

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_qsym(dyn_dim,dyn_dim))
      allocate(dyn_mat_q_full(dyn_dim,dyn_dim,qpts%nkptf))
      allocate(fft_grid(dyn_dim,dyn_dim,qpts%nkptf))
      fft_grid(:,:,:) = cmplx(0.0,0.0)

      do iqfull = 1, qpts%nkptf
         ! Get q vector index and that of its representative in the irreducible wedge
         iq = qpts%bkp(iqfull)
         ! Fold vector back to 1st BZ if necessary
         q_full = qpts%bkf(:,iqfull)
         isym = qpts%bksym(iqfull)
         call sym%get_sym_operation_int_coord(isym,mrot,invmrot,trans,l_inv)
         if (.not. all(trans == 0 )) call juDFT_error("dynMat interpolation with non symmorphic symmetries is currently not supported. & 
                                                      Please redo the calculation of the IBZ and interpolation with a symmorphic group.",calledby="dfpt_dynmat_sym.f90")
         if (l_inv) isym = isym - sym%nop ! the corresponding symmetry operation 
         dyn_mat_qsym(:,:) = cmplx(0.0,0.0)
         call rotate_dynmat(atoms,sym,isym,mrot,invmrot,l_inv,amat,qpts%bk(:,iq),dyn_mat_q(:,:,iq),dyn_mat_qsym)
         dyn_mat_qsym(:,:) = dyn_mat_qsym(:,:)
         dyn_mat_q_full(:,:,iqfull) = dyn_mat_qsym

         ! Perform the actual FT onto the lattice vector grid
         call ft_dyn_direct(ft_lim,1,q_full,dyn_mat_qsym,fft_grid)
      end do 

      fft_grid(:,:,:)=fft_grid(:,:,:)/qpts%nkptf 
      ! unroll the FCM on supercell index
      call unfold_grid(ft_lim,fft_grid,dyn_mat_r)

   END SUBROUTINE

   SUBROUTINE ft_dyn_direct(ft_lim,isn,bqpt,dyn_mat_q,dyn_mat_r)
      INTEGER, INTENT(IN) :: ft_lim(2,3), isn
      REAL,    INTENT(IN) :: bqpt(3)

      COMPLEX,INTENT(INOUT) :: dyn_mat_q(:,:)
      COMPLEX,INTENT(INOUT) :: dyn_mat_r(:,:,:)

      INTEGER :: iGrid, ix, iy, iz
      REAL    :: phas
      COMPLEX :: phase_fac
      iGrid=0
      DO iz=ft_lim(1,3),ft_lim(2,3)
         DO iy=ft_lim(1,2),ft_lim(2,2)
            DO ix=ft_lim(1,1),ft_lim(2,1)
               iGrid = iGrid+1
               phas=isn*tpi_const*(bqpt(1)*ix+bqpt(2)*iy+bqpt(3)*iz)
               phase_fac=cmplx(cos(phas),sin(phas))
               IF (isn==1) THEN
                  dyn_mat_r(:,:,iGrid) = dyn_mat_r(:,:,iGrid) + phase_fac*dyn_mat_q(:,:)
               ELSE IF (isn==-1) THEN
                  dyn_mat_q(:,:)    = dyn_mat_q(:,:)    + phase_fac*dyn_mat_r(:,:,iGrid)
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE

   subroutine build_ws_ft(ft_lim,bigBoxLim,weights,nNZ,Rvecs,indStored,weightNZ)
      ! Precompute the compact list of nonzero Wigner-Seitz weights 
      !   Rvecs(:,i) : big-box integer coordinates (ix,iy,iz) for the phase factor
      !   indStored(:,i) : the modulo-folded 0-based storage indices (nx,ny,nz)
      !   weightNZ(i)    : the (nonzero) weight value
      integer, intent(in)  :: ft_lim(2,3), bigBoxLim(2,3)
      real,    intent(in)  :: weights(:)
      integer, intent(out) :: nNZ
      integer, allocatable, intent(out) :: Rvecs(:,:), indStored(:,:)
      real,    allocatable, intent(out) :: weightNZ(:)

      integer :: iGrid, ix, iy, iz, nx, ny, nz, i

      nNZ = count(weights /= 0.0)
      allocate(Rvecs(3,nNZ), indStored(3,nNZ), weightNZ(nNZ))

      iGrid = 0
      i = 0
      do iz=bigBoxLim(1,3),bigBoxLim(2,3)
         do iy=bigBoxLim(1,2),bigBoxLim(2,2)
            do ix=bigBoxLim(1,1),bigBoxLim(2,1)
               iGrid = iGrid+1
               if (weights(iGrid) == 0.0) cycle
               i = i+1
               Rvecs(:,i) = (/ix,iy,iz/)
               ! map to smaller box using ft_lim bounds, shift to 0-based storage
               nx = ft_lim(1,1) + modulo(ix - ft_lim(1,1), ft_lim(2,1) - ft_lim(1,1) + 1)
               ny = ft_lim(1,2) + modulo(iy - ft_lim(1,2), ft_lim(2,2) - ft_lim(1,2) + 1)
               nz = ft_lim(1,3) + modulo(iz - ft_lim(1,3), ft_lim(2,3) - ft_lim(1,3) + 1)
               indStored(:,i) = (/nx - ft_lim(1,1), ny - ft_lim(1,2), nz - ft_lim(1,3)/)
               weightNZ(i) = weights(iGrid)
            end do
         end do
      end do
   end subroutine

   subroutine ft_fcm_weight_packed(isn,nNZ,Rvecs,indStored,weightNZ,bqpt,dyn_mat_q,dyn_mat_r)
      ! Backward FT from a bigger box, iterating the precomputed nonzero-weight
      ! points (see build_ws_ft). isn = -1 : r -> k (all callers), isn = 1 : k -> r.
      integer, intent(in) :: isn
      integer, intent(in) :: nNZ
      integer, intent(in) :: Rvecs(:,:), indStored(:,:)
      real,    intent(in) :: weightNZ(:)
      real,    intent(in) :: bqpt(3)

      complex, intent(inout) :: dyn_mat_q(:,:)
      complex, intent(in)    :: dyn_mat_r(:,:,0:,0:,0:)

      integer :: i
      real    :: phas
      complex :: phase_fac

      do i = 1, nNZ
         phas = isn*tpi_const*(bqpt(1)*Rvecs(1,i)+bqpt(2)*Rvecs(2,i)+bqpt(3)*Rvecs(3,i))
         phase_fac = cmplx(cos(phas),sin(phas))
         dyn_mat_q(:,:) = dyn_mat_q(:,:) + phase_fac*weightNZ(i)*dyn_mat_r(:,:,indStored(1,i),indStored(2,i),indStored(3,i))
      end do
   end subroutine

   subroutine ft_fcm_weight(isn,ft_lim,bigBoxLim,weights,bqpt,dyn_mat_q,dyn_mat_r)
      ! Perform a weighted fourier transform 
      integer, intent(in) :: isn
      integer, intent(in) :: ft_lim(2,3), bigBoxlim(2,3)
      real,    intent(in) :: bqpt(3)
      real,    intent(in) :: weights(:)

      complex,intent(inout) :: dyn_mat_q(:,:)
      complex,intent(in) :: dyn_mat_r(:,:,0:,0:,0:)

      integer :: nNZ
      integer, allocatable :: Rvecs(:,:), indStored(:,:)
      real,    allocatable :: weightNZ(:)

      call build_ws_ft(ft_lim,bigBoxLim,weights,nNZ,Rvecs,indStored,weightNZ)
      call ft_fcm_weight_packed(isn,nNZ,Rvecs,indStored,weightNZ,bqpt,dyn_mat_q,dyn_mat_r)
   end subroutine

   subroutine unfold_grid(ft_lim, grid, cube)
        !  Unfold a flat Fourier grid (as accumulated by the discrete fourier transform into
        !  a cube that starts at index 0, using the same iz-outer / iy / ix-inner ordering as ft_dyn_direct.
        integer, intent(in)  :: ft_lim(2,3)
        complex, intent(in)  :: grid(:,:,:)
        complex, intent(out) :: cube(:,:,0:,0:,0:)

        integer :: ix, iy, iz, nx, ny, nz, iGrid

        iGrid = 1
        do iz=ft_lim(1,3),ft_lim(2,3)
            do iy=ft_lim(1,2),ft_lim(2,2)
                do ix=ft_lim(1,1),ft_lim(2,1)
                    nx = ix - ft_lim(1,1)
                    ny = iy - ft_lim(1,2)
                    nz = iz - ft_lim(1,3)
                    cube(:,:,nx,ny,nz) = grid(:,:,iGrid)
                    iGrid = iGrid + 1
                end do !ix
            end do !iy
        end do !iz
    end subroutine unfold_grid

   SUBROUTINE rotate_dynmat(atoms,sym,isym,mrot,invmrot,l_inv,amat,bqpt,dyn,dyn_mat_qsym)
      !! Applies a symmetry operation to the dynamical matrix of an IBZ q vector
      !! to find the matrix of its mapped q vector in the full BZ. This is done
      !! by using the symmetry relation of the FCM when a symmetry operation
      !! \(\underline{B}\) maps atoms \(\beta',\alpha'\) onto \(\beta,\alpha\)
      !! in Cartesian coordinates 
      !! (\(\underline{B}\boldsymbol{\tau}_{\beta'}=\boldsymbol{\tau}_{\beta}\))
      !! $$\underline{\Phi}_{\alpha'+\boldsymbol{R},\beta'}=\underline{B}\underline{\Phi}_{\alpha+\underline{B}\boldsymbol{R},\beta}\underline{B}^{-1}$$
      !! Resulting (with the definition of the DM as the mass-scaled Fourier Transform
      !! of the FCM) in the corresponding relation:
      !! $$\underline{D}_{\alpha',\beta'}(\boldsymbol{q})=p(\alpha,\beta)*\underline{B}\underline{D}_{\alpha,\beta}(\boldsymbol{q}_{\mathrm{rep}})\underline{B}^{-1}$$
      !! with a phase factor
      !! $$f(\alpha,\beta)=exp(ix),x=\boldsymbol{q}_{\mathrm{red}}\cdot(\boldsymbol{\tau}_{\beta'}-\boldsymbol{\tau}_{\alpha'})-\boldsymbol{q}_{\mathrm{red}}\cdot(\boldsymbol{\tau}_{\beta}-\boldsymbol{\tau}_{\alpha})$$
      !! which can be written as a product
      !! $$f(\alpha,\beta) = f(\alpha)f^{*}(\beta), f(\alpha)=exp(i(\boldsymbol{q}_{\mathrm{red}}\cdot\boldsymbol{\tau}_{\alpha}-\boldsymbol{q}_{\mathrm{red}}\cdot\boldsymbol{\tau}_{\alpha'})).$$
      !! The real space rotation is related to the rotation matrix of the symmetry operation by the Bravais matrix of the system
      !! $$\underline{B}=\underline{A}\underline{S}^{-1}\underline{A}^{-1}$$
      USE m_inv3
      type(t_atoms), INTENT(IN)    :: atoms
      type(t_sym),   INTENT(IN)    :: sym
      INTEGER,       INTENT(IN)    :: isym
      INTEGER,          INTENT(IN)    :: mrot(3,3)
      INTEGER,          INTENT(IN)    :: invmrot(3,3) 
      LOGICAL,       INTENT(IN)    :: l_inv
      REAL,          INTENT(IN)    :: amat(3,3)
      REAL,          INTENT(IN)    :: bqpt(3)
      COMPLEX,       INTENT(IN)    :: dyn(:,:)
      COMPLEX,       INTENT(INOUT) :: dyn_mat_qsym(:,:)
      
      INTEGER :: iAtom, jAtom , iAtom_map, jAtom_map
      INTEGER :: iAlpha, iBeta, iAlpha_map, iBeta_map
      REAL    :: invamat(3,3), phas, det , pos_map(3)
      COMPLEX :: brot(3,3), temp_mat_1(3,3), temp_mat_2(3,3)
      COMPLEX :: phase_fac

      INTEGER :: map(atoms%nat)
      COMPLEX :: phase_map(atoms%nat)

      CALL inv3(amat,invamat,det)
      temp_mat_1 = MATMUL(invmrot,invamat)
      brot = MATMUL(amat,temp_mat_1)
      
      ! Find the q vector in the full BZ
      !q_full = MATMUL(transpose(mrot),bqpt)
      
      ! Calculate the array of phases
      DO iAtom = 1, atoms%nat
         jAtom = sym%mapped_atom(isym,iAtom)
         ! Create list to which atom we map 
         map(iAtom)=jAtom
         ! Find rotated atom position
         pos_map = matmul(mrot, atoms%taual(:,iAtom))
         IF (l_inv) pos_map = -1 * pos_map ! Inversion does not exist in real space, 
                                           ! we have to introduce a minus that negates the minus of mrot
         ! Calculate the phase factor f(alpha)
         phas= tpi_const*(dot_product(bqpt(:),atoms%taual(:,jAtom)) - dot_product(bqpt(:), pos_map(:))) 
         phase_fac=cmplx(cos(phas),sin(phas))
         phase_map(iAtom)=phase_fac
      END DO
 
      IF (l_inv) phase_map = conjg(phase_map) ! inversion maps q -> -q which results in a complex conjugate of the phase
      ! Transform the dynamical matrix from the representative atom and q vector to the unfolded ones
      DO iAtom=1, atoms%nat
         iAlpha = 3*(iAtom-1)
         iAlpha_map = 3*(map(iAtom)-1)
         iAtom_map = map(iAtom)
         DO jAtom=1, atoms%nat
            iBeta = 3*(jAtom-1)
            iBeta_map = 3*(map(jAtom)-1)
            jAtom_map = map(jAtom)
            temp_mat_1 = dyn(iAlpha+1:iAlpha+3,iBeta+1:iBeta+3)
            if (l_inv) temp_mat_1 = conjg(temp_mat_1) ! inversion maps q -> -q  , which results into a conjugate
            temp_mat_2 = MATMUL(brot,temp_mat_1)
            temp_mat_1 = MATMUL(temp_mat_2,TRANSPOSE(brot))
            phase_fac=phase_map(iAtom_map)*conjg(phase_map(jAtom_map))

            dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) &
          = dyn_mat_qsym(iAlpha_map+1:iAlpha_map+3,iBeta_map+1:iBeta_map+3) + phase_fac*temp_mat_1
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE ift_dyn(atoms,qpts,nNZ,Rvecs,indStored,weightNZ,bqpt,dyn_mat_r,dyn_mat_q)
      !! Transforms the dynamical matrix on a real space lattice vector grid
      !! (--> mass-normalized FCM) back onto a specific q vector provided as
      !! input (bqpt) by the inverse Fourier Transformation as compared to
      type(t_atoms), intent(in) :: atoms
      type(t_kpts),  intent(in) :: qpts
      integer,       intent(in) :: nNZ
      integer,       intent(in) :: Rvecs(:,:), indStored(:,:)
      real,          intent(in) :: weightNZ(:)
      real,          intent(in) :: bqpt(3)
      complex,       intent(in) :: dyn_mat_r(:,:,0:,0:,0:)
      complex, allocatable, intent(out) :: dyn_mat_q(:,:)

      integer :: dyn_dim

      dyn_dim = 3*atoms%nat

      allocate(dyn_mat_q(dyn_dim,dyn_dim))
      dyn_mat_q(:,:) = cmplx(0.0,0.0)

      call ft_fcm_weight_packed(-1,nNZ,Rvecs,indStored,weightNZ,bqpt,dyn_mat_q,dyn_mat_r)

   END SUBROUTINE

   SUBROUTINE interpolate_dynmat(atoms, sym, cell, qpts_coarse, dyn_mat_coarse,l_WSinterpol, q_target, dyn_mat_interp)
      !! Fourier interpolation of the dynamical matrix from a coarse q-mesh onto an
      !! arbitrary set of target q-points. The coarse (IBZ, "raw") dynamical matrices
      !! are unfolded and transformed to the real-space (mass-normalized) FCM once via
      !! ft_dyn, then de-normalized and inverse-transformed (ift_dyn) to every target q.
      !! The result is NOT diagonalized: it is returned in the convention that
      !! DiagonalizeDynMat(..., l_scalemass=.TRUE.) expects, so the caller does the
      !! diagonalization.
      type(t_atoms), intent(in) :: atoms
      type(t_sym),   intent(in) :: sym
      type(t_cell),  intent(in) :: cell
      type(t_kpts),  intent(in) :: qpts_coarse                  ! coarse q-mesh
      complex,       intent(in) :: dyn_mat_coarse(:,:,:)        ! (3nat,3nat,nq_coarse)
      logical,       intent(in) :: l_WSinterpol
      real,          intent(in) :: q_target(:,:)               ! (3,nTarget)
      complex, allocatable, intent(out) :: dyn_mat_interp(:,:,:) ! (3nat,3nat,nTarget), NOT diagonalized

      type(t_cell) :: cellLocal
      integer :: dyn_dim, iq, iDir, iDir2, nx, ny, nz, iGrid, boxSize
      integer :: ft_lim(2,3), bigBox_lim(2,3)
      real    :: mass_mat(3*atoms%nat, 3*atoms%nat)
      integer, allocatable :: supercellR(:,:)
      real,    allocatable :: FTweight(:)
      integer :: nNZ                              ! packed nonzero-weight list (built once)
      integer, allocatable :: Rvecs(:,:), indStored(:,:)
      real,    allocatable :: weightNZ(:)
      complex, allocatable :: dyn_mat_r(:,:,:,:,:), dyn_mat_q_full(:,:,:), dyn_mat_q(:,:)

      dyn_dim = 3*atoms%nat
      
      cellLocal = cell 

      ! Fourier box limits
      ft_lim(2,:) =  qpts_coarse%nkpt3(:)/2
      ft_lim(1,:) = ft_lim(2,:) - qpts_coarse%nkpt3(:) + 1

      if (l_WSinterpol) then
         ! create Wigner-Seitz cell and weights on a bigger fft mesh
         bigBox_lim(2,:) =   2*qpts_coarse%nkpt3(:)
         bigBox_lim(1,:) = - 2*qpts_coarse%nkpt3(:)
         boxSize = (4*qpts_coarse%nkpt3(1)+1) * (4*qpts_coarse%nkpt3(2)+1) * (4*qpts_coarse%nkpt3(3)+1)
         allocate(FTweight(boxSize))
         allocate(supercellR(3,boxSize))
         FTweight = 0.0
         supercellR = 0
         iGrid = 1
         do nz=bigBox_lim(1,3),bigBox_lim(2,3)
            do ny=bigBox_lim(1,2),bigBox_lim(2,2)
               do nx=bigBox_lim(1,1),bigBox_lim(2,1)
                  supercellR(:,iGrid) = (/nx,ny,nz/)
                  iGrid = iGrid+1
               end do
            end do
         end do
         call cellLocal%calculate_WSweight(supercellR,FTweight,scaleSupercell=qpts_coarse%nkpt3(:))
      else
         ! simple back-transformation
         bigBox_lim = ft_lim
         boxSize = (qpts_coarse%nkpt3(1)) * (qpts_coarse%nkpt3(2)) * (qpts_coarse%nkpt3(3))
         allocate(FTweight(boxSize))
         FTweight = 1.0
      end if

      call build_ws_ft(ft_lim, bigBox_lim, FTweight, nNZ, Rvecs, indStored, weightNZ)

      ! coarse q dynamical matrices --> real-space (mass-normalized) FCM
      allocate(dyn_mat_r(dyn_dim,dyn_dim,0:(qpts_coarse%nkpt3(1)-1),0:(qpts_coarse%nkpt3(2)-1),0:(qpts_coarse%nkpt3(3)-1)))
      call ft_dyn(atoms, qpts_coarse, sym, ft_lim, cell%amat, dyn_mat_coarse, dyn_mat_r, dyn_mat_q_full)

      ! De-normalize the FCM: the normal diagonalization routines re-apply the mass
      ! scaling (l_scalemass=.TRUE.), so the FCM must be non-normalized here.
      do iDir = 1, dyn_dim
         do iDir2 = 1, dyn_dim
            mass_mat(iDir,iDir2) = massInElectronMasses * &
               SQRT(atomicMasses_const(atoms%nz(atoms%itype(CEILING(iDir/3.0))))*atomicMasses_const(atoms%nz(atoms%itype(CEILING(iDir2/3.0)))))
         end do
      end do
      do nz = 0, qpts_coarse%nkpt3(3)-1
         do ny = 0, qpts_coarse%nkpt3(2)-1
            do nx = 0, qpts_coarse%nkpt3(1)-1
               dyn_mat_r(:,:,nx,ny,nz) = dyn_mat_r(:,:,nx,ny,nz) * mass_mat(:,:)
            end do
         end do
      end do

      ! TODO: NAC / polar long-range term (currently inactive in dfpt_interpolation)

      ! inverse Fourier transform onto every target q-point
      allocate(dyn_mat_interp(dyn_dim,dyn_dim,size(q_target,2)))
      do iq = 1, size(q_target,2)
         call ift_dyn(atoms, qpts_coarse, nNZ, Rvecs, indStored, weightNZ, q_target(:,iq), dyn_mat_r, dyn_mat_q)
         dyn_mat_interp(:,:,iq) = dyn_mat_q
         deallocate(dyn_mat_q)
      end do

   END SUBROUTINE interpolate_dynmat

END MODULE m_dfpt_dynmat_fourier
