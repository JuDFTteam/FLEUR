!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_lambda
   !! Symmetry bookkeeping for the unfolding of the electron-phonon matrix
   !! elements from the irreducible q wedge onto the full Brillouin zone.
   !!
   !! See README_symmetry.md in this directory for the derivation. Two things are
   !! provided:
   !!
   !! 1. `dfpt_read_fullsym`: the DFPT run itself is restricted to
   !!    \(N_{\mathrm{op}}=1\) (see dfpt_check.F90), so the symmetry operations
   !!    and the symmetry-reduced q mesh have to come from a second,
   !!    `fullsym_`-prefixed input set. This mirrors what dfpt_interpolation.F90
   !!    already does for the dynamical matrix.
   !!
   !! 2. `dfpt_build_lambda`: the gauge matrix
   !!    $$\Lambda_{pn}(\boldsymbol{k},\mathcal{S})
   !!      =\langle\psi_{p,S\boldsymbol{k}}|\hat{\mathcal{S}}|\psi_{n,\boldsymbol{k}}\rangle$$
   !!    relating the symmetry-rotated state to the independently diagonalized
   !!    eigenvector stored at \(S\boldsymbol{k}\). It is unitary and block
   !!    diagonal in the energy, and it depends only on
   !!    \((\boldsymbol{k},\mathcal{S})\) -- not on q and not on the perturbation
   !!    index -- so it is built once for every operation of the group and kept.
   !!
   !! The rotation itself is not reimplemented here. The muffin-tin coefficients
   !! are rotated by `waveftrafo_gen_cmt` (m_trafo), which handles the atom
   !! permutation, the lattice-translation phase, the Wigner d matrices at full
   !! lmax and time reversal, and which covers local orbitals through the stride
   !! `num_radfun_per_l`. The plane-wave block is rotated by
   !! `waveftrafo_gen_zmat`. Its local-orbital rows are left at zero, which is
   !! correct: local orbitals are confined to the muffin tins and contribute
   !! nothing to the interstitial.

#ifdef CPP_MPI
   use mpi
#endif
   use m_juDFT

#ifdef _OPENACC_later
   use cublas
#define CPP_zgemm cublaszgemm
#else
#define CPP_zgemm zgemm
#endif

   use m_types
   use m_constants

   implicit none

   private

   public :: dfpt_read_fullsym, dfpt_build_lambda, dfpt_check_lambda,dfpt_unfold_gmat, dfpt_qmesh_full_bz

contains

   subroutine dfpt_read_fullsym(fmpi, fi, sym_full, qvec_full)
      !! Reads the `fullsym_` input set and hands back the full symmetry group
      !! together with the symmetry-reduced q mesh.
      !!
      !! The convention is the one already established by dfpt_interpolation.F90:
      !! the *k*-point set of the `fullsym_` input is used as the *q* mesh, and
      !! its ordering must agree with `fi%dfpt%qvec`, because the
      !! `q<iQ>_b<b>_j<j>.hdf` and `dynMatq=<iQ>` file names are keyed on that
      !! index. That assumption is checked here rather than left implicit.

      use m_fleur_init

      type(t_mpi),        intent(in)  :: fmpi
      type(t_fleurinput), intent(in)  :: fi
      type(t_sym),        intent(out) :: sym_full
      type(t_kpts),       intent(out) :: qvec_full

      type(t_mpi)        :: fmpi_fullsym
      type(t_fleurinput) :: fi_fullsym
      type(t_sphhar)     :: sphhar_fullsym
      type(t_stars)      :: stars_fullsym
      type(t_nococonv)   :: nococonv_fullsym
      type(t_enpara)     :: enpara_fullsym
      type(t_results)    :: results_fullsym
      type(t_wann)       :: wann_fullsym
      type(t_hybdat)     :: hybdat_fullsym
      type(t_mpdata)     :: mpdata_fullsym

      class(t_xcpot),     allocatable :: xcpot_fullsym
      class(t_forcetheo), allocatable :: forcetheo_fullsym

      character(len=100) :: inp_pref
      integer            :: iQ, ik, isym
      integer            :: mrot(3, 3), invmrot(3, 3)
      logical            :: l_trs
      real               :: trans(3), rkpt(3)

      inp_pref = adjustl("fullsym_")

      fmpi_fullsym%l_mpi_multithreaded = fmpi%l_mpi_multithreaded
      fmpi_fullsym%mpi_comm            = fmpi%mpi_comm

      ! Skip setupMPI: this fleur_init only re-reads the fullsym_ input set; the
      ! parallel-solver setup it would otherwise do is never used here and fails
      ! when the fullsym q-mesh does not factor evenly onto the MPI ranks.
      call fleur_init(fmpi_fullsym, fi_fullsym, sphhar_fullsym, stars_fullsym, nococonv_fullsym, forcetheo_fullsym, &
                      enpara_fullsym, xcpot_fullsym, results_fullsym, wann_fullsym, hybdat_fullsym, mpdata_fullsym, inp_pref, l_skip_setupmpi=.true.)

      sym_full  = fi_fullsym%sym
      qvec_full = fi_fullsym%kpts

      if (sym_full%nop < 1) call juDFT_error("fullsym_ input carries no symmetry operations.", calledby="dfpt_read_fullsym")

      ! The full BZ expansion must be present; without it there is nothing to unfold.
      if (.not.allocated(qvec_full%bkf) .or. .not.allocated(qvec_full%bkp) .or. .not.allocated(qvec_full%bksym)) then
         call juDFT_error("fullsym_ q mesh has no full-BZ expansion (bkf/bkp/bksym missing).", calledby="dfpt_read_fullsym")
      end if

      ! Ordering consistency: the IBZ list must be the one the density-response
      ! files were written for.
      if (qvec_full%nkpt /= fi%dfpt%qvec%nkpt) then
         call juDFT_error("fullsym_ q mesh has a different number of irreducible q than dfpt%qvec.", calledby="dfpt_read_fullsym")
      end if

      do iQ = 1, qvec_full%nkpt
         if (any(abs(qvec_full%bk(:, iQ) - fi%dfpt%qvec%bk(:, iQ)) > 1e-7)) then
            call juDFT_error("fullsym_ q mesh ordering differs from dfpt%qvec; the q<iQ>_* and dynMatq=<iQ> file names would be mismatched.", calledby="dfpt_read_fullsym")
         end if
      end do

      ! The unfolding needs S k on the k mesh for every k and every operation.
      ! This is the analogue of EPW's check_q_points_sym and it is cheaper to
      ! fail here than downstream.
      do isym = 1, sym_full%nsym
         call sym_full%get_sym_operation_int_coord(isym, mrot, invmrot, trans, l_trs)
         do ik = 1, fi%kpts%nkptf
            rkpt = matmul(transpose(mrot), fi%kpts%bkf(:, ik))
            if (fi%kpts%get_nk(rkpt) < 1) then
               call juDFT_error("k mesh is not closed under the full symmetry group; regenerate it with the fullsym_ symmetry.", calledby="dfpt_read_fullsym")
            end if
         end do
      end do

   end subroutine dfpt_read_fullsym

   subroutine dfpt_build_lambda(fi, sym_full, fmpi, enpara, vTot, nococonv, stars, eig_id, jsp, bandWindow, lambda, ikrot)
      !! Builds \(\Lambda(\boldsymbol{k},\mathcal{S})\) for every k of the mesh and
      !! **every operation of the full group at once**, together with the index
      !! map \(\boldsymbol{k}\mapsto S\boldsymbol{k}\).
      !!
      !! \(\Lambda\) is assembled from its muffin-tin and interstitial parts,
      !! mirroring how eigen_hssetup builds the overlap matrix. No rotated
      !! eigenvector in the LAPW basis is required, which is what keeps the local
      !! orbitals out of the plane-wave rotation entirely.
      !!
      !! **Band range.** Only the bands of `bandWindow` are built -- the same ones
      !! the electron-phonon elements carry. \(\Lambda\) is block diagonal in the
      !! energy, so in
      !! $$g_{mn}(S\boldsymbol{k},S\boldsymbol{q})=\sum_{pr}\Lambda_{mp}\,
      !!   g_{pr}(\boldsymbol{k},\boldsymbol{q})\,\Lambda^{*}_{nr}$$
      !! the indices p and r are confined to the degenerate multiplets of m and n.
      !! If those multiplets lie inside the window, bands outside it cannot
      !! contribute and computing them would be pure waste. If a multiplet
      !! straddles an edge of the window, the missing contributions are the
      !! \(g_{pr}\) themselves, which the caller does not have either -- so a
      !! larger \(\Lambda\) would not rescue the unfolding, it would only hide the
      !! problem. `dfpt_check_lambda` is what detects it: the window block of a
      !! unitary matrix is itself unitary exactly when no multiplet straddles an
      !! edge.
      !!
      !! Memory: `lambda` is nb^2 * nkpt * nsym complex with nb the window size,
      !! i.e. 16 * nb^2 * nkpt * nsym bytes. Check this against the machine before
      !! calling with a wide window and a dense mesh.

      use m_eig66_io,   only: read_eig
      use m_trafo,      only: waveftrafo_gen_zmat, waveftrafo_gen_cmt
      use m_genMTBasis, only: genMTBasis
      use m_hs_int_direct

      type(t_fleurinput),   intent(in)  :: fi
      type(t_sym),          intent(in)  :: sym_full
      type(t_mpi),          intent(in)  :: fmpi
      type(t_enpara),       intent(in)  :: enpara
      type(t_potden),       intent(in)  :: vTot
      type(t_nococonv),     intent(in)  :: nococonv
      type(t_stars),        intent(in)  :: stars
      integer,              intent(in)  :: eig_id, jsp
      integer,              intent(in)  :: bandWindow(2)       !! first and last band to carry
      complex, allocatable, intent(out) :: lambda(:, :, :, :)  !! (nb,nb,nkpt,nsym), nb = window size
      integer, allocatable, intent(out) :: ikrot(:, :)         !! (nkpt,nsym): index of S k

      ! --- rotation scaffolding, local to this routine ------------------------
      type(t_hybinp)       :: hybinp   ! map / tvec / d_wgn2 on the full group
      type(t_mpdata)       :: mpdata   ! num_radfun_per_l only
      type(t_usdus)        :: usdus       ! radial overlaps, built below for jsp only
      complex, allocatable :: olapmt(:, :, :, :)
      integer, allocatable :: n_r(:, :)
      integer              :: maxlmindx, maxn_r

      ! scratch for genMTBasis; the radial functions themselves are not needed,
      ! only the overlap integrals they leave behind in usdus
      real, allocatable :: fRad(:, :, :), gRad(:, :, :), floRad(:, :, :)

      type(t_lapw) :: lapw_k, lapw_kPrime
      type(t_mat)  :: zMatK, zMatKPrime, zMatRot
      type(t_mat)  :: smat_is, hmat_dummy
      type(t_mpi)  :: fmpi_serial

      complex, allocatable :: cmt_all(:, :, :, :)   ! (nbands,maxlmindx,nat,nkpt)
      complex, allocatable :: cmt_rot(:, :, :)
      complex, allocatable :: c_phase(:), vpw_zero(:)
      complex, allocatable :: zKPrime_c(:, :), zRot_c(:, :), tmp_is(:, :)

      integer, allocatable :: ikrot_inv(:, :), ev_list(:)
      real,    allocatable :: eig_dummy(:)

      integer :: ik, ikPrime, isym, i, itype, nv_kPrime, neig, nbands
      integer :: mrot(3, 3), invmrot(3, 3)
      logical :: l_trs
      real    :: trans(3), rkpt(3)
#ifdef CPP_MPI
      integer :: ierr
#endif

      nbands = bandWindow(2) - bandWindow(1) + 1
      if (bandWindow(1) < 1 .or. nbands < 1) then
         call juDFT_error("Invalid band window handed to dfpt_build_lambda.", calledby="dfpt_build_lambda")
      end if

      if (fi%sym%nop /= 1) then
         call juDFT_error("dfpt_build_lambda needs the DFPT run to carry no symmetry (nop == 1); local MT frames would otherwise need rot_to_unrotated.", calledby="dfpt_build_lambda")
      end if

      call timestart("dfpt lambda")

      ! Radial basis for this spin. Only the overlap integrals in usdus survive;
      ! they are what the muffin-tin half of Lambda contracts with.
      allocate (fRad(fi%atoms%jmtd, 2, 0:fi%atoms%lmaxd))
      allocate (gRad(fi%atoms%jmtd, 2, 0:fi%atoms%lmaxd))
      allocate (floRad(fi%atoms%jmtd, 2, fi%atoms%nlod))
      call usdus%init(fi%atoms, fi%input%jspins)
      do itype = 1, fi%atoms%ntype
         call genMTBasis(fi%atoms, enpara, vTot, fmpi, itype, jsp, usdus, fRad, gRad, floRad)
      end do
      deallocate (fRad, gRad, floRad)

      call setup_rotation_data(fi, sym_full, usdus, jsp, hybinp, mpdata, olapmt, n_r, maxlmindx, maxn_r)

      allocate (lambda(nbands, nbands, fi%kpts%nkpt, sym_full%nsym), source=cmplx_0)
      allocate (ikrot(fi%kpts%nkpt, sym_full%nsym), source=0)
      allocate (ikrot_inv(fi%kpts%nkpt, sym_full%nsym), source=0)

      allocate (ev_list(nbands))
      ev_list = [(i, i=bandWindow(1), bandWindow(2))]
      allocate (eig_dummy(nbands), source=0.0)
      allocate (c_phase(nbands), source=cmplx_0)

      ! --- k -> S k map, and its inverse ---------------------------------------
      do isym = 1, sym_full%nsym
         call sym_full%get_sym_operation_int_coord(isym, mrot, invmrot, trans, l_trs)
         do ik = 1, fi%kpts%nkpt
            rkpt = matmul(transpose(mrot), fi%kpts%bkf(:, ik))
            ikPrime = fi%kpts%get_nk(rkpt)
            if (ikPrime < 1 .or. ikPrime > fi%kpts%nkpt) call juDFT_error("S k is not on the k mesh.", calledby="dfpt_build_lambda")
            ikrot(ik, isym)          = ikPrime
            ikrot_inv(ikPrime, isym) = ik
         end do
      end do

      ! --- muffin-tin coefficients, once per k ---------------------------------
      ! These do not depend on the symmetry operation, so they are built here and
      ! reused for the whole group.
      !
      ! Both this loop and the k' loop below are spread over the ranks and summed
      ! back. Every rank needs the complete result -- the symmetry operation maps
      ! a k out of whatever subset a rank happens to own, and eig66_mpi serves any
      ! k to any rank through one-sided RMA anyway, so there is no ownership to
      ! respect. Each slot is written by exactly one rank, so a plain sum over the
      ! zero-initialised array reproduces it without double counting.
      call timestart("dfpt lambda: cmt")
      allocate (cmt_all(nbands, maxlmindx, fi%atoms%nat, fi%kpts%nkpt), source=cmplx_0)
      do ik = fmpi%irank + 1, fi%kpts%nkpt, fmpi%isize
         call lapw_k%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ik, fi%cell)
         call zMatK%init(fi%input%l_real, lapw_k%nv(jsp) + fi%atoms%nlotot, nbands)
         call read_eig(eig_id, ik, jsp, list=ev_list, neig=neig, eig=eig_dummy, zmat=zMatK)
         call cmt_from_z(fi, usdus, nococonv, lapw_k, zMatK, jsp, nbands, maxlmindx, cmt_all(:, :, :, ik))
         call zMatK%free()
      end do
#ifdef CPP_MPI
      if (fmpi%isize > 1) then
         call MPI_ALLREDUCE(MPI_IN_PLACE, cmt_all, size(cmt_all), MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
      end if
#endif
      call timestop("dfpt lambda: cmt")

      ! hs_int_direct distributes its column loop over fmpi%n_size; we want the
      ! whole block on this rank, so hand it a serial descriptor.
      fmpi_serial%n_rank = 0
      fmpi_serial%n_size = 1

      allocate (cmt_rot(nbands, maxlmindx, fi%atoms%nat))
      allocate (vpw_zero(size(stars%ustep)), source=cmplx_0)

      ! --- k' outer, so the interstitial overlap is built only once per k' -----
      ! At fixed isym the map k' -> ikrot_inv(k',isym) is a bijection, so spreading
      ! k' over the ranks leaves every lambda(:,:,ik,isym) written exactly once.
      do ikPrime = fmpi%irank + 1, fi%kpts%nkpt, fmpi%isize
         call lapw_kPrime%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ikPrime, fi%cell)
         nv_kPrime = lapw_kPrime%nv(jsp)

         call zMatKPrime%init(fi%input%l_real, nv_kPrime + fi%atoms%nlotot, nbands)
         call read_eig(eig_id, ikPrime, jsp, list=ev_list, neig=neig, eig=eig_dummy, zmat=zMatKPrime)

         ! Both the stored and the rotated state live on the LAPW basis at k',
         ! so this is simply the ordinary interstitial overlap matrix there.
         call smat_is%init(.false., nv_kPrime, nv_kPrime)
         call hmat_dummy%init(.false., nv_kPrime, nv_kPrime)
         call hs_int_direct(fmpi_serial, stars, fi%cell%bbmat, lapw_kPrime%gvec(:, :, jsp), lapw_kPrime%gvec(:, :, jsp), &
                            fi%kpts%bkf(:, ikPrime), fi%kpts%bkf(:, ikPrime), nv_kPrime, nv_kPrime, 0, 1, .true., .true., vpw_zero, hmat_dummy, smat_is)

         allocate (zKPrime_c(nv_kPrime, nbands), zRot_c(nv_kPrime, nbands), tmp_is(nv_kPrime, nbands))
         if (zMatKPrime%l_real) then
            zKPrime_c = cmplx(zMatKPrime%data_r(1:nv_kPrime, 1:nbands), 0.0)
         else
            zKPrime_c = zMatKPrime%data_c(1:nv_kPrime, 1:nbands)
         end if

         do isym = 1, sym_full%nsym
            ik = ikrot_inv(ikPrime, isym)

            call lapw_k%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ik, fi%cell)
            call zMatK%init(fi%input%l_real, lapw_k%nv(jsp) + fi%atoms%nlotot, nbands)
            call zMatRot%init(fi%input%l_real, nv_kPrime + fi%atoms%nlotot, nbands)
            call read_eig(eig_id, ik, jsp, list=ev_list, neig=neig, eig=eig_dummy, zmat=zMatK)

            ! Plane waves first: in the real (inversion-symmetric) case
            ! waveftrafo_gen_zmat picks a per-band phase so the rotated z stays
            ! real, and gen_cmt must be handed the same phase or the two halves
            ! of Lambda end up in different gauges.
            call waveftrafo_gen_zmat(zMatK, ik, isym, fi%kpts, sym_full, jsp, nbands, lapw_k, lapw_kPrime, zMatRot, c_phase)

            cmt_rot = cmplx_0
            call waveftrafo_gen_cmt(cmt_all(:, :, :, ik), c_phase, fi%input%l_real, ik, isym, fi%atoms, mpdata, hybinp, fi%kpts, sym_full, nbands, cmt_rot)

            ! --- muffin-tin block ---
            call mt_overlap(fi, olapmt, n_r, maxlmindx, maxn_r, cmt_all(:, :, :, ikPrime), cmt_rot, nbands, lambda(:, :, ik, isym))

            ! --- interstitial block ---
            if (zMatRot%l_real) then
               zRot_c = cmplx(zMatRot%data_r(1:nv_kPrime, 1:nbands), 0.0)
            else
               zRot_c = zMatRot%data_c(1:nv_kPrime, 1:nbands)
            end if

            call CPP_zgemm('N', 'N', nv_kPrime, nbands, nv_kPrime, cmplx_1, smat_is%data_c, nv_kPrime, zRot_c, nv_kPrime, cmplx_0, tmp_is, nv_kPrime)
            call CPP_zgemm('C', 'N', nbands, nbands, nv_kPrime, cmplx_1, zKPrime_c, nv_kPrime, tmp_is, nv_kPrime, cmplx_1, lambda(:, :, ik, isym), nbands)

            call zMatK%free()
            call zMatRot%free()
         end do

         deallocate (zKPrime_c, zRot_c, tmp_is)
         call smat_is%free()
         call hmat_dummy%free()
         call zMatKPrime%free()
      end do

#ifdef CPP_MPI
      ! One reduction per operation rather than one for the whole array: same
      ! total bytes, but the buffer stays at nbands^2 * nkpt instead of nsym
      ! times that.
      !
      ! isym is the trailing dimension, so lambda(1,1,1,isym) starts a contiguous
      ! run of nbands^2 * nkpt elements. Slicing on any earlier dimension would
      ! not be contiguous and must not be reduced this way.
      if (fmpi%isize > 1) then
         do isym = 1, sym_full%nsym
            call MPI_ALLREDUCE(MPI_IN_PLACE, lambda(1, 1, 1, isym), nbands*nbands*fi%kpts%nkpt, MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
         end do
      end if
#endif

      call timestop("dfpt lambda")

   end subroutine dfpt_build_lambda

   subroutine dfpt_check_lambda(lambda, tol)
      !! Unitarity check on every \(\Lambda(\boldsymbol{k},\mathcal{S})\) block.
      !!
      !! This is the cheapest and sharpest test in the chain, and it catches two
      !! independent failures:
      !!
      !! 1. An error in the muffin-tin rotation, in the radial overlap table or
      !!    in the interstitial half -- in isolation, before any q is involved.
      !! 2. A degenerate multiplet straddling an edge of the band window. The
      !!    window block of a unitary matrix is unitary precisely when no
      !!    multiplet is cut, so a violation here is exactly the condition under
      !!    which the unfolding is not well defined on the stored bands. It is
      !!    the same condition Wannier90's outer window has to satisfy.

      complex, intent(in) :: lambda(:, :, :, :)
      real,    intent(in) :: tol

      integer :: ik, isym, i, nb
      real    :: devMax
      complex, allocatable :: prod(:, :)

      nb = size(lambda, 1)
      allocate (prod(nb, nb))
      devMax = 0.0

      do isym = 1, size(lambda, 4)
         do ik = 1, size(lambda, 3)
            prod = matmul(conjg(transpose(lambda(:, :, ik, isym))), lambda(:, :, ik, isym))
            do i = 1, nb
               prod(i, i) = prod(i, i) - cmplx_1
            end do
            devMax = max(devMax, maxval(abs(prod)))
         end do
      end do

      write (oUnit, '(a,es12.4)') "dfpt_check_lambda: max |Lambda^dag Lambda - 1| = ", devMax
      if (devMax > tol) then
         call juDFT_warn("Lambda is not unitary to tolerance; the unfolded el-ph elements will be wrong.", calledby="dfpt_check_lambda")
      end if

   end subroutine dfpt_check_lambda

   subroutine dfpt_qmesh_full_bz(qvec_full, qvec_bz)
      !! Repackages the symmetry-reduced q mesh as a plain, symmetry-free mesh
      !! whose `nkpt`/`bk` are the *full* BZ list.
      !!
      !! Everything downstream of the unfolding (the forward Fourier transform in
      !! matrix_interpolation, the round-trip check) sums over `%nkpt` and reads
      !! `%bk`. Handing those routines this object lets them consume the full zone
      !! without a single change: they simply see a larger mesh that happens to
      !! carry no symmetry, which is exactly the situation they were written for.

      type(t_kpts), intent(in)  :: qvec_full
      type(t_kpts), intent(out) :: qvec_bz

      qvec_bz = qvec_full
      qvec_bz%nkpt = qvec_full%nkptf

      if (allocated(qvec_bz%bk)) deallocate (qvec_bz%bk)
      allocate (qvec_bz%bk(3, qvec_full%nkptf))
      qvec_bz%bk = qvec_full%bkf(:, 1:qvec_full%nkptf)

      if (allocated(qvec_bz%wtkpt)) deallocate (qvec_bz%wtkpt)
      allocate (qvec_bz%wtkpt(qvec_full%nkptf))
      qvec_bz%wtkpt = 1.0/real(qvec_full%nkptf)

   end subroutine dfpt_qmesh_full_bz

   subroutine dfpt_unfold_gmat(fi, sym_full, qvec_full, jsp, lambda, ikrot, gmatIBZ, gmatFull)
      !! Fills the full q Brillouin zone from the electron-phonon matrix elements
      !! computed on the irreducible wedge, by
      !! $$\underline{g}^{\kappa'\beta}(S\boldsymbol{k},S\boldsymbol{q})
      !!   =\Lambda(\boldsymbol{k}+\boldsymbol{q},\mathcal{S})
      !!    \Big[e^{-i(S\boldsymbol{q})\cdot\boldsymbol{R}_{S\kappa}}
      !!         \sum_\alpha S_{\beta\alpha}\,
      !!         \underline{g}^{\kappa\alpha}(\boldsymbol{k},\boldsymbol{q})\Big]
      !!    \Lambda^\dagger(\boldsymbol{k},\mathcal{S})$$
      !! which is Eq. (9) of README_symmetry.md. No eigenvector work happens here;
      !! `lambda` is supplied by the caller and is not recomputed.
      !!
      !! Conventions, all fixed by `get_sym_operation_int_coord`
      !! (types_sym.f90) and `mapatom` (mapatom.F90):
      !!
      !! * `invmrot` returned for operation `iop` is `sym%mrot(:,:,iop)` itself,
      !!   and `trans` is `sym%tau(:,iop)` -- precisely the pair that
      !!   `mapped_atom(iop,:)` was tabulated with. So
      !!   \(\boldsymbol{R}_{S\kappa}\) closes exactly in these variables.
      !! * For the time-reversal half (`isym > nop`) `invmrot` carries an extra
      !!   minus sign that belongs to k space, not to real space, and
      !!   `mapped_atom` is only tabulated up to `nop`. Both are undone here:
      !!   the real-space operation is `-invmrot` at index `isym - nop`.
      !! * Time reversal itself enters as
      !!   \(g^{\kappa\alpha}(-\boldsymbol{k},-\boldsymbol{q})
      !!     =[g^{\kappa\alpha}(\boldsymbol{k},\boldsymbol{q})]^{*}\).
      !!   Written in terms of the stored \(\boldsymbol{q}_f=-S\boldsymbol{q}\)
      !!   the phase keeps the same form, so only the source block is conjugated.

      use m_inv3

      type(t_fleurinput), intent(in)    :: fi
      type(t_sym),        intent(in)    :: sym_full
      type(t_kpts),       intent(in)    :: qvec_full
      integer,            intent(in)    :: jsp
      complex,            intent(in)    :: lambda(:, :, :, :)        !! (nb,nb,nkpt,nsym), same bands as gmat
      integer,            intent(in)    :: ikrot(:, :)               !! (nkpt,nsym): index of S k
      complex,            intent(in)    :: gmatIBZ(:, :, :, :, :, :) !! (nb,nb,nkpt,jspins,3nat,nkptIBZ)
      complex,            intent(inout) :: gmatFull(:, :, :, :, :, :)!! (nb,nb,nkpt,jspins,3nat,nkptf)

      integer :: iqf, iq, isym, iopAtom, ik, ikPrime, ikq, iAtom, jAtom
      integer :: alpha, beta, iPerturb, iPerturbPrime, nb
      integer :: mrot(3, 3), invmrot(3, 3), rotReal(3, 3), unitMat(3, 3)
      logical :: l_trs
      real    :: trans(3), q_f(3), q_0(3), rlat(3), phas, det
      real    :: invamat(3, 3), brot(3, 3)
      integer :: atomMap(fi%atoms%nat)
      complex :: phaseAtom(fi%atoms%nat)
      complex, allocatable :: gtmp(:, :), gsrc(:, :), lamK(:, :), lamKQ(:, :)

      nb = size(gmatIBZ, 1)

      if (size(lambda, 1) /= nb) then
         call juDFT_error("Lambda and gmat do not carry the same bands.", calledby="dfpt_unfold_gmat")
      end if
      if (size(gmatFull, 6) /= qvec_full%nkptf) then
         call juDFT_error("gmatFull q dimension is not the full BZ size.", calledby="dfpt_unfold_gmat")
      end if

      call timestart("dfpt unfold gmat")

      unitMat = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
      call inv3(fi%cell%amat, invamat, det)

      allocate (gtmp(nb, nb), gsrc(nb, nb), lamK(nb, nb), lamKQ(nb, nb))

      do iqf = 1, qvec_full%nkptf
         iq   = qvec_full%bkp(iqf)
         isym = qvec_full%bksym(iqf)
         call sym_full%get_sym_operation_int_coord(isym, mrot, invmrot, trans, l_trs)

         iopAtom = isym
         rotReal = invmrot
         if (l_trs) then
            iopAtom = isym - sym_full%nop
            rotReal = -invmrot
         end if

         q_f = qvec_full%bkf(:, iqf)
         q_0 = qvec_full%bk(:, iq)

         ! The identity is copied rather than round-tripped through Lambda, so
         ! the irreducible blocks stay bit-exact.
         if ((.not.l_trs) .and. all(rotReal == unitMat) .and. all(abs(trans) < 1e-9)) then
            gmatFull(:, :, :, jsp, :, iqf) = gmatIBZ(:, :, :, jsp, :, iq)
            cycle
         end if

         brot = matmul(fi%cell%amat, matmul(real(rotReal), invamat))

         do iAtom = 1, fi%atoms%nat
            jAtom = sym_full%mapped_atom(iopAtom, iAtom)
            atomMap(iAtom) = jAtom
            rlat = matmul(real(rotReal), fi%atoms%taual(:, iAtom)) + trans - fi%atoms%taual(:, jAtom)
            phas = -tpi_const*dot_product(q_f, rlat)
            phaseAtom(iAtom) = cmplx(cos(phas), sin(phas))
         end do

         do ik = 1, fi%kpts%nkpt
            ikPrime = ikrot(ik, isym)
            ikq     = fi%kpts%get_nk(fi%kpts%bkf(:, ik) + q_0)
            if (ikq < 1 .or. ikq > fi%kpts%nkpt) then
               call juDFT_error("k+q is not on the k mesh; the q mesh must be commensurate with it.", calledby="dfpt_unfold_gmat")
            end if

            lamK  = lambda(:, :, ik, isym)
            lamKQ = lambda(:, :, ikq, isym)

            do iAtom = 1, fi%atoms%nat
               jAtom = atomMap(iAtom)
               do beta = 1, 3
                  gtmp = cmplx_0
                  do alpha = 1, 3
                     iPerturb = alpha + 3*(iAtom - 1)
                     gsrc = gmatIBZ(:, :, ik, jsp, iPerturb, iq)
                     if (l_trs) gsrc = conjg(gsrc)
                     gtmp = gtmp + brot(beta, alpha)*gsrc
                  end do
                  gtmp = phaseAtom(iAtom)*gtmp

                  iPerturbPrime = beta + 3*(jAtom - 1)
                  gmatFull(:, :, ikPrime, jsp, iPerturbPrime, iqf) = matmul(lamKQ, matmul(gtmp, conjg(transpose(lamK))))
               end do
            end do
         end do
      end do

      call timestop("dfpt unfold gmat")

   end subroutine dfpt_unfold_gmat

   subroutine setup_rotation_data(fi, sym_full, usdus, jsp, hybinp, mpdata, olapmt, n_r, maxlmindx, maxn_r)
      !! Everything that depends neither on k nor on the symmetry operation: the
      !! atom map, translation vectors and Wigner d matrices on the full group,
      !! plus the radial overlap table in the t_abc radial ordering.

      use m_dwigner

      type(t_fleurinput),   intent(in)  :: fi
      type(t_sym),          intent(in)  :: sym_full
      type(t_usdus),        intent(in)  :: usdus
      integer,              intent(in)  :: jsp
      type(t_hybinp),       intent(out) :: hybinp
      type(t_mpdata),       intent(out) :: mpdata
      complex, allocatable, intent(out) :: olapmt(:, :, :, :)
      integer, allocatable, intent(out) :: n_r(:, :)
      integer,              intent(out) :: maxlmindx, maxn_r

      integer :: itype, l, m1, m2, isym, iisym, ilo, jlo, iOrd, jOrd, lmaxd
      integer :: n_l(0:fi%atoms%lmaxd), lo_ord(fi%atoms%nlod)

      lmaxd = fi%atoms%lmaxd

      ! t_hybinp%init only fills map/tvec/d_wgn2 for hybrid functionals, so call
      ! the generators directly. Both are cheap and self-contained.
      call hybinp%gen_map(fi%atoms, sym_full)

      ! Wigner d matrices, l = 0..lmaxd, including the time-reversal half.
      ! Copied from init_hybinp (types_hybinp.f90); d_wigner fills l >= 1 only,
      ! and the l = 0 slice is unity by definition.
      allocate (hybinp%d_wgn2(-lmaxd:lmaxd, -lmaxd:lmaxd, 0:lmaxd, sym_full%nsym))
      call d_wigner(sym_full%nop, sym_full%mrot, fi%cell%bmat, lmaxd, hybinp%d_wgn2(:, :, 1:, :sym_full%nop))
      hybinp%d_wgn2(:, :, 0, :) = 1

      do isym = sym_full%nop + 1, sym_full%nsym
         iisym = isym - sym_full%nop
         do l = 0, lmaxd
            do m2 = -l, l
               do m1 = -l, -1
                  hybinp%d_wgn2(m1, m2, l, isym)  = hybinp%d_wgn2(-m1, m2, l, iisym)*(-1)**m1
                  hybinp%d_wgn2(-m1, m2, l, isym) = hybinp%d_wgn2(m1, m2, l, iisym)*(-1)**m1
               end do
               hybinp%d_wgn2(0, m2, l, isym) = hybinp%d_wgn2(0, m2, l, iisym)
            end do
         end do
      end do

      ! waveftrafo_gen_cmt strides through cmt with mpdata%num_radfun_per_l, so
      ! that component has to be present even though nothing else of t_mpdata is.
      call mpdata%set_num_radfun_per_l(fi%atoms)

      allocate (n_r(0:lmaxd, fi%atoms%ntype))
      do itype = 1, fi%atoms%ntype
         n_r(:, itype) = fi%atoms%num_radial_functions_per_l(itype)
      end do
      maxn_r = maxval(n_r)

      maxlmindx = 0
      do itype = 1, fi%atoms%ntype
         maxlmindx = max(maxlmindx, sum([((2*l + 1)*n_r(l, itype), l=0, fi%atoms%lmax(itype))]))
      end do

      ! Radial overlap table. t_usdus already holds every integral we need, so
      ! the hybrid-only bas1/bas2/gridf machinery is not required:
      !   <u|u>       = 1        (u is normalised)
      !   <u|udot>    = 0        (by construction)
      !   <udot|udot> = ddn
      !   <u|ulo>     = uulon
      !   <udot|ulo>  = dulon
      !   <ulo|ulo'>  = uloulopn
      allocate (olapmt(maxn_r, maxn_r, 0:lmaxd, fi%atoms%ntype), source=cmplx_0)

      do itype = 1, fi%atoms%ntype
         do l = 0, fi%atoms%lmax(itype)
            olapmt(1, 1, l, itype) = cmplx(1.0, 0.0)
            olapmt(2, 2, l, itype) = cmplx(usdus%ddn(l, itype, jsp), 0.0)
         end do

         ! Reproduce exactly the iOrd assignment of t_abc%calc_abc: the local
         ! orbitals of a given l are numbered in the order they appear in llo.
         n_l    = 2
         lo_ord = 0
         do ilo = 1, fi%atoms%nlo(itype)
            l           = fi%atoms%llo(ilo, itype)
            n_l(l)      = n_l(l) + 1
            lo_ord(ilo) = n_l(l)
         end do

         do ilo = 1, fi%atoms%nlo(itype)
            l    = fi%atoms%llo(ilo, itype)
            iOrd = lo_ord(ilo)

            olapmt(1, iOrd, l, itype) = cmplx(usdus%uulon(ilo, itype, jsp), 0.0)
            olapmt(iOrd, 1, l, itype) = olapmt(1, iOrd, l, itype)
            olapmt(2, iOrd, l, itype) = cmplx(usdus%dulon(ilo, itype, jsp), 0.0)
            olapmt(iOrd, 2, l, itype) = olapmt(2, iOrd, l, itype)

            do jlo = 1, fi%atoms%nlo(itype)
               if (fi%atoms%llo(jlo, itype) /= l) cycle
               jOrd = lo_ord(jlo)
               olapmt(iOrd, jOrd, l, itype) = cmplx(usdus%uloulopn(ilo, jlo, itype, jsp), 0.0)
            end do
         end do
      end do

   end subroutine setup_rotation_data

   subroutine cmt_from_z(fi, usdus, nococonv, lapw, zMat, jsp, nbands, maxlmindx, cmt)
      !! Expands an eigenvector in muffin-tin coefficients, in the layout used by
      !! the hybrid code: cmt(band, lmindx, atom) with lmindx running l, then m,
      !! then the radial-function index (u, udot, local orbitals).
      !!
      !! This is the core of calc_cmt (m_calc_cmt) with two deliberate
      !! differences: the rotation branch is dropped, because the caller applies
      !! `waveftrafo_gen_cmt` itself; and the inversion-pair symmetrization is
      !! dropped, because it is referenced to a specific k point and would not
      !! cancel between the two cmt sets entering Lambda.

      use m_types_abc

      type(t_fleurinput), intent(in)  :: fi
      type(t_usdus),      intent(in)  :: usdus
      type(t_nococonv),   intent(in)  :: nococonv
      type(t_lapw),       intent(in)  :: lapw
      type(t_mat),        intent(in)  :: zMat
      integer,            intent(in)  :: jsp, nbands, maxlmindx
      complex,            intent(out) :: cmt(:, :, :)

      type(t_abc) :: abc
      integer     :: itype, na, iatom, indx, l, ll, m, lm, i
      complex     :: cdum

      cmt = cmplx_0

      do itype = 1, fi%atoms%ntype
         call abc%init(fi%input, fi%atoms, nbands, itype)
         call abc%calc_abc(fi%input, fi%atoms, fi%sym, fi%cell, lapw, nbands, usdus, fi%noco, nococonv, jsp, itype, zMat)

         ! No rot_to_unrotated here: the DFPT run has nop == 1, so every atom is
         ! its own representative and the local-frame rotation is the identity.

         do na = 1, fi%atoms%neq(itype)
            iatom = fi%atoms%firstAtom(itype) + na - 1
            indx = 0
            do l = 0, fi%atoms%lmax(itype)
               ll   = l*(l + 1)
               cdum = ImagUnit**l
               do m = -l, l
                  lm = ll + m
                  do i = 1, abc%n_r(l)
                     indx = indx + 1
                     cmt(:, indx, iatom) = cdum*abc%cof(:, lm, i, na)
                  end do
               end do
            end do
         end do
      end do

   end subroutine cmt_from_z

   subroutine mt_overlap(fi, olapmt, n_r, maxlmindx, maxn_r, cmt_bra, cmt_ket, nbands, lam)
      !! Muffin-tin part of Lambda,
      !! $$\sum_{a}\sum_{lm}\sum_{ij}
      !!   c^{*}_{p,(lm\,i),a}\,\langle u_i|u_j\rangle_{l}\,c_{n,(lm\,j),a}.$$
      !!
      !! Same contraction as the wavefolap block of symm_hf.F90, but with an
      !! independent bra and ket instead of the same cmt twice.

      type(t_fleurinput), intent(in)    :: fi
      complex,            intent(in)    :: olapmt(:, :, 0:, :)
      integer,            intent(in)    :: n_r(0:, :)
      integer,            intent(in)    :: maxlmindx, maxn_r, nbands
      complex,            intent(in)    :: cmt_bra(:, :, :), cmt_ket(:, :, :)
      complex,            intent(inout) :: lam(:, :)

      complex, allocatable :: braT(:, :), ketT(:, :), tmp(:, :)
      integer :: iatom, itype, l, m, nn, lm

      allocate (braT(maxlmindx, nbands), ketT(maxlmindx, nbands), tmp(maxn_r, nbands))

      do iatom = 1, fi%atoms%nat
         itype = fi%atoms%itype(iatom)
         braT  = transpose(cmt_bra(:, :, iatom))
         ketT  = transpose(cmt_ket(:, :, iatom))

         lm = 0
         do l = 0, fi%atoms%lmax(itype)
            nn = n_r(l, itype)
            do m = -l, l
               call CPP_zgemm('N', 'N', nn, nbands, nn, cmplx_1, olapmt(1, 1, l, itype), maxn_r, ketT(lm + 1, 1), maxlmindx, cmplx_0, tmp, maxn_r)
               call CPP_zgemm('C', 'N', nbands, nbands, nn, cmplx_1, braT(lm + 1, 1), maxlmindx, tmp, maxn_r, cmplx_1, lam, nbands)
               lm = lm + nn
            end do
         end do
      end do

   end subroutine mt_overlap

end module m_dfpt_lambda
