!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

module m_dfpt_write_epw
   !! Write the real-space (Wannier) electron-phonon quantities computed by the
   !! FLEUR DFPT/Wannier interpolation as EPW (Quantum ESPRESSO) restart files, so
   !! that EPW can restart with `epwread = .true., epwwrite = .false.` and run the
   !! (isotropic/anisotropic) Eliashberg machinery without redoing the QE
   !! nscf + Wannierization + coarse e-ph steps.
   !!
   !! Files produced (all in the run directory):
   !!   crystal.fmt, epwdata.fmt, vmedata.fmt, dmedata.fmt, wigner.fmt,
   !!   <prefix>.ukk, <prefix>.epmatwp
   !!
   !! Conventions / unit conversions (FLEUR is Hartree/bohr native, EPW is Ry/bohr):
   !!   * energies  : E_Ry  = 2 E_Htr           -> chw *2, epmatwp *2
   !!   * dyn-matrix : rdw block(a,b) = 4 * amu_ry * sqrt(m_a m_b) * D_norm(R)
   !!                  so that EPW's rdw/sqrt(amass_a amass_b) (amass in amu_ry)
   !!                  yields eigenvalues omega^2 in Ry^2 (= 4 omega^2_Htr).
   !!   * masses    : amass = m[amu] * amu_ry   (QE internal Rydberg mass unit)
   !!   * FT sign   : EPW builds the real-space objects with e^{-i k.R} (see
   !!                 hambloch2wan/dynbloch2wan/ephbloch2wanp) and interpolates back
   !!                 with e^{+i k.R}. FLEUR uses the opposite sign (e^{+i k.R}
   !!                 forward, e^{-i k.R} backward). Because EPW re-interpolates from
   !!                 the data we write, we make EPW reproduce FLEUR's result exactly
   !!                 by *negating the WS lattice vectors* written to wigner.fmt while
   !!                 keeping the stored tensor values. Then, per WS vector i,
   !!                 EPW's  sum_i e^{+i k.(-R_i)} T_i  ==  FLEUR's e^{-i k.R_i} T_i.
   !!   * degeneracy: EPW stores integer ndegen and divides by it; FLEUR stores a
   !!                 real weight = 1/ndegen. So ndegen = nint(1/weight).
   !!
   !! velocity/position/dipole matrix elements (cvmew/crrw/cdmew) and the Wannier
   !! centers are written as zeros: they are read but unused by isotropic/standard
   !! anisotropic Eliashberg runs. (In the QE reference files crrw is itself zero.)

   use m_juDFT
   use m_types
   use m_constants
   use m_matrix_interpolation,  only: t_wann_ft
   use m_dfpt_dynmat_fourier,   only: ft_dyn, ft_dyn_direct, unfold_grid, build_ws_ft

   implicit none
   private

   public :: write_epw_restart_files

   ! QE internal Rydberg atomic mass unit: AMU_SI / ELECTRONMASS_SI / 2
   real, parameter :: amu_ry = 911.4442421
   ! Hartree -> Rydberg energy conversion
   real, parameter :: htr2ry = 2.0

contains

   subroutine write_epw_restart_files(fi, results, dynMats, U_full, ftRealspace)
      !! Orchestrate writing all EPW restart files. Call on rank 0 only.
      type(t_fleurinput), intent(in) :: fi
      type(t_results),    intent(in) :: results
      complex,            intent(in) :: dynMats(:,:,:)       ! (3nat,3nat,nq_coarse) mass-normalized D(q)
      complex,            intent(in) :: U_full(:,:,:,:)      ! (num_bands,num_wann,nkptf,jspins)
      type(t_wann_ft),    intent(in) :: ftRealspace(:,:)     ! (3*nat, jspins) elph real-space tensors

      integer :: nbndsub, nmodes, nat, ispin
      integer :: nrr_k, nrr_q, nrr_g
      integer :: ft_lim_k(2,3)
      complex, allocatable :: chw(:,:,:)          ! (nbndsub,nbndsub,nrr_k)   electron H, real space [Ry]
      complex, allocatable :: rdw(:,:,:)          ! (nmodes,nmodes,nrr_q)     dyn-matrix, real space
      complex, allocatable :: chw_cube(:,:,:,:,:) ! dense electron-H cube on the k-mesh box
      integer, allocatable :: irvec_q(:,:)        ! (3,nrr_q) phonon WS vectors (own q-mesh)
      real,    allocatable :: weight_q(:)         ! (nrr_q)   phonon WS weights (= 1/ndegen)
      character(len=:), allocatable :: prefix

      nat     = fi%atoms%nat
      nmodes  = 3*nat
      nbndsub = fi%wannierlib%num_wann
      ispin   = 1
      prefix  = trim(adjustl(fi%dfpt%epw_prefix))

      if (fi%input%jspins /= 1) &
         call juDFT_warn("EPW restart export only writes spin channel 1 (jspins>1 not supported yet).", &
                         calledby="write_epw_restart_files")

      ! Electron (R_e) and e-ph (R_g) WS lists come from the elph tensor. Note that
      ! wannier_matrixq_forward builds that tensor with fi%kpts as BOTH k- and q-mesh,
      ! so its R_p axis (ws_q) lives on the k-mesh box -> this is EPW's irvec_g / nrr_g.
      ! The dynamical matrix (rdw) is instead handled on its OWN coarse phonon q-mesh
      ! (fi%dfpt%qvec), which may differ from the k-mesh. EPW keeps irvec_q and irvec_g
      ! as independent lists, so nrr_q and nrr_g are allowed to differ.
      nrr_k = ftRealspace(1,ispin)%ws_k_nNZ
      nrr_g = ftRealspace(1,ispin)%ws_q_nNZ
      ft_lim_k = ftRealspace(1,ispin)%ft_lim_k

      ! --- electron Hamiltonian real space (R_e) ---
      call build_chw_cube(fi, results, U_full, ispin, ft_lim_k, nbndsub, chw_cube)
      call unpack_cube(chw_cube, ftRealspace(1,ispin)%ws_k_indStored, nrr_k, htr2ry, chw)

      ! --- dynamical matrix real space (R_q) on the coarse phonon q-mesh ---
      call build_rdw(fi, dynMats, nmodes, rdw, irvec_q, weight_q, nrr_q)

      write (oUnit,'(/,a)') "Writing EPW (epwread) restart files:"
      write (oUnit,'(a,5(1x,i0))') "   nbndsub, nmodes, nrr_k, nrr_q, nrr_g =", nbndsub, nmodes, nrr_k, nrr_q, nrr_g
      write (oUnit,'(a,a)')        "   prefix = ", prefix

      ! --- write files ---
      call write_crystal_fmt(fi, results, nbndsub, nmodes)
      call write_epwdata_fmt(results, nbndsub, nmodes, nrr_k, nrr_q, nrr_g, nat, chw, rdw)
      call write_vme_dme_fmt(nbndsub, nrr_k)
      call write_wigner_fmt(fi, ftRealspace(1,ispin), irvec_q, weight_q, nrr_k, nrr_q, nrr_g)
      call write_ukk_fmt(fi, prefix)
      call write_epmatwp(fi, ftRealspace(:,ispin), nbndsub, nmodes, nrr_k, nrr_g, prefix)

      write (oUnit,'(a)') "   ... EPW restart files written."

   end subroutine write_epw_restart_files

   !---------------------------------------------------------------------------

   subroutine unpack_cube(cube, indStored, nNZ, scal, list)
      !! Copy the dense real-space cube values at the WS folded indices into a
      !! flat (m,n,ir) list, scaled by `scal`.
      complex,              intent(in)  :: cube(:,:,0:,0:,0:)
      integer,              intent(in)  :: indStored(:,:)   ! (3,nNZ) 0-based folded indices
      integer,              intent(in)  :: nNZ
      real,                 intent(in)  :: scal
      complex, allocatable, intent(out) :: list(:,:,:)

      integer :: ir, m, n

      m = size(cube,1); n = size(cube,2)
      allocate(list(m,n,nNZ))
      do ir = 1, nNZ
         list(:,:,ir) = scal * cube(:,:, indStored(1,ir), indStored(2,ir), indStored(3,ir))
      end do
   end subroutine unpack_cube

   subroutine build_chw_cube(fi, results, U_full, ispin, ft_lim, nwann, cube)
      !! Real-space Wannier-gauge electron Hamiltonian on the coarse k-mesh box.
      !! Mirrors wannier_matrix_interpolate's forward transform (m_matrix_interpolation).
      type(t_fleurinput), intent(in)  :: fi
      type(t_results),    intent(in)  :: results
      complex,            intent(in)  :: U_full(:,:,:,:)
      integer,            intent(in)  :: ispin, ft_lim(2,3), nwann
      complex, allocatable, intent(out) :: cube(:,:,:,:,:)

      integer :: num_bands, ikpt, ib, nk1, nk2, nk3, ngrid
      complex, allocatable :: H_bloch(:,:), matRot(:,:), fft_grid(:,:,:)

      num_bands = fi%wannierlib%max_band - fi%wannierlib%min_band + 1
      ! size the box from ft_lim so the folding matches the WS indStored exactly
      nk1 = ft_lim(2,1) - ft_lim(1,1) + 1
      nk2 = ft_lim(2,2) - ft_lim(1,2) + 1
      nk3 = ft_lim(2,3) - ft_lim(1,3) + 1
      ngrid = nk1*nk2*nk3

      allocate(H_bloch(num_bands,num_bands), matRot(nwann,nwann))
      allocate(fft_grid(nwann,nwann,ngrid))
      allocate(cube(nwann,nwann,0:nk1-1,0:nk2-1,0:nk3-1))
      fft_grid = cmplx(0.0,0.0)

      do ikpt = 1, fi%kpts%nkpt
         H_bloch = cmplx(0.0,0.0)
         do ib = 1, num_bands
            H_bloch(ib,ib) = cmplx(results%eig(fi%wannierlib%min_band + ib - 1, fi%kpts%bkp(ikpt), ispin), 0.0)
         end do
         ! Wannier gauge: U^dagger H U
         matRot = matmul(conjg(transpose(U_full(:,:,ikpt,ispin))), matmul(H_bloch, U_full(:,:,ikpt,ispin)))
         call ft_dyn_direct(ft_lim, 1, fi%kpts%bk(:,ikpt), matRot, fft_grid)
      end do
      fft_grid = fft_grid / fi%kpts%nkpt
      call unfold_grid(ft_lim, fft_grid, cube)
   end subroutine build_chw_cube

   subroutine build_rdw(fi, dynMats, nmodes, rdw, Rvecs_q, weightNZ_q, nrr_q)
      !! Real-space dynamical matrix (rdw) on the coarse PHONON q-mesh (fi%dfpt%qvec),
      !! together with the phonon Wigner-Seitz list (irvec_q, ndegen_q = 1/weight).
      !! This mesh is independent of the electron/e-ph (k-mesh) WS lists.
      !! Scaling into EPW's convention (see module header):
      !!   rdw block(a,b) = 4 * amu_ry * sqrt(m_a m_b) * D_norm(R).
      type(t_fleurinput),   intent(in)  :: fi
      complex,              intent(in)  :: dynMats(:,:,:)
      integer,              intent(in)  :: nmodes
      complex, allocatable, intent(out) :: rdw(:,:,:)        ! (nmodes,nmodes,nrr_q)
      integer, allocatable, intent(out) :: Rvecs_q(:,:)      ! (3,nrr_q) WS vectors (un-negated)
      real,    allocatable, intent(out) :: weightNZ_q(:)     ! (nrr_q)   = 1/ndegen
      integer,              intent(out) :: nrr_q

      integer :: nq(3), ft_lim(2,3), bigBox(2,3), boxSize, ix, iy, iz, iGrid, ia, ja, i0, j0
      real    :: fac
      integer, allocatable :: supercellR(:,:), indStored_q(:,:)
      real,    allocatable :: FTweight(:)
      complex, allocatable :: cube(:,:,:,:,:), dyn_mat_q_full(:,:,:)
      type(t_cell) :: cellLocal

      nq = fi%dfpt%qvec%nkpt3(:)

      ! Fourier box and Wigner-Seitz weights of the coarse phonon q-mesh supercell
      ! (mirrors interpolate_dynmat / wannier_matrixq_forward: big box +-2N, always WS).
      ft_lim(2,:) = nq(:)/2
      ft_lim(1,:) = ft_lim(2,:) - nq(:) + 1
      bigBox(2,:) =  2*nq(:)
      bigBox(1,:) = -2*nq(:)
      boxSize = (4*nq(1)+1) * (4*nq(2)+1) * (4*nq(3)+1)
      allocate(supercellR(3,boxSize), FTweight(boxSize))
      FTweight = 0.0
      iGrid = 1
      do iz = bigBox(1,3), bigBox(2,3)
         do iy = bigBox(1,2), bigBox(2,2)
            do ix = bigBox(1,1), bigBox(2,1)
               supercellR(:,iGrid) = [ix,iy,iz]
               iGrid = iGrid + 1
            end do
         end do
      end do
      cellLocal = fi%cell
      call cellLocal%calculate_WSweight(supercellR, FTweight, scaleSupercell=nq)
      call build_ws_ft(ft_lim, bigBox, FTweight, nrr_q, Rvecs_q, indStored_q, weightNZ_q)

      ! D_norm(R) = (1/Nq) sum_q e^{+i q.R} D(q)   (symmetry-unfolded to the full BZ)
      allocate(cube(nmodes,nmodes,0:nq(1)-1,0:nq(2)-1,0:nq(3)-1))
      call ft_dyn(fi%atoms, fi%dfpt%qvec, fi%sym, ft_lim, fi%cell%amat, dynMats, cube, dyn_mat_q_full)

      ! per-(atom,atom) block mass/energy scaling (see module header)
      do ia = 1, fi%atoms%nat
         i0 = 3*(ia-1)
         do ja = 1, fi%atoms%nat
            j0 = 3*(ja-1)
            fac = htr2ry**2 * amu_ry * sqrt(atomicMasses_const(fi%atoms%nz(ia)) * atomicMasses_const(fi%atoms%nz(ja)))
            cube(i0+1:i0+3, j0+1:j0+3, :,:,:) = fac * cube(i0+1:i0+3, j0+1:j0+3, :,:,:)
         end do
      end do

      call unpack_cube(cube, indStored_q, nrr_q, 1.0, rdw)
   end subroutine build_rdw

   !---------------------------------------------------------------------------

   subroutine write_crystal_fmt(fi, results, nbndsub, nmodes)
      !! crystal.fmt -- exact record order of EPW's epw_write (io.f90:219-235).
      type(t_fleurinput), intent(in) :: fi
      type(t_results),    intent(in) :: results
      integer,            intent(in) :: nbndsub, nmodes

      integer, parameter :: ntypx = 10           ! QE Modules/parameters.f90
      integer :: iu, ia, it
      real    :: alat, at(3,3), bg(3,3), amass(ntypx)
      real, allocatable :: tau(:,:), w_centers(:,:)

      alat = norm2(fi%cell%amat(:,1))             ! any consistent scale works (w_centers=0)
      at   = fi%cell%amat / alat                  ! lattice vectors in columns, alat units
      do it = 1, 3
         bg(:,it) = fi%cell%bmat(it,:) * alat / tpi_const   ! recip vectors in columns, 2pi/alat units
      end do
      allocate(tau(3,fi%atoms%nat), w_centers(3,nbndsub))
      tau = fi%atoms%pos / alat
      w_centers = 0.0
      amass = 0.0
      do it = 1, fi%atoms%ntype
         if (it <= ntypx) amass(it) = atomicMasses_const(fi%atoms%nz(it)) * amu_ry
      end do

      open(newunit=iu, file="crystal.fmt", status='replace', action='write', form='formatted')
      write(iu,*) fi%atoms%nat
      write(iu,*) nmodes
      write(iu,*) fi%input%zelec, fi%wannierlib%min_band - 1     ! nelec, nbndskip
      write(iu,*) at
      write(iu,*) bg
      write(iu,*) fi%cell%omtil                                  ! omega [bohr^3]
      write(iu,*) alat                                           ! [bohr]
      write(iu,*) tau                                            ! (3,nat) alat units
      write(iu,*) amass                                          ! (ntypx) amu_ry
      write(iu,*) (fi%atoms%itype(ia), ia=1,fi%atoms%nat)        ! ityp
      write(iu,*) .false.                                        ! noncolin
      write(iu,*) .false.                                        ! do_cutoff_2D_epw
      write(iu,*) w_centers                                      ! (3,nbndsub) zeros
      write(iu,*) 0.0                                            ! L (2D cutoff length)
      write(iu,*) htr2ry * fi%input%tkb                          ! degauss [Ry]
      write(iu,*) 0                                              ! ngauss
      write(iu,*) .false.                                        ! lda_plus_u
      close(iu)
   end subroutine write_crystal_fmt

   subroutine write_epwdata_fmt(results, nbndsub, nmodes, nrr_k, nrr_q, nrr_g, nat, chw, rdw)
      !! epwdata.fmt -- EPW's epw_write (io.f90:237-273).
      type(t_results), intent(in) :: results
      integer,         intent(in) :: nbndsub, nmodes, nrr_k, nrr_q, nrr_g, nat
      complex,         intent(in) :: chw(:,:,:), rdw(:,:,:)

      integer :: iu, ibnd, jbnd, irk, imode, jmode, irq
      real    :: zstar(3,3,nat), epsi(3,3)   ! EPW stores these as REAL (18 reals in epwdata.fmt)

      zstar = 0.0     ! Born effective charges (0 for non-polar metal)
      epsi  = 0.0     ! dielectric tensor

      open(newunit=iu, file="epwdata.fmt", status='replace', action='write', form='formatted')
      write(iu,*) htr2ry * results%ef                            ! Fermi energy [Ry]
      write(iu,*) nbndsub, nrr_k, nmodes, nrr_q, nrr_g
      write(iu,*) zstar, epsi
      do ibnd = 1, nbndsub
         do jbnd = 1, nbndsub
            do irk = 1, nrr_k
               write(iu,*) chw(ibnd, jbnd, irk)
            end do
         end do
      end do
      do imode = 1, nmodes
         do jmode = 1, nmodes
            do irq = 1, nrr_q
               write(iu,*) rdw(imode, jmode, irq)
            end do
         end do
      end do
      close(iu)
   end subroutine write_epwdata_fmt

   subroutine write_vme_dme_fmt(nbndsub, nrr_k)
      !! vmedata.fmt / dmedata.fmt -- velocity/position/dipole matrix elements.
      !! Written as zeros (unused by Eliashberg); exact interleaving of io.f90:246-251.
      integer, intent(in) :: nbndsub, nrr_k

      integer :: iuv, iud, ibnd, jbnd, irk, ipol
      complex, parameter :: c0 = (0.0,0.0)

      open(newunit=iuv, file="vmedata.fmt", status='replace', action='write', form='formatted')
      open(newunit=iud, file="dmedata.fmt", status='replace', action='write', form='formatted')
      do ibnd = 1, nbndsub
         do jbnd = 1, nbndsub
            do irk = 1, nrr_k
               do ipol = 1, 3
                  write(iuv,*) c0        ! cvmew(ipol,...)
                  write(iuv,*) c0        ! crrw(ipol,...)
                  write(iud,*) c0        ! cdmew(ipol,...)
               end do
            end do
         end do
      end do
      close(iuv)
      close(iud)
   end subroutine write_vme_dme_fmt

   subroutine write_wigner_fmt(fi, ft, Rvecs_q, weight_q, nrr_k, nrr_q, nrr_g)
      !! wigner.fmt -- EPW's epw_write_ws_data (io.f90:669-737), dims=dims2=1 case.
      !! WS lattice vectors are NEGATED to convert FLEUR's FT sign to EPW's (see header).
      !!   R_e (electron) and R_g (e-ph) blocks use the elph-tensor WS lists (k-mesh box);
      !!   R_q (phonon)   block uses the independent dynamical-matrix WS list (q-mesh).
      type(t_fleurinput), intent(in) :: fi
      type(t_wann_ft),    intent(in) :: ft
      integer,            intent(in) :: Rvecs_q(:,:)     ! (3,nrr_q) phonon WS vectors (un-negated)
      real,               intent(in) :: weight_q(:)      ! (nrr_q)   = 1/ndegen
      integer,            intent(in) :: nrr_k, nrr_q, nrr_g

      integer :: iu, ir, irvec(3)
      real    :: alat

      alat = norm2(fi%cell%amat(:,1))

      open(newunit=iu, file="wigner.fmt", status='replace', action='write', form='formatted')
1     format(1000(I0, :, 1X))
2     format(3I6, ES26.17E3)
      write(iu,1) nrr_k, nrr_q, nrr_g, 1, 1
      ! electron block (R_e) -- elph-tensor k-mesh WS list
      do ir = 1, nrr_k
         irvec = -ft%ws_k_Rvecs(:,ir)
         write(iu,2) irvec, ws_len(fi%cell%amat, ft%ws_k_Rvecs(:,ir), alat)
         write(iu,1) nint(1.0/ft%ws_k_weightNZ(ir))
      end do
      ! phonon block (R_q) -- dynamical-matrix q-mesh WS list
      do ir = 1, nrr_q
         irvec = -Rvecs_q(:,ir)
         write(iu,2) irvec, ws_len(fi%cell%amat, Rvecs_q(:,ir), alat)
         write(iu,1) nint(1.0/weight_q(ir))
      end do
      ! e-ph block (R_g) -- elph-tensor k-mesh WS list (R_p axis of the tensor)
      do ir = 1, nrr_g
         irvec = -ft%ws_q_Rvecs(:,ir)
         write(iu,2) irvec, ws_len(fi%cell%amat, ft%ws_q_Rvecs(:,ir), alat)
         write(iu,1) nint(1.0/ft%ws_q_weightNZ(ir))
      end do
      close(iu)
   end subroutine write_wigner_fmt

   real function ws_len(amat, irvec, alat)
      !! |amat.irvec| / alat  == length of the WS vector in alat units (matches EPW wslen).
      real,    intent(in) :: amat(3,3)
      integer, intent(in) :: irvec(3)
      real,    intent(in) :: alat
      real :: rc(3)
      rc = matmul(amat, real(irvec))
      ws_len = norm2(rc) / alat
   end function ws_len

   subroutine write_ukk_fmt(fi, prefix)
      !! <prefix>.ukk -- band-manifold bookkeeping file. EPW's restart path
      !! (epwread=.true., wannierize=.false.) still calls loadbm() (io.f90:3196),
      !! which opens this file and reads ONLY:
      !!    READ(iu,*) nbndep, nbndskip
      !!    DO ibnd = 1, nbndep;  READ(iu,*) ibndkept(ibnd);  END DO
      !! and then closes it. The U(k) rotation, lwindow, excluded_band and Wannier
      !! centres that a full write_filukk stores are NOT read on restart (the coarse
      !! unfolding that would use them is skipped), so we write only the header +
      !! kept-band list. FLEUR's disentanglement window is contiguous, so:
      !!    nbndep   = max_band - min_band + 1
      !!    nbndskip = min_band - 1
      !!    ibndkept = [min_band, ..., max_band]
      type(t_fleurinput), intent(in) :: fi
      character(len=*),   intent(in) :: prefix

      integer :: iu, ibnd, nbndep, nbndskip

      nbndep   = fi%wannierlib%max_band - fi%wannierlib%min_band + 1
      nbndskip = fi%wannierlib%min_band - 1

      open(newunit=iu, file=trim(prefix)//".ukk", status='replace', action='write', form='formatted')
      write(iu,*) nbndep, nbndskip
      do ibnd = 1, nbndep
         write(iu,*) fi%wannierlib%min_band + ibnd - 1
      end do
      close(iu)
   end subroutine write_ukk_fmt

   subroutine write_epmatwp(fi, ft_spin, nbndsub, nmodes, nrr_k, nrr_g, prefix)
      !! <prefix>.epmatwp -- direct-access unformatted, one record per irg,
      !! record = epmatwp(nbndsub,nbndsub,nrr_k,nmodes) (EPW serial branch, io.f90:338-351).
      type(t_fleurinput), intent(in) :: fi
      type(t_wann_ft),    intent(in) :: ft_spin(:)          ! ftRealspace(:,ispin), size 3*nat
      integer,            intent(in) :: nbndsub, nmodes, nrr_k, nrr_g
      character(len=*),   intent(in) :: prefix

      integer :: iu, irg, imode, irk, lrepmatw, direct_io_factor, unf_recl
      integer :: ke(3), kg(3)
      real    :: dummy
      complex, allocatable :: rec4d(:,:,:,:)

      lrepmatw = 2 * nbndsub * nbndsub * nrr_k * nmodes
      inquire(iolength=direct_io_factor) dummy
      unf_recl = direct_io_factor * lrepmatw

      allocate(rec4d(nbndsub,nbndsub,nrr_k,nmodes))
      open(newunit=iu, file=trim(prefix)//".epmatwp", status='replace', action='write', &
           form='unformatted', access='direct', recl=unf_recl)
      do irg = 1, nrr_g
         kg = ft_spin(1)%ws_q_indStored(:,irg)
         do imode = 1, nmodes
            do irk = 1, nrr_k
               ke = ft_spin(1)%ws_k_indStored(:,irk)
               rec4d(:,:,irk,imode) = htr2ry * &
                  ft_spin(imode)%matWannier(:,:, ke(1),ke(2),ke(3), kg(1),kg(2),kg(3))
            end do
         end do
         write(iu, rec=irg) rec4d
      end do
      close(iu)
   end subroutine write_epmatwp

end module m_dfpt_write_epw
