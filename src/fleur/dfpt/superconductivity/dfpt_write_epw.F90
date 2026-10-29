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
   !!   * energies  : E_Ry  = 2 E_Htr           -> hamWannier *2, epmatwp *2
   !!   * dyn-matrix : fcm block(a,b) = 4 * amu_ry * sqrt(m_a m_b) * D_norm(R)
   !!                  so that EPW's rdw/sqrt(amass_a amass_b) (amass in amu_ry)
   !!                  yields eigenvalues omega^2 in Ry^2 (= 4 omega^2_Htr).
   !!   * masses    : amass = m[amu] * amu_ry   (QE internal Rydberg mass unit)
   !!   * WS lists  : on restart EPW does not read wigner.fmt; it regenerates its WS
   !!                 lists (wigner_seitzkq, sort=.true.) and reads the data in that
   !!                 order. Every R-indexed array is therefore written in EPW's order.
   !!                 FLEUR interpolates back with e^{-i k.R}, EPW with e^{+i k.R}, so
   !!                 the FLEUR value at R is stored at EPW's vector -R.
   !!                 wigner.fmt is currently not read by EPW; it is written in EPW's
   !!                 format and order (use_ws = .false.) for reference.
   !!   * degeneracy: EPW stores integer ndegen and divides by it; FLEUR stores a
   !!                 real weight = 1/ndegen. So ndegen = nint(1/weight).
   !!
   !! velocity/position/dipole matrix elements (cvmew/crrw/cdmew) and the Wannier
   !! centers are written as zeros: they are read but unused by isotropic/standard
   !! anisotropic Eliashberg runs. (In the QE reference files crrw is itself zero.)

   use m_juDFT
   use m_constants
   use m_matrix_interpolation,  only: t_wann_ft
   use m_dfpt_dynmat_fourier,   only: ft_dyn, ft_ws_forward, build_ws_ft
   use m_dfpt_sumrules,         only: dfpt_born_huang
   use m_sort,                  only: sort
   use m_types_cell
   use m_types_fleurinput
   use m_types_kpts
   use m_types_misc
   use m_types_sym

   implicit none
   private

   public :: write_epw_restart_files

   ! QE internal Rydberg atomic mass unit: AMU_SI / ELECTRONMASS_SI / 2
   real, parameter :: amu_ry = 911.4442421
   ! Hartree -> Rydberg energy conversion
   real, parameter :: htr2ry = 2.0

contains

   subroutine write_epw_restart_files(fi, results, dynMats, U_full, gmatRealspace, qpts_dyn, sym_dyn)
      !! Orchestrate writing all EPW restart files. Call on rank 0 only.
      !!
      !! qpts_dyn/sym_dyn describe the q mesh dynMats lives on: ft_dyn unfolds
      !! dynMats from qpts_dyn%bk onto qpts_dyn%bkf using sym_dyn, and
      !! qpts_dyn%nkpt3 sets the phonon Fourier box.
      type(t_fleurinput), intent(in) :: fi
      type(t_results),    intent(in) :: results
      complex,            intent(in) :: dynMats(:,:,:)       ! (3nat,3nat,nq_coarse) mass-normalized D(q)
      complex,            intent(in) :: U_full(:,:,:,:)      ! (num_bands,num_wann,nkptf,jspins)
      type(t_wann_ft),    intent(in) :: gmatRealspace(:,:)     ! (3*nat, jspins) elph real-space tensors
      type(t_kpts),       intent(in) :: qpts_dyn
      type(t_sym),        intent(in) :: sym_dyn

      integer :: nwann, nmodes, nat, ispin
      integer :: nNZ_k, nNZ_fcm, nNZ_g
      complex, allocatable :: hamWannier(:,:,:)   ! (nwann,nwann,nNZ_k)       electron H, real space [Ry]
      complex, allocatable :: fcm(:,:,:)          ! (nmodes,nmodes,nNZ_fcm)   force constants, real space
      integer, allocatable :: Rvecs_fcm(:,:)      ! (3,nNZ_fcm) phonon WS vectors (own q-mesh)
      real,    allocatable :: weightNZ_fcm(:)     ! (nNZ_fcm)   phonon WS weights (= 1/ndegen)
      integer, allocatable :: ord_k(:), ord_fcm(:), ord_g(:)   ! FLEUR index of EPW's i-th WS vector
      character(len=:), allocatable :: prefix

      nat     = fi%atoms%nat
      nmodes  = 3*nat
      nwann = fi%wannierlib%num_wann
      ispin   = 1
      prefix  = trim(adjustl(fi%dfpt%epw_prefix))

      if (fi%input%jspins /= 1) &
         call juDFT_warn("EPW restart export only writes spin channel 1 (jspins>1 not supported yet).", &
                         calledby="write_epw_restart_files")

      ! Electron (R_e) and e-ph (R_g) WS lists come from the elph tensor, whose R_p
      ! axis (ws_q) belongs to the coarse phonon q mesh. The force constants (fcm)
      ! get their own WS list from qpts_dyn; EPW keeps the two phonon lists independent.
      nNZ_k = gmatRealspace(1,ispin)%ws_k_nNZ
      nNZ_g = gmatRealspace(1,ispin)%ws_q_nNZ

      ! --- electron Hamiltonian real space (R_e) ---
      call build_ham_wannier(fi, results, U_full, ispin, gmatRealspace(1,ispin), nwann, hamWannier)

      ! --- force constants real space (R_q) on the coarse phonon q-mesh ---
      call build_fcm(fi, dynMats, nmodes, fcm, Rvecs_fcm, weightNZ_fcm, nNZ_fcm, qpts_dyn, sym_dyn)

      call epw_ws_order(gmatRealspace(1,ispin)%ws_k_Rvecs, ord_k)
      call epw_ws_order(Rvecs_fcm, ord_fcm)
      call epw_ws_order(gmatRealspace(1,ispin)%ws_q_Rvecs, ord_g)

      write (oUnit,'(/,a)') "Writing EPW (epwread) restart files:"
      write (oUnit,'(a,5(1x,i0))') "   nwann, nmodes, nNZ_k, nNZ_fcm, nNZ_g =", nwann, nmodes, nNZ_k, nNZ_fcm, nNZ_g
      write (oUnit,'(a,a)')        "   prefix = ", prefix
      write (oUnit,'(a,3(1x,i0),a,3(1x,i0))') "   EPW input must use nk1..3 =", fi%kpts%nkpt3, &
                                              ",  nq1..3 =", qpts_dyn%nkpt3

      ! --- write files ---
      call write_crystal_fmt(fi, results, nwann, nmodes)
      call write_epwdata_fmt(results, nwann, nmodes, nNZ_k, nNZ_fcm, nNZ_g, nat, hamWannier, fcm, ord_k, ord_fcm)
      call write_vme_dme_fmt(nwann, nNZ_k)
      call write_wigner_fmt(fi, gmatRealspace(1,ispin), Rvecs_fcm, weightNZ_fcm, nNZ_k, nNZ_fcm, nNZ_g, ord_k, ord_fcm, ord_g)
      call write_ukk_fmt(fi, prefix)
      call write_epmatwp(fi, gmatRealspace(:,ispin), nwann, nmodes, nNZ_k, nNZ_g, prefix, ord_k, ord_g)

      write (oUnit,'(a)') "   ... EPW restart files written."

   end subroutine write_epw_restart_files

   !---------------------------------------------------------------------------

   subroutine build_ham_wannier(fi, results, U_full, ispin, ft, nwann, hamWannier)
      !! Real-space Wannier-gauge electron Hamiltonian [Ry] on the WS vectors ft%ws_k_Rvecs.
      !! Mirrors wannier_matrix_interpolate's forward transform (m_matrix_interpolation).
      type(t_fleurinput), intent(in)  :: fi
      type(t_results),    intent(in)  :: results
      complex,            intent(in)  :: U_full(:,:,:,:)
      integer,            intent(in)  :: ispin, nwann
      type(t_wann_ft),    intent(in)  :: ft
      complex, allocatable, intent(out) :: hamWannier(:,:,:)

      integer :: num_bands, ikpt, ib
      complex, allocatable :: H_bloch(:,:), matRot(:,:)

      num_bands = fi%wannierlib%max_band - fi%wannierlib%min_band + 1

      allocate(H_bloch(num_bands,num_bands), matRot(nwann,nwann))
      allocate(hamWannier(nwann,nwann,ft%ws_k_nNZ))
      hamWannier = cmplx(0.0,0.0)

      do ikpt = 1, fi%kpts%nkpt
         H_bloch = cmplx(0.0,0.0)
         do ib = 1, num_bands
            H_bloch(ib,ib) = cmplx(results%eig(fi%wannierlib%min_band + ib - 1, fi%kpts%bkp(ikpt), ispin), 0.0)
         end do
         ! Wannier gauge: U^dagger H U
         matRot = matmul(conjg(transpose(U_full(:,:,ikpt,ispin))), matmul(H_bloch, U_full(:,:,ikpt,ispin)))
         call ft_ws_forward(ft%ws_k_nNZ, ft%ws_k_Rvecs, fi%kpts%bk(:,ikpt), matRot, hamWannier)
      end do
      hamWannier = htr2ry * hamWannier / fi%kpts%nkpt
   end subroutine build_ham_wannier

   subroutine build_fcm(fi, dynMats, nmodes, fcm, Rvecs_fcm, weightNZ_fcm, nNZ_fcm, qpts_dyn, sym_dyn)
      !! Real-space force constants (fcm) on the coarse PHONON q-mesh, together
      !! with the phonon Wigner-Seitz list (Rvecs_fcm, ndegen = 1/weightNZ_fcm).
      !! This mesh is independent of the electron/e-ph (k-mesh) WS lists.
      !!
      !! dynMats lives on qpts_dyn%bk; both the Fourier box (nkpt3) and the
      !! unfolding inside ft_dyn refer to the full zone qpts_dyn%bkf.
      !!
      !! Scaling into EPW's convention (see module header):
      !!   fcm block(a,b) = 4 * amu_ry * sqrt(m_a m_b) * D_norm(R).
      !! With l_bornhuang, fcm is projected onto the invariance conditions.
      type(t_fleurinput),   intent(in)  :: fi
      complex,              intent(in)  :: dynMats(:,:,:)
      integer,              intent(in)  :: nmodes
      complex, allocatable, intent(out) :: fcm(:,:,:)          ! (nmodes,nmodes,nNZ_fcm)
      integer, allocatable, intent(out) :: Rvecs_fcm(:,:)      ! (3,nNZ_fcm) WS vectors (un-negated)
      real,    allocatable, intent(out) :: weightNZ_fcm(:)     ! (nNZ_fcm)   = 1/ndegen
      integer,              intent(out) :: nNZ_fcm
      type(t_kpts),         intent(in)  :: qpts_dyn
      type(t_sym),          intent(in)  :: sym_dyn

      integer :: nq(3), bigBox(2,3), boxSize, ix, iy, iz, iGrid, ia, ja, i0, j0
      real    :: fac
      integer, allocatable :: supercellR(:,:)
      real,    allocatable :: FTweight(:)
      complex, allocatable :: dyn_mat_q_full(:,:,:)
      type(t_cell) :: cellLocal

      if (size(dynMats, 3) < qpts_dyn%nkpt) then
         call juDFT_error("Fewer dynamical matrices than irreducible q points.", calledby="build_fcm")
      end if

      nq = qpts_dyn%nkpt3(:)

      ! Wigner-Seitz weights of the coarse phonon q-mesh supercell
      ! (mirrors interpolate_dynmat / wannier_matrixq_forward: big box +-2N, always WS).
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
      call build_ws_ft(bigBox, FTweight, nNZ_fcm, Rvecs_fcm, weightNZ_fcm)

      ! D_norm(R) = (1/Nq) sum_q e^{+i q.R} D(q)   (symmetry-unfolded to the full BZ)
      allocate(fcm(nmodes,nmodes,nNZ_fcm))
      call ft_dyn(fi%atoms, qpts_dyn, sym_dyn, nNZ_fcm, Rvecs_fcm, fi%cell%amat, dynMats, fcm, dyn_mat_q_full)

      ! per-(atom,atom) block mass/energy scaling (see module header)
      do ia = 1, fi%atoms%nat
         i0 = 3*(ia-1)
         do ja = 1, fi%atoms%nat
            j0 = 3*(ja-1)
            fac = htr2ry**2 * amu_ry * sqrt(atomicMasses_const(fi%atoms%nz(fi%atoms%itype(ia))) * atomicMasses_const(fi%atoms%nz(fi%atoms%itype(ja))))
            fcm(i0+1:i0+3, j0+1:j0+3, :) = fac * fcm(i0+1:i0+3, j0+1:j0+3, :)
         end do
      end do

      ! fcm is now proportional to the de-normalized force constants
      if (fi%dfpt%l_bornhuang) then
         call timestart("Born-Huang projection")
         call dfpt_born_huang(fi%atoms, fi%cell, nq, nNZ_fcm, Rvecs_fcm, weightNZ_fcm, fcm)
         call timestop("Born-Huang projection")
      end if
   end subroutine build_fcm

   !---------------------------------------------------------------------------

   subroutine epw_ws_order(Rvecs, ord)
      !! ord(i) = index into Rvecs of EPW's i-th WS vector, where EPW's vector is -Rvecs(:,ord(i)).
      !! Same key as EPW's wigner_seitzkq/wigner_seitzg sort.
      integer,              intent(in)  :: Rvecs(:,:)
      integer, allocatable, intent(out) :: ord(:)

      real, allocatable :: key(:)

      allocate(ord(size(Rvecs,2)), key(size(Rvecs,2)))
      key = -(Rvecs(1,:) * 1.0e6 + Rvecs(2,:) * 1.0e3 + Rvecs(3,:) * 1.0)
      call sort(ord, key)
   end subroutine epw_ws_order

   subroutine write_crystal_fmt(fi, results, nwann, nmodes)
      !! crystal.fmt -- exact record order of EPW's epw_write (io.f90:219-235).
      type(t_fleurinput), intent(in) :: fi
      type(t_results),    intent(in) :: results
      integer,            intent(in) :: nwann, nmodes

      integer, parameter :: ntypx = 10           ! QE Modules/parameters.f90
      integer :: iu, ia, it
      real    :: alat, at(3,3), bg(3,3), amass(ntypx)
      real, allocatable :: tau(:,:), w_centers(:,:)

      alat = norm2(fi%cell%amat(:,1))             ! any consistent scale works (w_centers=0)
      at   = fi%cell%amat / alat                  ! lattice vectors in columns, alat units
      do it = 1, 3
         bg(:,it) = fi%cell%bmat(it,:) * alat / tpi_const   ! recip vectors in columns, 2pi/alat units
      end do
      allocate(tau(3,fi%atoms%nat), w_centers(3,nwann))
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
      write(iu,*) w_centers                                      ! (3,nwann) zeros
      write(iu,*) 0.0                                            ! L (2D cutoff length)
      write(iu,*) htr2ry * fi%input%tkb                          ! degauss [Ry]
      write(iu,*) 0                                              ! ngauss
      write(iu,*) .false.                                        ! lda_plus_u
      close(iu)
   end subroutine write_crystal_fmt

   subroutine write_epwdata_fmt(results, nwann, nmodes, nNZ_k, nNZ_fcm, nNZ_g, nat, hamWannier, fcm, ord_k, ord_fcm)
      !! epwdata.fmt -- EPW's epw_write (io.f90:237-273), R index in EPW's order.
      type(t_results), intent(in) :: results
      integer,         intent(in) :: nwann, nmodes, nNZ_k, nNZ_fcm, nNZ_g, nat
      complex,         intent(in) :: hamWannier(:,:,:), fcm(:,:,:)
      integer,         intent(in) :: ord_k(:), ord_fcm(:)

      integer :: iu, ibnd, jbnd, iRk, imode, jmode, iRfcm
      real    :: zstar(3,3,nat), epsi(3,3)   ! EPW stores these as REAL (18 reals in epwdata.fmt)

      zstar = 0.0     ! Born effective charges (0 for non-polar metal)
      epsi  = 0.0     ! dielectric tensor

      open(newunit=iu, file="epwdata.fmt", status='replace', action='write', form='formatted')
      write(iu,*) htr2ry * results%ef                            ! Fermi energy [Ry]
      write(iu,*) nwann, nNZ_k, nmodes, nNZ_fcm, nNZ_g
      write(iu,*) zstar, epsi
      do ibnd = 1, nwann
         do jbnd = 1, nwann
            do iRk = 1, nNZ_k
               write(iu,*) hamWannier(ibnd, jbnd, ord_k(iRk))
            end do
         end do
      end do
      do imode = 1, nmodes
         do jmode = 1, nmodes
            do iRfcm = 1, nNZ_fcm
               write(iu,*) fcm(imode, jmode, ord_fcm(iRfcm))
            end do
         end do
      end do
      close(iu)
   end subroutine write_epwdata_fmt

   subroutine write_vme_dme_fmt(nwann, nNZ_k)
      !! vmedata.fmt / dmedata.fmt -- velocity/position/dipole matrix elements.
      !! Written as zeros (unused by Eliashberg); exact interleaving of io.f90:246-251.
      integer, intent(in) :: nwann, nNZ_k

      integer :: iuv, iud, ibnd, jbnd, iRk, ipol
      complex, parameter :: c0 = (0.0,0.0)

      open(newunit=iuv, file="vmedata.fmt", status='replace', action='write', form='formatted')
      open(newunit=iud, file="dmedata.fmt", status='replace', action='write', form='formatted')
      do ibnd = 1, nwann
         do jbnd = 1, nwann
            do iRk = 1, nNZ_k
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

   subroutine write_wigner_fmt(fi, ft, Rvecs_fcm, weightNZ_fcm, nNZ_k, nNZ_fcm, nNZ_g, ord_k, ord_fcm, ord_g)
      !! wigner.fmt -- EPW's epw_write_ws_data (io.f90:669-737), dims=dims2=1 case.
      !! Lists EPW's WS vectors (negated FLEUR vectors, EPW order) as the data is stored.
      type(t_fleurinput), intent(in) :: fi
      type(t_wann_ft),    intent(in) :: ft
      integer,            intent(in) :: Rvecs_fcm(:,:)     ! (3,nNZ_fcm) phonon WS vectors (un-negated)
      real,               intent(in) :: weightNZ_fcm(:)    ! (nNZ_fcm)   = 1/ndegen
      integer,            intent(in) :: nNZ_k, nNZ_fcm, nNZ_g
      integer,            intent(in) :: ord_k(:), ord_fcm(:), ord_g(:)

      integer :: iu, ir, j, irvec(3)
      real    :: alat

      alat = norm2(fi%cell%amat(:,1))

      open(newunit=iu, file="wigner.fmt", status='replace', action='write', form='formatted')
1     format(1000(I0, :, 1X))
2     format(3I6, ES26.17E3)
      write(iu,1) nNZ_k, nNZ_fcm, nNZ_g, 1, 1
      ! electron block (R_e) -- elph-tensor k-mesh WS list
      do ir = 1, nNZ_k
         j = ord_k(ir)
         irvec = -ft%ws_k_Rvecs(:,j)
         write(iu,2) irvec, ws_len(fi%cell%amat, ft%ws_k_Rvecs(:,j), alat)
         write(iu,1) nint(1.0/ft%ws_k_weightNZ(j))
      end do
      ! phonon block (R_q) -- force-constant q-mesh WS list
      do ir = 1, nNZ_fcm
         j = ord_fcm(ir)
         irvec = -Rvecs_fcm(:,j)
         write(iu,2) irvec, ws_len(fi%cell%amat, Rvecs_fcm(:,j), alat)
         write(iu,1) nint(1.0/weightNZ_fcm(j))
      end do
      ! e-ph block (R_g) -- elph-tensor q-mesh WS list (R_p axis of the tensor)
      do ir = 1, nNZ_g
         j = ord_g(ir)
         irvec = -ft%ws_q_Rvecs(:,j)
         write(iu,2) irvec, ws_len(fi%cell%amat, ft%ws_q_Rvecs(:,j), alat)
         write(iu,1) nint(1.0/ft%ws_q_weightNZ(j))
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

   subroutine write_epmatwp(fi, gmatSpin, nwann, nmodes, nNZ_k, nNZ_g, prefix, ord_k, ord_g)
      !! <prefix>.epmatwp -- direct-access unformatted, one record per iRg,
      !! record = epmatwp(nwann,nwann,nNZ_k,nmodes) (EPW serial branch, io.f90:338-351).
      !! Both R indices in EPW's order.
      type(t_fleurinput), intent(in) :: fi
      type(t_wann_ft),    intent(in) :: gmatSpin(:)             ! gmatRealspace(:,ispin), size 3*nat
      integer,            intent(in) :: nwann, nmodes, nNZ_k, nNZ_g
      character(len=*),   intent(in) :: prefix
      integer,            intent(in) :: ord_k(:), ord_g(:)

      integer :: iu, iRg, imode, iRk, lrepmatw, direct_io_factor, unf_recl
      real    :: dummy
      complex, allocatable :: gmatRecord(:,:,:,:)

      lrepmatw = 2 * nwann * nwann * nNZ_k * nmodes
      inquire(iolength=direct_io_factor) dummy
      unf_recl = direct_io_factor * lrepmatw

      allocate(gmatRecord(nwann,nwann,nNZ_k,nmodes))
      open(newunit=iu, file=trim(prefix)//".epmatwp", status='replace', action='write', &
           form='unformatted', access='direct', recl=unf_recl)
      do iRg = 1, nNZ_g
         do imode = 1, nmodes
            gmatRecord(:,:,:,imode) = htr2ry * gmatSpin(imode)%matWannier(:,:,ord_k-1,ord_g(iRg)-1)
         end do
         write(iu, rec=iRg) gmatRecord
      end do
      close(iu)
   end subroutine write_epmatwp

end module m_dfpt_write_epw
