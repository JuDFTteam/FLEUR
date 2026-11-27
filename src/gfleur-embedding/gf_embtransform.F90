!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_embtransform
   !Embedding potential of a substrate on a plane S for an ad-layer whose
   !2D lattice differs from the substrate lattice (incommensurate or a
   !supercell). Controlled by the file incomm_inp (&INCOMM).
   !
   !In the ad-layer plane-wave basis kappa+G the substrate couples kappa+G_j
   !to kappa+G_i if (G_i-G_j) is a substrate reciprocal vector u:
   !   Sigma0(i,j) = exp(-i u.d) sigma_u((kappa+G_j)M),
   !M maps ad-layer to substrate fractional coordinates, d is the lateral
   !substrate shift. The remaining Umklapp couplings to kappa+u+G are folded
   !back with the embedding potential Sigma_up of the ad-layer region at the
   !Bloch vectors kappa+u (Loewdin partitioning, G_S^-1 = Sigma_L + Sigma_R):
   !   Sigma_eff = Sigma0 - Delta_0Q [D_Q + Delta_QQ]^-1 Delta_Q0,
   !   D_u = Sigma0(kappa+u) + Sigma_up(kappa+u).
#ifdef CPP_MPI
   USE mpi
#endif
   USE hdf5
   USE m_hdf_tools
   USE m_juDFT
   USE m_constants, ONLY: oUnit, pi_const
   USE m_gf_types
   USE m_gf_energies, ONLY: gf_noen, gf_Z
   USE m_gf_io2dmat, ONLY: gf_io2dmat_create, t_io2dmeta
   USE m_gf_sigtab
   USE m_gf_out, ONLY: gf_out_newheader
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: gf_embtransform

   INTEGER, PARAMETER :: ORDER_FULL = 99

   TYPE t_incomm
      CHARACTER(LEN=100) :: subfile = "", outfile = "embpot11", upfile = ""
      REAL    :: theta = 0.0, asub(2, 2) = 0.0, shift(2) = 0.0
      REAL    :: eps_comm = 1.0E-6, vvac = 0.0, eps_energy = 1.0E-5
      INTEGER :: nshell = 1, order = 0, interp = INTERP_CUBIC
      LOGICAL :: write_kaux = .FALSE., write_ksub = .FALSE., stop_after = .FALSE.
      LOGICAL :: l_asub = .FALSE., l_vacextrap = .FALSE.
   END TYPE t_incomm

   !geometry: ad-layer (o) and rotated substrate (b) lattices
   TYPE t_igeo
      REAL :: ao(2, 2), bo(2, 2), bbo(2, 2)
      REAL :: asub(2, 2), ab(2, 2), bb(2, 2), rot(2, 2)
      !fractional ad-layer -> fractional substrate: f_b = f_o m
      REAL :: m(2, 2), minv(2, 2)
      REAL :: rk
   END TYPE t_igeo

CONTAINS

   SUBROUTINE gf_embtransform(gmpi, ld, gkpts, l_vref, vref)
      !ld: the ad-layer region (layer 1), its lower plane is S; vref: energy
      !reference of this run (potential at the end of its vacuum), if known
      TYPE(t_gfmpi), INTENT(IN)   :: gmpi
      TYPE(t_gflayer), INTENT(IN) :: ld
      TYPE(t_kpts), INTENT(IN)    :: gkpts
      LOGICAL, INTENT(IN)         :: l_vref
      REAL, INTENT(IN)            :: vref

      TYPE(t_incomm) :: inp
      INTEGER :: ierr

      CALL priv_read_input(inp)
      CALL gf_out_newheader("Substrate embedding potential for the ad-layer cell (incomm_inp)")
      IF (gmpi%pe0) CALL priv_run(inp, ld, gkpts, l_vref, vref)
#ifdef CPP_MPI
      CALL MPI_BARRIER(gmpi%fmpi%mpi_comm, ierr)
#endif
      ierr = 0
      IF (inp%write_kaux .OR. inp%write_ksub) &
         CALL juDFT_end("incomm: k-point list written", gmpi%fmpi%irank)
      IF (inp%stop_after) CALL juDFT_end("incomm: embedding potential written", gmpi%fmpi%irank)
   END SUBROUTINE gf_embtransform

   SUBROUTINE priv_read_input(inp)
      TYPE(t_incomm), INTENT(OUT) :: inp

      CHARACTER(LEN=100) :: subfile, outfile, upfile
      CHARACTER(LEN=10)  :: order, interp, extrap
      REAL    :: theta, asub(2, 2), shift(2), eps_comm, vvac, eps_energy
      INTEGER :: nshell, ios
      LOGICAL :: write_kaux, write_ksub, stop_after
      NAMELIST /INCOMM/ subfile, outfile, upfile, theta, asub, shift, eps_comm, vvac, eps_energy, &
         nshell, order, interp, extrap, write_kaux, write_ksub, stop_after

      subfile = inp%subfile; outfile = inp%outfile; upfile = inp%upfile
      theta = inp%theta; asub = 0.0; shift = inp%shift
      eps_comm = inp%eps_comm; vvac = inp%vvac; eps_energy = inp%eps_energy
      nshell = inp%nshell; order = "0"; interp = "cubic"; extrap = "none"
      write_kaux = .FALSE.; write_ksub = .FALSE.; stop_after = .FALSE.
      OPEN (99, FILE="incomm_inp", STATUS="old", ACTION="read")
      READ (99, NML=INCOMM, IOSTAT=ios)
      CLOSE (99)
      IF (ios /= 0) CALL juDFT_error("Error reading &INCOMM in incomm_inp", calledby="gf_embtransform")
      IF (LEN_TRIM(subfile) == 0) CALL juDFT_error("&INCOMM: subfile missing", calledby="gf_embtransform")
      inp%subfile = subfile; inp%outfile = outfile; inp%upfile = upfile
      IF (TRIM(subfile) == TRIM(outfile)) CALL juDFT_error("&INCOMM: subfile and outfile must differ", &
                                                           calledby="gf_embtransform")
      inp%theta = theta; inp%shift = shift
      inp%l_asub = ANY(asub /= 0.0)
      inp%asub = asub
      inp%eps_comm = eps_comm; inp%vvac = vvac; inp%eps_energy = eps_energy
      inp%nshell = nshell
      inp%write_kaux = write_kaux; inp%write_ksub = write_ksub; inp%stop_after = stop_after
      SELECT CASE (TRIM(ADJUSTL(order)))
      CASE ("0")
         inp%order = 0
      CASE ("2")
         inp%order = 2
      CASE ("3")
         inp%order = 3
      CASE ("full")
         inp%order = ORDER_FULL
      CASE DEFAULT
         CALL juDFT_error("&INCOMM: order must be '0', '2', '3' or 'full'", calledby="gf_embtransform")
      END SELECT
      SELECT CASE (TRIM(ADJUSTL(interp)))
      CASE ("cubic")
         inp%interp = INTERP_CUBIC
      CASE ("linear")
         inp%interp = INTERP_LINEAR
      CASE ("exact")
         inp%interp = INTERP_EXACT
      CASE DEFAULT
         CALL juDFT_error("&INCOMM: interp must be 'cubic', 'linear' or 'exact'", calledby="gf_embtransform")
      END SELECT
      SELECT CASE (TRIM(ADJUSTL(extrap)))
      CASE ("none")
         inp%l_vacextrap = .FALSE.
      CASE ("vac")
         inp%l_vacextrap = .TRUE.
      CASE DEFAULT
         CALL juDFT_error("&INCOMM: extrap must be 'none' or 'vac'", calledby="gf_embtransform")
      END SELECT
      IF (inp%order > 0 .AND. inp%nshell < 1) CALL juDFT_error("&INCOMM: order>0 needs nshell>=1", &
                                                               calledby="gf_embtransform")
      IF (inp%order > 0 .AND. .NOT. inp%write_kaux .AND. LEN_TRIM(upfile) == 0) &
         CALL juDFT_error("&INCOMM: order>0 needs upfile", hint="generate it with write_kaux=t and an auxiliary run", &
                          calledby="gf_embtransform")
   END SUBROUTINE priv_read_input

   SUBROUTINE priv_run(inp, ld, gkpts, l_vref, vref)
      TYPE(t_incomm), INTENT(IN)  :: inp
      TYPE(t_gflayer), INTENT(IN) :: ld
      TYPE(t_kpts), INTENT(IN)    :: gkpts
      LOGICAL, INTENT(IN)         :: l_vref
      REAL, INTENT(IN)            :: vref

      TYPE(t_embfile)  :: f, up
      TYPE(t_sigtab)   :: tab
      TYPE(t_igeo)     :: geo
      TYPE(t_io2dmeta) :: fmeta
      INTEGER, ALLOCATABLE :: lc(:, :), q(:, :), ulist(:, :), dmap(:, :), hl(:), ie_sub(:), ie_up(:)
      COMPLEX, ALLOCATABLE :: s0(:, :, :), seff(:, :), sfull(:, :)
      REAL, ALLOCATABLE    :: a(:, :, :)
      INTEGER(HID_T) :: fid, varid
      INTEGER :: ng, nlc, nq, dmax, n, nk, ichan, nchan, jsp_out, ns, nh, ia, ib, iu, hdferr, ic, ispin
      INTEGER :: nmiss
      REAL    :: th, kap(2), p(2), uc(2), dnorm(4), interr
      COMPLEX :: val(4), ph
      LOGICAL :: ok, l_noco
      INTEGER, PARAMETER :: outspin(4) = (/1, 4, 3, 2/)

      l_noco = ld%fi%noco%l_noco
      ng = SIZE(ld%lapw_gf%global2Dlist, 2)

      CALL f%open(inp%subfile)
      !<-- geometry
      geo%ao = ld%fi%cell%amat(1:2, 1:2)
      geo%bo = ld%fi%cell%bmat(1:2, 1:2)
      geo%bbo = ld%fi%cell%bbmat(1:2, 1:2)
      geo%rk = ld%lapw_gf%rkmax
      IF (inp%l_asub) THEN
         geo%asub = inp%asub
      ELSE IF (f%l_amat) THEN
         geo%asub = f%amat2d
      ELSE
         CALL juDFT_error("The substrate file has no amat2d", hint="set asub in &INCOMM", calledby="gf_embtransform")
      END IF
      th = inp%theta*pi_const/180.0
      geo%rot = RESHAPE((/COS(th), SIN(th), -SIN(th), COS(th)/), (/2, 2/))
      geo%ab = MATMUL(geo%rot, geo%asub)
      geo%bb = 2.0*pi_const*priv_inv2(geo%ab)
      geo%m = MATMUL(geo%bo, priv_inv2(geo%bb))
      geo%minv = priv_inv2(geo%m)
      WRITE (oUnit, "(a,4f12.6)") "ad-layer amat (columns):           ", geo%ao
      WRITE (oUnit, "(a,4f12.6)") "substrate amat (columns, rotated): ", geo%ab
      WRITE (oUnit, "(a,f10.4,a,2f10.4)") "rotation theta (deg):", inp%theta, "  lateral shift (a.u.):", inp%shift
      WRITE (oUnit, "(a,4f12.6)") "M (f_sub = f_ad M):                ", geo%m
      IF (f%l_zplane) WRITE (oUnit, "(a,f10.4)") "substrate: distance plane - top atom plane (a.u.):", f%zplane
      !energy zero: the substrate and the ad-layer must use the same one
      IF (f%l_vref .AND. l_vref) THEN
         WRITE (oUnit, "(a,f10.5,a,f10.5)") "energy reference (vacuum potential, Htr): substrate", f%vref, &
            "  ad-layer", vref
         IF (ABS(f%vref - vref) > 1.0E-3) CALL priv_note("NOTE: the energy zeros of the substrate and the ad-layer differ, "// &
                                                         "check vshift of gfprep")
      ELSE
         WRITE (oUnit, "(a)") "energy reference: not available in both runs (no vacuum embedding)"
      END IF
      IF (f%l_rkmax) THEN
         WRITE (oUnit, "(a,f10.4,a,f10.4)") "Kmax substrate:", f%rkmax, "  ad-layer:", geo%rk
      ELSE
         WRITE (oUnit, "(a)") "substrate file without rkmax: basis taken from the nonzero diagonal"
      END IF

      !<-- common reciprocal vectors Lambda_c: G_i-G_j that are substrate vectors
      dmax = 2*MAXVAL(ABS(ld%lapw_gf%global2Dlist))
      ALLOCATE (dmap(-dmax:dmax, -dmax:dmax), lc(2, (2*dmax + 1)**2))
      CALL priv_common(ld%lapw_gf%global2Dlist, geo, inp%eps_comm, dmax, dmap, lc, nlc)
      WRITE (oUnit, "(a,i0)") "common reciprocal vectors (Lambda_c without 0) in the basis: ", nlc - 1
      IF (nlc == 1) WRITE (oUnit, "(a)") "incommensurate: Sigma0 is diagonal in the ad-layer basis"

      !<-- Umklapp vectors of the correction
      nq = 0
      ALLOCATE (q(2, 0))
      IF (inp%order > 0 .OR. inp%write_kaux) THEN
         IF (nlc > 1) CALL juDFT_error("The Umklapp correction needs an incommensurate ad-layer", &
                                       hint="use order='0' for commensurate cells (Sigma0 is exact there)", &
                                       calledby="gf_embtransform")
         CALL priv_shells(geo, inp%nshell, q)
         nq = SIZE(q, 2)
         WRITE (oUnit, "(a,i0,a,i0,a)") "Umklapp correction: ", nq, " vectors in ", inp%nshell, " shell(s)"
         DO n = 1, nq
            WRITE (oUnit, "(a,2i4,a,f10.5)") "   u=", q(:, n), "  |u|=", NORM2(MATMUL(REAL(q(:, n)), geo%bb))
         END DO
      END IF

      IF (inp%write_kaux) THEN
         CALL priv_write_kaux(gkpts, geo, q)
         CALL f%close()
         RETURN
      END IF
      IF (inp%write_ksub) THEN
         CALL priv_write_ksub(gkpts, geo, q, ld%lapw_gf%global2Dlist)
         CALL f%close()
         RETURN
      END IF

      !<-- symmetry of the ad-layer k-set must be a symmetry of the substrate
      CALL priv_symcheck(ld%fi%sym, geo, f, inp%shift, nlc > 1)

      !<-- table of the substrate amplitudes
      CALL priv_ulist(lc(:, :nlc), q, ulist)
      !a noco ad-layer needs 2x2 spin blocks, also from a collinear substrate
      CALL tab%init(f, ulist, inp%interp, geo%asub, .TRUE., l_spinor=l_noco)
      tab%l_vacextrap = inp%l_vacextrap
      tab%vvac = inp%vvac
      ns = tab%ns

      !<-- energies
      ALLOCATE (ie_sub(gf_noen()), ie_up(gf_noen()))
      ie_up = 0
      DO n = 1, gf_noen()
         ie_sub(n) = f%find_energy(gf_Z(n, 0), inp%eps_energy)
         IF (ie_sub(n) == 0) CALL juDFT_error("Energy of gf_en not in the substrate file", &
                                              hint="use the same energy grid for the substrate run", &
                                              calledby="gf_embtransform")
      END DO
      IF (inp%order > 0) THEN
         CALL up%open(inp%upfile)
         DO n = 1, gf_noen()
            ie_up(n) = up%find_energy(gf_Z(n, 0), inp%eps_energy)
            IF (ie_up(n) == 0) CALL juDFT_error("Energy of gf_en not in upfile", calledby="gf_embtransform")
         END DO
      END IF

      !<-- spin channels
      jsp_out = ld%fi%input%jspins
      IF (l_noco) jsp_out = 4
      IF (ns == 2) THEN
         IF (.NOT. l_noco) CALL juDFT_error("A noco substrate file needs a noco ad-layer", calledby="gf_embtransform")
         nchan = 1
      ELSE
         nchan = jsp_out
      END IF

      !<-- output file
      fmeta%l_set = .TRUE.
      fmeta%amat2d = geo%ao
      fmeta%rkmax = geo%rk
      fmeta%l_noco = l_noco
      fmeta%l_soc = ld%fi%noco%l_soc
      fmeta%l_timerev = .NOT. fmeta%l_soc
      fmeta%zplane = -1.0
      IF (f%l_zplane) fmeta%zplane(1) = f%zplane
      fmeta%l_vref = f%l_vref
      fmeta%vref = f%vref
      CALL gf_io2dmat_create(TRIM(inp%outfile), ld%lapw_gf, gkpts, ld%fi%sym, jsp_out, gkpts%nkpt, ng, &
                             fmeta, fid, varid, l_serial=.TRUE.)
      ALLOCATE (a(ng, ng, 2), hl(ng))

      DO n = 1, gf_noen()
         DO ichan = 1, nchan
            ispin = ichan
            CALL tab%fill(f, ie_sub(n), ispin)
            tab%energy = gf_Z(n, 0)
            CALL tab%reset_stat()
            dnorm = 0.0
            nmiss = 0
            DO nk = 1, gkpts%nkpt
               kap = gkpts%bk(1:2, nk)
               CALL priv_home(kap, ld%lapw_gf%global2Dlist, ld%fi%cell%bbmat, geo%rk, hl, nh)
               ALLOCATE (s0(nh, nh, ns**2))
               s0 = CMPLX(0.0, 0.0)
               DO ib = 1, nh
                  p = MATMUL(kap + ld%lapw_gf%global2Dlist(:, hl(ib)), geo%m)
                  DO ia = 1, nh
                     iu = dmap(ld%lapw_gf%global2Dlist(1, hl(ia)) - ld%lapw_gf%global2Dlist(1, hl(ib)), &
                               ld%lapw_gf%global2Dlist(2, hl(ia)) - ld%lapw_gf%global2Dlist(2, hl(ib)))
                     IF (iu == 0) CYCLE
                     CALL tab%eval(f, p, tab%uindex(lc(:, iu)), val(:ns**2), ok)
                     IF (.NOT. ok) CALL priv_fail(p, lc(:, iu), tab)
                     uc = MATMUL(REAL(lc(:, iu)), geo%bb)
                     ph = EXP(CMPLX(0.0, -DOT_PRODUCT(uc, inp%shift)))
                     s0(ia, ib, :) = ph*val(:ns**2)
                     IF (ns == 2) THEN
                        !spin frame rotated with the lattice: U sigma U^+, U=exp(-i theta sigma_z/2)
                        s0(ia, ib, 2) = s0(ia, ib, 2)*EXP(CMPLX(0.0, th))
                        s0(ia, ib, 3) = s0(ia, ib, 3)*EXP(CMPLX(0.0, -th))
                     END IF
                  END DO
               END DO
               IF (inp%order > 0) THEN
                  !full spinor matrix (index (s-1)*nh+G), components 11,21,12,22
                  ALLOCATE (seff(ns*nh, ns*nh), sfull(ns*nh, ns*nh))
                  DO ic = 1, ns**2
                     sfull(MOD(ic - 1, ns)*nh + 1:(MOD(ic - 1, ns) + 1)*nh, &
                           ((ic - 1)/ns)*nh + 1:((ic - 1)/ns + 1)*nh) = s0(:, :, ic)
                  END DO
                  CALL priv_umklapp(inp, tab, f, up, ie_up(n), ispin, geo, kap, ld%lapw_gf%global2Dlist, hl(:nh), &
                                    ld%fi%cell%bbmat, q, sfull, seff, dnorm, nmiss)
                  DO ic = 1, ns**2
                     s0(:, :, ic) = seff(MOD(ic - 1, ns)*nh + 1:(MOD(ic - 1, ns) + 1)*nh, &
                                         ((ic - 1)/ns)*nh + 1:((ic - 1)/ns + 1)*nh)
                  END DO
                  DEALLOCATE (seff, sfull)
               END IF
               !<-- write: noco spins 1..4 = 11,22,12,21
               DO ic = 1, ns**2
                  a = 0.0
                  a(hl(:nh), hl(:nh), 1) = REAL(s0(:, :, ic))
                  a(hl(:nh), hl(:nh), 2) = AIMAG(s0(:, :, ic))
                  IF (ns == 2) THEN
                     ispin = outspin(ic)
                  ELSE
                     ispin = ichan
                  END IF
                  CALL priv_write(varid, a, ispin, n, nk)
               END DO
               DEALLOCATE (s0)
            END DO
            !<-- diagnostics of this energy and spin channel
            WRITE (oUnit, "(a,i4,a,i1,a,2f10.5,a,i0,a,i0,a,i0,a,i0)") "incomm E(", n, ") spin ", ichan, &
               " z=", gf_Z(n, 0), ": exact ", tab%n_exact, ", cubic ", tab%n_cubic, ", linear ", tab%n_linear, &
               ", vacuum extrapolation ", tab%n_extrap
            IF (tab%n_cubic > 0) THEN
               interr = tab%interr/MAX(tab%valmax, 1.0E-30)
               WRITE (oUnit, "(a,es10.2)") "   interpolation estimate max|cubic-linear|/max|sigma|: ", interr
               IF (interr > 1.0E-3) CALL priv_note("NOTE: interpolation estimate above 1e-3, "// &
                                                   "a denser substrate k-mesh or interp='exact' is recommended")
            END IF
            IF (tab%n_extrap > 0) CALL priv_note("NOTE: substrate amplitudes beyond the cutoff "// &
                                                 "replaced by the vacuum form (extrap='vac')")
            IF (inp%order > 0) THEN
               WRITE (oUnit, "(a,es10.2,a,es10.2)") "   max ||Delta_u0||/||Sigma0||:", dnorm(1), &
                  "   max ||Sigma~2||:", dnorm(2)
               IF (inp%order == ORDER_FULL .OR. inp%order == 3) &
                  WRITE (oUnit, "(a,es10.2)") "   max ||Sigma~ - Sigma~2||:", dnorm(3)
               WRITE (oUnit, "(a,i0,a,i0)") "   near-singular D_u (rcond<1e-8): ", NINT(dnorm(4)), &
                  ",  couplings dropped at the cutoff: ", nmiss
               IF (dnorm(4) > 0.0 .AND. inp%order == 2) &
                  CALL priv_note("NOTE: near-singular D_u (resonant Umklapp coupling), order='full' is recommended")
            END IF
         END DO
      END DO
      CALL priv_copy_boundary(f, fid, ld%stars%ng2, geo%rot, inp%shift)
      CALL io_dclose(varid, hdferr)
      CALL io_hdfclose(fid, hdferr)
      CALL f%close()
      IF (inp%order > 0) CALL up%close()
      WRITE (oUnit, "(a,a,a)") "incomm: embedding potential written to ", TRIM(inp%outfile), ".hdf"
      WRITE (*, "(a,a,a)") "incomm: embedding potential written to ", TRIM(inp%outfile), ".hdf"
   END SUBROUTINE priv_run

   SUBROUTINE priv_note(text)
      CHARACTER(LEN=*), INTENT(IN) :: text
      WRITE (oUnit, "(3x,a)") text
      WRITE (*, "(a)") text
   END SUBROUTINE priv_note

   SUBROUTINE priv_write(varid, a, ispin, n, nk)
      INTEGER(HID_T), INTENT(IN) :: varid
      REAL, INTENT(IN)           :: a(:, :, :)
      INTEGER, INTENT(IN)        :: ispin, n, nk
      CALL io_write(varid, (/1, 1, 1, ispin, n, nk/), (/SIZE(a, 1), SIZE(a, 2), 1, 1, 1, 1/), "A", a(:, :, 1))
      CALL io_write(varid, (/1, 1, 2, ispin, n, nk/), (/SIZE(a, 1), SIZE(a, 2), 1, 1, 1, 1/), "A", a(:, :, 2))
   END SUBROUTINE priv_write

   SUBROUTINE priv_fail(p, u, tab)
      REAL, INTENT(IN)           :: p(2)
      INTEGER, INTENT(IN)        :: u(2)
      TYPE(t_sigtab), INTENT(IN) :: tab
      WRITE (oUnit, "(a,2f12.6,a,2i4,a,f10.4)") "missing substrate amplitude: p=", p, " u=", u, &
         " |p|=", NORM2(MATMUL(p, tab%bsub))
      CALL juDFT_error("Substrate amplitude not available (outside its cutoff or k-set)", &
                       hint="increase Kmax of the substrate run, use interp='exact' only with the k-list of "// &
                       "write_ksub, or set extrap='vac'", calledby="gf_embtransform")
   END SUBROUTINE priv_fail

   SUBROUTINE priv_home(kap, glist, bbmat, rk, hl, nh)
      !G of the global list with |kappa+G|<=rkmax (same test as gf_apws)
      REAL, INTENT(IN)     :: kap(2), bbmat(3, 3), rk
      INTEGER, INTENT(IN)  :: glist(:, :)
      INTEGER, INTENT(OUT) :: hl(:), nh
      INTEGER :: n
      REAL    :: s(3)
      nh = 0
      s(3) = 0.0
      DO n = 1, SIZE(glist, 2)
         s(1:2) = kap + glist(:, n)
         IF (DOT_PRODUCT(s, MATMUL(s, bbmat)) <= rk*rk) THEN
            nh = nh + 1
            hl(nh) = n
         END IF
      END DO
   END SUBROUTINE priv_home

   SUBROUTINE priv_common(glist, geo, eps, dmax, dmap, lc, nlc)
      !differences of basis vectors that are substrate reciprocal vectors:
      !dmap(G_i-G_j) = index into lc (0: no substrate vector), lc(:,1)=0
      INTEGER, INTENT(IN)        :: glist(:, :)
      TYPE(t_igeo), INTENT(IN)   :: geo
      REAL, INTENT(IN)           :: eps
      INTEGER, INTENT(IN)        :: dmax
      INTEGER, INTENT(OUT)       :: dmap(-dmax:, -dmax:)
      INTEGER, INTENT(OUT)       :: lc(:, :)
      INTEGER, INTENT(OUT)       :: nlc
      INTEGER :: i, j, d(2), u(2), n
      REAL    :: w(2)
      dmap = -1
      lc = 0
      nlc = 1
      dmap(0, 0) = 1
      DO i = 1, SIZE(glist, 2)
         DO j = 1, SIZE(glist, 2)
            d = glist(:, i) - glist(:, j)
            IF (dmap(d(1), d(2)) >= 0) CYCLE
            w = MATMUL(REAL(d), geo%m)
            dmap(d(1), d(2)) = 0
            IF (ALL(ABS(w - NINT(w)) < eps)) THEN
               u = NINT(w)
               DO n = 1, nlc
                  IF (ALL(lc(:, n) == u)) EXIT
               END DO
               IF (n > nlc) THEN
                  nlc = nlc + 1
                  lc(:, nlc) = u
               END IF
               dmap(d(1), d(2)) = n
            END IF
         END DO
      END DO
      WHERE (dmap < 0) dmap = 0
   END SUBROUTINE priv_common

   SUBROUTINE priv_shells(geo, nshell, q)
      !substrate vectors u/=0 of the first nshell shells (by length)
      TYPE(t_igeo), INTENT(IN)            :: geo
      INTEGER, INTENT(IN)                 :: nshell
      INTEGER, ALLOCATABLE, INTENT(INOUT) :: q(:, :)
      INTEGER, PARAMETER :: nb = 12
      INTEGER :: i, j, n, ns, k
      REAL    :: len((2*nb + 1)**2), rad(nshell), l
      INTEGER :: cand(2, (2*nb + 1)**2)
      n = 0
      DO i = -nb, nb
         DO j = -nb, nb
            IF (i == 0 .AND. j == 0) CYCLE
            n = n + 1
            cand(:, n) = (/i, j/)
            len(n) = NORM2(MATMUL(REAL(cand(:, n)), geo%bb))
         END DO
      END DO
      !radii of the first nshell shells
      ns = 0
      DO k = 1, nshell
         l = HUGE(1.0)
         DO i = 1, n
            IF (ns > 0) THEN
               IF (len(i) <= rad(ns)*(1.0 + 1.0E-6)) CYCLE
            END IF
            l = MIN(l, len(i))
         END DO
         ns = ns + 1
         rad(ns) = l
      END DO
      DEALLOCATE (q)
      ALLOCATE (q(2, COUNT(len(:n) <= rad(ns)*(1.0 + 1.0E-6))))
      k = 0
      DO i = 1, n
         IF (len(i) <= rad(ns)*(1.0 + 1.0E-6)) THEN
            k = k + 1
            q(:, k) = cand(:, i)
         END IF
      END DO
   END SUBROUTINE priv_shells

   SUBROUTINE priv_ulist(lc, q, ulist)
      !all Umklapp vectors needed: Lambda_c, the shells q, -q and q-q'
      INTEGER, INTENT(IN)               :: lc(:, :), q(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: ulist(:, :)
      INTEGER, ALLOCATABLE :: t(:, :)
      INTEGER :: n, i, j
      ALLOCATE (t(2, SIZE(lc, 2) + 2*SIZE(q, 2) + SIZE(q, 2)**2))
      n = 0
      DO i = 1, SIZE(lc, 2)
         CALL add(lc(:, i))
      END DO
      DO i = 1, SIZE(q, 2)
         CALL add(q(:, i))
         CALL add(-q(:, i))
         DO j = 1, SIZE(q, 2)
            CALL add(q(:, i) - q(:, j))
         END DO
      END DO
      ulist = t(:, :n)
   CONTAINS
      SUBROUTINE add(u)
         INTEGER, INTENT(IN) :: u(2)
         INTEGER :: k
         DO k = 1, n
            IF (ALL(t(:, k) == u)) RETURN
         END DO
         n = n + 1
         t(:, n) = u
      END SUBROUTINE add
   END SUBROUTINE priv_ulist

   SUBROUTINE priv_symcheck(sym, geo, f, shift, l_comm)
      !every in-plane operation of the ad-layer must be a substrate operation
      !(rotated by theta); with common reciprocal vectors also the translations
      TYPE(t_sym), INTENT(IN)     :: sym
      TYPE(t_igeo), INTENT(IN)    :: geo
      TYPE(t_embfile), INTENT(IN) :: f
      REAL, INTENT(IN)            :: shift(2)
      LOGICAL, INTENT(IN)         :: l_comm
      INTEGER :: io, ib
      REAL    :: ro(2, 2), rb(2, 2), t(2), x(2)
      LOGICAL :: found
      IF (f%nsym == 0) THEN
         CALL juDFT_warn("The substrate file has no symmetry operations: symmetry check skipped", &
                         calledby="gf_embtransform")
         RETURN
      END IF
      DO io = 1, sym%nop
         ro = MATMUL(MATMUL(geo%ao, REAL(sym%mrot(1:2, 1:2, io))), priv_inv2(geo%ao))
         found = .FALSE.
         DO ib = 1, f%nsym
            rb = MATMUL(MATMUL(geo%ab, REAL(f%symop(:, :, ib))), priv_inv2(geo%ab))
            IF (ANY(ABS(rb - ro) > 1.0E-5)) CYCLE
            IF (l_comm) THEN
               !{R|tau_o} must equal the substrate operation {R|tau_b+d-Rd} modulo the substrate lattice
               t = MATMUL(geo%ao, sym%tau(1:2, io)) - MATMUL(geo%ab, f%symtau(:, ib)) - shift + MATMUL(ro, shift)
               x = MATMUL(priv_inv2(geo%ab), t)
               IF (ANY(ABS(x - NINT(x)) > 1.0E-5)) CYCLE
            END IF
            found = .TRUE.
            EXIT
         END DO
         IF (.NOT. found) THEN
            WRITE (oUnit, "(a,i0,a,4i3)") "ad-layer operation ", io, " is no substrate operation: ", sym%mrot(1:2, 1:2, io)
            CALL juDFT_error("The ad-layer k-set uses an operation that is no symmetry of the substrate", &
                             hint="remove it from the ad-layer (gfprep zmirror=f, sym.xml) or use the full k-mesh", &
                             calledby="gf_embtransform")
         END IF
      END DO
      WRITE (oUnit, "(a,i0,a)") "symmetry check: all ", sym%nop, " ad-layer operations are substrate operations"
   END SUBROUTINE priv_symcheck

   SUBROUTINE priv_umklapp(inp, tab, f, up, ie_up, ispin, geo, kap, glist, hl, bbmat, q, s0, seff, dnorm, nmiss)
      !Sigma_eff = Sigma0 - Delta_0Q [D_Q + Delta_QQ]^-1 Delta_Q0 for one kappa. The matrices have
      !nsp=tab%ns spin components in the order of the noco basis: index (s-1)*n + G
      TYPE(t_incomm), INTENT(IN)    :: inp
      TYPE(t_sigtab), INTENT(INOUT) :: tab
      TYPE(t_embfile), INTENT(IN)   :: f, up
      INTEGER, INTENT(IN)           :: ie_up, ispin
      TYPE(t_igeo), INTENT(IN)      :: geo
      REAL, INTENT(IN)              :: kap(2), bbmat(3, 3)
      INTEGER, INTENT(IN)           :: glist(:, :), hl(:), q(:, :)
      COMPLEX, INTENT(IN)           :: s0(:, :)
      COMPLEX, INTENT(OUT)          :: seff(:, :)
      REAL, INTENT(INOUT)           :: dnorm(4)
      INTEGER, INTENT(INOUT)        :: nmiss

      TYPE t_block
         INTEGER :: nb = 0, off = 0, ik = 0
         REAL    :: q(2)
         INTEGER, ALLOCATABLE :: g(:, :), gf(:), ih(:)
         COMPLEX, ALLOCATABLE :: d(:, :), du0(:, :), d0u(:, :), x(:, :)
         INTEGER, ALLOCATABLE :: ipiv(:)
      END TYPE t_block
      TYPE(t_block), ALLOCATABLE :: b(:)
      COMPLEX, ALLOCATABLE :: s2(:, :), s3(:, :), big(:, :), rhs(:, :), d0q(:, :), z(:, :), duv(:, :)
      INTEGER, ALLOCATABLE :: ipiv(:)
      INTEGER :: nq, iq, jq, nh, nsp, mh, ik, ih, nshift(2), n, m, nn, ntot, info
      REAL    :: s(3), p(2), anorm, rcond, rwork(1), th
      COMPLEX :: sg(2, 2)
      COMPLEX, ALLOCATABLE :: work(:)
      REAL, ALLOCATABLE    :: rw(:)
      REAL, EXTERNAL :: zlange

      nq = SIZE(q, 2)
      nh = SIZE(hl)
      nsp = tab%ns
      mh = nsp*nh
      th = inp%theta*pi_const/180.0
      ALLOCATE (b(nq))
      s(3) = 0.0
      !<-- blocks u: basis kappa+u+G with |kappa+u+G|<=rkmax
      ntot = 0
      DO iq = 1, nq
         b(iq)%q = kap + MATMUL(REAL(q(:, iq)), geo%minv)
         ik = up%find_k(b(iq)%q, 1.0E-5, nshift)
         IF (ik == 0) THEN
            WRITE (oUnit, "(a,2f12.6)") "missing k-point in upfile: ", b(iq)%q - NINT(b(iq)%q)
            CALL juDFT_error("k-point kappa+u not in upfile", hint="generate the k-list with write_kaux=t", &
                             calledby="gf_embtransform")
         END IF
         b(iq)%ik = ik
         !relative G: kappa+u+G = k_file + (G - nshift)
         ALLOCATE (b(iq)%g(2, up%ng), b(iq)%gf(up%ng), b(iq)%ih(up%ng))
         n = 0
         DO m = 1, up%ng
            s(1:2) = up%kpts(:, ik) + up%gvecs(:, m)
            IF (DOT_PRODUCT(s, MATMUL(s, bbmat)) > geo%rk**2) CYCLE
            n = n + 1
            b(iq)%g(:, n) = up%gvecs(:, m) + nshift
            b(iq)%gf(n) = m
            !the same relative G in the home block
            b(iq)%ih(n) = 0
            DO ih = 1, nh
               IF (ALL(glist(:, hl(ih)) == b(iq)%g(:, n))) b(iq)%ih(n) = ih
            END DO
         END DO
         b(iq)%nb = n
         b(iq)%off = ntot
         ntot = ntot + nsp*n
         ALLOCATE (b(iq)%d(nsp*n, nsp*n), b(iq)%du0(nsp*n, mh), b(iq)%d0u(mh, nsp*n), b(iq)%ipiv(nsp*n))
         CALL priv_dblock(iq, b(iq)%d)
         b(iq)%du0 = 0.0
         b(iq)%d0u = 0.0
         DO n = 1, b(iq)%nb
            IF (b(iq)%ih(n) == 0) THEN
               nmiss = nmiss + 1
               CYCLE
            END IF
            !Delta_u0 = sigma_u(p), Delta_0u = sigma_-u(p+u), p = (kappa+G)M
            p = MATMUL(kap + b(iq)%g(:, n), geo%m)
            CALL priv_sig(p, q(:, iq), sg)
            CALL priv_put(b(iq)%du0, n, b(iq)%nb, b(iq)%ih(n), nh, sg)
            CALL priv_sig(p + q(:, iq), -q(:, iq), sg)
            CALL priv_put(b(iq)%d0u, b(iq)%ih(n), nh, n, b(iq)%nb, sg)
         END DO
         dnorm(1) = MAX(dnorm(1), zlange("F", nsp*b(iq)%nb, mh, b(iq)%du0, MAX(nsp*b(iq)%nb, 1), rwork)/ &
                        MAX(zlange("F", mh, mh, s0, mh, rwork), 1.0E-30))
      END DO

      !<-- second order: -sum_u Delta_0u D_u^-1 Delta_u0
      ALLOCATE (s2(mh, mh), work(2*nsp*MAXVAL(b%nb) + 1), rw(2*nsp*MAXVAL(b%nb) + 1))
      s2 = 0.0
      DO iq = 1, nq
         n = nsp*b(iq)%nb
         IF (n == 0) CYCLE
         anorm = zlange("1", n, n, b(iq)%d, n, rwork)
         CALL zgetrf(n, n, b(iq)%d, n, b(iq)%ipiv, info)
         IF (info /= 0) CALL juDFT_error("singular D_u in the Umklapp correction", calledby="gf_embtransform")
         CALL zgecon("1", n, b(iq)%d, n, anorm, rcond, work, rw, info)
         IF (rcond < 1.0E-8) dnorm(4) = dnorm(4) + 1.0
         ALLOCATE (b(iq)%x(n, mh))
         b(iq)%x = b(iq)%du0
         CALL zgetrs("N", n, mh, b(iq)%d, n, b(iq)%ipiv, b(iq)%x, n, info)
         s2 = s2 - MATMUL(b(iq)%d0u, b(iq)%x)
      END DO
      dnorm(2) = MAX(dnorm(2), zlange("F", mh, mh, s2, mh, rwork))
      seff = s0 + s2
      IF (inp%order == 2) RETURN

      IF (inp%order == 3) THEN
         !<-- third order: + sum_{u/=u'} Delta_0u D_u^-1 Delta_uu' D_u'^-1 Delta_u'0
         ALLOCATE (s3(mh, mh))
         s3 = 0.0
         DO iq = 1, nq
            n = nsp*b(iq)%nb
            IF (n == 0) CYCLE
            !z = Delta_0u D_u^-1 = (D_u^-T Delta_0u^T)^T
            ALLOCATE (z(n, mh))
            z = TRANSPOSE(b(iq)%d0u)
            CALL zgetrs("T", n, mh, b(iq)%d, n, b(iq)%ipiv, z, n, info)
            DO jq = 1, nq
               IF (jq == iq .OR. b(jq)%nb == 0) CYCLE
               CALL priv_duv(iq, jq, duv)
               s3 = s3 + MATMUL(TRANSPOSE(z), MATMUL(duv, b(jq)%x))
            END DO
            DEALLOCATE (z)
         END DO
         seff = seff + s3
         dnorm(3) = MAX(dnorm(3), zlange("F", mh, mh, s3, mh, rwork))
         RETURN
      END IF

      !<-- full: solve (D_Q + Delta_QQ) X = Delta_Q0
      ALLOCATE (big(ntot, ntot), rhs(ntot, mh), d0q(mh, ntot), ipiv(ntot))
      big = 0.0
      DO iq = 1, nq
         n = nsp*b(iq)%nb
         IF (n == 0) CYCLE
         !D_u again (the LU factor replaced it)
         CALL priv_dblock(iq, big(b(iq)%off + 1:b(iq)%off + n, b(iq)%off + 1:b(iq)%off + n))
         rhs(b(iq)%off + 1:b(iq)%off + n, :) = b(iq)%du0
         d0q(:, b(iq)%off + 1:b(iq)%off + n) = b(iq)%d0u
         DO jq = 1, nq
            IF (jq == iq .OR. b(jq)%nb == 0) CYCLE
            CALL priv_duv(iq, jq, duv)
            nn = nsp*b(jq)%nb
            big(b(iq)%off + 1:b(iq)%off + n, b(jq)%off + 1:b(jq)%off + nn) = duv
         END DO
      END DO
      CALL zgesv(ntot, mh, big, ntot, ipiv, rhs, ntot, info)
      IF (info /= 0) CALL juDFT_error("singular Umklapp system", calledby="gf_embtransform")
      seff = s0 - MATMUL(d0q, rhs)
      dnorm(3) = MAX(dnorm(3), zlange("F", mh, mh, seff - s0 - s2, mh, rwork))

   CONTAINS

      SUBROUTINE priv_sig(pp, u, sg)
         !sigma_u(p) as nsp x nsp block, spin frame rotated with the lattice
         REAL, INTENT(IN)     :: pp(2)
         INTEGER, INTENT(IN)  :: u(2)
         COMPLEX, INTENT(OUT) :: sg(2, 2)
         COMPLEX :: v(4)
         LOGICAL :: lok
         v = CMPLX(0.0, 0.0)
         CALL tab%eval(f, pp, tab%uindex(u), v(:nsp**2), lok)
         IF (.NOT. lok) CALL priv_fail(pp, u, tab)
         sg = RESHAPE(v, (/2, 2/))
         IF (nsp == 2) THEN
            sg(2, 1) = sg(2, 1)*EXP(CMPLX(0.0, th))
            sg(1, 2) = sg(1, 2)*EXP(CMPLX(0.0, -th))
         END IF
      END SUBROUTINE priv_sig

      SUBROUTINE priv_put(mat, ir, nr, ic, nc, sg)
         !G-diagonal block sg at plane wave ir (of nr) and ic (of nc)
         COMPLEX, INTENT(INOUT) :: mat(:, :)
         INTEGER, INTENT(IN)    :: ir, nr, ic, nc
         COMPLEX, INTENT(IN)    :: sg(2, 2)
         INTEGER :: s1, s2
         DO s2 = 1, nsp
            DO s1 = 1, nsp
               mat((s1 - 1)*nr + ir, (s2 - 1)*nc + ic) = mat((s1 - 1)*nr + ir, (s2 - 1)*nc + ic) + sg(s1, s2)
            END DO
         END DO
      END SUBROUTINE priv_put

      SUBROUTINE priv_dblock(iq, d)
         !D_u = Sigma_up(kappa+u) + Sigma0 of the block (G-diagonal, incommensurate)
         INTEGER, INTENT(IN)  :: iq
         COMPLEX, INTENT(OUT) :: d(:, :)
         COMPLEX, ALLOCATABLE :: um(:, :)
         INTEGER, PARAMETER   :: comp(2, 2) = RESHAPE((/1, 4, 3, 2/), (/2, 2/))
         INTEGER :: s1, s2, k, nb, js
         nb = b(iq)%nb
         ALLOCATE (um(up%ng, up%ng))
         d = 0.0
         DO s2 = 1, nsp
            DO s1 = 1, nsp
               !noco file: spins 11,22,12,21; collinear file: the diagonal
               IF (up%spins == 4) THEN
                  js = comp(s1, s2)
               ELSE IF (s1 /= s2) THEN
                  CYCLE
               ELSE IF (nsp == 2) THEN
                  js = MIN(s1, up%spins)
               ELSE
                  js = MIN(ispin, up%spins)
               END IF
               CALL up%read(ie_up, js, b(iq)%ik, um)
               d((s1 - 1)*nb + 1:s1*nb, (s2 - 1)*nb + 1:s2*nb) = um(b(iq)%gf(:nb), b(iq)%gf(:nb))
            END DO
         END DO
         DO k = 1, nb
            CALL priv_sig(MATMUL(b(iq)%q + b(iq)%g(:, k), geo%m), (/0, 0/), sg)
            CALL priv_put(d, k, nb, k, nb, sg)
         END DO
      END SUBROUTINE priv_dblock

      SUBROUTINE priv_duv(iq, jq, duv)
         !Delta_uu'(G,G) = sigma_{u-u'}((kappa+u'+G)M), same relative G
         INTEGER, INTENT(IN)                 :: iq, jq
         COMPLEX, ALLOCATABLE, INTENT(INOUT) :: duv(:, :)
         INTEGER :: n1, n2
         IF (ALLOCATED(duv)) DEALLOCATE (duv)
         ALLOCATE (duv(nsp*b(iq)%nb, nsp*b(jq)%nb))
         duv = 0.0
         DO n2 = 1, b(jq)%nb
            DO n1 = 1, b(iq)%nb
               IF (ANY(b(iq)%g(:, n1) /= b(jq)%g(:, n2))) CYCLE
               CALL priv_sig(MATMUL(b(jq)%q + b(jq)%g(:, n2), geo%m), q(:, iq) - q(:, jq), sg)
               CALL priv_put(duv, n1, b(iq)%nb, n2, b(jq)%nb, sg)
            END DO
         END DO
      END SUBROUTINE priv_duv
   END SUBROUTINE priv_umklapp

   SUBROUTINE priv_write_kaux(gkpts, geo, q)
      !auxiliary k-points kappa+u (ad-layer cell, reduced) and their map
      TYPE(t_kpts), INTENT(IN) :: gkpts
      TYPE(t_igeo), INTENT(IN) :: geo
      INTEGER, INTENT(IN)      :: q(:, :)
      REAL, ALLOCATABLE    :: kall(:, :), kl(:, :)
      INTEGER, ALLOCATABLE :: id(:), gs(:, :)
      INTEGER :: nk, iq, n
      REAL    :: qq(2)
      ALLOCATE (kall(2, gkpts%nkpt*SIZE(q, 2)), gs(2, gkpts%nkpt*SIZE(q, 2)))
      n = 0
      DO nk = 1, gkpts%nkpt
         DO iq = 1, SIZE(q, 2)
            n = n + 1
            qq = gkpts%bk(1:2, nk) + MATMUL(REAL(q(:, iq)), geo%minv)
            gs(:, n) = NINT(qq)
            kall(:, n) = qq - gs(:, n)
         END DO
      END DO
      CALL priv_unique(kall, kl, id)
      OPEN (98, FILE="kaux_map", STATUS="replace", ACTION="write")
      WRITE (98, "(a)") "# kappa-index  u(1) u(2)  aux-index  Gs(1) Gs(2)   kappa+u = k_aux + Gs"
      n = 0
      DO nk = 1, gkpts%nkpt
         DO iq = 1, SIZE(q, 2)
            n = n + 1
            WRITE (98, "(i6,2i5,i8,2i5)") nk, q(:, iq), id(n), NINT(kall(:, n) + gs(:, n) - kl(:, id(n)))
         END DO
      END DO
      CLOSE (98)
      CALL priv_klist_write("kpts_aux", kl)
      WRITE (oUnit, "(a,i0,a)") "incomm: ", SIZE(kl, 2), " auxiliary k-points written to kpts_aux (map in kaux_map)"
   END SUBROUTINE priv_write_kaux

   SUBROUTINE priv_write_ksub(gkpts, geo, q, glist)
      !substrate k-points [(kappa+G)M] needed for exact lookups
      TYPE(t_kpts), INTENT(IN) :: gkpts
      TYPE(t_igeo), INTENT(IN) :: geo
      INTEGER, INTENT(IN)      :: q(:, :), glist(:, :)
      REAL, ALLOCATABLE    :: kall(:, :), kl(:, :), t(:, :)
      INTEGER, ALLOCATABLE :: id(:)
      INTEGER :: nk, iq, n, ig
      REAL    :: p(2), kk(2), s(2)
      ALLOCATE (kall(2, 1024))
      n = 0
      DO nk = 1, gkpts%nkpt
         DO iq = 0, SIZE(q, 2)
            kk = gkpts%bk(1:2, nk)
            IF (iq > 0) kk = kk + MATMUL(REAL(q(:, iq)), geo%minv)
            kk = kk - NINT(kk)
            DO ig = 1, SIZE(glist, 2)
               s = kk + glist(:, ig)
               IF (NORM2(MATMUL(s, geo%bo)) > geo%rk) CYCLE
               p = MATMUL(s, geo%m)
               IF (n == SIZE(kall, 2)) THEN
                  ALLOCATE (t(2, 2*n))
                  t(:, :n) = kall
                  CALL MOVE_ALLOC(t, kall)
               END IF
               n = n + 1
               kall(:, n) = p - NINT(p)
            END DO
         END DO
      END DO
      CALL priv_unique(kall(:, :n), kl, id)
      CALL priv_klist_write("kpts_sub", kl)
      WRITE (oUnit, "(a,i0,a)") "incomm: ", SIZE(kl, 2), " substrate k-points written to kpts_sub (substrate cell)"
   END SUBROUTINE priv_write_ksub

   SUBROUTINE priv_unique(kin, kl, id)
      !distinct k-points modulo reciprocal vectors (canonical in [-1/2,1/2)),
      !kin(:,n) is equivalent to kl(:,id(n))
      REAL, INTENT(IN)                  :: kin(:, :)
      REAL, ALLOCATABLE, INTENT(OUT)    :: kl(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
      REAL, PARAMETER :: tol = 1.0E-8
      REAL, ALLOCATABLE    :: kc(:, :), ku(:, :)
      INTEGER, ALLOCATABLE :: idx(:)
      INTEGER :: n, i, j, nu
      n = SIZE(kin, 2)
      ALLOCATE (kc(2, n), ku(2, n), idx(n), id(n))
      DO i = 1, n
         kc(:, i) = gf_kcanon(kin(:, i), tol)
         idx(i) = i
      END DO
      CALL priv_mergesort(kc(1, :), idx)
      nu = 0
      DO i = 1, n
         id(idx(i)) = 0
         !earlier distinct points with k(1) within tol
         DO j = nu, 1, -1
            IF (kc(1, idx(i)) - ku(1, j) > tol) EXIT
            IF (ABS(kc(2, idx(i)) - ku(2, j)) < tol) THEN
               id(idx(i)) = j
               EXIT
            END IF
         END DO
         IF (id(idx(i)) == 0) THEN
            nu = nu + 1
            ku(:, nu) = kc(:, idx(i))
            id(idx(i)) = nu
         END IF
      END DO
      kl = ku(:, :nu)
   END SUBROUTINE priv_unique

   RECURSIVE SUBROUTINE priv_mergesort(key, idx)
      !sort idx by key(idx)
      REAL, INTENT(IN)       :: key(:)
      INTEGER, INTENT(INOUT) :: idx(:)
      INTEGER, ALLOCATABLE :: t(:)
      INTEGER :: n, m, i, j, k
      n = SIZE(idx)
      IF (n < 2) RETURN
      m = n/2
      CALL priv_mergesort(key, idx(:m))
      CALL priv_mergesort(key, idx(m + 1:))
      ALLOCATE (t(n))
      i = 1; j = m + 1; k = 0
      DO WHILE (i <= m .AND. j <= n)
         k = k + 1
         IF (key(idx(i)) <= key(idx(j))) THEN
            t(k) = idx(i); i = i + 1
         ELSE
            t(k) = idx(j); j = j + 1
         END IF
      END DO
      DO WHILE (i <= m)
         k = k + 1; t(k) = idx(i); i = i + 1
      END DO
      DO WHILE (j <= n)
         k = k + 1; t(k) = idx(j); j = j + 1
      END DO
      idx = t
   END SUBROUTINE priv_mergesort

   SUBROUTINE priv_klist_write(name, kl)
      !k-list file read by gfprep (kptsfile): count, then one k-point per line
      CHARACTER(LEN=*), INTENT(IN) :: name
      REAL, INTENT(IN)             :: kl(:, :)
      INTEGER :: n
      OPEN (98, FILE=name, STATUS="replace", ACTION="write")
      WRITE (98, "(i0)") SIZE(kl, 2)
      DO n = 1, SIZE(kl, 2)
         WRITE (98, "(2f20.14)") kl(:, n)
      END DO
      CLOSE (98)
   END SUBROUTINE priv_klist_write

   SUBROUTINE priv_copy_boundary(f, fid, ng2, rot, shift)
      !Coulomb boundary data of the substrate file: only the 2D star G=0 is
      !kept, the other stars are set to zero (screened by the metal)
      TYPE(t_embfile), INTENT(IN) :: f
      INTEGER(HID_T), INTENT(IN)  :: fid
      INTEGER, INTENT(IN)         :: ng2
      REAL, INTENT(IN)            :: rot(2, 2), shift(2)
      INTEGER(HID_T) :: gid, gout, did
      INTEGER :: n2, hdferr, dims(3), na, n
      REAL    :: pos, c, dt
      REAL, ALLOCATABLE :: br(:), bi(:), outr(:), outi(:), ps(:, :, :), pso(:, :, :), rmt(:), atpos(:, :)

      IF (io_groupexists(f%fid, "Boundary")) THEN
         CALL io_gopen(f%fid, "Boundary", gid, hdferr)
         ALLOCATE (outr(ng2), outi(ng2))
         outr = 0.0; outi = 0.0
         CALL io_gcreate(fid, "Boundary", gout, hdferr)
         IF (io_attexists(gid, "n2")) THEN
            CALL io_read_att(gid, "n2", n2)
            ALLOCATE (br(n2), bi(n2))
            CALL io_read_att(gid, "bulk_real", br)
            CALL io_read_att(gid, "bulk_imag", bi)
            CALL io_read_att(gid, "pos", pos)
            outr(1) = br(1); outi(1) = bi(1)
            CALL io_write_att(gout, "n2", ng2)
            CALL io_write_att(gout, "bulk_real", outr)
            CALL io_write_att(gout, "bulk_imag", outi)
            CALL io_write_att(gout, "pos", pos)
         END IF
         IF (io_dataexists(gid, "pseudo")) THEN
            CALL io_dopen(gid, "pseudo", did, hdferr)
            CALL io_datadim(did, dims)
            ALLOCATE (ps(dims(1), dims(2), dims(3)), pso(dims(1), dims(2), ng2))
            CALL io_read(did, (/1, 1, 1/), dims, "pseudo", ps)
            CALL io_dclose(did, hdferr)
            pso = 0.0
            pso(:, :, 1) = ps(:, :, 1)
            CALL io_write_var(gout, "pseudo", pso)
            CALL io_read_att(gid, "c", c)
            CALL io_read_att(gid, "dt", dt)
            CALL io_write_att(gout, "c", c)
            CALL io_write_att(gout, "dt", dt)
         END IF
         IF (io_dataexists(gid, "rmt")) THEN
            CALL io_dopen(gid, "rmt", did, hdferr)
            CALL io_datadim(did, dims(1:1))
            CALL io_dclose(did, hdferr)
            na = dims(1)
            ALLOCATE (rmt(na), atpos(3, na))
            IF (na > 0) THEN
               CALL io_read_var(gid, "rmt", rmt)
               CALL io_read_var(gid, "atpos", atpos)
               !substrate frame -> ad-layer frame
               DO n = 1, na
                  atpos(1:2, n) = MATMUL(rot, atpos(1:2, n)) + shift
               END DO
               CALL juDFT_warn("Substrate atoms cut the embedding plane: their lateral positions are not "// &
                               "commensurate with the ad-layer", calledby="gf_embtransform")
            END IF
            CALL io_write_var(gout, "rmt", rmt)
            CALL io_write_var(gout, "atpos", atpos)
         END IF
         CALL io_gclose(gout, hdferr)
         CALL io_gclose(gid, hdferr)
         WRITE (oUnit, "(a)") "incomm: Boundary data copied (2D star G=0 only)"
      END IF
      IF (io_groupexists(f%fid, "CoulBoundary")) THEN
         CALL io_gopen(f%fid, "CoulBoundary", gid, hdferr)
         CALL io_read_att(gid, "n", n2)
         IF (ALLOCATED(br)) DEALLOCATE (br, bi)
         ALLOCATE (br(n2), bi(n2))
         CALL io_read_var(gid, "v_coul_real", br)
         CALL io_read_var(gid, "v_coul_aimag", bi)
         CALL io_gclose(gid, hdferr)
         IF (.NOT. ALLOCATED(outr)) ALLOCATE (outr(ng2), outi(ng2))
         outr = 0.0; outi = 0.0
         outr(1) = br(1); outi(1) = bi(1)
         CALL io_gcreate(fid, "CoulBoundary", gout, hdferr)
         CALL io_write_att(gout, "n", ng2)
         CALL io_write_var(gout, "v_coul_real", outr)
         CALL io_write_var(gout, "v_coul_aimag", outi)
         CALL io_gclose(gout, hdferr)
      END IF
   END SUBROUTINE priv_copy_boundary

   FUNCTION priv_inv2(a) RESULT(b)
      REAL, INTENT(IN) :: a(2, 2)
      REAL             :: b(2, 2)
      REAL :: det
      det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
      IF (ABS(det) < 1.0E-12) CALL juDFT_error("singular 2x2 matrix", calledby="gf_embtransform")
      b(1, 1) = a(2, 2)/det; b(2, 2) = a(1, 1)/det
      b(1, 2) = -a(1, 2)/det; b(2, 1) = -a(2, 1)/det
   END FUNCTION priv_inv2

END MODULE m_gf_embtransform
