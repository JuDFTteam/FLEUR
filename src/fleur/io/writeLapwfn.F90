!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
! Contributed by Mike Towler (fleurqmc): writer of lapwfn.data/pwfn.data for CASINO.
! Copyright (c) 2026 Mike Towler. MIT licence; see src/contributed/CASINO/LICENSE.
!--------------------------------------------------------------------------------
MODULE m_writeLapwfn
   !! The all-electron LAPW orbitals of a converged ground state in the code-neutral
   !! format read by CASINO (lapwfn.data, format version 1, src/contributed/CASINO/SPEC.md)
   !! with the interstitial plane waves in pwfn.data. Switched on by output/lapwfn.
   !! Sphere coefficients are formed with FLEUR's matching in the global frame (a copy
   !! of sym with ngopr=1, invsat=0) times i^l; radial functions from t_radfun, core
   !! states from differ; interstitial coefficients are the eigenvector over sqrt(Omega).
   !! At k=G/2 bands are made real (one phase, or real combinations of a degenerate
   !! level). Output files are written as .part and renamed when complete; a refusal
   !! deletes them. With lapwfn_points.in FLEUR's own orbital values at its points are
   !! written to lapwfn_points.out, the reference for a reader.
   USE m_juDFT
   USE m_constants
   USE m_types_atoms
   USE m_types_input
   USE m_types_sym
   USE m_types_cell
   USE m_types_kpts
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_banddos
   USE m_types_enpara
   USE m_types_potden
   USE m_types_misc
   USE m_types_mpi
   USE m_types_lapw
   USE m_types_mat
   USE m_types_radfun
   USE m_types_abc
   USE m_types_xcpot
   USE m_types_fleurinput
   USE m_eig66_io, ONLY: read_eig
   USE m_cdn_io, ONLY: readPrevEFermi, getIOMode
   USE m_cdnpot_io_hdf
   USE m_abcoeff_store, ONLY: l_use_abcoeff_store
   USE m_compile_descr, ONLY: get_compile_desc
   USE m_differ
   USE m_ylm
#ifdef CPP_HDF
   USE hdf5
#endif
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_CHAR, C_NULL_CHAR
   USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int64
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: write_lapwfn

   CHARACTER(LEN=*), PARAMETER :: writer_version = "fleurqmc 0.3"
   INTEGER, PARAMETER :: nfmaxd = 10          ! radial functions per l at most
   REAL, PARAMETER :: c_nonrel = 1.0e5        ! a smaller speed of light is a relativistic calculation
   REAL, PARAMETER :: tail_floor = 1.0e-12    ! the core mesh ends where the tail is below this of max u
   REAL, PARAMETER :: real_tol = 1.0e-8       ! the two tests of a real orbital
   REAL, PARAMETER :: inv_tol = 1.0e-8        ! the inversion relation's residual
   REAL, PARAMETER :: lev_tol = 1.0e-6        ! consecutive eigenvalues within this are one degenerate level
   REAL, PARAMETER :: clos_tol = 1.0e-6       ! a level closes under conjugation when K V lies in its span to this
   REAL, PARAMETER :: ef_tol = 1.0e-6         ! a band at or below efermi + ef_tol counts as occupied
   CHARACTER(LEN=*), PARAMETER :: part = ".part"   ! output files are written under this suffix and renamed
   INTEGER, PARAMETER :: nfiles = 4
   CHARACTER(LEN=24), PARAMETER :: fname(nfiles) = [CHARACTER(LEN=24) :: "lapwfn.data", "pwfn.data", &
                                                    "lapwfn_measurements.out", "lapwfn_points.out"]
   INTEGER :: funit(nfiles) = 0
   LOGICAL :: fopen(nfiles) = .FALSE.
   LOGICAL :: l_realise = .TRUE.              ! output/lapwfn/@realise

   INTERFACE
      FUNCTION c_rename(old, new) BIND(C, NAME="rename") RESULT(r)
         IMPORT :: C_INT, C_CHAR
         IMPLICIT NONE
         INTEGER(C_INT) :: r
         CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) :: old, new
      END FUNCTION c_rename
   END INTERFACE

   TYPE corestate_t
      INTEGER :: n, l, nc
      REAL :: e, occ, kappa, rj, a, c2, c3, leak_tail, leak_fleur, norm, small
      REAL, ALLOCATABLE :: r(:), u(:), u_raw(:)
   END TYPE corestate_t

   INTEGER :: ulapw, upw, umeas, upts, nkdone
   INTEGER :: ngtot, gmin(3), gmax(3)
   INTEGER, ALLOCATABLE :: gall(:,:), gindex(:,:,:)
   INTEGER, ALLOCATABLE :: nfl(:,:), atype(:)
   REAL, ALLOCATABLE :: utab(:,:,:,:)
   REAL, ALLOCATABLE :: rjoin(:,:,:)
   INTEGER :: npts
   REAL, ALLOCATABLE :: pts(:,:)
   REAL :: efermi, etot_in
   CHARACTER(LEN=200) :: title_in, func_in

CONTAINS

   SUBROUTINE write_lapwfn(fi, xcpot, nococonv, enpara, vTot, results, fmpi, eig_id)
      !! Driver, called on rank 0 after eigen when output/lapwfn is set: every k point of the
      !! selected list, all bands computed.
      TYPE(t_fleurinput), INTENT(IN) :: fi
      CLASS(t_xcpot), INTENT(IN)     :: xcpot
      TYPE(t_nococonv), INTENT(IN)   :: nococonv
      TYPE(t_enpara), INTENT(IN)     :: enpara
      TYPE(t_potden), INTENT(IN)     :: vTot
      TYPE(t_results), INTENT(IN)    :: results
      TYPE(t_mpi), INTENT(IN)        :: fmpi
      INTEGER, INTENT(IN)            :: eig_id

      TYPE(t_radfun), ALLOCATABLE :: radfun(:)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat)  :: zMat
      INTEGER :: n, ik, nbands, nbasfcn
      LOGICAL :: l_real

      CALL timestart("write lapwfn")
      ! the cached matching coefficients belong to the local frames
      l_use_abcoeff_store = .FALSE.
      l_realise = fi%banddos%lapwfn_realise
      CALL lapwfn_refuse(fi%input, fi%atoms, fi%noco, fi%banddos, results, fmpi, fi%sym)
      ALLOCATE (radfun(fi%atoms%ntype))
      DO n = 1, fi%atoms%ntype
         CALL radfun(n)%generate_radial_functions(fi%atoms, fi%input, enpara, fmpi, vTot, n)
      END DO
      CALL lapwfn_functional(xcpot, func_in)
      CALL lapwfn_init(fi%input, fi%atoms, fi%sym, fi%cell, fi%kpts, fi%noco, nococonv, fi%banddos, enpara, &
                       radfun, vTot, fmpi, 1)
      WRITE (umeas, "(a,es16.8,a,es16.8)") "last_distance ", results%last_distance, " minDistance ", fi%input%minDistance
      nkdone = 0
      l_real = fi%sym%invs .AND. (.NOT.fi%noco%l_soc) .AND. (.NOT.fi%noco%l_noco) .AND. fi%atoms%n_hia == 0
      DO ik = 1, fi%kpts%nkpt
         CALL lapw%init(fi%input, fi%noco, nococonv, fi%kpts, fi%atoms, fi%sym, ik, fi%cell, fmpi)
         nbands = results%neig(ik, 1)
         nbasfcn = lapw%nv(1) + fi%atoms%nlotot
         CALL zMat%init(l_real, nbasfcn, nbands)
         CALL read_eig(eig_id, ik, 1, zmat=zMat)
         CALL lapwfn_kpoint(ik, fi%input, fi%atoms, fi%sym, fi%cell, lapw, zMat, results%eig(1:nbands, ik, 1), &
                            nbands, radfun, fi%noco, nococonv, 1)
      END DO
      CALL lapwfn_finish()
      CALL timestop("write lapwfn")
   END SUBROUTINE write_lapwfn

   SUBROUTINE lapwfn_kpoint(ikpt, input, atoms, sym, cell, lapw, zMat, eig, noccbd, radfun, noco, nococonv, jspin)
      INTEGER, INTENT(IN) :: ikpt, noccbd, jspin
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_lapw), INTENT(IN) :: lapw
      TYPE(t_mat), INTENT(IN) :: zMat
      REAL, INTENT(IN) :: eig(:)
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv

      TYPE(t_sym) :: sym_l
      TYPE(t_abc) :: abcg, abcl
      INTEGER :: n, na, ia, l, m, lm, j, ib, ig, nv, idx, nlmd
      INTEGER :: twok(3)
      INTEGER :: nlev, nw, ilv, i1, i2, nl, nok, nocl, nocl0
      INTEGER, ALLOCATABLE :: lev1(:), lev2(:)
      COMPLEX, ALLOCATABLE :: cof(:,:,:,:), cofl(:,:,:,:), cpw(:,:), tmat(:,:), phs(:)
      COMPLEX :: il
      LOGICAL :: l_halfg, ok
      LOGICAL, ALLOCATABLE :: lreal(:)
      REAL, ALLOCATABLE :: eigw(:), devp(:), devm(:)
      REAL :: kc(3), spread, closure, dmm, dms, sn, sn1, dcc, dlab
      CHARACTER(LEN=200) :: str
      CHARACTER(LEN=80) :: why

      IF (ANY(sym%invsat == 2) .AND. .NOT.zMat%l_real) then
        CALL lapwfn_error("lapwfn: inversion partners with complex eigenvectors")
      end if

      nv = lapw%nv(1)
      nlmd = atoms%lmaxd*(atoms%lmaxd+2)
      ALLOCATE (cof(noccbd, 0:nlmd, nfmaxd, atoms%nat), cpw(ngtot, noccbd), lreal(noccbd), phs(noccbd))
      phs = CMPLX(1.0, 0.0)
      cof = CMPLX(0.0, 0.0)
      cpw = CMPLX(0.0, 0.0)

      ! the sphere coefficients in the global frame: FLEUR's matching under the identity frame;
      ! inversion partners are computed directly (their LO data point to the partner's z entries)
      sym_l = sym
      sym_l%ngopr = 1
      sym_l%invsat = 0
      DO n = 1, atoms%ntype
         CALL abcg%init(input, atoms, noccbd, n)
         CALL abcg%calc_abc(input, atoms, sym_l, cell, lapw, noccbd, radfun(n), noco, nococonv, jspin, n, zMat)
         DO ia = 1, atoms%neq(n)
            na = atoms%firstAtom(n) + ia - 1
            DO l = 0, atoms%lmax(n)
               il = CMPLX(0.0, 1.0)**l
               DO m = -l, l
                  lm = l*(l+1) + m
                  DO j = 1, nfl(l, n)
                     cof(:, lm, j, na) = il*abcg%cof(:, lm, j, ia)
                  END DO
               END DO
            END DO
         END DO
         IF (ANY(sym%invsat(atoms%firstAtom(n):atoms%firstAtom(n)+atoms%neq(n)-1) == 2)) then
           CALL lapwfn_inversion(n, atoms, sym, lapw, noccbd, ikpt, cof)
         end if
      END DO

      ! the interstitial coefficients: the eigenvector over sqrt(Omega), on the union of G vectors
      DO ig = 1, nv
         idx = gindex(lapw%gvec(1, ig, 1), lapw%gvec(2, ig, 1), lapw%gvec(3, ig, 1))
         IF (idx <= 0) CALL lapwfn_error("lapwfn: a G vector is missing from the union")
         IF (zMat%l_real) THEN
            cpw(idx, :) = CMPLX(zMat%data_r(ig, 1:noccbd), 0.0)/SQRT(cell%omtil)
         ELSE
            cpw(idx, :) = zMat%data_c(ig, 1:noccbd)/SQRT(cell%omtil)
         END IF
      END DO

      ! the levels: consecutive eigenvalues within lev_tol are one degenerate level
      ALLOCATE (lev1(noccbd), lev2(noccbd), eigw(noccbd), devp(noccbd), devm(noccbd), tmat(noccbd, noccbd))
      nlev = 1
      lev1(1) = 1
      DO ib = 2, noccbd
         IF (eig(ib) - eig(ib - 1) > lev_tol) THEN
            lev2(nlev) = ib - 1
            nlev = nlev + 1
            lev1(nlev) = ib
         END IF
      END DO
      lev2(nlev) = noccbd
      ! the list ends at the last level boundary below the highest band computed, whose level may be incomplete
      IF (nlev < 2) THEN
         WRITE (str, "(a,i0,a,i0,a)") "lapwfn: the ", noccbd, " bands computed at k point ", ikpt, &
            " form one degenerate level, which cannot be known complete; raise numbands"
         CALL lapwfn_error(TRIM(str))
      END IF
      nw = lev2(nlev - 1)
      WRITE (umeas, "(a,i5,a,i0,a,i0,a,i0,a,i0,a)") "dropped ", ikpt, " bands ", nw + 1, " to ", noccbd, " of the ", &
         noccbd, " computed (the level of the highest computed band may be incomplete); ", nw, " bands written"
      IF (COUNT(eig(1:nw) <= efermi + ef_tol) == nw) THEN
         WRITE (str, "(a,i0,a,i0,a)") "lapwfn: all ", nw, " bands written at k point ", ikpt, &
            " lie at or below the Fermi energy; raise numbands so that a band above it remains at every k point"
         CALL lapwfn_error(TRIM(str))
      END IF
      IF (nw < NINT(input%zelec/2.0)) THEN
         WRITE (str, "(a,i0,a,i0,a)") "lapwfn: only ", nw, " bands written at k point ", ikpt, &
            ", fewer than the occupied ones; raise numbands"
         CALL lapwfn_error(TRIM(str))
      END IF

      ! real orbitals at k = G/2: one phase per band where it makes the band real, else the real orthonormal
      ! combinations of a degenerate level not divided by the Fermi energy; tmat holds each written band as a
      ! combination of FLEUR's, for the probe
      eigw(1:noccbd) = eig(1:noccbd)
      tmat = CMPLX(0.0, 0.0)
      DO ib = 1, noccbd
         tmat(ib, ib) = CMPLX(1.0, 0.0)
      END DO
      lreal = .FALSE.
      phs = CMPLX(1.0, 0.0)
      devp = 0.0
      devm = 0.0
      l_halfg = ALL(ABS(2.0*lapw%bkpt - NINT(2.0*lapw%bkpt)) < 1.0e-8)
      twok = NINT(2.0*lapw%bkpt)
      IF (l_halfg) THEN
         DO ilv = 1, nlev - 1
            i1 = lev1(ilv)
            i2 = lev2(ilv)
            nl = i2 - i1 + 1
            DO ib = i1, i2
               CALL lapwfn_phase(ib, twok, cof, cpw, lreal(ib), phs(ib), devp(ib), devm(ib))
               tmat(ib, ib) = phs(ib)
            END DO
            nok = COUNT(lreal(i1:i2))
            nocl = COUNT(eig(i1:i2) <= efermi + ef_tol)
            nocl0 = COUNT(eig(i1:i2) <= efermi)
            spread = eig(i2) - eig(i1)
            closure = 0.0; dmm = 0.0; dms = 0.0; sn = 0.0; sn1 = 0.0; dcc = 0.0; dlab = 0.0
            ok = .FALSE.
            IF (nl == 1) THEN
               why = "single"
            ELSE IF (nok == nl) THEN
               why = 'real band by band in FLEUR''s own vectors'
            ELSE IF (.NOT.l_realise) THEN
               why = 'combinations switched off (realise="F"): left as FLEUR gives it'
            ELSE IF (nocl > 0 .AND. nocl < nl) THEN
               why = "divided by the Fermi energy: left as FLEUR gives it"
            ELSE
               CALL lapwfn_realise(i1, i2, twok, cof, cpw, eig, eigw, tmat, phs, lreal, devp, devm, &
                                   closure, dmm, dms, sn, sn1, dcc, dlab, ok, why)
            END IF
            WRITE (umeas, "(a,i5,2i5,i3,es10.2,2i3,7es10.2,l3,1x,a)") "level ", ikpt, i1, i2, nl, spread, nocl, nocl0, &
               closure, dmm, dms, sn, sn1, dcc, dlab, ok, TRIM(why)
         END DO
         DO ib = 1, nw
            WRITE (umeas, "(a,i5,i5,2l3,2es12.3)") "phase ", ikpt, ib, l_halfg, lreal(ib), devp(ib), devm(ib)
         END DO
      ELSE
         WRITE (umeas, "(a,i5,a)") "phase ", ikpt, " not at G/2: every band complex"
      END IF

      ! the band blocks
      kc = MATMUL(TRANSPOSE(cell%bmat), lapw%bkpt)
      WRITE (ulapw, "(a)") "k-point # ; # of bands (up spin/down spin) ; k-point coords (au)"
      WRITE (ulapw, "(3i6,3es24.15)") ikpt, nw, 0, kc
      WRITE (upw, "(a)") " k-point # ; # of bands (up spin/down spin) ; k-point coords (au)"
      WRITE (upw, "(1x,3i6,3es24.15)") ikpt, nw, 0, kc
      DO ib = 1, nw
         WRITE (ulapw, "(a)") "Band, spin, eigenvalue (au), real (T/F)"
         WRITE (ulapw, "(2i6,es24.15,l3)") ib, 1, eigw(ib), lreal(ib)
         DO na = 1, atoms%nat
            n = atype(na)
            WRITE (ulapw, "(a,i0,a)") "Sphere coefficients of atom ", na, &
               ": l, m, then (re, im) for each radial function of that l"
            DO l = 0, atoms%lmax(n)
               DO m = -l, l
                  lm = l*(l+1) + m
                  WRITE (ulapw, "(2i4,*(1x,es23.15))") l, m, (REAL(cof(ib, lm, j, na)), AIMAG(cof(ib, lm, j, na)), j=1, nfl(l, n))
               END DO
            END DO
         END DO
         WRITE (upw, "(a)") " Band, spin, eigenvalue (au)"
         WRITE (upw, "(1x,2i6,es24.15)") ib, 1, eigw(ib)
         WRITE (upw, "(a)") " Eigenvector coefficients"
         DO ig = 1, ngtot
            WRITE (upw, '(" (",es23.15,",",es23.15,")")') REAL(cpw(ig, ib)), AIMAG(cpw(ig, ib))
         END DO
      END DO

      ! the probe: FLEUR's own local-frame coefficients at the listed points, combined as the written bands are
      IF (npts > 0) THEN
         ALLOCATE (cofl(noccbd, 0:nlmd, nfmaxd, atoms%nat))
         cofl = CMPLX(0.0, 0.0)
         DO n = 1, atoms%ntype
            CALL abcl%init(input, atoms, noccbd, n)
            CALL abcl%calc_abc(input, atoms, sym, cell, lapw, noccbd, radfun(n), noco, nococonv, jspin, n, zMat)
            DO ia = 1, atoms%neq(n)
               na = atoms%firstAtom(n) + ia - 1
               DO l = 0, atoms%lmax(n)
                  il = CMPLX(0.0, 1.0)**l
                  DO m = -l, l
                     lm = l*(l + 1) + m
                     DO j = 1, nfl(l, n)
                        cofl(:, lm, j, na) = il*abcl%cof(:, lm, j, ia)
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CALL lapwfn_probe(ikpt, atoms, sym, cell, lapw, noccbd, nw, cofl, cpw, tmat)
         DEALLOCATE (cofl)
      END IF

      nkdone = nkdone + 1
      DEALLOCATE (cof, cpw, lreal, phs, lev1, lev2, eigw, devp, devm, tmat)
   END SUBROUTINE lapwfn_kpoint

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_refuse(input, atoms, noco, banddos, results, fmpi, sym)
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_banddos), INTENT(IN) :: banddos
      TYPE(t_results), INTENT(IN) :: results
      TYPE(t_mpi), INTENT(IN) :: fmpi
      TYPE(t_sym), INTENT(IN) :: sym
      CHARACTER(LEN=300) :: str
      INTEGER :: n, l, lo, nsp
      LOGICAL :: l_err
      IF (input%jspins /= 1) CALL lapwfn_error("lapwfn: spin-polarised run; format version 1 is unpolarised")
      CALL lapwfn_stored_spins(nsp, l_err)
      IF (.NOT.l_err .AND. nsp /= input%jspins) THEN
         WRITE (str, "(a,i0,a,i0,a)") "lapwfn: the stored density (cdn.hdf) was written with ", nsp, &
            " spins where the input has ", input%jspins, "; the state as stored is what the writer takes"
         CALL lapwfn_error(TRIM(str))
      END IF
      IF (noco%l_noco .OR. noco%l_ss) CALL lapwfn_error("lapwfn: non-collinear or spin-spiral run; the eigenvector is a spinor")
      IF (noco%l_soc) CALL lapwfn_error("lapwfn: spin-orbit coupling; second-variational spinors are not in the format")
      IF (input%film) CALL lapwfn_error("lapwfn: film geometry; the vacuum region is not in the format")
      DO n = 1, atoms%ntype
         DO l = 0, atoms%lmax(n)
            IF (atoms%l_apw(l, n)) THEN
               WRITE (str, "(a,i0,a,i0,a)") "lapwfn: channel l=", l, " of type ", n, &
                  " is APW (an APW+lo local orbital or l > lmaxAPW); format version 1 is LAPW only"
               CALL lapwfn_error(TRIM(str))
            END IF
         END DO
         DO lo = 1, atoms%nlo(n)
            IF (atoms%l_relLO(lo, n)) THEN
               WRITE (str, "(a,i0,a,i0,a)") "lapwfn: local orbital ", lo, " of type ", n, &
                  " is a relativistic (Dirac) local orbital; the file carries non-relativistic orbitals"
               CALL lapwfn_error(TRIM(str))
            END IF
         END DO
      END DO
      IF (c_light(1.0) < c_nonrel) THEN
         WRITE (str, "(a,es10.3,a)") "lapwfn: this calculation is relativistic (c = ", c_light(1.0), &
            '); the file carries non-relativistic orbitals: set expertModes/@warp_factor="10000.0" for the SCF and this run'
         CALL lapwfn_error(TRIM(str))
      END IF
      IF (banddos%lapwfn_etot == HUGE(1.0)) CALL lapwfn_error("lapwfn: output/lapwfn/@totalEnergy is missing; "// &
         'the file carries the converged calculation''s total energy (totalEnergy in its out.xml)')
      IF (results%last_distance < 0.0) CALL lapwfn_error("lapwfn: no converged density is stored (no SCF distance)")
      IF (results%last_distance > input%minDistance) THEN
         WRITE (str, "(a,es10.3,a,es10.3,a)") "lapwfn: the stored density is not converged (distance ", results%last_distance, &
            " above minDistance ", input%minDistance, ")"
         CALL lapwfn_error(TRIM(str))
      END IF
      IF (fmpi%isize > 1) CALL lapwfn_error("lapwfn: one MPI process only")
      IF (ABS(input%zelec - NINT(input%zelec)) > 1.0e-8) then
        CALL lapwfn_error("lapwfn: a fractional number of valence electrons; the file carries an integer count")
      end if
      IF (input%frcor) CALL lapwfn_error('lapwfn: frozen core; the core would not be the converged potential''s')
      IF (input%kcrel /= 0) CALL lapwfn_error("lapwfn: kcrel must be 0")
      IF (atoms%n_hia > 0) CALL lapwfn_error("lapwfn: Hubbard-1 states are not handled")
      IF (ANY(sym%mrot(:, :, 1) /= RESHAPE([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])) .OR. ANY(ABS(sym%tau(:, 1)) > 1.0e-10)) then
        CALL lapwfn_error("lapwfn: operation 1 is not the identity")
      end if
   END SUBROUTINE lapwfn_refuse

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_init(input, atoms, sym, cell, kpts, noco, nococonv, banddos, enpara, radfun, vTot, fmpi, jspin)
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_noco), INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_banddos), INTENT(IN) :: banddos
      TYPE(t_enpara), INTENT(IN) :: enpara
      TYPE(t_radfun), INTENT(IN) :: radfun(:)
      TYPE(t_potden), INTENT(IN) :: vTot
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: jspin

      TYPE(t_lapw) :: lapw2
      TYPE(corestate_t), ALLOCATABLE :: core(:,:)
      INTEGER, ALLOCATABLE :: ncore(:)
      REAL, ALLOCATABLE :: fnew(:)
      REAL :: ncore_el, pos(3), a, c2, c3, rj, dev3, c23, jump, sraw, lraw, small, en
      INTEGER :: n, na, l, lo, j, ik, ig, i, nj, iu, jri, ic, nc
      LOGICAL :: l_err, l_exist
      CHARACTER(LEN=8) :: date
      CHARACTER(LEN=10) :: time
      CHARACTER(LEN=4) :: kind
      CHARACTER(LEN=40) :: what
      CHARACTER(LEN=300) :: generated
      CHARACTER(:), ALLOCATABLE :: gitdesc, githash, gitbranch, cdate, cuser, chost, cflags, lflags

      ALLOCATE (atype(atoms%nat))
      DO n = 1, atoms%ntype
         DO na = atoms%firstAtom(n), atoms%firstAtom(n) + atoms%neq(n) - 1
            atype(na) = n
         END DO
      END DO

      title_in = banddos%lapwfn_title
      etot_in = banddos%lapwfn_etot

      CALL readPrevEFermi(efermi, l_err)
      IF (l_err) CALL lapwfn_error("lapwfn: the stored Fermi energy could not be read (cdn.hdf needed)")

      umeas = lapwfn_open(3)
      WRITE (umeas, "(a)") "# lapwfn measurements: joins, cores, small components, inversion relation, phases"
      WRITE (umeas, "(a,es16.8)") "c_light ", c_light(1.0)
      WRITE (umeas, "(a,es24.15)") "efermi_stored ", efermi

      ! the radial functions of every type (slot j of l as in t_radfun and t_abc), divided by r
      ALLOCATE (nfl(0:atoms%lmaxd, atoms%ntype))
      ALLOCATE (utab(atoms%jmtd, nfmaxd, 0:atoms%lmaxd, atoms%ntype), rjoin(nfmaxd, 0:atoms%lmaxd, atoms%ntype))
      ALLOCATE (fnew(atoms%jmtd))
      nfl = 0
      utab = 0.0
      rjoin = 0.0

      ulapw = lapwfn_open(1)
      upw = lapwfn_open(2)
      CALL DATE_AND_TIME(date, time)
      CALL get_compile_desc(gitdesc, githash, gitbranch, cdate, cuser, chost, cflags, lflags)
      generated = writer_version//" from "//TRIM(version_const)//" ("//TRIM(gitdesc)//"), case "//TRIM(title_in)// &
         ", "//date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:6)

      ! core states first (written after the types, but their count is needed with the types)
      ALLOCATE (core(MAX(1, MAXVAL(atoms%econf(:)%num_core_states)), atoms%ntype), ncore(atoms%ntype))
      ncore_el = 0.0
      DO n = 1, atoms%ntype
         CALL lapwfn_core(n, input, atoms, vTot, jspin, core(:, n), ncore(n))
         DO ic = 1, ncore(n)
            ncore_el = ncore_el + atoms%neq(n)*core(ic, n)%occ
         END DO
      END DO

      ! lapwfn.data: BASIC INFO and GEOMETRY
      WRITE (ulapw, "(a)") "LAPWFN.DATA: LAPW orbitals, format version 1"
      WRITE (ulapw, "(a)") "Title"
      WRITE (ulapw, "(a)") TRIM(title_in)
      WRITE (ulapw, "(a)") ""
      WRITE (ulapw, "(a)") "BASIC INFO"
      WRITE (ulapw, "(a)") "----------"
      WRITE (ulapw, "(a)") "Generated by:"
      WRITE (ulapw, "(a)") TRIM(generated)
      WRITE (ulapw, "(a)") "Method:"
      WRITE (ulapw, "(a)") "DFT"
      WRITE (ulapw, "(a)") "DFT functional:"
      WRITE (ulapw, "(a)") TRIM(func_in)
      WRITE (ulapw, "(a)") "Relativistic treatment:"
      WRITE (ulapw, "(a)") "none"
      WRITE (ulapw, "(a)") "Spin polarized:"
      WRITE (ulapw, "(a)") "F"
      WRITE (ulapw, "(a)") "Units:"
      WRITE (ulapw, "(a)") "hartree, bohr"
      WRITE (ulapw, "(a)") "Total energy (au per primitive cell)"
      WRITE (ulapw, "(es24.15)") etot_in
      WRITE (ulapw, "(a)") "Fermi energy (au)"
      WRITE (ulapw, "(es24.15)") efermi
      WRITE (ulapw, "(a)") "Number of electrons per primitive cell (core included)"
      WRITE (ulapw, "(1x,i0)") NINT(input%zelec + ncore_el)
      WRITE (ulapw, "(a)") "Number of valence electrons per primitive cell (in the WAVE FUNCTION section)"
      WRITE (ulapw, "(1x,i0)") NINT(input%zelec)
      WRITE (ulapw, "(a)") "Interstitial part:"
      WRITE (ulapw, "(a)") "pwfn.data: same k points and bands in the same order; coefficients multiply exp(i(k+G).r) with no prefactor"
      WRITE (ulapw, "(a)") ""
      WRITE (ulapw, "(a)") "GEOMETRY"
      WRITE (ulapw, "(a)") "--------"
      WRITE (ulapw, "(a)") "Number of atoms per primitive cell"
      WRITE (ulapw, "(1x,i0)") atoms%nat
      WRITE (ulapw, "(a)") "Atomic number, type, position (au) of each atom"
      DO na = 1, atoms%nat
         pos = MATMUL(cell%amat, atoms%taual(:, na))
         WRITE (ulapw, "(2i5,3es24.15)") NINT(atoms%zatom(atype(na))), atype(na), pos
      END DO
      WRITE (ulapw, "(a)") "Primitive lattice vectors (au)"
      DO i = 1, 3
         WRITE (ulapw, "(3es24.15)") cell%amat(:, i)
      END DO
      WRITE (ulapw, "(a)") ""
      WRITE (ulapw, "(a)") "ATOM TYPES"
      WRITE (ulapw, "(a)") "----------"
      WRITE (ulapw, "(a)") "Number of atom types"
      WRITE (ulapw, "(1x,i0)") atoms%ntype

      DO n = 1, atoms%ntype
         jri = atoms%jri(n)
         DO l = 0, atoms%lmax(n)
            nfl(l, n) = radfun(n)%n_r(l)
            IF (nfl(l, n) > nfmaxd) CALL lapwfn_error("lapwfn: too many radial functions at one l")
            DO j = 1, nfl(l, n)
               utab(1:jri, j, l, n) = radfun(n)%r(1:jri, 1, j, l, jspin)/atoms%rmsh(1:jri, n)
            END DO
         END DO
         WRITE (ulapw, "(a)") "Type, label, atomic number, sphere radius (au), number of mesh points"
         WRITE (ulapw, "(i4,1x,a,i4,es24.15,i6)") n, TRIM(atoms%speciesName(n)), NINT(atoms%zatom(n)), atoms%rmt(n), jri
         WRITE (ulapw, "(a)") "Mesh points (au)"
         WRITE (ulapw, "(4es24.15)") atoms%rmsh(1:jri, n)
         WRITE (ulapw, "(a)") "Maximum l of the sphere expansion"
         WRITE (ulapw, "(1x,i0)") atoms%lmax(n)
         WRITE (ulapw, "(a)") "Number of radial functions for l = 0 to L"
         WRITE (ulapw, "(*(i3))") nfl(0:atoms%lmax(n), n)
         DO l = 0, atoms%lmax(n)
            DO j = 1, nfl(l, n)
               en = enpara%el0(l, n, jspin)
               IF (j == 1) THEN
                  kind = "u"
               ELSE IF (j == 2) THEN
                  kind = "udot"
               ELSE
                  kind = "lo"
                  DO lo = 1, atoms%nlo(n)
                     IF (atoms%llo(lo, n) == l .AND. atoms%slot_of_lo(lo, n) == j) THEN
                        IF (atoms%ulo_der(lo, n) < 1) en = enpara%ello0(lo, n, jspin)
                     END IF
                  END DO
               END IF
               small = MAXVAL(ABS(radfun(n)%r(1:jri, 2, j, l, jspin)))/ &
                       MAX(MAXVAL(ABS(radfun(n)%r(1:jri, 1, j, l, jspin))), 1.0e-300)
               WRITE (what, "(a,i0,a,i0,a,a)") "type ", n, " l ", l, " ", TRIM(kind)
               CALL innerseries(jri, atoms%rmsh(1:jri, n), utab(1:jri, j, l, n), l, atoms%zatom(n), nj, dev3, c23, iu, &
                                rj, a, c2, c3, jump, fnew(1:jri), sraw, lraw, what)
               rjoin(j, l, n) = rj
               WRITE (umeas, "(a,i3,i3,1x,a4,a,i5,3es12.3,a,i5,4es16.8,es10.2,a,2es14.6,a,es10.2)") "join ", n, l, kind, &
                  " | 3-term", nj, atoms%rmsh(nj, n), dev3, c23, " | used", iu, rj, a, c2, c3, jump, " | raw", sraw, lraw, &
                  " | small", small
               WRITE (ulapw, "(a)") "Radial function: l, j, kind, energy (au)"
               WRITE (ulapw, "(2i4,1x,a,es24.15)") l, j, TRIM(kind), en
               WRITE (ulapw, "(a)") "Inner series: join radius, A, c2, c3"
               WRITE (ulapw, "(4es24.15)") rj, a, c2, c3
               WRITE (ulapw, "(a)") "Values u(r) on the mesh"
               WRITE (ulapw, "(4es24.15)") fnew(1:jri)
            END DO
         END DO
      END DO

      WRITE (ulapw, "(a)") ""
      WRITE (ulapw, "(a)") "CORE STATES"
      WRITE (ulapw, "(a)") "-----------"
      WRITE (ulapw, "(a)") "Number of core states of each type"
      WRITE (ulapw, "(*(i4))") ncore(1:atoms%ntype)
      DO n = 1, atoms%ntype
         DO ic = 1, ncore(n)
            nc = core(ic, n)%nc
            WRITE (ulapw, "(a)") "Core state: type, label, n, l, energy (au), occupation, number of mesh points"
            IF (ABS(core(ic, n)%occ - NINT(core(ic, n)%occ)) > 1.0e-6) then
              CALL lapwfn_error("lapwfn: a core state with a fractional occupation")
            end if
            WRITE (ulapw, "(i4,1x,a,2i4,es24.15,i4,i6)") n, TRIM(atoms%speciesName(n)), core(ic, n)%n, core(ic, n)%l, &
               core(ic, n)%e, NINT(core(ic, n)%occ), nc
            WRITE (ulapw, "(a)") "Inner series: join radius, A, c2, c3"
            WRITE (ulapw, "(4es24.15)") core(ic, n)%rj, core(ic, n)%a, core(ic, n)%c2, core(ic, n)%c3
            WRITE (ulapw, "(a)") "Tail beyond the sphere: kappa"
            WRITE (ulapw, "(es24.15)") core(ic, n)%kappa
            WRITE (ulapw, "(a)") "Mesh points (au)"
            WRITE (ulapw, "(4es24.15)") core(ic, n)%r(1:nc)
            WRITE (ulapw, "(a)") "Values u(r) on the mesh"
            WRITE (ulapw, "(4es24.15)") core(ic, n)%u(1:nc)
            WRITE (umeas, "(a,i3,2i3,es20.12,f8.3,a,es14.6,a,es12.4,a,es12.4,a,es12.4,a,es14.8,a,es10.2,a,es12.4,a,i6)") &
               "core ", n, core(ic, n)%n, core(ic, n)%l, core(ic, n)%e, core(ic, n)%occ, " kappa ", core(ic, n)%kappa, &
               " leak_tail ", core(ic, n)%leak_tail, " leak_fleur ", core(ic, n)%leak_fleur, " ratio ", &
               core(ic, n)%leak_tail/MAX(core(ic, n)%leak_fleur, 1.0e-300), " norm ", core(ic, n)%norm, " small ", &
               core(ic, n)%small, " join ", core(ic, n)%rj, " points ", nc
         END DO
      END DO
      WRITE (ulapw, "(a)") ""
      WRITE (ulapw, "(a)") "WAVE FUNCTION"
      WRITE (ulapw, "(a)") "-------------"
      WRITE (ulapw, "(a)") "Number of k-points"
      WRITE (ulapw, "(1x,i0)") kpts%nkpt

      ! the union of the G vectors over the k points
      gmin = 0
      gmax = 0
      DO ik = 1, kpts%nkpt
         CALL lapw2%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, fmpi)
         DO i = 1, 3
            gmin(i) = MIN(gmin(i), MINVAL(lapw2%gvec(i, 1:lapw2%nv(1), 1)))
            gmax(i) = MAX(gmax(i), MAXVAL(lapw2%gvec(i, 1:lapw2%nv(1), 1)))
         END DO
      END DO
      ALLOCATE (gindex(gmin(1):gmax(1), gmin(2):gmax(2), gmin(3):gmax(3)))
      gindex = 0
      ngtot = 0
      DO ik = 1, kpts%nkpt
         CALL lapw2%init(input, noco, nococonv, kpts, atoms, sym, ik, cell, fmpi)
         DO ig = 1, lapw2%nv(1)
            IF (gindex(lapw2%gvec(1, ig, 1), lapw2%gvec(2, ig, 1), lapw2%gvec(3, ig, 1)) == 0) THEN
               ngtot = ngtot + 1
               gindex(lapw2%gvec(1, ig, 1), lapw2%gvec(2, ig, 1), lapw2%gvec(3, ig, 1)) = ngtot
            END IF
         END DO
      END DO
      ALLOCATE (gall(3, ngtot))
      CALL lapwfn_sort_g(cell)

      ! pwfn.data: the header and the G vectors
      WRITE (upw, "(1x,a)") TRIM(title_in)//" (interstitial plane-wave part; spheres and core in lapwfn.data)"
      WRITE (upw, "(a)") " "
      WRITE (upw, "(a)") " BASIC INFO"
      WRITE (upw, "(a)") " ----------"
      WRITE (upw, "(a)") " Generated by:"
      WRITE (upw, "(1x,a)") TRIM(generated)
      WRITE (upw, "(a)") " Method:"
      WRITE (upw, "(a)") " DFT"
      WRITE (upw, "(a)") " DFT Functional"
      WRITE (upw, "(1x,a)") TRIM(func_in)
      WRITE (upw, "(a)") " Pseudopotential"
      WRITE (upw, "(a)") " none (all-electron LAPW; the core and the spheres are in lapwfn.data)"
      WRITE (upw, "(a)") " Plane wave cutoff (au)"
      WRITE (upw, "(1x,es24.15)") 0.5*input%rkmax**2
      WRITE (upw, "(a)") " Spin polarized:"
      WRITE (upw, "(a)") " F"
      WRITE (upw, "(a)") " Total energy (au per primitive cell)"
      WRITE (upw, "(1x,es24.15)") etot_in
      WRITE (upw, "(a)") " Kinetic energy (au per primitive cell)"
      WRITE (upw, "(a)") " 0.0"
      WRITE (upw, "(a)") " Local potential energy (au per primitive cell)"
      WRITE (upw, "(a)") " 0.0"
      WRITE (upw, "(a)") " Non-local potential energy (au per primitive cell)"
      WRITE (upw, "(a)") " 0.0"
      WRITE (upw, "(a)") " Electron-electron energy (au per primitive cell)"
      WRITE (upw, "(a)") " 0.0"
      WRITE (upw, "(a)") " Ion-ion energy (au per primitive cell)"
      WRITE (upw, "(a)") " 0.0"
      WRITE (upw, "(a)") " Number of electrons per primitive cell"
      WRITE (upw, "(1x,i0)") NINT(input%zelec)
      WRITE (upw, "(a)") " "
      WRITE (upw, "(a)") " GEOMETRY"
      WRITE (upw, "(a)") " --------"
      WRITE (upw, "(a)") " Number of atoms per primitive cell"
      WRITE (upw, "(1x,i0)") atoms%nat
      WRITE (upw, "(a)") " Atomic numbers and positions of atoms (au)"
      DO na = 1, atoms%nat
         pos = MATMUL(cell%amat, atoms%taual(:, na))
         WRITE (upw, "(1x,i0,3es24.15)") NINT(atoms%zatom(atype(na))), pos
      END DO
      WRITE (upw, "(a)") " Primitive lattice vectors (au)"
      DO i = 1, 3
         WRITE (upw, "(1x,3es24.15)") cell%amat(:, i)
      END DO
      WRITE (upw, "(a)") " "
      WRITE (upw, "(a)") " G VECTORS"
      WRITE (upw, "(a)") " ---------"
      WRITE (upw, "(a)") " Number of G-vectors"
      WRITE (upw, "(1x,i0)") ngtot
      WRITE (upw, "(a)") " Gx Gy Gz (au)"
      DO ig = 1, ngtot
         WRITE (upw, "(1x,3es24.15)") MATMUL(TRANSPOSE(cell%bmat), REAL(gall(:, ig)))
      END DO
      WRITE (upw, "(a)") " "
      WRITE (upw, "(a)") " WAVE FUNCTION"
      WRITE (upw, "(a)") " -------------"
      WRITE (upw, "(a)") " Number of k-points"
      WRITE (upw, "(1x,i0)") kpts%nkpt

      ! the probe's points
      npts = 0
      INQUIRE (FILE="lapwfn_points.in", EXIST=l_exist)
      IF (l_exist) THEN
         OPEN (NEWUNIT=i, FILE="lapwfn_points.in", STATUS="old", ACTION="read")
         READ (i, *) npts
         ALLOCATE (pts(3, npts))
         DO j = 1, npts
            READ (i, *) pts(:, j)
         END DO
         CLOSE (i)
         upts = lapwfn_open(4)
         WRITE (upts, "(a)") "# ik ib ipt atom r | re im (FLEUR, local frame; atom 0: plane waves) | re im (frame transposed)"// &
            " | re im (plane-wave sum at the point) : orbital values at the points of lapwfn_points.in"
      END IF
      WRITE (*, "(a,i0,a,i0,a)") "lapwfn: writing lapwfn.data and pwfn.data (", kpts%nkpt, " k points, ", ngtot, " G vectors)"
      DEALLOCATE (fnew)
   END SUBROUTINE lapwfn_init

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_sort_g(cell)
      ! the union sorted by |G| then by components, gindex renumbered accordingly
      IMPLICIT NONE
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER :: i1, i2, i3, n, i, j
      INTEGER, ALLOCATABLE :: tmp(:,:), ord(:)
      REAL, ALLOCATABLE :: g2(:)
      REAL :: gc(3)
      LOGICAL :: swap
      ALLOCATE (tmp(3, ngtot), ord(ngtot), g2(ngtot))
      n = 0
      DO i3 = gmin(3), gmax(3)
         DO i2 = gmin(2), gmax(2)
            DO i1 = gmin(1), gmax(1)
               IF (gindex(i1, i2, i3) > 0) THEN
                  n = n + 1
                  tmp(:, n) = [i1, i2, i3]
                  gc = MATMUL(TRANSPOSE(cell%bmat), REAL(tmp(:, n)))
                  g2(n) = DOT_PRODUCT(gc, gc)
                  ord(n) = n
               END IF
            END DO
         END DO
      END DO
      ! insertion sort on (|G|^2 to 1e-10, then components); the sets are a few thousand at most
      DO i = 2, n
         j = i
         DO WHILE (j > 1)
            swap = .FALSE.
            IF (g2(ord(j)) < g2(ord(j-1)) - 1.0e-10) THEN
               swap = .TRUE.
            ELSE IF (ABS(g2(ord(j)) - g2(ord(j-1))) <= 1.0e-10) THEN
               IF (tmp(1, ord(j)) < tmp(1, ord(j-1))) THEN
                  swap = .TRUE.
               ELSE IF (tmp(1, ord(j)) == tmp(1, ord(j-1))) THEN
                  IF (tmp(2, ord(j)) < tmp(2, ord(j-1))) THEN
                     swap = .TRUE.
                  ELSE IF (tmp(2, ord(j)) == tmp(2, ord(j-1))) THEN
                     swap = tmp(3, ord(j)) < tmp(3, ord(j-1))
                  END IF
               END IF
            END IF
            IF (.NOT.swap) EXIT
            i1 = ord(j); ord(j) = ord(j-1); ord(j-1) = i1
            j = j - 1
         END DO
      END DO
      DO i = 1, n
         gall(:, i) = tmp(:, ord(i))
         gindex(gall(1, i), gall(2, i), gall(3, i)) = i
      END DO
      DEALLOCATE (tmp, ord, g2)
   END SUBROUTINE lapwfn_sort_g

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_inversion(n, atoms, sym, lapw, noccbd, ikpt, cof)
      ! Check of the inversion partners (sym%invsat == 2) of type n: their global-frame coefficients, computed
      ! directly, must satisfy A'(l,m) = p (-1)^(l+m) conj(A(l,-m)) with the representative's, p = exp(ik.T) the
      ! Bloch phase of the pair's lattice offset. p is read from the largest LAPW coefficient; the residual is
      ! measured separately on the LAPW (u, udot) and the local-orbital slots, and the writer stops if it fails.
      INTEGER, INTENT(IN) :: n, noccbd, ikpt
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_lapw), INTENT(IN) :: lapw
      COMPLEX, INTENT(IN) :: cof(:, 0:, :, :)
      INTEGER :: ia, na, rep, l, m, lm, lmm, j, ib, lb, mb, jb, ibb
      REAL :: amax, res, reslo, resalt, arg, d
      COMPLEX :: p, pgeom, ref
      CHARACTER(LEN=300) :: str

      DO ia = 1, atoms%neq(n)
         na = atoms%firstAtom(n) + ia - 1
         IF (sym%invsat(na) /= 2) CYCLE
         rep = sym%invsatnr(na)
         IF (atype(rep) /= n) CALL lapwfn_error("lapwfn: an inversion partner of another type")
         amax = 0.0
         ibb = 0; lb = 0; mb = 0; jb = 0
         DO ib = 1, noccbd
            DO l = 0, atoms%lmax(n)
               DO m = -l, l
                  DO j = 1, MIN(2, nfl(l, n))
                     IF (ABS(cof(ib, l*(l + 1) + m, j, na)) > amax) THEN
                        amax = ABS(cof(ib, l*(l + 1) + m, j, na))
                        ibb = ib; lb = l; mb = m; jb = j
                     END IF
                  END DO
               END DO
            END DO
         END DO
         IF (amax <= 0.0) CALL lapwfn_error("lapwfn: an inversion partner has no LAPW coefficients")
         ref = psign(lb + mb)*CONJG(cof(ibb, lb*(lb + 1) - mb, jb, rep))
         IF (ABS(ref) <= 0.0) CALL lapwfn_error("lapwfn: the inversion relation has no reference")
         p = cof(ibb, lb*(lb + 1) + mb, jb, na)/ref
         arg = tpi_const*DOT_PRODUCT(lapw%bkpt, atoms%taual(:, na) + atoms%taual(:, rep))
         pgeom = CMPLX(COS(arg), SIN(arg))
         res = 0.0
         reslo = 0.0
         resalt = 0.0
         DO ib = 1, noccbd
            DO l = 0, atoms%lmax(n)
               DO m = -l, l
                  lm = l*(l + 1) + m
                  lmm = l*(l + 1) - m
                  DO j = 1, nfl(l, n)
                     d = ABS(cof(ib, lm, j, na) - p*psign(l + m)*CONJG(cof(ib, lmm, j, rep)))
                     IF (j <= 2) THEN
                        res = MAX(res, d)
                        resalt = MAX(resalt, ABS(cof(ib, lm, j, na) - p*psign(m)*CONJG(cof(ib, lmm, j, rep))))
                     ELSE
                        reslo = MAX(reslo, d)
                     END IF
                  END DO
               END DO
            END DO
         END DO
         res = res/amax
         reslo = reslo/amax
         resalt = resalt/amax
         WRITE (umeas, "(a,i5,2i4,a,2es17.9,a,2es17.9,a,3es10.2,a,es10.2)") "inversion ", ikpt, na, rep, " p", p, &
            " p_geom", pgeom, " residual lapw, lo, lapw with (-1)^m", res, reslo, resalt, " |p|-1", ABS(p) - 1.0
         IF (ABS(ABS(p) - 1.0) > inv_tol .OR. res > inv_tol .OR. reslo > inv_tol) THEN
            WRITE (str, "(a,i0,a,es10.2,a,es10.2,a)") "lapwfn: the inversion relation fails for atom ", na, &
               " (residual ", res, ", local orbitals ", reslo, "); see lapwfn_measurements.out"
            CALL lapwfn_error(TRIM(str))
         END IF
      END DO
   END SUBROUTINE lapwfn_inversion

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_core(n, input, atoms, vTot, jspin, core, ncore)
      ! The core states of type n solved as cored does (differ on the mesh extended to about ten bohr
      ! beyond the sphere, the potential continued as the input asks), the two kappa solutions of one
      ! (n, l) merged, the tail beyond the sphere replaced by SPEC's exponential and the function
      ! normalised to one over all space.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n, jspin
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_potden), INTENT(IN) :: vTot
      TYPE(corestate_t), INTENT(INOUT) :: core(:)
      INTEGER, INTENT(OUT) :: ncore
      INTEGER :: nstd, nst, korb, jri, ncmsh, i, ierr, l, nc, ic, nj, iu, imax
      INTEGER, ALLOCATABLE :: nprnc(:), kappa(:)
      REAL, ALLOCATABLE :: occ_h(:, :), occ(:), vrd(:), a(:), b(:), r(:), u(:), fnew(:), rext(:), uext(:)
      REAL :: c, z, dxx, rnot, d, rn, t1, t2, rr, fn, fj, fl, weight, eig, rmt, small, kap, dur, ur, nin, tail, norm
      REAL :: dev3, c23, rj, aa, c2, c3, jump, sraw, lraw, leak_f, umax, rend, orig, sgn, r1
      LOGICAL :: found
      CHARACTER(LEN=*), PARAMETER :: spdf = "spdfghik"
      CHARACTER(LEN=40) :: what

      c = c_light(1.0)
      z = atoms%zatom(n)
      dxx = atoms%dx(n)
      rnot = atoms%rmsh(1, n)
      d = EXP(dxx)
      jri = atoms%jri(n)
      rmt = atoms%rmt(n)
      ncmsh = MIN(NINT(LOG((rmt + 10.0)/rnot)/dxx + 1), atoms%msh)
      rn = rnot*d**(ncmsh - 1)
      nstd = MAX(1, atoms%econf(n)%num_states)
      ALLOCATE (nprnc(nstd), kappa(nstd), occ_h(nstd, 2), occ(nstd), vrd(atoms%msh), a(atoms%msh), b(atoms%msh))
      ALLOCATE (r(ncmsh), u(ncmsh), fnew(ncmsh))
      ncore = 0
      IF (atoms%econf(n)%num_core_states == 0) RETURN
      CALL atoms%econf(n)%get_core(nst, nprnc, kappa, occ_h)
      occ(1:nst) = occ_h(1:nst, jspin)
      r(1:jri) = atoms%rmsh(1:jri, n)
      DO i = jri + 1, ncmsh
         r(i) = r(i - 1)*d
      END DO
      vrd = 0.0
      vrd(1:jri) = vTot%mt(1:jri, 0, n, jspin)
      IF (input%l_core_confpot) THEN
         t1 = 0.125
         t1 = MAX((vrd(jri) - vrd(jri - 1)*d)*d/(rmt**2*(d - 1)), t1)
         t2 = vrd(jri)/rmt - rmt*t1
         rr = rmt
      ELSE
         t2 = vrd(jri)/(jri - ncmsh)
      END IF
      DO i = jri + 1, ncmsh
         IF (input%l_core_confpot) THEN
            rr = d*rr
            vrd(i) = rr*(t2 + rr*t1)
         ELSE
            vrd(i) = vrd(jri) + t2*(i - jri)
         END IF
      END DO
      nst = atoms%econf(n)%num_core_states
      DO korb = 1, nst
         IF (occ(korb) == 0.0) CYCLE
         fn = nprnc(korb)
         fj = ABS(kappa(korb)) - 0.5
         weight = 2*occ(korb)
         fl = fj + 0.5*SIGN(1, kappa(korb))
         eig = -2*(z/(fn + fl))**2
         CALL differ(fn, fl, fj, c, z, dxx, rnot, rn, d, ncmsh, vrd, eig, a, b, ierr)
         IF (ierr /= 0) CALL lapwfn_error("lapwfn: the core solver (differ) failed")
         l = NINT(fl)
         small = MAXVAL(ABS(b(1:ncmsh)))/MAX(MAXVAL(ABS(a(1:ncmsh))), 1.0e-300)
         leak_f = lnrint(ncmsh - jri + 1, r(jri:ncmsh), a(jri:ncmsh)**2 + b(jri:ncmsh)**2)
         found = .FALSE.
         DO ic = 1, ncore
            IF (core(ic)%n == nprnc(korb) .AND. core(ic)%l == l) THEN
               found = .TRUE.
               WRITE (umeas, "(a,i3,2i3,a,i3,a,es10.2,a,es10.2)") "core-merge ", n, nprnc(korb), l, " kappa", kappa(korb), &
                  " energy difference", ABS(core(ic)%e - eig), " table difference", &
                  MAXVAL(ABS(a(1:jri)/r(1:jri) - core(ic)%u_raw(1:jri)))/MAXVAL(ABS(core(ic)%u_raw(1:jri)))
               IF (ABS(core(ic)%e - eig) > 1.0e-6) then
                 CALL lapwfn_error("lapwfn: the two kappa solutions of one core (n, l) differ by more than 1e-6 Ha")
               end if
               core(ic)%e = (core(ic)%occ*core(ic)%e + weight*eig)/(core(ic)%occ + weight)
               core(ic)%occ = core(ic)%occ + weight
               core(ic)%small = MAX(core(ic)%small, small)
            END IF
         END DO
         IF (found) CYCLE
         ncore = ncore + 1
         ic = ncore
         core(ic)%n = nprnc(korb)
         core(ic)%l = l
         core(ic)%e = eig
         core(ic)%occ = weight
         core(ic)%small = small
         core(ic)%leak_fleur = leak_f
         ALLOCATE (core(ic)%u_raw(ncmsh))
         u(1:ncmsh) = a(1:ncmsh)/r(1:ncmsh)
         imax = MAXLOC(ABS(u(1:MAX(1, jri/10))), 1)
         sgn = 1.0
         IF (u(imax) < 0.0) sgn = -1.0
         u(1:ncmsh) = sgn*u(1:ncmsh)
         core(ic)%u_raw(1:ncmsh) = u(1:ncmsh)
         ur = u(jri)
         dur = tabderiv(ncmsh, r, u, jri, l, 1)
         IF (ur == 0.0) CALL lapwfn_error("lapwfn: a core state vanishes at the sphere")
         kap = -dur/ur
         IF (kap <= 0.0) CALL lapwfn_error("lapwfn: a core state grows at the sphere")
         umax = MAXVAL(ABS(u(1:jri)))
         ! the tail is tabulated to where it is below tail_floor of the maximum; a shallow state needs
         ! the mesh extended beyond differ's, which is harmless, the tail being analytic there
         rend = rmt + LOG(MAX(ABS(ur)/(tail_floor*umax), 1.0))/kap
         nc = jri + MAX(1, CEILING(LOG(rend/rmt)/dxx))
         IF (nc > SIZE(r)) THEN
            ALLOCATE (rext(nc), uext(nc))
            rext(1:jri) = r(1:jri)
            DO i = jri + 1, nc
               rext(i) = rext(i - 1)*d
            END DO
            uext(1:jri) = u(1:jri)
            CALL MOVE_ALLOC(rext, r)
            CALL MOVE_ALLOC(uext, u)
            DEALLOCATE (fnew)
            ALLOCATE (fnew(nc))
         END IF
         DO i = jri + 1, nc
            u(i) = ur*EXP(-kap*(r(i) - rmt))
         END DO
         ALLOCATE (core(ic)%r(nc), core(ic)%u(nc))
         nin = lnrint(jri, r(1:jri), u(1:jri)**2*r(1:jri)**2)
         tail = ur**2*(rmt**2/(2*kap) + 2*rmt/(2*kap)**2 + 2/(2*kap)**3)
         norm = SQRT(nin + tail)
         u(1:nc) = u(1:nc)/norm
         ur = ur/norm
         core(ic)%leak_tail = tail/(nin + tail)
         core(ic)%kappa = kap
         core(ic)%nc = nc
         core(ic)%r(1:nc) = r(1:nc)
         WRITE (what, "(a,i0,a,i0,a)") "type ", n, " core ", nprnc(korb), spdf(l + 1:l + 1)
         CALL innerseries(nc, r(1:nc), u(1:nc), l, z, nj, dev3, c23, iu, rj, aa, c2, c3, jump, fnew(1:nc), sraw, lraw, what)
         ! the norm of what is written, the series from the origin to the first point (its leading terms,
         ! analytically), the table inside the sphere and the analytic tail, made one; then measured
         ! again with the tabulated tail, so that the measurement is not the construction
         r1 = r(1)
         orig = aa**2*r1**(2*l + 3)*(1.0/(2*l + 3) - 2*z/(l + 1)*r1/(2*l + 4) + (z**2/(l + 1)**2 + 2*c2)*r1**2/(2*l + 5))
         norm = SQRT(orig + lnrint(jri, r(1:jri), fnew(1:jri)**2*r(1:jri)**2) &
                     + ur**2*(rmt**2/(2*kap) + 2*rmt/(2*kap)**2 + 2/(2*kap)**3))
         fnew(1:nc) = fnew(1:nc)/norm
         aa = aa/norm
         core(ic)%norm = orig/norm**2 + lnrint(nc, r(1:nc), fnew(1:nc)**2*r(1:nc)**2)
         core(ic)%u(1:nc) = fnew(1:nc)
         core(ic)%rj = rj
         core(ic)%a = aa
         core(ic)%c2 = c2
         core(ic)%c3 = c3
         WRITE (umeas, "(a,i3,i3,1x,a4,a,i5,3es12.3,a,i5,4es16.8,es10.2,a,2es14.6,a,es10.2)") "join ", n, l, &
            "core", " | 3-term", nj, r(nj), dev3, c23, " | used", iu, rj, aa, c2, c3, jump, " | raw", sraw, lraw, &
            " | small", small
      END DO
      DEALLOCATE (nprnc, kappa, occ_h, occ, vrd, a, b, r, u, fnew)
   END SUBROUTINE lapwfn_core

   !-----------------------------------------------------------------------------
   INTEGER FUNCTION psign(k)
      ! (-1)^k for any integer k
      INTEGER, INTENT(IN) :: k
      psign = 1 - 2*MOD(ABS(k), 2)
   END FUNCTION psign

   REAL FUNCTION polyval6(md, xa, ya, x0)
      ! the md-th derivative (0 to 3) at x0 of the polynomial of degree 5 through six points, from the
      ! Vandermonde system in (x - x0)/h, h the first spacing, solved by dgesv
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: md
      REAL, INTENT(IN) :: xa(:), ya(:), x0
      INTEGER :: k, q, ipiv(6), info
      REAL :: mm(6, 6), bb(6), t, h, fac
      h = xa(2) - xa(1)
      DO k = 1, 6
         t = (xa(k) - x0)/h
         mm(k, 1) = 1.0
         DO q = 2, 6
            mm(k, q) = mm(k, q - 1)*t
         END DO
         bb(k) = ya(k)
      END DO
      CALL dgesv(6, 1, mm, 6, ipiv, bb, 6, info)
      IF (info /= 0) CALL lapwfn_error("lapwfn: a six-point polynomial fit is singular")
      fac = 1.0
      DO q = 2, md
         fac = fac*q
      END DO
      polyval6 = fac*bb(md + 1)/h**md
   END FUNCTION polyval6

   REAL FUNCTION tabderiv(np, r, f, i, l, md)
      ! the md-th derivative (0 to 3) at r(i) of a radial function behaving as r^l: g = f/r^l is fitted by
      ! a polynomial of degree 5 in ln r through the six mesh points about i, and the derivatives of
      ! f = g r^l follow from those of g
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np, i, l, md
      REAL, INTENT(IN) :: r(:), f(:)
      INTEGER :: i0, k
      REAL :: xa(6), ga(6), x0, r0, rl, g0, gx, gxx, gxxx, gr, grr, grrr
      i0 = MIN(MAX(i - 2, 1), np - 5)
      DO k = 1, 6
         xa(k) = LOG(r(i0 + k - 1))
         ga(k) = f(i0 + k - 1)/r(i0 + k - 1)**l
      END DO
      r0 = r(i)
      x0 = LOG(r0)
      rl = r0**l
      g0 = polyval6(0, xa, ga, x0)
      IF (md == 0) THEN
         tabderiv = g0*rl
         RETURN
      END IF
      gx = polyval6(1, xa, ga, x0)
      gr = gx/r0
      IF (md == 1) THEN
         tabderiv = gr*rl + l*g0*rl/r0
         RETURN
      END IF
      gxx = polyval6(2, xa, ga, x0)
      grr = (gxx - gx)/r0**2
      IF (md == 2) THEN
         tabderiv = grr*rl + 2*l*gr*rl/r0 + l*(l - 1)*g0*rl/r0**2
         RETURN
      END IF
      gxxx = polyval6(3, xa, ga, x0)
      grrr = (gxxx - 3*gxx + 2*gx)/r0**3
      tabderiv = grrr*rl + 3*l*grr*rl/r0 + 3*l*(l - 1)*gr*rl/r0**2 + l*(l - 1)*(l - 2)*g0*rl/r0**3
   END FUNCTION tabderiv

   REAL FUNCTION tabval(np, r, f, r0, l)
      ! the value at r0 of a tabulated radial function behaving as r^l, by the same six-point fit in ln r
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np, l
      REAL, INTENT(IN) :: r(:), f(:), r0
      INTEGER :: i, i0, k
      REAL :: xa(6), ga(6)
      i = 1
      DO WHILE (i < np .AND. r(i + 1) <= r0)
         i = i + 1
      END DO
      i0 = MIN(MAX(i - 2, 1), np - 5)
      DO k = 1, 6
         xa(k) = LOG(r(i0 + k - 1))
         ga(k) = f(i0 + k - 1)/r(i0 + k - 1)**l
      END DO
      tabval = polyval6(0, xa, ga, LOG(r0))*r0**l
   END FUNCTION tabval

   REAL FUNCTION lnrint(np, r, f)
      ! the integral of f dr on the logarithmic mesh: Simpson's rule in ln r of f r, with a trapezoid for
      ! the last interval when the number of points is even
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np
      REAL, INTENT(IN) :: r(:), f(:)
      REAL :: h, s
      INTEGER :: k, nq
      h = LOG(r(2)/r(1))
      nq = np
      IF (MOD(np, 2) == 0) nq = np - 1
      s = f(1)*r(1) + f(nq)*r(nq)
      DO k = 2, nq - 1
         IF (MOD(k, 2) == 0) THEN
            s = s + 4*f(k)*r(k)
         ELSE
            s = s + 2*f(k)*r(k)
         END IF
      END DO
      s = s*h/3
      IF (nq /= np) s = s + 0.5*h*(f(np - 1)*r(np - 1) + f(np)*r(np))
      lnrint = s
   END FUNCTION lnrint

   REAL FUNCTION sder(l, zn, a, c2, c3, r, md)
      ! the md-th derivative of the series a r^l (1 - zn r/(l+1) + c2 r^2 + c3 r^3)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l, md
      REAL, INTENT(IN) :: zn, a, c2, c3, r
      REAL :: b(0:3), s, ff
      INTEGER :: k, q, p
      b(0) = 1.0
      b(1) = -zn/(l + 1)
      b(2) = c2
      b(3) = c3
      s = 0.0
      DO k = 0, 3
         p = l + k - md
         IF (p < 0) CYCLE
         ff = 1.0
         DO q = l + k - md + 1, l + k
            ff = ff*q
         END DO
         s = s + b(k)*ff*r**p
      END DO
      sder = a*s
   END FUNCTION sder

   SUBROUTINE fit3(l, zn, rt, ut, dut, a, c2)
      ! the three-term series matched in value and slope at rt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l
      REAL, INTENT(IN) :: zn, rt, ut, dut
      REAL, INTENT(OUT) :: a, c2
      REAL :: g, gp, rl
      rl = rt**l
      g = ut/rl
      gp = (dut - l*ut/rt)/rl
      a = (g - 0.5*gp*rt)/(1.0 - zn*rt/(2*(l + 1)))
      c2 = (gp/a + zn/(l + 1))/(2*rt)
   END SUBROUTINE fit3

   SUBROUTINE fit4(l, zn, rt, ut, dut, d2ut, a, c2, c3)
      ! the four-term series matched in value, slope and curvature at rt
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: l
      REAL, INTENT(IN) :: zn, rt, ut, dut, d2ut
      REAL, INTENT(OUT) :: a, c2, c3
      REAL :: g, gp, gpp, rl, b1, mm(3, 3), bb(3)
      INTEGER :: ipiv(3), info
      rl = rt**l
      g = ut/rl
      gp = (dut - l*ut/rt)/rl
      gpp = (d2ut - 2*l*dut/rt + l*(l + 1)*ut/rt**2)/rl
      b1 = -zn/(l + 1)
      mm(1, 1) = 1.0 + b1*rt; mm(1, 2) = rt**2; mm(1, 3) = rt**3
      mm(2, 1) = b1; mm(2, 2) = 2*rt; mm(2, 3) = 3*rt**2
      mm(3, 1) = 0.0; mm(3, 2) = 2.0; mm(3, 3) = 6*rt
      bb(1) = g; bb(2) = gp; bb(3) = gpp
      CALL dgesv(3, 1, mm, 3, ipiv, bb, 3, info)
      IF (info /= 0) CALL lapwfn_error("lapwfn: the four-term series could not be fitted")
      a = bb(1)
      c2 = bb(2)/a
      c3 = bb(3)/a
   END SUBROUTINE fit4

   SUBROUTINE innerseries(np, r, f, l, zn, nj, dev3, c23, iu, rj, a, c2, c3, jump, fnew, sraw, lraw, what)
      ! The regular series inside the join radius, by the rule of the reference converter: the three-term
      ! join (the smallest mesh radius at which the series matched in value and slope reproduces the table to
      ! 1e-5 of its largest value over [r, 2r], or the radius of least deviation if none does and the least
      ! is below 1e-3); the join used at one and a half times it and at least 1e-4 bohr, moved outward for s
      ! and p until the three- and four-term fits agree on c2 to ten per cent (bounded by 0.03 bohr); the
      ! four-term series (s, p) or the three-term series (l >= 2) replaces the table inside. sraw and lraw
      ! record the raw table's log-log slope between its first two points and the logarithmic derivative at
      ! its second.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: np, l
      REAL, INTENT(IN) :: r(:), f(:), zn
      INTEGER, INTENT(OUT) :: nj, iu
      REAL, INTENT(OUT) :: dev3, c23, rj, a, c2, c3, jump, fnew(:), sraw, lraw
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER :: k, k2, kb, i0
      REAL, PARAMETER :: tol3 = 1.0e-5
      REAL :: rt, ut, dut, d2ut, a3, c2_3, a4, c2_4, c3_4, dev, devb, devj, fm, rbound, d3t, d3s
      LOGICAL :: l_settled
      CHARACTER(LEN=300) :: str
      sraw = LOG(ABS(f(2)/f(1)))/LOG(r(2)/r(1))
      lraw = tabderiv(np, r, f, 2, l, 1)/f(2)
      nj = 0
      kb = 0
      devb = HUGE(1.0)
      devj = HUGE(1.0)
      DO k = 3, np - 1
         rt = r(k)
         k2 = k
         DO WHILE (k2 < np .AND. r(MIN(k2 + 1, np)) <= 2*rt)
            k2 = k2 + 1
         END DO
         IF (k2 - k + 1 < 4) EXIT
         ut = tabderiv(np, r, f, k, l, 0)
         dut = tabderiv(np, r, f, k, l, 1)
         CALL fit3(l, zn, rt, ut, dut, a3, c2_3)
         fm = MAX(MAXVAL(ABS(f(k:k2))), 1.0e-300)
         dev = 0.0
         DO i0 = k, k2
            dev = MAX(dev, ABS(sder(l, zn, a3, c2_3, 0.0, r(i0), 0) - f(i0)))
         END DO
         dev = dev/fm
         IF (dev < devb) THEN
            devb = dev
            kb = k
         END IF
         IF (dev < tol3) THEN
            nj = k
            devj = dev
            EXIT
         END IF
      END DO
      IF (nj == 0) THEN
         IF (devb > 1.0e-3) THEN
            WRITE (str, "(a,a,a,es8.1,a)") "lapwfn: no join radius reproduces the radial function ", TRIM(what), &
               " to 1e-3 (best ", devb, ")"
            CALL lapwfn_error(TRIM(str))
         END IF
         WRITE (*, "(a,a,a,es8.1,a,es8.1,a)") "lapwfn: ", TRIM(what), ": three-term join radius of least deviation taken, ", &
            devb, " (over ", tol3, ")"
         WRITE (umeas, "(a,a,a,es8.1)") "join-fallback ", TRIM(what), " least deviation ", devb
         nj = kb
         devj = devb
      END IF
      dev3 = devj
      ut = tabderiv(np, r, f, nj, l, 0)
      dut = tabderiv(np, r, f, nj, l, 1)
      CALL fit3(l, zn, r(nj), ut, dut, a3, c23)
      rbound = MIN(0.05, 0.1*r(np))
      i0 = nj
      DO WHILE (i0 < np .AND. (r(i0) < 1.5*r(nj) .OR. r(i0) < 1.0e-4))
         i0 = i0 + 1
      END DO
      IF (r(i0) > rbound) THEN
         WRITE (str, "(a,es10.3,a,a,a)") "lapwfn: the join radius ", r(i0), " of the radial function ", TRIM(what), &
            " would exceed the bound of the format (0.05 bohr and a tenth of the sphere radius)"
         CALL lapwfn_error(TRIM(str))
      END IF
      iu = i0
      IF (l <= 1) THEN
         l_settled = .FALSE.
         DO k = i0, np - 1
            IF (r(k) >= MIN(0.03, rbound)) EXIT
            ut = tabderiv(np, r, f, k, l, 0)
            dut = tabderiv(np, r, f, k, l, 1)
            d2ut = tabderiv(np, r, f, k, l, 2)
            CALL fit3(l, zn, r(k), ut, dut, a3, c2_3)
            CALL fit4(l, zn, r(k), ut, dut, d2ut, a4, c2_4, c3_4)
            IF (ABS(c2_4/c2_3 - 1.0) < 0.1) THEN
               iu = k
               l_settled = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT.l_settled) WRITE (umeas, "(a,a,a,es10.3,a)") "Info: the four-term series of ", TRIM(what), &
            " did not settle (c2 of the four- and three-term fits within 10 per cent) below r = ", MIN(0.03, rbound), &
            "; the three-term radius is used"
      END IF
      rj = r(iu)
      ut = tabderiv(np, r, f, iu, l, 0)
      dut = tabderiv(np, r, f, iu, l, 1)
      d2ut = tabderiv(np, r, f, iu, l, 2)
      fnew(1:np) = f(1:np)
      IF (l <= 1) THEN
         CALL fit4(l, zn, rj, ut, dut, d2ut, a, c2, c3)
         d3t = tabderiv(np, r, f, iu, l, 3)
         d3s = sder(l, zn, a, c2, c3, rj, 3)
         jump = ABS(d3s - d3t)/MAX(ABS(d3t), 1.0e-300)
      ELSE
         CALL fit3(l, zn, rj, ut, dut, a, c2)
         c3 = 0.0
         d3s = sder(l, zn, a, c2, c3, rj, 2)
         jump = ABS(d3s - d2ut)/MAX(ABS(d2ut), 1.0e-300)
         IF (jump > 2.0e-3) THEN
            WRITE (str, "(a,es8.1,a,a,a)") "lapwfn: the curvature jump ", jump, " at the join of the radial function ", &
               TRIM(what), " exceeds the 2e-3 of the format"
            CALL lapwfn_error(TRIM(str))
         END IF
      END IF
      DO k = 1, iu - 1
         fnew(k) = sder(l, zn, a, c2, c3, r(k), 0)
      END DO
   END SUBROUTINE innerseries

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_probe(ikpt_i, atoms, sym, cell, lapw, nb, nw, cofl, cpw, tmat)
      ! FLEUR's own orbitals at the points of lapwfn_points.in, combined as the written bands are (tmat gives
      ! each written band as a combination of FLEUR's nb bands: the phases, and the real combinations of a
      ! degenerate level): inside a sphere from the local-frame coefficients (calc_abc with the real symmetry),
      ! the harmonics taken at R x, R the Cartesian form of the atom's operation (minus it for an inversion
      ! partner), with the Bloch phase of the image; beside it the same with R transposed, and the plane-wave
      ! sum at the point.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ikpt_i, nb, nw
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_sym), INTENT(IN) :: sym
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_lapw), INTENT(IN) :: lapw
      COMPLEX, INTENT(IN) :: cofl(:, 0:, :, :), cpw(:, :), tmat(:, :)
      INTEGER :: ip, na, n, l, m, lm, j, ib, ig, i1, i2, i3, op, jas
      REAL :: x(3), pos(3), dd(3), f0(3), tt(3), pp(3), dmin, tv(3), p(3), r, dir(3), kc(3), arg, ainv(3, 3)
      REAL :: rg(3, 3, sym%nop), rr(3, 3), ul(nfmaxd, 0:atoms%lmaxd), gc(3)
      COMPLEX :: ylm((atoms%lmaxd + 1)**2), ylmt((atoms%lmaxd + 1)**2), val, valt, valpw, bloch, vb(nb), vbt(nb)
      COMPLEX, ALLOCATABLE :: eikg(:)
      DO op = 1, sym%nop
         rg(:, :, op) = MATMUL(cell%amat, MATMUL(REAL(sym%mrot(:, :, op)), cell%bmat))/tpi_const
      END DO
      ainv = cell%bmat/tpi_const
      kc = MATMUL(TRANSPOSE(cell%bmat), lapw%bkpt)
      ALLOCATE (eikg(ngtot))
      DO ip = 1, npts
         x = pts(:, ip)
         DO ig = 1, ngtot
            gc = MATMUL(TRANSPOSE(cell%bmat), REAL(gall(:, ig)))
            arg = DOT_PRODUCT(kc + gc, x)
            eikg(ig) = CMPLX(COS(arg), SIN(arg))
         END DO
         jas = 0
         tv = 0.0
         p = 0.0
         dmin = HUGE(1.0)
         DO na = 1, atoms%nat
            pos = MATMUL(cell%amat, atoms%taual(:, na))
            dd = x - pos
            f0 = ANINT(MATMUL(ainv, dd))
            dmin = HUGE(1.0)
            DO i3 = -1, 1
               DO i2 = -1, 1
                  DO i1 = -1, 1
                     tt = MATMUL(cell%amat, f0 + [REAL(i1), REAL(i2), REAL(i3)])
                     pp = dd - tt
                     IF (SUM(pp**2) < dmin) THEN
                        dmin = SUM(pp**2)
                        tv = tt
                        p = pp
                     END IF
                  END DO
               END DO
            END DO
            IF (SQRT(dmin) < atoms%rmt(atype(na))) THEN
               jas = na
               EXIT
            END IF
         END DO
         IF (jas == 0) THEN
            DO ib = 1, nw
               valpw = SUM(cpw(:, ib)*eikg)
               WRITE (upts, "(3i6,i4,es16.8,6es24.15)") ikpt_i, ib, ip, 0, 0.0, valpw, valpw, valpw
            END DO
            CYCLE
         END IF
         na = jas
         n = atype(na)
         r = SQRT(dmin)
         IF (r < atoms%rmsh(1, n)) THEN     ! at the nucleus: the value at the first mesh point, along z
            dir = [0.0, 0.0, 1.0]
            IF (r > 0.0) dir = p/r
            r = atoms%rmsh(1, n)
         ELSE
            dir = p/r
         END IF
         arg = DOT_PRODUCT(kc, tv)
         bloch = CMPLX(COS(arg), SIN(arg))
         IF (sym%invsat(na) == 2) THEN
            rr = -rg(:, :, sym%ngopr(sym%invsatnr(na)))
         ELSE
            rr = rg(:, :, sym%ngopr(na))
         END IF
         CALL ylm4(atoms%lmax(n), MATMUL(rr, dir), ylm)
         CALL ylm4(atoms%lmax(n), MATMUL(TRANSPOSE(rr), dir), ylmt)
         ul = 0.0
         DO l = 0, atoms%lmax(n)
            DO j = 1, nfl(l, n)
               ul(j, l) = tabval(atoms%jri(n), atoms%rmsh(1:atoms%jri(n), n), utab(1:atoms%jri(n), j, l, n), r, l)
            END DO
         END DO
         DO ib = 1, nb
            val = CMPLX(0.0, 0.0)
            valt = CMPLX(0.0, 0.0)
            DO l = 0, atoms%lmax(n)
               DO m = -l, l
                  lm = l*(l + 1) + m
                  DO j = 1, nfl(l, n)
                     val = val + cofl(ib, lm, j, na)*ul(j, l)*ylm(lm + 1)
                     valt = valt + cofl(ib, lm, j, na)*ul(j, l)*ylmt(lm + 1)
                  END DO
               END DO
            END DO
            vb(ib) = bloch*val
            vbt(ib) = bloch*valt
         END DO
         DO ib = 1, nw
            val = SUM(tmat(1:nb, ib)*vb)
            valt = SUM(tmat(1:nb, ib)*vbt)
            valpw = SUM(cpw(:, ib)*eikg)
            WRITE (upts, "(3i6,i4,es16.8,6es24.15)") ikpt_i, ib, ip, na, r, val, valt, valpw
         END DO
      END DO
      DEALLOCATE (eikg)
   END SUBROUTINE lapwfn_probe

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_ylmtest()
      ! FLEUR's harmonics at fifty directions, for a reader to compare with its own; written by lapwfn_finish,
      ! when every k point is done, so that a refusal leaves no file behind
      IMPLICIT NONE
      INTEGER, PARAMETER :: lmx = 10, ndir = 50
      INTEGER :: id, l, m, u, k
      INTEGER(int64) :: s
      REAL :: v(3)
      COMPLEX :: ylm((lmx + 1)**2)
      OPEN (NEWUNIT=u, FILE="lapwfn_ylm.out"//part, STATUS="replace", ACTION="write")
      WRITE (u, "(a)") '# l m vx vy vz re im : FLEUR''s ylm4 (Y_lm) at fifty directions from an LCG, l <= 10'
      s = 12345
      DO id = 1, ndir
         DO k = 1, 3
            s = MOD(1103515245_int64*s + 12345_int64, 2147483648_int64)
            v(k) = 2.0*REAL(s)/2147483648.0 - 1.0
         END DO
         CALL ylm4(lmx, v, ylm)
         DO l = 0, lmx
            DO m = -l, l
               WRITE (u, "(2i4,3es20.12,2es24.15)") l, m, v, ylm(l*(l + 1) + 1 + m)
            END DO
         END DO
      END DO
      CLOSE (u)
      IF (c_rename("lapwfn_ylm.out"//part//C_NULL_CHAR, "lapwfn_ylm.out"//C_NULL_CHAR) /= 0) then
        CALL lapwfn_error("lapwfn: could not rename lapwfn_ylm.out"//part)
      end if
   END SUBROUTINE lapwfn_ylmtest

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_functional(xcpot, func)
      ! the functional's name for the file headers
      CLASS(t_xcpot), INTENT(IN) :: xcpot
      CHARACTER(LEN=*), INTENT(OUT) :: func
      IF (xcpot%l_libxc) THEN
         WRITE (func, "(a,i0,a,i0)") "LibXC exchange id ", xcpot%func_vxc_id_x, ", correlation id ", xcpot%func_vxc_id_c
      ELSE
         func = TRIM(xcpot%inbuild_name)//" (FLEUR built-in)"
      END IF
   END SUBROUTINE lapwfn_functional

   SUBROUTINE lapwfn_finish()
      ! every output file complete: closed and renamed from its .part name, lapwfn_ylm.out written first
      IMPLICIT NONE
      INTEGER :: i
      CALL lapwfn_ylmtest()
      WRITE (umeas, "(a,i0,a)") "done ", nkdone, " k points"
      DO i = 1, nfiles
         CALL lapwfn_close(i, .TRUE.)
      END DO
      WRITE (*, "(a,i0,a)") "lapwfn: lapwfn.data and pwfn.data written (", nkdone, " k points)"
   END SUBROUTINE lapwfn_finish

   !-----------------------------------------------------------------------------
   INTEGER FUNCTION lapwfn_open(i)
      ! output file i opened under its .part name
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      OPEN (NEWUNIT=lapwfn_open, FILE=TRIM(fname(i))//part, STATUS="replace", ACTION="write")
      funit(i) = lapwfn_open
      fopen(i) = .TRUE.
   END FUNCTION lapwfn_open

   SUBROUTINE lapwfn_close(i, l_keep)
      ! output file i closed: renamed from its .part name when complete (l_keep), deleted otherwise
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i
      LOGICAL, INTENT(IN) :: l_keep
      IF (.NOT.fopen(i)) RETURN
      fopen(i) = .FALSE.
      IF (l_keep) THEN
         CLOSE (funit(i))
         IF (c_rename(TRIM(fname(i))//part//C_NULL_CHAR, TRIM(fname(i))//C_NULL_CHAR) /= 0) then
           CALL juDFT_error("lapwfn: could not rename "//TRIM(fname(i))//part//" to "//TRIM(fname(i)), calledby="lapwfn")
         end if
      ELSE
         CLOSE (funit(i), STATUS="delete")
      END IF
   END SUBROUTINE lapwfn_close

   SUBROUTINE lapwfn_error(msg)
      ! a refusal, wherever it comes: this run's output files, all written under .part names, are deleted, so
      ! that the directory is left as it was (a previous good pair untouched), and FLEUR stops with the message
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: msg
      INTEGER :: i
      DO i = 1, nfiles
         CALL lapwfn_close(i, .FALSE.)
      END DO
      CALL juDFT_error(msg, calledby="lapwfn")
   END SUBROUTINE lapwfn_error

   SUBROUTINE lapwfn_stored_spins(nspins, l_err)
      ! the spin count the stored density was written with, read from cdn.hdf as readPrevEFermi reads the
      ! stored Fermi energy; l_err when the density is not stored in HDF5 form
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: nspins
      LOGICAL, INTENT(OUT) :: l_err
      INTEGER :: mode
#ifdef CPP_HDF
      INTEGER(HID_T) :: fileID
      INTEGER :: currentStarsIndex, currentLatharmsIndex, currentStructureIndex, currentStepfunctionIndex
      INTEGER :: readDensityIndex, lastDensityIndex, starsIndex, latharmsIndex, structureIndex, stepfunctionIndex
      INTEGER :: date, time, iter, previousDensityIndex
      REAL :: fermiEnergy, distance
      LOGICAL :: l_qfix
      CHARACTER(LEN=30) :: archiveName
#endif
      nspins = 0
      l_err = .TRUE.
      CALL getIOMode(mode)
      IF (mode /= 3) RETURN               ! CDN_HDF5_MODE of m_cdn_io, which keeps the constant private
#ifdef CPP_HDF
      CALL openCDN_HDF(fileID, currentStarsIndex, currentLatharmsIndex, currentStructureIndex, &
                       currentStepfunctionIndex, readDensityIndex, lastDensityIndex)
      WRITE (archiveName, "(a,i0)") "/cdn-", readDensityIndex
      CALL peekDensityEntryHDF(fileID, archiveName, DENSITY_TYPE_FFN_IN_const, iter, starsIndex, latharmsIndex, &
                               structureIndex, stepfunctionIndex, previousDensityIndex, nspins, date, time, &
                               distance, fermiEnergy, l_qfix)
      CALL closeCDNPOT_HDF(fileID)
      l_err = .FALSE.
#endif
   END SUBROUTINE lapwfn_stored_spins

   !-----------------------------------------------------------------------------
   SUBROUTINE lapwfn_conj(twok, cb, pb, kcb, kpb)
      ! The conjugate orbital's coefficients at k = G/2, the antiunitary map K: conj c_{-G-2k} on the plane
      ! waves (zero where the partner lies outside the union) and (-1)^m conj c_{l,-m} on every atom's sphere
      ! coefficients in the global frame.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: twok(3)
      COMPLEX, INTENT(IN) :: cb(0:, :, :), pb(:)
      COMPLEX, INTENT(OUT) :: kcb(0:, :, :), kpb(:)
      INTEGER :: na, n, l, m, j, ig, idx, g2(3)
      kcb = CMPLX(0.0, 0.0)
      kpb = CMPLX(0.0, 0.0)
      DO na = 1, SIZE(cb, 3)
         n = atype(na)
         DO l = 0, UBOUND(nfl, 1)
            DO m = -l, l
               DO j = 1, nfl(l, n)
                  kcb(l*(l + 1) + m, j, na) = psign(m)*CONJG(cb(l*(l + 1) - m, j, na))
               END DO
            END DO
         END DO
      END DO
      DO ig = 1, ngtot
         g2 = -gall(:, ig) - twok
         idx = 0
         IF (ALL(g2 >= gmin) .AND. ALL(g2 <= gmax)) idx = gindex(g2(1), g2(2), g2(3))
         IF (idx > 0) kpb(ig) = CONJG(pb(idx))
      END DO
   END SUBROUTINE lapwfn_conj

   COMPLEX FUNCTION lapwfn_dot(ca, pa, cb, pb)
      ! the Euclidean inner product of two full coefficient vectors
      IMPLICIT NONE
      COMPLEX, INTENT(IN) :: ca(0:, :, :), pa(:), cb(0:, :, :), pb(:)
      lapwfn_dot = SUM(CONJG(ca)*cb) + SUM(CONJG(pa)*pb)
   END FUNCTION lapwfn_dot

   SUBROUTINE lapwfn_phase(ib, twok, cof, cpw, lok, ph, devp, devm)
      ! A constant phase making band ib real at k = G/2, if one exists: from the band's largest coefficient
      ! and its partner under the relation of a real function (c_{-G-2k} = c_G^*, A_{l,-m} = (-1)^m A_lm^*),
      ! the partner of a band real up to a phase being exp(-2i phi) conj(largest); the phase is applied and
      ! the two tests of a real function made, each part against its own largest coefficient, so that large
      ! sphere coefficients cannot hide an imaginary plane-wave part; a band that fails is left as it was
      ! (ph = 1). devp and devm are the deviations relative to the band's largest coefficient, -1 when the
      ! band has no partner coefficient.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ib, twok(3)
      COMPLEX, INTENT(INOUT) :: cof(:, 0:, :, :), cpw(:, :)
      LOGICAL, INTENT(OUT) :: lok
      COMPLEX, INTENT(OUT) :: ph
      REAL, INTENT(OUT) :: devp, devm
      INTEGER :: na, n, l, m, j, ig, idx, g2(3)
      COMPLEX :: cmax, cpart, q2
      COMPLEX, ALLOCATABLE :: kcb(:, :, :), kpb(:)
      REAL :: amax, devpw, devmt
      lok = .FALSE.
      ph = CMPLX(1.0, 0.0)
      devp = -1.0
      devm = -1.0
      amax = MAX(MAXVAL(ABS(cof(ib, :, :, :))), MAXVAL(ABS(cpw(:, ib))))
      IF (amax <= 0.0) RETURN
      cmax = CMPLX(0.0, 0.0)
      cpart = CMPLX(0.0, 0.0)
      IF (MAXVAL(ABS(cof(ib, :, :, :))) >= MAXVAL(ABS(cpw(:, ib)))) THEN
         DO na = 1, SIZE(cof, 4)
            n = atype(na)
            DO l = 0, UBOUND(nfl, 1)
               DO m = -l, l
                  DO j = 1, nfl(l, n)
                     IF (ABS(cof(ib, l*(l+1)+m, j, na)) > ABS(cmax)) THEN
                        cmax = cof(ib, l*(l+1)+m, j, na)
                        cpart = psign(m)*cof(ib, l*(l+1)-m, j, na)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      ELSE
         DO ig = 1, ngtot
            IF (ABS(cpw(ig, ib)) > ABS(cmax)) THEN
               cmax = cpw(ig, ib)
               g2 = -gall(:, ig) - twok
               cpart = CMPLX(0.0, 0.0)
               IF (ALL(g2 >= gmin) .AND. ALL(g2 <= gmax)) THEN
                  idx = gindex(g2(1), g2(2), g2(3))
                  IF (idx > 0) cpart = cpw(idx, ib)
               END IF
            END IF
         END DO
      END IF
      IF (ABS(cpart) < 1.0e-3*ABS(cmax)) RETURN
      q2 = CONJG(cmax)/cpart            ! exp(2 i phi)
      ph = CMPLX(COS(0.5*ATAN2(AIMAG(q2), REAL(q2))), SIN(0.5*ATAN2(AIMAG(q2), REAL(q2))))
      cof(ib, :, :, :) = ph*cof(ib, :, :, :)
      cpw(:, ib) = ph*cpw(:, ib)
      ALLOCATE (kcb(0:UBOUND(cof, 2), SIZE(cof, 3), SIZE(cof, 4)), kpb(SIZE(cpw, 1)))
      CALL lapwfn_conj(twok, cof(ib, :, :, :), cpw(:, ib), kcb, kpb)
      devpw = MAXVAL(ABS(kpb - cpw(:, ib)))
      devmt = MAXVAL(ABS(kcb - cof(ib, :, :, :)))
      DEALLOCATE (kcb, kpb)
      devp = devpw/amax
      devm = devmt/amax
      lok = (devpw <= real_tol*MAXVAL(ABS(cpw(:, ib)))) .AND. (devmt <= real_tol*MAXVAL(ABS(cof(ib, :, :, :))))
      IF (.NOT.lok) THEN                ! no phase is applied to a band left complex
         cof(ib, :, :, :) = CONJG(ph)*cof(ib, :, :, :)
         cpw(:, ib) = CONJG(ph)*cpw(:, ib)
         ph = CMPLX(1.0, 0.0)
      END IF
   END SUBROUTINE lapwfn_phase

   SUBROUTINE lapwfn_realise(i1, i2, twok, cof, cpw, eig, eigw, tmat, phs, lreal, devp, devm, &
                             closure, dmm, dms, sn, sn1, dcc, dlab, ok, why)
      ! The real orthonormal combinations of the degenerate level i1..i2 at k = G/2 (the recipe of Elk's and
      ! exciting's writers). The conjugate of a member, K v (lapwfn_conj), lies in the level's span when the
      ! level is complete: K V = V M with M by least squares on the full coefficient vectors (exact for a
      ! closed level whatever the inner product), closure its residual relative to the level's largest
      ! coefficient; M is then unitary and symmetric (dmm, dms). The 2n candidates (I + M) e_j and
      ! i(I - M) e_j, coordinates in the level's own basis, are each K-invariant, and in that basis, FLEUR's
      ! eigenvectors being orthonormal, the inner product is the Euclidean one and real between K-invariant
      ! vectors: the real SVD of their real and imaginary parts stacked gives n orthonormal real combinations
      ! C (n singular values of order one, sn the smallest of them, and n at rounding, sn1 the largest of
      ! those), C unitary to dcc, and W = V C real. Each new band is labelled with its expectation energy
      ! sum_j |C_jk|^2 e_j (dlab the largest change from the eigenvalue it replaces) and the level is written
      ! in ascending order of these; the two tests of a real function are then applied to every new band. ok
      ! is false, and the level left as FLEUR gave it, when it does not close, the candidates do not span, or
      ! a new band fails the tests; why says which. tmat receives each written band as a combination of
      ! FLEUR's (the phases lapwfn_phase applied to members included), for the probe.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i1, i2, twok(3)
      COMPLEX, INTENT(INOUT) :: cof(:, 0:, :, :), cpw(:, :), tmat(:, :), phs(:)
      REAL, INTENT(IN) :: eig(:)
      REAL, INTENT(INOUT) :: eigw(:), devp(:), devm(:)
      LOGICAL, INTENT(INOUT) :: lreal(:)
      REAL, INTENT(OUT) :: closure, dmm, dms, sn, sn1, dcc, dlab
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=*), INTENT(OUT) :: why
      INTEGER :: n, j, k, kk, info, lwork, nlmd, nf, nat, ng
      INTEGER, ALLOCATABLE :: ipiv(:)
      COMPLEX, ALLOCATABLE :: kc(:, :, :, :), kp(:, :), g(:, :), b(:, :), m(:, :), x(:, :), c(:, :)
      COMPLEX, ALLOCATABLE :: wc(:, :, :, :), wp(:, :), rc(:, :, :), rp(:), kcb(:, :, :), kpb(:), ccol(:)
      REAL, ALLOCATABLE :: a(:, :), s(:), u(:, :), vt(:, :), work(:), en(:), dp(:), dm(:)
      REAL :: amax, wq(1), t, bmax
      n = i2 - i1 + 1
      nlmd = UBOUND(cof, 2)
      nf = SIZE(cof, 3)
      nat = SIZE(cof, 4)
      ng = SIZE(cpw, 1)
      closure = 0.0; dmm = 0.0; dms = 0.0; sn = 0.0; sn1 = 0.0; dcc = 0.0; dlab = 0.0
      ok = .FALSE.
      why = "recombined"
      ALLOCATE (kc(0:nlmd, nf, nat, n), kp(ng, n), g(n, n), b(n, n), m(n, n), ipiv(n), rc(0:nlmd, nf, nat), rp(ng))
      amax = 0.0
      DO j = 1, n
         CALL lapwfn_conj(twok, cof(i1 - 1 + j, :, :, :), cpw(:, i1 - 1 + j), kc(:, :, :, j), kp(:, j))
         amax = MAX(amax, MAXVAL(ABS(cof(i1 - 1 + j, :, :, :))), MAXVAL(ABS(cpw(:, i1 - 1 + j))))
      END DO
      DO j = 1, n
         DO k = 1, n
            g(k, j) = lapwfn_dot(cof(i1 - 1 + k, :, :, :), cpw(:, i1 - 1 + k), cof(i1 - 1 + j, :, :, :), cpw(:, i1 - 1 + j))
            b(k, j) = lapwfn_dot(cof(i1 - 1 + k, :, :, :), cpw(:, i1 - 1 + k), kc(:, :, :, j), kp(:, j))
         END DO
      END DO
      m = b
      CALL zgesv(n, n, g, n, ipiv, m, n, info)
      IF (info /= 0) THEN
         why = 'the level''s Gram matrix is singular: left as FLEUR gives it'
         RETURN
      END IF
      DO j = 1, n
         rc = kc(:, :, :, j)
         rp = kp(:, j)
         DO k = 1, n
            rc = rc - m(k, j)*cof(i1 - 1 + k, :, :, :)
            rp = rp - m(k, j)*cpw(:, i1 - 1 + k)
         END DO
         closure = MAX(closure, MAX(MAXVAL(ABS(rc)), MAXVAL(ABS(rp)))/amax)
      END DO
      IF (closure > clos_tol) THEN
         why = "does not close under conjugation: left as FLEUR gives it"
         RETURN
      END IF
      ALLOCATE (x(n, 2*n), c(n, n))
      x(:, 1:n) = MATMUL(m, CONJG(TRANSPOSE(m)))
      DO j = 1, n
         x(j, j) = x(j, j) - CMPLX(1.0, 0.0)
      END DO
      dmm = MAXVAL(ABS(x(:, 1:n)))
      dms = MAXVAL(ABS(m - TRANSPOSE(m)))
      ! the candidates (I + M) e_j and i(I - M) e_j in the level's own basis, their real and imaginary parts
      x(:, 1:n) = m
      x(:, n + 1:2*n) = -CMPLX(0.0, 1.0)*m
      DO j = 1, n
         x(j, j) = x(j, j) + CMPLX(1.0, 0.0)
         x(j, n + j) = x(j, n + j) + CMPLX(0.0, 1.0)
      END DO
      ALLOCATE (a(2*n, 2*n), s(2*n), u(2*n, 2*n), vt(1, 1))
      a(1:n, :) = REAL(x)
      a(n + 1:2*n, :) = AIMAG(x)
      CALL dgesvd("A", "N", 2*n, 2*n, a, 2*n, s, u, 2*n, vt, 1, wq, -1, info)
      lwork = MAX(NINT(wq(1)), 10*n + 10)
      ALLOCATE (work(lwork))
      CALL dgesvd("A", "N", 2*n, 2*n, a, 2*n, s, u, 2*n, vt, 1, work, lwork, info)
      IF (info /= 0 .OR. s(1) <= 0.0) THEN
         why = "the singular value decomposition failed: left as FLEUR gives it"
         RETURN
      END IF
      sn = s(n)/s(1)
      sn1 = s(n + 1)/s(1)
      IF (sn < 1.0e-3 .OR. sn1 > 1.0e-6) THEN
         why = "the real candidates do not span the level: left as FLEUR gives it"
         RETURN
      END IF
      DO k = 1, n
         c(:, k) = CMPLX(u(1:n, k), u(n + 1:2*n, k))
      END DO
      x(:, 1:n) = MATMUL(CONJG(TRANSPOSE(c)), c)
      DO j = 1, n
         x(j, j) = x(j, j) - CMPLX(1.0, 0.0)
      END DO
      dcc = MAXVAL(ABS(x(:, 1:n)))
      ! the new bands, their expectation energies, in ascending order of these
      ALLOCATE (en(n), ccol(n))
      DO k = 1, n
         en(k) = SUM(ABS(c(:, k))**2*eig(i1:i2))
      END DO
      DO k = 1, n - 1
         kk = k
         DO j = k + 1, n
            IF (en(j) < en(kk)) kk = j
         END DO
         IF (kk /= k) THEN
            t = en(k); en(k) = en(kk); en(kk) = t
            ccol = c(:, k); c(:, k) = c(:, kk); c(:, kk) = ccol
         END IF
      END DO
      ALLOCATE (wc(0:nlmd, nf, nat, n), wp(ng, n), kcb(0:nlmd, nf, nat), kpb(ng), dp(n), dm(n))
      DO k = 1, n
         wc(:, :, :, k) = CMPLX(0.0, 0.0)
         wp(:, k) = CMPLX(0.0, 0.0)
         DO j = 1, n
            wc(:, :, :, k) = wc(:, :, :, k) + c(j, k)*cof(i1 - 1 + j, :, :, :)
            wp(:, k) = wp(:, k) + c(j, k)*cpw(:, i1 - 1 + j)
         END DO
         ! the two tests of a real function, each part against its own largest coefficient
         CALL lapwfn_conj(twok, wc(:, :, :, k), wp(:, k), kcb, kpb)
         bmax = MAX(MAXVAL(ABS(wc(:, :, :, k))), MAXVAL(ABS(wp(:, k))))
         dp(k) = MAXVAL(ABS(kpb - wp(:, k)))
         dm(k) = MAXVAL(ABS(kcb - wc(:, :, :, k)))
         IF (dp(k) > real_tol*MAXVAL(ABS(wp(:, k))) .OR. dm(k) > real_tol*MAXVAL(ABS(wc(:, :, :, k)))) THEN
            why = "a combination fails the tests of a real function: left as FLEUR gives it"
            RETURN
         END IF
         dp(k) = dp(k)/bmax
         dm(k) = dm(k)/bmax
      END DO
      ! the level replaced
      dlab = MAXVAL(ABS(en - eig(i1:i2)))
      tmat(i1:i2, i1:i2) = CMPLX(0.0, 0.0)
      DO k = 1, n
         cof(i1 - 1 + k, :, :, :) = wc(:, :, :, k)
         cpw(:, i1 - 1 + k) = wp(:, k)
         eigw(i1 - 1 + k) = en(k)
         lreal(i1 - 1 + k) = .TRUE.
         devp(i1 - 1 + k) = dp(k)
         devm(i1 - 1 + k) = dm(k)
         DO j = 1, n
            tmat(i1 - 1 + j, i1 - 1 + k) = c(j, k)*phs(i1 - 1 + j)
         END DO
      END DO
      phs(i1:i2) = CMPLX(1.0, 0.0)
      ok = .TRUE.
   END SUBROUTINE lapwfn_realise

END MODULE m_writeLapwfn
