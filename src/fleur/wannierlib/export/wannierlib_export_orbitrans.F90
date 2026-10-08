!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The four structural files the orbitrans transport code reads, written from the
!>  quantities FLEUR already holds.
!>
!>  orbitrans takes its geometry from `input/lattice`, `input/atom_centers`,
!>  `input/wann_centers` and `input/N_wann_atom`. Those were recovered with regular
!>  expressions from a Wannier90 `.wout`, which the library mode never writes: it hands
!>  Wannier90 an internal seedname and redirects its output into `out`. A library-mode
!>  run therefore could not feed the transport code without borrowing the four files
!>  from a different calculation.
!>
!>  All four are in memory to full precision rather than the six decimals a printed
!>  table carries, which is the other half of the reason to write them here.
!>
!>  Only the four are written. The operator matrices orbitrans also reads are the same
!>  arrays <operators_r> already exports.
MODULE m_wannierlib_export_orbitrans
   USE m_juDFT
   USE m_constants, ONLY: oUnit, bohr_to_angstrom_const
   USE m_types_atoms
   USE m_types_cell
   USE m_types_mpi
   USE m_types_wgauge_bmesh, ONLY: t_wgauge_bmesh
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannierlib_export_orbitrans

   !> A directory of their own: the names orbitrans uses are generic enough to collide
   !> with anything else sitting in the run directory.
   CHARACTER(LEN=*), PARAMETER :: ORB_DIR = 'orbitrans_input'

CONTAINS

   !> Write this spin channel's share of the four files. Rank 0 only.
   !>
   !> A collinear two-channel run wannierises each channel separately and the transport
   !> code reads the two as one basis of 2*num_wann functions, so the centres are laid
   !> down one channel at a time and the per-atom counts are only known once the last
   !> channel has been seen. Nothing is carried between calls: which slice this channel
   !> owns follows from its index, and the total follows from the channel count.
   SUBROUTINE wannierlib_export_orbitrans(this, atoms, cell, bmesh, fmpi, channel, n_channels)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: channel      !< 1-based index of the spin channel being written
      INTEGER, INTENT(IN) :: n_channels   !< how many channels the run will write in total

      INTEGER, ALLOCATABLE :: nw_atom(:)

      IF (fmpi%irank /= 0) RETURN

      !> The centres are the only one of the four that is computed rather than copied,
      !> and they exist only once Wannier90 has wannierised.
      IF (.NOT. ALLOCATED(bmesh%centres)) THEN
         WRITE (oUnit, '(a)') 'wannierlib: <export orbitrans="T"> needs the Wannier centres '// &
            'and none were reported -- skipped'
         RETURN
      END IF

      !> The cell and the atoms do not depend on the spin channel, so the first one writes
      !> them and the rest leave them alone.
      IF (channel == 1) THEN
         CALL orbitrans_make_dir()
         CALL orbitrans_check_lattice_rows(atoms, cell)
         CALL orbitrans_write_lattice(cell)
         CALL orbitrans_write_atom_centers(atoms)
      END IF

      CALL orbitrans_write_wann_centers(bmesh, this%num_wann, channel)

      !> Last channel: now the total Wannier count is known.
      IF (channel == n_channels) THEN
         CALL orbitrans_count_wann_per_atom(this, n_channels, nw_atom)
         CALL orbitrans_write_n_wann_atom(nw_atom)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote '//ORB_DIR// &
            '/{lattice,atom_centers,wann_centers,N_wann_atom} for ', &
            n_channels*this%num_wann, ' Wannier functions'
      END IF
   END SUBROUTINE wannierlib_export_orbitrans

   !> Fortran has no portable mkdir, so this is the one shell call in the export and it
   !> is isolated here. A failure is raised rather than left to the opens below, which
   !> would otherwise report a missing file and hide the cause.
   SUBROUTINE orbitrans_make_dir()
      INTEGER :: stat
      CALL EXECUTE_COMMAND_LINE('mkdir -p '//ORB_DIR, exitstat=stat)
      IF (stat /= 0) CALL juDFT_error('wannierlib: could not create '//ORB_DIR, &
                                      calledby='wannierlib_export_orbitrans')
   END SUBROUTINE orbitrans_make_dir

   !> Does a row of the lattice file really hold a Bravais vector?
   !>
   !> orbitrans rebuilds a real-space vector as rpts = n1*a_vecs(1,:) + n2*a_vecs(2,:) +
   !> n3*a_vecs(3,:), so row i of the file has to be Bravais vector i. Here that row is
   !> amat(:,i), a COLUMN of amat, because amat is what turns a fractional coordinate
   !> cartesian: pos = MATMUL(amat, taual) makes a position a combination of amat's
   !> columns, so the columns are the vectors.
   !>
   !> Writing the rows instead would be a transpose, and a cubic cell cannot see it --
   !> its amat is symmetric, which is how the same mistake survived for months in the
   !> k-path code. So the check runs at every export rather than in a test case: rebuild
   !> each atom the way the consumer rebuilds a lattice vector and compare against the
   !> position FLEUR already holds. On a cubic cell it passes either way and costs
   !> nothing; on the first cell that is not, it fires.
   SUBROUTINE orbitrans_check_lattice_rows(atoms, cell)
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER :: ia, i
      REAL :: rebuilt(3), worst

      worst = 0.0
      DO ia = 1, atoms%nat
         rebuilt = 0.0
         DO i = 1, 3
            rebuilt = rebuilt + atoms%taual(i, ia)*cell%amat(:, i)
         END DO
         worst = MAX(worst, MAXVAL(ABS(rebuilt - atoms%pos(:, ia))))
      END DO

      IF (worst > 1.0e-8) &
         CALL juDFT_error('wannierlib: <export orbitrans="T"> would write a transposed lattice', &
                          hint='rebuilding the atomic positions from the rows this export writes '// &
                          'does not reproduce atoms%pos, so a row is not a Bravais vector', &
                          calledby='wannierlib_export_orbitrans')
   END SUBROUTINE orbitrans_check_lattice_rows

   !> One Bravais vector per line, cartesian, in BOHR: orbitrans multiplies the file by
   !> the Bohr radius on the way in.
   !>
   !> amat carries the vectors as COLUMNS -- pos = MATMUL(cell%amat, taual) is what turns
   !> a fractional coordinate cartesian, so amat(:,i) is Bravais vector i, and the file
   !> wants it as row i. A cubic cell cannot tell the two apart because its amat is
   !> symmetric, so no bulk case would catch this being transposed.
   SUBROUTINE orbitrans_write_lattice(cell)
      TYPE(t_cell), INTENT(IN) :: cell
      INTEGER :: iu, i
      OPEN (newunit=iu, file=ORB_DIR//'/lattice', status='replace')
      DO i = 1, 3
         WRITE (iu, '(3(2x,f18.10))') cell%amat(:, i)
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_lattice

   !> Index then cartesian position, in BOHR, one atom per line. atoms%pos is cartesian
   !> and already in Bohr; taual is the fractional one.
   SUBROUTINE orbitrans_write_atom_centers(atoms)
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER :: iu, ia
      OPEN (newunit=iu, file=ORB_DIR//'/atom_centers', status='replace')
      DO ia = 1, atoms%nat
         WRITE (iu, '(i5,3(2x,f18.10))') ia, atoms%pos(:, ia)
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_atom_centers

   !> Index then centre, one Wannier function per line. The first channel opens the file
   !> and the rest append after it, so the indices run over the whole basis.
   !>
   !> Written in BOHR, so all four structural files share one convention: orbitrans reads
   !> `lattice`, `atom_centers` and `wann_centers` as Bohr and multiplies each by the Bohr
   !> radius on the way in (io_module.f90, lines 401, 563 and 597). The first two already
   !> hold Bohr -- cell%amat and atoms%pos are FLEUR's own -- but w90_get_centres returns
   !> Angstrom, so the centres are the one quantity that has to be converted, and this is
   !> the only place it happens. Everything else, inside FLEUR and in what Wannier90
   !> prints, keeps Angstrom.
   !>
   !> This makes the file differ from the hand conversion it replaces, which passed the
   !> Angstrom values straight through and so came out short by the Bohr radius. That was
   !> invisible in every case run so far: the file feeds nothing but the assignment of
   !> Wannier functions to atoms (atom_projection_index.f90 ranks |r_atom - r_wann| and
   !> gives atom n its N_wann_atom(n) nearest), and with one atom in the cell every
   !> function is assigned to it whatever the scale. With two atoms the mismatch -- the
   !> atom centres correctly scaled, the Wannier centres short by 0.529 -- would decide
   !> the assignment, and with it everything resolved per atom.
   SUBROUTINE orbitrans_write_wann_centers(bmesh, num_wann, channel)
      TYPE(t_wgauge_bmesh), INTENT(IN) :: bmesh
      INTEGER, INTENT(IN) :: num_wann, channel
      INTEGER :: iu, n, offset

      IF (SIZE(bmesh%centres, 2) /= num_wann) &
         CALL juDFT_error('wannierlib: <export orbitrans="T"> got a centre count that is not num_wann', &
                          calledby='wannierlib_export_orbitrans')

      offset = (channel - 1)*num_wann
      IF (channel == 1) THEN
         OPEN (newunit=iu, file=ORB_DIR//'/wann_centers', status='replace')
      ELSE
         OPEN (newunit=iu, file=ORB_DIR//'/wann_centers', status='old', position='append')
      END IF
      DO n = 1, num_wann
         WRITE (iu, '(i5,3(2x,f18.10))') offset + n, bmesh%centres(:, n)/bohr_to_angstrom_const
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_wann_centers

   !> How many Wannier functions sit on each atom, one count per line, summed over the
   !> spin channels because that is the basis orbitrans sees.
   !>
   !> The projection list is already expanded one entry per atom, so each entry carries
   !> its atom. A spinor run turns one projection into more than one Wannier function;
   !> that ratio is read off num_wann instead of assumed, and a ratio that is not whole
   !> means the projections do not account for the Wannier functions -- which would make
   !> the file wrong in a way orbitrans has no way to detect.
   SUBROUTINE orbitrans_count_wann_per_atom(this, n_channels, nw_atom)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      INTEGER, INTENT(IN) :: n_channels
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nw_atom(:)
      INTEGER :: nproj, nat, ip, ia, per_proj

      nproj = SIZE(this%proj)
      IF (nproj == 0) CALL juDFT_error('wannierlib: <export orbitrans="T"> found no projections', &
                                       calledby='wannierlib_export_orbitrans')
      IF (MOD(this%num_wann, nproj) /= 0) &
         CALL juDFT_error('wannierlib: <export orbitrans="T"> cannot split num_wann over the projections', &
                          calledby='wannierlib_export_orbitrans')
      per_proj = this%num_wann/nproj
      nat = MAXVAL(this%proj(:)%atom)

      ALLOCATE (nw_atom(nat))
      nw_atom = 0
      DO ip = 1, nproj
         ia = this%proj(ip)%atom
         nw_atom(ia) = nw_atom(ia) + per_proj*n_channels
      END DO
   END SUBROUTINE orbitrans_count_wann_per_atom

   SUBROUTINE orbitrans_write_n_wann_atom(nw_atom)
      INTEGER, INTENT(IN) :: nw_atom(:)
      INTEGER :: iu, ia
      OPEN (newunit=iu, file=ORB_DIR//'/N_wann_atom', status='replace')
      DO ia = 1, SIZE(nw_atom)
         WRITE (iu, '(i5)') nw_atom(ia)
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_n_wann_atom

END MODULE m_wannierlib_export_orbitrans
