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
!>  arrays <operatorsRs> already exports.
!>
!>  The 2N assembly a collinear run owes the same consumer lives here too. It is a
!>  different kind of work -- it reads back the rendered files rather than FLEUR's own
!>  arrays, and needs none of the types above -- but it exists only because orbitrans
!>  reads one 18x18 where a collinear run rendered two 9x9, so it belongs to this export
!>  and to no other. Its own account is above its entry point below.
MODULE m_wannierlib_export_orbitrans
   USE m_juDFT
   USE m_constants, ONLY: oUnit, bohr_to_angstrom_const, hartree_to_ev_const
   USE m_types_atoms
   USE m_types_cell
   USE m_types_mpi
   USE m_types_wgauge_bmesh, ONLY: t_wgauge_bmesh
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   USE m_wannierlib_export_paths, ONLY: wannierlib_export_dir, wannierlib_export_file
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannierlib_export_orbitrans, orbitrans_assemble_2n

   !> The same directory the operator matrices were written into, which for this flag is
   !> `input/`: orbitrans opens all of them by a relative path and expects them together.
   !> Taken from the resolver rather than repeated here, so there is one name for it.

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
         CALL orbitrans_check_lattice_rows(atoms, cell)
         CALL orbitrans_write_lattice(cell)
         CALL orbitrans_write_atom_centers(atoms)
         CALL orbitrans_write_kpts(this)
      END IF

      CALL orbitrans_write_wann_centers(bmesh, this%num_wann, channel)

      !> Last channel: now the total Wannier count is known.
      IF (channel == n_channels) THEN
         CALL orbitrans_count_wann_per_atom(this, atoms, n_channels, nw_atom)
         CALL orbitrans_write_n_wann_atom(nw_atom)
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote '//wannierlib_export_dir()// &
            '/{lattice,atom_centers,wann_centers,N_wann_atom} for ', &
            n_channels*this%num_wann, ' Wannier functions'
      END IF
   END SUBROUTINE wannierlib_export_orbitrans

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
      OPEN (newunit=iu, file=wannierlib_export_dir()//'/lattice', status='replace')
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
      OPEN (newunit=iu, file=wannierlib_export_dir()//'/atom_centers', status='replace')
      DO ia = 1, atoms%nat
         WRITE (iu, '(i5,3(2x,f18.10))') ia, atoms%pos(:, ia)
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_atom_centers

   !> The k path orbitrans resolves on, in the format its own kptgen expects.
   !>
   !> kptgen reads `input/kpts` as `N_kpts scale` and then one k per line, and rebuilds
   !> each one as kpt(1)*b_vecs(1,:) + kpt(2)*b_vecs(2,:) + kpt(3)*b_vecs(3,:) over scale.
   !> Those are internal coordinates, which is what bk already holds, so the list goes out
   !> as it stands and scale is 1.
   !>
   !> The path is not invented here. It is the <domain> list the run already names, after
   !> @refinementPoints has subdivided it, and that is the same set wannierlib interpolates
   !> its own bands on -- so the transport code and FLEUR sample the same k, and N_kpts on
   !> line 1 cannot disagree with what was interpolated. Choosing a path by lattice instead
   !> would be a choice made for the caller, and a wrong one still runs to completion.
   !>
   !> Without a <domain> there is no path to take and the file is left alone. Only the kres
   !> drivers read it; a mesh run never opens it.
   !>
   !> dom_kset is filled on rank 0 alone, which is the only rank that reaches this.
   SUBROUTINE orbitrans_write_kpts(this)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      INTEGER :: iu, ik, idom, ipick

      IF (this%n_domains < 1) THEN
         WRITE (oUnit, '(a)') 'wannierlib: no '//wannierlib_export_dir()// &
            '/kpts written -- the run declares no <domain> to take a k path from'
         RETURN
      END IF

      !> @suffix may be left off for at most one domain, so a domain without one is the
      !> distinguished set of the run, and the natural source for the single file orbitrans
      !> reads. With a suffix on every domain there is no distinguished one and the first
      !> is taken, which the log names so the choice is never silent.
      ipick = 1
      DO idom = 1, this%n_domains
         IF (LEN_TRIM(this%dom_suffix(idom)) == 0) THEN
            ipick = idom
            EXIT
         END IF
      END DO

      OPEN (newunit=iu, file=wannierlib_export_dir()//'/kpts', status='replace')
      WRITE (iu, '(2x,i0,2x,f10.5)') this%dom_kset(ipick)%nkpt, 1.0
      DO ik = 1, this%dom_kset(ipick)%nkpt
         WRITE (iu, '(3(2x,f18.12))') this%dom_kset(ipick)%bk(:, ik)
      END DO
      CLOSE (iu)

      IF (LEN_TRIM(this%dom_suffix(ipick)) == 0) THEN
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote '//wannierlib_export_dir()// &
            '/kpts with ', this%dom_kset(ipick)%nkpt, ' k-points, from the <domain> with no suffix'
      ELSE
         WRITE (oUnit, '(a,i0,a)') 'wannierlib: wrote '//wannierlib_export_dir()// &
            '/kpts with ', this%dom_kset(ipick)%nkpt, ' k-points, from <domain> suffix "'// &
            TRIM(this%dom_suffix(ipick))//'"'
      END IF
   END SUBROUTINE orbitrans_write_kpts

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
         OPEN (newunit=iu, file=wannierlib_export_dir()//'/wann_centers', status='replace')
      ELSE
         OPEN (newunit=iu, file=wannierlib_export_dir()//'/wann_centers', status='old', position='append')
      END IF
      DO n = 1, num_wann
         WRITE (iu, '(i5,3(2x,f18.10))') offset + n, bmesh%centres(:, n)/bohr_to_angstrom_const
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_wann_centers

   !> How many Wannier functions sit on each atom, one count per line, summed over the
   !> spin channels because that is the basis orbitrans sees.
   !>
   !> A projection carries its atom as the pair (ntype, atom), where `atom` counts within
   !> the equivalence class of its type -- which is what builds A_mn, since the type is
   !> what reaches the radial functions. The index of the atom in the cell is
   !> firstAtom(ntype) + atom - 1, and this file is per atom of the cell, so that is what
   !> indexes it and atoms%nat is what sizes it.
   !>
   !> Reading `atom` as that index instead agrees with it whenever every class holds one
   !> atom, which is every cell of one species -- so a monolayer cannot tell the two
   !> apart. A film of five layers has classes of two, and the file came out with one
   !> line per distinct class size: two counts for five atoms, and the consumer read off
   !> the end of it. The assertion below is what makes that visible here instead.
   !>
   !> A spinor run turns one projection into more than one Wannier function; that ratio is
   !> read off num_wann instead of assumed, and a ratio that is not whole means the
   !> projections do not account for the Wannier functions.
   SUBROUTINE orbitrans_count_wann_per_atom(this, atoms, n_channels, nw_atom)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: n_channels
      INTEGER, ALLOCATABLE, INTENT(OUT) :: nw_atom(:)
      INTEGER :: nproj, ip, ia, per_proj

      nproj = SIZE(this%proj)
      IF (nproj == 0) CALL juDFT_error('wannierlib: <export orbitrans="T"> found no projections', &
                                       calledby='wannierlib_export_orbitrans')
      IF (MOD(this%num_wann, nproj) /= 0) &
         CALL juDFT_error('wannierlib: <export orbitrans="T"> cannot split num_wann over the projections', &
                          calledby='wannierlib_export_orbitrans')
      per_proj = this%num_wann/nproj

      ALLOCATE (nw_atom(atoms%nat))
      nw_atom = 0
      DO ip = 1, nproj
         ia = atoms%firstAtom(this%proj(ip)%ntype) + this%proj(ip)%atom - 1
         IF ((ia < 1) .OR. (ia > atoms%nat)) &
            CALL juDFT_error('wannierlib: <export orbitrans="T"> built an atom index '// &
                             'outside the cell', calledby='wannierlib_export_orbitrans')
         nw_atom(ia) = nw_atom(ia) + per_proj*n_channels
      END DO

      IF (SUM(nw_atom) /= n_channels*this%num_wann) &
         CALL juDFT_error('wannierlib: <export orbitrans="T"> per-atom counts do not sum '// &
                          'to the Wannier functions there are', &
                          hint='the projections do not account for every Wannier function', &
                          calledby='wannierlib_export_orbitrans')
   END SUBROUTINE orbitrans_count_wann_per_atom

   !> One line per atom, in atom order, with how many Wannier functions sit on it. It is
   !> what lets orbitrans cut the Wannier index into per-atom blocks, and the file has to
   !> exist even where every atom carries the same number.
   SUBROUTINE orbitrans_write_n_wann_atom(nw_atom)
      INTEGER, INTENT(IN) :: nw_atom(:)
      INTEGER :: iu, ia
      OPEN (newunit=iu, file=wannierlib_export_dir()//'/N_wann_atom', status='replace')
      DO ia = 1, SIZE(nw_atom)
         WRITE (iu, '(i5)') nw_atom(ia)
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_n_wann_atom


   !===============================================================================
   ! The 2N assembly: what a collinear run owes orbitrans that a spinorial one does not.
   !===============================================================================

   !>  Build the 2N Hamiltonian and position operator a collinear run owes orbitrans.
   !>
   !>  A collinear run wannierises each spin channel on its own, so the operator layer
   !>  renders two 9x9 objects, WF1_* and WF2_*, while orbitrans reads one 18x18. That
   !>  combination is physics, not a rename, so it cannot live in the renderer: the same
   !>  rendering also feeds JuWIP and the standalone Wannier90 files, and those want the
   !>  channels apart. It is assembled here instead, in the orbitrans consumer itself,
   !>  from the files the renderer has already written.
   !>
   !>    hoppings(R) = diag( H_up(R), H_down(R) ) + Htr * SOC(R)
   !>    position(R) = diag( A_up(R), A_down(R) )
   !>    bmn(R), cmn(R), fmn(R) = the same diagonal stack, with no extra term
   !>
   !>  in the order [up x num_wann, down x num_wann], which is the one the diagonal of
   !>  S_z(R=0) shows and the one the spin indices of rssocmat.1 follow.
   !>
   !>  The spin-orbit term belongs INSIDE hoppings because orbitrans subtracts rather
   !>  than adds: read_spin_orbit ends with H_nosoc = H_tot - H_soc, so it takes the
   !>  hopping file for the full Hamiltonian and uses spin_orbit only to split the two
   !>  apart. Handing it a spin-diagonal H instead leaves no matrix element between
   !>  states of equal spin and sigma^Sx and sigma^Sy come out identically zero, with no
   !>  error and with the orbital channel still looking healthy.
   !>
   !>  rssocmat.1 is in Hartree and the hopping file in eV, hence the conversion. Its
   !>  orbital indices already run over all 2N functions and the (s1,s2) pair sits on
   !>  top of them, so the element that counts is the one whose spin labels match those
   !>  of the two functions. Treating the pair as a block offset instead leaves the
   !>  up-down block untouched, which is the same silent zero again.

   !> Replace the channel-one files in the export directory by their 2N assembly.
   !> Rank 0 only, and only once the last channel has been rendered.
   !>
   !> Everything below this point reads back what the per-channel export rendered and
   !> writes the result in the same format. It takes nothing from FLEUR's own arrays, and
   !> it must not: sharing a module with the structural writers puts their type modules in
   !> scope, so the compiler no longer refuses a reach for atoms or cell.
   SUBROUTINE orbitrans_assemble_2n(num_wann)
      INTEGER, INTENT(IN) :: num_wann          !< functions per spin channel

      COMPLEX, ALLOCATABLE :: oup(:, :, :, :), odn(:, :, :, :), o2n(:, :, :, :)
      INTEGER, ALLOCATABLE :: rv(:, :), rv2(:, :)
      INTEGER :: nrpts, nrpts2, nw2
      LOGICAL :: l_soc, l_h, l_a

      nw2 = 2*num_wann

      !> Only what the run actually rendered. The consumer flag does not decide the
      !> operator list -- <operatorsRs> does -- so a run that asked for the Hamiltonian
      !> and not for the position leaves nothing to stack for the second one. Absent
      !> means not asked for; present on one channel and missing on the other is a
      !> broken run and says so below.
      INQUIRE (file=wannierlib_export_file('WF1_hr.dat'), exist=l_h)
      INQUIRE (file=wannierlib_export_file('WF1_r.dat'), exist=l_a)

      !> --- H(R). The channel-one file already carries the orbitrans name.
      IF (l_h) THEN
         CALL orbitrans_read_op(wannierlib_export_file('WF1_hr.dat'), num_wann, 1, .FALSE., .FALSE., &
                                nrpts, rv, oup)
         CALL orbitrans_read_op(wannierlib_export_file('WF2_hr.dat'), num_wann, 1, .TRUE., .TRUE., &
                                nrpts2, rv2, odn)
         CALL orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, 'WF2_hr.dat')
         CALL orbitrans_stack(oup, odn, num_wann, nrpts, 1, o2n)

         INQUIRE (file=wannierlib_export_file('rssocmat.1'), exist=l_soc)
         IF (l_soc) THEN
            CALL orbitrans_add_soc(o2n, num_wann, nrpts, rv)
         ELSE
            WRITE (oUnit, '(a)') 'wannierlib: WARNING -- no spin-orbit operator to fold into '// &
               'hoppings. The Hamiltonian stays diagonal in spin and the transverse spin '// &
               'conductivities will come out identically zero.'
         END IF
         CALL orbitrans_write_op(wannierlib_export_file('WF1_hr.dat'), nw2, nrpts, 1, rv, o2n)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: assembled hoppings as ', nw2, 'x', nw2, &
            ' from the two spin channels'
      END IF

      !> --- A(R), the same stacking and no extra term. Its R mesh is checked against the
      !> file itself rather than against the Hamiltonian's, which a run without one has not
      !> read.
      IF (l_a) THEN
         CALL orbitrans_read_op(wannierlib_export_file('WF1_r.dat'), num_wann, 3, .FALSE., .FALSE., &
                                nrpts, rv, oup)
         CALL orbitrans_read_op(wannierlib_export_file('WF2_r.dat'), num_wann, 3, .TRUE., .FALSE., &
                                nrpts2, rv2, odn)
         CALL orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, 'WF2_r.dat')
         CALL orbitrans_stack(oup, odn, num_wann, nrpts, 3, o2n)
         CALL orbitrans_write_op(wannierlib_export_file('WF1_r.dat'), nw2, nrpts, 3, rv, o2n)
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: assembled position as ', nw2, 'x', nw2, &
            ' from the two spin channels'
      END IF

      !> --- B(R), C(R) and F(R), the three the modern theory of orbital magnetisation
      !> reads. Same diagonal stack and no extra term.
      !>
      !> F = <0n|r_a r_b|Rm> is purely spatial and the two channels are orthogonal
      !> spinors, so its up-down block vanishes by algebra: the stack is the whole
      !> operator.
      !>
      !> B = <0n|H (r-R)|Rm> and C = <0n|r_a H r_b|Rm> carry an H, and in a collinear
      !> run that H is the channel's own, without the spin-orbit term -- which is only
      !> known after the channel loop, and is why the Hamiltonian is folded here rather
      !> than rendered that way. So their up-down block stays empty. Measured on bcc Fe
      !> against the spinorial route, that block holds 0.41 % of the norm of B and
      !> 0.56 % of C, against 0.37 % for H itself. Filling it asks for B and C built
      !> from the full Hamiltonian, which belongs in the matrix-element layer and not
      !> in this consumer.
      !>
      !> Channel two keeps the layout its name promises: WF2_bmn.dat still carries the
      !> Wannier90 header, WF2_cmn.dat and WF2_fmn.dat never had one. Getting that wrong
      !> aborts in the reader rather than passing silently.
      CALL orbitrans_stack_one('WF1_bmn.dat', 'WF2_bmn.dat', 'bmn', num_wann, 3, .TRUE.)
      CALL orbitrans_stack_one('WF1_cmn.dat', 'WF2_cmn.dat', 'cmn', num_wann, 9, .FALSE.)
      CALL orbitrans_stack_one('WF1_fmn.dat', 'WF2_fmn.dat', 'fmn', num_wann, 9, .FALSE.)
   END SUBROUTINE orbitrans_assemble_2n

   !> Stack one per-channel operator into the 2N file its consumer reads, and do nothing
   !> when the run did not render it: <operatorsRs> decides the list, not this flag.
   SUBROUTINE orbitrans_stack_one(ch1, ch2, label, nw, ncomp, l_header2)
      CHARACTER(LEN=*), INTENT(IN) :: ch1, ch2, label
      INTEGER, INTENT(IN) :: nw, ncomp
      LOGICAL, INTENT(IN) :: l_header2   !< channel two still carries the Wannier90 header

      COMPLEX, ALLOCATABLE :: oup(:, :, :, :), odn(:, :, :, :), o2n(:, :, :, :)
      INTEGER, ALLOCATABLE :: rv(:, :), rv2(:, :)
      INTEGER :: nrpts, nrpts2
      LOGICAL :: l_have
      CHARACTER(LEN=24) :: fmt

      INQUIRE (file=wannierlib_export_file(ch1), exist=l_have)
      IF (.NOT. l_have) RETURN

      CALL orbitrans_read_op(wannierlib_export_file(ch1), nw, ncomp, .FALSE., .FALSE., &
                             nrpts, rv, oup)
      CALL orbitrans_read_op(wannierlib_export_file(ch2), nw, ncomp, l_header2, .FALSE., &
                             nrpts2, rv2, odn)
      CALL orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, ch2)
      CALL orbitrans_stack(oup, odn, nw, nrpts, ncomp, o2n)

      !> Wide enough for the components these carry, and for the digits the renderer
      !> wrote: the consumer reads them list-directed, so only the values matter.
      WRITE (fmt, '(a,i0,a)') '(5i5,', 2*ncomp, 'f22.8)'
      CALL orbitrans_write_op(wannierlib_export_file(ch1), 2*nw, nrpts, ncomp, rv, o2n, fmt)
      WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: assembled '//TRIM(label)//' as ', 2*nw, &
         'x', 2*nw, ' from the two spin channels'
   END SUBROUTINE orbitrans_stack_one

   !> One real-space operator file into (nw,nw,nrpts,ncomp), with its R vectors.
   !>
   !> The channel-one files have had their header taken off and their degeneracies moved
   !> to `ndegen`, because that is what the consumer reads; WF2_* keeps the Wannier90
   !> header, so the number of R vectors is read from it rather than guessed, and the
   !> degeneracy block that follows is ceil(nrpts/15) lines wide.
   SUBROUTINE orbitrans_read_op(fname, nw, ncomp, l_header, l_ndegen, nrpts, rv, o)
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: nw, ncomp
      LOGICAL, INTENT(IN) :: l_header    !< the Wannier90 three-line header is still on
      LOGICAL, INTENT(IN) :: l_ndegen    !< and behind it the degeneracy block (hr only)
      INTEGER, INTENT(OUT) :: nrpts
      INTEGER, ALLOCATABLE, INTENT(OUT) :: rv(:, :)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: o(:, :, :, :)

      INTEGER :: iu, ios, irpt, i, j, c, k, nw_file, ndlines, ir(3), ii, jj
      REAL :: val(2*ncomp)
      CHARACTER(LEN=256) :: line

      OPEN (newunit=iu, file=TRIM(fname), status='old', action='read', iostat=ios)
      IF (ios /= 0) CALL juDFT_error('<export orbitrans="T"> cannot open '//TRIM(fname)// &
                                     ': a collinear run has to leave both spin channels behind', &
                                     calledby='orbitrans_read_op')

      IF (l_header) THEN
         READ (iu, '(a)') line
         READ (iu, *) nw_file
         READ (iu, *) nrpts
         IF (nw_file /= nw) CALL juDFT_error('channel two of '//TRIM(fname)//' has a '// &
                                             'different number of Wannier functions', calledby='orbitrans_read_op')
         IF (l_ndegen) THEN
            ndlines = (nrpts + 14)/15
            DO k = 1, ndlines
               READ (iu, '(a)') line
            END DO
         END IF
      ELSE
         CALL orbitrans_count_rpts(TRIM(fname), nw, nrpts)
      END IF

      ALLOCATE (o(nw, nw, nrpts, ncomp), rv(3, nrpts))
      DO irpt = 1, nrpts
         DO j = 1, nw
            DO i = 1, nw
               READ (iu, *, iostat=ios) ir, ii, jj, val
               IF (ios /= 0) CALL juDFT_error(TRIM(fname)//' ended early: it does not hold '// &
                                              'nrpts*nw*nw records', calledby='orbitrans_read_op')
               IF (ii /= i .OR. jj /= j) &
                  CALL juDFT_error(TRIM(fname)//' is not in the order this reader assumes '// &
                                   '(R outermost, then column, then row)', calledby='orbitrans_read_op')
               IF (i == 1 .AND. j == 1) rv(:, irpt) = ir
               DO c = 1, ncomp
                  o(i, j, irpt, c) = CMPLX(val(2*c - 1), val(2*c))
               END DO
            END DO
         END DO
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_read_op

   !> How many R vectors a headerless file holds, from its record count.
   SUBROUTINE orbitrans_count_rpts(fname, nw, nrpts)
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: nw
      INTEGER, INTENT(OUT) :: nrpts

      INTEGER :: iu, ios, n
      CHARACTER(LEN=256) :: line

      n = 0
      OPEN (newunit=iu, file=TRIM(fname), status='old', action='read')
      DO
         READ (iu, '(a)', iostat=ios) line
         IF (ios /= 0) EXIT
         IF (LEN_TRIM(line) > 0) n = n + 1
      END DO
      CLOSE (iu)
      IF (MOD(n, nw*nw) /= 0) CALL juDFT_error(TRIM(fname)//' does not hold a whole number '// &
                                               'of nw x nw blocks', calledby='orbitrans_count_rpts')
      nrpts = n/(nw*nw)
   END SUBROUTINE orbitrans_count_rpts

   !> Both channels have to live on the same Wigner-Seitz mesh, or the blocks being put
   !> on one diagonal belong to different R vectors. Same mesh, same order, every vector.
   SUBROUTINE orbitrans_same_mesh(nrpts, nrpts2, rv, rv2, who)
      INTEGER, INTENT(IN) :: nrpts, nrpts2, rv(:, :), rv2(:, :)
      CHARACTER(LEN=*), INTENT(IN) :: who

      INTEGER :: irpt

      IF (nrpts /= nrpts2) CALL juDFT_error('the two spin channels disagree on the number '// &
                                            'of R vectors ('//TRIM(who)//')', calledby='orbitrans_same_mesh')
      DO irpt = 1, nrpts
         IF (ANY(rv(:, irpt) /= rv2(:, irpt))) &
            CALL juDFT_error('the two spin channels list the R vectors in a different order ('// &
                             TRIM(who)//')', calledby='orbitrans_same_mesh')
      END DO
   END SUBROUTINE orbitrans_same_mesh

   !> diag(up, down) in the order [up x nw, down x nw]. The off-diagonal blocks are the
   !> ones the spin-orbit term later fills; without it they stay zero, which is what a
   !> Hamiltonian that does not mix the channels means.
   SUBROUTINE orbitrans_stack(oup, odn, nw, nrpts, ncomp, o2n)
      COMPLEX, INTENT(IN) :: oup(:, :, :, :), odn(:, :, :, :)
      INTEGER, INTENT(IN) :: nw, nrpts, ncomp
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: o2n(:, :, :, :)

      IF (ALLOCATED(o2n)) DEALLOCATE (o2n)
      ALLOCATE (o2n(2*nw, 2*nw, nrpts, ncomp), source=CMPLX(0.0, 0.0))
      o2n(1:nw, 1:nw, :, :) = oup
      o2n(nw + 1:2*nw, nw + 1:2*nw, :, :) = odn
   END SUBROUTINE orbitrans_stack

   !> Fold rssocmat.1 into the stacked Hamiltonian, in eV.
   !>
   !> The file is nrpts x 2N x 2N x 2 x 2: the orbital indices already run over the whole
   !> 2N basis and carry the spinor pair on top. Function i belongs to the up channel when
   !> i <= nw, so the element that is physically its own is the one with s1 and s2 equal to
   !> those labels; the other three of each quartet belong to a different pair of functions.
   SUBROUTINE orbitrans_add_soc(h2n, nw, nrpts, rv)
      COMPLEX, INTENT(INOUT) :: h2n(:, :, :, :)
      INTEGER, INTENT(IN) :: nw, nrpts, rv(:, :)

      INTEGER :: iu, ios, irpt, k, nw2, ir(3), i, j, s1, s2, si, sj
      REAL :: re, im, cross
      CHARACTER(LEN=:), ALLOCATABLE :: fname

      nw2 = 2*nw
      fname = wannierlib_export_file('rssocmat.1')
      OPEN (newunit=iu, file=fname, status='old', action='read', iostat=ios)
      IF (ios /= 0) CALL juDFT_error('cannot open '//fname, calledby='orbitrans_add_soc')

      DO irpt = 1, nrpts
         DO k = 1, nw2*nw2*4
            READ (iu, *, iostat=ios) ir, i, j, s1, s2, re, im
            IF (ios /= 0) CALL juDFT_error(fname//' ended early: it does not hold '// &
                                           'nrpts*2N*2N*4 records', calledby='orbitrans_add_soc')
            IF (ANY(ir /= rv(:, irpt))) &
               CALL juDFT_error(fname//' runs over the R vectors in a different order than '// &
                                'the hopping file', calledby='orbitrans_add_soc')
            si = MERGE(1, 2, i <= nw)
            sj = MERGE(1, 2, j <= nw)
            IF (s1 /= si .OR. s2 /= sj) CYCLE
            h2n(i, j, irpt, 1) = h2n(i, j, irpt, 1) + hartree_to_ev_const*CMPLX(re, im)
         END DO
      END DO
      CLOSE (iu)

      !> The whole point of the fold is the block that connects the two channels. If it is
      !> still empty the spin indices did not line up, and the run would end in a transverse
      !> spin conductivity of exactly zero with nothing to show for it.
      cross = 0.0
      DO irpt = 1, nrpts
         DO j = 1, nw2
            DO i = 1, nw2
               IF ((i <= nw) .NEQV. (j <= nw)) cross = MAX(cross, ABS(h2n(i, j, irpt, 1)))
            END DO
         END DO
      END DO
      IF (cross <= 0.0) CALL juDFT_error('the spin-orbit term left the up-down block of '// &
                                         'hoppings empty, so it did not enter the Hamiltonian', &
                                         calledby='orbitrans_add_soc')
      WRITE (oUnit, '(a,f12.6,a)') 'wannierlib: spin-orbit folded into hoppings, '// &
         'largest up-down element ', cross, ' eV'
   END SUBROUTINE orbitrans_add_soc

   !> Write it back in the layout the file came in: R outermost, then column, then row,
   !> and no header, which is what the consumer reads from the first line.
   SUBROUTINE orbitrans_write_op(fname, nw, nrpts, ncomp, rv, o, fmt)
      CHARACTER(LEN=*), INTENT(IN) :: fname
      INTEGER, INTENT(IN) :: nw, nrpts, ncomp, rv(:, :)
      COMPLEX, INTENT(IN) :: o(:, :, :, :)
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: fmt   !< absent: the six-real layout

      INTEGER :: iu, irpt, i, j, c
      REAL :: val(2*ncomp)
      CHARACTER(LEN=24) :: ffmt

      ffmt = '(5i5,6f12.6)'
      IF (PRESENT(fmt)) ffmt = fmt

      OPEN (newunit=iu, file=TRIM(fname), status='replace')
      DO irpt = 1, nrpts
         DO j = 1, nw
            DO i = 1, nw
               DO c = 1, ncomp
                  val(2*c - 1) = REAL(o(i, j, irpt, c))
                  val(2*c) = AIMAG(o(i, j, irpt, c))
               END DO
               WRITE (iu, TRIM(ffmt)) rv(:, irpt), i, j, val
            END DO
         END DO
      END DO
      CLOSE (iu)
   END SUBROUTINE orbitrans_write_op

END MODULE m_wannierlib_export_orbitrans
