!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The Wannier90 input file, WF<spin>.win, written beside the .amn, .mmn and .eig so
!> that the three of them can be replayed by a standalone wannier90.x.
!>
!> Wannier90 cannot write this file for us. Its own dumper, w90_readwrite_write_win
!> (readwrite.F90:4948), walks the option store and prints what it finds there, which
!> leaves out exactly what a reader cannot do without: the atoms never enter the store
!> -- atoms_frac is a block of character data and w90_set_option takes none, which is
!> why the call that would set it is commented out in the adapter -- and the dump ends
!> at the first "dump_inputs" or "distk" entry, because the loop leaves on them with
!> EXIT rather than CYCLE. It also names the file <seed>.win_dump. So the writer is
!> ours, and it writes the same values the adapter hands the library, in the same
!> units, so that a standalone run and the library run see one setup and not two.
!>
!> The content follows the classic path, wann_write_win in wann_wan90prep.F of the
!> MaX5.1 tree, which has written this file since 2006. What is new here is that the
!> four energy windows are written as numbers rather than as commented-out placeholders:
!> they are known by the time this runs, because wannierlib_default_windows has already
!> resolved the ones the inp.xml left out.
MODULE m_wannierlib_export_win
   USE m_juDFT
   USE m_constants, ONLY: namat_const, hartree_to_ev_const, bohr_to_angstrom_const
   USE m_types_atoms
   USE m_types_cell
   USE m_types_kpts
   USE m_wgauge_mesh, ONLY: wgauge_mp_grid
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannierlib_write_win

CONTAINS

   !> Write WF<spin>.win. Called on one rank only: everything below is already whole.
   !>
   !> search_shells comes from the adapter rather than from a number repeated here, so
   !> that a standalone run searches the same list of shells as the library did and
   !> therefore settles on the same b vectors -- without that, the .mmn beside this file
   !> would be a table of overlaps for neighbours the reader never chose.
   SUBROUTINE wannierlib_write_win(this, atoms, cell, kpts, filename, l_spinors, search_shells)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_cell), INTENT(IN)  :: cell
      TYPE(t_kpts), INTENT(IN)  :: kpts
      CHARACTER(LEN=*), INTENT(IN) :: filename
      LOGICAL, INTENT(IN) :: l_spinors
      INTEGER, INTENT(IN) :: search_shells

      INTEGER :: iunit, io, na, itype, idir, ikpt
      INTEGER :: mp_grid(3)
      CHARACTER(LEN=2) :: atom_symbol

      OPEN (NEWUNIT=iunit, FILE=filename, FORM='formatted', STATUS='replace', IOSTAT=io)
      IF (io /= 0) CALL juDFT_error('wannierlib: cannot write '//TRIM(filename), &
                                    calledby='wannierlib_write_win')

      WRITE (iunit, '(a)') '! Written by FLEUR/wannierlib, to go with the .amn, .mmn and'
      WRITE (iunit, '(a)') '! .eig beside it. Energies in eV, lengths in Angstrom.'
      WRITE (iunit, '(a)') ''
      WRITE (iunit, '(a,i0)') 'num_wann  = ', this%num_wann
      WRITE (iunit, '(a,i0)') 'num_bands = ', this%num_bands
      IF (l_spinors) WRITE (iunit, '(a)') 'spinors   = .true.'
      WRITE (iunit, '(a,i0)') 'search_shells = ', search_shells
      IF (this%kmesh_tol >= 0.0) WRITE (iunit, '(a,es22.14)') 'kmesh_tol = ', this%kmesh_tol
      WRITE (iunit, '(a)') ''

      !> The minimisation, in both its stages. A negative value is the sentinel for "not
      !> given", and leaving the keyword out is what hands the reader Wannier90's default
      !> -- the same thing silence does on our side.
      IF (this%num_iter >= 0)     WRITE (iunit, '(a,i0)') 'num_iter      = ', this%num_iter
      IF (this%dis_num_iter >= 0) WRITE (iunit, '(a,i0)') 'dis_num_iter  = ', this%dis_num_iter
      IF (this%conv_tol > 0.0)     WRITE (iunit, '(a,es22.14)') 'conv_tol      = ', this%conv_tol
      IF (this%dis_conv_tol > 0.0) WRITE (iunit, '(a,es22.14)') 'dis_conv_tol  = ', this%dis_conv_tol
      WRITE (iunit, '(a)') ''

      !> The four windows, as numbers. They are the effective ones: whatever the inp.xml
      !> left out was filled in from the bands before this runs, so the file says what the
      !> wannierisation actually used and not what somebody typed.
      WRITE (iunit, '(a,es22.14)') 'dis_win_min  = ', hartree_to_ev_const*this%dis_win_min
      WRITE (iunit, '(a,es22.14)') 'dis_win_max  = ', hartree_to_ev_const*this%dis_win_max
      WRITE (iunit, '(a,es22.14)') 'dis_froz_min = ', hartree_to_ev_const*this%dis_froz_min
      WRITE (iunit, '(a,es22.14)') 'dis_froz_max = ', hartree_to_ev_const*this%dis_froz_max
      WRITE (iunit, '(a)') ''

      !> Angstrom, because that is the unit the adapter converts the cell into before
      !> handing it to the library. Writing Bohr here would be correct on its own and
      !> wrong beside the .mmn.
      WRITE (iunit, '(a)') 'begin unit_cell_cart'
      WRITE (iunit, '(a)') 'Ang'
      DO idir = 1, 3
         WRITE (iunit, '(3f20.12)') bohr_to_angstrom_const*cell%amat(:, idir)
      END DO
      WRITE (iunit, '(a)') 'end unit_cell_cart'
      WRITE (iunit, '(a)') ''

      WRITE (iunit, '(a)') 'begin atoms_frac'
      DO na = 1, atoms%nat
         itype = atoms%itype(na)
         atom_symbol = ADJUSTL(namat_const(atoms%nz(itype)))
         WRITE (iunit, '(a2,1x,3f20.12)') atom_symbol, atoms%taual(:, na)
      END DO
      WRITE (iunit, '(a)') 'end atoms_frac'
      WRITE (iunit, '(a)') ''

      mp_grid = wgauge_mp_grid(kpts)
      WRITE (iunit, '(a,3(1x,i0))') 'mp_grid =', mp_grid(1), mp_grid(2), mp_grid(3)
      WRITE (iunit, '(a)') ''

      !> The full-zone list, in the order the .amn, .mmn and .eig are written in, because
      !> all four files index k the same way and nothing in them says so.
      WRITE (iunit, '(a)') 'begin kpoints'
      DO ikpt = 1, kpts%nkptf
         WRITE (iunit, '(3f20.14)') kpts%bkf(:, ikpt)
      END DO
      WRITE (iunit, '(a)') 'end kpoints'

      CLOSE (iunit)
   END SUBROUTINE wannierlib_write_win

END MODULE m_wannierlib_export_win
