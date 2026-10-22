!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Where the real-space operator files land, and under which name.
!>
!>  <operatorsRs> renders every operator ONCE. A consumer flag in <export> does not add
!>  a second rendering next to the first: it decides the directory the one rendering is
!>  written into and the name it is written under. Nothing is copied and nothing is
!>  moved after the fact -- the writers below already take their filename as an
!>  argument, so the choice is made before the file is opened.
!>
!>  The two consumers hardcode what they read:
!>
!>    orbitrans  opens a relative `input/` and INQUIREs for `spin`, `orbital`,
!>               `spin_orbit`, `hoppings`, `rpts`, `bmn`, `cmn`, `fmn` -- generic names
!>               that carry no spin channel and no seedname.
!>    JuWIP      opens `w90files/` (a hardcoded character(len=9) in its source) and reads
!>               OUR names unchanged, because it was written against the files a FLEUR
!>               run already leaves behind.
!>
!>  So JuWIP is a move and orbitrans is a move plus a rename, and both are this table.
!>
!>  What does NOT get renamed: a second spin channel. A collinear two-channel run
!>  wannierises each channel on its own and orbitrans wants one 2N basis, which is a
!>  physical combination (col2ot.py) and not a name. WF2_* therefore keeps its own name
!>  inside the directory rather than colliding with WF1_* over `hoppings`.
MODULE m_wannierlib_export_paths
   USE m_juDFT
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: wannierlib_export_dir_set, wannierlib_export_file, wannierlib_export_dir

   !> Empty until a run sets it, which is also what a run without a consumer flag keeps:
   !> the canonical names in the run directory, exactly as before this table existed.
   CHARACTER(LEN=:), ALLOCATABLE :: out_dir
   LOGICAL :: rename_for_orbitrans = .FALSE.

CONTAINS

   !> Decide the directory once per run and create it. Rank 0 creates, every rank records,
   !> because the name is needed wherever a filename is built and only rank 0 opens.
   SUBROUTINE wannierlib_export_dir_set(wl, irank)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: wl
      INTEGER, INTENT(IN) :: irank
      INTEGER :: stat

      !> Two consumers want the same rendering in two different places, and the whole
      !> point of this table is that there is only one. Saying so is better than
      !> silently serving one of them.
      IF (wl%export%orbitrans .AND. wl%export%juwip) &
         CALL juDFT_error('wannierlib: <export orbitrans="T"/> and juwip="T" both move the '// &
                          'operator files, and they move them to different directories. '// &
                          'Ask for one of them.', calledby='wannierlib_export_dir_set')

      rename_for_orbitrans = wl%export%orbitrans
      IF (wl%export%orbitrans) THEN
         out_dir = 'input'
      ELSE IF (wl%export%juwip) THEN
         out_dir = 'w90files'
      ELSE
         out_dir = ''
         RETURN
      END IF

      IF (irank == 0) THEN
         !> Fortran has no portable mkdir, so this is a shell call, as in the orbitrans
         !> export. Raising here beats letting every OPEN below fail on a missing path.
         CALL EXECUTE_COMMAND_LINE('mkdir -p '//out_dir, exitstat=stat)
         IF (stat /= 0) CALL juDFT_error('wannierlib: could not create '//out_dir//'/', &
                                         calledby='wannierlib_export_dir_set')
      END IF
   END SUBROUTINE wannierlib_export_dir_set

   !> The directory in force, '' when the files stay in the run directory. For the export
   !> that writes structural files of its own and has to land in the same place.
   FUNCTION wannierlib_export_dir() RESULT(d)
      CHARACTER(LEN=:), ALLOCATABLE :: d
      d = ''
      IF (ALLOCATED(out_dir)) d = out_dir
   END FUNCTION wannierlib_export_dir

   !> Canonical name in, path to open out.
   FUNCTION wannierlib_export_file(canonical) RESULT(path)
      CHARACTER(LEN=*), INTENT(IN) :: canonical
      CHARACTER(LEN=:), ALLOCATABLE :: path

      path = TRIM(canonical)
      IF (.NOT. ALLOCATED(out_dir)) RETURN
      IF (LEN(out_dir) == 0) RETURN
      IF (rename_for_orbitrans) path = orbitrans_name(path)
      path = out_dir//'/'//path
   END FUNCTION wannierlib_export_file

   !> The names orbitrans has in its OPEN statements. Anything absent from this list --
   !> a second channel, a diagnostic -- keeps its own name and simply moves.
   FUNCTION orbitrans_name(canonical) RESULT(name)
      CHARACTER(LEN=*), INTENT(IN) :: canonical
      CHARACTER(LEN=:), ALLOCATABLE :: name

      SELECT CASE (canonical)
      CASE ('rspauli.1');    name = 'spin'
      CASE ('anglmomrs.1');  name = 'orbital'
      CASE ('rssocmat.1');   name = 'spin_orbit'
      CASE ('WF1_hr.dat');   name = 'hoppings'
      !> It reads `dum, R(1:3)` from this one and our fifth column, the degeneracy, is
      !> simply not in its read list -- list-directed input takes what it asks for.
      CASE ('wig_vectors');  name = 'rpts'
      CASE ('WF1_bmn.dat');  name = 'bmn'
      CASE ('WF1_cmn.dat');  name = 'cmn'
      CASE ('WF1_fmn.dat');  name = 'fmn'
      CASE DEFAULT;          name = TRIM(canonical)
      END SELECT
   END FUNCTION orbitrans_name

END MODULE m_wannierlib_export_paths
