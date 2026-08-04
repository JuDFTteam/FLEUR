!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  Output-domain plumbing for the matrix-element interpolation.
!>
!>  An <interpolation> block may declare several output domains (<path>/<plane>/<grid>);
!>  each one interpolates the requested operators on its own k-set and gets its output
!>  files suffixed so the domains do not overwrite each other. This module owns the
!>  k-set generation (kpts_interpol) and the post-run file renaming. Rank-0 file IO only.
!>
!>  Extracted verbatim from m_wannierlib_w90_adapter: it is matrix-element output
!>  bookkeeping and has nothing to do with the Wannier90 library API.
MODULE m_melem_domains
  USE m_juDFT
  USE m_types_wannierlib
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: melem_write_domain_kpts, melem_rename_domain_outputs, melem_shell
CONTAINS

  ! Write kpts_interpol for a generated output domain (plane/grid). For an explicit
  ! <path file="..">, copy that file to kpts_interpol; legacy/default path uses the
  ! existing kpts_interpol as-is. Rank-0 file I/O only.
  SUBROUTINE melem_write_domain_kpts(this, kind)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: kind
    INTEGER :: i, j, k, iu, np
    REAL :: t1, t2, kf(3)
    LOGICAL :: lex2
    np = 0
    SELECT CASE (TRIM(kind))
    CASE ('path', 'legacy')
      ! restore the user's original kpts_interpol if a generated (plane/grid) domain overwrote it
      INQUIRE(file='.kpts_interpol_userbak', exist=lex2)
      IF (lex2) CALL melem_shell('cp -f .kpts_interpol_userbak kpts_interpol')
      ! <path listName=".."> : write the (optionally subdivided) named list as kpts_interpol.
      ! Takes precedence over @file; path_kpts was filled from kpts.xml in read_xml.
      IF (TRIM(kind) == 'path' .AND. this%path_kset%nkpt > 0) THEN
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') this%path_kset%nkpt
        DO i = 1, this%path_kset%nkpt
          WRITE(iu,'(3(f18.12,1x))') this%path_kset%bk(:, i)
        END DO
        CLOSE(iu)
        RETURN
      END IF
      ! explicit <path file="..">: use that file as the k-list (default 'kpts_interpol' -> no-op)
      IF (TRIM(kind) == 'path' .AND. TRIM(this%path_file) /= 'kpts_interpol') THEN
        INQUIRE(file=TRIM(this%path_file), exist=lex2)
        IF (.NOT. lex2) CALL juDFT_error('wannierlib: <path>/@file "'//TRIM(this%path_file)//'" not found', &
                                         calledby='melem_write_domain_kpts')
        CALL melem_shell('cp -f '//TRIM(this%path_file)//' kpts_interpol')
      END IF
      RETURN
    CASE ('plane')
      IF (this%plane_kset%nkpt > 0) THEN     ! <plane listName=".."/>: use the named list as-is
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') this%plane_kset%nkpt
        DO i = 1, this%plane_kset%nkpt
          WRITE(iu,'(3(f18.12,1x))') this%plane_kset%bk(:, i)
        END DO
        CLOSE(iu)
      END IF
    CASE ('grid')
      IF (this%grid_kset%nkpt > 0) THEN      ! <grid listName=".."/>: use the named list as-is
        OPEN(newunit=iu, file='kpts_interpol', status='replace')
        WRITE(iu,'(i0)') this%grid_kset%nkpt
        DO i = 1, this%grid_kset%nkpt
          WRITE(iu,'(3(f18.12,1x))') this%grid_kset%bk(:, i)
        END DO
        CLOSE(iu)
      END IF
    END SELECT
  END SUBROUTINE melem_write_domain_kpts

  ! Rename this domain's operator output files bands_wann_<x>.dat -> bands_wann_<x><suffix>.dat
  SUBROUTINE melem_rename_domain_outputs(this, suffix)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    INTEGER :: iop
    DO iop = 1, this%n_ops
      SELECT CASE (TRIM(this%op_name(iop)))
      CASE ('hamiltonian')
        CALL ren('bands_wann_interpol', suffix)
        CALL ren('bands_wann_interpol_ev', suffix)
      CASE ('spin')
        IF (this%op_total(iop) == 1)   CALL ren('bands_wann_spin', suffix)
      CASE ('orbital')
        IF (this%op_total(iop) == 1)   CALL ren('bands_wann_orbmom', suffix)
      CASE ('soc')
        CALL ren('bands_wann_soc', suffix)
      CASE ('velocity')
        CALL ren('bands_wann_velocity', suffix)
        CALL ren('bands_wann_berrycurv', suffix)
      CASE ('spinCurrent')
        CALL ren('bands_wann_spincurrent', suffix)
      CASE ('orbitalCurrent')
        CALL ren('bands_wann_orbcurrent', suffix)
      CASE ('eigenstates')
        CALL ren('bands_wann_eigenstates', suffix)
      END SELECT
    END DO
  CONTAINS
    SUBROUTINE ren(base, suf)
      CHARACTER(LEN=*), INTENT(IN) :: base, suf
      LOGICAL :: lexr
      INQUIRE(file=TRIM(base)//'.dat', exist=lexr)
      IF (lexr) CALL melem_shell('mv -f '//TRIM(base)//'.dat '//TRIM(base)//TRIM(suf)//'.dat')
    END SUBROUTINE ren
  END SUBROUTINE melem_rename_domain_outputs

  ! Run a shell command (synchronous) and abort with a clear message if it fails,
  ! so a failed cp/mv in the domain file-shuffling never passes silently.
  SUBROUTINE melem_shell(cmd)
    CHARACTER(LEN=*), INTENT(IN) :: cmd
    INTEGER :: cs, es
    cs = 0; es = 0
    CALL EXECUTE_COMMAND_LINE(cmd, wait=.TRUE., cmdstat=cs, exitstat=es)
    IF (cs /= 0 .OR. es /= 0) CALL juDFT_error('wannierlib: shell command failed: '//TRIM(cmd), &
                                               calledby='melem_shell')
  END SUBROUTINE melem_shell

END MODULE m_melem_domains
