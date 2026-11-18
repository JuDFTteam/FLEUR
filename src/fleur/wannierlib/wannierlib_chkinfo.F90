!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_wannierlib_chkinfo
  !> The witness written beside a Wannier90 checkpoint, and the comparison that reads it.
  !>
  !> This lives apart from the Wannier90 adapter on purpose. The adapter's job is to put
  !> FLEUR's quantities into the shapes the library expects and to hand back what it
  !> produced; the .chkinfo is neither written nor read by Wannier90 -- it is ours, and what
  !> handles it is a text parser and a set of comparisons, which is a different job. The
  !> seedname arrives as an argument for the same reason, and that is the whole of what the
  !> two had in common.
  !>
  !> It also means this file does not need CPP_WANNLIB_API: the checks compile and can be
  !> exercised whether or not the library is there.
  USE m_constants, ONLY : oUnit
  USE m_juDFT
  USE m_judft_string, ONLY : int2str
  USE m_types_wannierlib
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: wl_chkinfo_write, wl_chkinfo_check

  !> How far apart two window energies may be and still be the same window. It is not a
  !> physical tolerance: both numbers come out of the same parser reading the same inp.xml
  !> text, so they agree bit for bit unless that text changed, and this is margin for the
  !> last one of them. Anything anyone actually edits moves a bound by millihartree.
  REAL, PARAMETER :: WANNIERLIB_CHKINFO_ETOL = 1.0e-10

  !> The fields of one projection, in the order they are written and read back. They are
  !> named so that a refusal can say which one moved rather than only that something did.
  CHARACTER(LEN=6), PARAMETER :: WL_PROJ_I(6) = &
    [CHARACTER(LEN=6) :: 'ntype', 'atom', 'l', 'm', 'spin', 'rwf']
  CHARACTER(LEN=6), PARAMETER :: WL_PROJ_R(11) = &
    [CHARACTER(LEN=6) :: 'alpha', 'beta', 'gamma', 'zona', 'regio', 'j', 'mj', 'weight', &
     'shiftX', 'shiftY', 'shiftZ']

CONTAINS

  !> The name of the witness: the checkpoint's own stem with another extension, so the two
  !> travel together and the per-channel suffix the seedname already carries is inherited
  !> for free. Not '.window', which would sit one letter away from Wannier90's '.win' input.
  FUNCTION wl_chkinfo_file(seedname) RESULT(f)
    CHARACTER(LEN=*), INTENT(IN) :: seedname
    CHARACTER(LEN=:), ALLOCATABLE :: f
    f = TRIM(seedname)//'.chkinfo'
  END FUNCTION wl_chkinfo_file

  !> What a Wannier90 checkpoint cannot carry, written beside it.
  !>
  !> The .chk holds lwindow, which is the CONSEQUENCE of the energy bounds, and imposes it
  !> when it is loaded -- the bounds themselves are not in the file. So three edits survive
  !> all ten of Wannier90's checks and do not do what they look like they do: disWinMin and
  !> disWinMax have no effect at all, disFrozMin and disFrozMax are ignored by a minimisation
  !> that does not redo the disentanglement, and minBand with maxBand shifted together leave
  !> numBands unchanged and hand the run a gauge belonging to other bands. This file is what
  !> lets the next run refuse rather than ignore in silence.
  !>
  !> The projections are in here for the same reason, and they are the fourth such edit: a
  !> restart loads the gauge and skips w90_project_overlap, so an edited <wannierproj> is
  !> simply never used -- measured as a bit-for-bit identical result. They are written as
  !> this%proj stands HERE, which after init is one entry per Wannier function rather than
  !> one per <wannierproj>, so what is compared is the list that was really used.
  SUBROUTINE wl_chkinfo_write(this, route, seedname, irank)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: route
    CHARACTER(LEN=*), INTENT(IN) :: seedname
    INTEGER, INTENT(IN) :: irank
    INTEGER :: iu, io, ip, np

    !> Rank zero alone. Every rank calls the checkpoint writer, because that one gathers the
    !> distributed matrices itself, so without this guard every rank opens the same path for
    !> replacement at the same moment.
    IF (irank /= 0) RETURN

    OPEN (NEWUNIT=iu, FILE=wl_chkinfo_file(seedname), FORM='formatted', STATUS='replace', IOSTAT=io)
    IF (io /= 0) CALL juDFT_error('wannierlib: cannot write '//wl_chkinfo_file(seedname), &
                                  calledby='wl_chkinfo_write')
    WRITE (iu, '(a)')         '# wannierlib: what the checkpoint beside this file cannot carry.'
    WRITE (iu, '(a)')         '# Read back when that checkpoint is read. Energies in Hartree.'
    WRITE (iu, '(a,i0)')      'version    ', 2
    WRITE (iu, '(a,a)')       'route      ', TRIM(route)
    WRITE (iu, '(a,i0)')      'numBands   ', this%num_bands
    WRITE (iu, '(a,i0)')      'minBand    ', this%min_band
    WRITE (iu, '(a,i0)')      'maxBand    ', this%max_band
    WRITE (iu, '(a,i0)')      'numWann    ', this%num_wann
    WRITE (iu, '(a,es22.14)') 'disWinMin  ', this%dis_win_min
    WRITE (iu, '(a,es22.14)') 'disWinMax  ', this%dis_win_max
    WRITE (iu, '(a,es22.14)') 'disFrozMin ', this%dis_froz_min
    WRITE (iu, '(a,es22.14)') 'disFrozMax ', this%dis_froz_max

    !> One record per projection, each carrying every field of it. The species is quoted so
    !> that an empty one survives the list-directed read that takes the record back.
    np = 0
    IF (ALLOCATED(this%proj)) np = SIZE(this%proj)
    WRITE (iu, '(a,i0)') 'nproj      ', np
    DO ip = 1, np
      WRITE (iu, '(a,i0,6(1x,i0),11(1x,es22.14),1x,a)') 'proj ', ip, &
        this%proj(ip)%ntype, this%proj(ip)%atom, this%proj(ip)%l, this%proj(ip)%m, &
        this%proj(ip)%spin, this%proj(ip)%rwf, &
        this%proj(ip)%alpha, this%proj(ip)%beta, this%proj(ip)%gamma, this%proj(ip)%zona, &
        this%proj(ip)%regio, this%proj(ip)%j, this%proj(ip)%mj, this%proj(ip)%weight, &
        this%proj(ip)%shift(1), this%proj(ip)%shift(2), this%proj(ip)%shift(3), &
        "'"//TRIM(this%proj(ip)%species)//"'"
    END DO
    CLOSE (iu)
  END SUBROUTINE wl_chkinfo_write

  !> Judge the witness. Every rank reads it and reaches its own verdict, rather than rank
  !> zero deciding and broadcasting: the verdict is a function of a few hundred bytes, so
  !> all ranks reach the same one and abort together, whereas a verdict taken on rank zero
  !> alone would need a broadcast and without one leaves the others in the next collective.
  SUBROUTINE wl_chkinfo_check(this, route, seedname)
    TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
    CHARACTER(LEN=*), INTENT(IN) :: route
    CHARACTER(LEN=*), INTENT(IN) :: seedname
    INTEGER :: iu, io, ver, i_f, ip, k, np_f, np_now, idx
    INTEGER :: pi_f(6), pi_now(6)
    REAL :: r_f, pr_f(11), pr_now(11)
    CHARACTER(LEN=32) :: a_f, sp_f
    CHARACTER(LEN=512) :: rec
    LOGICAL :: l_there

    INQUIRE (FILE=wl_chkinfo_file(seedname), EXIST=l_there)
    !> A checkpoint from before this file existed. Its gauge is still usable and all ten of
    !> Wannier90's checks still run on it, so the run goes on -- but the one thing those
    !> checks cannot see is now going unchecked, and that is worth a line rather than being
    !> the quietest part of the run.
    IF (.NOT. l_there) THEN
      WRITE (oUnit, '(a)') ' wannierlib: no '//wl_chkinfo_file(seedname)//' beside the checkpoint: '// &
                           'THE WINDOW AND THE BAND INDICES WERE NOT CHECKED'
      RETURN
    END IF

    OPEN (NEWUNIT=iu, FILE=wl_chkinfo_file(seedname), FORM='formatted', STATUS='old', IOSTAT=io)
    IF (io /= 0) CALL juDFT_error('wannierlib: cannot read '//wl_chkinfo_file(seedname), &
                                  calledby='wl_chkinfo_check')
    READ (iu, '(a)', IOSTAT=io)      ! the two comment lines
    READ (iu, '(a)', IOSTAT=io)

    CALL rd_i('version', ver)
    !> Version 2 added the projections. A version 1 file still answers for everything it does
    !> carry, so it is read rather than refused -- and the one check it cannot serve is said
    !> out loud, the same way the operator file treats its own version 1.
    IF (ver < 1 .OR. ver > 2) &
      CALL juDFT_error('wannierlib: '//wl_chkinfo_file(seedname)//' is layout version '// &
                       int2str(ver)//' and this build reads 1 and 2', &
                       hint='delete it to go on without the check, or write the '// &
                            'checkpoint again', calledby='wl_chkinfo_check')

    CALL rd_a('route', a_f);      CALL same_a('route',      a_f, route)
    CALL rd_i('numBands', i_f);   CALL same_i('numBands',   i_f, this%num_bands)
    CALL rd_i('minBand', i_f);    CALL same_i('minBand',    i_f, this%min_band)
    CALL rd_i('maxBand', i_f);    CALL same_i('maxBand',    i_f, this%max_band)
    CALL rd_i('numWann', i_f);    CALL same_i('numWann',    i_f, this%num_wann)
    CALL rd_r('disWinMin', r_f);  CALL same_r('disWinMin',  r_f, this%dis_win_min, ' Htr')
    CALL rd_r('disWinMax', r_f);  CALL same_r('disWinMax',  r_f, this%dis_win_max, ' Htr')
    CALL rd_r('disFrozMin', r_f); CALL same_r('disFrozMin', r_f, this%dis_froz_min, ' Htr')
    CALL rd_r('disFrozMax', r_f); CALL same_r('disFrozMax', r_f, this%dis_froz_max, ' Htr')

    !> The projections, record by record. The count goes first so that a list of another
    !> length is refused as a count rather than as a read running off the end of the file.
    np_now = 0
    IF (ALLOCATED(this%proj)) np_now = SIZE(this%proj)
    IF (ver < 2) THEN
      WRITE (oUnit, '(a)') ' wannierlib: '//wl_chkinfo_file(seedname)//' is layout version 1 and '// &
                           'carries no projections: THE PROJECTIONS WERE NOT CHECKED'
      CLOSE (iu)
      RETURN
    END IF
    CALL rd_i('nproj', np_f)
    CALL same_i('nproj', np_f, np_now)
    DO ip = 1, np_now
      CALL rd_key('proj', rec)
      READ (rec, *, IOSTAT=io) idx, pi_f, pr_f, sp_f
      IF (io /= 0) CALL juDFT_error('wannierlib: projection '//int2str(ip)//' of '// &
                                    wl_chkinfo_file(seedname)//' is not readable', &
                                    calledby='wl_chkinfo_check')
      CALL same_i('the index of projection '//int2str(ip), idx, ip)
      pi_now = [this%proj(ip)%ntype, this%proj(ip)%atom, this%proj(ip)%l, &
                this%proj(ip)%m, this%proj(ip)%spin, this%proj(ip)%rwf]
      pr_now = [this%proj(ip)%alpha, this%proj(ip)%beta, this%proj(ip)%gamma, &
                this%proj(ip)%zona, this%proj(ip)%regio, this%proj(ip)%j, &
                this%proj(ip)%mj, this%proj(ip)%weight, this%proj(ip)%shift(1), &
                this%proj(ip)%shift(2), this%proj(ip)%shift(3)]
      DO k = 1, 6
        CALL same_i('projection '//int2str(ip)//' '//TRIM(WL_PROJ_I(k)), pi_f(k), pi_now(k))
      END DO
      DO k = 1, 11
        !> No unit: the eleven fields of a projection do not share one, and the field
        !> name is what says which it is.
        CALL same_r('projection '//int2str(ip)//' '//TRIM(WL_PROJ_R(k)), pr_f(k), pr_now(k), '')
      END DO
      CALL same_s('the species of projection '//int2str(ip), sp_f, this%proj(ip)%species)
    END DO
    CLOSE (iu)

    WRITE (oUnit, '(a)') ' wannierlib: the window, the band indices, the route and '// &
                         int2str(np_now)//' projections agree with the checkpoint'

  CONTAINS

    !> Each read states which key it expects, so a file written by something else, or cut
    !> short, is refused where it goes wrong instead of being read as the next field.
    SUBROUTINE rd_key(want, line)
      CHARACTER(LEN=*), INTENT(IN) :: want
      CHARACTER(LEN=*), INTENT(OUT) :: line
      INTEGER :: jo
      READ (iu, '(a)', IOSTAT=jo) line
      IF (jo /= 0 .OR. INDEX(ADJUSTL(line), want) /= 1) CALL juDFT_error( &
        'wannierlib: '//wl_chkinfo_file(seedname)//' has no '//want//' where one is due', &
        hint='it was written by another version of this code, or truncated: delete it to '// &
             'go on without the check, or write the checkpoint again', &
        calledby='wl_chkinfo_check')
      line = ADJUSTL(line)
      line = line(LEN(want) + 1:)
    END SUBROUTINE rd_key

    SUBROUTINE rd_i(want, val)
      CHARACTER(LEN=*), INTENT(IN) :: want
      INTEGER, INTENT(OUT) :: val
      CHARACTER(LEN=200) :: line
      CALL rd_key(want, line)
      READ (line, *) val
    END SUBROUTINE rd_i

    SUBROUTINE rd_r(want, val)
      CHARACTER(LEN=*), INTENT(IN) :: want
      REAL, INTENT(OUT) :: val
      CHARACTER(LEN=200) :: line
      CALL rd_key(want, line)
      READ (line, *) val
    END SUBROUTINE rd_r

    SUBROUTINE rd_a(want, val)
      CHARACTER(LEN=*), INTENT(IN) :: want
      CHARACTER(LEN=*), INTENT(OUT) :: val
      CHARACTER(LEN=200) :: line
      CALL rd_key(want, line)
      val = TRIM(ADJUSTL(line))
    END SUBROUTINE rd_a

    !> The two values always go in the message. One that says only that something does not
    !> match sends the reader back to the file to work out what.
    SUBROUTINE same_i(what, in_file, in_run)
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN) :: in_file, in_run
      IF (in_file == in_run) RETURN
      CALL juDFT_error('wannierlib: the checkpoint was written with '//what//' = '// &
                       int2str(in_file)//', this run has '//int2str(in_run), &
                       hint=wl_chkinfo_hint(), calledby='wl_chkinfo_check')
    END SUBROUTINE same_i

    SUBROUTINE same_a(what, in_file, in_run)
      CHARACTER(LEN=*), INTENT(IN) :: what, in_file, in_run
      IF (TRIM(in_file) == TRIM(in_run)) RETURN
      CALL juDFT_error('wannierlib: the checkpoint was written on the '//TRIM(in_file)// &
                       ' route and this is the '//TRIM(in_run)//' one', &
                       hint=wl_chkinfo_hint(), calledby='wl_chkinfo_check')
    END SUBROUTINE same_a

    !> For anything else held as text, which says what it is rather than reading as a route.
    SUBROUTINE same_s(what, in_file, in_run)
      CHARACTER(LEN=*), INTENT(IN) :: what, in_file, in_run
      IF (TRIM(in_file) == TRIM(in_run)) RETURN
      CALL juDFT_error('wannierlib: the checkpoint has '//TRIM(in_file)//' as '//what// &
                       ', this run has '//TRIM(in_run), &
                       hint=wl_chkinfo_hint(), calledby='wl_chkinfo_check')
    END SUBROUTINE same_s

    SUBROUTINE same_r(what, in_file, in_run, unit)
      CHARACTER(LEN=*), INTENT(IN) :: what
      REAL, INTENT(IN) :: in_file, in_run
      !> Stated by the caller, because these are not all energies: the window bounds are in
      !> Hartree, a projection's shift is a position and its angles are angles. A message
      !> that stamps Htr on all of them names the wrong unit to whoever has to act on it.
      CHARACTER(LEN=*), INTENT(IN) :: unit
      CHARACTER(LEN=220) :: msg
      IF (ABS(in_file - in_run) <= WANNIERLIB_CHKINFO_ETOL) RETURN
      WRITE (msg, '(a,a,es20.12,a,es20.12,a)') 'wannierlib: the checkpoint was written with ', &
        what//' = ', in_file, unit//', this run has ', in_run, unit
      CALL juDFT_error(TRIM(msg), hint=wl_chkinfo_hint(), calledby='wl_chkinfo_check')
    END SUBROUTINE same_r

    !> One hint for all of them, because the reason is the same every time: the gauge in the
    !> file was built for those bands inside that window, and a restart reuses it rather
    !> than rebuilding it. The minimisation settings are deliberately not in the witness --
    !> numIter, convTol, convWindow and precond stay free, which is the whole point of
    !> resuming.
    FUNCTION wl_chkinfo_hint() RESULT(h)
      CHARACTER(LEN=:), ALLOCATABLE :: h
      h = 'a restart reuses the gauge in the checkpoint and cannot rebuild it for another '// &
          'window or another set of bands: put the value back, or drop restart= and '// &
          'wannierise from the start'
    END FUNCTION wl_chkinfo_hint

  END SUBROUTINE wl_chkinfo_check

END MODULE m_wannierlib_chkinfo
