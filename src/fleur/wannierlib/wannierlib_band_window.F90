!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> What the band window selects out of the SCF results: the energy bounds the input
!> left out, and the eigenvalues of the selected bands over the wannierization mesh.
MODULE m_wannierlib_band_window
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   USE m_types_kpts
   USE m_types_input
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   USE m_types_misc
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_default_windows, wannierlib_create_eig, wannierlib_check_bands
CONTAINS

   !> That the bands the wannierization is about to read were actually computed.
   !>
   !> `numbands` in <cutoffs> and `maxBand` in <wannierlib><bands> are two different knobs
   !> with nearly the same name: the first is how many bands FLEUR KEEPS, the second how
   !> many the wannierization READS. Asking for a band that is not there reports nothing --
   !> the read runs past the end of the record and the run dies with
   !>
   !>      Signal 11 detected on PE: 0
   !>      Last timer: IO (read)  /  Timer: wannierlib
   !>
   !> which names neither bands nor windows, and which is insensitive to num_wann, to the
   !> disentanglement and to the operators, so the search goes everywhere except here. On
   !> Si in October 2026 that cost five runs: numbands was left at 0, FLEUR chose seven
   !> bands for the eight valence electrons, and maxBand asked for the eighth.
   !>
   !> Counted from results%neig, which is how many eigenvalues each k-point actually has.
   !> That is the one count free of convention: `numbands` means spinor bands with SOC and
   !> no noco, and per-channel bands under noco, where FLEUR doubles it -- so comparing
   !> against the input value would be right in one case and wrong in the other.
   SUBROUTINE wannierlib_check_bands(this, results, kpts, input, l_spinors)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_results), INTENT(IN) :: results
      TYPE(t_kpts),    INTENT(IN) :: kpts
      TYPE(t_input),   INTENT(IN) :: input
      LOGICAL,         INTENT(IN) :: l_spinors

      INTEGER :: ikpt, jsp, nsp, navail, kworst, jworst

      IF (this%min_band < 1 .OR. this%max_band < this%min_band) CALL juDFT_error( &
         'wannierlib: <bands> has minBand and maxBand the wrong way round or below one', &
         calledby='wannierlib_check_bands')

      IF (.NOT. ALLOCATED(results%eig) .OR. .NOT. ALLOCATED(results%neig)) CALL juDFT_error( &
         'wannierlib: the eigenvalues of the SCF are not available to the wannierization', &
         hint='run the wannierization as a separate step on a converged density', &
         calledby='wannierlib_check_bands')

      !> The scarcest k-point of the wannierization mesh, not the first: the count may
      !> differ between k, and the one that crashes is the one with fewest bands.
      nsp = MERGE(1, input%jspins, l_spinors)
      navail = HUGE(navail)
      kworst = 0; jworst = 0
      DO jsp = 1, nsp
         DO ikpt = 1, kpts%nkptf
            IF (results%neig(kpts%bkp(ikpt), jsp) < navail) THEN
               navail = results%neig(kpts%bkp(ikpt), jsp)
               kworst = ikpt; jworst = jsp
            END IF
         END DO
      END DO

      IF (this%max_band > navail .OR. this%max_band > SIZE(results%eig, 1)) THEN
         WRITE (oUnit, '(/,a)')   'wannierlib: maxBand asks for bands that were not computed'
         WRITE (oUnit, '(a,2i8)') '   bands asked for (min, max) :', this%min_band, this%max_band
         WRITE (oUnit, '(a,i8)')  '   bands available, worst k    :', navail
         WRITE (oUnit, '(a,2i8)') '   at k-point / spin           :', kworst, jworst
         WRITE (oUnit, '(a,i8)')  '   room in results%eig         :', SIZE(results%eig, 1)
         CALL juDFT_error( &
            'wannierlib: maxBand is above the number of bands the calculation has', &
            hint='raise numbands in <cutoffs> to at least maxBand and rerun; '// &
                 'numbands="0" lets FLEUR pick, and it picks what the SCF needs, not what '// &
                 'the wannierization reads', &
            calledby='wannierlib_check_bands')
      END IF
   END SUBROUTINE wannierlib_check_bands

   ! Fill in the energy windows the input did not state, from the range the selected bands
   ! actually span ON THE WANNIERIZATION MESH. That is what Wannier90 falls back to when the
   ! classic route leaves dis_win_min/max out of the .win, and it is the only definition that
   ! cannot cut the manifold.
   !
   ! Measuring the range anywhere else is a silent trap: taken from the SCF mesh, which need
   ! not contain Gamma, the minimum comes out ABOVE the true one, states fall outside the
   ! window at some k, and Wannier90 drops them without a word.
   !
   ! A margin of one microhartree is added on each side so that a band sitting exactly on the
   ! edge is not dropped by rounding.
   !
   ! Both spin channels are scanned when there are two, so a single window holds every band
   ! that will be wannierized in either.
   SUBROUTINE wannierlib_default_windows(this, results, kpts, input, l_spinors)
      TYPE(t_wannierlib_wannierize), INTENT(INOUT) :: this
      TYPE(t_results), INTENT(IN) :: results
      TYPE(t_kpts), INTENT(IN) :: kpts
      TYPE(t_input), INTENT(IN) :: input
      LOGICAL, INTENT(IN) :: l_spinors

      REAL, PARAMETER :: margin = 1.0e-6      ! Hartree
      INTEGER :: ikpt, jsp, nsp
      REAL :: emin, emax

      !> With as many bands as Wannier functions there is no subspace to choose, so there is
      !> no window either: w90 runs no disentanglement and the options are never read. Filling
      !> them in anyway hands it bounds for a step it does not take, and it refuses the run.
      IF (this%num_bands <= this%num_wann) RETURN

      !> Each bound is handled on its own: a window left at its 0.0 default was not given.
      !> Testing the pair instead would break a half-stated window -- give only disWinMin
      !> and disWinMax would stay at zero.
      IF (this%dis_win_min /= 0.0 .AND. this%dis_win_max /= 0.0 .AND. this%dis_froz_min /= 0.0) RETURN

      nsp = MERGE(1, input%jspins, l_spinors)
      emin = HUGE(emin)
      emax = -HUGE(emax)
      DO jsp = 1, nsp
         DO ikpt = 1, kpts%nkptf
            emin = MIN(emin, results%eig(this%min_band, kpts%bkp(ikpt), jsp))
            emax = MAX(emax, results%eig(this%max_band, kpts%bkp(ikpt), jsp))
         END DO
      END DO

      !> An empty or degenerate range means the eigenvalues are not in results%eig at the
      !> point this runs, and the window would come out [0,0]: Wannier90 then reports
      !> "energy window contains fewer states than target WFs" at some k, or worse,
      !> disentangles an empty subspace. Refusing here turns a silent wrong window into
      !> a message that names the remedy.
      IF (.NOT. (emax > emin)) THEN
         !> Say what was scanned before giving up. An empty range means either a band pair
         !> that selects nothing or eigenvalues that are not here yet, and the two need
         !> opposite fixes, so the numbers go out rather than being left to be guessed.
         WRITE (oUnit, '(a)') 'wannierlib: the outer energy window could not be derived'
         WRITE (oUnit, '(a,2i8)')  '   bands scanned (min, max)  :', this%min_band, this%max_band
         WRITE (oUnit, '(a,i8)')   '   k-points on the full mesh :', kpts%nkptf
         WRITE (oUnit, '(a,2es16.6)') '   range over those bands   :', emin, emax
         !> The whole array, not just the slice above: if it holds real eigenvalues while the
         !> slice came out flat, the band range is at fault; if the array itself is flat or
         !> absent, the eigenvalues are not here yet. One line separates the two.
         IF (ALLOCATED(results%eig)) THEN
            WRITE (oUnit, '(a,3i8)')     '   shape of results%eig      :', SHAPE(results%eig)
            WRITE (oUnit, '(a,2es16.6)') '   min/max of the whole array:', &
               MINVAL(results%eig), MAXVAL(results%eig)
         ELSE
            WRITE (oUnit, '(a)')         '   results%eig is NOT allocated'
         END IF
         CALL juDFT_error( &
            'wannierlib: the outer energy window cannot be derived from the bands', &
            hint='state disWinMin and disWinMax explicitly in <disentanglement>; '// &
                 'the numbers behind this are in the out file', &
            calledby='wannierlib_default_windows')
      END IF

      IF (this%dis_win_min == 0.0) this%dis_win_min = emin - margin
      IF (this%dis_win_max == 0.0) this%dis_win_max = emax + margin
      !> The frozen window starts where the outer one does unless the input says
      !> otherwise, which is what the established recipes write by hand.
      IF (this%dis_froz_min == 0.0) this%dis_froz_min = this%dis_win_min
   END SUBROUTINE wannierlib_default_windows

   subroutine wannierlib_create_eig(this, results, kpts, jspin, eig)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_results), INTENT(IN) :: results
      TYPE(t_kpts), INTENT(IN) :: kpts
      INTEGER, INTENT(IN) :: jspin
      REAL, ALLOCATABLE, INTENT(OUT) :: eig(:, :)

      INTEGER :: nkpt
      allocate (eig(this%num_bands, kpts%nkptf))

      DO nkpt = 1, kpts%nkptf
         eig(:this%num_bands, nkpt) = results%eig(this%min_band:this%max_band, kpts%bkp(nkpt), jspin)
      END DO

   END SUBROUTINE wannierlib_create_eig

END MODULE m_wannierlib_band_window
