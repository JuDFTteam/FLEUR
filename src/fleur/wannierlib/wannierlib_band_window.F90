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
   USE m_types, ONLY: t_results
   USE m_types_kpts
   USE m_types_input
   USE m_types_wannierlib, ONLY: t_wannierlib_wannierize
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_default_windows, wannierlib_create_eig
CONTAINS

   ! Fill in the energy windows the input did not state, from the range the selected bands
   ! actually span ON THE WANNIERIZATION MESH. That is what Wannier90 falls back to when the
   ! classic route leaves dis_win_min/max out of the .win, and it is the only definition that
   ! cannot cut the manifold.
   !
   ! Measuring the range anywhere else is a silent trap: taken from the SCF mesh, which need
   ! not contain Gamma, the minimum comes out ABOVE the true one, states fall outside the
   ! window at some k, and Wannier90 drops them without a word -- in Pt that cost 7 of 512
   ! k-points two of their 36 states.
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
