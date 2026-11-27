!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> The one place that knows which of the selected bands lie inside the outer window.
!>
!> Wannier90 does not keep the band index. Its disentanglement works on the states inside
!> the outer window and PACKS them: the rows of u_opt run over those states in increasing
!> band order, rows past their count are zeroed (disentangle.F90:287), and the overlaps it
!> was handed are packed to match in place (internal_slim_m, disentangle.F90:206). That
!> packing is an internal convention of the library and nothing in the returned arrays
!> says it happened.
!>
!> Everything on the FLEUR side indexes by band: the Bloch overlaps, the matrix elements,
!> the eigenvalues. Contracting a packed u_opt against any of them pairs row 1 with the
!> lowest SELECTED band when it means the lowest band IN THE WINDOW, so the two differ by
!> however many bands the window cuts off below -- and agree exactly when it cuts off
!> none, which is what the default window does. The result is a defect that stays hidden
!> through every ordinary run and appears the first time someone narrows the window by
!> hand.
!>
!> So the gauge is unpacked once, here, as soon as it comes back from the library and
!> before anything reads it. From that point on u_opt(b, :, k) means band b of the
!> selected set, with zero rows for the bands the window excludes, and every consumer --
!> the operators, the interpolation, the overlap cross-check, the electron-phonon
!> interpolation that reads results%U_dis -- is correct without knowing any of this.
MODULE m_wgauge_window
   USE m_juDFT
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wgauge_unpack_gauge

CONTAINS

   !> eig is in Hartree and holds the SELECTED bands, the same slice that was handed to
   !> the library, so a row of eig and a row of u_opt are the same band before the
   !> packing and have to be made so again after it.
   SUBROUTINE wgauge_unpack_gauge(u_opt, eig, win_min, win_max, l_froz_proj)
      COMPLEX, INTENT(INOUT) :: u_opt(:, :, :)   !< (nb,nw,nk) packed in, by band out
      REAL,    INTENT(IN)    :: eig(:, :)        !< (nb,nk) Hartree, the selected bands
      REAL,    INTENT(IN)    :: win_min, win_max !< the outer window, as resolved
      LOGICAL, INTENT(IN)    :: l_froz_proj      !< projectability route: cannot be mapped

      LOGICAL, ALLOCATABLE :: inwin(:)
      INTEGER :: nb, nw, nk, k, b, p, nin

      nb = SIZE(u_opt, 1); nw = SIZE(u_opt, 2); nk = SIZE(u_opt, 3)

      !> Without disentanglement there is no window and no packing: u_opt is the
      !> projection itself, already one row per band.
      IF (nb <= nw .OR. nk == 0) RETURN

      !> The projectability route selects states by how well the trial orbitals span
      !> them, not by energy, and keeps them through an index list the library does not
      !> hand back. The mapping below would then be wrong rather than missing, so it is
      !> refused instead of guessed.
      IF (l_froz_proj) &
         CALL juDFT_error("wannierlib: the Wannier gauge cannot be put back on the band "// &
                          "index for a projectability disentanglement", &
                          hint='disFrozProj="T" selects states by projectability, and which '// &
                          'ones it kept is not among what Wannier90 returns', &
                          calledby="wgauge_unpack_gauge")

      ALLOCATE (inwin(nb))
      DO k = 1, nk
         inwin = (eig(:, k) >= win_min .AND. eig(:, k) <= win_max)
         nin = COUNT(inwin)

         !> Fewer states in the window than Wannier functions is not something that could
         !> have produced the gauge in hand, so reaching here means the window used to
         !> build it is not the window being read now.
         IF (nin < nw) &
            CALL juDFT_error("wannierlib: the outer window holds fewer states than there "// &
                             "are Wannier functions", hint='the window the gauge was built '// &
                             'with and the one read here are not the same', &
                             calledby="wgauge_unpack_gauge")

         IF (nin == nb) CYCLE   ! the window cuts nothing off: the packing is the identity

         !> THE PREMISE, MEASURED RATHER THAN ASSUMED. Wannier90 zeroes the rows past the
         !> state count it packed into, so those rows are where its count ends. If ours
         !> ended elsewhere -- a band exactly on the window edge rounded the other way,
         !> say -- there would be something non-zero here, and the rows below would be
         !> moved to the wrong bands in silence.
         IF (ANY(u_opt(nin + 1:, :, k) /= CMPLX(0.0, 0.0))) &
            CALL juDFT_error("wannierlib: the gauge does not end where the outer window "// &
                             "says it should", hint='Wannier90 packed a different number '// &
                             'of states into u_opt than the window counts here', &
                             calledby="wgauge_unpack_gauge")

         !> Downwards, so a row is read before it is written: the packed index of a band
         !> is never above the band itself, and the row it is read from is only reached
         !> later in this direction.
         p = nin
         DO b = nb, 1, -1
            IF (inwin(b)) THEN
               u_opt(b, :, k) = u_opt(p, :, k)
               p = p - 1
            END IF
         END DO
         DO b = 1, nb
            IF (.NOT. inwin(b)) u_opt(b, :, k) = CMPLX(0.0, 0.0)
         END DO
      END DO
   END SUBROUTINE wgauge_unpack_gauge

END MODULE m_wgauge_window
