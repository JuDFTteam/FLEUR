!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

!> Disentanglement that keeps the same number of states from each spin channel at every
!> k-point, for a non-collinear run whose moments are collinear in the global frame.
!>
!> WHAT THIS FIXES, and how small it turned out to be. Without spin-orbit coupling the two
!> spin channels are decoupled: the cross-spin block of M(k,b) is zero to machine precision,
!> so Z(k) = sum_b w_b M U U^dag M^dag is already block-diagonal and its eigenvectors already
!> come out spin-pure. The ordinary disentanglement then keeps the num_wann largest
!> eigenvalues GLOBALLY, and at an occasional k that means ten states of one channel and
!> eight of the other. Measured on bcc Fe with a 4x4x4 mesh: one k-point out of 64.
!>
!> That single point is not a small error. The Wannier index is the same at every k, so if
!> one k cannot be labelled nine-and-nine the labelling is inconsistent globally and the
!> gauge is forced to mix spin EVERYWHERE, not only there. Measured damage: the on-site
!> weight of S(R) falls to 52 %, which leaves the spin operator not interpolable at all.
!> Forcing the count costs 0.0000 % in Omega_I -- free and balanced both reach 11.441731 --
!> because the balanced subspace is an equally good optimum, just not the one a global sort
!> happens to pick.
!>
!> WHY HERE AND NOT IN WANNIER90. The library passes u_opt in and out, so a caller may fill
!> it and skip w90_disentangle; the wannierisation that follows neither knows nor cares who
!> chose the subspace. Wannier90 is read, never patched: a fork would not survive its own
!> upstream.
!>
!> SCOPE, narrow and stated. This assumes spin is a good quantum number, true for a
!> non-collinear run only when the moments are parallel or antiparallel in the global frame
!> -- the usual case of pointing a ferromagnet somewhere other than z. With spin-orbit
!> coupling, or with a genuinely non-collinear texture, no direction commutes with the
!> Hamiltonian and the labels are not exact. Every entry point below refuses in that case
!> rather than labelling on a coin toss, and the caller stops the run: asking for
!> spinBalanced and silently getting the ordinary disentanglement would hide exactly the
!> thing worth knowing.
MODULE m_wannierlib_disentangle_spin
   USE m_juDFT
   USE m_constants, ONLY: oUnit
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_disentangle_spin

   REAL, PARAMETER :: PURITY_MIN = 0.90     !> below this, spin is not a good quantum number
   INTEGER, PARAMETER :: MAXIT = 1000
   REAL, PARAMETER :: MIXING = 0.5, CONVTOL = 1.0e-10

CONTAINS

   !> Fills u_opt with a spin-balanced subspace. ok = .FALSE. means nothing was written and
   !> the caller should run the ordinary disentanglement instead.
   SUBROUTINE wannierlib_disentangle_spin(win_min, win_max, froz_min, froz_max, &
                                          nntot, nnkp, wb, eig_ev, mmn, amn, proj_spin, &
                                          u_opt, ok, s0)
      COMPLEX, INTENT(IN)    :: mmn(:, :, :, :)     !> (nb, nb, nntot, nk)
      COMPLEX, INTENT(IN)    :: amn(:, :, :)        !> (nb, nw, nk)
      COMPLEX, INTENT(INOUT) :: u_opt(:, :, :)      !> (nb, nw, nk)
      REAL,    INTENT(IN)    :: eig_ev(:, :)        !> (nb, nk), eV
      REAL,    INTENT(IN)    :: wb(:)               !> (nntot)
      INTEGER, INTENT(IN)    :: nnkp(:, :)          !> (nk, nntot)
      INTEGER, INTENT(IN)    :: nntot
      INTEGER, INTENT(IN)    :: proj_spin(:)     !> (nw) +1/-1: the channel of each projection
      REAL,    INTENT(IN)    :: win_min, win_max, froz_min, froz_max
      LOGICAL, INTENT(OUT)   :: ok
      !> The spin operator in the Bloch basis, (nb,nb,3,nk). When it is here the bands are
      !> labelled by their own spin; without it there is only the orbital criterion, which is
      !> the same thing in a ferromagnet and wrong in an antiferromagnet.
      COMPLEX, INTENT(IN), OPTIONAL :: s0(:, :, :, :)

      INTEGER :: nb, nw, nk, ik, ib, ich, nper
      INTEGER, ALLOCATABLE :: lab(:, :), col_ch(:)
      LOGICAL, ALLOCATABLE :: allow(:, :), froz(:, :)
      INTEGER :: na, nf

      nb = SIZE(amn, 1); nw = SIZE(amn, 2); nk = SIZE(amn, 3)
      ok = .FALSE.
      IF (MOD(nw, 2) /= 0 .OR. nb <= nw) RETURN
      nper = nw/2

      ALLOCATE (col_ch(nw), lab(nb, nk))
      !> The channel of each projection is read from proj_spin, not inferred. It cannot be
      !> inferred from the .amn: that array's row index is the BAND index, ordered by energy,
      !> and every band is a whole spinor, so those rows carry no spin structure to read.
      IF (ANY(proj_spin(1:nw) == 0)) THEN
         WRITE (oUnit, '(a)') 'wannierlib: the spin-balanced disentanglement needs '// &
            'projections with a definite spin; spinBalanced does not apply here.'
         RETURN
      END IF
      col_ch = MERGE(1, 2, proj_spin(1:nw) > 0)
      IF (COUNT(col_ch == 1) /= nper) THEN
         WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: spin-balanced disentanglement needs the '// &
            'same number of projections per channel; got ', COUNT(col_ch == 1), ' and ', &
            nw - COUNT(col_ch == 1), '. Falling back.'
         RETURN
      END IF
      IF (PRESENT(s0)) THEN
         CALL bands_by_spin(s0, lab, ok)
      ELSE
         CALL bands_by_channel(amn, col_ch, lab, ok)
      END IF
      IF (.NOT. ok) RETURN
      !> The labels say which channel each BAND belongs to; col_ch says which channel each
      !> Wannier column belongs to. They have to agree on which channel is called 1, and
      !> nothing so far forces that: the spin criterion calls 1 the one along +m, the columns
      !> call 1 the one with proj_spin > 0. If they disagree the two are swapped here rather
      !> than left to produce a subspace that is orthogonal to its own projections.
      IF (PRESENT(s0)) CALL match_labels(amn, col_ch, lab)

      ALLOCATE (allow(nb, nk), froz(nb, nk))
      DO ik = 1, nk
         DO ib = 1, nb
            allow(ib, ik) = eig_ev(ib, ik) >= win_min .AND. eig_ev(ib, ik) <= win_max
            froz(ib, ik) = allow(ib, ik) .AND. eig_ev(ib, ik) >= froz_min &
                           .AND. eig_ev(ib, ik) <= froz_max
         END DO
      END DO

      !> Two ways this cannot be done, both reported rather than worked around: a channel
      !> with more frozen states than it has Wannier functions, and a channel with too few
      !> states inside the outer window to fill it. The second is the one that bites, because
      !> a window chosen by looking at the Hamiltonian alone knows nothing about the split.
      DO ik = 1, nk
         DO ich = 1, 2
            na = COUNT(allow(:, ik) .AND. lab(:, ik) == ich)
            nf = COUNT(froz(:, ik) .AND. lab(:, ik) == ich)
            IF (nf > nper .OR. na < nper) THEN
               WRITE (oUnit, '(a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
                  'wannierlib: spin-balanced disentanglement impossible at k = ', ik, &
                  ': channel ', ich, ' has ', na, ' states in the window and ', nf, &
                  ' frozen, for ', nper, ' Wannier functions.'
               WRITE (oUnit, '(a)') '  widen disWinMax until every channel has enough, '// &
                  'or drop spinBalanced.'
               ok = .FALSE.; RETURN
            END IF
         END DO
      END DO

      CALL smv_balanced(mmn, amn, lab, col_ch, allow, froz, nnkp, wb, nntot, nper, u_opt)
      WRITE (oUnit, '(a,i0,a,i0,a)') 'wannierlib: spin-balanced disentanglement, ', nk, &
         ' k-points, ', nper, ' states per channel at every one'
      ok = .TRUE.
   END SUBROUTINE wannierlib_disentangle_spin

   !> The channel of each band, from which projections it overlaps. Refuses when the
   !> separation is poor: a label assigned there would be a coin toss dressed as data.
   SUBROUTINE bands_by_channel(amn, col_ch, lab, ok)
      COMPLEX, INTENT(IN) :: amn(:, :, :)
      INTEGER, INTENT(IN) :: col_ch(:)
      INTEGER, INTENT(OUT) :: lab(:, :)
      LOGICAL, INTENT(OUT) :: ok
      INTEGER :: ib, iw, ik
      REAL :: w1, w2, worst
      worst = 1.0
      DO ik = 1, SIZE(amn, 3)
         DO ib = 1, SIZE(amn, 1)
            w1 = 0.0; w2 = 0.0
            DO iw = 1, SIZE(amn, 2)
               IF (col_ch(iw) == 1) THEN
                  w1 = w1 + ABS(amn(ib, iw, ik))**2
               ELSE
                  w2 = w2 + ABS(amn(ib, iw, ik))**2
               END IF
            END DO
            lab(ib, ik) = MERGE(1, 2, w1 >= w2)
            !> A band with no projection weight says nothing about its own spin, and there
            !> are always some: only bands that DO project are asked to be clean.
            IF (w1 + w2 > 1.0e-6) worst = MIN(worst, MAX(w1, w2)/(w1 + w2))
         END DO
      END DO
      ok = worst >= PURITY_MIN
      IF (.NOT. ok) WRITE (oUnit, '(a,f6.3,a)') &
         'wannierlib: spin is not a good quantum number here (worst projection purity ', &
         worst, '); spinBalanced does not apply. Drop it to use the ordinary one.'
   END SUBROUTINE bands_by_channel

   !> The channel of each band, from its OWN spin. The quantisation axis is taken from the
   !> data -- the leading eigenvector of sum_{n,k} <sigma><sigma>^T -- and not built from the
   !> noco angles, because in an antiferromagnet those are per atom and no single pair of them
   !> describes the cell, while the global axis is still perfectly well defined.
   SUBROUTINE bands_by_spin(s0, lab, ok)
      COMPLEX, INTENT(IN) :: s0(:, :, :, :)      !> (nb, nb, 3, nk)
      INTEGER, INTENT(OUT) :: lab(:, :)
      LOGICAL, INTENT(OUT) :: ok
      INTEGER :: ib, ik, ic, jc, nb, nk, n_pol, n_col
      REAL :: t(3, 3), ev(3), vec(3, 3), p(3), s, worst, pol, worst_pol
      COMPLEX :: tc(3, 3), vc(3, 3)
      nb = SIZE(s0, 1); nk = SIZE(s0, 4)

      t = 0.0
      DO ik = 1, nk
         DO ib = 1, nb
            DO ic = 1, 3
               p(ic) = REAL(s0(ib, ib, ic, ik))
            END DO
            DO ic = 1, 3
               DO jc = 1, 3
                  t(ic, jc) = t(ic, jc) + p(ic)*p(jc)
               END DO
            END DO
         END DO
      END DO
      tc = CMPLX(t, 0.0)
      CALL herm_eigen(tc, ev, vc)
      vec = REAL(vc)
      !> zheev returns ascending eigenvalues; the axis is the last one.
      p = vec(:, 3)

      !> Two independent questions, and reporting only one of them gave a false diagnosis
      !> for a year of this work. POLARISATION asks whether a state is a spinor at all,
      !> |<sigma>|, which is 1 whichever way it points. COLLINEARITY asks whether it points
      !> along the common axis. A ferromagnet answers 1 and 1. A 3Q antiferromagnet like
      !> Mn3Ir answers 0.90 and 0: its states have perfectly good spin, there is simply no
      !> single axis -- measured, 80.6 % of its bands are pure to better than 0.90, along the
      !> three <110> directions of its sublattices, carrying 36.5 %, 36.5 % and 10.4 %.
      !> Saying that is what points at the generalisation this mode would need (one channel
      !> per sublattice); saying "spin is not a good quantum number" is false and points
      !> nowhere.
      worst = 1.0; worst_pol = 1.0; n_pol = 0; n_col = 0
      DO ik = 1, nk
         DO ib = 1, nb
            s = 0.0; pol = 0.0
            DO ic = 1, 3
               s = s + p(ic)*REAL(s0(ib, ib, ic, ik))
               pol = pol + REAL(s0(ib, ib, ic, ik))**2
            END DO
            pol = SQRT(pol)
            lab(ib, ik) = MERGE(1, 2, s >= 0.0)
            worst = MIN(worst, ABS(s))
            worst_pol = MIN(worst_pol, pol)
            !> Counted as well as minimised: a worst-case over hundreds of bands is set by a
            !> handful of outliers and hides what is actually in the way. In Mn3Ir the worst
            !> polarisation is 0.0000 and twelve non-magnetic Ir bands put it there, while
            !> four bands in five are pure spinors -- along three different axes, which is the
            !> real obstacle and the one a minimum never shows.
            IF (pol >= PURITY_MIN) THEN
               n_pol = n_pol + 1
               IF (ABS(s) >= PURITY_MIN*pol) n_col = n_col + 1
            END IF
         END DO
      END DO
      ok = worst >= PURITY_MIN
      WRITE (oUnit, '(a,3f8.4,a)') &
         'wannierlib: spin axis for the balanced disentanglement (', p, ')'
      WRITE (oUnit, '(a,f6.1,a,f6.1,a)') &
         '  bands with a definite spin ', 100.0*n_pol/REAL(nb*nk), &
         ' %, of those collinear with that axis ', &
         MERGE(100.0*n_col/REAL(MAX(n_pol, 1)), 0.0, n_pol > 0), ' %'
      WRITE (oUnit, '(a,f7.4,a,f7.4)') &
         '  worst band polarisation |<sigma>| ', worst_pol, &
         ', worst projection on that axis ', worst
      IF (.NOT. ok) THEN
         IF (n_pol > 0 .AND. n_col < n_pol) THEN
            WRITE (oUnit, '(a)') &
               '  the bands ARE spin eigenstates, but not along one common axis: this is a '// &
               'genuinely non-collinear texture and it does not fit in two channels. It '// &
               'would need one channel per magnetic sublattice, which this mode does not do.'
         ELSE
            WRITE (oUnit, '(a)') &
               '  some bands carry no spin polarisation at all, so no labelling of them is '// &
               'possible. Drop spinBalanced to use the ordinary disentanglement.'
         END IF
      END IF
   END SUBROUTINE bands_by_spin

   !> Make the band labels and the column labels agree on which channel is number 1. The two
   !> come from different places -- the bands from the sign of their spin, the columns from
   !> proj_spin -- and if they disagree every band is matched to the projections of the other
   !> channel: the subspace comes out orthogonal to what it is supposed to represent, and
   !> nothing reports it. Decided by which pairing carries more projection weight.
   SUBROUTINE match_labels(amn, col_ch, lab)
      COMPLEX, INTENT(IN) :: amn(:, :, :)
      INTEGER, INTENT(IN) :: col_ch(:)
      INTEGER, INTENT(INOUT) :: lab(:, :)
      INTEGER :: ib, iw, ik
      REAL :: same, cross, w
      same = 0.0; cross = 0.0
      DO ik = 1, SIZE(amn, 3)
         DO iw = 1, SIZE(amn, 2)
            DO ib = 1, SIZE(amn, 1)
               w = ABS(amn(ib, iw, ik))**2
               IF (lab(ib, ik) == col_ch(iw)) THEN
                  same = same + w
               ELSE
                  cross = cross + w
               END IF
            END DO
         END DO
      END DO
      IF (cross > same) THEN
         WRITE (oUnit, '(a)') 'wannierlib: swapping the spin labels so that channel 1 of the '// &
            'bands is channel 1 of the projections'
         lab = 3 - lab
      END IF
   END SUBROUTINE match_labels

   !> Souza-Marzari-Vanderbilt, with the eigenvalue selection made once per channel instead
   !> of once globally. That single change is the whole module: without spin-orbit coupling
   !> Z is already block-diagonal, so its eigenvectors need no help -- only the counting does.
   SUBROUTINE smv_balanced(mmn, amn, lab, col_ch, allow, froz, nnkp, wb, nntot, nper, u_opt)
      COMPLEX, INTENT(IN) :: mmn(:, :, :, :), amn(:, :, :)
      COMPLEX, INTENT(INOUT) :: u_opt(:, :, :)
      INTEGER, INTENT(IN) :: lab(:, :), col_ch(:), nnkp(:, :), nntot, nper
      LOGICAL, INTENT(IN) :: allow(:, :), froz(:, :)
      REAL, INTENT(IN) :: wb(:)
      COMPLEX, ALLOCATABLE :: z(:, :, :), zold(:, :, :), v(:, :)
      INTEGER :: nb, nw, nk, ik, nn, jk, it
      REAL :: oi, oi_old
      nb = SIZE(u_opt, 1); nw = SIZE(u_opt, 2); nk = SIZE(u_opt, 3)
      ALLOCATE (z(nb, nb, nk), zold(nb, nb, nk), v(nb, nw))
      zold = CMPLX(0.0, 0.0); oi_old = HUGE(1.0); oi = 0.0

      !> Start from the projections, restricted to the blocks: the same point Wannier90
      !> would start from, minus the columns that cross channels.
      DO ik = 1, nk
         CALL pick_from_matrix(amn(:, :, ik), lab(:, ik), col_ch, allow(:, ik), froz(:, ik), &
                               nper, u_opt(:, :, ik))
      END DO

      DO it = 1, MAXIT
         z = CMPLX(0.0, 0.0)
         DO ik = 1, nk
            DO nn = 1, nntot
               jk = nnkp(ik, nn)
               v = MATMUL(mmn(:, :, nn, ik), u_opt(:, :, jk))
               z(:, :, ik) = z(:, :, ik) + wb(nn)*MATMUL(v, CONJG(TRANSPOSE(v)))
            END DO
         END DO
         !> The mixing goes on Z and not on u_opt: two subspaces from different iterations
         !> have no meaningful average, their Z matrices do. This is what Wannier90 does too.
         IF (it > 1) z = MIXING*z + (1.0 - MIXING)*zold
         zold = z
         DO ik = 1, nk
            CALL pick_from_matrix(z(:, :, ik), lab(:, ik), col_ch, allow(:, ik), froz(:, ik), &
                                  nper, u_opt(:, :, ik), l_eigen=.TRUE.)
         END DO
         oi = omega_i(mmn, u_opt, nnkp, wb, nntot)
         IF (it > 20 .AND. ABS(oi - oi_old) < CONVTOL) EXIT
         oi_old = oi
      END DO
      WRITE (oUnit, '(a,i0,a,f14.9)') '  converged in ', MIN(it, MAXIT), &
         ' iterations, Omega_I = ', oi
   END SUBROUTINE smv_balanced

   !> The columns of u_opt at one k. Frozen bands enter as unit vectors; the rest are either
   !> the leading eigenvectors of Z inside the channel (l_eigen) or the orthonormalised
   !> projections (the starting point). Columns are laid out in the order col_ch gives, so
   !> the Wannier index keeps the meaning the projections gave it.
   SUBROUTINE pick_from_matrix(src, lab, col_ch, allow, froz, nper, u, l_eigen)
      COMPLEX, INTENT(IN) :: src(:, :)
      INTEGER, INTENT(IN) :: lab(:), col_ch(:), nper
      LOGICAL, INTENT(IN) :: allow(:), froz(:)
      COMPLEX, INTENT(OUT) :: u(:, :)
      LOGICAL, INTENT(IN), OPTIONAL :: l_eigen
      INTEGER :: nb, nw, ich, ib, iw, nfr, nfe, n, i, j, icol
      INTEGER, ALLOCATABLE :: idf(:), idr(:), cols(:)
      COMPLEX, ALLOCATABLE :: sub(:, :), vec(:, :)
      REAL, ALLOCATABLE :: ev(:)
      LOGICAL :: eig
      eig = .FALSE.; IF (PRESENT(l_eigen)) eig = l_eigen
      nb = SIZE(src, 1); nw = SIZE(u, 2)
      u = CMPLX(0.0, 0.0)

      DO ich = 1, 2
         nfr = COUNT(froz .AND. lab == ich)
         nfe = nper - nfr
         ALLOCATE (idf(MAX(nfr, 1)), idr(nb), cols(nper))
         n = 0
         DO ib = 1, nb
            IF (froz(ib) .AND. lab(ib) == ich) THEN
               n = n + 1; idf(n) = ib
            END IF
         END DO
         n = 0
         DO ib = 1, nb
            IF (allow(ib) .AND. .NOT. froz(ib) .AND. lab(ib) == ich) THEN
               n = n + 1; idr(n) = ib
            END IF
         END DO
         i = 0
         DO iw = 1, nw
            IF (col_ch(iw) == ich) THEN
               i = i + 1; cols(i) = iw
            END IF
         END DO

         DO i = 1, nfr
            u(idf(i), cols(i)) = CMPLX(1.0, 0.0)
         END DO
         IF (nfe > 0) THEN
            ALLOCATE (sub(n, n), vec(n, n), ev(n))
            IF (eig) THEN
               DO i = 1, n
                  DO j = 1, n
                     sub(i, j) = src(idr(i), idr(j))
                  END DO
               END DO
               CALL herm_eigen(sub, ev, vec)
               !> zheev returns ascending eigenvalues; the subspace wants the largest.
               DO i = 1, nfe
                  DO j = 1, n
                     u(idr(j), cols(nfr + i)) = vec(j, n - i + 1)
                  END DO
               END DO
            ELSE
               CALL lowdin_columns(src, idr(1:n), cols(nfr + 1:nfr + nfe), u)
            END IF
            DEALLOCATE (sub, vec, ev)
         END IF
         DEALLOCATE (idf, idr, cols)
      END DO
   END SUBROUTINE pick_from_matrix

   !> Orthonormalised projections restricted to a set of rows: the Lowdin of the submatrix,
   !> done through its singular vectors so a rank-deficient block fails loudly rather than
   !> silently returning something that is not orthonormal.
   SUBROUTINE lowdin_columns(a, rows, cols, u)
      COMPLEX, INTENT(IN) :: a(:, :)
      INTEGER, INTENT(IN) :: rows(:), cols(:)
      COMPLEX, INTENT(INOUT) :: u(:, :)
      COMPLEX, ALLOCATABLE :: x(:, :), uu(:, :), vt(:, :), work(:)
      REAL, ALLOCATABLE :: s(:), rwork(:)
      INTEGER :: m, n, i, j, info, lwork
      m = SIZE(rows); n = SIZE(cols)
      ALLOCATE (x(m, n), uu(m, m), vt(n, n), s(MIN(m, n)), rwork(5*MIN(m, n)))
      DO j = 1, n
         DO i = 1, m
            x(i, j) = a(rows(i), cols(j))
         END DO
      END DO
      lwork = 4*MAX(m, n) + 64*MAX(m, n)
      ALLOCATE (work(lwork))
      CALL zgesvd('A', 'A', m, n, x, m, s, uu, m, vt, n, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('wannierlib: SVD failed in the balanced disentanglement', &
                                      calledby='lowdin_columns')
      DO j = 1, n
         DO i = 1, m
            u(rows(i), cols(j)) = SUM(uu(i, 1:n)*vt(1:n, j))
         END DO
      END DO
   END SUBROUTINE lowdin_columns

   SUBROUTINE herm_eigen(a, ev, vec)
      COMPLEX, INTENT(IN) :: a(:, :)
      REAL, INTENT(OUT) :: ev(:)
      COMPLEX, INTENT(OUT) :: vec(:, :)
      COMPLEX, ALLOCATABLE :: work(:)
      REAL, ALLOCATABLE :: rwork(:)
      INTEGER :: n, lwork, info
      n = SIZE(a, 1)
      vec = a
      lwork = 2*n + 64*n
      ALLOCATE (work(lwork), rwork(MAX(1, 3*n - 2)))
      CALL zheev('V', 'U', n, vec, n, ev, work, lwork, rwork, info)
      IF (info /= 0) CALL juDFT_error('wannierlib: zheev failed in the balanced disentanglement', &
                                      calledby='herm_eigen')
   END SUBROUTINE herm_eigen

   !> The invariant spread of a subspace, which is what this iteration maximises the negative
   !> of. Reported so the balanced result can be compared with the ordinary one directly.
   REAL FUNCTION omega_i(mmn, u, nnkp, wb, nntot)
      COMPLEX, INTENT(IN) :: mmn(:, :, :, :), u(:, :, :)
      INTEGER, INTENT(IN) :: nnkp(:, :), nntot
      REAL, INTENT(IN) :: wb(:)
      COMPLEX, ALLOCATABLE :: mw(:, :)
      INTEGER :: nk, nw, ik, nn, jk
      REAL :: acc
      nk = SIZE(u, 3); nw = SIZE(u, 2)
      ALLOCATE (mw(nw, nw))
      acc = 0.0
      DO ik = 1, nk
         DO nn = 1, nntot
            jk = nnkp(ik, nn)
            mw = MATMUL(CONJG(TRANSPOSE(u(:, :, ik))), &
                        MATMUL(mmn(:, :, nn, ik), u(:, :, jk)))
            acc = acc + wb(nn)*(REAL(nw) - SUM(ABS(mw)**2))
         END DO
      END DO
      omega_i = acc/REAL(nk)
   END FUNCTION omega_i

END MODULE m_wannierlib_disentangle_spin
