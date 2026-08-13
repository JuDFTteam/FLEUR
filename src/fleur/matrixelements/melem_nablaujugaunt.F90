!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The muffin-tin half of <u_a| e^{-i b.r} nabla |u_b>, tabulated the way the plain
!>  overlap is: radial integrals against a Bessel function, times an angular factor,
!>  once per distinct b so that the state loop only has matrices to multiply.
!>
!>  The gradient of a radial function times a spherical harmonic reaches two angular
!>  momenta, l-1 and l+1, and the two carry different radial combinations. Both are
!>  therefore built here and kept apart until the angular factors, which differ, are
!>  applied. In the notation of momentum_matrix, from which the two channels and their
!>  coefficients are taken, l+1 is the f channel and l-1 the g channel.
!>
!>  Radial functions arrive as r*u(r), which is what makes the integrand carry no extra
!>  power of r and what puts the -p/r term into the derivative.
!>
!>  The radial derivative is taken from the equations the functions solve, not from the mesh:
!>  the scalar-relativistic solver advances r*dp/dr = rm*q + p, so the derivative is already
!>  known wherever the function is, and its accuracy is that of the function rather than of a
!>  difference formula. The small component is stored divided by c, which is why it comes back
!>  multiplied by it.
!>
!>  The energy derivative solves the same pair with a source, and the source is not the
!>  function itself: it is what differentiating the pair with respect to the energy brings
!>  down, and it comes out crossed, the equation for the large component fed by the small
!>  one and the other way round. radsrd writes it as t1 = cin*r*q and t2 = r*p.
MODULE m_melem_nablaujugaunt
   USE m_juDFT
   USE m_types_atoms
   USE m_types_cell
   USE m_types_radfun
   USE m_sphbes
   USE m_ylm
   USE m_intgr, ONLY: intgr3
   USE m_gaunt, ONLY: gaunt1
   USE m_types_enpara
   USE m_constants, ONLY: c_light
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_nablaujugaunt
CONTAINS

   !> The table, one entry per vector of kdiff and per spherical component of the
   !> gradient. The component index runs -1, 0, +1 and is turned into a cartesian one by
   !> the caller, which is where the two conventions meet:
   !>
   !>     x = (m_{-1} - m_{+1})/sqrt(2),  y = -i(m_{-1} + m_{+1})/sqrt(2),  z = m_0
   !>
   !> and the whole thing carries the -i that makes nabla into the momentum.
   !>
   !> The gradient acts on the SECOND side, the one whose coefficients multiply from the
   !> right, so it is that side's angular momentum that is displaced by one.
   SUBROUTINE melem_nablaujugaunt(atoms, cell, enpara, nntot, kdiff, radfun, radfun_b, &
                                  jspin, jspin_b, nujug)
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)   :: cell
      !> Carries the potential the radial functions were generated with, and their energy
      !> parameters. Both are needed by the equations and neither may come from anywhere else.
      TYPE(t_enpara), INTENT(IN) :: enpara
      INTEGER, INTENT(IN)        :: nntot
      REAL, INTENT(IN)           :: kdiff(:, :)
      TYPE(t_radfun), INTENT(IN) :: radfun(atoms%ntype)
      TYPE(t_radfun), INTENT(IN) :: radfun_b(atoms%ntype)
      INTEGER, INTENT(IN)        :: jspin, jspin_b
      !> (lm, lm', r1, r2, type, b, component) with the component running -1:1
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: nujug(:, :, :, :, :, :, :)

      !> The radial halves, one per channel: nuju_p reaches lp+1, nuju_m reaches lp-1.
      REAL, ALLOCATABLE :: nuju_p(:, :, :, :, :, :), nuju_m(:, :, :, :, :, :)
      REAL, ALLOCATABLE :: rk_uni(:), gplus(:, :), gminus(:, :), x(:), jj(:, :)
      INTEGER, ALLOCATABLE :: imod(:)
      REAL    :: bpt(3), bkrot(3), rk, res
      COMPLEX :: ylmpp((atoms%lmaxd + 1)**2), factor, fac_ch
      REAL    :: fco(-1:1), gco(-1:1)
      REAL    :: cc, cin2, esl, rr, rm_sra, pv, qv, src_p, src_q, rve
      INTEGER :: n, lwn, total_nr, nmod, q, ikpt_b, i, ch, lpe, mpe, ilo, nlo_seen
      INTEGER :: l, m, lm, lp, mp, lpmp, lpp, mpp, llpp, lmpp, lmini, lmaxi, r1, r2
      REAL    :: jlpp(0:atoms%lmaxd + 1)

      CALL timestart("melem_nablaujugaunt")

      cc = c_light(1.0)
      cin2 = 1.0/(cc*cc)

      total_nr = 0
      DO n = 1, atoms%ntype
         total_nr = MAX(total_nr, MAXVAL(radfun(n)%n_r(:)))
      END DO

      !> Same grouping as the plain overlap: the modulus of b enters only through
      !> j_l(|b| r), so vectors that share it share the whole radial table. Grouped on the
      !> exact bit pattern, so the result does not depend on the grouping.
      ALLOCATE (rk_uni(nntot), imod(nntot))
      nmod = 0
      DO ikpt_b = 1, nntot
         bpt(:) = kdiff(:, ikpt_b)
         rk = SQRT(DOT_PRODUCT(bpt, MATMUL(cell%bbmat, bpt)))
         imod(ikpt_b) = 0
         DO q = 1, nmod
            IF (rk == rk_uni(q)) THEN
               imod(ikpt_b) = q
               EXIT
            END IF
         END DO
         IF (imod(ikpt_b) == 0) THEN
            nmod = nmod + 1
            rk_uni(nmod) = rk
            imod(ikpt_b) = nmod
         END IF
      END DO

      ALLOCATE (nuju_p(total_nr, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:atoms%lmaxd, nmod), source=0.0)
      ALLOCATE (nuju_m(total_nr, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:atoms%lmaxd, nmod), source=0.0)
      ALLOCATE (nujug(0:atoms%lmaxd*(atoms%lmaxd + 2), 0:atoms%lmaxd*(atoms%lmaxd + 2), &
                      total_nr, total_nr, atoms%ntype, nntot, -1:1), source=CMPLX(0.0, 0.0))
      ALLOCATE (jj(0:atoms%lmaxd, atoms%jmtd), x(atoms%jmtd))
      ALLOCATE (gplus(atoms%jmtd, 2), gminus(atoms%jmtd, 2))

      DO n = 1, atoms%ntype
         lwn = atoms%lmax(n)
         nuju_p = 0.0
         nuju_m = 0.0

         DO q = 1, nmod
            DO i = 1, atoms%jri(n)
               CALL sphbes(lwn, rk_uni(q)*atoms%rmsh(i, n), jlpp)
               jj(:, i) = jlpp(0:lwn)
            END DO

            DO lp = 0, lwn
               DO r2 = 1, radfun_b(n)%n_r(lp)
                  !> The two channels of the gradient. With r*u(r) stored, d/dr of the
                  !> function itself brings the -p/r, and the two angular channels then
                  !> differ by which multiple of p/r is added.
                  !> Which energy the slot was solved at: the first two share the one of the
                  !> l channel, and every one after them is a local orbital with its own.
                  esl = enpara%el0(lp, n, jspin_b)
                  IF (r2 > 2) THEN
                     nlo_seen = 2
                     DO ilo = 1, atoms%nlo(n)
                        IF (atoms%llo(ilo, n) /= lp) CYCLE
                        nlo_seen = nlo_seen + 1
                        IF (nlo_seen == r2) THEN
                           esl = enpara%ello0(ilo, n, jspin_b)
                           EXIT
                        END IF
                     END DO
                  END IF
                  DO i = 1, atoms%jri(n)
                     rr = atoms%rmsh(i, n)
                     rve = enpara%vr(i, n, jspin_b) - rr*esl
                     rm_sra = 2.0*rr - cin2*rve
                     pv = radfun_b(n)%r(i, 1, r2, lp, jspin_b)
                     qv = cc*radfun_b(n)%r(i, 2, r2, lp, jspin_b)
                     !> Orthogonalising the energy derivative against u adds a homogeneous
                     !> piece, which leaves its source unchanged.
                     src_p = 0.0; src_q = 0.0
                     IF (r2 == 2) THEN
                        !> The two sources are CROSSED: the equation for the large
                        !> component is fed by the small one and the other way round.
                        !> That is what differentiating the equations with respect to
                        !> the energy brings down, since rm carries e through the large
                        !> equation and rve through the small one. radsrd writes them as
                        !> t1 = cin*r*q and t2 = r*p, with q the STORED small component,
                        !> already divided by c; dividing by r leaves what follows.
                        src_p = radfun_b(n)%r(i, 2, 1, lp, jspin_b)/cc
                        src_q = -radfun_b(n)%r(i, 1, 1, lp, jspin_b)/cc
                     END IF
                     gplus(i, 1) = (rm_sra*qv - REAL(lp)*pv)/rr + src_p
                     gminus(i, 1) = (rm_sra*qv + REAL(lp + 1)*pv)/rr + src_p
                     !> and the small component, brought back to the scaling it is stored in.
                     gplus(i, 2) = (((REAL(lp*(lp + 1))/rm_sra + rve)*pv - qv) &
                                    - REAL(lp + 1)*qv)/(rr*cc) + src_q
                     gminus(i, 2) = (((REAL(lp*(lp + 1))/rm_sra + rve)*pv - qv) &
                                     + REAL(lp)*qv)/(rr*cc) + src_q
                  END DO

                  DO l = 0, lwn
                     DO r1 = 1, radfun(n)%n_r(l)
                        DO lpp = 0, lwn
                           !> The l+1 channel, kept only where that momentum exists. The
                           !> Bessel order has to lie between the two angular momenta it
                           !> couples, and the one on the ket side is the DISPLACED one.
                           IF (lp + 1 <= lwn) THEN
                              lmini = ABS(l - (lp + 1))
                              lmaxi = l + (lp + 1)
                              IF (.NOT. ((MOD(l + lp + 1 + lpp, 2) == 1) .OR. (lpp < lmini) .OR. (lpp > lmaxi))) THEN
                                 DO i = 1, atoms%jri(n)
                                    x(i) = (radfun(n)%r(i, 1, r1, l, jspin)*gplus(i, 1) + &
                                            radfun(n)%r(i, 2, r1, l, jspin)*gplus(i, 2))*jj(lpp, i)
                                 END DO
                                 CALL intgr3(x, atoms%rmsh(1:, n), atoms%dx(n), atoms%jri(n), res)
                                 nuju_p(r1, r2, l, lp, lpp, q) = res
                              END IF
                           END IF
                           !> and the l-1 channel, with its own displaced momentum.
                           IF (lp - 1 >= 0) THEN
                              lmini = ABS(l - (lp - 1))
                              lmaxi = l + (lp - 1)
                              IF (.NOT. ((MOD(l + lp - 1 + lpp, 2) == 1) .OR. (lpp < lmini) .OR. (lpp > lmaxi))) THEN
                                 DO i = 1, atoms%jri(n)
                                    x(i) = (radfun(n)%r(i, 1, r1, l, jspin)*gminus(i, 1) + &
                                            radfun(n)%r(i, 2, r1, l, jspin)*gminus(i, 2))*jj(lpp, i)
                                 END DO
                                 CALL intgr3(x, atoms%rmsh(1:, n), atoms%dx(n), atoms%jri(n), res)
                                 nuju_m(r1, r2, l, lp, lpp, q) = res
                              END IF
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO

         !> The angular half: once per vector of kdiff, and now also per channel, because
         !> the displaced angular momentum changes both the Gaunt coefficient and the
         !> harmonic the coefficients belong to.
         DO ikpt_b = 1, nntot
            bpt(:) = kdiff(:, ikpt_b)
            bkrot = MATMUL(bpt, cell%bmat)
            CALL ylm4(lwn, bkrot, ylmpp)

            DO l = 0, lwn
               DO m = -l, l
                  lm = l*(l + 1) + m
                  DO lp = 0, lwn
                     DO mp = -lp, lp
                        lpmp = lp*(lp + 1) + mp
                        CALL nabla_coeff(lp, mp, fco, gco)
                        DO ch = 1, 2
                           lpe = MERGE(lp + 1, lp - 1, ch == 1)
                           IF (lpe < 0 .OR. lpe > lwn) CYCLE
                           DO i = -1, 1
                              !> The channel index is not the spherical component of the
                              !> gradient: it labels the three vectors the cartesian matrix
                              !> below turns into x, y and z. Its coefficient and its
                              !> displacement of m both belong to that labelling.
                              mpe = mp - i
                              IF (ABS(mpe) > lpe) CYCLE
                              fac_ch = MERGE(fco(i), gco(i), ch == 1)
                              IF (fac_ch == 0.0) CYCLE
                              DO lpp = 0, lwn
                                 llpp = lpp*(lpp + 1)
                                 mpp = mpe - m
                                 lmpp = llpp + mpp
                                 lmini = ABS(l - lpp)
                                 lmaxi = l + lpp
                                 IF ((lmini <= lpe) .AND. (lpe <= lmaxi) .AND. &
                                     (MOD(l + lpe + lpp, 2) == 0) .AND. (ABS(mpp) <= lpp)) THEN
                                    !> The power of i belongs to the basis functions, not to
                                    !> the gradient: it is their own i^l convention, and the
                                    !> coefficients this table is contracted against are the
                                    !> ones of lp, not of the momentum the gradient displaces
                                    !> it to. So the undisplaced lp goes here, exactly as in
                                    !> the plain overlap. The Gaunt and the harmonic do take
                                    !> the displaced pair, because those are the integral.
                                    factor = CONJG(ylmpp(lmpp + 1))*(CMPLX(0.0, 1.0)**(l + lpp - lp)) &
                                             *gaunt1(lpe, lpp, l, mpe, mpp, m, atoms%lmaxd)*fac_ch
                                    IF (ch == 1) THEN
                                       nujug(lpmp, lm, :, :, n, ikpt_b, i) = nujug(lpmp, lm, :, :, n, ikpt_b, i) &
                                                                             + factor*nuju_p(:, :, l, lp, lpp, imod(ikpt_b))
                                    ELSE
                                       nujug(lpmp, lm, :, :, n, ikpt_b, i) = nujug(lpmp, lm, :, :, n, ikpt_b, i) &
                                                                             + factor*nuju_m(:, :, l, lp, lpp, imod(ikpt_b))
                                    END IF
                                 END IF
                              END DO
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (nuju_p, nuju_m, jj, x, gplus, gminus, rk_uni, imod)
      CALL timestop("melem_nablaujugaunt")
   END SUBROUTINE melem_nablaujugaunt

   !> The angular coefficients of the two channels, taken from momentum_matrix: f for
   !> l -> l+1 and g for l -> l-1, one per spherical component of the gradient. They are
   !> the Clebsch-Gordan factors of coupling a vector to a spherical harmonic, and they
   !> are what the plain overlap has no counterpart for.
   PURE SUBROUTINE nabla_coeff(l, m, fco, gco)
      INTEGER, INTENT(IN) :: l, m
      REAL, INTENT(OUT) :: fco(-1:1), gco(-1:1)

      REAL :: rl, rm

      REAL :: fb(-1:1), gb(-1:1)
      INTEGER :: mu

      rl = REAL(l); rm = REAL(m)
      fb = 0.0; gb = 0.0

      fb(-1) = -SQRT((rl + rm + 1)*(rl + rm + 2)/(2*(2*rl + 1)*(2*rl + 3)))
      fb(0) = SQRT((rl - rm + 1)*(rl + rm + 1)/((2*rl + 1)*(2*rl + 3)))
      fb(1) = -SQRT((rl - rm + 1)*(rl - rm + 2)/(2*(2*rl + 1)*(2*rl + 3)))

      IF (l >= 1) THEN
         IF ((rl - rm)*(rl - rm - 1) > 0.0) gb(-1) = SQRT((rl - rm)*(rl - rm - 1)/(2*(2*rl - 1)*(2*rl + 1)))
         IF ((rl - rm)*(rl + rm) > 0.0) gb(0) = SQRT((rl - rm)*(rl + rm)/((2*rl - 1)*(2*rl + 1)))
         IF ((rl + rm)*(rl + rm - 1) > 0.0) gb(1) = SQRT((rl + rm)*(rl + rm - 1)/(2*(2*rl - 1)*(2*rl + 1)))
      END IF

      !> The muffin-tin half is not hermitian on its own: integrating by parts leaves a
      !> surface term which only the interstitial half cancels, so neither half can be
      !> checked against hermiticity and only their sum can.
      DO mu = -1, 1
         fco(mu) = fb(mu)
         gco(mu) = gb(mu)
      END DO
   END SUBROUTINE nabla_coeff

END MODULE m_melem_nablaujugaunt
