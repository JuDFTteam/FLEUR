!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The table the muffin-tin <u_{k+b1}|H_k|u_{k+b2}> contracts against, in the two halves
!>  m_melem_ujugaunt is built in: a radial pass per pair of b moduli, and an angular pass per
!>  pair of b vectors.
!>
!>  It is m_melem_ujugaunt with a Hamiltonian in the middle, and the middle adds a dimension.
!>  H_sph is diagonal in angular momentum, so one side reaches an intermediate channel through
!>  its own Bessel function and that channel reaches the other side through the other Bessel
!>  function:
!>
!>      lB,mB  --(j of bB)-->  lt,mt  --(j of bA)-->  lA,mA
!>
!>  Two Gaunt coefficients where the overlap needs one, and the radial integral is taken at
!>  the intermediate channel lt rather than at either end. That index is easy to leave out
!>  and nothing downstream would notice it missing.
!>
!>  The radial pass runs per pair of moduli, the angular pass per pair of vectors.
MODULE m_melem_uhu_gaunt
   USE m_juDFT
   USE m_types_atoms
   USE m_ylm
   USE m_gaunt, ONLY: gaunt1
   USE m_constants, ONLY: ImagUnit
   USE m_melem_uhu_radint, ONLY: melem_uhu_radint
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: melem_uhu_radpass, melem_uhu_angpass
CONTAINS

   !> The radial half: < r u_lA j_jA | H_sph(lt) | r u_lB j_jB > for every pair of radial
   !> slots and every angular combination the two Gaunt coefficients can reach. Both sides
   !> arrive already dressed by m_melem_uhu_dress, one per modulus, so nothing is
   !> differentiated here.
   !>
   !> WHICH ENERGY THE MASS ENHANCEMENT SEES: the average of the two linearisation energies,
   !> epar(lA) and epar(lB). The reference leaves the intermediate channel's own epar(lt)
   !> commented beside it, so the choice was made and not inherited; it matters only through
   !> (e - V)/2c^2 and is recorded here because nothing in the result shows which was used.
   !>
   !> The selection rules are the ones the Gaunt coefficients will impose later, applied here
   !> so the integral is not run for combinations that get multiplied by zero. Skipping them
   !> would be correct and roughly an order of magnitude slower.
   !>
   !> THE SYMMETRY IT HAS TO HAVE: m_melem_uhu_radint applies H to each side and averages,
   !> so exchanging the two sides together with their slots and channels returns the same
   !> number. asym_max carries the largest R_MT surface term met instead, which is the part
   !> that exchanging the sides does not cancel.
   SUBROUTINE melem_uhu_radpass(atoms, n, lwn, epar_a, epar_b, vr0, &
                                uj_a, duj_a, uj_b, duj_b, uhu, asym_max)
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: n                      !> atom type
      INTEGER, INTENT(IN) :: lwn                    !> lmax of this type
      REAL, INTENT(IN)    :: epar_a(0:), epar_b(0:) !> linearisation energies, per l
      REAL, INTENT(IN)    :: vr0(:)                 !> r*V, spherical channel, Y_00 taken out
      REAL, INTENT(IN)    :: uj_a(:, :, :, 0:, 0:)  !> (jri,2,slot,l,L) side A, at |bA|
      REAL, INTENT(IN)    :: duj_a(:, :, :, 0:, 0:)
      REAL, INTENT(IN)    :: uj_b(:, :, :, 0:, 0:)  !> side B, at |bB|
      REAL, INTENT(IN)    :: duj_b(:, :, :, 0:, 0:)
      REAL, INTENT(INOUT) :: uhu(:, :, 0:, 0:, 0:, 0:, 0:)  !> (r1,r2,lA,lB,jA,jB,lt)
      REAL, INTENT(INOUT) :: asym_max

      INTEGER :: lt, la, lb, ja, jb, r1, r2, jri, nra, nrb
      REAL    :: e, val, asym

      CALL timestart("melem_uhu_radpass")
      jri = atoms%jri(n)
      nra = SIZE(uj_a, 3)
      nrb = SIZE(uj_b, 3)

      DO lt = 0, lwn
         DO la = 0, lwn
            DO ja = 0, lwn
               !> gaunt1(lt, ja, la, ...) vanishes outside the triangle and for odd parity
               IF (la < ABS(lt - ja) .OR. la > lt + ja) CYCLE
               IF (MOD(lt + ja + la, 2) /= 0) CYCLE
               DO lb = 0, lwn
                  e = 0.5*(epar_a(la) + epar_b(lb))
                  DO jb = 0, lwn
                     !> and gaunt1(lb, jb, lt, ...) the same way
                     IF (lt < ABS(lb - jb) .OR. lt > lb + jb) CYCLE
                     IF (MOD(lb + jb + lt, 2) /= 0) CYCLE
                     DO r1 = 1, nra
                        DO r2 = 1, nrb
                           CALL melem_uhu_radint(jri, atoms%rmsh(1:, n), atoms%dx(n), e, vr0, &
                                                 uj_a(:, :, r1, la, ja), uj_b(:, :, r2, lb, jb), &
                                                 duj_a(:, :, r1, la, ja), duj_b(:, :, r2, lb, jb), &
                                                 lt, val, asym)
                           uhu(r1, r2, la, lb, ja, jb, lt) = val
                           asym_max = MAX(asym_max, ABS(asym))
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL timestop("melem_uhu_radpass")
   END SUBROUTINE melem_uhu_radpass

   !> The angular half: the two Gaunt coefficients, the two spherical harmonics of the two b
   !> directions, and the phase, accumulated onto the (lm, lm) table the contraction reads.
   !>
   !> The loop bounds and the selection rules are transcribed rather than rederived. Each pair
   !> of nested limits is a Gaunt condition written as a range, and rewriting them from the
   !> conditions is where an off-by-one hides that no test downstream would catch:
   !>
   !>      G1 = G( (lt,mt), (ja,mja), (la,mla) )      side A, at bA, reaching the intermediate
   !>      G2 = G( (lb,mb), (jb,mjb), (lt,mt) )       side B, at bB, leaving it
   !>
   !> with m1 = m2 + m3 on each, which is what fixes mt and mb rather than looping over them.
   !>
   !> NO EXCHANGE AVERAGE HERE, and that is a choice the reference makes too. It carries
   !> weights sy1, sy2 that combine this ordering with its mirror, and sets them to 1 and 0
   !> whenever the radial integral is the symmetrised one; the averaging is the alternative to
   !> symmetrising inside the integral, not an addition to it. m_melem_uhu_radint symmetrises,
   !> so the mirror term would double-count.
   !>
   !> INDEX ORDER, WHICH NO IDENTITY HERE CAN SETTLE: the table goes out with the KET side
   !> first and the bra second, while the radial slots stay in bra-then-ket order. That
   !> asymmetry is not ours, it is what m_melem_ujugaunt builds and what m_melem_mmkb_sph
   !> contracts, where conjg of the ket's coefficients closes on the first lm index.
   !> Hermiticity at bA = bB survives transposing the pair, so it holds either way round
   !> and only the shape of the contraction decides which way is right.
   SUBROUTINE melem_uhu_angpass(atoms, n, lwn, bkrot_a, bkrot_b, uhu, ipair, uhug)
      TYPE(t_atoms), INTENT(IN) :: atoms
      INTEGER, INTENT(IN) :: n, lwn
      REAL, INTENT(IN)    :: bkrot_a(3), bkrot_b(3)   !> the two b, in cartesian
      REAL, INTENT(IN)    :: uhu(:, :, 0:, 0:, 0:, 0:, 0:)  !> (r1,r2,la,lb,ja,jb,lt)
      INTEGER, INTENT(IN) :: ipair                     !> which pair of vectors this is
      COMPLEX, INTENT(INOUT) :: uhug(0:, 0:, :, :, :, :)   !> (lm_ket,lm_bra,r_bra,r_ket,type,pair)

      COMPLEX :: yla((atoms%lmaxd + 1)**2), ylb((atoms%lmaxd + 1)**2)
      COMPLEX :: cci, fac1, cil1, cil
      REAL    :: gc1, gc2
      INTEGER :: ja, mja, la, mla, lt, mt, jb, mjb, lb, mb
      INTEGER :: lmja, lmjb, lma, lmb, lo, hi, hic, lo2, hi2, hic2

      CALL timestart("melem_uhu_angpass")
      cci = CMPLX(0.0, -1.0)
      CALL ylm4(lwn, bkrot_a, yla)
      CALL ylm4(lwn, bkrot_b, ylb)

      DO ja = 0, lwn
         DO mja = -ja, ja
            lmja = ja*(ja + 1) + mja
            fac1 = (cci**ja)*CONJG(yla(lmja + 1))
            IF (fac1 == CMPLX(0.0, 0.0)) CYCLE

            DO la = 0, lwn
               lo = ABS(ja - la)
               hi = ja + la
               hic = MIN(hi, lwn)
               DO mla = -la, la
                  lma = la*(la + 1) + mla
                  mt = mja + mla
                  IF (ABS(mt) > hic) CYCLE

                  DO lt = lo, hic
                     IF (MOD(hi + lt, 2) == 1 .OR. ABS(mt) > lt) CYCLE
                     gc1 = gaunt1(lt, ja, la, mt, mja, mla, atoms%lmaxd)
                     cil1 = (ImagUnit**la)*fac1*gc1

                     DO jb = 0, lwn
                        lo2 = ABS(lt - jb)
                        hi2 = lt + jb
                        hic2 = MIN(hi2, lwn)
                        DO mjb = -jb, jb
                           lmjb = jb*(jb + 1) + mjb
                           mb = mjb + mt
                           IF (ABS(mb) > hic2) CYCLE

                           DO lb = lo2, hic2
                              IF (MOD(hi2 + lb, 2) == 1 .OR. ABS(mb) > lb) CYCLE
                              lmb = lb*(lb + 1) + mb
                              gc2 = gaunt1(lb, jb, lt, mb, mjb, mt, atoms%lmaxd)
                              cil = cil1*(ImagUnit**(jb - lb))*gc2*CONJG(ylb(lmjb + 1))
                              IF (cil == CMPLX(0.0, 0.0)) CYCLE

                              uhug(lmb, lma, :, :, n, ipair) = uhug(lmb, lma, :, :, n, ipair) &
                                 + cil*uhu(:, :, la, lb, ja, jb, lt)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      CALL timestop("melem_uhu_angpass")
   END SUBROUTINE melem_uhu_angpass

END MODULE m_melem_uhu_gaunt
