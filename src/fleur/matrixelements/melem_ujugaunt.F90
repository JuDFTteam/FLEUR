!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_melem_ujugaunt
   USE m_juDFT
   USE m_types_atoms
   USE m_types_cell
   USE m_types_radfun
   USE m_sphbes
   USE m_ylm
   USE m_intgr, ONLY: intgr3
   USE m_gaunt, ONLY: gaunt1
   IMPLICIT NONE
  !> One routine, and the radial tables it builds are its only product.
   PRIVATE
   PUBLIC :: melem_ujugaunt
CONTAINS
!   melem_ujugaunt calculates integrals of radial wave functions with
!   bessel functions and multiplies them with an angular factor.
!   Calculating them only once gives some speed-up of melem_mmkb_sph.
!   Frank Freimuth, October 2006, refactored DW in 2026
   SUBROUTINE melem_ujugaunt(atoms, cell, nntot, kdiff, radfun, radfun_b, jspin, jspin_b, &
                                  l_q, sign_q, ujug)
      TYPE(t_atoms), INTENT(IN)  :: atoms
      TYPE(t_cell), INTENT(IN)   :: cell
      INTEGER, INTENT(IN)        :: nntot
      REAL, INTENT(IN)           :: kdiff(:, :)
      TYPE(t_radfun), INTENT(IN) :: radfun(atoms%ntype)
      TYPE(t_radfun), INTENT(IN) :: radfun_b(atoms%ntype)
      INTEGER, INTENT(IN)        :: jspin
      INTEGER, INTENT(IN)        :: jspin_b
      LOGICAL, INTENT(IN)        :: l_q
      INTEGER, INTENT(IN)        :: sign_q

      COMPLEX, ALLOCATABLE, intent(out) :: ujug(:, :, :, :, :, :)

      REAL, ALLOCATABLE :: uju(:, :, :, :, :, :)
      REAL, ALLOCATABLE :: rk_uni(:)
      INTEGER, ALLOCATABLE :: imod(:)
      INTEGER :: ikpt_b, i, lwn, n, lpp, l, lp, q, nmod
      INTEGER :: lmini, lmaxi, m, mp, llpp, mpp, r1, r2
      INTEGER :: lmpp, lm, lpmp, total_nr
      REAL :: rk, bpt(3), gs, jlpp(0:atoms%lmaxd)
      REAL :: jj(0:atoms%lmaxd, atoms%jmtd), x(atoms%jmtd)
      REAL :: bkrot(3)
      COMPLEX :: ylmpp((atoms%lmaxd + 1)**2), factor

      CALL timestart("melem_ujugaunt")
      total_nr = 0
      DO n = 1, atoms%ntype
         total_nr = max(total_nr, maxval(radfun(n)%n_r(:)))
      end do
      !> The modulus of b enters only through j_l(|b| r), inside the radial integral; the
      !> direction enters through Y_lm(b) and the Gaunt coefficients, which are products.
      !> Entries of kdiff that share a modulus therefore share the whole radial table.
      !> Grouped on the exact bit pattern, not on a tolerance: two moduli that differ in
      !> the last bits stay apart, so the result does not depend on the grouping.
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
      ALLOCATE (uju(total_nr, total_nr, 0:atoms%lmaxd, 0:atoms%lmaxd, 0:atoms%lmaxd, nmod), source=0.0)
      allocate (ujug( 0:atoms%lmaxd*(atoms%lmaxd+2), 0:atoms%lmaxd*(atoms%lmaxd+2), total_nr, total_nr,atoms%ntype, nntot), source=CMPLX(0.0, 0.0))

      DO n = 1, atoms%ntype
         lwn = atoms%lmax(n)

         !> The radial half: once per distinct modulus, not once per vector. Zeroed whole,
         !> because the branch that survives the selection rule fills only the block that
         !> this type's n_r reaches and the rest is read back as part of the slice.
         uju = 0.0
         DO q = 1, nmod
            DO i = 1, atoms%jri(n)
               gs = rk_uni(q)*atoms%rmsh(i, n)
               CALL sphbes(lwn, gs, jlpp)
               jj(:, i) = jlpp(:)
            END DO
            do lpp = 0, lwn   ! lpp is the ang. momentum of the bessel function
            DO lp = 0, lwn
               DO l = 0, lwn
                  lmini = ABS(lp - l)
                  lmaxi = lp + l
                  IF (.NOT. ((MOD(l + lp + lpp, 2) == 1) .OR. (lpp < lmini) .OR. (lpp > lmaxi))) THEN
                     do r1 = 1, radfun(n)%n_r(l)
                        do r2 = 1, radfun_b(n)%n_r(lp)
                           DO i = 1, atoms%jri(n)
                              x(i) = (radfun(n)%r(i, 1, r1, l, jspin)*radfun_b(n)%r(i, 1, r2, lp, jspin_b) + &
                                      radfun(n)%r(i, 2, r1, l, jspin)*radfun_b(n)%r(i, 2, r2, lp, jspin_b))* &
                                     jj(lpp, i)
                           END DO
                           CALL intgr3(x, atoms%rmsh(1:, n), atoms%dx(n), atoms%jri(n), uju(r1, r2, l, lp, lpp, q))
                        end do
                     end do
                  END IF
               END DO
            END DO
            end do
         END DO

         !> The angular half: once per vector of kdiff.
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
                        DO lpp = 0, lwn
                           llpp = lpp*(lpp + 1)
                           mpp = mp - m
                           lmpp = llpp + mpp
                           lmini = ABS(l - lpp)
                           lmaxi = l + lpp
                           IF ((lmini <= lp) .AND. (lp <= lmaxi) .AND. (MOD(l + lp + lpp, 2) == 0) .AND. (ABS(mpp) <= lpp)) THEN
                       factor = CONJG(ylmpp(lmpp + 1))*(cmplx(0.0, 1.0)**(l + lpp - lp))*gaunt1(lp, lpp, l, mp, mpp, m, atoms%lmaxd)
                              IF (l_q) factor = (sign_q**lpp)*factor
                              ujug(lpmp, lm, :,:, n, ikpt_b) = ujug(lpmp, lm, :,:, n, ikpt_b) &
                                                              + factor*uju(:, :, l, lp, lpp, imod(ikpt_b))
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         end do
      end do

      CALL timestop("melem_ujugaunt")
   END SUBROUTINE melem_ujugaunt

END MODULE m_melem_ujugaunt
