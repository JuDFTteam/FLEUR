!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_subvxc

CONTAINS

   SUBROUTINE subvxc(lapw, bk, DIMENSION, input, jsp, vr0, atoms, usdus, mpbasis, hybrid, el, ello, sym, &
                     cell, sphhar, stars, xcpot, mpi, oneD, hmat, vx)

      USE m_types
      USE m_judft
      USE m_intgr, ONLY: intgr3
      USE m_constants
      USE m_gaunt, ONLY: gaunt1
      USE m_wrapper
      USE m_loddop
      USE m_radflo
      USE m_radfun
      USE m_abcof3

      IMPLICIT NONE

      CLASS(t_xcpot), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: mpi
      TYPE(t_dimension), INTENT(IN)    :: dimension
      TYPE(t_oneD), INTENT(IN)    :: oneD
      TYPE(t_mpbasis), intent(inout) :: mpbasis
      TYPE(t_hybrid), INTENT(INOUT) :: hybrid
      TYPE(t_input), INTENT(IN)    :: input
      TYPE(t_sym), INTENT(IN)    :: sym
      TYPE(t_stars), INTENT(IN)    :: stars
      TYPE(t_cell), INTENT(IN)    :: cell
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_atoms), INTENT(IN)    :: atoms
      TYPE(t_lapw), INTENT(IN)    :: lapw
      TYPE(t_usdus), INTENT(INOUT) :: usdus
      TYPE(t_potden), INTENT(IN)    :: vx
      TYPE(t_mat), INTENT(INOUT) :: hmat

      ! Scalar Arguments
      INTEGER, INTENT(IN) :: jsp

      ! Array Arguments
      REAL, INTENT(IN) :: vr0(:,:,:)               ! just for radial functions
      REAL, INTENT(IN) :: el(0:atoms%lmaxd, atoms%ntype, input%jspins)
      REAL, INTENT(IN) :: ello(:,:,:)
      REAL, INTENT(IN) :: bk(:)

      ! Local Scalars
      INTEGER               ::  ic, indx, m, ig1, ig2, n, nn
      INTEGER               ::  nlharm, nnbas, typsym, lm
      INTEGER               ::  noded, nodeu
      INTEGER               ::  nbasf0
      INTEGER               ::  i, j, l, ll, l1, l2, m1, m2, j1, j2
      INTEGER               ::  ok, p1, p2, lh, mh, pp1, pp2
      INTEGER               ::  igrid, itype, ilharm, istar
      INTEGER               ::  ineq, iatom, ilo, ilop, ieq, icentry
      INTEGER               ::  ikvecat, ikvecprevat, invsfct, ikvec, ikvecp
      INTEGER               ::  lp, mp, pp
      REAL                  ::  a_ex
      REAL                  ::  wronk
      COMPLEX               ::  rc, rr

      ! Local Arrays
      INTEGER               ::  gg(3)
      INTEGER               ::  pointer_lo(atoms%nlod, atoms%ntype)

      REAL                  ::  integ(0:sphhar%nlhd, maxval(mpbasis%num_radfun_per_l), 0:atoms%lmaxd, maxval(mpbasis%num_radfun_per_l), 0:atoms%lmaxd)
      REAL                  ::  grid(atoms%jmtd)
      REAL                  ::  vr(atoms%jmtd, 0:sphhar%nlhd)
      REAL                  ::  f(atoms%jmtd, 2, 0:atoms%lmaxd), g(atoms%jmtd, 2, 0:atoms%lmaxd)
      REAL                  ::  flo(atoms%jmtd, 2, atoms%nlod)
      REAL                  ::  uuilon(atoms%nlod, atoms%ntype), duilon(atoms%nlod, atoms%ntype)
      REAL                  ::  ulouilopn(atoms%nlod, atoms%nlod, atoms%ntype)

      REAL                  ::  bas1(atoms%jmtd, maxval(mpbasis%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype)
      REAL                  ::  bas2(atoms%jmtd, maxval(mpbasis%num_radfun_per_l), 0:atoms%lmaxd, atoms%ntype)

      COMPLEX               ::  vpw(stars%ng3)
      COMPLEX               ::  vxc(hmat%matsize1*(hmat%matsize1 + 1)/2)
      COMPLEX               ::  vrmat(hybrid%maxlmindx, hybrid%maxlmindx)
      COMPLEX               ::  carr(hybrid%maxlmindx, DIMENSION%nvd), carr1(DIMENSION%nvd, DIMENSION%nvd)
      COMPLEX, ALLOCATABLE  ::  ahlp(:, :, :), bhlp(:, :, :)
      COMPLEX, ALLOCATABLE  ::  bascof(:, :, :)
#ifndef CPP_OLDINTEL
      COMPLEX               ::  bascof_lo(3, -atoms%llod:atoms%llod, 4*atoms%llod + 2, atoms%nlod, atoms%nat)
#endif

      CALL timestart("subvxc")

      integ = 0.0
      bas1 = 0.0
      bas2 = 0.0

      vxc = 0

      ! Calculate radial functions
      mpbasis%num_radfun_per_l = 2
      DO itype = 1, atoms%ntype

         ! Generate the radial basis-functions for each l
         WRITE (6, '(a,i3,a)') new_LINE('n')//new_LINE('n')//' wavefunction parameters for atom type', itype, ':'
         WRITE (6, '(31x,a,32x,a)') 'radial function', 'energy derivative'
         WRITE (6, '(a)') '  l    energy            value        '// &
            'derivative    nodes          value        derivative    nodes       norm        wronskian'
         DO l = 0, atoms%lmax(itype)
            CALL radfun(l, itype, jsp, el(l, itype, jsp), vr0(:, itype, jsp), atoms, f(1, 1, l), g(1, 1, l), usdus, nodeu, noded, wronk)
            WRITE (6, FMT=8010) l, el(l, itype, jsp), usdus%us(l, itype, jsp), &
               usdus%dus(l, itype, jsp), nodeu, usdus%uds(l, itype, jsp), usdus%duds(l, itype, jsp), noded, &
               usdus%ddn(l, itype, jsp), wronk
         END DO
8010     FORMAT(i3, f10.5, 2(5x, 1p, 2e16.7, i5), 1p, 2e16.7)

         bas1(:, 1, :, itype) = f(:, 1, :)
         bas1(:, 2, :, itype) = g(:, 1, :)
         bas2(:, 1, :, itype) = f(:, 2, :)
         bas2(:, 2, :, itype) = g(:, 2, :)

         ! Generate the extra radial basis-functions for the local orbitals, if there are any.
         IF (atoms%nlo(itype) >= 1) THEN
            CALL radflo(atoms, itype, jsp, ello(:,:, jsp), vr0(:, itype, jsp), f, g, mpi, &
                        usdus, uuilon, duilon, ulouilopn, flo, .TRUE.)

            DO i = 1, atoms%nlo(itype)
               mpbasis%num_radfun_per_l(atoms%llo(i, itype), itype) = mpbasis%num_radfun_per_l(atoms%llo(i, itype), itype) + 1
               pointer_lo(i, itype) = mpbasis%num_radfun_per_l(atoms%llo(i, itype), itype)
               bas1(:, mpbasis%num_radfun_per_l(atoms%llo(i, itype), itype), atoms%llo(i, itype), itype) = flo(:, 1, i)
               bas2(:, mpbasis%num_radfun_per_l(atoms%llo(i, itype), itype), atoms%llo(i, itype), itype) = flo(:, 2, i)
            END DO
         END IF
      END DO

      ! Compute APW coefficients

      ! Calculate bascof
      allocate(ahlp(DIMENSION%nvd, 0:DIMENSION%lmd, atoms%nat), bhlp(DIMENSION%nvd, 0:DIMENSION%lmd, atoms%nat), stat=ok)
      IF (ok /= 0) call judft_error('subvxc: error in allocation of ahlp/bhlp')
#ifndef CPP_OLDINTEL
      CALL abcof3(input, atoms, sym, jsp, cell, bk, lapw, usdus, oneD, ahlp, bhlp, bascof_lo)
#endif
      allocate(bascof(DIMENSION%nvd, 2*(DIMENSION%lmd + 1), atoms%nat), stat=ok)
      IF (ok /= 0) call judft_error('subvxc: error in allocation of bascof')

      bascof = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            indx = 0
            DO l = 0, atoms%lmax(itype)
               ll = l*(l + 1)
               DO M = -l, l
                  lm = ll + M
                  DO i = 1, 2
                     indx = indx + 1
                     IF (i == 1) THEN
                        bascof(:, indx, ic) = ahlp(:, lm, ic)
                     ELSE IF (i == 2) THEN
                        bascof(:, indx, ic) = bhlp(:, lm, ic)
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      deallocate(ahlp, bhlp)

      ! Loop over atom types
      iatom = 0
      DO itype = 1, atoms%ntype

         typsym = atoms%ntypsy(SUM(atoms%neq(:itype - 1)) + 1)
         nlharm = sphhar%nlh(typsym)

         ! Calculate vxc = vtot - vcoul
         DO l = 0, nlharm
            DO i = 1, atoms%jri(itype)
               IF (l == 0) THEN
                  ! vr(i,0)= vrtot(i,0,itype)*sfp/rmsh(i,itype) - vrcou(i,0,itype,jsp)
                  vr(i, 0) = vx%mt(i, 0, itype, jsp)! * sfp_const / atoms%rmsh(i,itype)
               ELSE ! vxc = vtot - vcoul
                  ! vr(i,l)= vrtot(i,l,itype)-vrcou(i,l,itype,jsp)
                  vr(i, l) = vx%mt(i, l, itype, jsp)
               END IF
            END DO
         END DO

         ! Calculate MT contribution to vxc matrix elements
         ! Precompute auxiliary radial integrals
         DO ilharm = 0, nlharm
            i = 0
            DO l1 = 0, atoms%lmax(itype)
               DO p1 = 1, 2
                  i = i + 1
                  j = 0
                  DO l2 = 0, atoms%lmax(itype)
                     DO p2 = 1, 2
                        j = j + 1
                        IF (j <= i) THEN
                           DO igrid = 1, atoms%jri(itype)
                              grid(igrid) = vr(igrid, ilharm)*(bas1(igrid, p1, l1, itype)*bas1(igrid, p2, l2, itype) + &
                                                               bas2(igrid, p1, l1, itype)*bas2(igrid, p2, l2, itype))
                           END DO
                           CALL intgr3(grid, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), integ(ilharm, p1, l1, p2, l2))
                           integ(ilharm, p2, l2, p1, l1) = integ(ilharm, p1, l1, p2, l2)
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END DO

         ! Calculate muffin tin contribution to vxc matrix
         vrmat = 0

         j1 = 0
         DO l1 = 0, atoms%lmax(itype) ! loop: left basis function
            DO m1 = -l1, l1
               DO p1 = 1, 2
                  j1 = j1 + 1
                  j2 = 0
                  DO l2 = 0, atoms%lmax(itype) ! loop: right basis function
                     DO m2 = -l2, l2
                        DO p2 = 1, 2
                           j2 = j2 + 1
                           rr = 0
                           DO ilharm = 0, nlharm ! loop: lattice harmonics of vxc
                              l = sphhar%llh(ilharm, typsym)
                              DO i = 1, sphhar%nmem(ilharm, typsym)
                                 M = sphhar%mlh(i, ilharm, typsym)
                                 rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, l, l2, m1, M, m2, atoms%lmaxd)
                                 rr = rr + integ(ilharm, p1, l1, p2, l2)*rc
                              END DO
                           END DO
                           rc = CMPLX(0, 1)**(l2 - l1) ! adjusts to a/b/ccof-scaling
                           vrmat(j1, j2) = rr*rc
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         nnbas = j1

         ! Project on bascof
         DO ineq = 1, atoms%neq(itype)
            iatom = iatom + 1
            carr(:nnbas, :lapw%nv(jsp)) = CONJG(MATMUL(vrmat(:nnbas, :nnbas), &
                                                       TRANSPOSE(bascof(:lapw%nv(jsp), :nnbas, iatom))))
            carr1(:lapw%nv(jsp), :lapw%nv(jsp)) = MATMUL(bascof(:lapw%nv(jsp), :nnbas, iatom), carr(:nnbas, :lapw%nv(jsp)))
            ic = 0
            DO j = 1, lapw%nv(jsp)
               ! carr(:nnbas) =  matmul(vrmat(:nnbas,:nnbas),bascof(j,:nnbas,iatom) )
               DO i = 1, j
                  ic = ic + 1
                  vxc(ic) = vxc(ic) + carr1(i, j)
                  ! vxc(ic) = vxc(ic) + conjg(dot_product ( bascof(i,:nnbas,iatom),carr(:nnbas) ))
               END DO
            END DO
         END DO
      END DO ! End loop over atom types

      ! Calculate plane wave contribution
      DO i = 1, stars%ng3
         vpw(i) = vx%pw_w(i, jsp)
         ! vpw(i)=vpwtot(i)-vpwcou(i,jsp)
      END DO

      ! Calculate vxc-matrix,  left basis function (ig1)
      !                        right basis function (ig2)
      ic = 0
      DO ig1 = 1, lapw%nv(jsp)
         DO ig2 = 1, ig1
            ic = ic + 1
            gg = lapw%gvec(:,ig1,jsp) - lapw%gvec(:,ig2,jsp)
            istar = stars%ig(gg(1), gg(2), gg(3))
            IF (istar /= 0) THEN
               vxc(ic) = vxc(ic) + stars%rgphs(gg(1), gg(2), gg(3))*vpw(istar)
            ELSE
               IF (mpi%irank == 0) THEN
                  WRITE (6, '(A,/6I5)') 'Warning: Gi-Gj not in any star:', &
                     lapw%gvec(:,ig1, jsp), lapw%gvec(:,ig2, jsp)
               END IF
            END IF
         END DO
      END DO

      ! Calculate local orbital contribution
      IF (ANY(atoms%nlo /= 0)) THEN

         nbasf0 = lapw%nv(jsp)*(lapw%nv(jsp) + 1)/2  ! number of pure APW contributions
         icentry = nbasf0                          ! icentry counts the entry in the matrix vxc
         iatom = 0
         ikvecat = 0
         ikvecprevat = 0

         DO itype = 1, atoms%ntype

            typsym = atoms%ntypsy(SUM(atoms%neq(:itype - 1)) + 1)
            nlharm = sphhar%nlh(typsym)

            ! Calculate vxc = vtot - vcoul
            DO l = 0, nlharm
               DO i = 1, atoms%jri(itype)
                  IF (l == 0) THEN
                     ! vr(i,0)= vrtot(i,0,itype)*sfp/rmsh(i,itype) -  vrcou(i,0,itype,jsp)
                     vr(i, 0) = vx%mt(i, 0, itype, jsp) !*sfp_const/atoms%rmsh(i, itype)
                  ELSE ! vxc = vtot - vcoul
                     vr(i, l) = vx%mt(i, l, itype, jsp)                    !
                     ! vr(i,l)=  vrtot(i,l,itype)-vrcou(i,l,itype,jsp)
                  END IF
               END DO
            END DO

            ! Precompute auxiliary radial integrals
            DO ilharm = 0, nlharm
               i = 0
               DO l1 = 0, atoms%lmax(itype)
                  DO p1 = 1, mpbasis%num_radfun_per_l(l1, itype)
                     i = i + 1
                     j = 0
                     DO l2 = 0, atoms%lmax(itype)
                        DO p2 = 1, mpbasis%num_radfun_per_l(l2, itype)
                           j = j + 1
                           IF (j <= i) THEN
                              DO igrid = 1, atoms%jri(itype)
                                 grid(igrid) = vr(igrid, ilharm)*(bas1(igrid, p1, l1, itype)*bas1(igrid, p2, l2, itype) + &
                                                                  bas2(igrid, p1, l1, itype)*bas2(igrid, p2, l2, itype))
                              END DO
                              CALL intgr3(grid, atoms%rmsh(:, itype), atoms%dx(itype), atoms%jri(itype), integ(ilharm, p1, l1, p2, l2))
                              integ(ilharm, p2, l2, p1, l1) = integ(ilharm, p1, l1, p2, l2)
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               IF ((atoms%invsat(iatom) == 0) .OR. (atoms%invsat(iatom) == 1)) THEN

                  IF (atoms%invsat(iatom) == 0) invsfct = 1
                  IF (atoms%invsat(iatom) == 1) invsfct = 2

                  DO ilo = 1, atoms%nlo(itype)
#ifdef CPP_OLDINTEL
                     CALL judft_error("no LOs & hybrid with old intel compiler!", calledby="subvxc.F90")
#else
                     l1 = atoms%llo(ilo, itype)
                     DO ikvec = 1, invsfct*(2*l1 + 1)
                        DO m1 = -l1, l1
                           DO p1 = 1, 3
                              IF (p1 == 3) THEN
                                 pp1 = pointer_lo(ilo, itype)
                              ELSE
                                 pp1 = p1
                              END IF

                              IF (mpbasis%num_radfun_per_l(l1, itype) <= 2) call judft_error('subvxc: error mpbasis%num_radfun_per_l')

                              lm = 0

                              !loop over APW
                              DO l2 = 0, atoms%lmax(itype)
                                 DO m2 = -l2, l2
                                    DO p2 = 1, 2
                                       lm = lm + 1

                                       rr = 0
                                       DO ilharm = 0, nlharm
                                          lh = sphhar%llh(ilharm, typsym)
                                          DO i = 1, sphhar%nmem(ilharm, typsym)
                                             mh = sphhar%mlh(i, ilharm, typsym)
                                             rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, lh, l2, m1, mh, m2, atoms%lmaxd)
                                             rr = rr + integ(ilharm, p2, l2, pp1, l1)*rc
                                          END DO
                                       END DO

                                       rc = CMPLX(0.0, 1.0)**(l2 - l1) ! adjusts to a/b/ccof-scaling

                                       ! ic counts the entry in vxc
                                       ic = icentry
                                       DO i = 1, lapw%nv(jsp)
                                          ic = ic + 1
                                          IF (hmat%l_real) THEN
                                             vxc(ic) = vxc(ic) + invsfct*REAL(rr*rc*bascof(i, lm, iatom)* &
                                                                              CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom)))
                                          ELSE
                                             vxc(ic) = vxc(ic) + rr*rc*bascof(i, lm, iatom)* &
                                                       CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom))
                                          END IF
                                       END DO
                                    END DO  !p2
                                 END DO  ! m2
                              END DO ! l2 ->  loop over APW

                              ! calcualte matrix-elements with local orbitals at the same atom
                              IF (ic /= icentry + lapw%nv(jsp)) call judft_error('subvxc: error counting ic')

                              ic = ic + ikvecprevat

                              DO ilop = 1, ilo - 1
                                 lp = atoms%llo(ilop, itype)
                                 DO ikvecp = 1, invsfct*(2*lp + 1)
                                    ic = ic + 1
                                    DO mp = -lp, lp
                                       DO pp = 1, 3
                                          IF (pp == 3) THEN
                                             pp2 = pointer_lo(ilop, itype)
                                          ELSE
                                             pp2 = pp
                                          END IF

                                          rr = 0
                                          DO ilharm = 0, nlharm
                                             lh = sphhar%llh(ilharm, typsym)
                                             DO i = 1, sphhar%nmem(ilharm, typsym)
                                                mh = sphhar%mlh(i, ilharm, typsym)
                                                rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, lh, lp, m1, mh, mp, atoms%lmaxd)
                                                rr = rr + integ(ilharm, pp2, lp, pp1, l1)*rc
                                             END DO
                                          END DO

                                          rc = CMPLX(0.0, 1.0)**(lp - l1) ! adjusts to a/b/ccof-scaling

                                          IF (hmat%l_real) THEN
                                             vxc(ic) = vxc(ic) + invsfct*REAL(rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                                              CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom)))
                                          ELSE
                                             vxc(ic) = vxc(ic) + rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                       CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom))
                                          END IF
                                       END DO ! pp
                                    END DO ! mp

                                 END DO !ikvecp
                              END DO ! ilop

                              ! calculate matrix-elements of one local orbital with itself
                              DO ikvecp = 1, ikvec
                                 ic = ic + 1

                                 lp = l1
                                 ilop = ilo
                                 DO mp = -lp, lp
                                    DO pp = 1, 3
                                       IF (pp == 3) THEN
                                          pp2 = pointer_lo(ilop, itype)
                                       ELSE
                                          pp2 = pp
                                       END IF

                                       rr = 0
                                       DO ilharm = 0, nlharm
                                          lh = sphhar%llh(ilharm, typsym)
                                          DO i = 1, sphhar%nmem(ilharm, typsym)
                                             mh = sphhar%mlh(i, ilharm, typsym)
                                             rc = sphhar%clnu(i, ilharm, typsym)*gaunt1(l1, lh, lp, m1, mh, mp, atoms%lmaxd)
                                             rr = rr + integ(ilharm, pp2, lp, pp1, l1)*rc
                                          END DO
                                       END DO

                                       rc = CMPLX(0.0, 1.0)**(lp - l1) ! adjusts to a/b/ccof-scaling

                                       IF (hmat%l_real) THEN
                                          vxc(ic) = vxc(ic) + invsfct*REAL(rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                                           CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom)))
                                       ELSE
                                          vxc(ic) = vxc(ic) + rr*rc*bascof_lo(pp, mp, ikvecp, ilop, iatom)* &
                                                    CONJG(bascof_lo(p1, m1, ikvec, ilo, iatom))
                                       END IF
                                    END DO ! pp
                                 END DO ! mp
                              END DO ! ikvecp
                           END DO  ! p1
                        END DO  ! m1
                        icentry = ic
                     END DO !ikvec
                     ikvecat = ikvecat + invsfct*(2*l1 + 1)
#endif
                  END DO  ! ilo
                  ikvecprevat = ikvecprevat + ikvecat
                  ikvecat = 0
               END IF ! atoms%invsat(iatom)
            END DO ! ieq
         END DO !itype
      END IF ! if any atoms%llo

      !initialize weighting factor
      a_ex = xcpot%get_exchange_weight()

      i = 0
      DO n = 1, hmat%matsize1
         DO nn = 1, n
            i = i + 1
            IF (hmat%l_real) THEN
               hmat%data_r(nn, n) = hmat%data_r(nn, n) - a_ex*REAL(vxc(i))
            ELSE
               hmat%data_c(nn, n) = hmat%data_c(nn, n) - a_ex*vxc(i)
            ENDIF
         END DO
      END DO

      CALL timestop("subvxc")

      deallocate(bascof)

   END SUBROUTINE subvxc
END MODULE m_subvxc
