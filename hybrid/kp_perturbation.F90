MODULE m_kp_perturbation

CONTAINS

   SUBROUTINE ibs_correction( &
      nk, atoms, &
      dimension, input, jsp, &
      hybdat, hybrid, &
      lapw, kpts, nkpti, &
      cell, mnobd, &
      sym, &
      proj_ibsc, olap_ibsc)

      USE m_sphbes
      USE m_dsphbs
      USE m_constants
      USE m_ylm
      USE m_gaunt
      USE m_util
      USE m_types
      USE m_io_hybrid
      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(IN)   :: hybdat
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_hybrid), INTENT(INOUT)   :: hybrid
      TYPE(t_input), INTENT(IN)   :: input
      TYPE(t_sym), INTENT(IN)   :: sym
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_lapw), INTENT(IN)   :: lapw

      ! - scalars -
      INTEGER, INTENT(IN)   :: jsp
      INTEGER, INTENT(IN)   ::  mnobd
      INTEGER, INTENT(IN)   ::  nk, nkpti

      ! - arrays -

      COMPLEX, INTENT(INOUT)::  olap_ibsc(:,:,:,:)
      COMPLEX, INTENT(INOUT)::  proj_ibsc(:, :, :)!(3,mnobd,hybrid%nbands(nk))
      ! - local scalars -
      INTEGER               ::  i, itype, ieq, iatom, iatom1, iband, iband1
      INTEGER               ::  iband2, ilo, ibas, ic, ikpt, ikvec, invsfct
      INTEGER               ::  irecl_cmt, irecl_z
      INTEGER               ::  j, m
      INTEGER               ::  l1, m1, p1, l2, m2, p2, l, p, lm, &
                               lmp, lmp1, lmp2, lm1, lm2
      INTEGER               ::  ok, ig
      INTEGER               ::  idum
      REAL                  ::  const
      REAL                  ::  ka, kb
      REAL                  ::  kvecn
      REAL                  ::  olap_udot, olap_uulo, olap_udotulo
      REAL                  ::  rdum
      REAL                  ::  ws
      COMPLEX               ::  phase
      COMPLEX               ::  cj, cdj
      COMPLEX               ::  denom, enum
      COMPLEX               ::  cdum, cdum1, cdum2
      COMPLEX, PARAMETER    ::  img = (0.0, 1.0)
      ! - local arrays -
      INTEGER               ::  lmp_start(atoms%ntype)
      REAL                  ::  alo(atoms%nlod, atoms%ntype), blo(atoms%nlod, atoms%ntype), &
                               clo(atoms%nlod, atoms%ntype)
      REAL                  ::  u1_lo(atoms%jmtd, atoms%nlod, atoms%ntype), &
                               u2_lo(atoms%jmtd, atoms%nlod, atoms%ntype)
      REAL                  ::  kvec(3), qvec(3)
      REAL                  ::  sbes(0:atoms%lmaxd + 1), dsbes(0:atoms%lmaxd + 1)
      REAL                  ::  bas1_tmp(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd + 1, atoms%ntype), &
                               bas2_tmp(atoms%jmtd, hybrid%maxindx, 0:atoms%lmaxd + 1, atoms%ntype)
      REAL                  ::  bas1_MT_tmp(hybrid%maxindx, 0:atoms%lmaxd + 1, atoms%ntype), &
                               drbas1_MT_tmp(hybrid%maxindx, 0:atoms%lmaxd + 1, atoms%ntype)
      REAL                  ::  ru1(atoms%jmtd, 3, mnobd), ru2(atoms%jmtd, 3, mnobd)
      REAL                  ::  iu1(atoms%jmtd, 3, mnobd), iu2(atoms%jmtd, 3, mnobd)
      REAL                  ::  rintegrand(atoms%jmtd), iintegrand(atoms%jmtd), &
                               integrand(atoms%jmtd)

      COMPLEX               ::  f(atoms%jmtd, mnobd)
      COMPLEX               ::  carr(3), carr2(3, hybrid%nbands(nk))
      COMPLEX               ::  ylm((atoms%lmaxd + 2)**2)
      COMPLEX, ALLOCATABLE   ::  u1(:, :, :, :, :), u2(:, :, :, :, :)
      COMPLEX, ALLOCATABLE   ::  cmt_lo(:, :, :, :)
      COMPLEX, ALLOCATABLE   ::  cmt_apw(:, :, :)
      TYPE(t_mat)           ::  z
      REAL                  ::  work_r(dimension%neigd)
      COMPLEX               ::  work_c(dimension%neigd)

      !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)

      bas1_tmp(:, :, 0:atoms%lmaxd, :) = hybdat%bas1(:, :, 0:atoms%lmaxd, :)
      bas2_tmp(:, :, 0:atoms%lmaxd, :) = hybdat%bas2(:, :, 0:atoms%lmaxd, :)

      bas1_MT_tmp(:, 0:atoms%lmaxd, :) = hybdat%bas1_MT(:, 0:atoms%lmaxd, :)
      drbas1_MT_tmp(:, 0:atoms%lmaxd, :) = hybdat%drbas1_MT(:, 0:atoms%lmaxd, :)

      bas1_tmp(:, :, atoms%lmaxd + 1, :) = hybdat%bas1(:, :, atoms%lmaxd, :)
      bas2_tmp(:, :, atoms%lmaxd + 1, :) = hybdat%bas2(:, :, atoms%lmaxd, :)

      bas1_MT_tmp(:, atoms%lmaxd + 1, :) = hybdat%bas1_MT(:, atoms%lmaxd, :)
      drbas1_MT_tmp(:, atoms%lmaxd + 1, :) = hybdat%drbas1_MT(:, atoms%lmaxd, :)

      ! read in z coefficient from direct access file z at k-point nk

      call read_z(z, nk)

      ! construct local orbital consisting of radial function times spherical harmonic
      ! where the radial function vanishes on the MT sphere boundary
      ! with this the local orbitals have a trivial k-dependence

      ! compute radial lo matching coefficients
      hybrid%num_radfun_per_l = 2
      DO itype = 1, atoms%ntype
         DO ilo = 1, atoms%nlo(itype)
            l = atoms%llo(ilo, itype)
            hybrid%num_radfun_per_l(l, itype) = hybrid%num_radfun_per_l(l, itype) + 1
            p = hybrid%num_radfun_per_l(l, itype)

            ws = -wronskian(hybdat%bas1_MT(1, l, itype), hybdat%drbas1_MT(1, l, itype), hybdat%bas1_MT(2, l, itype), hybdat%drbas1_MT(2, l, itype))

            ka = 1.0/ws*wronskian(hybdat%bas1_MT(p, l, itype), hybdat%drbas1_MT(p, l, itype), hybdat%bas1_MT(2, l, itype), hybdat%drbas1_MT(2, l, itype))

            kb = 1.0/ws*wronskian(hybdat%bas1_MT(1, l, itype), hybdat%drbas1_MT(1, l, itype), hybdat%bas1_MT(p, l, itype), hybdat%drbas1_MT(p, l, itype))

            integrand = hybdat%bas1(:, 2, l, itype)*hybdat%bas1(:, 2, l, itype) + hybdat%bas2(:, 2, l, itype)*hybdat%bas2(:, 2, l, itype)
            olap_udot = intgrf(integrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)

            integrand = hybdat%bas1(:, 1, l, itype)*hybdat%bas1(:, p, l, itype) + hybdat%bas2(:, 1, l, itype)*hybdat%bas2(:, p, l, itype)
            olap_uulo = intgrf(integrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)

            integrand = hybdat%bas1(:, 2, l, itype)*hybdat%bas1(:, p, l, itype) + hybdat%bas2(:, 2, l, itype)*hybdat%bas2(:, p, l, itype)
            olap_udotulo = intgrf(integrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)

            rdum = ka**2 + (kb**2)*olap_udot + 1.0 + 2.0*ka*olap_uulo + 2.0*kb*olap_udotulo
            clo(ilo, itype) = 1.0/sqrt(rdum)
            alo(ilo, itype) = ka*clo(ilo, itype)
            blo(ilo, itype) = kb*clo(ilo, itype)

            u1_lo(:, ilo, itype) = alo(ilo, itype)*hybdat%bas1(:, 1, l, itype) + blo(ilo, itype)*hybdat%bas1(:, 2, l, itype) + clo(ilo, itype)*hybdat%bas1(:, p, l, itype)

            u2_lo(:, ilo, itype) = alo(ilo, itype)*hybdat%bas2(:, 1, l, itype) + blo(ilo, itype)*hybdat%bas2(:, 2, l, itype) + clo(ilo, itype)*hybdat%bas2(:, p, l, itype)
         END DO
      END DO

      ! calculate lo wavefunction coefficients
      allocate(cmt_lo(dimension%neigd, -atoms%llod:atoms%llod, atoms%nlod, atoms%nat))
      cmt_lo = 0
      iatom = 0
      ic = 0
      ibas = lapw%nv(jsp)
      DO itype = 1, atoms%ntype
         ! the program is in hartree units, therefore 1/wronskian is
         ! (rmt**2)/2.
         const = fpi_const*(atoms%rmt(itype)**2)/2/sqrt(cell%omtil)
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            IF ((atoms%invsat(iatom) == 0) .or. (atoms%invsat(iatom) == 1)) THEN
               IF (atoms%invsat(iatom) == 0) invsfct = 1
               IF (atoms%invsat(iatom) == 1) THEN
                  invsfct = 2
                  iatom1 = sym%invsatnr(iatom)
               END IF

               DO ilo = 1, atoms%nlo(itype)
                  l = atoms%llo(ilo, itype)
                  cdum = img**l*const
                  DO ikvec = 1, invsfct*(2*l + 1)
                     ic = ic + 1
                     ibas = ibas + 1
                     kvec = kpts%bk(:, nk) + lapw%gvec(:, hybdat%kveclo_eig(ic, nk), jsp)

                     phase = exp(img*tpi_const*dot_product(atoms%taual(:, iatom), kvec))
                     cdum1 = cdum*phase

                     CALL ylm4(l, matmul(kvec, cell%bmat), ylm)

                     lm = l**2
                     DO M = -l, l
                        lm = lm + 1
                        cdum2 = cdum1*conjg(ylm(lm))
                        if (z%l_real) THEN
                           work_r = z%data_r(ibas, :)
                           DO iband = 1, hybrid%nbands(nk)
                              cmt_lo(iband, M, ilo, iatom) = cmt_lo(iband, M, ilo, iatom) + cdum2*work_r(iband)
                              IF (invsfct == 2) THEN
                                 ! the factor (-1)**l is necessary as we do not calculate
                                 ! the cmt_lo in the local coordinate system of the atom
                                 cmt_lo(iband, -M, ilo, iatom1) = cmt_lo(iband, -M, ilo, iatom1) + (-1)**(l + M)*conjg(cdum2)*work_r(iband)
                              END IF
                           END DO
                        else
                           work_c = z%data_c(ibas, :)
                           DO iband = 1, hybrid%nbands(nk)
                              cmt_lo(iband, M, ilo, iatom) = cmt_lo(iband, M, ilo, iatom) + cdum2*work_c(iband)
                              IF (invsfct == 2) THEN
                                 ! the factor (-1)**l is necessary as we do not calculate
                                 ! the cmt_lo in the local coordinate system of the atom
                                 cmt_lo(iband, -M, ilo, iatom1) = cmt_lo(iband, -M, ilo, iatom1) + (-1)**(l + M)*conjg(cdum2)*work_c(iband)
                              END IF
                           END DO
                        end if

                     END DO

                  END DO  !ikvec
               END DO  ! ilo

            END IF

         END DO  !ieq
      END DO  !itype

      !
      ! calculate apw wavefunction coefficients up to lmax + 1
      ! note that the lo contribution is separated in cmt_lo
      !

      DO itype = 1, atoms%ntype
         lmp_start(itype) = sum((/(2*(2*l + 1), l=0, atoms%lmax(itype) + 1)/))
      END DO
      idum = maxval(lmp_start)

      allocate(cmt_apw(dimension%neigd, idum, atoms%nat))
      cmt_apw = 0
      DO i = 1, lapw%nv(jsp)
         kvec = kpts%bk(:, nk) + lapw%gvec(:, i, jsp)
         kvecn = sqrt(dot_product(matmul(kvec, cell%bmat), matmul(kvec, cell%bmat)))

         iatom = 0
         DO itype = 1, atoms%ntype
            !calculate spherical sperical harmonics
            CALL ylm4(atoms%lmax(itype) + 1, matmul(kvec, cell%bmat), ylm)

            !calculate spherical bessel function at |kvec|*R_MT(itype)
            CALL sphbes(atoms%lmax(itype) + 1, kvecn*atoms%rmt(itype), sbes)

            !calculate radial derivative of spherical bessel function at |kvec|*R_MT(itype)
            CALL dsphbs(atoms%lmax(itype) + 1, kvecn*atoms%rmt(itype), sbes, dsbes)
            dsbes = kvecn*dsbes

            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1

               phase = exp(img*tpi_const*dot_product(kvec, atoms%taual(:, iatom)))

               lm = 0
               lmp = 0
               DO l = 0, atoms%lmax(itype) + 1
                  denom = wronskian(bas1_MT_tmp(2, l, itype), drbas1_MT_tmp(2, l, itype), &
                                    bas1_MT_tmp(1, l, itype), drbas1_MT_tmp(1, l, itype))
                  cdum1 = fpi_const*img**l*sbes(l)*phase/sqrt(cell%omtil)
                  cdum2 = fpi_const*img**l*dsbes(l)*phase/sqrt(cell%omtil)
                  DO M = -l, l
                     lm = lm + 1
                     cj = cdum1*conjg(ylm(lm))
                     cdj = cdum2*conjg(ylm(lm))
                     DO p = 1, 2
                        lmp = lmp + 1
                        p1 = p + (-1)**(p - 1)

                        enum = CMPLX(wronskian(bas1_MT_tmp(p1, l, itype), drbas1_MT_tmp(p1, l, itype), REAL(cj), REAL(cdj)), &
                                     wronskian(bas1_MT_tmp(p1, l, itype), drbas1_MT_tmp(p1, l, itype), AIMAG(cj), AIMAG(cdj)))

                        cdum = (-1)**(p + 1)*enum/denom
                        if (z%l_real) THEN
                           work_r = z%data_r(i, :)
                           DO iband = 1, hybrid%nbands(nk)
                              cmt_apw(iband, lmp, iatom) = cmt_apw(iband, lmp, iatom) + cdum*work_r(iband)
                           END DO
                        else
                           work_c = z%data_c(i, :)
                           DO iband = 1, hybrid%nbands(nk)
                              cmt_apw(iband, lmp, iatom) = cmt_apw(iband, lmp, iatom) + cdum*work_c(iband)
                           END DO
                        end if
                     END DO  !p
                  END DO  !M
               END DO  !l

            END DO  !iatom
         END DO  ! itype
      END DO  ! i

      ! construct radial functions (complex) for the first order
      ! incomplete basis set correction

      allocate(u1(atoms%jmtd, 3, mnobd, (atoms%lmaxd + 1)**2, atoms%nat), stat=ok)!hybrid%nbands
      IF (ok /= 0) call judft_error('kp_perturbation: failure allocation u1')
      allocate(u2(atoms%jmtd, 3, mnobd, (atoms%lmaxd + 1)**2, atoms%nat), stat=ok)!hybrid%nbands
      IF (ok /= 0) call judft_error('kp_perturbation: failure allocation u2')
      u1 = 0; u2 = 0

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1

            lm1 = 0
            DO l1 = 0, atoms%lmax(itype)! + 1
               DO m1 = -l1, l1
                  lm1 = lm1 + 1

                  DO p1 = 1, 2

                     carr2 = 0
                     l2 = l1 + 1
                     lmp2 = 2*l2**2 + p1
                     DO m2 = -l2, l2
                        carr = gauntvec(l1, m1, l2, m2, atoms)
                        DO iband = 1, mnobd! hybrid%nbands
                           carr2(1:3, iband) = carr2(1:3, iband) + carr*cmt_apw(iband, lmp2, iatom)
                        END DO
                        lmp2 = lmp2 + 2
                     END DO

                     DO iband = 1, mnobd! hybrid%nbands
                        DO i = 1, 3
                           DO ig = 1, atoms%jri(itype)
                              ! the r factor is already included in bas1
                              u1(ig, i, iband, lm1, iatom) = u1(ig, i, iband, lm1, iatom) - img*bas1_tmp(ig, p1, l2, itype)*carr2(i, iband)
                              u2(ig, i, iband, lm1, iatom) = u2(ig, i, iband, lm1, iatom) - img*bas2_tmp(ig, p1, l2, itype)*carr2(i, iband)
                           END DO
                        END DO
                     END DO

                     l2 = l1 - 1
                     IF (l2 >= 0) THEN
                        carr2 = 0
                        lmp2 = 2*l2**2 + p1
                        DO m2 = -l2, l2
                           carr = gauntvec(l1, m1, l2, m2, atoms)
                           DO iband = 1, mnobd! hybrid%nbands
                              carr2(1:3, iband) = carr2(1:3, iband) + carr*cmt_apw(iband, lmp2, iatom)
                           END DO
                           lmp2 = lmp2 + 2
                        END DO

                        DO iband = 1, mnobd! hybrid%nbands
                           DO i = 1, 3
                              DO ig = 1, atoms%jri(itype)
                                 ! the r factor is already included in bas1
                                 u1(ig, i, iband, lm1, iatom) = u1(ig, i, iband, lm1, iatom) - img*bas1_tmp(ig, p1, l2, itype)*carr2(i, iband)
                                 u2(ig, i, iband, lm1, iatom) = u2(ig, i, iband, lm1, iatom) - img*bas2_tmp(ig, p1, l2, itype)*carr2(i, iband)
                              END DO
                           END DO
                        END DO

                     END IF

                     carr2 = 0
                     l2 = l1 + 1
                     lmp2 = 2*l2**2
                     DO m2 = -l2, l2
                        carr = gauntvec(l1, m1, l2, m2, atoms)
                        DO p2 = 1, 2
                           lmp2 = lmp2 + 1
                           rdum = w(p1, l1, p2, l2, itype, bas1_MT_tmp, drbas1_MT_tmp, atoms%rmt)
                           DO iband = 1, mnobd! hybrid%nbands
                              carr2(1:3, iband) = carr2(1:3, iband) + img*carr*rdum*cmt_apw(iband, lmp2, iatom)
                           END DO
                        END DO
                     END DO

                     DO iband = 1, mnobd! hybrid%nbands
                        DO i = 1, 3
                           DO ig = 1, atoms%jri(itype)
                              u1(ig, i, iband, lm1, iatom) = u1(ig, i, iband, lm1, iatom) + bas1_tmp(ig, p1, l1, itype)*carr2(i, iband)/atoms%rmsh(ig, itype)
                              u2(ig, i, iband, lm1, iatom) = u2(ig, i, iband, lm1, iatom) + bas2_tmp(ig, p1, l1, itype)*carr2(i, iband)/atoms%rmsh(ig, itype)
                           END DO
                        END DO
                     END DO

                     l2 = l1 - 1
                     IF (l2 >= 0) THEN
                        carr2 = 0
                        lmp2 = 2*l2**2
                        DO m2 = -l2, l2
                           carr = gauntvec(l1, m1, l2, m2, atoms)
                           DO p2 = 1, 2
                              lmp2 = lmp2 + 1
                              rdum = w(p1, l1, p2, l2, itype, bas1_MT_tmp, drbas1_MT_tmp, atoms%rmt)
                              DO iband = 1, mnobd! hybrid%nbands
                                 carr2(1:3, iband) = carr2(1:3, iband) + img*carr*rdum*cmt_apw(iband, lmp2, iatom)
                              END DO
                           END DO
                        END DO

                        DO iband = 1, mnobd! hybrid%nbands
                           DO i = 1, 3
                              DO ig = 1, atoms%jri(itype)
                                 u1(ig, i, iband, lm1, iatom) = u1(ig, i, iband, lm1, iatom) &
                                                                + bas1_tmp(ig, p1, l1, itype)*carr2(i, iband)/atoms%rmsh(ig, itype)
                                 u2(ig, i, iband, lm1, iatom) = u2(ig, i, iband, lm1, iatom) &
                                                                + bas2_tmp(ig, p1, l1, itype)*carr2(i, iband)/atoms%rmsh(ig, itype)
                              END DO
                           END DO
                        END DO
                     END IF

                  END DO  ! p1
               END DO  !m1
            END DO  !l1
         END DO  !ieq
      END DO  !iatom

      ! construct lo contribtution
      IF (any(atoms%llo == atoms%lmaxd)) call judft_error('ibs_correction: atoms%llo=atoms%lmaxd is not implemented')

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            hybrid%num_radfun_per_l = 2
            iatom = iatom + 1
            DO ilo = 1, atoms%nlo(itype)
               l1 = atoms%llo(ilo, itype)
               hybrid%num_radfun_per_l(l1, itype) = hybrid%num_radfun_per_l(l1, itype) + 1
               p1 = hybrid%num_radfun_per_l(l1, itype)

               l2 = l1 + 1
               lm2 = l2**2
               DO m2 = -l2, l2
                  lm2 = lm2 + 1
                  carr2 = 0

                  DO m1 = -l1, l1
                     carr = gauntvec(l2, m2, l1, m1, atoms)
                     DO iband = 1, mnobd
                        carr2(1:3, iband) = carr2(1:3, iband) + cmt_lo(iband, m1, ilo, iatom)*carr
                     END DO
                  END DO

                  DO iband = 1, mnobd
                     DO i = 1, 3
                        DO ig = 1, atoms%jri(itype)
                           ! the r factor is already included in
                           u1(ig, i, iband, lm2, iatom) = u1(ig, i, iband, lm2, iatom) - img*u1_lo(ig, ilo, itype)*carr2(i, iband)
                           u2(ig, i, iband, lm2, iatom) = u2(ig, i, iband, lm2, iatom) - img*u2_lo(ig, ilo, itype)*carr2(i, iband)
                        END DO
                     END DO
                  END DO

               END DO

               l2 = l1 - 1
               IF (l2 >= 0) THEN
                  lm2 = l2**2
                  DO m2 = -l2, l2
                     lm2 = lm2 + 1
                     carr2 = 0

                     DO m1 = -l1, l1
                        carr = gauntvec(l2, m2, l1, m1, atoms)
                        DO iband = 1, mnobd
                           carr2(1:3, iband) = carr2(1:3, iband) + cmt_lo(iband, m1, ilo, iatom)*carr
                        END DO
                     END DO

                     DO iband = 1, mnobd
                        DO i = 1, 3
                           DO ig = 1, atoms%jri(itype)
                              ! the r factor is already included in
                              u1(ig, i, iband, lm2, iatom) = u1(ig, i, iband, lm2, iatom) - img*u1_lo(ig, ilo, itype)*carr2(i, iband)
                              u2(ig, i, iband, lm2, iatom) = u2(ig, i, iband, lm2, iatom) - img*u2_lo(ig, ilo, itype)*carr2(i, iband)
                           END DO
                        END DO
                     END DO

                  END DO
               END IF

            END DO
         END DO
      END DO

      !
      ! calculate projection < phi(n',k)|phi^1(n,k)>
      ! n' = iband1 , n= iband2
      !
      iatom = 0
      proj_ibsc = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            lmp = 0
            DO l = 0, atoms%lmax(itype)
               DO M = -l, l
                  lm = lm + 1
                  ru1 = real(u1(:, :, :, lm, iatom))
                  iu1 = aimag(u1(:, :, :, lm, iatom))
                  ru2 = real(u2(:, :, :, lm, iatom))
                  iu2 = aimag(u2(:, :, :, lm, iatom))
                  DO p = 1, 2
                     lmp = lmp + 1

                     DO iband = 1, mnobd! hybrid%nbands
                        DO i = 1, 3

                           rintegrand = atoms%rmsh(:, itype)*(hybdat%bas1(:, p, l, itype)*ru1(:, i, iband) + hybdat%bas2(:, p, l, itype)*ru2(:, i, iband))

                           iintegrand = atoms%rmsh(:, itype)*(hybdat%bas1(:, p, l, itype)*iu1(:, i, iband) + hybdat%bas2(:, p, l, itype)*iu2(:, i, iband))

                           carr2(i, iband) = intgrf(rintegrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf) &
                                             + img*intgrf(iintegrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)

                        END DO
                     END DO

                     DO iband1 = 1, hybrid%nbands(nk)
                        cdum = conjg(cmt_apw(iband1, lmp, iatom))
                        DO iband2 = 1, mnobd! hybrid%nbands
                           proj_ibsc(1:3, iband2, iband1) = proj_ibsc(1:3, iband2, iband1) + cdum*carr2(1:3, iband2)
                        END DO
                     END DO

                  END DO!p
               END DO!M
            END DO!l

         END DO!ieq
      END DO!itype

      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO ilo = 1, atoms%nlo(itype)
               l = atoms%llo(ilo, itype)
               lm = l**2
               DO M = -l, l
                  lm = lm + 1
                  ru1 = real(u1(:, :, :, lm, iatom))
                  iu1 = aimag(u1(:, :, :, lm, iatom))
                  ru2 = real(u2(:, :, :, lm, iatom))
                  iu2 = aimag(u2(:, :, :, lm, iatom))

                  DO iband = 1, mnobd! hybrid%nbands
                     DO i = 1, 3

                        rintegrand = atoms%rmsh(:, itype)*(u1_lo(:, ilo, itype)*ru1(:, i, iband) + u2_lo(:, ilo, itype)*ru2(:, i, iband))

                        iintegrand = atoms%rmsh(:, itype)*(u1_lo(:, ilo, itype)*iu1(:, i, iband) + u2_lo(:, ilo, itype)*iu2(:, i, iband))

                        carr2(i, iband) = intgrf(rintegrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf) &
                                          + img*intgrf(iintegrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)

                     END DO
                  END DO

                  DO iband1 = 1, hybrid%nbands(nk)
                     cdum = conjg(cmt_lo(iband1, M, ilo, iatom))
                     DO iband2 = 1, mnobd! hybrid%nbands
                        proj_ibsc(1:3, iband2, iband1) = proj_ibsc(1:3, iband2, iband1) + cdum*carr2(1:3, iband2)
                     END DO
                  END DO

               END DO

            END DO
         END DO
      END DO

      !
      ! calculate <phi^1(n1,k)|phi^1(n2,k)>
      ! n1 and n2 occupied
      !
      iatom = 0
      olap_ibsc = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            lm = 0
            DO l = 0, atoms%lmax(itype)!+1
               DO M = -l, l
                  lm = lm + 1
                  ru1 = real(u1(:, :, :, lm, iatom))
                  iu1 = aimag(u1(:, :, :, lm, iatom))
                  ru2 = real(u2(:, :, :, lm, iatom))
                  iu2 = aimag(u2(:, :, :, lm, iatom))

                  DO iband1 = 1, mnobd ! hybrid%nbands
                     DO iband2 = 1, mnobd!iband1
                        DO i = 1, 3
                           DO j = 1, 3

                              rintegrand = atoms%rmsh(:, itype)**2*(ru1(:, i, iband1)*ru1(:, j, iband2) + ru2(:, i, iband1)*ru2(:, j, iband2) &
                                                                    + iu1(:, i, iband1)*iu1(:, j, iband2) + iu2(:, i, iband1)*iu2(:, j, iband2))

                              iintegrand = atoms%rmsh(:, itype)**2*(ru1(:, i, iband1)*iu1(:, j, iband2) + ru2(:, i, iband1)*iu2(:, j, iband2) &
                                                                    - iu1(:, i, iband1)*ru1(:, j, iband2) - iu2(:, i, iband1)*ru2(:, j, iband2))

                              olap_ibsc(i, j, iband2, iband1) = olap_ibsc(i, j, iband2, iband1) &
                                                                + intgrf(rintegrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf) &
                                                                + img*intgrf(iintegrand, atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)

                           END DO
                        END DO

                     END DO
                  END DO

               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE ibs_correction

   FUNCTION gauntvec(l1, m1, l2, m2, atoms)

      USE m_constants
      USE m_gaunt
      USE m_juDFT

      USE m_types
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN)   :: atoms

      INTEGER, INTENT(IN)  ::  l1, m1, l2, m2

      COMPLEX               ::  gauntvec(-1:1)

      INTEGER               ::  j, mj
      INTEGER               ::  point(-1:1)
      REAL                  ::  rfac
      COMPLEX               ::  cfac
      COMPLEX, PARAMETER   ::  img = (0.0, 1.0)

      rfac = sqrt(tpi_const/3)
      DO j = -1, 1
         ! j = -1 corresponds to x-direction
         ! j = 0  corresponds to z-direction
         ! j = 1  corresponds to y-direction
         mj = abs(j)
         cfac = img**((j + mj)/2.)*sqrt(2.)**(1 - mj)*rfac
         gauntvec(j) = cfac*(gaunt1(1, l1, l2, -mj, m1, m2, atoms%lmaxd + 1) + j*gaunt1(1, l1, l2, mj, m1, m2, atoms%lmaxd + 1))
      END DO

      ! transform onto cartesian coordinates
      point(-1) = -1
      point(0) = 1
      point(1) = 0

      gauntvec = gauntvec(point)

   END FUNCTION gauntvec

   FUNCTION w(p1, l1, p2, l2, itype, bas1_mt, drbas1_mt, &
              rmt)
      USE m_types
      USE m_juDFT
      IMPLICIT NONE

      INTEGER, INTENT(IN)       ::  p1, l1, p2, l2
      INTEGER, INTENT(IN)       ::  itype
      REAL, INTENT(IN)            :: rmt(:), bas1_mt(:, 0:, :), drbas1_mt(:, 0:, :)

      REAL                  ::  w

      INTEGER               ::  p
      REAL                  ::  denom, enum

      IF (p1 > 2 .or. p2 > 2) call judft_error('w: the formalism is only valid for p<=2')

      denom = wronskian(bas1_MT(2, l1, itype), drbas1_MT(2, l1, itype), bas1_MT(1, l1, itype), drbas1_MT(1, l1, itype))

      p = p1 + (-1)**(p1 - 1)

      enum = bas1_MT(p, l1, itype)*bas1_MT(p2, l2, itype) + rmt(itype)*wronskian(bas1_MT(p, l1, itype), &
                                                                                 drbas1_MT(p, l1, itype), bas1_MT(p2, l2, itype), drbas1_MT(p2, l2, itype))

      w = (-1)**(p1 + 1)*enum/denom

   END FUNCTION

   PURE FUNCTION wronskian(f, df, g, dg)
      IMPLICIT NONE
      REAL, INTENT(IN) ::  f, df, g, dg
      REAL              ::  wronskian

      wronskian = f*dg - df*g

   END FUNCTION

!     Calculates the derivative
!                                  ikr    s       s
!     dcprod(n',n,k,xyz) = d    < e    phi   | phi       > / sqrt(vol)
!                           xyz           qn      q+k,n'
!
!                                s             s
!                           < phi  | d    | phi    >
!                                nq   xyz      n'q
!                      = -i ------------------------ / sqrt(vol)  ,   s = ispin
!                                  s     s                        ,   n = occ.     ,   n' = unocc.
!                                 e    - e                        ,   bandi1 <= n <= bandf1 , bandi2 <= n' <= bandf2
!                                  n'q    nq
!
!     with kp perturbation theory and
!
!               d     d     d
!     d    =   ---,  ---,  --- .
!      xyz     dk    dk    dk
!                x     y     z
!
   SUBROUTINE dwavefproducts( &
      dcprod, nk, bandi1, bandf1, bandi2, bandf2, lwrite, &
      atoms, hybrid, &
      cell, &
      hybdat, kpts, nkpti, lapw, &
      dimension, jsp, &
      eig_irr)

      USE m_wrapper
      USE m_types
      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(IN)   :: hybdat

      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_lapw), INTENT(IN)   :: lapw

!     - scalars -
      INTEGER, INTENT(IN)      ::  nk, bandi1, bandf1, bandi2, bandf2
      INTEGER, INTENT(IN)      :: nkpti
      INTEGER, INTENT(IN)      :: jsp

!     - arrays -

      REAL, INTENT(IN)         ::  eig_irr(:,:)
      COMPLEX, INTENT(OUT)     ::  dcprod(bandi2:bandf2, bandi1:bandf1, 3)

!     - local scalars -
      INTEGER                 ::  ikpt, ikpt1, iband1, iband2
      REAL                    ::  rdum
      LOGICAL                 ::  lwrite

      !                                       __
      ! Get momentum-matrix elements -i < uj | \/ | ui >
      !
      CALL momentum_matrix( &
         dcprod, nk, bandi1, bandf1, bandi2, bandf2, &
         atoms, hybrid, &
         cell, &
         hybdat, kpts, lapw, &
         dimension, jsp)

      !                                                __
      !  Calculate expansion coefficients -i < uj | \/ | ui > / ( ei - ej ) for periodic function ui
      !
      DO iband1 = bandi1, bandf1
         DO iband2 = bandi2, bandf2
            rdum = eig_irr(iband2, nk) - eig_irr(iband1, nk)
            IF (abs(rdum) > 1e-6) THEN !10.0**-6
               dcprod(iband2, iband1, :) = dcprod(iband2, iband1, :)/rdum
            ELSE
               dcprod(iband2, iband1, :) = 0.0
            END IF
         END DO
      END DO

      dcprod = dcprod/sqrt(cell%omtil)

   END SUBROUTINE dwavefproducts

!     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!     Calculates the momentum matrix elements
!
!                                    s             s
!     momentum(n',n,q,xyz) = -i < phi  | d    | phi    >   ,   s  = ispin
!                                    nq   xyz      n'q     ,   n  = occ.    ( bandi1 <= n  <= bandf1 )  ,
!                                                              n' = unocc.  ( bandi2 <= n' <= bandf2 )
!
!
   SUBROUTINE momentum_matrix( &
      momentum, nk, bandi1, bandf1, bandi2, bandf2, &
      atoms, hybrid, &
      cell, &
      hybdat, kpts, lapw, &
      dimension, jsp)

      USE m_olap
      USE m_wrapper
      USE m_util, only: derivative, intgrf_init, intgrf
      USE m_dr2fdr
      USE m_constants
      USE m_types
      USE m_io_hybrid
      IMPLICIT NONE

      TYPE(t_hybdat), INTENT(IN)   :: hybdat
      TYPE(t_dimension), INTENT(IN)   :: dimension
      TYPE(t_hybrid), INTENT(IN)   :: hybrid
      TYPE(t_cell), INTENT(IN)   :: cell
      TYPE(t_kpts), INTENT(IN)   :: kpts
      TYPE(t_atoms), INTENT(IN)   :: atoms
      TYPE(t_lapw), INTENT(IN)   :: lapw

!     - scalars -
      INTEGER, INTENT(IN)      ::  bandi1, bandf1, bandi2, bandf2
      INTEGER, INTENT(IN)      :: nk
      INTEGER, INTENT(IN)      :: jsp

!     - arrays -
      TYPE(t_mat):: z
      COMPLEX, INTENT(OUT)     :: momentum(bandi2:bandf2, bandi1:bandf1, 3)

!     - local scalars -
      INTEGER                 ::  itype, ieq, ic, i, j, l, lm, n1, n2, ikpt, iband1, iband2, ll, mm
      INTEGER                 ::  lm_0, lm_1, lm0, lm1, lm2, lm3, n0, nn, n, l1, l2, m1, m2, ikpt1
      INTEGER                 ::  irecl_cmt, irecl_z, m
      COMPLEX                 ::  cdum
      COMPLEX                 ::  img = (0.0, 1.0)

!     - local arrays -
      INTEGER                 ::  gpt(3, lapw%nv(jsp))

      REAL                    ::  fcoeff((atoms%lmaxd + 1)**2, -1:1), gcoeff((atoms%lmaxd + 1)**2, -1:1)
      REAL                    ::  qmat1(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype), dbas1(atoms%jmtd)
      REAL                    ::  qmat2(hybrid%maxindx, hybrid%maxindx, 0:atoms%lmaxd, atoms%ntype), dbas2(atoms%jmtd)
      REAL                    ::  qg(lapw%nv(jsp), 3)

      COMPLEX                 ::  hlp(3, 3)
      COMPLEX                 ::  cvec1(hybrid%maxlmindx), cvec2(hybrid%maxlmindx), cvec3(hybrid%maxlmindx)
      COMPLEX                 ::  cmt1(hybrid%maxlmindx, bandi1:bandf1), cmt2(hybrid%maxlmindx, bandi2:bandf2)
      COMPLEX                 ::  carr1(3), carr2(3)
      COMPLEX                 ::  cmt(dimension%neigd, hybrid%maxlmindx, atoms%nat)
      REAL                    ::  olap_r(lapw%nv(jsp)*(lapw%nv(jsp) + 1)/2)
      COMPLEX                 ::  olap_c(lapw%nv(jsp)*(lapw%nv(jsp) + 1)/2)
      REAL                    ::  vec1_r(lapw%nv(jsp)), vec2_r(lapw%nv(jsp)), vec3_r(lapw%nv(jsp))
      COMPLEX                 ::  vec1_c(lapw%nv(jsp)), vec2_c(lapw%nv(jsp)), vec3_c(lapw%nv(jsp))

      ! read in cmt coefficients from direct access file cmt at kpoint nk

      call read_cmt(cmt, nk)

      ! read in z coefficients from direct access file z at kpoint nk

      call read_z(z, nk)

      !CALL intgrf_init(atoms%ntype,atoms%jmtd,atoms%jri,atoms%dx,atoms%rmsh,hybdat%gridf)
      gpt(:, 1:lapw%nv(jsp)) = lapw%gvec(:, 1:lapw%nv(jsp), jsp)

!     Define coefficients F and G
      lm = 0
      DO l = 0, atoms%lmaxd
         DO M = -l, l
            lm = lm + 1
            fcoeff(lm, -1) = -sqrt(1.0*(l + M + 1)*(l + M + 2)/(2*(2*l + 1)*(2*l + 3)))
            fcoeff(lm, 0) = sqrt(1.0*(l - M + 1)*(l + M + 1)/((2*l + 1)*(2*l + 3)))
            fcoeff(lm, 1) = -sqrt(1.0*(l - M + 1)*(l - M + 2)/(2*(2*l + 1)*(2*l + 3)))
            gcoeff(lm, -1) = sqrt(1.0*(l - M)*(l - M - 1)/(2*(2*l - 1)*(2*l + 1)))
            gcoeff(lm, 0) = sqrt(1.0*(l - M)*(l + M)/((2*l - 1)*(2*l + 1)))
            gcoeff(lm, 1) = sqrt(1.0*(l + M)*(l + M - 1)/(2*(2*l - 1)*(2*l + 1)))
         END DO
      END DO

!     Calculate olap int r**2*u*u' + w * int r*u*u, w = -l,l+1 ( -> qmat1/2 )
      qmat1 = 0
      qmat2 = 0
      ic = 0
      DO itype = 1, atoms%ntype
         DO l = 0, atoms%lmax(itype)
            DO n2 = 1, hybrid%num_radfun_per_l(l, itype)
               !ic = ic + 1
               CALL derivative(dbas1, hybdat%bas1(:, n2, l, itype), atoms, itype)
               dbas1 = dbas1 - hybdat%bas1(:, n2, l, itype)/atoms%rmsh(:, itype)

               CALL derivative(dbas2, hybdat%bas2(:, n2, l, itype), atoms, itype)
               dbas2 = dbas2 - hybdat%bas2(:, n2, l, itype)/atoms%rmsh(:, itype)

               IF (l /= 0) THEN
                  DO n1 = 1, hybrid%num_radfun_per_l(l - 1, itype)
                     ic = ic + 1
                     qmat1(n1, n2, l, itype) = intgrf(dbas1(:)*hybdat%bas1(:, n1, l - 1, itype) + &
                                                      dbas2(:)*hybdat%bas2(:, n1, l - 1, itype), atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf) &
                                               + intgrf((hybdat%bas1(:, n2, l, itype)*hybdat%bas1(:, n1, l - 1, itype) + hybdat%bas2(:, n2, l, itype)*hybdat%bas2(:, n1, l - 1, itype)) &
                                                        /atoms%rmsh(:, itype), atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)*(l + 1)

                  END DO
               END IF
               IF (l /= atoms%lmax(itype)) THEN
                  DO n1 = 1, hybrid%num_radfun_per_l(l + 1, itype)

                     qmat2(n1, n2, l, itype) = intgrf(dbas1(:)*hybdat%bas1(:, n1, l + 1, itype) + dbas2(:)*hybdat%bas2(:, n1, l + 1, itype), &
                                                      atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf) &
                                               - intgrf((hybdat%bas1(:, n2, l, itype)*hybdat%bas1(:, n1, l + 1, itype) + hybdat%bas2(:, n2, l, itype)*hybdat%bas2(:, n1, l + 1, itype)) &
                                                        /atoms%rmsh(:, itype), atoms%jri, atoms%jmtd, atoms%rmsh, atoms%dx, atoms%ntype, itype, hybdat%gridf)*l

                  END DO
               END IF

            END DO
         END DO
      END DO

      !                                                  __
      ! Calculate momentum matrix elements -i < uj | \/ | ui > wrt wave functions u (->momentum)
      !

      momentum = 0

      ! MT contribution

      ic = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            ic = ic + 1
            nn = sum((/((2*l + 1)*hybrid%num_radfun_per_l(l, itype), l=0, atoms%lmax(itype))/))
            DO iband1 = bandi1, bandf1
               cmt1(:nn, iband1) = cmt(iband1, :nn, ic)
            ENDDO
            DO iband2 = bandi2, bandf2
               cmt2(:nn, iband2) = cmt(iband2, :nn, ic)
            ENDDO
            DO iband1 = bandi1, bandf1

               cvec1 = 0; cvec2 = 0; cvec3 = 0
               ! build up left vector(s) ( -> cvec1/2/3 )
               lm_0 = 0              ! we start with s-functions (l=0)
               lm_1 = hybrid%num_radfun_per_l(0, itype) ! we start with p-functions (l=0+1)
               lm = 0
               DO l = 0, atoms%lmax(itype) - 1
                  n0 = hybrid%num_radfun_per_l(l, itype)
                  n1 = hybrid%num_radfun_per_l(l + 1, itype)
                  DO M = -l, l
                     lm = lm + 1
                     lm0 = lm_0 + (M + l)*n0
                     lm1 = lm_1 + (M + 1 + l + 1)*n1
                     lm2 = lm_1 + (M + l + 1)*n1
                     lm3 = lm_1 + (M - 1 + l + 1)*n1
                     cvec1(lm0 + 1:lm0 + n0) = fcoeff(lm, -1)*matmul(cmt1(lm1 + 1:lm1 + n1, iband1), qmat2(:n1, :n0, l, itype))
                     cvec2(lm0 + 1:lm0 + n0) = fcoeff(lm, 0)*matmul(cmt1(lm2 + 1:lm2 + n1, iband1), qmat2(:n1, :n0, l, itype))
                     cvec3(lm0 + 1:lm0 + n0) = fcoeff(lm, 1)*matmul(cmt1(lm3 + 1:lm3 + n1, iband1), qmat2(:n1, :n0, l, itype))
                  END DO
                  lm_0 = lm_0 + (2*l + 1)*n0
                  lm_1 = lm_1 + (2*l + 3)*n1
               END DO

               lm_0 = hybrid%num_radfun_per_l(0, itype) ! we start with p-functions (l=1)
               lm_1 = 0              ! we start with s-functions (l=1-1)
               lm = 1
               DO l = 1, atoms%lmax(itype)
                  n0 = hybrid%num_radfun_per_l(l, itype)
                  n1 = hybrid%num_radfun_per_l(l - 1, itype)
                  DO M = -l, l
                     lm = lm + 1
                     lm0 = lm_0 + (M + l)*n0
                     lm1 = lm_1 + (M + 1 + l - 1)*n1
                     lm2 = lm_1 + (M + l - 1)*n1
                     lm3 = lm_1 + (M - 1 + l - 1)*n1
                     IF (abs(M + 1) <= l - 1) cvec1(lm0 + 1:lm0 + n0) = cvec1(lm0 + 1:lm0 + n0) + gcoeff(lm, -1)*matmul(cmt1(lm1 + 1:lm1 + n1, iband1), qmat1(:n1, :n0, l, itype))
                     IF (abs(M) <= l - 1) cvec2(lm0 + 1:lm0 + n0) = cvec2(lm0 + 1:lm0 + n0) + gcoeff(lm, 0)*matmul(cmt1(lm2 + 1:lm2 + n1, iband1), qmat1(:n1, :n0, l, itype))
                     IF (abs(M - 1) <= l - 1) cvec3(lm0 + 1:lm0 + n0) = cvec3(lm0 + 1:lm0 + n0) + gcoeff(lm, 1)*matmul(cmt1(lm3 + 1:lm3 + n1, iband1), qmat1(:n1, :n0, l, itype))
                  END DO
                  lm_0 = lm_0 + (2*l + 1)*n0
                  lm_1 = lm_1 + (2*l - 1)*n1
               END DO
               ! multiply with right vector
               DO iband2 = bandi2, bandf2
                  momentum(iband2, iband1, 1) = momentum(iband2, iband1, 1) + dot_product(cvec1(:nn), cmt2(:nn, iband2))
                  momentum(iband2, iband1, 2) = momentum(iband2, iband1, 2) + dot_product(cvec2(:nn), cmt2(:nn, iband2))
                  momentum(iband2, iband1, 3) = momentum(iband2, iband1, 3) + dot_product(cvec3(:nn), cmt2(:nn, iband2))
               END DO ! iband2
            END DO ! iband1

         END DO ! ieq
      END DO ! itype

      ! Transform to cartesian coordinates
      hlp = 0
      hlp(1, 1) = 1.0/sqrt(2.0)
      hlp(1, 3) = -1.0/sqrt(2.0)
      hlp(2, 1) = -img/sqrt(2.0)
      hlp(2, 3) = -img/sqrt(2.0)
      hlp(3, 2) = 1.0
      DO iband1 = bandi1, bandf1
         DO iband2 = bandi2, bandf2
            momentum(iband2, iband1, :) = -img*matmul(momentum(iband2, iband1, :), transpose(hlp))
         END DO
      END DO

      ! plane-wave contribution

      CALL olap_pwp(z%l_real, olap_r, olap_c, gpt, lapw%nv(jsp), atoms, cell)

      DO nn = 1, lapw%nv(jsp)
         qg(nn, :) = matmul(kpts%bk(:, nk) + gpt(:, nn), cell%bmat)
      END DO

      if (z%l_real) THEN
      DO iband2 = bandi2, bandf2
         vec1_r = matvec(olap_r, z%data_r(:lapw%nv(jsp), iband2)*qg(:, 1))
         vec2_r = matvec(olap_r, z%data_r(:lapw%nv(jsp), iband2)*qg(:, 2))
         vec3_r = matvec(olap_r, z%data_r(:lapw%nv(jsp), iband2)*qg(:, 3))
         DO iband1 = bandi1, bandf1
            momentum(iband2, iband1, 1) = momentum(iband2, iband1, 1) + dot_product(z%data_r(:lapw%nv(jsp), iband1), vec1_r)
            momentum(iband2, iband1, 2) = momentum(iband2, iband1, 2) + dot_product(z%data_r(:lapw%nv(jsp), iband1), vec2_r)
            momentum(iband2, iband1, 3) = momentum(iband2, iband1, 3) + dot_product(z%data_r(:lapw%nv(jsp), iband1), vec3_r)
         END DO
      END DO
      else
      DO iband2 = bandi2, bandf2
         vec1_c = matvec(olap_c, z%data_c(:lapw%nv(jsp), iband2)*qg(:, 1))
         vec2_c = matvec(olap_c, z%data_c(:lapw%nv(jsp), iband2)*qg(:, 2))
         vec3_c = matvec(olap_c, z%data_c(:lapw%nv(jsp), iband2)*qg(:, 3))
         DO iband1 = bandi1, bandf1
            momentum(iband2, iband1, 1) = momentum(iband2, iband1, 1) + dot_product(z%data_c(:lapw%nv(jsp), iband1), vec1_c)
            momentum(iband2, iband1, 2) = momentum(iband2, iband1, 2) + dot_product(z%data_c(:lapw%nv(jsp), iband1), vec2_c)
            momentum(iband2, iband1, 3) = momentum(iband2, iband1, 3) + dot_product(z%data_c(:lapw%nv(jsp), iband1), vec3_c)
         END DO
      END DO
      end if

   END SUBROUTINE momentum_matrix

END MODULE m_kp_perturbation
