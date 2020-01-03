MODULE m_spmvec

CONTAINS
   !Note this module contains a real/complex version of spmvec

   SUBROUTINE spmvec_invs(&
              atoms, mpdata, hybinp,&
              ikpt, &
              coulomb_mt1, coulomb_mt2, coulomb_mt3,&
              coulomb_mtir, vecin,&
              vecout)

      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_mpdata), intent(in)  :: mpdata
      TYPE(t_atoms), INTENT(IN)    :: atoms

      ! - scalars -
      INTEGER, INTENT(IN) ::  ikpt

      ! - arrays -
      REAL, INTENT(IN) ::  coulomb_mt1(maxval(mpdata%num_radbasfn) - 1, maxval(mpdata%num_radbasfn) - 1,&
                                          0:maxval(hybinp%lcutm1), atoms%ntype)
      REAL, INTENT(IN) ::  coulomb_mt2(maxval(mpdata%num_radbasfn) - 1, -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1),&
                                          0:maxval(hybinp%lcutm1) + 1, atoms%nat)
      REAL, INTENT(IN) ::  coulomb_mt3(:, :, :)
      REAL, INTENT(IN) ::  coulomb_mtir(:)
      REAL, INTENT(IN) ::  vecin(:)!(hybinp%nbasm)
      REAL, INTENT(INOUT)::  vecout(:)!(hybinp%nbasm)

      ! - local scalars -
      INTEGER             ::  itype, ieq, iatom, ishift
      INTEGER             ::  itype1, ieq1, iatom1, ishift1
      INTEGER             ::  indx0, indx1, indx2, indx3, indx4
      INTEGER             ::  ibasm
      INTEGER             ::  l, n, m
      ! - local arrays -

      REAL                ::  vecinhlp(hybinp%nbasm(ikpt))

      call timestart("spmvec_invs")
      vecinhlp = vecin

      CALL reorder(hybinp%nbasm(ikpt), atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), mpdata%num_radbasfn, 1, vecinhlp)

      ibasm = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybinp%lcutm1(itype)
               DO m = -l, l
                  ibasm = ibasm + mpdata%num_radbasfn(l, itype) - 1
               END DO
            END DO
         END DO
      END DO

      ! compute vecout for the indices from 0:ibasm
      iatom = 0
      indx1 = 0; indx2 = 0; indx3 = ibasm
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybinp%lcutm1(itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + mpdata%num_radbasfn(l, itype) - 1
                  indx3 = indx3 + 1

                  vecout(indx1:indx2) = matmul(coulomb_mt1(:mpdata%num_radbasfn(l, itype) - 1, :mpdata%num_radbasfn(l, itype) - 1, l, itype),&
                          vecinhlp(indx1:indx2))

                 vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, m, l, iatom)*vecinhlp(indx3)

                  indx1 = indx2
               END DO

            END DO
         END DO
      END DO

      IF (indx2 /= ibasm) call judft_error('spmvec: error counting basis functions')

      IF (ikpt == 1) THEN
         iatom = 0
         indx0 = 0
         DO itype = 1, atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, hybinp%lcutm1(itype))])
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               l = 0
               m = 0

               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2

               iatom1 = 0
               indx3 = ibasm
               DO itype1 = 1, atoms%ntype
                  ishift1 = (hybinp%lcutm1(itype1) + 1)**2
                  DO ieq1 = 1, atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     indx4 = indx3 + (ieq1 - 1)*ishift1 + 1
                     IF (iatom == iatom1) CYCLE

               vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt3(:mpdata%num_radbasfn(l, itype) - 1, iatom1, iatom)*vecinhlp(indx4)

                  END DO
                  indx3 = indx3 + atoms%neq(itype1)*ishift1
               END DO

               IF (indx3 /= hybinp%nbasp) call judft_error('spmvec: error counting index indx3')

                     vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom)*vecinhlp(indx3 + 1)

               indx0 = indx0 + ishift
            END DO

         END DO
      END IF

      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybinp%lcutm1(itype)),&
                                            itype=1, atoms%ntype)/)) + mpdata%n_g(ikpt)
      CALL dspmv('U', indx1, 1.0, coulomb_mtir, vecinhlp(ibasm + 1:), 1, 0.0, vecout(ibasm + 1:), 1)

      iatom = 0
      indx1 = ibasm; indx2 = 0; indx3 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybinp%lcutm1(itype)
               n = mpdata%num_radbasfn(l, itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + 1
                  indx3 = indx3 + n - 1

                  vecout(indx1) = vecout(indx1) + dot_product(coulomb_mt2(:n - 1, m, l, iatom), vecinhlp(indx2:indx3))
                  indx2 = indx3
               END DO

            END DO
         END DO
      END DO

      IF (ikpt == 1) THEN
         iatom = 0
         indx0 = 0
         DO itype = 1, atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, hybinp%lcutm1(itype))])
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(0, itype) - 2
                     vecout(hybinp%nbasp + 1) = vecout(hybinp%nbasp + 1) + dot_product(coulomb_mt2(:mpdata%num_radbasfn(0, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom), vecinhlp(indx1:indx2))

               indx0 = indx0 + ishift
            END DO
         END DO

         iatom = 0
         indx0 = ibasm
         DO itype = 1, atoms%ntype
            ishift = (hybinp%lcutm1(itype) + 1)**2
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1

               iatom1 = 0
               indx2 = 0
               DO itype1 = 1, atoms%ntype
                  ishift1 = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype1) - 1), l=0, hybinp%lcutm1(itype1))])
                  DO ieq1 = 1, atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     IF (iatom1 == iatom) CYCLE

                     indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                     indx4 = indx3 + mpdata%num_radbasfn(0, itype1) - 2

      vecout(indx1) = vecout(indx1) + dot_product(coulomb_mt3(:mpdata%num_radbasfn(0, itype1) - 1, iatom, iatom1), vecinhlp(indx3:indx4))

                  END DO
                  indx2 = indx2 + atoms%neq(itype1)*ishift1
               END DO
               indx0 = indx0 + ishift
            END DO
         END DO
         IF (indx0 /= hybinp%nbasp) call judft_error('spmvec: error index counting (indx0)')
      END IF

      CALL reorder(hybinp%nbasm(ikpt), atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), &
                   mpdata%num_radbasfn,2, vecout)
     call timestop("spmvec_invs")
   END SUBROUTINE spmvec_invs

   SUBROUTINE spmvec_noinvs(&
                atoms, mpdata, hybinp,&
                ikpt, &
                coulomb_mt1, coulomb_mt2, coulomb_mt3,&
                coulomb_mtir, vecin,&
                vecout)

      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(IN)  :: mpdata
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_atoms), INTENT(IN)    :: atoms

      ! - scalars -
      INTEGER, INTENT(IN) ::  ikpt

      ! - arrays -
      REAL, INTENT(IN) ::  coulomb_mt1(maxval(mpdata%num_radbasfn) - 1, maxval(mpdata%num_radbasfn) - 1,&
                                          0:maxval(hybinp%lcutm1), atoms%ntype)
      COMPLEX, INTENT(IN) ::  coulomb_mt2(maxval(mpdata%num_radbasfn) - 1, -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1),&
                                          0:maxval(hybinp%lcutm1) + 1, atoms%nat)
      COMPLEX, INTENT(IN) ::  coulomb_mt3(:, :, :)
      COMPLEX, INTENT(IN) ::  coulomb_mtir(:)
      COMPLEX, INTENT(IN) ::  vecin(:)!(hybinp%nbasm)
      COMPLEX, INTENT(OUT)::  vecout(:)!(hybinp%nbasm)

      ! - local scalars -
      INTEGER             ::  itype, ieq, iatom, ishift
      INTEGER             ::  itype1, ieq1, iatom1, ishift1
      INTEGER             ::  indx0, indx1, indx2, indx3, indx4
      INTEGER             ::  ibasm
      INTEGER             ::  l
      INTEGER             ::  n, m

      ! - local arrays -

      COMPLEX             ::  vecinhlp(hybinp%nbasm(ikpt))

      call timestart("spmvec_noinvs")
      vecinhlp = vecin

      CALL reorder(hybinp%nbasm(ikpt), atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), mpdata%num_radbasfn, 1, vec_c=vecinhlp)

      ibasm = 0
      iatom = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybinp%lcutm1(itype)
               DO m = -l, l
                  ibasm = ibasm + mpdata%num_radbasfn(l, itype) - 1
               END DO
            END DO
         END DO
      END DO

      ! compute vecout for the indices from 0:ibasm
      iatom = 0
      indx1 = 0; indx2 = 0; indx3 = ibasm
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybinp%lcutm1(itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + mpdata%num_radbasfn(l, itype) - 1
                  indx3 = indx3 + 1

                  vecout(indx1:indx2) = matmul(coulomb_mt1(:mpdata%num_radbasfn(l, itype) - 1, :mpdata%num_radbasfn(l, itype) - 1, l, itype),&
                          vecinhlp(indx1:indx2))

                 vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, m, l, iatom)*vecinhlp(indx3)

                  indx1 = indx2
               END DO

            END DO
         END DO
      END DO

      IF (indx2 /= ibasm) call judft_error('spmvec: error counting basis functions')

      IF (ikpt == 1) THEN
         iatom = 0
         indx0 = 0
         DO itype = 1, atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, hybinp%lcutm1(itype))])
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               l = 0
               m = 0

               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(l, itype) - 2

               iatom1 = 0
               indx3 = ibasm
               DO itype1 = 1, atoms%ntype
                  ishift1 = (hybinp%lcutm1(itype1) + 1)**2
                  DO ieq1 = 1, atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     indx4 = indx3 + (ieq1 - 1)*ishift1 + 1
                     IF (iatom == iatom1) CYCLE

                     vecout(indx1:indx2) = vecout(indx1:indx2)&
                                         + coulomb_mt3(:mpdata%num_radbasfn(l, itype) - 1, iatom1, iatom)*vecinhlp(indx4)

                  END DO
                  indx3 = indx3 + atoms%neq(itype1)*ishift1
               END DO

               IF (indx3 /= hybinp%nbasp) call judft_error('spmvec: error counting index indx3')

                     vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom)*vecinhlp(indx3 + 1)

               indx0 = indx0 + ishift
            END DO

         END DO
      END IF

      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybinp%lcutm1(itype)),&
                                            itype=1, atoms%ntype)/)) + mpdata%n_g(ikpt)
      call zhpmv('U', indx1, (1.0, 0.0), coulomb_mtir, vecinhlp(ibasm + 1), 1, (0.0, 0.0), vecout(ibasm + 1), 1)

      iatom = 0
      indx1 = ibasm; indx2 = 0; indx3 = 0
      DO itype = 1, atoms%ntype
         DO ieq = 1, atoms%neq(itype)
            iatom = iatom + 1
            DO l = 0, hybinp%lcutm1(itype)
               n = mpdata%num_radbasfn(l, itype)
               DO m = -l, l
                  indx1 = indx1 + 1
                  indx2 = indx2 + 1
                  indx3 = indx3 + n - 1

                  vecout(indx1) = vecout(indx1) + dot_product(coulomb_mt2(:n - 1, m, l, iatom), vecinhlp(indx2:indx3))
                  indx2 = indx3
               END DO

            END DO
         END DO
      END DO

      IF (ikpt == 1) THEN
         iatom = 0
         indx0 = 0
         DO itype = 1, atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, hybinp%lcutm1(itype))])
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(0, itype) - 2
                     vecout(hybinp%nbasp + 1) = vecout(hybinp%nbasp + 1) + dot_product(coulomb_mt2(:mpdata%num_radbasfn(0, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom), vecinhlp(indx1:indx2))

               indx0 = indx0 + ishift
            END DO
         END DO

         iatom = 0
         indx0 = ibasm
         DO itype = 1, atoms%ntype
            ishift = (hybinp%lcutm1(itype) + 1)**2
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1

               iatom1 = 0
               indx2 = 0
               DO itype1 = 1, atoms%ntype
                  ishift1 = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype1) - 1), l=0, hybinp%lcutm1(itype1))])
                  DO ieq1 = 1, atoms%neq(itype1)
                     iatom1 = iatom1 + 1
                     IF (iatom1 == iatom) CYCLE

                     indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                     indx4 = indx3 + mpdata%num_radbasfn(0, itype1) - 2

      vecout(indx1) = vecout(indx1) + dot_product(coulomb_mt3(:mpdata%num_radbasfn(0, itype1) - 1, iatom, iatom1), vecinhlp(indx3:indx4))

                  END DO
                  indx2 = indx2 + atoms%neq(itype1)*ishift1
               END DO
               indx0 = indx0 + ishift
            END DO
         END DO
         IF (indx0 /= hybinp%nbasp) call judft_error('spmvec: error index counting (indx0)')
      END IF

      CALL reorder(hybinp%nbasm(ikpt), atoms, hybinp%lcutm1, maxval(hybinp%lcutm1), mpdata%num_radbasfn,&
                   2,&
                   vec_c=vecout)
     call timestop("spmvec_noinvs")
   END SUBROUTINE spmvec_noinvs

   SUBROUTINE reorder(nbasm, atoms, lcutm, maxlcutm, nindxm, imode, vec_r, vec_c)
      USE m_types
      USE m_juDFT
      IMPLICIT NONE
      TYPE(t_atoms), INTENT(IN)   :: atoms

      ! - scalars -
      INTEGER, INTENT(IN)   ::  maxlcutm
      INTEGER, INTENT(IN)   ::  nbasm
      INTEGER, INTENT(IN)   ::  imode

      ! - arrays -
      INTEGER, INTENT(IN)   ::  lcutm(:)
      INTEGER, INTENT(IN)   ::  nindxm(0:maxlcutm, atoms%ntype)
      REAL, INTENT(INOUT), OPTIONAL::  vec_r(nbasm)
      COMPLEX, INTENT(INOUT), OPTIONAL::  vec_c(nbasm)
      ! - local scalars -
      INTEGER               ::  itype, ieq
      INTEGER               ::  indx1, indx2
      INTEGER               ::  l
      INTEGER               ::  n, m
      LOGICAL               :: l_real
      ! - local arrays -
      REAL                  ::  vechlp_r(nbasm)
      COMPLEX               ::  vechlp_c(nbasm)

      call timestart("reorder")
      l_real = PRESENT(vec_r)

      IF (imode /= 1 .and. imode /= 2) call judft_error('reorder: imode equals neither 1 nor 2')

      if (l_real) THEN
         vechlp_r = vec_r
      else
         vechlp_c = vec_c
      end if

      IF (imode == 1) THEN
         indx1 = 0
         indx2 = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               DO l = 0, lcutm(itype)
                  DO m = -l, l
                     DO n = 1, nindxm(l, itype) - 1
                        indx1 = indx1 + 1
                        indx2 = indx2 + 1
                        if (l_real) THEN
                           vec_r(indx1) = vechlp_r(indx2)
                        else
                           vec_c(indx1) = vechlp_c(indx2)
                        endif
                     END DO
                     indx2 = indx2 + 1
                  END DO
               END DO
            END DO
         END DO

         indx2 = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               DO l = 0, lcutm(itype)
                  DO m = -l, l
                     indx1 = indx1 + 1
                     indx2 = indx2 + nindxm(l, itype)
                     if (l_real) THEN
                        vec_r(indx1) = vechlp_r(indx2)
                     else
                        vec_c(indx1) = vechlp_c(indx2)
                     endif

                  END DO
               END DO
            END DO
         END DO
      ELSE IF (imode == 2) THEN
         indx1 = 0
         indx2 = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               DO l = 0, lcutm(itype)
                  DO m = -l, l
                     DO n = 1, nindxm(l, itype) - 1
                        indx1 = indx1 + 1
                        indx2 = indx2 + 1
                        if (l_real) THEN
                           vec_r(indx2) = vechlp_r(indx1)
                        else
                           vec_c(indx2) = vechlp_c(indx1)
                        endif
                     END DO
                     indx2 = indx2 + 1
                  END DO
               END DO
            END DO
         END DO

         indx2 = 0
         DO itype = 1, atoms%ntype
            DO ieq = 1, atoms%neq(itype)
               DO l = 0, lcutm(itype)
                  DO m = -l, l
                     indx1 = indx1 + 1
                     indx2 = indx2 + nindxm(l, itype)
                     if (l_real) THEN
                        vec_r(indx2) = vechlp_r(indx1)
                     else
                        vec_c(indx2) = vechlp_c(indx1)
                     endif
                  END DO
               END DO
            END DO
         END DO
      END IF
      !IR must not be rearranged
      call timestop("reorder")
   END SUBROUTINE reorder
END MODULE
