MODULE m_spmvec

CONTAINS
   !Note this module contains a real/complex version of spmvec

   SUBROUTINE spmvec_invs(&
              atoms, mpdata, hybinp, hybdat,&
              ikpt, &
              coulomb_mt1, coulomb_mt2, coulomb_mt3,&
              coulomb_mtir, vecin,&
              vecout)

      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_juDFT
      use m_reorder
      IMPLICIT NONE
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_hybdat), INTENT(IN)   :: hybdat
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
      REAL, INTENT(IN) ::  coulomb_mtir(:,:)
      REAL, INTENT(IN) ::  vecin(:)!(hybdat%nbasm)
      REAL, INTENT(INOUT)::  vecout(:)!(hybdat%nbasm)

      ! - local scalars -
      INTEGER             ::  itype, ieq, iatom, ishift
      INTEGER             ::  itype1, ieq1, iatom1, ishift1
      INTEGER             ::  indx0, indx1, indx2, indx3, indx4
      INTEGER             ::  ibasm
      INTEGER             ::  l, n, m
      ! - local arrays -

      REAL                ::  vecinhlp(hybdat%nbasm(ikpt))

      call timestart("spmvec_invs")
      vecinhlp = vecin

      call reorder_forw(hybdat%nbasm(ikpt), atoms, hybinp%lcutm1, mpdata%num_radbasfn, vecinhlp)

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

               IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

                     vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom)*vecinhlp(indx3 + 1)

               indx0 = indx0 + ishift
            END DO

         END DO
      END IF

      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum([(((2*l + 1)*atoms%neq(itype), l=0, hybinp%lcutm1(itype)),&
                                            itype=1, atoms%ntype)]) + mpdata%n_g(ikpt)
      call timestart("ibasm+1 -> dgemv")
      call dgemv("N", indx1, indx1, 1.0, coulomb_mtir, indx1, vecinhlp(ibasm + 1), 1, 1.0, vecout(ibasm + 1), 1 )
      call timestop("ibasm+1 -> dgemv")
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
                     vecout(hybdat%nbasp + 1) = vecout(hybdat%nbasp + 1) + dot_product(coulomb_mt2(:mpdata%num_radbasfn(0, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom), vecinhlp(indx1:indx2))

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
         IF (indx0 /= hybdat%nbasp) call judft_error('spmvec: error index counting (indx0)')
      END IF

      call reorder_back(hybdat%nbasm(ikpt), atoms, hybinp%lcutm1, mpdata%num_radbasfn, vecout)
     call timestop("spmvec_invs")
   END SUBROUTINE spmvec_invs

   SUBROUTINE spmvec_noinvs(&
                atoms, mpdata, hybinp, hybdat,&
                ikpt, &
                coulomb_mt1, coulomb_mt2, coulomb_mt3,&
                coulomb_mtir, vecin,&
                vecout)

      USE m_wrapper
      USE m_constants
      USE m_types
      USE m_juDFT
      use m_calc_l_m_from_lm
      use m_reorder
      IMPLICIT NONE
      TYPE(t_mpdata), INTENT(IN)  :: mpdata
      TYPE(t_hybinp), INTENT(IN)   :: hybinp
      TYPE(t_hybdat), INTENT(IN)   :: hybdat
      TYPE(t_atoms), INTENT(IN)    :: atoms

      ! - scalars -
      INTEGER, INTENT(IN) ::  ikpt

      ! - arrays -
      REAL, INTENT(IN) ::  coulomb_mt1(maxval(mpdata%num_radbasfn) - 1, maxval(mpdata%num_radbasfn) - 1,&
                                          0:maxval(hybinp%lcutm1), atoms%ntype)
      COMPLEX, INTENT(IN) ::  coulomb_mt2(maxval(mpdata%num_radbasfn) - 1, -maxval(hybinp%lcutm1):maxval(hybinp%lcutm1),&
                                          0:maxval(hybinp%lcutm1) + 1, atoms%nat)
      COMPLEX, INTENT(IN) ::  coulomb_mt3(:, :, :)
      COMPLEX, INTENT(IN) ::  coulomb_mtir(:,:)
      COMPLEX, INTENT(IN) ::  vecin(:)!(hybdat%nbasm)
      COMPLEX, INTENT(INOUT)::  vecout(:)!(hybdat%nbasm)

      ! - local scalars -
      INTEGER             ::  itype, ieq, iatom, ishift
      INTEGER             ::  itype1, ieq1, iatom1, ishift1
      INTEGER             ::  indx0, indx1, indx2, indx3, indx4
      INTEGER             ::  ibasm
      INTEGER             ::  l, lm
      INTEGER             ::  n, m, iatom2, l1, idx_start, idx_stop

      ! - local arrays -

      COMPLEX             ::  vecinhlp(hybdat%nbasm(ikpt))

      call timestart("spmvec_noinvs")
      vecout = CMPLX_NOT_INITALIZED
      vecinhlp = vecin

      call reorder_forw(hybdat%nbasm(ikpt), atoms, hybinp%lcutm1, mpdata%num_radbasfn, vecinhlp)

      ibasm = 0
      iatom = 0
      do iatom = 1, atoms%nat
         itype = atoms%itype(iatom)
         DO l = 0, hybinp%lcutm1(itype)
            ibasm = ibasm + (2*l+1) * (mpdata%num_radbasfn(l, itype) - 1)
         END DO
      END DO

      ! compute vecout for the indices from 0:ibasm

      call timestart("0->ibasm: matmul")
      iatom = 0
      indx3 = ibasm
      do iatom = 1, atoms%nat 
         itype = atoms%itype(iatom)
         do lm =1,(hybinp%lcutm1(itype)+1)**2
            call calc_l_m_from_lm(lm, l, m)
            idx_start = 0
            ! go through previous atoms
            do iatom2 =1,iatom-1
               do l1=0, hybinp%lcutm1(itype)
                  idx_start = idx_start + (2*l1+1)*(mpdata%num_radbasfn(l1, atoms%itype(iatom2)) - 1)
               enddo
            enddo
            ! current atom
            do l1=0, l-1
               idx_start = idx_start + (2*l1+1)*(mpdata%num_radbasfn(l1, itype) - 1)
            enddo
            !current l
            idx_start = 1 + idx_start + ((m + l) * (mpdata%num_radbasfn(l, itype)-1))


            idx_stop = idx_start + mpdata%num_radbasfn(l, itype) - 2
            indx3 = indx3 + 1

            vecout(idx_start:idx_stop) = matmul(coulomb_mt1(:mpdata%num_radbasfn(l, itype) - 1, :mpdata%num_radbasfn(l, itype) - 1, l, itype),&
                  vecinhlp(idx_start:idx_stop))

            vecout(idx_start:idx_stop) = vecout(idx_start:idx_stop) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, m, l, iatom)*vecinhlp(indx3)
         END DO
      END DO
      call timestop("0->ibasm: matmul")

      IF (idx_stop /= ibasm) call judft_error('spmvec: error counting basis functions')


      IF (ikpt == 1) THEN
         call timestart("gamma point 1")
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

               IF (indx3 /= hybdat%nbasp) call judft_error('spmvec: error counting index indx3')

                     vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:mpdata%num_radbasfn(l, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom)*vecinhlp(indx3 + 1)

               indx0 = indx0 + ishift
            END DO

         END DO
         call timestop("gamma point 1")
      END IF
      ! compute vecout for the index-range from ibasm+1:nbasm

      indx1 = sum([(((2*l + 1)*atoms%neq(itype), l=0, hybinp%lcutm1(itype)),&
                                            itype=1, atoms%ntype)]) + mpdata%n_g(ikpt)
      call timestart("ibasm+1->nbasm: zgemv")
      call zgemv("N", indx1, indx1, cmplx_1, coulomb_mtir, indx1, vecinhlp(ibasm + 1), 1, cmplx_0, vecout(ibasm + 1), 1 )
      call timestop("ibasm+1->nbasm: zgemv")

      call timestart("dot prod")
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
      call timestop("dot prod")


      IF (ikpt == 1) THEN
         call timestart("gamma point 2")
         iatom = 0
         indx0 = 0
         DO itype = 1, atoms%ntype
            ishift = sum([((2*l + 1)*(mpdata%num_radbasfn(l, itype) - 1), l=0, hybinp%lcutm1(itype))])
            DO ieq = 1, atoms%neq(itype)
               iatom = iatom + 1
               indx1 = indx0 + 1
               indx2 = indx1 + mpdata%num_radbasfn(0, itype) - 2
                     vecout(hybdat%nbasp + 1) = vecout(hybdat%nbasp + 1) + dot_product(coulomb_mt2(:mpdata%num_radbasfn(0, itype) - 1, 0, maxval(hybinp%lcutm1) + 1, iatom), vecinhlp(indx1:indx2))

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
         IF (indx0 /= hybdat%nbasp) call judft_error('spmvec: error index counting (indx0)')
         call timestop("gamma point 2")
      END IF

      call reorder_back(hybdat%nbasm(ikpt), atoms, hybinp%lcutm1, mpdata%num_radbasfn, vecout)
     call timestop("spmvec_noinvs")
   END SUBROUTINE spmvec_noinvs
END MODULE
