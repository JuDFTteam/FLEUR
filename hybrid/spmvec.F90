MODULE m_spmvec

      CONTAINS
         !Note this module contains a real/complex version of spmvec

         SUBROUTINE spmvec_invs(&
        &           atoms, hybrid,&
        &           hybdat, ikpt, kpts, cell,&
        &           coulomb_mt1, coulomb_mt2, coulomb_mt3,&
        &           coulomb_mtir, vecin,&
        &           vecout)

            USE m_wrapper
            USE m_constants
            USE m_types
            IMPLICIT NONE
            TYPE(t_hybdat), INTENT(IN)   :: hybdat
            TYPE(t_hybrid), INTENT(IN)   :: hybrid
            TYPE(t_cell), INTENT(IN)     :: cell
            TYPE(t_kpts), INTENT(IN)     :: kpts
            TYPE(t_atoms), INTENT(IN)    :: atoms

            ! - scalars -
            INTEGER, INTENT(IN) ::  ikpt

            ! - arrays -
            REAL, INTENT(IN) ::  coulomb_mt1(hybrid%maxindxm1 - 1, hybrid%maxindxm1 - 1,&
           &                                    0:hybrid%maxlcutm1, atoms%ntype)
            REAL, INTENT(IN) ::  coulomb_mt2(hybrid%maxindxm1 - 1, -hybrid%maxlcutm1:hybrid%maxlcutm1,&
           &                                    0:hybrid%maxlcutm1 + 1, atoms%nat)
            REAL, INTENT(IN) ::  coulomb_mt3(hybrid%maxindxm1 - 1, atoms%nat, atoms%nat)
#ifdef CPP_IRCOULOMBAPPROX
            REAL, INTENT(IN) ::  coulomb_mtir(:, :)
#else
            REAL, INTENT(IN) ::  coulomb_mtir(:)
#endif
            REAL, INTENT(IN) ::  vecin(:)!(hybrid%nbasm)
            REAL, INTENT(INOUT)::  vecout(:)!(hybrid%nbasm)

            ! - local scalars -
            INTEGER             ::  itype, ieq, iatom, ishift
            INTEGER             ::  itype1, ieq1, iatom1, ishift1
            INTEGER             ::  indx0, indx1, indx2, indx3, indx4
            INTEGER             ::  i, ibasm, igptm
            INTEGER             ::  l
            INTEGER             ::  n, m

            REAL                ::  gnorm

            COMPLEX, PARAMETER  ::  img = (0.0, 1.0)
            ! - local arrays -

            REAL                ::  vecinhlp(hybrid%nbasm(ikpt))
            REAL, ALLOCATABLE ::  coulhlp(:, :)

            vecinhlp = vecin

            CALL reorder(hybrid%nbasm(ikpt), hybrid%nbasp, atoms, hybrid%lcutm1, hybrid%maxlcutm1, hybrid%nindxm1, 1, vecinhlp)

            ibasm = 0
            iatom = 0
            DO itype = 1, atoms%ntype
               DO ieq = 1, atoms%neq(itype)
                  iatom = iatom + 1
                  DO l = 0, hybrid%lcutm1(itype)
                     DO m = -l, l
                        ibasm = ibasm + hybrid%nindxm1(l, itype) - 1
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
                  DO l = 0, hybrid%lcutm1(itype)
                     DO m = -l, l
                        indx1 = indx1 + 1
                        indx2 = indx2 + hybrid%nindxm1(l, itype) - 1
                        indx3 = indx3 + 1

                        vecout(indx1:indx2) = matmul(coulomb_mt1(:hybrid%nindxm1(l, itype) - 1, :hybrid%nindxm1(l, itype) - 1, l, itype),&
               &                vecinhlp(indx1:indx2))

                        vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:hybrid%nindxm1(l, itype) - 1, m, l, iatom)*vecinhlp(indx3)

                        indx1 = indx2
                     END DO

                  END DO
               END DO
            END DO

            IF (indx2 /= ibasm) STOP 'spmvec: error counting basis functions'

            IF (ikpt == 1) THEN
               iatom = 0
               indx0 = 0
               DO itype = 1, atoms%ntype
                  ishift = sum((/((2*l + 1)*(hybrid%nindxm1(l, itype) - 1), l=0, hybrid%lcutm1(itype))/))
                  DO ieq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     l = 0
                     m = 0

                     indx1 = indx0 + 1
                     indx2 = indx1 + hybrid%nindxm1(l, itype) - 2

                     iatom1 = 0
                     indx3 = ibasm
                     DO itype1 = 1, atoms%ntype
                        ishift1 = (hybrid%lcutm1(itype1) + 1)**2
                        DO ieq1 = 1, atoms%neq(itype1)
                           iatom1 = iatom1 + 1
                           indx4 = indx3 + (ieq1 - 1)*ishift1 + 1
                           IF (iatom == iatom1) CYCLE

                           vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt3(:hybrid%nindxm1(l, itype) - 1, iatom1, iatom)*vecinhlp(indx4)

                        END DO
                        indx3 = indx3 + atoms%neq(itype1)*ishift1
                     END DO

                     IF (indx3 /= hybrid%nbasp) STOP 'spmvec: error counting index indx3'

                     vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:hybrid%nindxm1(l, itype) - 1, 0, hybrid%maxlcutm1 + 1, iatom)*vecinhlp(indx3 + 1)

                     indx0 = indx0 + ishift
                  END DO

               END DO
            END IF

            ! compute vecout for the index-range from ibasm+1:nbasm

#ifdef CPP_IRCOULOMBAPPROX

            indx0 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybrid%lcutm1(itype)), itype=1, atoms%ntype)/)) + hybrid%ngptm
            indx1 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybrid%lcutm1(itype)), itype=1, atoms%ntype)/))

            CALL DGEMV('N', indx1, indx0, 1.0, coulomb_mtir, (hybrid%maxlcutm1 + 1)**2*atoms,&
           &          vecinhlp(ibasm + 1:), 1, 0.0, vecout(ibasm + 1:), 1)

            CALL DGEMV('T', indx1, hybrid, 1.0, coulomb_mtir(:indx1, indx1 + 1:),&
           &          indx1, vecinhlp(ibasm + 1:), 1, 0.0, vecout(ibasm + indx1 + 1:), 1)

!       vecout(ibasm+1:ibasm+indx1) = matmul( coulomb_mtir(:indx1,:indx0),vecinhlp(ibasm+1:ibasm+indx0) )
!       vecout(ibasm+indx1+1:ibasm+indx0) = matmul( conjg(transpose(coulomb_mtir(:indx1,indx1+1:indx0))),
!      &                                            vecinhlp(ibasm+1:ibasm+indx1) )

            indx0 = ibasm + indx1
            IF (indx0 /= hybrid%nbasp) STOP 'spmvec: error indx0'
            DO i = 1, hybrid%ngptm
               indx0 = indx0 + 1
               igptm = hybrid%pgptm(i)
               gnorm = sqrt(sum(matmul(kpts%bk(:) + hybrid%gptm(:, igptm), cell%bmat)**2))
               IF (gnorm == 0) CYCLE
               vecout(indx0) = vecout(indx0) + fpi*vecinhlp(indx0)/gnorm
            END DO

#else

            indx1 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybrid%lcutm1(itype)),&
           &                                      itype=1, atoms%ntype)/)) + hybrid%ngptm(ikpt)
            CALL dspmv('U', indx1, 1.0, coulomb_mtir, vecinhlp(ibasm + 1:), 1, 0.0, vecout(ibasm + 1:), 1)

#endif
            iatom = 0
            indx1 = ibasm; indx2 = 0; indx3 = 0
            DO itype = 1, atoms%ntype
               DO ieq = 1, atoms%neq(itype)
                  iatom = iatom + 1
                  DO l = 0, hybrid%lcutm1(itype)
                     n = hybrid%nindxm1(l, itype)
                     DO m = -l, l
                        indx1 = indx1 + 1
                        indx2 = indx2 + 1
                        indx3 = indx3 + n - 1

                        vecout(indx1) = vecout(indx1) + dotprod(coulomb_mt2(:n - 1, m, l, iatom), vecinhlp(indx2:indx3))
                        indx2 = indx3
                     END DO

                  END DO
               END DO
            END DO

            IF (ikpt == 1) THEN
               iatom = 0
               indx0 = 0
               DO itype = 1, atoms%ntype
                  ishift = sum((/((2*l + 1)*(hybrid%nindxm1(l, itype) - 1), l=0, hybrid%lcutm1(itype))/))
                  DO ieq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     indx1 = indx0 + 1
                     indx2 = indx1 + hybrid%nindxm1(0, itype) - 2
                     vecout(hybrid%nbasp + 1) = vecout(hybrid%nbasp + 1) + dotprod(coulomb_mt2(:hybrid%nindxm1(0, itype) - 1, 0, hybrid%maxlcutm1 + 1, iatom), vecinhlp(indx1:indx2))

                     indx0 = indx0 + ishift
                  END DO
               END DO

               iatom = 0
               indx0 = ibasm
               DO itype = 1, atoms%ntype
                  ishift = (hybrid%lcutm1(itype) + 1)**2
                  DO ieq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     indx1 = indx0 + 1

                     iatom1 = 0
                     indx2 = 0
                     DO itype1 = 1, atoms%ntype
                        ishift1 = sum((/((2*l + 1)*(hybrid%nindxm1(l, itype1) - 1), l=0, hybrid%lcutm1(itype1))/))
                        DO ieq1 = 1, atoms%neq(itype1)
                           iatom1 = iatom1 + 1
                           IF (iatom1 == iatom) CYCLE

                           indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                           indx4 = indx3 + hybrid%nindxm1(0, itype1) - 2

                           vecout(indx1) = vecout(indx1) + dotprod(coulomb_mt3(:hybrid%nindxm1(0, itype1) - 1, iatom, iatom1), vecinhlp(indx3:indx4))

                        END DO
                        indx2 = indx2 + atoms%neq(itype1)*ishift1
                     END DO
                     indx0 = indx0 + ishift
                  END DO
               END DO
               IF (indx0 /= hybrid%nbasp) STOP 'spmvec: error index counting (indx0)'
            END IF

            CALL reorder(hybrid%nbasm(ikpt), hybrid%nbasp, atoms, hybrid%lcutm1, hybrid%maxlcutm1, hybrid%nindxm1,&
           &             2,&
           &             vecout)

         END SUBROUTINE spmvec_invs

         SUBROUTINE spmvec_noinvs(&
          &           atoms, hybrid,&
          &           hybdat, ikpt, kpts, cell,&
          &           coulomb_mt1, coulomb_mt2, coulomb_mt3,&
          &           coulomb_mtir, vecin,&
          &           vecout)

            USE m_wrapper
            USE m_constants
            USE m_types
            IMPLICIT NONE
            TYPE(t_hybdat), INTENT(IN)   :: hybdat
            TYPE(t_hybrid), INTENT(IN)   :: hybrid
            TYPE(t_cell), INTENT(IN)     :: cell
            TYPE(t_kpts), INTENT(IN)     :: kpts
            TYPE(t_atoms), INTENT(IN)    :: atoms

            ! - scalars -
            INTEGER, INTENT(IN) ::  ikpt

            ! - arrays -
            REAL, INTENT(IN) ::  coulomb_mt1(hybrid%maxindxm1 - 1, hybrid%maxindxm1 - 1,&
           &                                    0:hybrid%maxlcutm1, atoms%ntype)
            COMPLEX, INTENT(IN) ::  coulomb_mt2(hybrid%maxindxm1 - 1, -hybrid%maxlcutm1:hybrid%maxlcutm1,&
           &                                    0:hybrid%maxlcutm1 + 1, atoms%nat)
            COMPLEX, INTENT(IN) ::  coulomb_mt3(hybrid%maxindxm1 - 1, atoms%nat, atoms%nat)
#ifdef CPP_IRCOULOMBAPPROX
            COMPLEX, INTENT(IN) ::  coulomb_mtir(:, :)
#else
            COMPLEX, INTENT(IN) ::  coulomb_mtir(:)
#endif
            COMPLEX, INTENT(IN) ::  vecin(:)!(hybrid%nbasm)
            COMPLEX, INTENT(OUT)::  vecout(:)!(hybrid%nbasm)

            ! - local scalars -
            INTEGER             ::  itype, ieq, iatom, ishift
            INTEGER             ::  itype1, ieq1, iatom1, ishift1
            INTEGER             ::  indx0, indx1, indx2, indx3, indx4
            INTEGER             ::  i, ibasm, igptm
            INTEGER             ::  l
            INTEGER             ::  n, m

            REAL                ::  gnorm

            COMPLEX, PARAMETER  ::  img = (0.0, 1.0)
            ! - local arrays -

            REAL                ::  vecr(hybrid%maxindxm1 - 1), veci(hybrid%maxindxm1 - 1)
            COMPLEX             ::  vecinhlp(hybrid%nbasm(ikpt))
            COMPLEX, ALLOCATABLE ::  coulhlp(:, :)

            vecinhlp = vecin

            CALL reorder(hybrid%nbasm(ikpt), hybrid%nbasp, atoms, hybrid%lcutm1, hybrid%maxlcutm1, hybrid%nindxm1, 1, vec_c=vecinhlp)

            ibasm = 0
            iatom = 0
            DO itype = 1, atoms%ntype
               DO ieq = 1, atoms%neq(itype)
                  iatom = iatom + 1
                  DO l = 0, hybrid%lcutm1(itype)
                     DO m = -l, l
                        ibasm = ibasm + hybrid%nindxm1(l, itype) - 1
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
                  DO l = 0, hybrid%lcutm1(itype)
                     DO m = -l, l
                        indx1 = indx1 + 1
                        indx2 = indx2 + hybrid%nindxm1(l, itype) - 1
                        indx3 = indx3 + 1

                        vecout(indx1:indx2) = matmul(coulomb_mt1(:hybrid%nindxm1(l, itype) - 1, :hybrid%nindxm1(l, itype) - 1, l, itype),&
               &                vecinhlp(indx1:indx2))

                        vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:hybrid%nindxm1(l, itype) - 1, m, l, iatom)*vecinhlp(indx3)

                        indx1 = indx2
                     END DO

                  END DO
               END DO
            END DO

            IF (indx2 /= ibasm) &
           &  STOP 'spmvec: error counting basis functions'

            IF (ikpt == 1) THEN
               iatom = 0
               indx0 = 0
               DO itype = 1, atoms%ntype
                  ishift = sum((/((2*l + 1)*(hybrid%nindxm1(l, itype) - 1), l=0, hybrid%lcutm1(itype))/))
                  DO ieq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     l = 0
                     m = 0

                     indx1 = indx0 + 1
                     indx2 = indx1 + hybrid%nindxm1(l, itype) - 2

                     iatom1 = 0
                     indx3 = ibasm
                     DO itype1 = 1, atoms%ntype
                        ishift1 = (hybrid%lcutm1(itype1) + 1)**2
                        DO ieq1 = 1, atoms%neq(itype1)
                           iatom1 = iatom1 + 1
                           indx4 = indx3 + (ieq1 - 1)*ishift1 + 1
                           IF (iatom == iatom1) CYCLE

                           vecout(indx1:indx2) = vecout(indx1:indx2)&
                &                              + coulomb_mt3(:hybrid%nindxm1(l, itype) - 1, iatom1, iatom)*vecinhlp(indx4)

                        END DO
                        indx3 = indx3 + atoms%neq(itype1)*ishift1
                     END DO

                     IF (indx3 /= hybrid%nbasp) STOP 'spmvec: error counting index indx3'

                     vecout(indx1:indx2) = vecout(indx1:indx2) + coulomb_mt2(:hybrid%nindxm1(l, itype) - 1, 0, hybrid%maxlcutm1 + 1, iatom)*vecinhlp(indx3 + 1)

                     indx0 = indx0 + ishift
                  END DO

               END DO
            END IF

            ! compute vecout for the index-range from ibasm+1:nbasm

#ifdef CPP_IRCOULOMBAPPROX

            indx0 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybrid%lcutm1(itype)), itype=1, atoms%ntype)/)) + hybrid%ngptm
            indx1 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybrid%lcutm1(itype)), itype=1, atoms%ntype)/))

            CALL ZGEMV('N', indx1, indx0, (1.0, 0.0), coulomb_mtir,&
           &          (hybrid%maxlcutm1 + 1)**2*atoms, vecinhlp(ibasm + 1:),&
           &          1, (0.0, 0.0), vecout(ibasm + 1:), 1)

            CALL ZGEMV('C', indx1, hybrid, (1.0, 0.0), coulomb_mtir(:indx1, indx1 + 1:)&
           &          , indx1, vecinhlp(ibasm + 1:), 1, (0.0, 0.0),&
           &          vecout(ibasm + indx1 + 1:), 1)

!       vecout(ibasm+1:ibasm+indx1) = matmul( coulomb_mtir(:indx1,:indx0),vecinhlp(ibasm+1:ibasm+indx0) )
!       vecout(ibasm+indx1+1:ibasm+indx0) = matmul( conjg(transpose(coulomb_mtir(:indx1,indx1+1:indx0))),
!      &                                            vecinhlp(ibasm+1:ibasm+indx1) )

            indx0 = ibasm + indx1
            IF (indx0 /= hybrid%nbasp) STOP 'spmvec: error indx0'
            DO i = 1, hybrid%ngptm
               indx0 = indx0 + 1
               igptm = hybrid%pgptm(i)
               gnorm = sqrt(sum(matmul(kpts%bk(:) + hybrid%gptm(:, igptm), cell%bmat)**2))
               IF (gnorm == 0) CYCLE
               vecout(indx0) = vecout(indx0) + fpi*vecinhlp(indx0)/gnorm
            END DO

#else

            indx1 = sum((/(((2*l + 1)*atoms%neq(itype), l=0, hybrid%lcutm1(itype)),&
           &                                      itype=1, atoms%ntype)/)) + hybrid%ngptm(ikpt)
            call zhpmv('U', indx1, (1.0, 0.0), coulomb_mtir, vecinhlp(ibasm + 1), 1, (0.0, 0.0), vecout(ibasm + 1), 1)

#endif
            iatom = 0
            indx1 = ibasm; indx2 = 0; indx3 = 0
            DO itype = 1, atoms%ntype
               DO ieq = 1, atoms%neq(itype)
                  iatom = iatom + 1
                  DO l = 0, hybrid%lcutm1(itype)
                     n = hybrid%nindxm1(l, itype)
                     DO m = -l, l
                        indx1 = indx1 + 1
                        indx2 = indx2 + 1
                        indx3 = indx3 + n - 1

                        vecout(indx1) = vecout(indx1) + dotprod(coulomb_mt2(:n - 1, m, l, iatom), vecinhlp(indx2:indx3))
                        indx2 = indx3
                     END DO

                  END DO
               END DO
            END DO

            IF (ikpt == 1) THEN
               iatom = 0
               indx0 = 0
               DO itype = 1, atoms%ntype
                  ishift = sum((/((2*l + 1)*(hybrid%nindxm1(l, itype) - 1), l=0, hybrid%lcutm1(itype))/))
                  DO ieq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     indx1 = indx0 + 1
                     indx2 = indx1 + hybrid%nindxm1(0, itype) - 2
                     vecout(hybrid%nbasp + 1) = vecout(hybrid%nbasp + 1) + dotprod(coulomb_mt2(:hybrid%nindxm1(0, itype) - 1, 0, hybrid%maxlcutm1 + 1, iatom), vecinhlp(indx1:indx2))

                     indx0 = indx0 + ishift
                  END DO
               END DO

               iatom = 0
               indx0 = ibasm
               DO itype = 1, atoms%ntype
                  ishift = (hybrid%lcutm1(itype) + 1)**2
                  DO ieq = 1, atoms%neq(itype)
                     iatom = iatom + 1
                     indx1 = indx0 + 1

                     iatom1 = 0
                     indx2 = 0
                     DO itype1 = 1, atoms%ntype
                        ishift1 = sum((/((2*l + 1)*(hybrid%nindxm1(l, itype1) - 1), l=0, hybrid%lcutm1(itype1))/))
                        DO ieq1 = 1, atoms%neq(itype1)
                           iatom1 = iatom1 + 1
                           IF (iatom1 == iatom) CYCLE

                           indx3 = indx2 + (ieq1 - 1)*ishift1 + 1
                           indx4 = indx3 + hybrid%nindxm1(0, itype1) - 2

                           vecout(indx1) = vecout(indx1) + dotprod(coulomb_mt3(:hybrid%nindxm1(0, itype1) - 1, iatom, iatom1), vecinhlp(indx3:indx4))

                        END DO
                        indx2 = indx2 + atoms%neq(itype1)*ishift1
                     END DO
                     indx0 = indx0 + ishift
                  END DO
               END DO
               IF (indx0 /= hybrid%nbasp) STOP 'spmvec: error index counting (indx0)'
            END IF

            CALL reorder(hybrid%nbasm(ikpt), hybrid%nbasp, atoms, hybrid%lcutm1, hybrid%maxlcutm1, hybrid%nindxm1,&
           &             2,&
           &             vec_c=vecout)

         END SUBROUTINE spmvec_noinvs

         SUBROUTINE reorder(nbasm, nbasp, atoms, lcutm, maxlcutm, nindxm, imode, vec_r, vec_c)
            USE m_types
            IMPLICIT NONE
            TYPE(t_atoms), INTENT(IN)   :: atoms

            ! - scalars -
            INTEGER, INTENT(IN)   ::  maxlcutm
            INTEGER, INTENT(IN)   ::  nbasm, nbasp
            INTEGER, INTENT(IN)   ::  imode

            ! - arrays -
            INTEGER, INTENT(IN)   ::  lcutm(atoms%ntype)
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

            l_real = PRESENT(vec_r)

            IF (imode /= 1 .and. imode /= 2) STOP 'reorder: imode equals neither 1 nor 2'

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

         END SUBROUTINE reorder

      END MODULE
