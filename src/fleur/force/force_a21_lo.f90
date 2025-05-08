!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_forcea21lo
CONTAINS
   SUBROUTINE force_a21_lo(atoms, isp, itype, we, eig, ne, abc, &
                           aveccof, bveccof, cveccof, tlmplm, usdus, a21)
      !--------------------------------------------------------------------------
      ! This subroutine calculates the local orbital contribution to A21,
      ! which is the combination of the terms A17 and A20 according to the
      ! paper of R.Yu et al. (PRB vol.43 no.8 p.64111991).
      ! p.kurz nov. 1997
      !--------------------------------------------------------------------------

      USE m_types_setup
      USE m_types_usdus
      USE m_types_tlmplm
      USE m_types_cdnval
      USE m_types_abc

      IMPLICIT NONE

      TYPE(t_usdus), INTENT(IN) :: usdus
      TYPE(t_tlmplm), INTENT(IN) :: tlmplm
      TYPE(t_atoms), INTENT(IN) :: atoms
      TYPE(t_abc), INTENT(IN) :: abc

      INTEGER, INTENT(IN) :: itype, ne, isp

      REAL, INTENT(IN)    :: we(ne), eig(:) !(input%neig)
      REAL, INTENT(INOUT) :: a21(3, atoms%nat)
      COMPLEX, INTENT(IN)    :: aveccof(3, ne, 0:atoms%lmaxd*(atoms%lmaxd + 2), atoms%nat)
      COMPLEX, INTENT(IN)    :: bveccof(3, ne, 0:atoms%lmaxd*(atoms%lmaxd + 2), atoms%nat)
      COMPLEX, INTENT(IN)    :: cveccof(3, -atoms%llod:atoms%llod, ne, atoms%nlod, atoms%nat)

      COMPLEX tuulo, tdulo, tuloulo
      INTEGER lo, lop, l, lp, mp, lm, lmp, iatom, ie, i, lolop, loplo, m, lo1, s, n_lo, n_lop, iatom_l

      !--- ABBREVIATIONS --------------------------------------------------------
      ! ccof       : coefficient of the local orbital function (u_lo*Y_lm)
      ! cveccof    : is defined equivalently to aveccof, but with the LO-fct.
      ! tuulo,tdulo and tuloulo are the MT hamiltonian matrix elements of the
      ! local orbitals with the flapw basisfct. and with themselves.
      ! for information on nlo,llo,nlol,lo1l,uulon,dulon, and uloulopn see
      ! comments in setlomap.
      !--------------------------------------------------------------------------

      DO lo = 1, atoms%nlo(itype)
         lo1 = SUM(atoms%nlo(:itype - 1)) + lo
         l = atoms%llo(lo, itype)
         n_lo = 2 + count(atoms%llo(:lo, itype) == l)
         DO m = -l, l
            lm = l*(l + 1) + m
            DO lp = 0, atoms%lnonsph(itype)
               s = tlmplm%h_loc2_nonsph(itype)
               DO mp = -lp, lp
                  lmp = lp*(lp + 1) + mp
                  DO iatom_l = 1, atoms%neq(itype)
                     ! Check whether the t-matrixelement is 0
                     ! (indmat.EQ.-9999)
                     iatom = iatom_l + atoms%firstAtom(itype) - 1
                     tuulo = tlmplm%h_LO(lmp, m, lo1, isp, isp)
                     tdulo = tlmplm%h_LO(lmp + s, m, lo1, isp, isp)

                     DO ie = 1, ne
                        DO i = 1, 3
                           a21(i, iatom) = a21(i, iatom) + 2.0*AIMAG( &
                                           CONJG(abc%cof(ie, lmp, 1, iatom_l))*tuulo &
                                           *cveccof(i, m, ie, lo, iatom) &
                                           + CONJG(abc%cof(ie, lmp, 2, iatom_l))*tdulo &
                                           *cveccof(i, m, ie, lo, iatom) &
                                           + CONJG(abc%cof(ie, lm, n_lo, iatom_l)) &
                                           *conjg(tuulo)*aveccof(i, ie, lmp, iatom) &
                                           + CONJG(abc%cof(ie, lm, n_lo, iatom_l)) &
                                           *conjg(tdulo)*bveccof(i, ie, lmp, iatom) &
                                           )*we(ie)/atoms%neq(itype)
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            DO lop = 1, atoms%nlo(itype)
               lp = atoms%llo(lop, itype)
               n_lop = 2 + count(atoms%llo(:lop, itype) == lp)
               DO mp = -lp, lp
                  lmp = lp*(lp + 1) + mp
                  DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
                     iatom_l = iatom - atoms%firstAtom(itype) + 1
                     tuloulo = tlmplm%tuloulo_newer(m, mp, lo, lop, itype, isp, isp)
                     DO ie = 1, ne
                        DO i = 1, 3
                           a21(i, iatom) = a21(i, iatom) + 2.0*AIMAG( &
                                           +CONJG(abc%cof(ie, lm, n_lop, iatom_l)) &
                                           *tuloulo*cveccof(i, mp, ie, lop, iatom) &
                                           )*we(ie)/atoms%neq(itype)
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
               iatom_l = iatom - atoms%firstAtom(itype) + 1
               DO ie = 1, ne
                  DO i = 1, 3
                     a21(i, iatom) = a21(i, iatom) - 2.0*AIMAG( &
                                     (CONJG(abc%cof(ie, lm, 1, iatom_l))*cveccof(i, m, ie, lo, iatom) + &
                                    CONJG(abc%cof(ie, lm, n_lo, iatom_l))*aveccof(i, ie, lm, iatom))*usdus%uulon(lo, itype, isp) + &
                                     (CONJG(abc%cof(ie, lm, 2, iatom_l))*cveccof(i, m, ie, lo, iatom) + &
                                      CONJG(abc%cof(ie, lm, n_lo, iatom_l))*bveccof(i, ie, lm, iatom))* &
                                     usdus%dulon(lo, itype, isp))*eig(ie)*we(ie)/atoms%neq(itype)
                  END DO
               END DO
            END DO

            ! Consider only the lop with l_lop = l_lo
            DO lop = atoms%lo1l(l, itype), (atoms%lo1l(l, itype) + atoms%nlol(l, itype) - 1)
               DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
                  iatom_l = iatom - atoms%firstAtom(itype) + 1
                  DO ie = 1, ne
                     DO i = 1, 3
                        a21(i, iatom) = a21(i, iatom) - 2.0*AIMAG( &
                                        CONJG(abc%cof(ie, lm, n_lo, iatom_l))* &
                                        cveccof(i, m, ie, lop, iatom)* &
                                        usdus%uloulopn(lo, lop, itype, isp))* &
                                        eig(ie)*we(ie)/atoms%neq(itype)

                     END DO
                  END DO
               END DO
            END DO
         END DO! m
      END DO ! lo

   END SUBROUTINE force_a21_lo
END MODULE m_forcea21lo
