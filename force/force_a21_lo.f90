!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_forcea21lo
CONTAINS
   SUBROUTINE force_a21_lo(atoms,isp,itype,we,eig,ne,eigVecCoeffs,&
                           aveccof,bveccof,cveccof,tlmplm,usdus,a21)
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

      IMPLICIT NONE

      TYPE(t_usdus),        INTENT(IN) :: usdus
      TYPE(t_tlmplm),       INTENT(IN) :: tlmplm
      TYPE(t_atoms),        INTENT(IN) :: atoms
      TYPE(t_eigVecCoeffs), INTENT(IN) :: eigVecCoeffs

      INTEGER, INTENT (IN) :: itype, ne, isp

      REAL,    INTENT(IN)    :: we(ne),eig(:) !(input%neig)
      REAL,    INTENT(INOUT) :: a21(3,atoms%nat)
      COMPLEX, INTENT(IN)    :: aveccof(3,ne,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
      COMPLEX, INTENT(IN)    :: bveccof(3,ne,0:atoms%lmaxd*(atoms%lmaxd+2),atoms%nat)
      COMPLEX, INTENT(IN)    :: cveccof(3,-atoms%llod:atoms%llod,ne,atoms%nlod,atoms%nat)

      COMPLEX tuulo, tdulo,  tuloulo
      INTEGER lo, lop, l, lp , mp, lm, lmp, iatom, ie, i, lolop, loplo, m, lo1,s

      !--- ABBREVIATIONS --------------------------------------------------------
      ! ccof       : coefficient of the local orbital function (u_lo*Y_lm)
      ! cveccof    : is defined equivalently to aveccof, but with the LO-fct.
      ! tuulo,tdulo and tuloulo are the MT hamiltonian matrix elements of the
      ! local orbitals with the flapw basisfct. and with themselves.
      ! for information on nlo,llo,nlol,lo1l,uulon,dulon, and uloulopn see
      ! comments in setlomap.
      !--------------------------------------------------------------------------

      DO lo = 1,atoms%nlo(itype)
         lo1=SUM(atoms%nlo(:itype-1))+lo
         l = atoms%llo(lo,itype)
         DO m = -l,l
            lm = l* (l+1) + m
            DO lp = 0,atoms%lnonsph(itype)
               s=tlmplm%h_loc2_nonsph(itype)
               DO mp = -lp,lp
                  lmp = lp* (lp+1) + mp
                  DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
                     ! Check whether the t-matrixelement is 0
                     ! (indmat.EQ.-9999)

                     tuulo = tlmplm%h_LO(lmp,m,lo1,isp,isp)
                     tdulo = tlmplm%h_LO(lmp+s,m,lo1,isp,isp)
                 
                     DO ie = 1,ne
                        DO i = 1,3
                           a21(i,iatom)=a21(i,iatom)+2.0*AIMAG(&
                                 CONJG(eigVecCoeffs%abcof(ie,lmp,0,iatom,isp))*tuulo&
                                 *cveccof(i,m,ie,lo,iatom)&
                                 + CONJG(eigVecCoeffs%abcof(ie,lmp,1,iatom,isp))*tdulo&
                                 *cveccof(i,m,ie,lo,iatom)&
                                 + CONJG(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))&
                                 *conjg(tuulo)*aveccof(i,ie,lmp,iatom)&
                                 + CONJG(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))&
                                 *conjg(tdulo)*bveccof(i,ie,lmp,iatom)&
                                 )*we(ie)/atoms%neq(itype)
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            DO lop = 1, atoms%nlo(itype)
               lp = atoms%llo(lop,itype)
               DO mp = -lp, lp
                  lmp = lp* (lp+1) + mp
                  DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
                     tuloulo = tlmplm%tuloulo_newer(m,mp,lo,lop,itype,isp,isp)
                     DO ie = 1,ne
                        DO i = 1,3
                           a21(i,iatom)=a21(i,iatom)+2.0*AIMAG(&
                              + CONJG(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))&
                              *tuloulo*cveccof(i,mp,ie,lop,iatom)&
                              )*we(ie)/atoms%neq(itype)
                        END DO
                     END DO
                  END DO
               END DO
            END DO

            DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
               DO ie = 1,ne
                  DO i = 1,3
                     a21(i,iatom)=a21(i,iatom)-2.0*AIMAG(&
                        (CONJG(eigVecCoeffs%abcof(ie,lm,0,iatom,isp))*cveccof(i,m,ie,lo,iatom)+&
                        CONJG(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))*aveccof(i,ie,lm,iatom))*usdus%uulon(lo,itype,isp)+&
                        (CONJG(eigVecCoeffs%abcof(ie,lm,1,iatom,isp))*cveccof(i,m,ie,lo,iatom)+&
                        CONJG(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))*bveccof(i,ie,lm,iatom))*&
                        usdus%dulon(lo,itype,isp))*eig(ie)*we(ie)/atoms%neq(itype)
                  END DO
               END DO
            END DO

            ! Consider only the lop with l_lop = l_lo
            DO lop = atoms%lo1l(l,itype),(atoms%lo1l(l,itype)+atoms%nlol(l,itype)-1)
               DO iatom = atoms%firstAtom(itype), atoms%firstAtom(itype) + atoms%neq(itype) - 1
                  DO ie = 1,ne
                     DO i = 1,3
                        a21(i,iatom)=a21(i,iatom)-2.0*AIMAG(&
                           CONJG(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))*&
                           cveccof(i,m,ie,lop,iatom)*&
                           usdus%uloulopn(lo,lop,itype,isp))*&
                           eig(ie)*we(ie)/atoms%neq(itype)

                     END DO
                  END DO
               END DO
            END DO
         END DO! m
      END DO ! lo

   END SUBROUTINE force_a21_lo
END MODULE m_forcea21lo
