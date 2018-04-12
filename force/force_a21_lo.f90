!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_forcea21lo
CONTAINS
  SUBROUTINE force_a21_lo(nobd,atoms,isp,itype,we,eig,ne,eigVecCoeffs,force,tlmplm,usdus,a21)
    !
    !***********************************************************************
    ! This subroutine calculates the local orbital contribution to A21,
    ! which is the combination of the terms A17 and A20 according to the
    ! paper of R.Yu et al. (PRB vol.43 no.8 p.64111991).
    ! p.kurz nov. 1997
    !***********************************************************************
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(IN)        :: usdus
    TYPE(t_tlmplm),INTENT(IN)       :: tlmplm
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_eigVecCoeffs),INTENT(IN) :: eigVecCoeffs
    TYPE(t_force),INTENT(IN)        :: force
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nobd     
    INTEGER, INTENT (IN) :: itype,ne,isp
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(nobd),eig(:)!(dimension%neigd)
    REAL, INTENT (INOUT) :: a21(3,atoms%nat)
    !     ..
    !     .. Local Scalars ..
    COMPLEX utulo,dtulo,cutulo,cdtulo,ulotulo
    INTEGER lo,lop,l,lp ,mp,lm,lmp,iatom,ie,i,lolop,loplo,in,m,lo1
    !     ..
    !     ..
    !*************** ABBREVIATIONS *****************************************
    ! ccof       : coefficient of the local orbital function (u_lo*Y_lm)
    ! cveccof    : is defined equivalently to aveccof, but with the LO-fct.
    ! tuulo,tdulo and tuloulo are the MT hamiltonian matrix elements of the
    ! local orbitals with the flapw basisfct. and with themselves.
    ! for information on nlo,llo,nlol,lo1l,uulon,dulon, and uloulopn see
    ! comments in setlomap.
    !***********************************************************************

    DO lo = 1,atoms%nlo(itype)
       lo1=SUM(atoms%nlo(:itype-1))+lo
       l = atoms%llo(lo,itype)
       DO m = -l,l
          lm = l* (l+1) + m
          DO lp = 0,atoms%lmax(itype)
             DO mp = -lp,lp
                lmp = lp* (lp+1) + mp
                DO iatom = sum(atoms%neq(:itype-1))+1,sum(atoms%neq(:itype))
                   !
                   !--->             check whether the t-matrixelement is 0
                   !--->             (indmat.EQ.-9999)
                   !
                   in = tlmplm%ind(lmp,lm,itype,1)
                   IF ((in.NE.-9999).OR.(lmp.EQ.lm)) THEN
                      utulo = tlmplm%tuulo(lmp,m,lo1,1)
                      dtulo = tlmplm%tdulo(lmp,m,lo1,1)
                      cutulo = conjg(tlmplm%tuulo(lmp,m,lo1,1))
                      cdtulo = conjg(tlmplm%tdulo(lmp,m,lo1,1))
                      DO ie = 1,ne
                         DO i = 1,3
                            a21(i,iatom)=a21(i,iatom)+2.0*aimag(&
                                 conjg(eigVecCoeffs%acof(ie,lmp,iatom,isp))*utulo&
                                 *force%cveccof(i,m,ie,lo,iatom)&
                                 + conjg(eigVecCoeffs%bcof(ie,lmp,iatom,isp))*dtulo&
                                 *force%cveccof(i,m,ie,lo,iatom)&
                                 + conjg(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))&
                                 *cutulo*force%aveccof(i,ie,lmp,iatom)&
                                 + conjg(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))&
                                 *cdtulo*force%bveccof(i,ie,lmp,iatom)&
                                 )*we(ie)/atoms%neq(itype)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO

             ENDDO
          ENDDO
          DO lop = 1,atoms%nlo(itype)
             lp = atoms%llo(lop,itype)
             DO mp = -lp,lp
                lmp = lp* (lp+1) + mp
                DO iatom = sum(atoms%neq(:itype-1))+1,sum(atoms%neq(:itype))
                   in = tlmplm%ind(lmp,lm,itype,1)
                   IF ((in.NE.-9999).OR.(lmp.EQ.lm)) THEN
                      lolop=DOT_PRODUCT(atoms%nlo(:itype-1),atoms%nlo(:itype-1)+1)/2
                      IF (lo.GE.lop) THEN
                         lolop = (lo-1)*lo/2 + lop + lolop
                         ulotulo = tlmplm%tuloulo(m,mp,lolop,1)
                      ELSE
                         loplo = (lop-1)*lop/2 + lo +lolop
                         ulotulo = conjg(tlmplm%tuloulo(mp,m,loplo,1))
                      ENDIF
                      DO ie = 1,ne
                         DO i = 1,3
                            a21(i,iatom)=a21(i,iatom)+2.0*aimag(&
                                 + conjg(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))&
                                 *ulotulo*force%cveccof(i,mp,ie,lop,iatom)&
                                 )*we(ie)/atoms%neq(itype)
                         ENDDO
                      ENDDO
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          DO iatom = sum(atoms%neq(:itype-1))+1,sum(atoms%neq(:itype))
             DO ie = 1,ne
                DO i = 1,3
                   a21(i,iatom)=a21(i,iatom)&
                        -2.0*aimag(&
                        (conjg(eigVecCoeffs%acof(ie,lm,iatom,isp))*force%cveccof(i,m,ie,lo,iatom)+&
                        conjg(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))*force%aveccof(i,ie,lm,iatom))*usdus%uulon(lo,itype,isp)+&
                        (conjg(eigVecCoeffs%bcof(ie,lm,iatom,isp))*force%cveccof(i,m,ie,lo,iatom)+&
                        conjg(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))*force%bveccof(i,ie,lm,iatom))*&
                        usdus%dulon(lo,itype,isp))*eig(ie)*we(ie)/atoms%neq(itype)
                ENDDO
             ENDDO
          ENDDO
          !--->       consider only the lop with l_lop = l_lo
          DO lop = atoms%lo1l(l,itype),(atoms%lo1l(l,itype)+atoms%nlol(l,itype)-1)
             DO iatom = sum(atoms%neq(:itype-1))+1,sum(atoms%neq(:itype))
                DO ie = 1,ne
                   DO i = 1,3
                      a21(i,iatom)=a21(i,iatom)-2.0*aimag(&
                           conjg(eigVecCoeffs%ccof(m,ie,lo,iatom,isp))*&
                           force%cveccof(i,m,ie,lop,iatom)*&
                           usdus%uloulopn(lo,lop,itype,isp))*&
                           eig(ie)*we(ie)/atoms%neq(itype)

                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          !--->    end of m loop
       ENDDO
       !---> end of lo loop
    ENDDO

  END SUBROUTINE force_a21_lo
END MODULE m_forcea21lo
