!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eparas
  !***********************************************************************
  ! Calculates qlo, enerlo and sqlo, which are needed to determine the 
  ! new energy parameters.
  ! Philipp Kurz 99/04
  !***********************************************************************
  ! also the 'normal' energy parameters are now included...
  !
  ! if (l_mcd) then mcd contains mcd spectrum: first index = polarization
  ! second = core level ; third = band index                  gb.2001
  ! corrected to work also for multiple LO's of same l at the same atom
  !                                                           gb.2005
  !*************** ABBREVIATIONS *****************************************
  ! qlo     : charge density of one local orbital at the current k-point
  ! sqlo    : qlo integrated over the Brillouin zone
  ! enerlo  : qlo*energy integrated over the Brillouin zone
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE eparas(jsp,atoms,noccbd, mpi,ikpt,ne,we,eig,skip_t,l_evp,eigVecCoeffs,&
       usdus, ncore,l_mcd,m_mcd, enerlo,sqlo,ener,sqal,qal,mcd)
    USE m_types
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(IN)        :: usdus
    TYPE(t_mpi),INTENT(IN)          :: mpi
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_eigVecCoeffs),INTENT(IN) :: eigVecCoeffs
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: noccbd,jsp     
    INTEGER, INTENT (IN) :: ne,ikpt  ,skip_t
    LOGICAL, INTENT (IN) :: l_mcd,l_evp
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN)  :: ncore(atoms%ntype)
    REAL,    INTENT (IN)  :: eig(:)!(dimension%neigd),
    REAL,    INTENT (IN)  :: we(noccbd) 
    COMPLEX, INTENT (IN)  :: m_mcd(:,:,:,:)!(dimension%nstd,(3+1)**2,3*ntypd ,2)
    REAL,    INTENT (INOUT) :: enerlo(atoms%nlod,atoms%ntype),sqlo(atoms%nlod,atoms%ntype)
    REAL,    INTENT (INOUT) :: ener(0:3,atoms%ntype),sqal(0:3,atoms%ntype)
    REAL,    INTENT (INOUT) :: qal(0:,:,:)!(0:3,atoms%ntype,dimension%neigd)
    REAL,    INTENT (INOUT) :: mcd(:,:,:)!(3*atoms%ntype,dimension%nstd,dimension%neigd)

    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lo,lop ,natom,nn,ntyp,m
    INTEGER nt1,nt2,lm,n,ll1,ipol,icore,index
    REAL fac
    COMPLEX suma,sumb,sumab,sumba
    !     ..
    !     .. Local Arrays ..
    REAL qlo(noccbd,atoms%nlod,atoms%nlod,atoms%ntype)
    REAL qaclo(noccbd,atoms%nlod,atoms%ntype),qbclo(noccbd,atoms%nlod,atoms%ntype)
    !     ..
    !
    !---> initialize ener, sqal, enerlo and sqlo on first call
    !

    IF ((ikpt.LE.mpi%isize).AND..NOT.l_evp) THEN
       IF (l_mcd) THEN
          mcd(:,:,:) = 0.0
       ENDIF
       ener(:,:) = 0.0
       sqal(:,:) = 0.0
       qal(:,:,:) = 0.0
       enerlo(:,:) = 0.0
       sqlo(:,:) = 0.0
    END IF
    !
    !--->    l-decomposed density for each occupied state
    !
    !         DO 140 i = (skip_t+1),ne    ! this I need for all states
    DO i = 1,ne              ! skip in next loop
       nt1 = 1
       DO n = 1,atoms%ntype
          fac = 1./atoms%neq(n)
          nt2 = nt1 + atoms%neq(n) - 1
          DO l = 0,3
             suma = CMPLX(0.,0.)
             sumb = CMPLX(0.,0.)
             ll1 = l* (l+1)
             DO m = -l,l
                lm = ll1 + m
                IF ( .NOT.l_mcd ) THEN
                   DO natom = nt1,nt2
                      suma = suma + eigVecCoeffs%acof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%acof(i,lm,natom,jsp))
                      sumb = sumb + eigVecCoeffs%bcof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%bcof(i,lm,natom,jsp))
                   ENDDO
                ELSE
                   suma = CMPLX(0.,0.) ; sumab = CMPLX(0.,0.) 
                   sumb = CMPLX(0.,0.) ; sumba = CMPLX(0.,0.)
                   DO natom = nt1,nt2
                      suma = suma + eigVecCoeffs%acof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%acof(i,lm,natom,jsp))
                      sumb = sumb + eigVecCoeffs%bcof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%bcof(i,lm,natom,jsp))
                      sumab= sumab + eigVecCoeffs%acof(i,lm,natom,jsp) *CONJG(eigVecCoeffs%bcof(i,lm,natom,jsp))
                      sumba= sumba + eigVecCoeffs%bcof(i,lm,natom,jsp) *CONJG(eigVecCoeffs%acof(i,lm,natom,jsp))
                   ENDDO
                   DO icore = 1, ncore(n)
                      DO ipol = 1, 3
                         index = 3*(n-1) + ipol
                         mcd(index,icore,i)=mcd(index,icore,i) + fac*(&
                              suma * CONJG(m_mcd(icore,lm+1,index,1))*m_mcd(icore,lm+1,index,1)  +&
                              sumb * CONJG(m_mcd(icore,lm+1,index,2))*m_mcd(icore,lm+1,index,2)  +&
                              sumab* CONJG(m_mcd(icore,lm+1,index,2))*m_mcd(icore,lm+1,index,1)  +&
                              sumba* CONJG(m_mcd(icore,lm+1,index,1))*m_mcd(icore,lm+1,index,2)  ) 
                      ENDDO
                   ENDDO
                ENDIF     ! end MCD
             ENDDO
             qal(l,n,i) = (suma+sumb*usdus%ddn(l,n,jsp))/atoms%neq(n)
          ENDDO
          nt1 = nt1 + atoms%neq(n)
       ENDDO
    ENDDO
    !
    !--->    perform Brillouin zone integration and summation over the
    !--->    bands in order to determine the energy parameters for each
    !--->    atom and angular momentum
    !
    DO l = 0,3
       DO n = 1,atoms%ntype
          DO i = (skip_t+1),noccbd
             ener(l,n) = ener(l,n) + qal(l,n,i)*we(i)*eig(i)
             sqal(l,n) = sqal(l,n) + qal(l,n,i)*we(i)
          ENDDO
       ENDDO
    ENDDO

    !---> initialize qlo

    qlo(:,:,:,:) = 0.0
    qaclo(:,:,:) = 0.0
    qbclo(:,:,:) = 0.0

    !---> density for each local orbital and occupied state

    natom = 0
    DO ntyp = 1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             ll1 = l* (l+1)
             DO m = -l,l
                lm = ll1 + m
                DO i = 1,ne
                   qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +REAL(&
                        eigVecCoeffs%bcof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,jsp))+&
                        eigVecCoeffs%ccof(m,i,lo,natom,jsp)*CONJG(eigVecCoeffs%bcof(i,lm,natom,jsp)) )
                   qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) + REAL(&
                        eigVecCoeffs%acof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,jsp))+&
                        eigVecCoeffs%ccof(m,i,lo,natom,jsp)*CONJG(eigVecCoeffs%acof(i,lm,natom,jsp)) )
                ENDDO
             ENDDO
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1,ne
                         qlo(i,lop,lo,ntyp) = qlo(i,lop,lo,ntyp) +  REAL(&
                              CONJG(eigVecCoeffs%ccof(m,i,lop,natom,jsp))*eigVecCoeffs%ccof(m,i,lo,natom,jsp))
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !---> perform brillouin zone integration and sum over bands

    DO ntyp = 1,atoms%ntype
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          ! llo > 3 used for unoccupied states only
          IF( l .GT. 3 ) CYCLE
          DO i = 1,ne
             qal(l,ntyp,i)= qal(l,ntyp,i)  + ( 1.0/atoms%neq(ntyp) )* (&
                  qaclo(i,lo,ntyp)*usdus%uulon(lo,ntyp,jsp)+qbclo(i,lo,ntyp)*usdus%dulon(lo,ntyp,jsp)     )
          END DO
          DO lop = 1,atoms%nlo(ntyp)
             IF (atoms%llo(lop,ntyp).EQ.l) THEN
                DO i = 1,ne
                   enerlo(lo,ntyp) = enerlo(lo,ntyp) +qlo(i,lop,lo,ntyp)*we(i)*eig(i)
                   sqlo(lo,ntyp) = sqlo(lo,ntyp) + qlo(i,lop,lo,ntyp)*we(i)
                   qal(l,ntyp,i)= qal(l,ntyp,i)  + ( 1.0/atoms%neq(ntyp) ) *&
                        qlo(i,lop,lo,ntyp)*usdus%uloulopn(lop,lo,ntyp,jsp)
                ENDDO
             ENDIF
          ENDDO
       END DO
    END DO

  END SUBROUTINE eparas
END MODULE m_eparas
