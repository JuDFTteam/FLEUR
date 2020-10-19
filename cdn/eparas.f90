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
  SUBROUTINE eparas(jsp,atoms,banddos,noccbd,ev_list,fmpi,ikpt,ne,we,eig,skip_t,l_evp,eigVecCoeffs,&
                    usdus,regCharges,dos,mcd)
    USE m_types
    use m_types_dos
    use m_types_mcd
    IMPLICIT NONE
    TYPE(t_usdus),         INTENT(IN)    :: usdus
    TYPE(t_mpi),           INTENT(IN)    :: fmpi
    TYPE(t_atoms),         INTENT(IN)    :: atoms
    TYPE(t_banddos),       INTENT(IN)    :: banddos
    TYPE(t_eigVecCoeffs),  INTENT(IN)    :: eigVecCoeffs
    TYPE(t_regionCharges), INTENT(INOUT) :: regCharges
    TYPE(t_dos),           INTENT(INOUT) :: dos
    TYPE(t_mcd), OPTIONAL, INTENT(INOUT) :: mcd
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: noccbd,jsp
    INTEGER, INTENT (IN) :: ne,ikpt  ,skip_t
    LOGICAL, INTENT (IN) :: l_evp
    INTEGER, INTENT (IN) :: ev_list(noccbd)
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN)  :: eig(:)!(input%neig),
    REAL,    INTENT (IN)  :: we(noccbd)

    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lo,lop ,natom,nn,ntyp,m,n_dos
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

    IF ((ikpt.LE.fmpi%isize).AND..NOT.l_evp) THEN
       regCharges%ener(:,:,jsp) = 0.0
       regCharges%sqal(:,:,jsp) = 0.0
       regCharges%enerlo(:,:,jsp) = 0.0
       regCharges%sqlo(:,:,jsp) = 0.0
       dos%qal(:,:,:,ikpt,jsp) = 0.0
    END IF

    !--->    l-decomposed density for each occupied state
    !
    !         DO 140 i = (skip_t+1),ne    ! this I need for all states
    DO i = 1,ne              ! skip in next loop
       nt1 = 1
       DO n_dos = 1,size(banddos%dos_typelist)
         n=banddos%dos_typelist(n_dos)
         fac = 1./atoms%neq(n)
         nt2 = nt1 + atoms%neq(n) - 1
         DO l = 0,3
           suma = CMPLX(0.,0.)
           sumb = CMPLX(0.,0.)
           ll1 = l* (l+1)
           DO m = -l,l
             lm = ll1 + m
             IF ( .NOT.banddos%l_mcd ) THEN
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
               DO icore = 1, mcd%ncore(n)
                 DO ipol = 1, 3
                   index = 3*(n_dos-1) + ipol
                   mcd%mcd(index,icore,ev_list(i),ikpt,jsp)=mcd%mcd(index,icore,ev_list(i),ikpt,jsp) + fac*(&
                   suma * CONJG(mcd%m_mcd(icore,lm+1,index,1))*mcd%m_mcd(icore,lm+1,index,1)  +&
                   sumb * CONJG(mcd%m_mcd(icore,lm+1,index,2))*mcd%m_mcd(icore,lm+1,index,2)  +&
                   sumab* CONJG(mcd%m_mcd(icore,lm+1,index,2))*mcd%m_mcd(icore,lm+1,index,1)  +&
                   sumba* CONJG(mcd%m_mcd(icore,lm+1,index,1))*mcd%m_mcd(icore,lm+1,index,2)  )
                 ENDDO
               ENDDO
             ENDIF     ! end MCD
           ENDDO
           dos%qal(l,n_dos,ev_list(i),ikpt,jsp) = (suma+sumb*usdus%ddn(l,n,jsp))/atoms%neq(n)
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
        if (banddos%map_atomtype(n)==0) then
          DO i = (skip_t+1),noccbd
            DO natom = sum(atoms%neq(:n-1))+1,sum(atoms%neq(:n))
              suma=suma+dot_product(eigVecCoeffs%acof(i,l* (l+1)-l:l* (l+1)+l,natom,jsp),eigVecCoeffs%acof(i,l* (l+1)-l:l* (l+1)+l,natom,jsp))
              sumb=sumb+dot_product(eigVecCoeffs%bcof(i,l* (l+1)-l:l* (l+1)+l,natom,jsp),eigVecCoeffs%bcof(i,l* (l+1)-l:l* (l+1)+l,natom,jsp))
            ENDDO
            regCharges%ener(l,n,jsp) = regCharges%ener(l,n,jsp) + (suma+sumb*usdus%ddn(l,n,jsp))/atoms%neq(n)*we(i)*eig(i)
            regCharges%sqal(l,n,jsp) = regCharges%sqal(l,n,jsp) + (suma+sumb*usdus%ddn(l,n,jsp))/atoms%neq(n)*we(i)
          ENDDO
        ELSE
          n_dos=banddos%map_atomtype(n)
          !data present in dos-type
          DO i = (skip_t+1),noccbd
            regCharges%ener(l,n,jsp) = regCharges%ener(l,n,jsp) + dos%qal(l,n_dos,ev_list(i),ikpt,jsp)*we(i)*eig(i)
            regCharges%sqal(l,n,jsp) = regCharges%sqal(l,n,jsp) + dos%qal(l,n_dos,ev_list(i),ikpt,jsp)*we(i)
          ENDDO
        ENDIF
      ENDDO
    ENDDO

    !---> initialize qlo

    qlo(:,:,:,:) = 0.0
    qaclo(:,:,:) = 0.0
    qbclo(:,:,:) = 0.0

    !---> density for each local orbital and occupied state

    DO ntyp=1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
         natom=sum(atoms%neq(:ntyp-1))
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
          n_dos=banddos%map_atomtype(ntyp)
          if (n_dos>0)THEN
            DO i = 1,ne
              dos%qal(l,n_dos,ev_list(i),ikpt,jsp)= dos%qal(l,n_dos,ev_list(i),ikpt,jsp)  + ( 1.0/atoms%neq(ntyp) )* (&
              qaclo(i,lo,ntyp)*usdus%uulon(lo,ntyp,jsp)+qbclo(i,lo,ntyp)*usdus%dulon(lo,ntyp,jsp)     )
            END DO
          ENDIF
          DO lop = 1,atoms%nlo(ntyp)
             IF (atoms%llo(lop,ntyp).EQ.l) THEN
                DO i = 1,ne
                   regCharges%enerlo(lo,ntyp,jsp) = regCharges%enerlo(lo,ntyp,jsp) +qlo(i,lop,lo,ntyp)*we(i)*eig(i)
                   regCharges%sqlo(lo,ntyp,jsp) = regCharges%sqlo(lo,ntyp,jsp) + qlo(i,lop,lo,ntyp)*we(i)
                   if (n_dos>0) dos%qal(l,n_dos,ev_list(i),ikpt,jsp)= dos%qal(l,n_dos,ev_list(i),ikpt,jsp)  + ( 1.0/atoms%neq(ntyp) ) *&
                        qlo(i,lop,lo,ntyp)*usdus%uloulopn(lop,lo,ntyp,jsp)
                ENDDO
             ENDIF
          ENDDO
       END DO
    END DO

  END SUBROUTINE eparas
END MODULE m_eparas
