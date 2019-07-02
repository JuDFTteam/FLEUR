MODULE m_qmtsl
CONTAINS
  !***********************************************************************
  ! Calculates the mt-spheres contribution to the layer charge for states 
  !  {En} at the current k-point. 
  !                                      Yury Koroteev 2003
  !                     from eparas.F  by  Philipp Kurz 99/04
  !
  !***********************************************************************
  !
  SUBROUTINE q_mt_sl(jsp,atoms,sym,nobd,ikpt,ne,skip_t,noccbd,eigVecCoeffs,usdus,slab)
    USE m_types_setup
    USE m_types_usdus
    USE m_types_cdnval, ONLY: t_eigVecCoeffs, t_slab
    IMPLICIT NONE
    TYPE(t_usdus),INTENT(IN)        :: usdus
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_sym),INTENT(IN)          :: sym
    TYPE(t_eigVecCoeffs),INTENT(IN) :: eigVecCoeffs
    TYPE(t_slab), INTENT(INOUT)     :: slab
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: nobd,jsp      
    INTEGER, INTENT (IN) :: ne,ikpt ,skip_t,noccbd
    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lo ,natom,nn,ntyp,nt1,nt2,m
    INTEGER lm,n,ll1,ipol,icore,index,nl
    REAL fac,sabd,ss,qq
    COMPLEX suma,sumb,sumab,sumba
    !     ..
    !     .. Local Arrays ..
    REAL, ALLOCATABLE :: qlo(:,:,:),qmt(:,:),qmtlo(:,:)
    REAL, ALLOCATABLE :: qaclo(:,:,:),qbclo(:,:,:),qmttot(:,:)
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC conjg,cmplx


    ALLOCATE ( qlo(nobd,atoms%nlod,atoms%ntype),qmt(atoms%ntype,SIZE(slab%qmtsl,2)) )
    ALLOCATE ( qaclo(nobd,atoms%nlod,atoms%ntype),qbclo(nobd,atoms%nlod,atoms%ntype) )
    ALLOCATE ( qmttot(atoms%ntype,SIZE(slab%qmtsl,2)),qmtlo(atoms%ntype,SIZE(slab%qmtsl,2)) )
    !
    !--->    l-decomposed density for each valence state
    !
    !         DO 140 i = (skip_t+1),ne    ! this I need for all states
    DO i = 1,ne              ! skip in next loop
       nt1 = 1
       DO n = 1,atoms%ntype
          fac = 1./atoms%neq(n)
          nt2 = nt1 + atoms%neq(n) - 1
          sabd = 0.0
          DO l = 0,atoms%lmax(n)
             suma = CMPLX(0.,0.)
             sumb = CMPLX(0.,0.)
             ll1 = l* (l+1)
             DO m = -l,l
                lm = ll1 + m
                DO natom = nt1,nt2
                   suma = suma + eigVecCoeffs%acof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%acof(i,lm,natom,jsp))
                   sumb = sumb + eigVecCoeffs%bcof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%bcof(i,lm,natom,jsp))
                ENDDO
             enddo
             ss = suma + sumb*usdus%ddn(l,n,jsp)
             sabd = sabd + ss
          enddo
          qmt(n,i) = sabd*fac
          nt1 = nt1 + atoms%neq(n)
       enddo
    enddo
    !                  
    !---> initialize qlo
    !
    qlo=0.0
    qaclo=0.0
    qbclo=0.0
    !
    !---> density for each local orbital and valence state
    !
    natom = 0
    DO ntyp = 1,atoms%ntype
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             ll1 = l* (l+1)
             DO i = 1,ne
                DO m = -l,l
                   lm = ll1 + m
                   qlo(i,lo,ntyp) = qlo(i,lo,ntyp) +&
                        eigVecCoeffs%ccof(m,i,lo,natom,jsp)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,jsp))
                   qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +&
                        eigVecCoeffs%bcof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,jsp)) +&
                        eigVecCoeffs%ccof(m,i,lo,natom,jsp)*CONJG(eigVecCoeffs%bcof(i,lm,natom,jsp))
                   qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +&
                        eigVecCoeffs%acof(i,lm,natom,jsp)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,jsp)) +&
                        eigVecCoeffs%ccof(m,i,lo,natom,jsp)*CONJG(eigVecCoeffs%acof(i,lm,natom,jsp))
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    natom = 1
    DO ntyp = 1,atoms%ntype
       IF (sym%invsat(natom).EQ.1) THEN
          DO lo = 1,atoms%nlo(ntyp)
             DO i = 1,ne
                qlo(i,lo,ntyp) = 2*qlo(i,lo,ntyp)
             ENDDO
          ENDDO
       ENDIF
       natom = natom + atoms%neq(ntyp)
    ENDDO
    !
    !--->  l-decomposed density for each valence state
    !--->      ( a contribution from local orbitals)
    !--->                       and
    !--->  total  l-decomposed density for each valence state
    !
    DO i = 1,ne
       DO ntyp = 1,atoms%ntype
          fac = 1.0/atoms%neq(ntyp)
          qq = 0.0
          DO lo = 1,atoms%nlo(ntyp)
             qq = qq + qlo(i,lo,ntyp)*usdus%uloulopn(lo,lo,ntyp,jsp) +&
                  qaclo(i,lo,ntyp)*usdus%uulon(lo,ntyp,jsp)     +&
                  qbclo(i,lo,ntyp)*usdus%dulon(lo,ntyp,jsp)    
          ENDDO
          qmtlo(ntyp,i) = qq*fac
          qmttot(ntyp,i) = qmt(ntyp,i) + qmtlo(ntyp,i)
       ENDDO
    ENDDO
    !
    DO i = 1,ne
       DO nl = 1,slab%nsl
          qq = 0.0
          DO ntyp = 1,atoms%ntype
             qq = qq + qmttot(ntyp,i)*slab%nmtsl(ntyp,nl)
          ENDDO
          slab%qmtsl(nl,i,ikpt,jsp) = qq
       ENDDO
    ENDDO
    !        DO ntyp = 1,ntype
    !        write(*,*) qmttot(ntyp,1)
    !        write(*,*) (nmtsl(ntyp,nl),nl=1,nsl)
    !        ENDDO
    !
    DEALLOCATE ( qlo,qmt,qmtlo,qaclo,qbclo,qmttot )

  END SUBROUTINE q_mt_sl
END MODULE m_qmtsl
