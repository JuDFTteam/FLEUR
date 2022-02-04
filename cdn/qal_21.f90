MODULE m_qal21
  !***********************************************************************
  ! Calculates qal21  needed to determine the off-diagonal parts of the
  ! DOS
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE qal_21(atoms,banddos,input,noccbd,ev_list,nococonv,eigVecCoeffs,denCoeffsOffdiag,ikpt,dos)
    use m_types_nococonv
    USE m_types_setup
    USE m_types_dos
    USE m_types_cdnval, ONLY: t_eigVecCoeffs
    USE m_types_denCoeffsOffdiag
    USE m_rotdenmat
    use m_constants
    IMPLICIT NONE

    TYPE(t_input),             INTENT(IN)    :: input
    TYPE(t_nococonv),          INTENT(IN)    :: nococonv
    TYPE(t_atoms),             INTENT(IN)    :: atoms
    TYPE(t_banddos),           INTENT(IN)    :: banddos
    TYPE(t_eigVecCoeffs),      INTENT(IN)    :: eigVecCoeffs
    TYPE(t_denCoeffsOffdiag),  INTENT(IN)    :: denCoeffsOffdiag
    TYPE(t_dos),               INTENT(INOUT) :: dos

    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: noccbd,ikpt

    INTEGER, INTENT (IN) :: ev_list(noccbd)

    !     .. Local Scalars ..
    INTEGER i,l,lo,lop ,natom,nn,ntyp
    INTEGER nt1,nt2,lm,ll1,ipol,icore,index,m,n_dos
    REAL fac
    COMPLEX sumaa,sumbb,sumab,sumba

    !     .. Local Arrays ..
    COMPLEX qlo(noccbd,atoms%nlod,atoms%nlod,atoms%ntype)
    COMPLEX qaclo(noccbd,atoms%nlod,atoms%ntype),qbclo(noccbd,atoms%nlod,atoms%ntype)
    COMPLEX qcloa(noccbd,atoms%nlod,atoms%ntype),qclob(noccbd,atoms%nlod,atoms%ntype)
    COMPLEX qal21(0:3,size(banddos%dos_typelist),input%neig)
    COMPLEX q_loc(2,2),q_hlp(2,2),chi(2,2)
    REAL    qmat(0:3,atoms%ntype,input%neig,4)

    !     .. Intrinsic Functions ..
    INTRINSIC conjg
    qal21=0.0
      !---> initialize qlo

    qlo(:,:,:,:) = CMPLX(0.,0.)
    qaclo(:,:,:) = CMPLX(0.,0.)
    qcloa(:,:,:) = CMPLX(0.,0.)
    qclob(:,:,:) = CMPLX(0.,0.)
    qbclo(:,:,:) = CMPLX(0.,0.)
    !--->    l-decomposed density for each occupied state
    states : DO i = 1, noccbd
       DO n_dos=1,size(banddos%dos_typelist)
         ntyp=banddos%dos_typelist(n_dos)
         nt1 = sum(atoms%neq(:ntyp-1))+1
         nt2 = nt1 + atoms%neq(ntyp) - 1
          ls : DO l = 0,3
             IF (i==1) THEN
             ENDIF
             sumaa = CMPLX(0.,0.) ; sumab = CMPLX(0.,0.)
             sumbb = CMPLX(0.,0.) ; sumba = CMPLX(0.,0.)
             ll1 = l* (l+1)
             ms : DO m = -l,l
                lm = ll1 + m
                atoms_loop : DO natom = nt1,nt2
                   sumaa = sumaa + eigVecCoeffs%acof(i,lm,natom,1)* CONJG(eigVecCoeffs%acof(i,lm,natom,input%jspins))
                   sumbb = sumbb + eigVecCoeffs%bcof(i,lm,natom,1)* CONJG(eigVecCoeffs%bcof(i,lm,natom,input%jspins))
                   sumba = sumba + eigVecCoeffs%acof(i,lm,natom,1) * CONJG(eigVecCoeffs%bcof(i,lm,natom,input%jspins))
                   sumab = sumab + eigVecCoeffs%bcof(i,lm,natom,1) * CONJG(eigVecCoeffs%acof(i,lm,natom,input%jspins))
                ENDDO atoms_loop
             ENDDO ms
             qal21(l,n_dos,i) = sumaa * denCoeffsOffdiag%uu21n(l,ntyp) + sumbb * denCoeffsOffdiag%dd21n(l,ntyp) +&
                            sumba * denCoeffsOffdiag%du21n(l,ntyp) + sumab * denCoeffsOffdiag%ud21n(l,ntyp)
          ENDDO ls
       ENDDO
    ENDDO states



    !---> density for each local orbital and occupied state

    DO n_dos=1,SIZE(banddos%dos_typelist)
    ntyp = banddos%dos_typelist(n_dos)
       natom=sum(atoms%neq(:ntyp-1))
       DO nn = 1,atoms%neq(ntyp)
          natom = natom + 1
          DO lo = 1,atoms%nlo(ntyp)
             l = atoms%llo(lo,ntyp)
             ll1 = l* (l+1)
             DO m = -l,l
                lm = ll1 + m
                DO i = 1, noccbd
                   qbclo(i,lo,n_dos) = qbclo(i,lo,n_dos) +      &
                        eigVecCoeffs%bcof(i,lm,natom,1)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,input%jspins))
                   qbclo(i,lo,n_dos) = qbclo(i,lo,n_dos) +      &
                        eigVecCoeffs%ccof(m,i,lo,natom,1)*CONJG(eigVecCoeffs%bcof(i,lm,natom,input%jspins))
                   qaclo(i,lo,n_dos) = qaclo(i,lo,n_dos) +       &
                        eigVecCoeffs%acof(i,lm,natom,1)*CONJG(eigVecCoeffs%ccof(m,i,lo,natom,input%jspins))
                   qaclo(i,lo,n_dos) = qaclo(i,lo,n_dos) +       &
                        eigVecCoeffs%ccof(m,i,lo,natom,1)*CONJG(eigVecCoeffs%acof(i,lm,natom,input%jspins))
                ENDDO
             ENDDO
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1, noccbd
                         qlo(i,lop,lo,n_dos) = qlo(i,lop,lo,n_dos) +  &
                              CONJG(eigVecCoeffs%ccof(m,i,lop,natom,input%jspins))*eigVecCoeffs%ccof(m,i,lo,natom,1) +&
                              CONJG(eigVecCoeffs%ccof(m,i,lo,natom,input%jspins))*eigVecCoeffs%ccof(m,i,lop,natom,1)
                      ENDDO
                   ENDDO
                ENDIF
             ENDDO
          ENDDO
       ENDDO

       !---> perform brillouin zone integration and sum over bands

       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          DO i = 1, noccbd
             qal21(l,n_dos,i)= qal21(l,n_dos,i)  + &
                  qaclo(i,lo,n_dos)*denCoeffsOffdiag%uulo21n(lo,ntyp) +&
                  qcloa(i,lo,n_dos)*denCoeffsOffdiag%ulou21n(lo,ntyp) +&
                  qclob(i,lo,n_dos)*denCoeffsOffdiag%ulod21n(lo,ntyp) +&
                  qbclo(i,lo,n_dos)*denCoeffsOffdiag%dulo21n(lo,ntyp)
          END DO
          DO lop = 1,atoms%nlo(ntyp)
             IF (atoms%llo(lop,ntyp).EQ.l) THEN
                DO i = 1, noccbd
                   qal21(l,n_dos,i)= qal21(l,n_dos,i)  + &
                        qlo(i,lop,lo,n_dos)*denCoeffsOffdiag%uloulop21n(lop,lo,ntyp)
                ENDDO
             ENDIF
          ENDDO
       END DO
       qal21(:,n_dos,:) = qal21(:,n_dos,:)/atoms%neq(ntyp)
       !
       ! rotate into global frame
       !
       chi(1,1) =  EXP(-ImagUnit*nococonv%alph(ntyp)/2)*COS(nococonv%beta(ntyp)/2)
       chi(1,2) = -EXP(-ImagUnit*nococonv%alph(ntyp)/2)*SIN(nococonv%beta(ntyp)/2)
       chi(2,1) =  EXP( ImagUnit*nococonv%alph(ntyp)/2)*SIN(nococonv%beta(ntyp)/2)
       chi(2,2) =  EXP( ImagUnit*nococonv%alph(ntyp)/2)*COS(nococonv%beta(ntyp)/2)
       state : DO i = 1, noccbd
          lls : DO l = 0,3
             CALL rot_den_mat(nococonv%alph(ntyp),nococonv%beta(ntyp),&
                  dos%qal(l,n_dos,ev_list(i),ikpt,1),dos%qal(l,n_dos,ev_list(i),ikpt,2),qal21(l,n_dos,i))
          ENDDO lls
       ENDDO state
     ENDDO

  END SUBROUTINE qal_21
END MODULE m_qal21
