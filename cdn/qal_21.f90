MODULE m_qal21 
  !***********************************************************************
  ! Calculates qal21  needed to determine the off-diagonal parts of the 
  ! DOS
  !***********************************************************************
  !
CONTAINS
  SUBROUTINE qal_21(atoms, input,noccbd,we,ccof, noco,acof,bcof,mt21,lo21,uloulopn21, qal,qmat)

    USE m_rotdenmat
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: noccbd 
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (INout)  :: we(noccbd),qal(0:,:,:,:)!(0:3,atoms%ntype,DIMENSION%neigd,input%jspins)
    REAL,    INTENT (IN)  :: uloulopn21(atoms%nlod,atoms%nlod,atoms%ntype)
    COMPLEX, INTENT (IN)  :: ccof(-atoms%llod:atoms%llod,noccbd,atoms%nlod,atoms%nat,input%jspins)
    COMPLEX, INTENT (IN)  :: acof(:,0:,:,:)!(noccbd,0:DIMENSION%lmd,atoms%nat,input%jspins)
    COMPLEX, INTENT (IN)  :: bcof(:,0:,:,:)!(noccbd,0:DIMENSION%lmd,atoms%nat,input%jspins)
    REAL,    INTENT (OUT) :: qmat(0:,:,:,:)!(0:3,atoms%ntype,DIMENSION%neigd,4)
    TYPE (t_mt21), INTENT (IN) :: mt21(0:atoms%lmaxd,atoms%ntype)
    TYPE (t_lo21), INTENT (IN) :: lo21(0:atoms%lmaxd,atoms%ntype)

    !     ..
    !     .. Local Scalars ..
    INTEGER i,l,lo,lop ,natom,nn,ntyp
    INTEGER nt1,nt2,lm,n,ll1,ipol,icore,index,m
    REAL fac
    COMPLEX sumaa,sumbb,sumab,sumba
    COMPLEX, PARAMETER :: ci = (0.0,1.0)

    !     ..
    !     .. Local Arrays ..
    COMPLEX qlo(noccbd,atoms%nlod,atoms%nlod,atoms%ntype)
    COMPLEX qaclo(noccbd,atoms%nlod,atoms%ntype),qbclo(noccbd,atoms%nlod,atoms%ntype)
    COMPLEX qcloa(noccbd,atoms%nlod,atoms%ntype),qclob(noccbd,atoms%nlod,atoms%ntype)
    COMPLEX qal21(0:3,atoms%ntype,size(qmat,3))
    COMPLEX q_loc(2,2),q_hlp(2,2),chi(2,2)
    !     ..
    !     .. Intrinsic Functions ..
    INTRINSIC conjg
    !
    !--->    l-decomposed density for each occupied state
    !
    states : DO i = 1, noccbd
       nt1 = 1
       types_loop : DO n = 1 ,atoms%ntype
          nt2 = nt1 + atoms%neq(n) - 1
          ls : DO l = 0,3
             IF (i==1) THEN
             ENDIF
             sumaa = CMPLX(0.,0.) ; sumab = CMPLX(0.,0.) 
             sumbb = CMPLX(0.,0.) ; sumba = CMPLX(0.,0.)
             ll1 = l* (l+1)
             ms : DO m = -l,l
                lm = ll1 + m
                atoms_loop : DO natom = nt1,nt2
                   sumaa = sumaa + acof(i,lm,natom,1)* CONJG(acof(i,lm,natom,input%jspins))
                   sumbb = sumbb + bcof(i,lm,natom,1)* CONJG(bcof(i,lm,natom,input%jspins))
                   sumba = sumba + acof(i,lm,natom,1) * CONJG(bcof(i,lm,natom,input%jspins))
                   sumab = sumab + bcof(i,lm,natom,1) * CONJG(acof(i,lm,natom,input%jspins))
                ENDDO atoms_loop
             ENDDO ms
             qal21(l,n,i) = sumaa * mt21(l,n)%uun + sumbb * mt21(l,n)%ddn +&
                  sumba * mt21(l,n)%dun + sumab * mt21(l,n)%udn 
          ENDDO ls
          nt1 = nt1 + atoms%neq(n)
       ENDDO types_loop
    ENDDO states

    !---> initialize qlo

    qlo(:,:,:,:) = CMPLX(0.,0.)
    qaclo(:,:,:) = CMPLX(0.,0.)
    qcloa(:,:,:) = CMPLX(0.,0.)
    qclob(:,:,:) = CMPLX(0.,0.)
    qbclo(:,:,:) = CMPLX(0.,0.)

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
                DO i = 1, noccbd
                   qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +      &
                        bcof(i,lm,natom,1)*CONJG(ccof(m,i,lo,natom,input%jspins)) 
                   qbclo(i,lo,ntyp) = qbclo(i,lo,ntyp) +      &
                        ccof(m,i,lo,natom,1)*CONJG(bcof(i,lm,natom,input%jspins)) 
                   qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       &
                        acof(i,lm,natom,1)*CONJG(ccof(m,i,lo,natom,input%jspins)) 
                   qaclo(i,lo,ntyp) = qaclo(i,lo,ntyp) +       &
                        ccof(m,i,lo,natom,1)*CONJG(acof(i,lm,natom,input%jspins)) 
                ENDDO
             ENDDO
             DO lop = 1,atoms%nlo(ntyp)
                IF (atoms%llo(lop,ntyp).EQ.l) THEN
                   DO m = -l,l
                      DO i = 1, noccbd
                         qlo(i,lop,lo,ntyp) = qlo(i,lop,lo,ntyp) +  &
                              CONJG(ccof(m,i,lop,natom,input%jspins))*ccof(m,i,lo,natom,1) +&
                              CONJG(ccof(m,i,lo,natom,input%jspins))*ccof(m,i,lop,natom,1)
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
          DO i = 1, noccbd
             qal21(l,ntyp,i)= qal21(l,ntyp,i)  + &
                  qaclo(i,lo,ntyp)*lo21(lo,ntyp)%uulon +&
                  qcloa(i,lo,ntyp)*lo21(lo,ntyp)%uloun +&
                  qclob(i,lo,ntyp)*lo21(lo,ntyp)%ulodn +&
                  qbclo(i,lo,ntyp)*lo21(lo,ntyp)%dulon 
          END DO
          DO lop = 1,atoms%nlo(ntyp)
             IF (atoms%llo(lop,ntyp).EQ.l) THEN
                DO i = 1, noccbd
                   qal21(l,ntyp,i)= qal21(l,ntyp,i)  + &
                        qlo(i,lop,lo,ntyp)*uloulopn21(lop,lo,ntyp)
                ENDDO
             ENDIF
          ENDDO
       END DO
    END DO

    DO n = 1,atoms%ntype
       fac = 1./atoms%neq(n)
       qal21(:,n,:) = qal21(:,n,:) * fac
    ENDDO
    !
    ! rotate into global frame
    !
    TYPE_loop : DO n = 1,atoms%ntype 
       chi(1,1) =  EXP(-ci*noco%alph(n)/2)*COS(noco%beta(n)/2)
       chi(1,2) = -EXP(-ci*noco%alph(n)/2)*SIN(noco%beta(n)/2)
       chi(2,1) =  EXP( ci*noco%alph(n)/2)*SIN(noco%beta(n)/2)
       chi(2,2) =  EXP( ci*noco%alph(n)/2)*COS(noco%beta(n)/2)
       state : DO i = 1, noccbd
          lls : DO l = 0,3
             CALL rot_den_mat(noco%alph(n),noco%beta(n),&
                  qal(l,n,i,1),qal(l,n,i,2),qal21(l,n,i))
             IF (.FALSE.) THEN
                IF (n==1) WRITE(*,'(3i3,4f10.5)') l,n,i,qal21(l,n,i),&
                     qal(l,n,i,:)
                q_loc(1,1) = qal(l,n,i,1); q_loc(2,2) = qal(l,n,i,2)
                q_loc(1,2) = qal21(l,n,i); q_loc(2,1) = CONJG(q_loc(1,2))
                q_hlp = MATMUL( TRANSPOSE( CONJG(chi) ) ,q_loc)
                q_loc = MATMUL(q_hlp,chi)
                qmat(l,n,i,1) = REAL(q_loc(1,1))
                qmat(l,n,i,2) = REAL(q_loc(1,2))
                qmat(l,n,i,3) = AIMAG(q_loc(1,2))
                qmat(l,n,i,4) = REAL(q_loc(2,2))
                IF (n==1) WRITE(*,'(3i3,4f10.5)') l,n,i,qmat(l,n,i,:)
             ENDIF
          ENDDO lls
       ENDDO state
    ENDDO TYPE_loop

  END SUBROUTINE qal_21
END MODULE m_qal21
