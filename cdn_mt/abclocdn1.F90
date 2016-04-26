MODULE m_abclocdn1
  !*********************************************************************
  ! Calculates the basis coeffcients (bascof_lo) for the local
  ! orbitals. 
  !*********************************************************************
  !*************** ABBREVIATIONS ***************************************
  ! nkvec   : stores the number of G-vectors that have been found and
  !           accepted during the construction of the local orbitals.
  ! kvec    : k-vector used in hssphn to attach the local orbital 'lo'
  !           of atom 'na' to it.
  !*********************************************************************
CONTAINS
  SUBROUTINE abclocdn1(atoms,sym, con1,phase,ylm,ntyp,na,k,s,nv,&
       nbasf0,alo1,blo1,clo1,kvec, nkvec,enough,bascof_lo )
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: k,na,ntyp,nv
    REAL,    INTENT (IN) :: con1 ,s
    COMPLEX, INTENT (IN) :: phase
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nbasf0(atoms%nlod,atoms%natd) 
    REAL,    INTENT (IN) :: alo1(atoms%nlod,atoms%ntypd),blo1(atoms%nlod,atoms%ntypd)
    REAL,    INTENT (IN) :: clo1(atoms%nlod,atoms%ntypd)
    COMPLEX, INTENT (IN) :: ylm( (atoms%lmaxd+1)**2 )
    INTEGER, INTENT (IN) :: kvec(2*(2*atoms%llod+1) ,atoms%natd)
    LOGICAL, INTENT (OUT) :: enough(atoms%natd)
    COMPLEX, INTENT (INOUT) :: bascof_lo(3,-atoms%llod:atoms%llod,4*atoms%llod+2,atoms%nlod,atoms%natd)
    INTEGER, INTENT (INOUT) :: nkvec(atoms%nlod,atoms%natd)

    !     ..
    !     .. Local Scalars ..
    COMPLEX ctmp,term1
    REAL,PARAMETER:: linindq=1.0e-4,eps=1.0e-30
    INTEGER i,l,ll1,lm,lo ,mind,nbasf,na2,lmp,m
    LOGICAL linind
    !     ..
    !     .. Local Arrays ..
    COMPLEX clotmp(-atoms%llod:atoms%llod)
    !     ..
    enough(na) = .TRUE.
    term1 = con1 * ((atoms%rmt(ntyp)**2)/2) * phase

    !---> the whole program is in hartree units, therefore 1/wronskian is
    !---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
    !---> and c coefficients, is included in the t-matrices. thus, it does
    !---> not show up in the formula above.
    DO lo = 1,atoms%nlo(ntyp)
       l = atoms%llo(lo,ntyp)
       IF (.NOT.((s.LE.eps).AND.(l.GE.1))) THEN
          IF (atoms%invsat(na).EQ.0) THEN

             IF ((nkvec(lo,na)).LT. (2*atoms%llo(lo,ntyp)+1)) THEN
                enough(na) = .FALSE.
                nkvec(lo,na) = nkvec(lo,na) + 1
                nbasf = nbasf0(lo,na) + nkvec(lo,na)
                l = atoms%llo(lo,ntyp)
                ll1 = l* (l+1)
                DO m = -l,l
                   clotmp(m) = term1*CONJG(ylm(ll1+m+1))
                END DO
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                   DO m = -l,l
                      lm = ll1 + m
                      !WRITE(*,*) 'nkvec(lo,na)',nkvec(lo,na)
                      bascof_lo(1,m,nkvec(lo,na),lo,na) =&
                           &                                           clotmp(m)*alo1(lo,ntyp)
                      bascof_lo(2,m,nkvec(lo,na),lo,na) =&
                           &                                           clotmp(m)*blo1(lo,ntyp)
                      bascof_lo(3,m,nkvec(lo,na),lo,na) =&
                           &                                           clotmp(m)*clo1(lo,ntyp)
                   END DO
                ELSE
                   nkvec(lo,na) = nkvec(lo,na) - 1
                ENDIF ! linind
             ENDIF   ! nkvec < 2*atoms%llo

          ELSEIF (atoms%invsat(na).EQ.1) THEN
             IF ((nkvec(lo,na)).LT. (2* (2*atoms%llo(lo,ntyp)+1))) THEN
                enough(na) = .FALSE.
                nkvec(lo,na) = nkvec(lo,na) + 1
                nbasf = nbasf0(lo,na) + nkvec(lo,na)
                l = atoms%llo(lo,ntyp)
                ll1 = l* (l+1)
                DO m = -l,l
                   clotmp(m) = term1*CONJG(ylm(ll1+m+1))
                END DO
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                   DO m = -l,l
                      lm = ll1 + m
                      bascof_lo(1,m,nkvec(lo,na),lo,na) =&
                           &                                           clotmp(m)*alo1(lo,ntyp)
                      bascof_lo(2,m,nkvec(lo,na),lo,na) =&
                           &                                           clotmp(m)*blo1(lo,ntyp)
                      bascof_lo(3,m,nkvec(lo,na),lo,na) =&
                           &                                           clotmp(m)*clo1(lo,ntyp)
                   ENDDO  ! m
                ELSE       
                   nkvec(lo,na) = nkvec(lo,na) - 1
                ENDIF       ! linind
             ENDIF         ! nkvec < 2*atoms%llo
          ELSE
             STOP 'invsat =/= 0 or 1'
          ENDIF
       ELSE
          enough(na) = .FALSE.
       ENDIF  ! s > eps  & l >= 1
    END DO
    IF ((k.EQ.nv) .AND. (.NOT.enough(na))) THEN
       WRITE (6,FMT=*)&
            &     'abclocdn did not find enough linearly independent'
       WRITE (6,FMT=*)&
            &     'ccof coefficient-vectors. the linear independence'
       WRITE (6,FMT=*) 'quality, linindq, is set to: ',linindq,'.'
       WRITE (6,FMT=*) 'this value might be to large.'
       STOP 'abclocdn: did not find enough lin. ind. ccof-vectors'
    END IF

  END SUBROUTINE abclocdn1
END MODULE m_abclocdn1
