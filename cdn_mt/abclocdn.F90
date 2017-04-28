!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abclocdn
  USE m_juDFT
  !*********************************************************************
  ! Calculates the (upper case) A, B and C coefficients for the local
  ! orbitals. The difference to abccoflo is, that a summation over the
  ! Gs ist performed. The A, B and C coeff. are set up for each eigen-
  ! state.
  ! Philipp Kurz 99/04
  !*********************************************************************
  !*************** ABBREVIATIONS ***************************************
  ! nkvec   : stores the number of G-vectors that have been found and
  !           accepted during the construction of the local orbitals.
  ! kvec    : k-vector used in hssphn to attach the local orbital 'lo'
  !           of atom 'na' to it.
  !*********************************************************************
CONTAINS
  SUBROUTINE abclocdn(atoms, sym, noco,ccchi,kspin,iintsp,con1,phase,ylm,&
       ntyp,na,k,s,nv,ne,nbasf0,alo1,blo1,clo1,kvec,nkvec,enough,acof,bcof,ccof,zMat)
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_zMat),INTENT(IN)   :: zMat
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: kspin,iintsp
    INTEGER, INTENT (IN) :: k,na,ne,ntyp,nv
    REAL,    INTENT (IN) :: con1 ,s
    COMPLEX, INTENT (IN) :: phase
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nbasf0(atoms%nlod,atoms%nat) 
    REAL,    INTENT (IN) :: alo1(atoms%nlod,atoms%ntype),blo1(atoms%nlod,atoms%ntype)
    REAL,    INTENT (IN) :: clo1(atoms%nlod,atoms%ntype)
    COMPLEX, INTENT (IN) :: ylm( (atoms%lmaxd+1)**2 )
    COMPLEX, INTENT (IN) :: ccchi(2)
    INTEGER, INTENT (IN) :: kvec(2*(2*atoms%llod+1),atoms%nlod )
    LOGICAL, INTENT (OUT) :: enough ! enough(na)
    COMPLEX, INTENT (INOUT) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: ccof(-atoms%llod:,:,:,:)!(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%nat)
    INTEGER, INTENT (INOUT) :: nkvec(atoms%nlod,atoms%nat)
    !     ..
    !     .. Local Scalars ..
    COMPLEX ctmp,term1
    REAL,PARAMETER:: eps=1.0e-30
    INTEGER i,l,ll1,lm,lo ,mind,nbasf,na2,lmp,m
    !     ..
    !     .. Local Arrays ..
    COMPLEX clotmp(-atoms%llod:atoms%llod)
    !     ..
    LOGICAL :: l_real
    l_real=zMat%l_real
    !     ..
    enough = .TRUE.
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
                enough = .FALSE.
                nkvec(lo,na) = nkvec(lo,na) + 1
                nbasf = nbasf0(lo,na) + nkvec(lo,na)
                l = atoms%llo(lo,ntyp)
                ll1 = l* (l+1)
                DO m = -l,l
                   clotmp(m) = term1*CONJG(ylm(ll1+m+1))
                END DO
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                   !                   write(*,'(i3,5(2f10.5,2x))')k,(z(nbasf,i),i=11,15)
                   DO i = 1,ne
                      DO m = -l,l
                         lm = ll1 + m
                         !+gu_con
                         IF (noco%l_noco) THEN
                            IF (noco%l_ss) THEN
                               ctmp = clotmp(m)*ccchi(iintsp)*zMat%z_c(kspin+nbasf,i)
                            ELSE
                               ctmp = clotmp(m)*( ccchi(1)*zMat%z_c(nbasf,i)+ccchi(2)*zMat%z_c(kspin+nbasf,i) )
                            ENDIF
                         ELSE
                            IF (l_real) THEN
                               ctmp = zMat%z_r(nbasf,i)*clotmp(m)
                            ELSE
                               ctmp = zMat%z_c(nbasf,i)*clotmp(m)
                            ENDIF
                         ENDIF
                         acof(i,lm,na) = acof(i,lm,na) +ctmp*alo1(lo,ntyp)
                         bcof(i,lm,na) = bcof(i,lm,na) +ctmp*blo1(lo,ntyp)
                         ccof(m,i,lo,na) = ccof(m,i,lo,na) +ctmp*clo1(lo,ntyp)
                      END DO
                   END DO
                   !                  write(6,9000) nbasf,k,lo,na,
                   !     +                          (clo1(lo,ntyp)*clotmp(m),m=-l,l)
                   ! 9000             format(2i4,2i2,7(' (',e9.3,',',e9.3,')'))
                ELSE
                   nkvec(lo,na) = nkvec(lo,na) - 1
                ENDIF ! kvec = k
             ENDIF   ! nkvec < 2*atoms%llo

          ELSEIF (atoms%invsat(na).EQ.1) THEN
             IF ((nkvec(lo,na)).LT. (2* (2*atoms%llo(lo,ntyp)+1))) THEN
                enough = .FALSE.
                nkvec(lo,na) = nkvec(lo,na) + 1
                nbasf = nbasf0(lo,na) + nkvec(lo,na)
                l = atoms%llo(lo,ntyp)
                ll1 = l* (l+1)
                DO m = -l,l
                   clotmp(m) = term1*CONJG(ylm(ll1+m+1))
                END DO
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                   !                  write(*,*)'k vector nr ',k,' has been accepted'
                   !                  write(*,'(i3,5(2f10.5,2x))')k,(z(nbasf,i),i=11,15)
                   DO i = 1,ne
                      DO m = -l,l
                         lm = ll1 + m
                         !                        if(i.eq.1 .and. l.eq.1) then
                         !              write(*,*)'k=',k,' z=',z(nbasf,i),' clotmp=',clotmp(m)
                         !              write(*,*)'clo1=',clo1(lo,ntyp),' term1=',term1
                         !                         endif
                         !+gu_con
                         IF (noco%l_noco) THEN
                            IF (noco%l_ss) THEN
                               ctmp = clotmp(m)*ccchi(iintsp)*zMat%z_c(kspin+nbasf,i)
                            ELSE
                               ctmp = clotmp(m)*( ccchi(1)*zMat%z_c(nbasf,i)+ ccchi(2)*zMat%z_c(kspin+nbasf,i) )
                            ENDIF
                         ELSE
                            IF (l_real) THEN
                               ctmp = zMat%z_r(nbasf,i)*clotmp(m)
                            ELSE
                               ctmp = zMat%z_c(nbasf,i)*clotmp(m)
                            ENDIF
                         ENDIF
                         acof(i,lm,na) = acof(i,lm,na) +ctmp*alo1(lo,ntyp)
                         bcof(i,lm,na) = bcof(i,lm,na) +ctmp*blo1(lo,ntyp)
                         ccof(m,i,lo,na) = ccof(m,i,lo,na) +ctmp*clo1(lo,ntyp)
                         IF (noco%l_soc.AND.sym%invs) THEN
                            ctmp = zMat%z_c(nbasf,i)*CONJG(clotmp(m))*(-1)**(l-m)
                            na2 = sym%invsatnr(na)
                            lmp = ll1 - m
                            acof(i,lmp,na2) = acof(i,lmp,na2) +ctmp*alo1(lo,ntyp)
                            bcof(i,lmp,na2) = bcof(i,lmp,na2) +ctmp*blo1(lo,ntyp)
                            ccof(-m,i,lo,na2) = ccof(-m,i,lo,na2) +ctmp*clo1(lo,ntyp)
                         ENDIF
                      ENDDO  ! m
                   ENDDO     ! i = 1,ne
                ELSE       
                   nkvec(lo,na) = nkvec(lo,na) - 1
                ENDIF       ! kvec = k
             ENDIF         ! nkvec < 2*atoms%llo
          ELSE
             CALL juDFT_error("invsat =/= 0 or 1",calledby ="abclocdn")
          ENDIF
       ELSE
          enough = .FALSE.
       ENDIF  ! s > eps  & l >= 1
    END DO
    IF ((k.EQ.nv) .AND. (.NOT.enough)) THEN
       WRITE (6,FMT=*) 'abclocdn did not find enough linearly independent'
       WRITE (6,FMT=*) 'ccof coefficient-vectors.'
       CALL juDFT_error("did not find enough lin. ind. ccof-vectors" ,calledby ="abclocdn")
    END IF

  END SUBROUTINE abclocdn
END MODULE m_abclocdn
