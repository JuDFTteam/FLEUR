!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_abclocdnpulay
  USE m_juDFT
CONTAINS
  SUBROUTINE abclocdn_pulay(&
       &                          atoms,sym,&
       &                          noco,ccchi,kspin,iintsp,&
       &                          con1,phase,ylm,ntyp,na,k,fgp,&
       &                          s,nv,ne,nbasf0,alo1,blo1,clo1,&
       &                          kvec,nkvec,enough,acof,bcof,ccof,&
       &                          acoflo,bcoflo,aveccof,bveccof,cveccof,zMat,realdata)
    !
    !*********************************************************************
    ! for details see abclocdn; calles by to_pulay
    !*********************************************************************
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_zMat),INTENT(IN)   :: zMat
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: iintsp
    INTEGER, INTENT (IN) :: k,na,ne,ntyp,nv,kspin
    REAL,    INTENT (IN) :: con1 ,s
    COMPLEX, INTENT (IN) :: phase
    !     ..
    !     .. Array Arguments ..
    INTEGER, INTENT (IN) :: nbasf0(atoms%nlod,atoms%nat) 
    INTEGER, INTENT (IN) :: kvec(2*(2*atoms%llod+1),atoms%nlod)
    REAL,    INTENT (IN) :: alo1(atoms%nlod,atoms%ntype),blo1(atoms%nlod,atoms%ntype)
    REAL,    INTENT (IN) :: clo1(atoms%nlod,atoms%ntype)
    REAL,    INTENT (IN) :: fgp(3)
    COMPLEX, INTENT (IN) :: ylm( (atoms%lmaxd+1)**2 ),ccchi(2)
    COMPLEX, INTENT (INOUT) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: ccof(-atoms%llod:,:,:,:)!(-llod:llod,nobd,atoms%nlod,atoms%nat)
    COMPLEX, INTENT (INOUT) :: acoflo(-atoms%llod:,:,:,:)
    COMPLEX, INTENT (INOUT) :: bcoflo(-atoms%llod:,:,:,:)
    COMPLEX, INTENT (INOUT) :: aveccof(:,:,0:,:)!(3,nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: bveccof(:,:,0:,:)!(3,nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (INOUT) :: cveccof(:,-atoms%llod:,:,:,:)!(3,-atoms%llod:llod,nobd,atoms%nlod,atoms%nat)
    LOGICAL, INTENT (OUT) :: enough(atoms%nat)
    INTEGER, INTENT (INOUT) :: nkvec(atoms%nlod,atoms%nat)
    LOGICAL,OPTIONAL,INTENT(IN) ::realdata
    !     ..
    !     .. Local Scalars ..
    COMPLEX ctmp,term1
    REAL,PARAMETER:: linindq=1.0e-4,eps=1.e-30
    INTEGER i,ie,l,ll1,lm,lo ,mind,nbasf,na2,lmp,m
    LOGICAL linind
    !     ..
    !     .. Local Arrays ..
    COMPLEX clotmp(-atoms%llod:atoms%llod)
    !     ..
    LOGICAL:: l_real
    l_real=zMat%l_real
    IF (PRESENT(realdata)) l_real=realdata
    enough(na) = .TRUE.
    term1 = con1* ((atoms%rmt(ntyp)**2)/2)*phase
    !
    !---> the whole program is in hartree units, therefore 1/wronskian is
    !---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
    !---> and c coefficients, is included in the t-matrices. thus, it does
    !---> not show up in the formula above.
    !
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
                !
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                   DO ie = 1,ne
                      DO m = -l,l
                         lm = ll1 + m
                         IF (noco%l_noco) THEN
                            IF (noco%l_ss) THEN
                               ctmp = clotmp(m)* ccchi(iintsp)*zMat%z_c(kspin+nbasf,ie)
                            ELSE
                               ctmp = clotmp(m)*( ccchi(1)*zMat%z_c(nbasf,ie)+ccchi(2)*zMat%z_c(kspin+nbasf,ie) )
                            ENDIF
                         ELSE
                            IF (l_real) THEN
                               ctmp = zMat%z_r(nbasf,ie)*clotmp(m)
                            ELSE
                               ctmp = zMat%z_c(nbasf,ie)*clotmp(m)
                            ENDIF
                         ENDIF
                         acof(ie,lm,na)     = acof(ie,lm,na) +ctmp*alo1(lo,ntyp)
                         bcof(ie,lm,na)     = bcof(ie,lm,na) +ctmp*blo1(lo,ntyp)
                         ccof(m,ie,lo,na)   = ccof(m,ie,lo,na) +ctmp*clo1(lo,ntyp)
                         acoflo(m,ie,lo,na) = acoflo(m,ie,lo,na) +ctmp*alo1(lo,ntyp)
                         bcoflo(m,ie,lo,na) = bcoflo(m,ie,lo,na) +ctmp*blo1(lo,ntyp)
                         DO i = 1,3
                            aveccof(i,ie,lm,na)=aveccof(i,ie,lm,na) +fgp(i)*ctmp*alo1(lo,ntyp)
                            bveccof(i,ie,lm,na)=bveccof(i,ie,lm,na) +fgp(i)*ctmp*blo1(lo,ntyp)
                            cveccof(i,m,ie,lo,na) =cveccof(i,m,ie,lo,na) +fgp(i)*ctmp*clo1(lo,ntyp)
                         ENDDO
                      END DO
                   END DO
                   !                    write(6,9000) nbasf,k,lo,na,
                   !     +                          (clo1(lo,ntyp)*clotmp(m),m=-l,l)
                   ! 9000               format(2i4,2i2,7(' (',e9.3,',',e9.3,')'))
                ELSE
                   nkvec(lo,na) = nkvec(lo,na) - 1
                END IF
             END IF
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
                !
                IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                   !                     write(*,*)'k vector nr ',k,' has been accepted'
                   DO ie = 1,ne
                      DO m = -l,l
                         lm = ll1 + m
                         IF (noco%l_noco) THEN
                            IF (noco%l_ss) THEN
                               ctmp = clotmp(m)*ccchi(iintsp)*zMat%z_c(kspin+nbasf,ie)
                            ELSE
                               ctmp = clotmp(m)*( ccchi(1)*zMat%z_c(nbasf,ie)+ccchi(2)*zMat%z_c(kspin+nbasf,ie) )
                            ENDIF
                         ELSE
                            IF (l_real) THEN
                               ctmp = zMat%z_r(nbasf,ie)*clotmp(m)
                            ELSE
                               ctmp = zMat%z_c(nbasf,ie)*clotmp(m)
                            END IF
                         ENDIF
                         acof(ie,lm,na) = acof(ie,lm,na) +ctmp*alo1(lo,ntyp)
                         bcof(ie,lm,na) = bcof(ie,lm,na) +ctmp*blo1(lo,ntyp)
                         ccof(m,ie,lo,na) = ccof(m,ie,lo,na) +ctmp*clo1(lo,ntyp)
                         acoflo(m,ie,lo,na) = acoflo(m,ie,lo,na) +ctmp*alo1(lo,ntyp)
                         bcoflo(m,ie,lo,na) = bcoflo(m,ie,lo,na) +ctmp*blo1(lo,ntyp)
                         DO i = 1,3
                            aveccof(i,ie,lm,na)=aveccof(i,ie,lm,na) +fgp(i)*ctmp*alo1(lo,ntyp)
                            bveccof(i,ie,lm,na)=bveccof(i,ie,lm,na) +fgp(i)*ctmp*blo1(lo,ntyp)
                            cveccof(i,m,ie,lo,na)=cveccof(i,m,ie,lo,na)+fgp(i)*ctmp*clo1(lo,ntyp)
                         ENDDO
                         IF (noco%l_soc.AND.sym%invs) THEN
                            ctmp = zMat%z_c(nbasf,ie) * CONJG(clotmp(m))*(-1)**(l-m)
                            na2 = sym%invsatnr(na)
                            lmp = ll1 - m
                            acof(ie,lmp,na2) = acof(ie,lmp,na2) +ctmp*alo1(lo,ntyp)
                            bcof(ie,lmp,na2) = bcof(ie,lmp,na2) +ctmp*blo1(lo,ntyp)
                            ccof(-m,ie,lo,na2) = ccof(-m,ie,lo,na2) + ctmp*clo1(lo,ntyp)
                            acoflo(-m,ie,lo,na2) = acoflo(-m,ie,lo,na2) +ctmp*alo1(lo,ntyp)
                            bcoflo(-m,ie,lo,na2) = bcoflo(-m,ie,lo,na2) +ctmp*blo1(lo,ntyp)
                            DO i = 1,3
                               aveccof(i,ie,lmp,na2)=aveccof(i,ie,lmp,na2)-fgp(i)*ctmp*alo1(lo,ntyp)
                               bveccof(i,ie,lmp,na2)=bveccof(i,ie,lmp,na2)-fgp(i)*ctmp*blo1(lo,ntyp)
                               cveccof(i,-m,ie,lo,na2) =cveccof(i,-m,ie,lo,na2) -fgp(i)*ctmp*clo1(lo,ntyp)
                            ENDDO
                         ENDIF
                      ENDDO ! loop over m
                   ENDDO    ! loop over eigenstates (ie)
                ELSE
                   nkvec(lo,na) = nkvec(lo,na) - 1
                END IF   ! linind
             END IF      ! nkvec(lo,na) < 2*(2*atoms%llo + 1)
          ELSE
             CALL juDFT_error("invsat =/= 0 or 1",calledby ="abclocdn_pulay")
          ENDIF
       ELSE
          enough(na) = .FALSE.
       ENDIF ! s > eps  & l >= 1  
    END DO
    IF ((k.EQ.nv) .AND. (.NOT.enough(na))) THEN
       WRITE (6,FMT=*) 'abclocdn did not find enough linearly independent'
       WRITE (6,FMT=*) 'ccof coefficient-vectors. the linear independence'
       WRITE (6,FMT=*) 'quality, linindq, is set to: ',linindq,'.'
       WRITE (6,FMT=*) 'this value might be to large.'
       CALL juDFT_error("did not find enough lin. ind. ccof-vectors" ,calledby ="abclocdn_pulay")
    END IF

  END SUBROUTINE abclocdn_pulay
END MODULE m_abclocdnpulay
