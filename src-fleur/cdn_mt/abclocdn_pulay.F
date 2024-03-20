      MODULE m_abclocdnpulay
      USE m_fleurenv
      CONTAINS
      SUBROUTINE abclocdn_pulay(
     >                          nobd,natd,nlod,llod,lmaxd,neigd,ntypd,
     >                          lmd,nbasfcn,nlo,llo,invsat,invsatnr,
     >                          l_noco,ccchi,kspin,l_ss,iintsp,
     >                          con1,rmt,phase,ylm,ntyp,na,k,fgp,
     >                          s,nv,ne,z,nbasf0,alo1,blo1,clo1,
     <                          kvec,nkvec,enough,acof,bcof,ccof,
     X                          acoflo,bcoflo,aveccof,bveccof,cveccof)
c
c*********************************************************************
c for details see abclocdn; called by to_pulay
c*********************************************************************
c
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: nobd,natd,nlod,llod,lmaxd,neigd,ntypd
      INTEGER, INTENT (IN) :: lmd,nbasfcn,iintsp
      INTEGER, INTENT (IN) :: k,na,ne,ntyp,nv,kspin
      REAL,    INTENT (IN) :: con1,rmt,s
      COMPLEX, INTENT (IN) :: phase
      LOGICAL, INTENT (IN) :: l_noco,l_ss
C     ..
C     .. Array Arguments ..
      INTEGER, INTENT (IN) :: nbasf0(nlod,natd),invsat(natd)
      INTEGER, INTENT (IN) :: nlo(ntypd),llo(nlod,ntypd),invsatnr(natd)
      INTEGER, INTENT (IN) :: kvec(2*(2*llod+1),nlod)
      REAL,    INTENT (IN) :: alo1(nlod,ntypd),blo1(nlod,ntypd)
      REAL,    INTENT (IN) :: clo1(nlod,ntypd)
      REAL,    INTENT (IN) :: fgp(3)
      COMPLEX, INTENT (IN) :: ylm( (lmaxd+1)**2 ),ccchi(2)
      COMPLEX, INTENT (INOUT) :: acof(nobd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: bcof(nobd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: ccof(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (INOUT) :: acoflo(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (INOUT) :: bcoflo(-llod:llod,nobd,nlod,natd)
      COMPLEX, INTENT (INOUT) :: aveccof(3,nobd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: bveccof(3,nobd,0:lmd,natd)
      COMPLEX, INTENT (INOUT) :: cveccof(3,-llod:llod,nobd,nlod,natd)
      LOGICAL, INTENT (OUT) :: enough(natd)
      INTEGER, INTENT (INOUT) :: nkvec(nlod,natd)
#if ( defined(CPP_INVERSION) && !defined(CPP_SOC) )
      REAL,    INTENT (IN) :: z(nbasfcn,neigd)
#else
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#endif
C     ..
C     .. Local Scalars ..
      COMPLEX ctmp,term1
      REAL linindq,eps
      INTEGER i,ie,l,ll1,lm,lo,m,mind,nbasf,na2,lmp
      LOGICAL linind
C     ..
C     .. Local Arrays ..
      COMPLEX clotmp(-llod:llod)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC conjg
C     ..
C     .. Data statements ..

      DATA linindq/1.0e-4/,eps/1.0e-30/
C     ..

      enough(na) = .true.
      term1 = con1* ((rmt**2)/2)*phase
c
c---> the whole program is in hartree units, therefore 1/wronskian is
c---> (rmt**2)/2. the factor i**l, which usually appears in the a, b
c---> and c coefficients, is included in the t-matrices. thus, it does
c---> not show up in the formula above.
c
      DO lo = 1,nlo(ntyp)
         l = llo(lo,ntyp)
         IF (.NOT.((s.LE.eps).AND.(l.GE.1))) THEN
            IF (invsat(na).EQ.0) THEN
               IF ((nkvec(lo,na)).LT. (2*llo(lo,ntyp)+1)) THEN
                  enough(na) = .false.
                  nkvec(lo,na) = nkvec(lo,na) + 1
                  nbasf = nbasf0(lo,na) + nkvec(lo,na)
                  l = llo(lo,ntyp)
                  ll1 = l* (l+1)
                  DO m = -l,l
                     clotmp(m) = term1*conjg(ylm(ll1+m+1))
                  END DO
c
                  IF ( kvec(nkvec(lo,na),lo) == k ) THEN
                     DO ie = 1,ne
                        DO m = -l,l
                           lm = ll1 + m
                           IF (l_noco) THEN
                             IF (l_ss) THEN
                               ctmp = clotmp(m)*
     +                                ccchi(iintsp)*z(kspin+nbasf,ie)
                             ELSE
                              ctmp = clotmp(m)*( ccchi(1)*z(nbasf,ie)
     +                                 +   ccchi(2)*z(kspin+nbasf,ie) )
                             ENDIF
                           ELSE
                              ctmp = z(nbasf,ie)*clotmp(m)
                           ENDIF
                           acof(ie,lm,na)     = acof(ie,lm,na) +
     +                                         ctmp*alo1(lo,ntyp)
                           bcof(ie,lm,na)     = bcof(ie,lm,na) +
     +                                         ctmp*blo1(lo,ntyp)
                           ccof(m,ie,lo,na)   = ccof(m,ie,lo,na) +
     +                                         ctmp*clo1(lo,ntyp)
                           acoflo(m,ie,lo,na) = acoflo(m,ie,lo,na) +
     +                                         ctmp*alo1(lo,ntyp)
                           bcoflo(m,ie,lo,na) = bcoflo(m,ie,lo,na) +
     +                                         ctmp*blo1(lo,ntyp)
                           DO i = 1,3
                              aveccof(i,ie,lm,na)=aveccof(i,ie,lm,na) +
     +                                        fgp(i)*ctmp*alo1(lo,ntyp)
                              bveccof(i,ie,lm,na)=bveccof(i,ie,lm,na) +
     +                                        fgp(i)*ctmp*blo1(lo,ntyp)
                              cveccof(i,m,ie,lo,na) =
     +                        cveccof(i,m,ie,lo,na) +
     +                                        fgp(i)*ctmp*clo1(lo,ntyp)
                           ENDDO
                        END DO
                     END DO
c                    write(6,9000) nbasf,k,lo,na,
c     +                          (clo1(lo,ntyp)*clotmp(m),m=-l,l)
c 9000               format(2i4,2i2,7(' (',e9.3,',',e9.3,')'))
                  ELSE
                     nkvec(lo,na) = nkvec(lo,na) - 1
                  END IF
               END IF
            ELSEIF (invsat(na).EQ.1) THEN
               IF ((nkvec(lo,na)).LT. (2* (2*llo(lo,ntyp)+1))) THEN
                  enough(na) = .false.
                  nkvec(lo,na) = nkvec(lo,na) + 1
                  nbasf = nbasf0(lo,na) + nkvec(lo,na)
                  l = llo(lo,ntyp)
                  ll1 = l* (l+1)
                  DO m = -l,l
                     clotmp(m) = term1*conjg(ylm(ll1+m+1))
                  END DO
c
                  IF ( kvec(nkvec(lo,na),lo) == k ) THEN
c                     write(*,*)'k vector nr ',k,' has been accepted'
                     DO ie = 1,ne
                        DO m = -l,l
                           lm = ll1 + m
                           IF (l_noco) THEN
                            IF (l_ss) THEN
                              ctmp = clotmp(m)*
     +                               ccchi(iintsp)*z(kspin+nbasf,ie)
                            ELSE
                              ctmp = clotmp(m)*( ccchi(1)*z(nbasf,ie)
     +                                 +   ccchi(2)*z(kspin+nbasf,ie) )
                            ENDIF
                           ELSE
                              ctmp = z(nbasf,ie)*clotmp(m)
                           ENDIF
                           acof(ie,lm,na) = acof(ie,lm,na) +
     +                                      ctmp*alo1(lo,ntyp)
                           bcof(ie,lm,na) = bcof(ie,lm,na) +
     +                                      ctmp*blo1(lo,ntyp)
                           ccof(m,ie,lo,na) = ccof(m,ie,lo,na) +
     +                                      ctmp*clo1(lo,ntyp)
                           acoflo(m,ie,lo,na) = acoflo(m,ie,lo,na) +
     +                                         ctmp*alo1(lo,ntyp)
                           bcoflo(m,ie,lo,na) = bcoflo(m,ie,lo,na) +
     +                                         ctmp*blo1(lo,ntyp)
                           DO i = 1,3
                            aveccof(i,ie,lm,na)=aveccof(i,ie,lm,na) +
     +                                        fgp(i)*ctmp*alo1(lo,ntyp)
                            bveccof(i,ie,lm,na)=bveccof(i,ie,lm,na) +
     +                                        fgp(i)*ctmp*blo1(lo,ntyp)
                            cveccof(i,m,ie,lo,na)=cveccof(i,m,ie,lo,na)+
     +                                        fgp(i)*ctmp*clo1(lo,ntyp)
                           ENDDO
#if ( defined(CPP_SOC) && defined(CPP_INVERSION) )
                           ctmp = z(nbasf,ie) * 
     +                                     conjg(clotmp(m))*(-1)**(l-m)
                           na2 = invsatnr(na)
                           lmp = ll1 - m
                           acof(ie,lmp,na2) = acof(ie,lmp,na2) +
     +                                               ctmp*alo1(lo,ntyp)
                           bcof(ie,lmp,na2) = bcof(ie,lmp,na2) +
     +                                               ctmp*blo1(lo,ntyp)
                           ccof(-m,ie,lo,na2) = ccof(-m,ie,lo,na2) +
     +                                               ctmp*clo1(lo,ntyp)
                           acoflo(-m,ie,lo,na2) = acoflo(-m,ie,lo,na2) +
     +                                                ctmp*alo1(lo,ntyp)
                           bcoflo(-m,ie,lo,na2) = bcoflo(-m,ie,lo,na2) +
     +                                                ctmp*blo1(lo,ntyp)
                           DO i = 1,3
                            aveccof(i,ie,lmp,na2)=aveccof(i,ie,lmp,na2)-
     +                                        fgp(i)*ctmp*alo1(lo,ntyp)
                            bveccof(i,ie,lmp,na2)=bveccof(i,ie,lmp,na2)-
     +                                        fgp(i)*ctmp*blo1(lo,ntyp)
                            cveccof(i,-m,ie,lo,na2) =
     +                      cveccof(i,-m,ie,lo,na2) -
     +                                        fgp(i)*ctmp*clo1(lo,ntyp)
                           ENDDO
#endif
                        ENDDO ! loop over m
                     ENDDO    ! loop over eigenstates (ie)
                  ELSE
                     nkvec(lo,na) = nkvec(lo,na) - 1
                  END IF   ! linind
               END IF      ! nkvec(lo,na) < 2*(2*llo + 1)
            ELSE
               CALL fleur_err("invsat =/= 0 or 1",calledby
     +              ="abclocdn_pulay")
            ENDIF
         ELSE
            enough(na) = .false.
         ENDIF ! s > eps  & l >= 1  
      END DO
      IF ((k.EQ.nv) .AND. (.NOT.enough(na))) THEN
         WRITE (6,FMT=*)
     +     'abclocdn did not find enough linearly independent'
         WRITE (6,FMT=*)
     +     'ccof coefficient-vectors. the linear independence'
         WRITE (6,FMT=*) 'quality, linindq, is set to: ',linindq,'.'
         WRITE (6,FMT=*) 'this value might be to large.'
         CALL fleur_err("did not find enough lin. ind. ccof-vectors"
     +        ,calledby ="abclocdn_pulay")
      END IF

      END SUBROUTINE abclocdn_pulay
      END MODULE m_abclocdnpulay
