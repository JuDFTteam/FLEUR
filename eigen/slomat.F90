!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_slomat
  !***********************************************************************
  ! updates the overlap matrix with the contributions from the local
  ! orbitals.
  !                                                p.kurz sept. 1996
  !***********************************************************************
CONTAINS
  SUBROUTINE slomat(&
       input,atoms, ntyp,na,lapw,con1,n_size,n_rank,&
       gk,rph,cph,f,g,kvec,usp,ud, alo1,blo1,clo1,noco,&
       ab_dim,iintsp,jintsp,chi11,chi22,chi21, iilo,locol,nkvecprevat,l_real,bb_r,bb_c)
    !***********************************************************************
    ! iilo is the starting position for the LO-contribution of atom 'na'
    !      in the Hamilton-matrix minus 1 (starting from (nv+1)*nv/2)
    ! nkvecprevat gives the number of G-vectors already processed on other
    !             atoms (each l-LO gives (2l+1)*invsfact G-vectors)
    ! nkvecprevlo gives the number of G-vectors already processed on this
    !             atoms due to other lo's (appears in hlomat only)
    ! locol stores the number of columns already processed; on parallel
    !       computers this decides, whether the LO-contribution is
    !       done on this node                                          gb00
    !
    ! function legpol() at end of module
    !***********************************************************************
#include"cpp_double.h" 
    USE m_constants,ONLY: fpi_const
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN)   :: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,ntyp,n_size,n_rank 
    INTEGER, INTENT (IN) :: ab_dim,iintsp,jintsp
    COMPLEX, INTENT (IN) :: chi11,chi22,chi21
    REAL,    INTENT (IN) :: con1
    INTEGER, INTENT (INOUT) :: iilo,nkvecprevat,locol
    INTEGER, INTENT(IN)  :: usp
    !     ..
    !     .. Array Arguments ..
    INTEGER,INTENT (IN) :: kvec(:,:)!(2*llod_1),nlod  )
    REAL,   INTENT (IN) :: alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod) 
    REAL,   INTENT (IN) :: gk(:,:,:)!(dimension%nvd,3,ab_dim)
    REAL,   INTENT (IN) :: rph(:,:),cph(:,:)!(nvd,ab_dim)
    REAL,   INTENT (IN) :: f(:,0:,:,:)!(dimension%nvd,0:atoms%lmaxd,atoms%ntypd,ab_dim)
    REAL,   INTENT (IN) :: g(:,0:,:,:)!(dimension%nvd,0:atoms%lmaxd,atoms%ntypd,ab_dim)
    TYPE(t_usdus),INTENT(IN):: ud
    REAL,  ALLOCATABLE,  OPTIONAL,INTENT (INOUT) :: bb_r(:)!(matsize)
    COMPLEX,ALLOCATABLE, OPTIONAL,INTENT (INOUT) :: bb_c(:)!(matsize)
    LOGICAL,INTENT(IN)                  :: l_real
    !     ..
    !     .. Local Scalars ..
    REAL con,dotp,fact1,fact2,fact3,fl2p1
    INTEGER invsfct,k ,l,l2p1,lo,lop,lp,nkvec,nkvecat,nkvecp,kp
    INTEGER iilo_s,locol_s,nkvecprevat_s,ic,ii,ij,n
    COMPLEX chihlp
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: bhelp(:)
    !     ..
    !$OMP MASTER

    IF ((atoms%invsat(na).EQ.0) .OR. (atoms%invsat(na).EQ.1)) THEN
       !--->    if this atom is the first of two atoms related by inversion,
       !--->    the contributions to the overlap matrix of both atoms are added
       !--->    at once. where it is made use of the fact, that the sum of
       !--->    these contributions is twice the real part of the contribution
       !--->    of each atom. note, that in this case there are twice as many
       !--->    (2*(2*l+1)) k-vectors (compare abccoflo and comments there).
       IF (atoms%invsat(na).EQ.0) invsfct = 1
       IF (atoms%invsat(na).EQ.1) invsfct = 2
       !+noco
       IF (noco%l_ss) THEN
          iilo_s = iilo
          locol_s = locol
          nkvecprevat_s = nkvecprevat ! save for other interstitial spin loops
          ic = 0                                        ! update b-matrix
          DO lo = 1,atoms%nlo(ntyp)
             ic = ic + invsfct* (2*atoms%llo(lo,ntyp)+1)
          ENDDO
          k = ic*(lapw%nv(jintsp)+nkvecprevat) + (ic+1)*ic/2
          ALLOCATE ( bhelp(k) )       ! initialize help-array 
          bhelp=CMPLX(0.,0.)
          iilo = 0
       ENDIF
       !-noco
       con = con1* ((atoms%rmt(ntyp))**2)/2.0
       nkvecat = 0
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          l2p1 = 2*l + 1
          fl2p1 = l2p1/fpi_const
          IF (input%l_useapw) THEN
             fact1 = (con**2)* fl2p1 * (&
                  alo1(lo)* (  alo1(lo) + &
                  2*clo1(lo)*   ud%uulon(lo,ntyp,usp) ) +&
                  blo1(lo)* (  blo1(lo)*     ud%ddn(l,ntyp,usp) +&
                  2*clo1(lo)*   ud%dulon(lo,ntyp,usp) ) +&
                  clo1(lo)*    clo1(lo)*ud%uloulopn(lo,lo,ntyp,usp) )
          ELSE
             fact1 = (con**2)* fl2p1 * (&
                  alo1(lo)* (  alo1(lo) + &
                  2*clo1(lo) * ud%uulon(lo,ntyp,usp) ) +&
                  blo1(lo)* (  blo1(lo) * ud%ddn(l, ntyp,usp) +&
                  2*clo1(lo) * ud%dulon(lo,ntyp,usp) ) +&
                  clo1(lo)*    clo1(lo) )
          ENDIF
          DO nkvec = 1,invsfct* (2*l+1)
             !+t3e
             locol = locol + 1
             IF (MOD(locol-1,n_size).EQ.n_rank) THEN
                !-t3e
                k = kvec(nkvec,lo)
                !--->          calculate the overlap matrix elements with the regular
                !--->          flapw basis-functions
                DO kp = 1,lapw%nv(jintsp)
                   iilo = iilo + 1
                   fact2 = con * fl2p1 * (&
                        f(kp,l,ntyp,jintsp)* ( alo1(lo) + &
                        clo1(lo)*ud%uulon(lo,ntyp,usp))+&
                        g(kp,l,ntyp,jintsp)* ( blo1(lo) * ud%ddn(l,ntyp,usp)+&
                        clo1(lo)*ud%dulon(lo,ntyp,usp)))
                   dotp = gk(k,1,iintsp) * gk(kp,1,jintsp) + &
                        gk(k,2,iintsp) * gk(kp,2,jintsp) +&
                        gk(k,3,iintsp) * gk(kp,3,jintsp)
                   IF (l_real) THEN
                      bb_r(iilo) = bb_r(iilo) + invsfct*fact2 * legpol(l,dotp) *&
                           ( rph(k,iintsp)*rph(kp,jintsp) +&
                           cph(k,iintsp)*cph(kp,jintsp) )
                   ELSE
                      IF (.NOT.noco%l_ss) THEN
                         bb_c(iilo) = bb_c(iilo) + invsfct*fact2 * legpol(l,dotp) *&
                              CMPLX( ( rph(k,iintsp)*rph(kp,jintsp) +&
                              cph(k,iintsp)*cph(kp,jintsp) ) ,&
                              ( cph(k,iintsp)*rph(kp,jintsp) -&
                              rph(k,iintsp)*cph(kp,jintsp) ) )
                      ELSE 
                         bhelp(iilo) = invsfct*fact2*legpol(l,dotp) *&
                              CMPLX( ( rph(k,iintsp)*rph(kp,jintsp) +&
                              cph(k,iintsp)*cph(kp,jintsp) ) ,&
                              ( cph(k,iintsp)*rph(kp,jintsp) -&
                              rph(k,iintsp)*cph(kp,jintsp) ) )
                      ENDIF
                   ENDIF
                END DO
                !--->          calculate the overlap matrix elements with other local
                !--->          orbitals at the same atom, if they have the same l
                iilo = iilo + nkvecprevat
                DO lop = 1, (lo-1)
                   lp = atoms%llo(lop,ntyp)
                   IF (l.EQ.lp) THEN
                      fact3 = con**2 * fl2p1 * (&
                           alo1(lop)*(alo1(lo) + &
                           clo1(lo)*ud%uulon(lo,ntyp,usp))+&
                           blo1(lop)*(blo1(lo)*ud%ddn(l,ntyp,usp) +&
                           clo1(lo)*ud%dulon(lo,ntyp,usp))+&
                           clo1(lop)*(alo1(lo)*ud%uulon(lop,ntyp,usp)+&
                           blo1(lo)*ud%dulon(lop,ntyp,usp)+&
                           clo1(lo)*ud%uloulopn(lop,lo,ntyp,usp)))
                      DO nkvecp = 1,invsfct* (2*lp+1)
                         iilo = iilo + 1
                         kp = kvec(nkvecp,lop) !
                         dotp = gk(k,1,iintsp) * gk(kp,1,jintsp) +&
                              gk(k,2,iintsp) * gk(kp,2,jintsp) +&
                              gk(k,3,iintsp) * gk(kp,3,jintsp)
                         IF (l_real) THEN
                            bb_r(iilo) =bb_r(iilo)+invsfct*fact3*legpol(l,dotp)* &
                                 ( rph(k,iintsp)*rph(kp,jintsp) +&
                                 cph(k,iintsp)*cph(kp,jintsp) )
                         ELSE
                            IF (.NOT.noco%l_ss) THEN
                               bb_c(iilo) =bb_c(iilo)+invsfct*fact3*legpol(l,dotp)*&
                                    CMPLX( ( rph(k,iintsp)*rph(kp,jintsp) +&
                                    cph(k,iintsp)*cph(kp,jintsp) ) ,&
                                    ( cph(k,iintsp)*rph(kp,jintsp) -&
                                    rph(k,iintsp)*cph(kp,jintsp) ) )
                            ELSE 
                               bhelp(iilo) = invsfct*fact3*legpol(l,dotp)*&
                                    CMPLX( ( rph(k,iintsp)*rph(kp,jintsp) +&
                                    cph(k,iintsp)*cph(kp,jintsp) ) ,&
                                    ( cph(k,iintsp)*rph(kp,jintsp) -&
                                    rph(k,iintsp)*cph(kp,jintsp) ) )
                            ENDIF
                         ENDIF
                      END DO
                   ELSE
                      iilo = iilo + invsfct* (2*lp+1)
                   END IF
                END DO
                !--->          calculate the overlap matrix elements of one local
                !--->          orbital with itself
                DO nkvecp = 1,nkvec
                   iilo = iilo + 1
                   kp = kvec(nkvecp,lo)
                   dotp = gk(k,1,iintsp) * gk(kp,1,jintsp) +&
                        gk(k,2,iintsp) * gk(kp,2,jintsp) +&
                        gk(k,3,iintsp) * gk(kp,3,jintsp)
                   IF (l_real) THEN
                      bb_r(iilo) = bb_r(iilo) + invsfct*fact1*legpol(l,dotp) *&
                           ( rph(k,iintsp)*rph(kp,jintsp) +&
                           cph(k,iintsp)*cph(kp,jintsp) )
                   ELSE
                      IF (.NOT.noco%l_ss) THEN
                         bb_c(iilo) = bb_c(iilo) + invsfct*fact1*legpol(l,dotp)*&
                              CMPLX( ( rph(k,iintsp)*rph(kp,jintsp) +&
                              cph(k,iintsp)*cph(kp,jintsp) ) ,&
                              ( cph(k,iintsp)*rph(kp,jintsp) -&
                              rph(k,iintsp)*cph(kp,jintsp) ) )
                      ELSE 
                         bhelp(iilo) = invsfct*fact1*legpol(l,dotp)*&
                              CMPLX( ( rph(k,iintsp)*rph(kp,jintsp) +&
                              cph(k,iintsp)*cph(kp,jintsp) ) ,&
                              ( cph(k,iintsp)*rph(kp,jintsp) -&
                              rph(k,iintsp)*cph(kp,jintsp) ) )
                      ENDIF
                   ENDIF
                END DO
             ENDIF ! mod(locol-1,n_size) = nrank 
             !-t3e
          END DO
          nkvecat = nkvecat + invsfct* (2*l+1)
       END DO
       nkvecprevat = nkvecprevat + nkvecat
       !+noco
       IF (noco%l_ss) THEN
          IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN     !---> spin-up spin-up part
             chihlp = chi11
             n = lapw%nv(1) + nkvecprevat_s
             ij = (n+1)*n/2
          ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN !---> spin-down spin-down part
             chihlp = chi22
             n = lapw%nv(1) + lapw%nv(iintsp) + atoms%nlotot + nkvecprevat_s
             ij = (n+1)*n/2 + lapw%nv(1) + atoms%nlotot
          ELSE                                          !---> spin-down spin-up part
             chihlp = chi21
             n = lapw%nv(1) + lapw%nv(iintsp) + atoms%nlotot + nkvecprevat_s
             ij = (n+1)*n/2 
          ENDIF

          ic = 0                                        ! update b-matrix
          DO lo = 1,atoms%nlo(ntyp)
             ic = ic + invsfct* (2*atoms%llo(lo,ntyp)+1)
          ENDDO
          IF (.NOT.( iintsp.EQ.1 .AND. jintsp.EQ.2 )) THEN
             ii = 0
             DO k = 1, ic
                n = k + lapw%nv(jintsp) + nkvecprevat_s
                bb_c(ij+1:ij+n-1)=bb_c(ij+1:ij+n-1)+chihlp*bhelp(ii+1:ii+n-1)
                IF (.NOT.(iintsp.EQ.2 .AND. jintsp.EQ.1 )) THEN
                   bb_c(ij+n)=bb_c(ij+n)+chihlp*bhelp(ii+n)
                ENDIF
                ii = ii + n
                ij = ij + n + (lapw%nv(3-jintsp)+atoms%nlotot)*(iintsp-1)
             ENDDO
          ELSE                                         ! special treatment for up-down:
             n = lapw%nv(1) + atoms%nlotot
             ii = 0
             DO k = 1, ic
                ij = (n+1)*n/2 + lapw%nv(1) + k + nkvecprevat_s
                DO kp = 1, k + lapw%nv(jintsp) + nkvecprevat_s
                   ii = ii + 1
                   bb_c(ij) = bb_c(ij) +  chihlp *CONJG( bhelp(ii) )
                   ij = ij + lapw%nv(1) + kp + atoms%nlotot
                ENDDO
             ENDDO
          ENDIF
          !          n1=nv(jintsp)+1
          !          n2=n1+nv(jintsp)+2
          !          n3=n2+nv(jintsp)+3
          !          do n = 1,nv(jintsp)+(nlotot/2)
          !           write(4,'(8f15.8)') bhelp(n),bhelp(n1+n),bhelp(n2+n),bhelp(n3+n)
          !          enddo
          DEALLOCATE ( bhelp )
          IF (.NOT.( iintsp.EQ.2 .AND. jintsp.EQ.2 )) THEN
             iilo = iilo_s
             locol = locol_s             ! restore for other loops
             nkvecprevat = nkvecprevat_s
          ENDIF
       ENDIF ! noco%l_ss
       !-noco
    END IF
    !$OMP END MASTER
  END SUBROUTINE slomat
  !===========================================================================
  REAL FUNCTION legpol(l,arg)
    !
    IMPLICIT NONE
    !     ..
    !     .. Scalar Arguments ..
    REAL arg
    INTEGER l
    !     ..
    !     .. Local Scalars ..
    INTEGER lp
    !     ..
    !     .. Local Arrays ..
    REAL plegend(0:l)
    !     ..
    plegend(0) = 1.0
    IF (l.GE.1) THEN
       plegend(1) = arg
       DO lp = 1,l - 1
          plegend(lp+1) = (lp+lp+1)*arg*plegend(lp)/ (lp+1) -lp*plegend(lp-1)/ (lp+1)
       END DO
    END IF
    legpol = plegend(l)
  END FUNCTION legpol
  !===========================================================================
END MODULE m_slomat
