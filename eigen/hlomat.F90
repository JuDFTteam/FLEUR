!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hlomat
!***********************************************************************
! updates the hamiltonian  matrix with the contributions from the local
! orbitals.
! p.kurz sept. 1996
!***********************************************************************
CONTAINS
  SUBROUTINE hlomat(input,atoms,mpi,lapw,ud,tlmplm,sym,cell,noco,isp,&
       ntyp,na,fj,gj,alo1,blo1,clo1, iintsp,jintsp,chi,hmat)
    !
    USE m_hsmt_ab
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN)   :: lapw
    TYPE(t_mpi),INTENT(IN)    :: mpi
    TYPE(t_usdus),INTENT(IN)  :: ud
    TYPE(t_tlmplm),INTENT(IN) :: tlmplm
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_cell),INTENT(IN)   :: cell
    TYPE(t_noco),INTENT(IN)   :: noco
    
    
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: na,ntyp  
    INTEGER, INTENT (IN) :: isp !spin for usdus and tlmplm
    INTEGER, INTENT (IN) :: iintsp,jintsp
    COMPLEX, INTENT (IN) :: chi
    !     ..
    !     .. Array Arguments ..
    REAL, INTENT (IN) :: alo1(:),blo1(:),clo1(:)
    REAL,INTENT(IN)      :: fj(:,0:,:),gj(:,0:,:)

    TYPE(t_lapwmat),INTENT (INOUT) :: hmat
    !     ..
    !     .. Local Scalars ..
    COMPLEX axx,bxx,cxx,dtd,dtu,dtulo,ulotd,ulotu,ulotulo,utd,utu, utulo
    INTEGER im,in,invsfct,l,lm,lmp,lo,lolo,lolop,lop,lp,i  
    INTEGER mp,nkvec,nkvecp,lmplm,loplo,kp,m,mlo,mlolo
    INTEGER locol,lorow,ii,ij,n,k,ab_size
    !     ..
    !     .. Local Arrays ..
    COMPLEX, ALLOCATABLE :: ab(:,:),ax(:),bx(:),cx(:)
    COMPLEX,ALLOCATABLE  :: abclo(:,:,:,:,:)
    !     ..

    
    !-->              synthesize the complex conjugates of a and b
    ALLOCATE(ab(MAXVAL(lapw%nv),0:2*atoms%lmaxd*(atoms%lmaxd+2)+1))
    ALLOCATE(ax(MAXVAL(lapw%nv)),bx(MAXVAL(lapw%nv)),cx(MAXVAL(lapw%nv)))
    ALLOCATE(abclo(3,-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod,2))
    DO i=MIN(iintsp,jintsp),MAX(iintsp,jintsp)
       CALL hsmt_ab(sym,atoms,noco,isp,iintsp,ntyp,na,cell,lapw,fj,gj,ab(:,:),ab_size,.TRUE.,abclo(:,:,:,:,i),alo1,blo1,clo1)
    ENDDO

    
    mlo=0;mlolo=0
    DO m=1,ntyp-1
       mlo=mlo+atoms%nlo(m)
       mlolo=mlolo+atoms%nlo(m)*(atoms%nlo(m)+1)/2
    ENDDO


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
       !
    
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          !--->       calculate the hamiltonian matrix elements with the regular
          !--->       flapw basis-functions
          DO m = -l,l
             lm = l* (l+1) + m
             DO kp = 1,lapw%nv(jintsp)
                ax(kp) = CMPLX(0.0,0.0)
                bx(kp) = CMPLX(0.0,0.0)
                cx(kp) = CMPLX(0.0,0.0)
             END DO
             DO lp = 0,atoms%lnonsph(ntyp)
                DO mp = -lp,lp
                   lmp = lp* (lp+1) + mp
                   in = tlmplm%ind(lmp,lm,ntyp,isp)
                   IF (in.NE.-9999) THEN
                      IF (in.GE.0) THEN
                         utu = tlmplm%tuu(in,ntyp,isp)
                         dtu = tlmplm%tdu(in,ntyp,isp)
                         utd = tlmplm%tud(in,ntyp,isp)
                         dtd = tlmplm%tdd(in,ntyp,isp)
                      ELSE
                         im = -in
                         utu = CONJG(tlmplm%tuu(im,ntyp,isp))
                         dtu = CONJG(tlmplm%tud(im,ntyp,isp))
                         utd = CONJG(tlmplm%tdu(im,ntyp,isp))
                         dtd = CONJG(tlmplm%tdd(im,ntyp,isp))
                      END IF
                      utulo = tlmplm%tuulo(lmp,m,lo+mlo,isp)
                      dtulo = tlmplm%tdulo(lmp,m,lo+mlo,isp)
                      !--->                   note, that utu,dtu... are the t-matrices and
                      !--->                   not their complex conjugates as in hssphn
                      !--->                   and that a,b,alo... are the complex
                      !--->                   conjugates of the a,b...-coefficients
                      DO kp = 1,lapw%nv(jintsp)
                         ax(kp) = ax(kp) + ab(kp,lmp)*utu + ab(kp,ab_size/2+lmp)*dtu
                         bx(kp) = bx(kp) + ab(kp,lmp)*utd + ab(kp,ab_size/2+lmp)*dtd
                         cx(kp) = cx(kp) + ab(kp,lmp)*utulo + ab(kp,ab_size/2+lmp)*dtulo
                      END DO
                   END IF
                END DO
             END DO
             !+t3e
             DO nkvec = 1,invsfct* (2*l+1)
               locol= lapw%nv(iintsp)+lapw%index_lo(lo,ntyp)+nkvec !this is the column of the matrix
                IF (MOD(locol-1,mpi%n_size).EQ.mpi%n_rank) THEN
                   locol=(locol-1)/mpi%n_size+1 !this is the column in local storage
                   !-t3e
                   IF (hmat%l_real) THEN
                      DO kp = 1,lapw%nv(jintsp)
                         hmat%data_r(kp,locol) = hmat%data_r(kp,locol) + chi*invsfct * (&
                              REAL(abclo(1,m,nkvec,lo,iintsp))* REAL(ax(kp)) +&
                              AIMAG(abclo(1,m,nkvec,lo,iintsp))*AIMAG(ax(kp)) +&
                              REAL(abclo(2,m,nkvec,lo,iintsp))* REAL(bx(kp)) +&
                              AIMAG(abclo(2,m,nkvec,lo,iintsp))*AIMAG(bx(kp)) +&
                              REAL(abclo(3,m,nkvec,lo,iintsp))* REAL(cx(kp)) +&
                              AIMAG(abclo(3,m,nkvec,lo,iintsp))*AIMAG(cx(kp)) )
                         IF (input%l_useapw) THEN
                            !---> APWlo
                            hmat%data_r(kp,locol) = hmat%data_r(kp,locol) + 0.25 * atoms%rmt(ntyp)**2 * chi*invsfct * (&
                                 (CONJG(ab(kp,lm))* ud%us(l,ntyp,isp)+&
                                 CONJG(ab(kp,ab_size/2+lm))*ud%uds(l,ntyp,isp))*&
                                 (abclo(1,m,nkvec,lo,iintsp)*  ud%dus(l,ntyp,isp)&
                                 +abclo(2,m,nkvec,lo,iintsp)* ud%duds(l,ntyp,isp)&
                                 +abclo(3,m,nkvec,lo,iintsp)*ud%dulos(lo,ntyp,isp) ))
                         ENDIF
                      ENDDO
                   ELSE
                      DO kp = 1,lapw%nv(jintsp)
                            hmat%data_c(kp,locol) = hmat%data_c(kp,locol) + chi*invsfct * (&
                                 abclo(1,m,nkvec,lo,iintsp) * CONJG( ax(kp) ) +&
                                 abclo(2,m,nkvec,lo,iintsp) * CONJG( bx(kp) ) +&
                                 abclo(3,m,nkvec,lo,iintsp) * CONJG( cx(kp) ) )
                            IF (input%l_useapw) THEN
                               !---> APWlo
                               hmat%data_c(kp,locol)=hmat%data_c(kp,locol) + 0.25 * atoms%rmt(ntyp)**2 * chi*invsfct*(&
                                    (CONJG(ab(kp,lm))* ud%us(l,ntyp,isp)+&
                                    CONJG(ab(kp,ab_size/2+lm))*ud%uds(l,ntyp,isp))*&
                                    (abclo(1,m,nkvec,lo,iintsp)*  ud%dus(l,ntyp,isp)&
                                    +abclo(2,m,nkvec,lo,iintsp)* ud%duds(l,ntyp,isp)&
                                    +abclo(3,m,nkvec,lo,iintsp)*ud%dulos(lo,ntyp,isp) ))
                            ENDIF
                      ENDDO
                   ENDIF
                   !--->             jump to the last matrixelement of the current row
                ENDIF
             END DO
          END DO
          !--->       calculate the hamiltonian matrix elements with other
          !--->       local orbitals at the same atom and with itself
          DO nkvec = 1,invsfct* (2*l+1)
             locol = lapw%nv(iintsp)+lapw%index_lo(lo,ntyp)+nkvec !this is the column of the matrix
             IF (MOD(locol-1,mpi%n_size).EQ.mpi%n_rank) THEN
                locol=(locol-1)/mpi%n_size+1 !this is the column in local storage
                !-t3e
                !--->          calculate the hamiltonian matrix elements with other
                !--->          local orbitals at the same atom, if they have the same l
                DO lop = 1, (lo-1)
                   lp = atoms%llo(lop,ntyp)
                   DO nkvecp = 1,invsfct* (2*lp+1)
                      lorow=lapw%nv(jintsp)+lapw%index_lo(lop,ntyp)+nkvecp
                      DO m = -l,l
                         lm = l* (l+1) + m
                         DO mp = -lp,lp
                            lmp = lp* (lp+1) + mp
                            in = tlmplm%ind(lmp,lm,ntyp,isp)
                            IF (in.NE.-9999) THEN
                               IF (in.GE.0) THEN
                                  utu = tlmplm%tuu(in,ntyp,isp)
                                  dtu = tlmplm%tdu(in,ntyp,isp)
                                  utd = tlmplm%tud(in,ntyp,isp)
                                  dtd = tlmplm%tdd(in,ntyp,isp)
                               ELSE
                                  im = -in
                                  utu = CONJG(tlmplm%tuu(im,ntyp,isp))
                                  dtu = CONJG(tlmplm%tud(im,ntyp,isp))
                                  utd = CONJG(tlmplm%tdu(im,ntyp,isp))
                                  dtd = CONJG(tlmplm%tdd(im,ntyp,isp))
                               END IF
                               utulo = tlmplm%tuulo(lmp,m,lo+mlo,isp)
                               dtulo = tlmplm%tdulo(lmp,m,lo+mlo,isp)
                               ulotu=CONJG(tlmplm%tuulo(lm,mp,lop+mlo,isp))
                               ulotd=CONJG(tlmplm%tdulo(lm,mp,lop+mlo,isp))
                               !--->                         note that lo > lop
                               lolop = ((lo-1)*lo)/2 + lop
                               ulotulo = CONJG(tlmplm%tuloulo (m,mp,lolop+mlolo,isp))
                               axx=CONJG(abclo(1,m,nkvec,lo,iintsp))*utu +&
                                    CONJG(abclo(2,m,nkvec,lo,iintsp))*utd +&
                                    CONJG(abclo(3,m,nkvec,lo,iintsp))*utulo
                               bxx=CONJG(abclo(1,m,nkvec,lo,iintsp))*dtu +&
                                    CONJG(abclo(2,m,nkvec,lo,iintsp))*dtd +&
                                    CONJG(abclo(3,m,nkvec,lo,iintsp))*dtulo
                               cxx = &
                                    CONJG(abclo(1,m,nkvec,lo,iintsp))*ulotu +&
                                    CONJG(abclo(2,m,nkvec,lo,iintsp))*ulotd +&
                                    CONJG(abclo(3,m,nkvec,lo,iintsp))*ulotulo
                               IF (hmat%l_real) THEN
                                  hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) + chi*invsfct * (&
                                       REAL(abclo(1,mp,nkvecp,lop,jintsp))* REAL(axx) -&
                                       AIMAG(abclo(1,mp,nkvecp,lop,jintsp))*AIMAG(axx) +&
                                       REAL(abclo(2,mp,nkvecp,lop,jintsp))* REAL(bxx) -&
                                       AIMAG(abclo(2,mp,nkvecp,lop,jintsp))*AIMAG(bxx) +&
                                       REAL(abclo(3,mp,nkvecp,lop,jintsp))* REAL(cxx) -&
                                       AIMAG(abclo(3,mp,nkvecp,lop,jintsp))*AIMAG(cxx) )
                               ELSE
                                  hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) + chi*invsfct * CONJG(&
                                       abclo(1,mp,nkvecp,lop,jintsp) * axx +&
                                       abclo(2,mp,nkvecp,lop,jintsp) * bxx +&
                                       abclo(3,mp,nkvecp,lop,jintsp) * cxx )
                               ENDIF
                            END IF
                         END DO
                      END DO
                   END DO
                END DO
                !--->          calculate the hamiltonian matrix elements of one local
                !--->          orbital with itself
                DO nkvecp = 1,nkvec
                   lorow=lapw%nv(jintsp)+lapw%index_lo(lop,ntyp)+nkvecp
                   DO m = -l,l
                      lm = l* (l+1) + m
                      DO mp = -l,l
                         lmp = l* (l+1) + mp
                         in = tlmplm%ind(lmp,lm,ntyp,isp)
                         IF (in.NE.-9999) THEN
                            IF (in.GE.0) THEN
                               utu = tlmplm%tuu(in,ntyp,isp)
                               dtu = tlmplm%tdu(in,ntyp,isp)
                               utd = tlmplm%tud(in,ntyp,isp)
                               dtd = tlmplm%tdd(in,ntyp,isp)
                            ELSE
                               im = -in
                               utu = CONJG(tlmplm%tuu(im,ntyp,isp))
                               dtu = CONJG(tlmplm%tud(im,ntyp,isp))
                               utd = CONJG(tlmplm%tdu(im,ntyp,isp))
                               dtd = CONJG(tlmplm%tdd(im,ntyp,isp))
                            END IF
                            utulo = tlmplm%tuulo(lmp,m,lo+mlo,isp)
                            dtulo = tlmplm%tdulo(lmp,m,lo+mlo,isp)
                            ulotu = CONJG(tlmplm%tuulo(lm,mp,lo+mlo,isp))
                            ulotd = CONJG(tlmplm%tdulo(lm,mp,lo+mlo,isp))
                            lolo = ((lo-1)*lo)/2 + lo
                            ulotulo =CONJG(tlmplm%tuloulo(m,mp,lolo+mlolo,isp))
                            axx = CONJG(abclo(1,m,nkvec,lo,iintsp))*utu +&
                                 CONJG(abclo(2,m,nkvec,lo,iintsp))*utd +&
                                 CONJG(abclo(3,m,nkvec,lo,iintsp))*utulo
                            bxx = CONJG(abclo(1,m,nkvec,lo,iintsp))*dtu +&
                                 CONJG(abclo(2,m,nkvec,lo,iintsp))*dtd +&
                                 CONJG(abclo(3,m,nkvec,lo,iintsp))*dtulo
                            cxx = CONJG(abclo(1,m,nkvec,lo,iintsp))*ulotu +&
                                 CONJG(abclo(2,m,nkvec,lo,iintsp))*ulotd +&
                                 CONJG(abclo(3,m,nkvec,lo,iintsp))*ulotulo
                            IF (hmat%l_real) THEN
                               hmat%data_r(lorow,locol) = hmat%data_r(lorow,locol) + chi*invsfct* (&
                                    REAL(abclo(1,mp,nkvecp,lo,jintsp))* REAL(axx) -&
                                    AIMAG(abclo(1,mp,nkvecp,lo,jintsp))*AIMAG(axx) +&
                                    REAL(abclo(2,mp,nkvecp,lo,jintsp))* REAL(bxx) -&
                                    AIMAG(abclo(2,mp,nkvecp,lo,jintsp))*AIMAG(bxx) +&
                                    REAL(abclo(3,mp,nkvecp,lo,jintsp))* REAL(cxx) -&
                                    AIMAG(abclo(3,mp,nkvecp,lo,jintsp))*AIMAG(cxx) )
                            ELSE
                               hmat%data_c(lorow,locol) = hmat%data_c(lorow,locol) + chi*invsfct* CONJG(&
                                    abclo(1,mp,nkvecp,lo,jintsp)*axx +&
                                    abclo(2,mp,nkvecp,lo,jintsp)*bxx +&
                                    abclo(3,mp,nkvecp,lo,jintsp)*cxx )
                            ENDIF
                         END IF
                      END DO
                   END DO
                END DO
             ENDIF
             !-t3e
          END DO
       END DO ! end of lo = 1,atoms%nlo loop
   
    END IF
    !$OMP END MASTER
    !$OMP barrier
  END SUBROUTINE hlomat
      END MODULE m_hlomat
