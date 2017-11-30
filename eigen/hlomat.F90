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
  SUBROUTINE hlomat(input,atoms,usp,tsp,&
       n_size,n_rank, ntyp,na,lapw,ar,br,ai,bi, el,alo,blo,clo,ud,&
       noco,iintsp,jintsp,chi11,chi22,chi21, iilo,locol,nkvecprevat,tlmplm,l_real,aa_r,aa_c)
#include"cpp_double.h"
    !
    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)  :: input
    TYPE(t_noco),INTENT(IN)   :: noco
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_lapw),INTENT(IN)   :: lapw
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: n_size,n_rank
    INTEGER, INTENT (IN) :: na,ntyp  
    INTEGER, INTENT (IN) :: usp,tsp !spin for usdus and tlmplm
    INTEGER, INTENT (IN) :: iintsp,jintsp
    COMPLEX, INTENT (IN) :: chi11,chi22,chi21
    INTEGER, INTENT (INOUT) :: iilo,nkvecprevat,locol
    !     ..
    !     .. Array Arguments ..
    COMPLEX, INTENT (IN) :: alo(-atoms%llod:,:,:,:)
    COMPLEX, INTENT (IN) :: blo(-atoms%llod:,:,:,:)
    COMPLEX, INTENT (IN) :: clo(-atoms%llod:,:,:,:)
    REAL,    INTENT (IN) :: ar(:,0:),br(:,0:) !(nvd,0:lmd)
    REAL,    INTENT (IN) :: ai(:,0:),bi(:,0:)
    REAL,    INTENT (IN) :: el(0:)!(0:atoms%lmax)
    TYPE(t_usdus),INTENT(IN)     ::ud
    TYPE(t_tlmplm),INTENT(INOUT) :: tlmplm

    LOGICAL,INTENT(IN)      :: l_real
    REAL,  ALLOCATABLE,  OPTIONAL,INTENT (INOUT) :: aa_r(:)!(matsize)
    COMPLEX,ALLOCATABLE, OPTIONAL,INTENT (INOUT) :: aa_c(:)
    !     ..
    !     .. Local Scalars ..
    COMPLEX axx,bxx,cxx,dtd,dtu,dtulo,ulotd,ulotu,ulotulo,utd,utu, utulo,chihlp
    INTEGER im,in,invsfct ,l,lm,lmp,lo,lolo,lolop,lop,lp ,matel0 
    INTEGER locol0,mp,nkvec,nkvecp,nkvecprevlo,lmplm , loplo,kp,m
    INTEGER iilo_s,locol_s,nkvecprevat_s,ic,ii,ij,n,k,mlo,mlolo
    !     ..
    !     .. Local Arrays ..
    COMPLEX ax(SIZE(ar,1)),bx(SIZE(ar,1)),cx(SIZE(ar,1))
    COMPLEX, ALLOCATABLE :: ahelp(:)
    INTEGER indt(0:SIZE(ar,2)-1)
    !     ..

    !     .. External Subroutines ..
    EXTERNAL CPP_BLAS_caxpy
    !     ..

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
       nkvecprevlo = 0
       !-noco
       IF (noco%l_ss) THEN
          iilo_s = iilo
          locol_s = locol
          nkvecprevat_s = nkvecprevat ! save for other interstitial spin loops
          ic = 0                                        ! update b-matrix
          DO lo = 1,atoms%nlo(ntyp)
             ic = ic + invsfct* (2*atoms%llo(lo,ntyp)+1)
          ENDDO
          k = ic*(lapw%nv(jintsp)+nkvecprevat) + (ic+1)*ic/2
          ALLOCATE ( ahelp(k) )       ! initialize help-array
          ahelp=CMPLX(0.,0.)
          iilo = 0
       ENDIF
       !-noco
       !
       ! temporarily update the diagonal elements
       !
       DO l = 0,atoms%lnonsph(ntyp)
          DO  m = -l,l
             lm = l* (l+1) + m
             lmplm = (lm* (lm+3))/2
             tlmplm%tuu(lmplm,ntyp,tsp)=tlmplm%tuu(lmplm,ntyp,tsp) + el(l)
             tlmplm%tdd(lmplm,ntyp,tsp)=tlmplm%tdd(lmplm,ntyp,tsp) + el(l)*ud%ddn(l,ntyp,usp)
             tlmplm%tud(lmplm,ntyp,tsp)=tlmplm%tud(lmplm,ntyp,tsp) + 0.5
             tlmplm%tdu(lmplm,ntyp,tsp)=tlmplm%tdu(lmplm,ntyp,tsp) + 0.5
             indt(lm) = tlmplm%ind(lm,lm,ntyp,tsp)
             tlmplm%ind(lm,lm,ntyp,tsp) = lmplm
          ENDDO
       ENDDO
       !
       DO lo = 1,atoms%nlo(ntyp)
          l = atoms%llo(lo,ntyp)
          matel0 = iilo
          locol0 = locol
          !--->       calculate the hamiltonian matrix elements with the regular
          !--->       flapw basis-functions
          DO m = -l,l
             lm = l* (l+1) + m
             iilo = matel0
             locol = locol0
             DO kp = 1,lapw%nv(jintsp)
                ax(kp) = CMPLX(0.0,0.0)
                bx(kp) = CMPLX(0.0,0.0)
                cx(kp) = CMPLX(0.0,0.0)
             END DO
             DO lp = 0,atoms%lnonsph(ntyp)
                DO mp = -lp,lp
                   lmp = lp* (lp+1) + mp
                   in = tlmplm%ind(lmp,lm,ntyp,tsp)
                   IF (in.NE.-9999) THEN
                      IF (in.GE.0) THEN
                         utu = tlmplm%tuu(in,ntyp,tsp)
                         dtu = tlmplm%tdu(in,ntyp,tsp)
                         utd = tlmplm%tud(in,ntyp,tsp)
                         dtd = tlmplm%tdd(in,ntyp,tsp)
                      ELSE
                         im = -in
                         utu = CONJG(tlmplm%tuu(im,ntyp,tsp))
                         dtu = CONJG(tlmplm%tud(im,ntyp,tsp))
                         utd = CONJG(tlmplm%tdu(im,ntyp,tsp))
                         dtd = CONJG(tlmplm%tdd(im,ntyp,tsp))
                      END IF
                      utulo = tlmplm%tuulo(lmp,m,lo+mlo,tsp)
                      dtulo = tlmplm%tdulo(lmp,m,lo+mlo,tsp)
                      !--->                   note, that utu,dtu... are the t-matrices and
                      !--->                   not their complex conjugates as in hssphn
                      !--->                   and that a,b,alo... are the complex
                      !--->                   conjugates of the a,b...-coefficients
                      DO kp = 1,lapw%nv(jintsp)
                         ax(kp) = ax(kp) + &
                              CMPLX(ar(kp,lmp),ai(kp,lmp))*utu +&
                              CMPLX(br(kp,lmp),bi(kp,lmp))*dtu
                         bx(kp) = bx(kp) + &
                              CMPLX(ar(kp,lmp),ai(kp,lmp))*utd +&
                              CMPLX(br(kp,lmp),bi(kp,lmp))*dtd
                         cx(kp) = cx(kp) + &
                              CMPLX(ar(kp,lmp),ai(kp,lmp))*utulo +&
                              CMPLX(br(kp,lmp),bi(kp,lmp))*dtulo
                      END DO
                   END IF
                END DO
             END DO
             !+t3e
             DO nkvec = 1,invsfct* (2*l+1)
                locol = locol + 1
                IF (MOD(locol-1,n_size).EQ.n_rank) THEN
                   !-t3e
                   IF (l_real) THEN
                      DO kp = 1,lapw%nv(jintsp)
                         iilo = iilo + 1
                         aa_r(iilo) = aa_r(iilo) + invsfct * (&
                              REAL(alo(m,nkvec,lo,iintsp))* REAL(ax(kp)) +&
                              AIMAG(alo(m,nkvec,lo,iintsp))*AIMAG(ax(kp)) +&
                              REAL(blo(m,nkvec,lo,iintsp))* REAL(bx(kp)) +&
                              AIMAG(blo(m,nkvec,lo,iintsp))*AIMAG(bx(kp)) +&
                              REAL(clo(m,nkvec,lo,iintsp))* REAL(cx(kp)) +&
                              AIMAG(clo(m,nkvec,lo,iintsp))*AIMAG(cx(kp)) )
                         IF (input%l_useapw) THEN
                            !---> APWlo
                            aa_r(iilo) = aa_r(iilo) + 0.25 * atoms%rmt(ntyp)**2 * invsfct * (&
                                 (ar(kp,lm)*ud%us(l,ntyp,usp)+&
                                 br(kp,lm)*ud%uds(l,ntyp,usp))*&
                                 REAL( alo(m,nkvec,lo,iintsp)*  ud%dus(l,ntyp,usp) +&
                                 blo(m,nkvec,lo,iintsp)* ud%duds(l,ntyp,usp)  +&
                                 clo(m,nkvec,lo,iintsp)*ud%dulos(lo,ntyp,usp) ) +&
                                 (ai(kp,lm)*ud%us(l,ntyp,usp)+bi(kp,lm)*ud%uds(l,ntyp,usp))*&
                                 AIMAG( alo(m,nkvec,lo,iintsp)*  ud%dus(l,ntyp,usp) +&
                                 blo(m,nkvec,lo,iintsp)* ud%duds(l,ntyp,usp) +&
                                 clo(m,nkvec,lo,iintsp)*ud%dulos(lo,ntyp,usp) ) )
                         ENDIF
                      ENDDO
                   ELSE
                      DO kp = 1,lapw%nv(jintsp)
                         iilo = iilo + 1
                         IF (.NOT.noco%l_ss) THEN
                            aa_c(iilo) = aa_c(iilo) + invsfct * (&
                                 alo(m,nkvec,lo,iintsp) * CONJG( ax(kp) ) +&
                                 blo(m,nkvec,lo,iintsp) * CONJG( bx(kp) ) +&
                                 clo(m,nkvec,lo,iintsp) * CONJG( cx(kp) ) )
                            IF (input%l_useapw) THEN
                               !---> APWlo
                               aa_c(iilo)=aa_c(iilo) + 0.25 * atoms%rmt(ntyp)**2 * invsfct*(&
                                    (CMPLX(ar(kp,lm),-ai(kp,lm))* ud%us(l,ntyp,usp)+&
                                    CMPLX(br(kp,lm),-bi(kp,lm))*ud%uds(l,ntyp,usp))*&
                                    (alo(m,nkvec,lo,iintsp)*  ud%dus(l,ntyp,usp)&
                                    +blo(m,nkvec,lo,iintsp)* ud%duds(l,ntyp,usp)&
                                    +clo(m,nkvec,lo,iintsp)*ud%dulos(lo,ntyp,usp) ))
                            ENDIF
                         ELSE
                            ahelp(iilo) = ahelp(iilo) + invsfct * (&
                                 alo(m,nkvec,lo,iintsp) * CONJG( ax(kp) ) +&
                                 blo(m,nkvec,lo,iintsp) * CONJG( bx(kp) ) +&
                                 clo(m,nkvec,lo,iintsp) * CONJG( cx(kp) ) )
                            IF (input%l_useapw) THEN

                               !---> APWlo
                               ahelp(iilo)=ahelp(iilo) + 0.25 * atoms%rmt(ntyp)**2 * invsfct*(&
                                    (CMPLX(ar(kp,lm),-ai(kp,lm))* ud%us(l,ntyp,usp)+&
                                    CMPLX(br(kp,lm),-bi(kp,lm))*ud%uds(l,ntyp,usp))*&
                                    (alo(m,nkvec,lo,iintsp)*  ud%dus(l,ntyp,usp)&
                                    +blo(m,nkvec,lo,iintsp)* ud%duds(l,ntyp,usp)&
                                    +clo(m,nkvec,lo,iintsp)*ud%dulos(lo,ntyp,usp) ))
                            ENDIF
                         ENDIF
                      ENDDO
                   ENDIF
                   !--->             jump to the last matrixelement of the current row
                   iilo = iilo + nkvecprevat + nkvecprevlo + nkvec
                ENDIF
             END DO
          END DO
          !--->       calculate the hamiltonian matrix elements with other
          !--->       local orbitals at the same atom and with itself
          iilo = matel0
          !+t3e
          locol = locol0
          DO nkvec = 1,invsfct* (2*l+1)
             locol = locol + 1
             IF (MOD(locol-1,n_size).EQ.n_rank) THEN
                !-t3e
                !--->          skip the matrixelements with regular flapw-fcn. and
                !--->          with local orbitals at other atoms
                iilo = iilo + lapw%nv(jintsp) + nkvecprevat
                !--->          calculate the hamiltonian matrix elements with other
                !--->          local orbitals at the same atom, if they have the same l
                DO lop = 1, (lo-1)
                   lp = atoms%llo(lop,ntyp)
                   DO nkvecp = 1,invsfct* (2*lp+1)
                      iilo = iilo + 1
                      DO m = -l,l
                         lm = l* (l+1) + m
                         DO mp = -lp,lp
                            lmp = lp* (lp+1) + mp
                            in = tlmplm%ind(lmp,lm,ntyp,tsp)
                            IF (in.NE.-9999) THEN
                               IF (in.GE.0) THEN
                                  utu = tlmplm%tuu(in,ntyp,tsp)
                                  dtu = tlmplm%tdu(in,ntyp,tsp)
                                  utd = tlmplm%tud(in,ntyp,tsp)
                                  dtd = tlmplm%tdd(in,ntyp,tsp)
                               ELSE
                                  im = -in
                                  utu = CONJG(tlmplm%tuu(im,ntyp,tsp))
                                  dtu = CONJG(tlmplm%tud(im,ntyp,tsp))
                                  utd = CONJG(tlmplm%tdu(im,ntyp,tsp))
                                  dtd = CONJG(tlmplm%tdd(im,ntyp,tsp))
                               END IF
                               utulo = tlmplm%tuulo(lmp,m,lo+mlo,tsp)
                               dtulo = tlmplm%tdulo(lmp,m,lo+mlo,tsp)
                               ulotu=CONJG(tlmplm%tuulo(lm,mp,lop+mlo,tsp))
                               ulotd=CONJG(tlmplm%tdulo(lm,mp,lop+mlo,tsp))
                               !--->                         note that lo > lop
                               lolop = ((lo-1)*lo)/2 + lop
                               ulotulo = CONJG(tlmplm%tuloulo (m,mp,lolop+mlolo,tsp))
                               axx=CONJG(alo(m,nkvec,lo,iintsp))*utu +&
                                    CONJG(blo(m,nkvec,lo,iintsp))*utd +&
                                    CONJG(clo(m,nkvec,lo,iintsp))*utulo
                               bxx=CONJG(alo(m,nkvec,lo,iintsp))*dtu +&
                                    CONJG(blo(m,nkvec,lo,iintsp))*dtd +&
                                    CONJG(clo(m,nkvec,lo,iintsp))*dtulo
                               cxx = &
                                    CONJG(alo(m,nkvec,lo,iintsp))*ulotu +&
                                    CONJG(blo(m,nkvec,lo,iintsp))*ulotd +&
                                    CONJG(clo(m,nkvec,lo,iintsp))*ulotulo
                               IF (l_real) THEN
                                  aa_r(iilo) = aa_r(iilo) + invsfct * (&
                                       REAL(alo(mp,nkvecp,lop,jintsp))* REAL(axx) -&
                                       AIMAG(alo(mp,nkvecp,lop,jintsp))*AIMAG(axx) +&
                                       REAL(blo(mp,nkvecp,lop,jintsp))* REAL(bxx) -&
                                       AIMAG(blo(mp,nkvecp,lop,jintsp))*AIMAG(bxx) +&
                                       REAL(clo(mp,nkvecp,lop,jintsp))* REAL(cxx) -&
                                       AIMAG(clo(mp,nkvecp,lop,jintsp))*AIMAG(cxx) )
                               ELSE
                                  IF (.NOT.noco%l_ss) THEN
                                     aa_c(iilo) = aa_c(iilo) + invsfct * CONJG(&
                                          alo(mp,nkvecp,lop,jintsp) * axx +&
                                          blo(mp,nkvecp,lop,jintsp) * bxx +&
                                          clo(mp,nkvecp,lop,jintsp) * cxx )
                                  ELSE
                                     ahelp(iilo)=ahelp(iilo)+invsfct*CONJG(&
                                          alo(mp,nkvecp,lop,jintsp) * axx +&
                                          blo(mp,nkvecp,lop,jintsp) * bxx +&
                                          clo(mp,nkvecp,lop,jintsp) * cxx )
                                  ENDIF
                               ENDIF
                            END IF
                         END DO
                      END DO
                   END DO
                END DO
                !--->          calculate the hamiltonian matrix elements of one local
                !--->          orbital with itself
                DO nkvecp = 1,nkvec
                   iilo = iilo + 1
                   DO m = -l,l
                      lm = l* (l+1) + m
                      DO mp = -l,l
                         lmp = l* (l+1) + mp
                         in = tlmplm%ind(lmp,lm,ntyp,tsp)
                         IF (in.NE.-9999) THEN
                            IF (in.GE.0) THEN
                               utu = tlmplm%tuu(in,ntyp,tsp)
                               dtu = tlmplm%tdu(in,ntyp,tsp)
                               utd = tlmplm%tud(in,ntyp,tsp)
                               dtd = tlmplm%tdd(in,ntyp,tsp)
                            ELSE
                               im = -in
                               utu = CONJG(tlmplm%tuu(im,ntyp,tsp))
                               dtu = CONJG(tlmplm%tud(im,ntyp,tsp))
                               utd = CONJG(tlmplm%tdu(im,ntyp,tsp))
                               dtd = CONJG(tlmplm%tdd(im,ntyp,tsp))
                            END IF
                            utulo = tlmplm%tuulo(lmp,m,lo+mlo,tsp)
                            dtulo = tlmplm%tdulo(lmp,m,lo+mlo,tsp)
                            ulotu = CONJG(tlmplm%tuulo(lm,mp,lo+mlo,tsp))
                            ulotd = CONJG(tlmplm%tdulo(lm,mp,lo+mlo,tsp))
                            lolo = ((lo-1)*lo)/2 + lo
                            ulotulo =CONJG(tlmplm%tuloulo(m,mp,lolo+mlolo,tsp))
                            axx = CONJG(alo(m,nkvec,lo,iintsp))*utu +&
                                 CONJG(blo(m,nkvec,lo,iintsp))*utd +&
                                 CONJG(clo(m,nkvec,lo,iintsp))*utulo
                            bxx = CONJG(alo(m,nkvec,lo,iintsp))*dtu +&
                                 CONJG(blo(m,nkvec,lo,iintsp))*dtd +&
                                 CONJG(clo(m,nkvec,lo,iintsp))*dtulo
                            cxx = CONJG(alo(m,nkvec,lo,iintsp))*ulotu +&
                                 CONJG(blo(m,nkvec,lo,iintsp))*ulotd +&
                                 CONJG(clo(m,nkvec,lo,iintsp))*ulotulo
                            IF (l_real) THEN
                               aa_r(iilo) = aa_r(iilo) + invsfct* (&
                                    REAL(alo(mp,nkvecp,lo,jintsp))* REAL(axx) -&
                                    AIMAG(alo(mp,nkvecp,lo,jintsp))*AIMAG(axx) +&
                                    REAL(blo(mp,nkvecp,lo,jintsp))* REAL(bxx) -&
                                    AIMAG(blo(mp,nkvecp,lo,jintsp))*AIMAG(bxx) +&
                                    REAL(clo(mp,nkvecp,lo,jintsp))* REAL(cxx) -&
                                    AIMAG(clo(mp,nkvecp,lo,jintsp))*AIMAG(cxx) )
                            ELSE
                               IF (.NOT.noco%l_ss) THEN
                                  aa_c(iilo) = aa_c(iilo) + invsfct* CONJG(&
                                       alo(mp,nkvecp,lo,jintsp)*axx +&
                                       blo(mp,nkvecp,lo,jintsp)*bxx +&
                                       clo(mp,nkvecp,lo,jintsp)*cxx )
                               ELSE
                                  ahelp(iilo) = ahelp(iilo) + invsfct* CONJG(&
                                       alo(mp,nkvecp,lo,jintsp)*axx +&
                                       blo(mp,nkvecp,lo,jintsp)*bxx +&
                                       clo(mp,nkvecp,lo,jintsp)*cxx )
                               ENDIF
                            ENDIF
                         END IF
                      END DO
                   END DO
                END DO
             ENDIF
             !-t3e
          END DO
          nkvecprevlo = nkvecprevlo + invsfct* (2*l+1)
       END DO ! end of lo = 1,atoms%nlo loop
       !
       ! remove the temporary update of the diagonal elements
       !
       DO l = 0,atoms%lnonsph(ntyp)
          DO  m = -l,l
             lm = l* (l+1) + m
             lmplm = (lm* (lm+3))/2
             tlmplm%tuu(lmplm,ntyp,tsp)=tlmplm%tuu(lmplm,ntyp,tsp) - el(l)
             tlmplm%tdd(lmplm,ntyp,tsp)=tlmplm%tdd(lmplm,ntyp,tsp) - el(l)*ud%ddn(l,ntyp,usp)
             tlmplm%tud(lmplm,ntyp,tsp)=tlmplm%tud(lmplm,ntyp,tsp) - 0.5
             tlmplm%tdu(lmplm,ntyp,tsp)=tlmplm%tdu(lmplm,ntyp,tsp) - 0.5
             tlmplm%ind(lm,lm,ntyp,tsp) = indt(lm)
          ENDDO
       ENDDO
       !
       nkvecprevat = nkvecprevat + nkvecprevlo
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
                IF (.NOT.l_real) THEN
                   IF (iintsp.EQ.2 .AND. jintsp.EQ.1 ) THEN
                      CALL CPP_BLAS_caxpy(n-1,chihlp,ahelp(ii+1),1,aa_c(ij+1),1)
                   ELSE
                      CALL CPP_BLAS_caxpy(n,chihlp,ahelp(ii+1),1,aa_c(ij+1),1)
                   ENDIF
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
                   aa_c(ij) = aa_c(ij) + chihlp * CONJG( ahelp(ii) )
                   ij = ij + lapw%nv(1) + kp + atoms%nlotot
                ENDDO
             ENDDO
          ENDIF
          DEALLOCATE ( ahelp )
          IF (.NOT.( iintsp.EQ.2 .AND. jintsp.EQ.2 )) THEN
             iilo = iilo_s
             locol = locol_s             ! restore for other loops
             nkvecprevat = nkvecprevat_s
          ENDIF
       ELSE
          k = lapw%nv(1) * ( lapw%nv(1) + 1 ) / 2
       ENDIF
       !-noco
    END IF
    !$OMP END MASTER
    !$OMP barrier
  END SUBROUTINE hlomat
      END MODULE m_hlomat
