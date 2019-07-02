!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsmt_extra
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE hsmt_extra(DIMENSION,atoms,sym,isp,n_size,n_rank,input,nintsp,sub_comm,&
       hlpmsize,lmaxb,noco,l_socfirst, lapw,cell,el, fj,gj,gk,vk,tlmplm,usdus, vs_mmp,oneD,& !in
       kveclo,l_real,aa_r,bb_r,aa_c,bb_c) !out/inout
    USE m_constants, ONLY : tpi_const,fpi_const
    USE m_uham
    USE m_ylm
    USE m_abccoflo
    USE m_hlomat
    USE m_slomat
    USE m_vecforlo
    USE m_setabc1lo
    USE m_hsmt_spinor
    USE m_hsmt_hlptomat
    USE m_types
    IMPLICIT NONE
    TYPE(t_dimension),INTENT(IN):: DIMENSION
    TYPE(t_oneD),INTENT(IN)     :: oneD
    TYPE(t_input),INTENT(IN)    :: input
    TYPE(t_noco),INTENT(IN)     :: noco
    TYPE(t_sym),INTENT(IN)      :: sym
    TYPE(t_cell),INTENT(IN)     :: cell
    TYPE(t_atoms),INTENT(IN)    :: atoms
    TYPE(t_lapw),INTENT(INOUT)  :: lapw !lapw%nmat is updated with lo-basisfcts
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp
    INTEGER, INTENT (IN) :: n_size,n_rank ,nintsp,sub_comm
    INTEGER, INTENT (IN) :: hlpmsize
    INTEGER, INTENT (IN) :: lmaxb
    LOGICAL, INTENT (IN) :: l_socfirst 


    !     ..
    !     .. Array Arguments ..

    REAL,    INTENT (IN) :: el(0:atoms%lmaxd,atoms%ntype,DIMENSION%jspd)
    COMPLEX,INTENT(IN)   :: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins)
    REAL,INTENT(IN),ALLOCATABLE :: fj(:,:,:,:),gj(:,:,:,:),gk(:,:,:),vk(:,:,:)


    TYPE(t_usdus),INTENT(IN)   :: usdus

    INTEGER, INTENT (OUT)    :: kveclo(atoms%nlotot)
    TYPE(t_tlmplm),INTENT(INOUT)::tlmplm


    REAL,    OPTIONAL,ALLOCATABLE,INTENT (INOUT) :: aa_r(:),bb_r(:)!(matsize)
    COMPLEX, OPTIONAL,ALLOCATABLE,INTENT (INOUT) :: aa_c(:),bb_c(:)!(matsize)
    LOGICAL,INTENT(IN)      :: l_real
    !     ..
    !     .. Local Scalars ..
    REAL con1,termi,termr,th,invsfct


    COMPLEX chi11,chi21,chi22
    INTEGER k,i,spin2,l,ll1,lo,jd
    INTEGER m,n,na,nn,np,i_u
    INTEGER iiloh,iilos,nkvecprevath,nkvecprevats,iintsp,jintsp
    INTEGER nc,locolh,locols,nkvecprevatu,iilou,locolu
    INTEGER nkvecprevatuTemp,iilouTemp,locoluTemp
    INTEGER ab_dim,nkvec_sv,fjstart
    LOGICAL enough,l_lo1
    !     ..
    !     .. Local Arrays ..
    INTEGER kvec(2*(2*atoms%llod+1),atoms%nlod )
    REAL alo1(atoms%nlod),blo1(atoms%nlod),clo1(atoms%nlod)
    REAL bmrot(3,3),gkrot(DIMENSION%nvd,3),vmult(3),v(3)
    REAL qssbti(3)
    REAL, ALLOCATABLE :: ar(:,:,:),ai(:,:,:),br(:,:,:),bi(:,:,:)
    REAL, ALLOCATABLE :: rph(:,:),cph(:,:)
    COMPLEX, ALLOCATABLE :: alo(:,:,:,:),blo(:,:,:,:),clo(:,:,:,:)
    INTEGER, ALLOCATABLE :: nkvec(:,:)
    COMPLEX,ALLOCATABLE :: aahlp(:),bbhlp(:)

    COMPLEX ylm( (atoms%lmaxd+1)**2 )

    COMPLEX chi(2,2)


    REAL, PARAMETER :: eps = 1.0E-30

    con1=fpi_const/SQRT(cell%omtil)
    ab_dim=1
    IF ( noco%l_noco .AND. (.NOT. noco%l_ss) )  ALLOCATE ( aahlp(hlpmsize),bbhlp(hlpmsize) )
    IF (noco%l_ss) ab_dim=2
    ALLOCATE(ar(DIMENSION%nvd,0:DIMENSION%lmd,ab_dim),ai(DIMENSION%nvd,0:DIMENSION%lmd,ab_dim))
    ALLOCATE(br(DIMENSION%nvd,0:DIMENSION%lmd,ab_dim),bi(DIMENSION%nvd,0:DIMENSION%lmd,ab_dim))
    ALLOCATE(alo(-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod,ab_dim))
    ALLOCATE(blo(-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod,ab_dim))
    ALLOCATE(clo(-atoms%llod:atoms%llod,2*(2*atoms%llod+1),atoms%nlod,ab_dim))
    ALLOCATE(rph(DIMENSION%nvd,ab_dim),cph(DIMENSION%nvd,ab_dim))
    ALLOCATE(nkvec(atoms%nlod,ab_dim))
    na = 0
    nkvecprevats = 0
    nkvecprevath = 0
    nkvecprevatu = 0
    nkvec_sv = 0
    !Determine index of first LO
    locols = lapw%nv(1)
    locolh = lapw%nv(1)
#ifdef CPP_MPI
    nc = 0
    k = 0
    DO  i=1+n_rank, lapw%nv(1), n_size
       nc = nc + 1
       k = k + n_size*(nc-1) + n_rank + 1
    ENDDO
    iilos = k
    iiloh = k
#else
    iilos = lapw%nv(1)* (lapw%nv(1)+1)/2 
    iiloh = lapw%nv(1)* (lapw%nv(1)+1)/2
#endif

    iilou = iilos
    locolu = locols

    i_u = 1
    ntype_loop: DO n=1,atoms%ntype

       IF (noco%l_noco) THEN
          IF (.NOT.noco%l_ss) aahlp=CMPLX(0.,0.)
          IF (.NOT.noco%l_ss) bbhlp=CMPLX(0.,0.)
          CALL hsmt_spinor(isp,n,noco,input,chi,chi11,chi21,chi22)
       ENDIF
       DO nn = 1,atoms%neq(n)
          na = na + 1
          IF (atoms%lnonsph(n).LT.0) CYCLE ntype_loop
          IF ((sym%invsat(na).EQ.0) .OR. (sym%invsat(na).EQ.1)) THEN
             IF (sym%invsat(na).EQ.0) invsfct = 1
             IF (sym%invsat(na).EQ.1) invsfct = 2
             np = sym%invtab(sym%ngopr(na))
             IF (oneD%odi%d1) THEN
                np = oneD%ods%ngopr(na)
             END IF
             CALL vec_for_lo(atoms,nintsp,sym,na, n,np,noco, lapw,cell, gk,vk, nkvec,kvec)
             DO lo = 1,atoms%nlo(n)
                kveclo(nkvec_sv+1:nkvec_sv+nkvec(lo,1)) = kvec(1:nkvec(lo,1),lo)
                nkvec_sv = nkvec_sv+nkvec(lo,1)
                nkvec(lo,:) = 0
             ENDDO

             !--->       loop over interstitial spins
             DO iintsp = 1,nintsp
                IF (noco%l_constr.OR.l_socfirst) THEN
                   spin2=isp
                ELSE
                   spin2=iintsp
                ENDIF
                IF (iintsp==1) THEN
                   qssbti(:) = - noco%qss(:)/2
                ELSE
                   qssbti(:) = + noco%qss(:)/2
                ENDIF
                !--->          set up phase factors
                DO k = 1,lapw%nv(iintsp)
                   th= DOT_PRODUCT((/lapw%k1(k,iintsp),lapw%k2(k,iintsp),lapw%k3(k,iintsp)/)+qssbti,atoms%taual(:,na))
                   rph(k,iintsp) = COS(tpi_const*th)
                   cph(k,iintsp) = -SIN(tpi_const*th)
                END DO

                !--->          set up the a,b and c  coefficients
                !--->          for the local orbitals, if necessary.
                IF (atoms%nlo(n).GE.1) THEN
                   enough = .FALSE.
                   CALL setabc1lo(atoms, n,usdus,isp, alo1,blo1,clo1) 
                ELSE
                   enough = .TRUE.
                END IF

                !--->          set up the a and b coefficients

                IF (np==1) THEN
                   gkrot( 1:lapw%nv(iintsp),:) = gk( 1:lapw%nv(iintsp),:,iintsp)
                ELSE 
                   IF (oneD%odi%d1) THEN
                      bmrot=MATMUL(oneD%ods%mrot(:,:,np),cell%bmat)
                   ELSE
                      bmrot=MATMUL(1.*sym%mrot(:,:,np),cell%bmat)
                   END IF
                   DO k = 1,lapw%nv(iintsp)
                      !-->                 apply the rotation that brings this atom into the
                      !-->                 representative (this is the definition of ngopr(na)
                      !-->                 and transform to cartesian coordinates
                      v(:) = vk(k,:,iintsp)
                      gkrot(k,:) = MATMUL(TRANSPOSE(bmrot),v)
                   END DO
                END IF
                DO k = 1,lapw%nv(iintsp)
                   !-->    generate spherical harmonics
                   vmult(:) =  gkrot(k,:)
                   CALL ylm4(atoms%lnonsph(n),vmult,ylm)
                   IF (.NOT.enough) THEN
                      l_lo1=.FALSE. 
                      IF ((lapw%rk(k,iintsp).LT.eps).AND.(.NOT.noco%l_ss)) THEN
                         l_lo1=.TRUE.
                      ELSE
                         l_lo1=.FALSE. 
                      ENDIF
                      CALL abccoflo(&
                           atoms,con1, rph(k,iintsp),cph(k,iintsp),ylm,n,&
                           na,k,lapw%nv(iintsp),l_lo1,alo1,blo1,clo1,&
                           nkvec(1,iintsp),enough, alo(-atoms%llod:,:,:,iintsp),blo(-atoms%llod:,:,:,iintsp),&
                           clo(-atoms%llod:,:,:,iintsp),kvec)
                      !-lo
                   ENDIF
                   !-->              synthesize the complex conjugates of a and b
                   DO l = 0,atoms%lnonsph(n)
                      ll1 = l* (l+1)
                      DO m = -l,l
                         termr = rph(k,iintsp)* REAL(ylm(ll1+m+1)) -&
                              cph(k,iintsp)*AIMAG(ylm(ll1+m+1))
                         termi = rph(k,iintsp)*AIMAG(ylm(ll1+m+1)) +&
                              cph(k,iintsp)* REAL(ylm(ll1+m+1))
                         ar(k,ll1+m,iintsp) = fj(k,l,n,spin2)*termr
                         ai(k,ll1+m,iintsp) = fj(k,l,n,spin2)*termi
                         br(k,ll1+m,iintsp) = gj(k,l,n,spin2)*termr
                         bi(k,ll1+m,iintsp) = gj(k,l,n,spin2)*termi
                      END DO
                   END DO

                END DO
                !--->       end loop over interstitial spin
             ENDDO

             !--->       add the local orbital contribution to the overlap and
             !--->       hamiltonian matrix, if they are used for this atom.

             IF (atoms%nlo(n).GE.1) THEN
                IF (atoms%n_u.GT.0) THEN
                   nkvecprevatu = nkvecprevats
                   iilou = iilos ; locolu = locols
                ENDIF
                IF ( noco%l_noco .AND. (.NOT.noco%l_ss) ) THEN
                   fjstart = 1
                   IF(l_socfirst) fjstart = isp
                   CALL slomat(&
                        input,atoms, n,na,lapw,con1,n_size,n_rank, gk,rph,cph,&
                        fj(:,0:,:,fjstart:), gj(:,0:,:,fjstart:),&
                        kvec,isp,usdus,alo1,blo1,clo1,noco, ab_dim,1,1,chi11,chi22,chi21,&
                        iilos,locols,nkvecprevats,.false.,bb_c=bbhlp)
                   CALL hlomat(input,atoms,isp,isp,n_size,n_rank,&
                        n,na,lapw,ar(:,0:,1),br(:,0:,1),ai(:,0:,1),bi(:,0:,1),&
                        el(:,n,isp),alo,blo,clo,usdus, noco,1,1,chi11,chi22,chi21,&
                        iiloh,locolh,nkvecprevath,tlmplm,.false.,aa_c=aahlp)
                ELSE
                   jd = 1 ; IF (noco%l_noco) jd = isp
                   DO iintsp = 1,nintsp
                      DO jintsp = 1,nintsp
                         CALL slomat(input,atoms,n,na,lapw,con1,n_size,n_rank,&
                              gk,rph,cph,fj,gj, kvec,isp,usdus,alo1,blo1,clo1,noco,&
                              ab_dim,iintsp,jintsp,chi11,chi22,chi21,&
                              iilos,locols,nkvecprevats,l_real,bb_r,bb_c)
                         CALL hlomat(input,atoms,isp,jd,n_size,n_rank,&
                              n,na,lapw,ar(:,0:,jintsp),br(:,0:,jintsp),ai(:,0:,jintsp),bi(:,0:,jintsp),&
                              el(:,n,isp),alo,blo,clo,usdus, noco,iintsp,jintsp,chi11,chi22,chi21,&
                              iiloh,locolh,nkvecprevath,tlmplm,l_real,aa_r,aa_c)
                      ENDDO
                   ENDDO
                ENDIF
             END IF


             IF (atoms%n_u.GT.0) THEN
                nkvecprevatuTemp = nkvecprevatu
                iilouTemp = iilou
                locoluTemp = locolu
                DO WHILE (i_u.LE.atoms%n_u)
                   IF (atoms%lda_u(i_u)%atomType.GT.n) EXIT
                   nkvecprevatuTemp = nkvecprevatu
                   iilouTemp = iilou
                   locoluTemp = locolu
                   IF (atoms%lda_u(i_u)%atomType.EQ.n) THEN
                      IF ((noco%l_noco).AND.(.NOT.noco%l_ss)) THEN
                         CALL u_ham(atoms,input,lapw,isp,n,i_u,invsfct,&
                                    ar,ai,br,bi,vs_mmp,lmaxb,&
                                    alo,blo,clo,&
                                    n_size,n_rank,isp,usdus,noco,&
                                    1,1,chi11,chi22,chi21,&
                                    nkvecprevatuTemp,iilouTemp,locoluTemp,.false.,aa_c=aahlp)
                      ELSE
                         DO iintsp = 1,nintsp
                            DO jintsp = 1,iintsp
                               CALL u_ham(atoms,input,lapw,isp,n,i_u,invsfct,&
                                          ar,ai,br,bi,vs_mmp,lmaxb,&
                                          alo,blo,clo,&
                                          n_size,n_rank,isp,usdus,noco,&
                                          iintsp,jintsp,chi11,chi22,chi21,&
                                          nkvecprevatuTemp,iilouTemp,locoluTemp,l_real,aa_r,aa_c)
                            END DO
                         END DO
                      END IF
                   END IF
                   i_u = i_u + 1
                END DO
                nkvecprevatu = nkvecprevatuTemp
                iilou = iilouTemp
                locolu = locoluTemp
             END IF

          ENDIF                                ! sym%invsat(na) = 0 or 1
          !--->    end loop over equivalent atoms
       END DO
       IF ( noco%l_noco .AND. (.NOT. noco%l_ss) ) CALL hsmt_hlptomat(atoms%nlotot,lapw%nv,sub_comm,chi11,chi21,chi22,aahlp,aa_c,bbhlp,bb_c)
       !---> end loop over atom types (ntype)
    ENDDO ntype_loop
    lapw%nmat = lapw%nmat + nkvecprevats

    RETURN
  END SUBROUTINE hsmt_extra
END MODULE m_hsmt_extra
