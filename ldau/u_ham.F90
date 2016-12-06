!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_uham
  USE m_juDFT
  !-------------------------------------------------------------------+
  !     For details see Eq.(21) of Shick et al. PRB 60, 10765 (1999)  |
  !     Sets up the LDA + U Hamilton matrix with eigenvalues in the   |
  !     diagonal:                                                     |
  !  s     ---       *s,j'  s,j    ..   *s,j'  s,j    s            s |
  ! H    = >      [ A      A    + <u|u> B      B    ] V     + d    e  |
  !  j',j  -- m,m'   l,m'   l,m          l,m'   l,m    m,m'    j,j' j |
  !                                                                   |
  !                                                  G.B. Oct. 2000   |
  !                                                                   |
  ! Version for 1st variation to be called from hssphn.               |
  ! now the Hamiltonmatrix is close packed aa()  and ar,ai and br,bi  |
  ! are provided instead of ahelp and bhelp.              Feb. 2001   |
  !-------------------------------------------------------------------+
CONTAINS
  SUBROUTINE u_ham(atoms, input,lapw,jsp,n,invsfct,&
       ar,ai,br,bi,vs_mmp,lmaxb, alo,blo,clo, n_size,n_rank,usp,ud,&
       noco,iintsp,jintsp,chi11,chi22,chi21, nkvecprevat,iilo,locol,l_real,aa_r,aa_c)

    USE m_types
    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_lapw),INTENT(IN)    :: lapw
    !     ..
    !     .. Scalar Arguments ..    
    INTEGER, INTENT (IN) :: lmaxb  
    INTEGER, INTENT (IN) :: jsp,n
    INTEGER, INTENT (IN) :: n_size,n_rank ,usp
    INTEGER, INTENT (IN) :: iintsp,jintsp
    INTEGER, INTENT (INOUT) :: nkvecprevat,iilo,locol
    REAL,    INTENT (IN) :: invsfct
    COMPLEX, INTENT (IN) :: chi11,chi22,chi21
    !     ..
    !     .. Array Arguments ..
    TYPE(t_usdus),INTENT(IN):: ud
    REAL,    INTENT (IN) :: ar(:,0:,:),ai(:,0:,:) !(dimension%nvd,0:dimension%lmd,ab_dim)
    REAL,    INTENT (IN) :: br(:,0:,:),bi(:,0:,:) !(dimension%nvd,0:dimension%lmd,ab_dim)
    COMPLEX, INTENT (IN) :: alo(-atoms%llod:,:,:,:)!(-llod:llod,2*(2*llod+1),nlod,ab_dim)
    COMPLEX, INTENT (IN) :: blo(-atoms%llod:,:,:,:)!(-llod:llod,2*(2*llod+1),nlod,ab_dim)
    COMPLEX, INTENT (IN) :: clo(-atoms%llod:,:,:,:)!(-llod:llod,2*(2*llod+1),nlod,ab_dim)
    COMPLEX,INTENT(IN):: vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins)

    LOGICAL, INTENT(IN)     :: l_real
    REAL,    OPTIONAL,INTENT (INOUT) :: aa_r(:)!(matsize)
    COMPLEX, OPTIONAL,INTENT (INOUT) :: aa_c(:)!(matsize)

    !     ..
    !     .. Local Scalars ..
    INTEGER  m,mp,ispin,itype,l ,igp,ll,llm,nvv,ii,n_l,lp,lop,ig
    INTEGER matel0,locol0,lo ,nkvec,nkvecp,nkvecprevlo,nkend
    INTEGER igp_m,i_invsf,nc
    COMPLEX chihlp
    REAL    dd,uulo,dulo,ulolo
    !     ..
    !     .. Local Arrays
    REAL, ALLOCATABLE :: ar_help(:,:),ai_help(:,:),va_help(:)
    REAL, ALLOCATABLE :: br_help(:,:),bi_help(:,:),vb_help(:)
    REAL, ALLOCATABLE :: cr_help(:,:),ci_help(:,:)
    COMPLEX, ALLOCATABLE :: ax(:,:),bx(:,:),cx(:,:)
    COMPLEX, ALLOCATABLE :: wa_help(:),wb_help(:),wc_help(:)
    COMPLEX, ALLOCATABLE :: wd_help(:),we_help(:)


    INTRINSIC cmplx,conjg

    i_invsf = NINT(invsfct)
    n_l = 0 
    DO itype = 1,n
       l = atoms%lda_u(itype)%l
       IF (l.GE.0) n_l = n_l + 1
    ENDDO
    nvv = lapw%nv(jsp)
    IF (noco%l_ss)  nvv = lapw%nv(jintsp)

    ALLOCATE ( ar_help(nvv,-l:l),ai_help(nvv,-l:l),va_help(-l:l),&
         br_help(nvv,-l:l),bi_help(nvv,-l:l),vb_help(-l:l) )
    !
    !-----------------------------------------------------------------------
    !
    !                __      mm'    *lm'                    __      mm'    *lm' . .
    ! A_help(G',m) = >     V       A     and B_help(G',m) = >     V       B    <u|u>
    !                --m'   l,n,s   G'                      --m'   l,n,s   G'
    !
    !------------------------------------------------------------------------
    ll = l * (l+1)
    DO m = -l,l

       DO mp = -l,l
          va_help(mp) = REAL( vs_mmp(m,mp,n_l,jsp) ) * invsfct
          vb_help(mp) = va_help(mp) * ud%ddn(l,n,jsp)
       ENDDO

       DO ig = 1, nvv
          ar_help(ig,m) = 0.0
          ai_help(ig,m) = 0.0
          br_help(ig,m) = 0.0
          bi_help(ig,m) = 0.0

          DO mp = -l,l
             ar_help(ig,m) = ar_help(ig,m) + va_help(mp) * ar(ig,ll+mp,jintsp)
             ai_help(ig,m) = ai_help(ig,m) - va_help(mp) * ai(ig,ll+mp,jintsp)
             br_help(ig,m) = br_help(ig,m) + vb_help(mp) * br(ig,ll+mp,jintsp)
             bi_help(ig,m) = bi_help(ig,m) - vb_help(mp) * bi(ig,ll+mp,jintsp)
          ENDDO

       ENDDO

       DO mp = -l,l
          va_help(mp) = AIMAG( vs_mmp(m,mp,n_l,jsp) ) * invsfct
          vb_help(mp) = va_help(mp) * ud%ddn(l,n,jsp)
       ENDDO
       DO ig = 1, nvv
          DO mp = -l,l
             ar_help(ig,m) = ar_help(ig,m) + va_help(mp) * ai(ig,ll+mp,jintsp)
             ai_help(ig,m) = ai_help(ig,m) + va_help(mp) * ar(ig,ll+mp,jintsp)
             br_help(ig,m) = br_help(ig,m) + vb_help(mp) * bi(ig,ll+mp,jintsp)
             bi_help(ig,m) = bi_help(ig,m) + vb_help(mp) * br(ig,ll+mp,jintsp)
          ENDDO
       ENDDO

    ENDDO
    !
    IF (noco%l_ss) nvv = lapw%nv(iintsp)
    !
    !--------------------------------------------
    !
    !  l,n,s    --        m   lm          m   lm
    ! H      =  >   A_help   A    + B_help   B
    !  G G'     --m       G'  G           G'  G 
    !
    !--------------------------------------------
    nc = 0
    DO ig = n_rank+1, nvv, n_size
       nc = nc + 1
       !
       IF (noco%l_ss) THEN
          IF (n_size > 1)  CALL juDFT_error("U+SS+EV-PAR",calledby ="u_ham")
          IF ( iintsp.EQ.1 .AND. jintsp.EQ.1 ) THEN     !---> spin-up spin-up part
             igp_m = ig
             chihlp = chi11
             ii = (ig-1)*(ig)/2
          ELSEIF ( iintsp.EQ.2 .AND. jintsp.EQ.2 ) THEN !---> spin-down spin-down part
             igp_m = ig
             chihlp = chi22
             ii = (lapw%nv(1)+atoms%nlotot+ig-1)*(lapw%nv(1)+atoms%nlotot+ig)/2 + &
                  lapw%nv(1)+atoms%nlotot
          ELSE                                          !---> spin-down spin-up part
             igp_m = lapw%nv(1)
             chihlp = chi21
             ii = (lapw%nv(1)+atoms%nlotot+ig-1)*(lapw%nv(1)+atoms%nlotot+ig)/2
          ENDIF
       ELSE                                             ! no spin-spiral
          chihlp = CMPLX(1.0,0.0)
          igp_m = ig
          !          ii = ig * (ig - 1) / 2
          ii = nc*(nc-1)/2*n_size - (nc-1)*(n_size-n_rank-1)
       ENDIF

       DO m = -l,l
          llm = ll + m

          if (l_real) THEN
             DO igp = 1, igp_m
                aa_r(ii+igp) = aa_r(ii+igp) + &
                     ar_help(igp,m) * ar(ig,llm,iintsp) + br_help(igp,m) * br(ig,llm,iintsp) -&
                     ai_help(igp,m) * ai(ig,llm,iintsp) - bi_help(igp,m) * bi(ig,llm,iintsp) 
             ENDDO ! igp 
          ELSE
             DO igp = 1, igp_m
                aa_c(ii+igp) = aa_c(ii+igp) + chihlp * CMPLX(&
                     ar_help(igp,m) * ar(ig,llm,iintsp) + br_help(igp,m) * br(ig,llm,iintsp) -&
                     ai_help(igp,m) * ai(ig,llm,iintsp) - bi_help(igp,m) * bi(ig,llm,iintsp) ,&
                     ar_help(igp,m) * ai(ig,llm,iintsp) + br_help(igp,m) * bi(ig,llm,iintsp) +&
                     ai_help(igp,m) * ar(ig,llm,iintsp) + bi_help(igp,m) * br(ig,llm,iintsp) )
             enddo
          ENDIF
       ENDDO   ! m
    ENDDO     ! ig = 1, nvv
    !--------------------------------------------------------------------------
    !
    !     local orbitals (if there are any with same l on this atom) 
    !                    for the meaning of iilo etc. see slomat
    !
    !--------------------------------------------------------------------------
    nkvecprevlo = 0
    !
    DO lo = 1, atoms%nlo(n)
       IF (l.EQ.atoms%llo(lo,n)) THEN
          IF (n_size > 1)  CALL juDFT_error("U+LO+EV-PAR",calledby="u_ham")
          !
          !--->    calculate the contribution to H (G_lo,G): help arrays
          !
          ALLOCATE ( cr_help(nvv,-l:l),ci_help(nvv,-l:l) )
          uulo =  ud%uulon(lo,n,usp) ! * invsfct
          dulo =  ud%dulon(lo,n,usp) ! * invsfct
          matel0 = iilo
          locol0 = locol
          ll = l * (l+1)
          DO m = -l,l
             DO mp = -l,l
                va_help(mp) = REAL( vs_mmp(m,mp,n_l,jsp) ) * uulo
                vb_help(mp) = REAL( vs_mmp(m,mp,n_l,jsp) ) * dulo
             ENDDO
             DO ig = 1, nvv
                cr_help(ig,m) = 0.0
                ci_help(ig,m) = 0.0
                DO mp = -l,l
                   cr_help(ig,m) = cr_help(ig,m) + va_help(mp) * ar(ig,ll+mp,jintsp)
                   ci_help(ig,m) = ci_help(ig,m) - va_help(mp) * ai(ig,ll+mp,jintsp)
                   cr_help(ig,m) = cr_help(ig,m) + vb_help(mp) * br(ig,ll+mp,jintsp)
                   ci_help(ig,m) = ci_help(ig,m) - vb_help(mp) * bi(ig,ll+mp,jintsp)
                ENDDO
             ENDDO
             DO mp = -l,l
                va_help(mp) = AIMAG( vs_mmp(m,mp,n_l,jsp) ) * uulo
                vb_help(mp) = AIMAG( vs_mmp(m,mp,n_l,jsp) ) * dulo
             ENDDO
             DO ig = 1, nvv
                DO mp = -l,l
                   cr_help(ig,m) = cr_help(ig,m) + va_help(mp) * ai(ig,ll+mp,jintsp)
                   ci_help(ig,m) = ci_help(ig,m) + va_help(mp) * ar(ig,ll+mp,jintsp)
                   cr_help(ig,m) = cr_help(ig,m) + vb_help(mp) * bi(ig,ll+mp,jintsp)
                   ci_help(ig,m) = ci_help(ig,m) + vb_help(mp) * br(ig,ll+mp,jintsp)
                ENDDO
             ENDDO

          ENDDO
          !
          !--->    now add to H (G_lo,G)
          !
          DO m = -l,l
             iilo = matel0
             locol = locol0
             DO nkvec = 1,i_invsf* (2*l+1)
                locol = locol + 1
                IF (MOD(locol-1,n_size).EQ.n_rank) THEN
                   !-t3e
                   IF (l_real) THEN
                      DO ig = 1,nvv
                         iilo = iilo + 1
                         aa_r(iilo) = aa_r(iilo) + &
                              ar_help(ig,m) * REAL(alo(m,nkvec,lo,iintsp)) + br_help(ig,m) * REAL(blo(m,nkvec,lo,iintsp)) +&
                              cr_help(ig,m) * REAL(clo(m,nkvec,lo,iintsp)) - ai_help(ig,m) *AIMAG(alo(m,nkvec,lo,iintsp)) -&
                              bi_help(ig,m) *AIMAG(blo(m,nkvec,lo,iintsp)) - ci_help(ig,m) *AIMAG(clo(m,nkvec,lo,iintsp))
                      ENDDO                     ! jump to end of row:
                   ELSE
                      DO ig = 1,nvv
                         iilo = iilo + 1
                         aa_c(iilo) = aa_c(iilo) + CMPLX(&
                              ar_help(ig,m) * REAL(alo(m,nkvec,lo,iintsp)) + br_help(ig,m) * REAL(blo(m,nkvec,lo,iintsp)) +&
                              cr_help(ig,m) * REAL(clo(m,nkvec,lo,iintsp)) - ai_help(ig,m) *AIMAG(alo(m,nkvec,lo,iintsp)) -&
                              bi_help(ig,m) *AIMAG(blo(m,nkvec,lo,iintsp)) - ci_help(ig,m) *AIMAG(clo(m,nkvec,lo,iintsp)) ,&
                              ai_help(ig,m) * REAL(alo(m,nkvec,lo,iintsp)) + bi_help(ig,m) * REAL(blo(m,nkvec,lo,iintsp)) +&
                              ci_help(ig,m) * REAL(clo(m,nkvec,lo,iintsp)) + ar_help(ig,m) *AIMAG(alo(m,nkvec,lo,iintsp)) +&
                              br_help(ig,m) *AIMAG(blo(m,nkvec,lo,iintsp)) + cr_help(ig,m) *AIMAG(clo(m,nkvec,lo,iintsp)) )
                      ENDDO                     ! jump to end of row:
                   ENDIF
                   iilo = iilo + nkvecprevat + nkvecprevlo + nkvec 
                ENDIF
             ENDDO
          ENDDO
          DEALLOCATE ( cr_help,ci_help )
          !
          !--->    calculate the contribution to H (G_lo,G_lo)
          !
          ALLOCATE ( ax(-l:l,2*(2*atoms%llod+1)),wc_help(-l:l), bx(-l:l,2*(2*atoms%llod+1)),wd_help(-l:l),&
               cx(-l:l,2*(2*atoms%llod+1)),we_help(-l:l), wa_help(-l:l),wb_help(-l:l))

          iilo = matel0
          locol = locol0
          DO nkvec = 1,i_invsf* (2*l+1)
             locol = locol + 1
             IF (MOD(locol-1,n_size).EQ.n_rank) THEN
                iilo = iilo + nvv + nkvecprevat
                !
                !---> with other LOs
                !
                DO lop = 1, lo
                   lp = atoms%llo(lop,n)
                   IF (l.EQ.lp) THEN
                      dd    = ud%ddn(l,n,jsp)       ! * invsfct
                      ulolo = ud%uloulopn(lo,lop,n,usp) ! * invsfct
                      nkend = (2*lp+1) * i_invsf
                      IF (lop.EQ.lo) nkend = nkvec
                      DO m = -l,l
                         DO mp = -l,l
                            wa_help(mp) = vs_mmp(m,mp,n_l,jsp) ! * invsfct
                            wb_help(mp) = vs_mmp(m,mp,n_l,jsp) * dd
                            wc_help(mp) = vs_mmp(m,mp,n_l,jsp) * uulo
                            wd_help(mp) = vs_mmp(m,mp,n_l,jsp) * dulo
                            we_help(mp) = vs_mmp(m,mp,n_l,jsp) * ulolo
                         ENDDO
                         DO nkvecp = 1, nkend
                            ax(m,nkvecp) = CMPLX(0.0,0.0)
                            bx(m,nkvecp) = CMPLX(0.0,0.0)
                            cx(m,nkvecp) = CMPLX(0.0,0.0)
                            DO mp = -l,l
                               ax(m,nkvecp) = ax(m,nkvecp) +&
                                    CONJG(alo(mp,nkvecp,lop,jintsp)) * wa_help(mp) +&
                                    CONJG(clo(mp,nkvecp,lop,jintsp)) * wc_help(mp)
                               bx(m,nkvecp) = bx(m,nkvecp) +&
                                    CONJG(blo(mp,nkvecp,lop,jintsp)) * wb_help(mp) +&
                                    CONJG(clo(mp,nkvecp,lop,jintsp)) * wd_help(mp)
                               cx(m,nkvecp) = cx(m,nkvecp) +&
                                    CONJG(alo(mp,nkvecp,lop,jintsp)) * wc_help(mp) +&
                                    CONJG(blo(mp,nkvecp,lop,jintsp)) * wd_help(mp) +&
                                    CONJG(clo(mp,nkvecp,lop,jintsp)) * we_help(mp)
                            ENDDO
                         ENDDO
                      ENDDO
                      !
                      !---> add to H, special bound for contribution of lo with with itself
                      !
                      DO nkvecp = 1,nkend  
                         iilo = iilo + 1
                         IF (l_real) THEN
                         DO m = -l,l
                            aa_r(iilo) = aa_r(iilo) + &
                                 REAL(alo(m,nkvec,lo,iintsp))* REAL(ax(m,nkvecp)) +&
                                 REAL(blo(m,nkvec,lo,iintsp))* REAL(bx(m,nkvecp)) +&
                                 REAL(clo(m,nkvec,lo,iintsp))* REAL(cx(m,nkvecp)) -&
                                 AIMAG(alo(m,nkvec,lo,iintsp))*AIMAG(ax(m,nkvecp)) -&
                                 AIMAG(blo(m,nkvec,lo,iintsp))*AIMAG(bx(m,nkvecp)) -&
                                 AIMAG(clo(m,nkvec,lo,iintsp))*AIMAG(cx(m,nkvecp)) 
                         ENDDO
                      ELSE
                         DO m = -l,l
                            aa_c(iilo) = aa_c(iilo) + alo(m,nkvec,lo,iintsp) * ax(m,nkvecp) +&
                                 blo(m,nkvec,lo,iintsp) * bx(m,nkvecp) + clo(m,nkvec,lo,iintsp) * cx(m,nkvecp) 
                         ENDDO
                      ENDIF
                      ENDDO
                   ELSE
                      iilo = iilo + i_invsf * (2*lp+1) ! = nkvecprevlo
                   ENDIF                               ! l = lp
                ENDDO
                !
             ENDIF ! T3e
             !-t3e
          END DO
          DEALLOCATE ( ax,bx,cx,wa_help,wb_help,wc_help,wd_help,we_help)
          nkvecprevlo = nkvecprevlo + i_invsf * (2*l+1)

       ELSE
          lp = i_invsf * (2*atoms%llo(lo,n)+1)
          iilo = iilo + (nvv + nkvecprevlo + nkvecprevat) * lp
          iilo = iilo + (lp+1) * lp / 2
          nkvecprevlo = nkvecprevlo + lp
          locol = locol + lp
       ENDIF  ! end of if(l=atoms%llo(lo,n))
    ENDDO     ! end of lo = 1,atoms%nlo loop
    nkvecprevat = nkvecprevat + nkvecprevlo

    DEALLOCATE ( ar_help,ai_help,br_help,bi_help,va_help,vb_help )
    !
  END SUBROUTINE u_ham
END MODULE m_uham
