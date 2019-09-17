MODULE m_uham21
!-------------------------------------------------------------------+
!     Sets up the spin-off-diagonal part of  the LDA + U Hamilton   |
!                                                                   +---------+
!  off   ---     s1,j'  s2,j   s1  s2      s1,j'  s2,j  s1 .s2           off  |
! H    = >    [ A      A     <u  |u  > +  A      B    <u  |u  > + ... ] V     |
!  j',j  -- m,m' l,m'   l,m                 l,m'  l,m                    m',m |
!                                                                   +---------+
!                                                  G.B. May 2019    |
!                                                                   |
! Version for 1st variation to be called from hssphn.               |
! now the Hamiltonmatrix is close packed aa()  and ar,ai and br,bi  |
! are provided instead of ahelp and bhelp.              Feb. 2001   |
!                                                                   |
! in case of ev || we need uu to store spin2 columns that would     | 
! appear in another MPI rank                                  gb`18 |
!-------------------------------------------------------------------+
! Adapted for new fleur version (Henning Janßen 19)                 !
! - no local orbitals                                               !
!-------------------------------------------------------------------!
      CONTAINS

      SUBROUTINE u_ham21(atoms,noco,input,mpi,lapw,denCoeffsOffdiag,&
                         i_u,n,iintsp,jintsp,jsp,invsfct,&
                         chi11,chi22,chi21,chi12,&
                         ar,br,ai,bi,alo,blo,clo,v_mmp21,hmat)

      !SUBROUTINE u_ham21(
      !                jspd,nvd,lmaxd,ntypd,lmd,matsize,llod,&
      !                nlod,jspins,nv,lda_u,jsp,n,invsfct,&
      !                ar,ai,br,bi,v_mmp21,n_u,lmaxb,&
      !                llo,nlo,alo,blo,clo,uun21,udn21,dun21,ddn21,&
      !                uulon21,dulon21,uloun21,ulodn21,uloulopn21,&
      !                n_size,n_rank,nlotot,l_noco,&
      !                l_ss,ab_dim,iintsp,jintsp,chi11,chi22,chi21,&
      !                nkvecprevat,iilo,locol,aa)
      USE m_types
      USE m_constants
      USE m_juDFT

      IMPLICIT NONE
!     ..
!     .. Type Arguments
      TYPE(t_atoms), INTENT(IN) :: atoms 
      TYPE(t_noco),  INTENT(IN) :: noco 
      TYPE(t_input), INTENT(IN) :: input 
      TYPE(t_lapw),  INTENT(IN) :: lapw
      TYPE(t_mpi),   INTENT(IN) :: mpi
      TYPE(t_denCoeffsOffDiag), INTENT(IN) :: denCoeffsOffdiag

!     ..
!     .. Scalar Arguments ..
      INTEGER, INTENT(IN) :: i_u, n

      REAL,    INTENT (IN) :: invsfct
      COMPLEX, INTENT (IN) :: chi11,chi22,chi21,chi12
!     ..
!     .. Array Arguments ..
      REAL,    INTENT (IN) :: ar(:,0:,:),ai(:,0:,:) !(dimension%nvd,0:dimension%lmd,ab_dim)
      REAL,    INTENT (IN) :: br(:,0:,:),bi(:,0:,:) !(dimension%nvd,0:dimension%lmd,ab_dim)
      COMPLEX, INTENT (IN) :: alo(-atoms%llod:,:,:,:)!(-llod:llod,2*(2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT (IN) :: blo(-atoms%llod:,:,:,:)!(-llod:llod,2*(2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT (IN) :: clo(-atoms%llod:,:,:,:)!(-llod:llod,2*(2*llod+1),nlod,ab_dim)
      COMPLEX, INTENT(IN)  :: v_mmp21(-lmaxU_const:lmaxU_const,-lmaxU_const:lmaxU_const,n_u)

      CLASS(t_mat), INTENT(INOUT) :: hmat
!     ..
!     .. Local Scalars ..
      INTEGER m,mp,ispin,itype,l,ig,igp,ll,llm,nvv,ii,lp,lop,ki,kj
      INTEGER matel0,locol0,lo,kp,nkvec,nkvecp,nkvecprevlo,nkend
      INTEGER igp_m,i_invsf,nc,nc2,ig2,i2,nr2
      COMPLEX chihlp,ab1,ab2
      REAL    dd,uulo,dulo,ulolo
      LOGICAL l_1sp,l_2sp
!     ..
!     .. Local Arrays
      COMPLEX, ALLOCATABLE :: a_help(:,:),v_help(:)
      COMPLEX, ALLOCATABLE :: b_help(:,:),c_help(:,:)
      COMPLEX, ALLOCATABLE :: ax(:,:),bx(:,:),cx(:,:),aax(:)
      COMPLEX, ALLOCATABLE :: wa_help(:),wb_help(:),wc_help(:)
      COMPLEX, ALLOCATABLE :: wd_help(:),we_help(:)

      i_invsf = nint(invsfct)
      l = atoms%lda_u(i_u)%l
      nvv = lapw%nv(jsp)
      IF (noco%l_ss)  nvv = lapw%nv(jintsp)
      
      ALLOCATE (a_help(nvv,-l:l),v_help(-l:l),b_help(nvv,-l:l),aax(nvv))
!
!-----------------------------------------------------------------------
!
!                __     off    * s2,G'                 __     off    * s2,G'
! A_help(G',m) = >     V      A     and B_help(G',m) = >     V      B   
!                --m'   m,m'   lm'                     --m'   m,m'   lm'
!
!------------------------------------------------------------------------
      ll = l * (l+1)
      DO m = -l,l
        DO mp = -l,l
          v_help(mp) = v_mmp21(m,mp,i_u) ! *invsfct 
        ENDDO
        DO ig = 1, nvv
          a_help(ig,m) = cmplx(0.0,0.0)
          b_help(ig,m) = cmplx(0.0,0.0)
          DO mp = -l,l
            a_help(ig,m) = a_help(ig,m) + v_help(mp) * cmplx(ar(ig,ll+mp,2),ai(ig,ll+mp,2)) ! 2
            b_help(ig,m) = b_help(ig,m) + v_help(mp) * cmplx(br(ig,ll+mp,2),bi(ig,ll+mp,2))
!            a_help(ig,m) = a_help(ig,m) + v_help(mp) * 
!     +                             cmplx(ar(ig,ll+mp,1),ai(ig,ll+mp,1)) ! 2
!            b_help(ig,m) = b_help(ig,m) + v_help(mp) *
!     +                             cmplx(br(ig,ll+mp,1),bi(ig,ll+mp,1))
          ENDDO
        ENDDO
      ENDDO
!
      IF (noco%l_ss) nvv = lapw%nv(iintsp)
!
!--------------------------------------------
!
!  l,n,s    --        m   lm   s2  s1     lm   s2  s1           m   lm
! H      =  >   A_help  [A   <u  |u  > + B   <u  |u  >] + B_help  [A   < ...
!  G G'     --m       G'  G,s1            G,s1                  G'  G,s1
!
!--------------------------------------------
      DO ki = 1, nvv
        IF(mpi%n_size > 1) CALL juDFT_error("EV-PAR uham21",calledby="u_ham21")
        IF(noco%l_ss) CALL juDFT_error("l_ss uham21",calledby="u_ham21")

        aax = 0.0
        DO m = -l,l
          llm = ll + m

!          ab1 = cmplx(ar(ig,llm,2),ai(ig,llm,2)) * denCoeffsOffDiag%uun21(l,n) +  ! 1 (see swap of spins in hssphn)
!     +          cmplx(br(ig,llm,2),bi(ig,llm,2)) * denCoeffsOffDiag%udn21(l,n) 
!          ab2 = cmplx(ar(ig,llm,2),ai(ig,llm,2)) * denCoeffsOffDiag%dun21(l,n) +
!     +          cmplx(br(ig,llm,2),bi(ig,llm,2)) * denCoeffsOffDiag%ddn21(l,n) 

          ab1 = cmplx(ar(ki,llm,1),ai(ki,llm,1)) * denCoeffsOffDiag%uun21(l,n) + & ! 1 (see swap of spins in hssphn)
                cmplx(br(ki,llm,1),bi(ki,llm,1)) * denCoeffsOffDiag%udn21(l,n) 
          ab2 = cmplx(ar(ki,llm,1),ai(ki,llm,1)) * denCoeffsOffDiag%dun21(l,n) + &
                cmplx(br(ki,llm,1),bi(ki,llm,1)) * denCoeffsOffDiag%ddn21(l,n) 

          DO kj = 1, lapw%nv(1) 
            aax(kj) = aax(kj)+ a_help(kj,m)*ab1 + b_help(kj,m)*ab2
          ENDDO ! igp 
        ENDDO   ! m

        DO kj = 1, ki
          IF(.NOT.hmat%l_real) THEN
            !spin-up/spin-up
            hmat(1,1)%data_c(kj,ki) = hmat(1,1)%data_c(kj,ki) + aax(kj) * chi11
            !spin-down/spin-down
            hmat(2,2)%data_c(kj,ki) = hmat(2,2)%data_c(kj,ki) + aax(kj) * chi22
            !spin-down/spin-up
            hmat(2,1)%data_c(kj,ki) = hmat(2,1)%data_c(kj,ki) + aax(kj) * chi21
            !spin-up/spin-down
            hmat(1,2)%data_c(kj,ki) = hmat(1,2)%data_c(kj,ki) + aax(kj) * chi12
          ENDIF
        ENDDO
      ENDDO

!-----------------------------------------------------------------------
! The Contributions from local orbitals start here
! THey are ignored for now (same as normal LDA+U)
!-----------------------------------------------------------------------

!!--------------------------------------------------------------------------
!!
!!     local orbitals (if there are any with same l on this atom) 
!!                    for the meaning of iilo etc. see slomat
!!
!!--------------------------------------------------------------------------
!      nkvecprevlo = 0
!!
!      DO lo = 1, nlo(n)
!         IF (l.EQ.llo(lo,n)) THEN
!         IF (n_size > 1) STOP 'U+LO+EV-PAR'
!!
!!--->    calculate the contribution to H (G_lo,G): help arrays
!!
!         ALLOCATE ( cr_help(nvv,-l:l),ci_help(nvv,-l:l) )
!         uulo =  uulon(lo,n) ! * invsfct
!         dulo =  dulon(lo,n) ! * invsfct
!         matel0 = iilo
!         locol0 = locol
!         ll = l * (l+1)
!         DO m = -l,l
!           DO mp = -l,l
!             va_help(mp) = real( vs_mmp(m,mp,i_u,jsp) ) * uulo
!             vb_help(mp) = real( vs_mmp(m,mp,i_u,jsp) ) * dulo
!           ENDDO
!           DO ig = 1, nvv
!             cr_help(ig,m) = 0.0
!             ci_help(ig,m) = 0.0
!             DO mp = -l,l
!              cr_help(ig,m) = cr_help(ig,m) + va_help(mp) *
!     +                                      ar(ig,ll+mp,jintsp)
!              ci_help(ig,m) = ci_help(ig,m) - va_help(mp) *
!     +                                      ai(ig,ll+mp,jintsp)
!              cr_help(ig,m) = cr_help(ig,m) + vb_help(mp) *
!     +                                      br(ig,ll+mp,jintsp)
!              ci_help(ig,m) = ci_help(ig,m) - vb_help(mp) *
!     +                                      bi(ig,ll+mp,jintsp)
!             ENDDO
!           ENDDO
!           DO mp = -l,l
!             va_help(mp) = aimag( vs_mmp(m,mp,i_u,jsp) ) * uulo
!             vb_help(mp) = aimag( vs_mmp(m,mp,i_u,jsp) ) * dulo
!           ENDDO
!           DO ig = 1, nvv
!             DO mp = -l,l
!              cr_help(ig,m) = cr_help(ig,m) + va_help(mp) *
!     +                                      ai(ig,ll+mp,jintsp)
!              ci_help(ig,m) = ci_help(ig,m) + va_help(mp) *
!     +                                      ar(ig,ll+mp,jintsp)
!              cr_help(ig,m) = cr_help(ig,m) + vb_help(mp) *
!     +                                      bi(ig,ll+mp,jintsp)
!              ci_help(ig,m) = ci_help(ig,m) + vb_help(mp) *
!     +                                      br(ig,ll+mp,jintsp)
!             ENDDO
!           ENDDO
!
!         ENDDO
!!
!!--->    now add to H (G_lo,G)
!!
!         DO m = -l,l
!            iilo = matel0
!            locol = locol0
!            DO nkvec = 1,i_invsf* (2*l+1)
!               locol = locol + 1
!               IF (mod(locol-1,n_size).EQ.n_rank) THEN
!!-t3e
!                  DO ig = 1,nvv
!                     iilo = iilo + 1
!#ifdef CPP_INVERSION
!                     aa(iilo) = aa(iilo) + 
!     +                    ar_help(ig,m) * real(alo(m,nkvec,lo,iintsp)) +
!     +                    br_help(ig,m) * real(blo(m,nkvec,lo,iintsp)) +
!     +                    cr_help(ig,m) * real(clo(m,nkvec,lo,iintsp)) -
!     +                    ai_help(ig,m) *aimag(alo(m,nkvec,lo,iintsp)) -
!     +                    bi_help(ig,m) *aimag(blo(m,nkvec,lo,iintsp)) -
!     +                    ci_help(ig,m) *aimag(clo(m,nkvec,lo,iintsp))
!#else
!                     aa(iilo) = aa(iilo) + cmplx(
!     +                    ar_help(ig,m) * real(alo(m,nkvec,lo,iintsp)) +
!     +                    br_help(ig,m) * real(blo(m,nkvec,lo,iintsp)) +
!     +                    cr_help(ig,m) * real(clo(m,nkvec,lo,iintsp)) -
!     +                    ai_help(ig,m) *aimag(alo(m,nkvec,lo,iintsp)) -
!     +                    bi_help(ig,m) *aimag(blo(m,nkvec,lo,iintsp)) -
!     +                    ci_help(ig,m) *aimag(clo(m,nkvec,lo,iintsp)) ,
!     +                    ai_help(ig,m) * real(alo(m,nkvec,lo,iintsp)) +
!     +                    bi_help(ig,m) * real(blo(m,nkvec,lo,iintsp)) +
!     +                    ci_help(ig,m) * real(clo(m,nkvec,lo,iintsp)) +
!     +                    ar_help(ig,m) *aimag(alo(m,nkvec,lo,iintsp)) +
!     +                    br_help(ig,m) *aimag(blo(m,nkvec,lo,iintsp)) +
!     +                    cr_help(ig,m) *aimag(clo(m,nkvec,lo,iintsp)) )
!#endif
!                  ENDDO                     ! jump to end of row:
!                  iilo = iilo + nkvecprevat + nkvecprevlo + nkvec 
!               ENDIF
!            ENDDO
!         ENDDO
!         DEALLOCATE ( cr_help,ci_help )
!!
!!--->    calculate the contribution to H (G_lo,G_lo)
!!
!         ALLOCATE ( ax(-l:l,2*(2*llod+1)),wc_help(-l:l),
!     +              bx(-l:l,2*(2*llod+1)),wd_help(-l:l),
!     +              cx(-l:l,2*(2*llod+1)),we_help(-l:l),
!     +              wa_help(-l:l),wb_help(-l:l))
!
!         iilo = matel0
!         locol = locol0
!         DO nkvec = 1,i_invsf* (2*l+1)
!            locol = locol + 1
!            IF (mod(locol-1,n_size).EQ.n_rank) THEN
!               iilo = iilo + nvv + nkvecprevat
!!
!!---> with other LOs
!!
!               DO lop = 1, lo
!                  lp = llo(lop,n)
!                  IF (l.EQ.lp) THEN
!                     dd    = ddn(l,n,jsp)       ! * invsfct
!                     ulolo = uloulopn(lo,lop,n) ! * invsfct
!                     nkend = (2*lp+1) * i_invsf
!                     IF (lop.EQ.lo) nkend = nkvec
!                     DO m = -l,l
!                       DO mp = -l,l
!                         wa_help(mp) = vs_mmp(m,mp,i_u,jsp) ! * invsfct
!                         wb_help(mp) = vs_mmp(m,mp,i_u,jsp) * dd
!                         wc_help(mp) = vs_mmp(m,mp,i_u,jsp) * uulo
!                         wd_help(mp) = vs_mmp(m,mp,i_u,jsp) * dulo
!                         we_help(mp) = vs_mmp(m,mp,i_u,jsp) * ulolo
!                       ENDDO
!                       DO nkvecp = 1, nkend
!                         ax(m,nkvecp) = cmplx(0.0,0.0)
!                         bx(m,nkvecp) = cmplx(0.0,0.0)
!                         cx(m,nkvecp) = cmplx(0.0,0.0)
!                         DO mp = -l,l
!                           ax(m,nkvecp) = ax(m,nkvecp) +
!     +                  conjg(alo(mp,nkvecp,lop,jintsp)) * wa_help(mp) +
!     +                  conjg(clo(mp,nkvecp,lop,jintsp)) * wc_help(mp)
!                           bx(m,nkvecp) = bx(m,nkvecp) +
!     +                  conjg(blo(mp,nkvecp,lop,jintsp)) * wb_help(mp) +
!     +                  conjg(clo(mp,nkvecp,lop,jintsp)) * wd_help(mp)
!                           cx(m,nkvecp) = cx(m,nkvecp) +
!     +                  conjg(alo(mp,nkvecp,lop,jintsp)) * wc_help(mp) +
!     +                  conjg(blo(mp,nkvecp,lop,jintsp)) * wd_help(mp) +
!     +                  conjg(clo(mp,nkvecp,lop,jintsp)) * we_help(mp)
!                         ENDDO
!                       ENDDO
!                     ENDDO
!!
!!---> add to H, special bound for contribution of lo with with itself
!!
!                     DO nkvecp = 1,nkend  
!                        iilo = iilo + 1
!                        DO m = -l,l
!#ifdef CPP_INVERSION
!                           aa(iilo) = aa(iilo) + 
!     +                real(alo(m,nkvec,lo,iintsp))* real(ax(m,nkvecp)) +
!     +                real(blo(m,nkvec,lo,iintsp))* real(bx(m,nkvecp)) +
!     +                real(clo(m,nkvec,lo,iintsp))* real(cx(m,nkvecp)) -
!     +               aimag(alo(m,nkvec,lo,iintsp))*aimag(ax(m,nkvecp)) -
!     +               aimag(blo(m,nkvec,lo,iintsp))*aimag(bx(m,nkvecp)) -
!     +               aimag(clo(m,nkvec,lo,iintsp))*aimag(cx(m,nkvecp)) 
!#else
!                           aa(iilo) = aa(iilo) + 
!     +                           alo(m,nkvec,lo,iintsp) * ax(m,nkvecp) +
!     +                           blo(m,nkvec,lo,iintsp) * bx(m,nkvecp) +
!     +                           clo(m,nkvec,lo,iintsp) * cx(m,nkvecp) 
!#endif
!                        ENDDO
!                     ENDDO
!                  ELSE
!                     iilo = iilo + i_invsf * (2*lp+1) ! = nkvecprevlo
!                  ENDIF                               ! l = lp
!               ENDDO
!!
!            ENDIF ! T3e
!!-t3e
!         END DO
!         DEALLOCATE ( ax,bx,cx,wa_help,wb_help,wc_help,wd_help,we_help)
!         nkvecprevlo = nkvecprevlo + i_invsf * (2*l+1)
!
!         ELSE
!           lp = i_invsf * (2*llo(lo,n)+1)
!           iilo = iilo + (nvv + nkvecprevlo + nkvecprevat) * lp
!           iilo = iilo + (lp+1) * lp / 2
!           nkvecprevlo = nkvecprevlo + lp
!           locol = locol + lp
!         ENDIF  ! end of if(l=llo(lo,n))
!      ENDDO     ! end of lo = 1,nlo loop
!      nkvecprevat = nkvecprevat + nkvecprevlo

      DEALLOCATE ( a_help,b_help,v_help )
!
      END SUBROUTINE u_ham21
      END MODULE m_uham21
