c****************************************c
c   Muffin tin contribution to uHu       c
c  < u_{k+b1} | H_{k}^{mt} | u_{k+b2} >  c
c****************************************c
c   Routine to set up T(lm,lmp) for all  c
c   pairs (b1,b2) and every atom n       c
c   due to local orbitals                c
c                                        c
c****************************************c
c               J.-P. Hanke, Dec. 2015   c
c****************************************c
      MODULE m_wann_uHu_tlo
      CONTAINS
      SUBROUTINE wann_uHu_tlo(
     >             memd,nlhd,ntypsd,jmtd,lmaxd,
     >             dx,rmsh,jri,lmax,ntypsy,natd,
     >             lnonsph,lmd,lmplmd,clnu,mlh,nmem,llh,nlh,
     >             irank,vr,nlod,llod,loplod,ello,llo,nlo,lo1l,l_dulo,
     >             ulo_der,f,g,flo,f_b,g_b,flo_b,vr0,epar,
     >             l_skip_sph,l_nns,
     >             rk1,rk2,invsfct,yl1,yl2,jj1,jj2,p,dp,p_b,dp_b,
     >             q,dq,q_b,dq_b,
     >             tuulo,tdulo,tulou,tulod,tuloulo)

      USE m_intgr, ONLY : intgr3
      USE m_radflo
      USE m_radfun
      USE m_tlo
      USE m_sphbes
      USE m_ylm
      USE m_gaunt, ONLY: gaunt1
      USE m_wann_uHu_radintsra3
      USE m_wann_uHu_radintsra5
      USE m_constants
      USE m_dujdr

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      REAL,    INTENT (IN) :: rk1,rk2,invsfct
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,jmtd,lmaxd
      INTEGER, INTENT (IN) :: lmd,lmplmd,irank
      INTEGER, INTENT (IN) :: nlod,llod,loplod,natd
      LOGICAL, INTENT (IN) :: l_skip_sph,l_nns
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (OUT):: tdulo(0:lmd,nlod,-llod:llod)
      COMPLEX, INTENT (OUT):: tuulo(0:lmd,nlod,-llod:llod)
      COMPLEX, INTENT (OUT):: tulod(0:lmd,nlod,-llod:llod)
      COMPLEX, INTENT (OUT):: tulou(0:lmd,nlod,-llod:llod)
      COMPLEX, INTENT (OUT):: tuloulo(nlod,-llod:llod,nlod,-llod:llod)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      COMPLEX, INTENT (IN) :: yl1((lmaxd+1)**2),yl2((lmaxd+1)**2)
      INTEGER, INTENT (IN) :: llo(nlod),nlo
      INTEGER, INTENT (IN) :: lo1l(0:llod)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri,lmax,ntypsy(natd)
      INTEGER, INTENT (IN) :: lnonsph,ulo_der(nlod)
      REAL,    INTENT (IN) :: dx,rmsh(jmtd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd)
      REAL,    INTENT (IN) :: ello(nlod)
      REAL,    INTENT (IN) :: epar(0:lmaxd)
      REAL,    INTENT (IN) :: vr0(jmtd)
      REAL,    INTENT (IN) :: f  (jmtd,2,0:lmaxd),
     >                        g  (jmtd,2,0:lmaxd),
     >                        f_b(jmtd,2,0:lmaxd),
     >                        g_b(jmtd,2,0:lmaxd),
     >                        flo  (jmtd,2,nlod),
     >                        flo_b(jmtd,2,nlod)
      REAL,    INTENT (IN) :: p(jmtd,2,0:lmaxd,0:lmaxd),
     >                        q(jmtd,2,0:lmaxd,0:lmaxd),
     >                        dp(jmtd,2,0:lmaxd,0:lmaxd),
     >                        dq(jmtd,2,0:lmaxd,0:lmaxd),
     >                        p_b(jmtd,2,0:lmaxd,0:lmaxd),
     >                        q_b(jmtd,2,0:lmaxd,0:lmaxd),
     >                        dp_b(jmtd,2,0:lmaxd,0:lmaxd),
     >                        dq_b(jmtd,2,0:lmaxd,0:lmaxd)
      REAL,    INTENT (IN) :: jj1(0:lmaxd,jmtd),jj2(0:lmaxd,jmtd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod)
C     ..
C     .. Local Scalars ..
      COMPLEX cil,ci,cil2
      INTEGER i,l,l2,lamda,lh,lm,lmin,lmp,lp
      INTEGER lp1,lpl,m,mem,mems,mp,mu,n,nh,noded,nodeu,nrec,nsym,na,ne
      INTEGER ljj1,ljj2,mjj1,mjj2,lo,lop
      INTEGER ltil1,ltil2,mtil1,mtil2
      INTEGER mlo, mlolo, iu
      INTEGER :: lmini,lmini2,lmini3
      INTEGER :: lmaxi,lmaxi2,lmaxi3
      INTEGER :: ll,llp,lljj1,lljj2,lmjj1,lmjj2
      INTEGER :: ikpt_b,ikpt_b2
      INTEGER :: lwn
      REAL gc1,gc2,gc3,e,sy1,sy2
      REAL t1,t0,t_sphbes,t_ylm,t_radint,t_gggint,t_sphint
      LOGICAL :: L_oldsym
C     ..
C     .. Local Arrays .. 
      REAL, ALLOCATABLE :: ulovulo_sph(:,:,:,:,:)
      REAL, ALLOCATABLE :: ulovulo_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: ulovu_sph(:,:,:,:,:),ulovu_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: ulovd_sph(:,:,:,:,:),ulovd_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: uvulo_sph(:,:,:,:,:),uvulo_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: dvulo_sph(:,:,:,:,:),dvulo_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: plo(:,:,:,:),dplo(:,:,:,:)
      REAL, ALLOCATABLE :: plo_b(:,:,:,:),dplo_b(:,:,:,:)
      REAL, ALLOCATABLE :: x(:)


#if (defined(CPP_MPI) && !defined(CPP_T90))
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,max,mod
C     ..

      l_oldsym=.false.
      sy1=1.0
      sy2=0.0
      if(l_oldsym) then
       sy1=0.5
       sy2=0.5
      endif

      t_sphbes = 0.
      t_ylm    = 0.
      t_radint = 0.
      t_sphint = 0.
      t_gggint = 0.

      lwn = lmax
      ci = cmplx(0.0,1.0)
      tdulo(:,:,:) = cmplx(0.0,0.0)
      tuulo(:,:,:) = cmplx(0.0,0.0)
      tulod(:,:,:) = cmplx(0.0,0.0)
      tulou(:,:,:) = cmplx(0.0,0.0)
      tuloulo(:,:,:,:) = cmplx(0.0,0.0)

      allocate( ulovulo_sph(0:lmaxd,0:lmaxd,0:lmaxd,nlod,nlod) )
      allocate( ulovu_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( ulovd_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( uvulo_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( dvulo_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( ulovulo_non(1:nlhd,0:lmaxd,0:lmaxd,nlod,nlod) )
      allocate( ulovu_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( ulovd_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( uvulo_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( dvulo_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( x(jmtd) )
      allocate( plo(jmtd,2,0:lmaxd,nlod), plo_b(jmtd,2,0:lmaxd,nlod) )
      allocate( dplo(jmtd,2,0:lmaxd,nlod) )
      allocate( dplo_b(jmtd,2,0:lmaxd,nlod) )

        DO lo = 1,nlo
         DO ljj1=0, lwn
          DO i=1,jri
           plo(i,:,ljj1,lo) = flo(i,:,lo)*jj1(ljj1,i)
           plo_b(i,:,ljj1,lo) = flo_b(i,:,lo)*jj2(ljj1,i)
          ENDDO
          CALL dujdr(jmtd,jri,rmsh(1),dx,flo(:,:,lo),
     >               jj1,rk1,ljj1,lmaxd,dplo(:,:,ljj1,lo))
          CALL dujdr(jmtd,jri,rmsh(1),dx,flo_b(:,:,lo),
     >               jj2,rk2,ljj1,lmaxd,dplo_b(:,:,ljj1,lo))
         ENDDO
        ENDDO

c****************************************************c
c   compute radial integrals                         c
c   <u(l')jj1(ljj1) | v(lamda,nu) | u(l)jj2(ljj2)>   c
c   where v(lamda,nu) ~ V^{sph} + V^{non-sph}        c
c****************************************************c
! lo - lo
        do lop = 1,nlo
         lp = llo(lop)
         do lo = 1,nlo
          l = llo(lo)
          do ljj1=0, lwn
           do ljj2=0, lwn

            if(l_skip_sph) GOTO 1232
            do lh=0,lwn
              lmini = abs(lp-ljj1)
              lmini2= abs(lh-ljj2)
              lmaxi = lp+ljj1
              lmaxi2= lh+ljj2
              if((lmini.gt.lh) .or. (lmaxi.lt.lh).or.
     >           (lmini2.gt.l) .or. (lmaxi2.lt.l).or.
     >           (mod(lp+lh+ljj1,2).eq.1).or.
     >           (mod( l+lh+ljj2,2).eq.1)) then
               ulovulo_sph(lh,ljj2,ljj1,lo,lop) = 0.
              else
                e = (ello(lo)+ello(lop))/2.0!epar(lh)
                if(.not.l_oldsym) then
                call wann_uHu_radintsra5(jmtd,jri,rmsh(1),dx,
     >            e,vr0(1),plo(:,:,ljj1,lop),plo_b(:,:,ljj2,lo),
     >            dplo(:,:,ljj1,lop),dplo_b(:,:,ljj2,lo),lmaxd,lh,
     >            ulovulo_sph(lh,ljj2,ljj1,lo,lop),irank)
                else
                call wann_uHu_radintsra3(jmtd,jri,rmsh(1),dx,
     >            ello(lo),vr0(1),plo(:,:,ljj1,lop),plo_b(:,:,ljj2,lo),
     >            dplo_b(:,:,ljj2,lo),lmaxd,lh,
     >            ulovulo_sph(lh,ljj2,ljj1,lo,lop))
                endif
              endif 
            enddo
 1232       continue

            if(l_nns) GOTO 1230
            do lh = 1, nh
             do i = 1,jri
              x(i) = (  plo(i,1,ljj1,lop)*plo_b(i,1,ljj2,lo)
     >                + plo(i,2,ljj1,lop)*plo_b(i,2,ljj2,lo) )
     >             * vr(i,lh)
             enddo
             CALL intgr3(x(1:jri),rmsh(1),dx,jri,
     >                   ulovulo_non(lh,ljj2,ljj1,lo,lop))
            enddo
 1230       continue
           enddo
          enddo
         enddo
        enddo

! apw - lo
         do lop = 1,nlo
          lp = llo(lop)
          do ljj1 = 0,lwn
           do ljj2 = 0,lwn
            do l = 0,lwn

            if(l_skip_sph) GOTO 1233
            do lh=0,lwn
              lmini = abs(lp-ljj1)
              lmini2= abs(lh-ljj2)
              lmaxi = lp+ljj1
              lmaxi2= lh+ljj2
              if((lmini.gt.lh) .or. (lmaxi.lt.lh).or.
     >           (lmini2.gt.l) .or. (lmaxi2.lt.l).or.
     >           (mod(lp+lh+ljj1,2).eq.1).or.
     >           (mod( l+lh+ljj2,2).eq.1)) then
               ulovu_sph(lh,l,ljj2,ljj1,lop) = 0.
               ulovd_sph(lh,l,ljj2,ljj1,lop) = 0.
              else
                e = (ello(lop)+epar(l))/2.0!epar(lh)
                if(.not.l_oldsym) then
                call wann_uHu_radintsra5(jmtd,jri,rmsh(1),dx,
     >            e,vr0(1),plo(:,:,ljj1,lop),p_b(:,:,l,ljj2),
     >            dplo(:,:,ljj1,lop),dp_b(:,:,l,ljj2),lmaxd,lh,
     >            ulovu_sph(lh,l,ljj2,ljj1,lop),irank)
                call wann_uHu_radintsra5(jmtd,jri,rmsh(1),dx,
     >            e,vr0(1),plo(:,:,ljj1,lop),q_b(:,:,l,ljj2),
     >            dplo(:,:,ljj1,lop),dq_b(:,:,l,ljj2),lmaxd,lh,
     >            ulovd_sph(lh,l,ljj2,ljj1,lop),irank)
                else
                call wann_uHu_radintsra3(jmtd,jri,rmsh(1),dx,
     >            epar(l),vr0(1),plo(:,:,ljj1,lop),p_b(:,:,l,ljj2),
     >            dp_b(:,:,l,ljj2),lmaxd,lh,
     >            ulovu_sph(lh,l,ljj2,ljj1,lop))
                call wann_uHu_radintsra3(jmtd,jri,rmsh(1),dx,
     >            epar(l),vr0(1),plo(:,:,ljj1,lop),q_b(:,:,l,ljj2),
     >            dq_b(:,:,l,ljj2),lmaxd,lh,
     >            ulovd_sph(lh,l,ljj2,ljj1,lop))
                endif
              endif 

              lmini = abs(l-ljj1)
              lmini2= abs(lh-ljj2)
              lmaxi = l+ljj1
              lmaxi2= lh+ljj2
              if((lmini.gt.lh) .or. (lmaxi.lt.lh).or.
     >           (lmini2.gt.lp) .or. (lmaxi2.lt.lp).or.
     >           (mod(l+lh+ljj1,2).eq.1).or.
     >           (mod(lp+lh+ljj2,2).eq.1)) then
               uvulo_sph(lh,l,ljj2,ljj1,lop) = 0.
               dvulo_sph(lh,l,ljj2,ljj1,lop) = 0.
              else
                e=(ello(lop)+epar(l))/2.0!epar(lh)
                if(.not.l_oldsym) then
                call wann_uHu_radintsra5(jmtd,jri,rmsh(1),dx,
     >            e,vr0(1),p(:,:,l,ljj1),plo_b(:,:,ljj2,lop),
     >            dp(:,:,l,ljj1),dplo_b(:,:,ljj2,lop),lmaxd,lh,
     >            uvulo_sph(lh,l,ljj2,ljj1,lop),irank)
                call wann_uHu_radintsra5(jmtd,jri,rmsh(1),dx,
     >            e,vr0(1),q(:,:,l,ljj1),plo_b(:,:,ljj2,lop),
     >            dq(:,:,l,ljj1),dplo_b(:,:,ljj2,lop),lmaxd,lh,
     >            dvulo_sph(lh,l,ljj2,ljj1,lop),irank)
                else
                call wann_uHu_radintsra3(jmtd,jri,rmsh(1),dx,
     >            ello(lop),vr0(1),p(:,:,l,ljj1),plo_b(:,:,ljj2,lop),
     >            dplo_b(:,:,ljj2,lop),lmaxd,lh,
     >            uvulo_sph(lh,l,ljj2,ljj1,lop))
                call wann_uHu_radintsra3(jmtd,jri,rmsh(1),dx,
     >            ello(lop),vr0(1),q(:,:,l,ljj1),plo_b(:,:,ljj2,lop),
     >            dplo_b(:,:,ljj2,lop),lmaxd,lh,
     >            dvulo_sph(lh,l,ljj2,ljj1,lop))
                endif
              endif 
            enddo
 1233       continue

            if(l_nns) GOTO 1231
            do lh = 1, nh
             do i = 1,jri
              x(i) = (  plo(i,1,ljj1,lop)*p_b(i,1,l,ljj2)
     >                + plo(i,2,ljj1,lop)*p_b(i,2,l,ljj2) )
     >             * vr(i,lh)
             enddo
             CALL intgr3(x(1:jri),rmsh(1),dx,jri,
     >                   ulovu_non(lh,l,ljj2,ljj1,lop))
             do i = 1,jri
              x(i) = (  plo(i,1,ljj1,lop)*q_b(i,1,l,ljj2)
     >                + plo(i,2,ljj1,lop)*q_b(i,2,l,ljj2) )
     >             * vr(i,lh)
             enddo
             CALL intgr3(x(1:jri),rmsh(1),dx,jri,
     >                   ulovd_non(lh,l,ljj2,ljj1,lop))
             do i = 1,jri
              x(i) = (  p(i,1,l,ljj1)*plo_b(i,1,ljj2,lop)
     >                + p(i,2,l,ljj1)*plo_b(i,2,ljj2,lop) )
     >             * vr(i,lh)
             enddo
             CALL intgr3(x(1:jri),rmsh(1),dx,jri,
     >                   uvulo_non(lh,l,ljj2,ljj1,lop))
             do i = 1,jri
              x(i) = (  q(i,1,l,ljj1)*plo_b(i,1,ljj2,lop)
     >                + q(i,2,l,ljj1)*plo_b(i,2,ljj2,lop) )
     >             * vr(i,lh)
             enddo
             CALL intgr3(x(1:jri),rmsh(1),dx,jri,
     >                   dvulo_non(lh,l,ljj2,ljj1,lop))
            enddo
 1231       continue
           enddo
          enddo
         enddo
        enddo

       if(l_skip_sph) GOTO 444
c************** SPHERICAL CONTRIBUTION *******************c
c   compute product of the two Gaunt coefficients         c
c   with the radial integrals (+prefactors)               c
c*********************************************************c
c   We deal with two Gaunt coefficients:                  c
c  G1 = G( (ltil1,mtil1), (ljj1,mjj1) , ( lp, mp) )       c
c  G2 = G( (l , m)      , (ljj2, mjj2), (ltil1, mtil1) )  c
c*********************************************************c
c   use Gaunt conditions to reduce number of operations.  c
c   coefficient G(L1,L2,L3) only nonzero if               c
c       a)  l1 + l2 + l3 = even                           c
c       b)  |l2-l3| <= l1 <= l2+l3                        c
c       c)  m1 = m2 + m3                                  c
c*********************************************************c

! lo-lo
       DO lop=1,nlo
        lp = llo(lop)
        DO mp=-lp,lp
         DO ltil1=0,lwn
          DO mtil1=-ltil1,ltil1
           mjj1=mtil1-mp
           DO ljj1=0,lwn
            lljj1 = ljj1 * (ljj1 + 1)
            lmjj1 = lljj1 + mjj1
            lmini = abs(ljj1-lp)
            lmaxi = ljj1+lp
            ! Gaunt conditions (G1):
            if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >         (mod(lp+ltil1+ljj1,2).ne.0).or.
     >         (abs(mjj1).gt.ljj1)) cycle

             gc1 = gaunt1(ltil1, ljj1, lp, 
     >                    mtil1, mjj1, mp, lmaxd)

         DO lo=1,nlo
          l = llo(lo)
          DO m=-l,l
           mjj2=m-mtil1
           DO ljj2=0,lwn
            lljj2 = ljj2 * (ljj2 + 1)
            lmjj2 = lljj2 + mjj2
            lmini2 = abs(ltil1-ljj2)
            lmaxi2 = ltil1+ljj2
             ! Gaunt conditions (G2):
             if((lmini2.gt.l).or.(lmaxi2.lt.l).or.
     >          (mod(l+ltil1+ljj2,2).ne.0).or.
     >          (abs(mjj2).gt.ljj2)) cycle

             gc2 = gaunt1( l, ljj2, ltil1,
     >                     m, mjj2, mtil1, lmaxd)

             ! set up prefactor
             cil =   ( ci ** (lp - l + ljj1 + ljj2) ) 
     +              *( (-1.0) **ljj1 )
     +              * gc1 * gc2
     +              * conjg( yl1(lmjj1 + 1) )
     +              * conjg( yl2(lmjj2 + 1) )

             ! additional factor from symmetrization (below)
c             cil = 0.5 * cil

             ! compute T-coefficients
             ! using symmetrized uvu,uvd,dvu,dvd integrals
             tuloulo(lo, m, lop, mp)
     >         = tuloulo(lo, m, lop, mp)
     >         + cil * ( sy1*ulovulo_sph(ltil1, ljj2, ljj1, lo, lop) 
     >                 + sy2*ulovulo_sph(ltil1, ljj1, ljj2, lop, lo) ) 

           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1

! apw-lo
       DO lop=1,nlo
        lp = llo(lop)
        DO mp=-lp,lp
         DO ltil1=0,lwn
          DO mtil1=-ltil1,ltil1
           mjj1=mtil1-mp
           DO ljj1=0,lwn
            lljj1 = ljj1 * (ljj1 + 1)
            lmjj1 = lljj1 + mjj1
            lmini = abs(ljj1-lp)
            lmaxi = ljj1+lp
            ! Gaunt conditions (G1):
            if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >         (mod(lp+ltil1+ljj1,2).ne.0).or.
     >         (abs(mjj1).gt.ljj1)) cycle

             gc1 = gaunt1(ltil1, ljj1, lp, 
     >                    mtil1, mjj1, mp, lmaxd)

         DO l=0,lwn
          ll = l*(l+1)
          DO m=-l,l
           lm = ll+m
           mjj2=m-mtil1
           DO ljj2=0,lwn
            lljj2 = ljj2 * (ljj2 + 1)
            lmjj2 = lljj2 + mjj2
            lmini2 = abs(ltil1-ljj2)
            lmaxi2 = ltil1+ljj2
             ! Gaunt conditions (G2):
             if((lmini2.gt.l).or.(lmaxi2.lt.l).or.
     >          (mod(l+ltil1+ljj2,2).ne.0).or.
     >          (abs(mjj2).gt.ljj2)) cycle

             gc2 = gaunt1( l, ljj2, ltil1,
     >                     m, mjj2, mtil1, lmaxd)

             ! set up prefactor
             cil =   ( ci ** (lp - l + ljj1 + ljj2) ) 
     +              *( (-1.0) **ljj1 )
     +              * gc1 * gc2
     +              * conjg( yl1(lmjj1 + 1) )
     +              * conjg( yl2(lmjj2 + 1) )

             ! additional factor from symmetrization (below)
c             cil = 0.5 * cil

             ! compute T-coefficients
             ! using symmetrized uvu,uvd,dvu,dvd integrals
             tulou(lm, lop, mp)
     >         = tulou(lm, lop, mp)
     >         + cil * ( sy1*ulovu_sph(ltil1, l, ljj2, ljj1, lop) 
     >                 + sy2*uvulo_sph(ltil1, l, ljj1, ljj2, lop) ) 

             tulod(lm, lop, mp)
     >         = tulod(lm, lop, mp)
     >         + cil * ( sy1*ulovd_sph(ltil1, l, ljj2, ljj1, lop) 
     >                 + sy2*dvulo_sph(ltil1, l, ljj1, ljj2, lop) ) 

           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1

       DO ltil1=0,lwn
        DO mtil1=-ltil1,ltil1
         DO l=0,lwn
          ll = l*(l+1)
          DO m=-l,l
           lm = ll+m
           mjj1=mtil1-m
           DO ljj1=0,lwn
            lljj1 = ljj1 * (ljj1 + 1)
            lmjj1 = lljj1 + mjj1
            lmini = abs(ljj1-l)
            lmaxi = ljj1+l
            ! Gaunt conditions (G1):
            if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >         (mod(l+ltil1+ljj1,2).ne.0).or.
     >         (abs(mjj1).gt.ljj1)) cycle

             gc1 = gaunt1(ltil1, ljj1, l, 
     >                    mtil1, mjj1, m, lmaxd)

         DO lop=1,nlo
          lp = llo(lop)
          DO mp=-lp,lp
           mjj2=mp-mtil1
           DO ljj2=0,lwn
            lljj2 = ljj2 * (ljj2 + 1)
            lmjj2 = lljj2 + mjj2
            lmini2 = abs(ltil1-ljj2)
            lmaxi2 = ltil1+ljj2
             ! Gaunt conditions (G2):
             if((lmini2.gt.lp).or.(lmaxi2.lt.lp).or.
     >          (mod(lp+ltil1+ljj2,2).ne.0).or.
     >          (abs(mjj2).gt.ljj2)) cycle

             gc2 = gaunt1( lp, ljj2, ltil1,
     >                     mp, mjj2, mtil1, lmaxd)

             ! set up prefactor
             cil =   ( ci ** (l - lp + ljj1 + ljj2) )
     +              *( (-1.0) **ljj1 )
     +              * gc1 * gc2
     +              * conjg( yl1(lmjj1 + 1) )
     +              * conjg( yl2(lmjj2 + 1) )

             ! additional factor from symmetrization (below)
c             cil = 0.5 * cil

             ! compute T-coefficients
             ! using symmetrized uvu,uvd,dvu,dvd integrals
             tuulo(lm, lop, mp)
     >         = tuulo(lm, lop, mp)
     >         + cil * ( sy1*uvulo_sph(ltil1, l, ljj2, ljj1, lop) 
     >                 + sy2*ulovu_sph(ltil1, l, ljj1, ljj2, lop) ) 

             tdulo(lm, lop, mp)
     >         = tdulo(lm, lop, mp)
     >         + cil * ( sy1*dvulo_sph(ltil1, l, ljj2, ljj1, lop) 
     >                 + sy2*ulovd_sph(ltil1, l, ljj1, ljj2, lop) ) 

           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1

  444  CONTINUE
       call cpu_time(t1)
       t_sphint = t_sphint + t1-t0
 
       IF( l_nns ) GOTO 555
c************** NON-SPHERICAL CONTRIBUTION ***************c
c   compute product of the three Gaunt coefficients       c
c   with the radial integrals (+prefactors)               c
c*********************************************************c
c   We deal with three Gaunt coefficients:                c
c  G1 = G( (ltil1,mtil1) , (lamda,mu)  , (ltil2,mtil2) )  c
c  G2 = G( (ltil1,mtil1) , (ljj1,mjj1) , (lp,mp)       )  c
c  G3 = G( (l,m)         , (ljj2,mjj2) , (ltil2,mtil2) )  c
c*********************************************************c
c   use Gaunt conditions to reduce number of operations.  c
c   coefficient G(L1,L2,L3) only nonzero if               c
c       a)  l1 + l2 + l3 = even                           c
c       b)  |l2-l3| <= l1 <= l2+l3                        c
c       c)  m1 = m2 + m3                                  c
c*********************************************************c

! lo-lo
       DO ltil1 = 0, lwn
        DO mtil1 = -ltil1, ltil1
         DO 120 lh = 1, nh
            lamda = llh(lh,nsym)
            mems = nmem(lh,nsym)
            DO 110 mem = 1, mems
               mu = mlh(mem,lh,nsym)
               mtil2 = mtil1-mu
               DO 1111 ltil2 = 0, lwn         
               lmini = abs(ltil2 - lamda)
               lmaxi = ltil2 + lamda

               ! Gaunt conditions (G1):
               if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >            (mod(ltil1+ltil2+lamda,2).ne.0).or.
     >            (abs(mtil2).gt.ltil2)) cycle

                  gc1 = gaunt1(ltil1, lamda, ltil2, 
     >                         mtil1,    mu, mtil2, lmaxd)

         DO 1001 lop = 1,nlo
            lp = llo(lop)
            if(lp.GT.lnonsph) cycle
            DO mp = -lp, lp
               mjj1 = mtil1 - mp
               DO 140 ljj1 = 0, lwn
                  lljj1 = ljj1 * (ljj1 + 1)
                  lmjj1 = lljj1 + mjj1
                  lmini2 = abs(lp-ljj1)
                  lmaxi2 = lp+ljj1

                  ! Gaunt conditions (G2):
                  if((lmini2.gt.ltil1).or.(lmaxi2.lt.ltil1).or.
     >               (mod(lp+ltil1+ljj1,2).ne.0).or.
     >               (abs(mjj1).gt.ljj1)) cycle

                  gc2 = gaunt1(ltil1,  ljj1,  lp,
     >                         mtil1,  mjj1,  mp, lmaxd)

         DO 100 lo = 1,nlo
            l = llo(lo)
            if(l.GT.lnonsph) cycle
            DO m = -l, l
               mjj2 = m - mtil2 
               DO 901 ljj2 = 0, lwn
                  lljj2 = ljj2 * (ljj2 + 1)
                  lmjj2 = lljj2 + mjj2
                  lmini3 = abs(ltil2 - ljj2)
                  lmaxi3 = ltil2 + ljj2

                  ! Gaunt conditions (G3):
                  if((lmini3.gt.l).or.(lmaxi3.lt.l).or.
     >               (mod(l+ltil2+ljj2,2).ne.0).or.
     >               (abs(mjj2).gt.ljj2)) cycle

                  gc3 = gaunt1(l,  ljj2,  ltil2,
     >                         m,  mjj2,  mtil2, lmaxd)

                  ! set up prefactor
                  cil =   ( ci ** (lp - l + ljj1 + ljj2) )
     +                   * conjg(clnu(mem,lh,nsym))
     +                   *( (-1.0) ** ljj1 )
     +                   * gc1 * gc2 * gc3
     +                   * conjg( yl1(lmjj1 + 1) )
     +                   * conjg( yl2(lmjj2 + 1) )
     +                   * invsfct

                  ! compute T-coefficients
                  tuloulo(lo, m, lop, mp)
     >              = tuloulo(lo, m, lop, mp)
     >              + cil * ulovulo_non(lh, ljj2, ljj1, lo, lop)

  901          CONTINUE !l
            ENDDO !mjj2
  100    CONTINUE !ljj2
  140          CONTINUE !lp
            ENDDO !mjj1
 1001    CONTINUE !mjj1
 1111          CONTINUE !ltil1
  110       CONTINUE !mem
  120    CONTINUE !lh
        ENDDO !mtil2
       ENDDO !ltil2

! apw-lo
       DO ltil1 = 0, lwn
        DO mtil1 = -ltil1, ltil1
         DO lh = 1, nh
            lamda = llh(lh,nsym)
            mems = nmem(lh,nsym)
            DO mem = 1, mems
               mu = mlh(mem,lh,nsym)
               mtil2 = mtil1-mu
               DO ltil2 = 0, lwn         
               lmini = abs(ltil2 - lamda)
               lmaxi = ltil2 + lamda

               ! Gaunt conditions (G1):
               if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >            (mod(ltil1+ltil2+lamda,2).ne.0).or.
     >            (abs(mtil2).gt.ltil2)) cycle

                  gc1 = gaunt1(ltil1, lamda, ltil2, 
     >                         mtil1,    mu, mtil2, lmaxd)

         DO lop = 1,nlo
            lp = llo(lop)
            if(lp.GT.lnonsph) cycle
            DO mp = -lp, lp
               mjj1 = mtil1 - mp
               DO ljj1 = 0, lwn
                  lljj1 = ljj1 * (ljj1 + 1)
                  lmjj1 = lljj1 + mjj1
                  lmini2 = abs(lp-ljj1)
                  lmaxi2 = lp+ljj1

                  ! Gaunt conditions (G2):
                  if((lmini2.gt.ltil1).or.(lmaxi2.lt.ltil1).or.
     >               (mod(lp+ltil1+ljj1,2).ne.0).or.
     >               (abs(mjj1).gt.ljj1)) cycle

                  gc2 = gaunt1(ltil1,  ljj1,  lp,
     >                         mtil1,  mjj1,  mp, lmaxd)

         DO l = 0,lnonsph
            ll = l*(l+1)
            DO m = -l, l
               lm = ll+m
               mjj2 = m - mtil2 
               DO ljj2 = 0, lwn
                  lljj2 = ljj2 * (ljj2 + 1)
                  lmjj2 = lljj2 + mjj2
                  lmini3 = abs(ltil2 - ljj2)
                  lmaxi3 = ltil2 + ljj2

                  ! Gaunt conditions (G3):
                  if((lmini3.gt.l).or.(lmaxi3.lt.l).or.
     >               (mod(l+ltil2+ljj2,2).ne.0).or.
     >               (abs(mjj2).gt.ljj2)) cycle

                  gc3 = gaunt1(l,  ljj2,  ltil2,
     >                         m,  mjj2,  mtil2, lmaxd)

                  ! set up prefactor
                  cil =   ( ci ** (lp - l + ljj1 + ljj2) )
     +                   * conjg(clnu(mem,lh,nsym))
     +                   *( (-1.0) ** ljj1 )
     +                   * gc1 * gc2 * gc3
     +                   * conjg( yl1(lmjj1 + 1) )
     +                   * conjg( yl2(lmjj2 + 1) )
     +                   * invsfct

                  ! compute T-coefficients
                  tulou(lm, lop, mp)
     >              = tulou(lm, lop, mp)
     >              + cil * ulovu_non(lh, l, ljj2, ljj1, lop)

                  tulod(lm, lop, mp)
     >              = tulod(lm, lop, mp)
     >              + cil * ulovd_non(lh, l, ljj2, ljj1, lop)

               ENDDO !l
            ENDDO !mjj2
         ENDDO !ljj2
               ENDDO !lp
            ENDDO !mjj1
         ENDDO !mjj1

         DO l = 0,lnonsph
            ll = l*(l+1)
            DO m = -l, l
               lm = ll+m
               mjj1 = mtil1 - m
               DO ljj1 = 0, lwn
                  lljj1 = ljj1 * (ljj1 + 1)
                  lmjj1 = lljj1 + mjj1
                  lmini2 = abs(l-ljj1)
                  lmaxi2 = l+ljj1

                  ! Gaunt conditions (G2):
                  if((lmini2.gt.ltil1).or.(lmaxi2.lt.ltil1).or.
     >               (mod(l+ltil1+ljj1,2).ne.0).or.
     >               (abs(mjj1).gt.ljj1)) cycle

                  gc2 = gaunt1(ltil1,  ljj1,  l,
     >                         mtil1,  mjj1,  m, lmaxd)

         DO lop = 1,nlo
            lp = llo(lop)
            if(lp.GT.lnonsph) cycle
            DO mp = -lp, lp
               mjj2 = mp - mtil2 
               DO ljj2 = 0, lwn
                  lljj2 = ljj2 * (ljj2 + 1)
                  lmjj2 = lljj2 + mjj2
                  lmini3 = abs(ltil2 - ljj2)
                  lmaxi3 = ltil2 + ljj2

                  ! Gaunt conditions (G3):
                  if((lmini3.gt.lp).or.(lmaxi3.lt.lp).or.
     >               (mod(lp+ltil2+ljj2,2).ne.0).or.
     >               (abs(mjj2).gt.ljj2)) cycle

                  gc3 = gaunt1(lp,  ljj2,  ltil2,
     >                         mp,  mjj2,  mtil2, lmaxd)

                  ! set up prefactor
                  cil =   ( ci ** (l - lp + ljj1 + ljj2) )
     +                   * conjg(clnu(mem,lh,nsym))
     +                   *( (-1.0) ** ljj1 )
     +                   * gc1 * gc2 * gc3
     +                   * conjg( yl1(lmjj1 + 1) )
     +                   * conjg( yl2(lmjj2 + 1) )
     +                   * invsfct

                  ! compute T-coefficients
                  tuulo(lm, lop, mp)
     >              = tuulo(lm, lop, mp)
     >              + cil * uvulo_non(lh, l, ljj2, ljj1, lop)

                  tdulo(lm, lop, mp)
     >              = tdulo(lm, lop, mp)
     >              + cil * dvulo_non(lh, l, ljj2, ljj1, lop)

               ENDDO !l
            ENDDO !mjj2
         ENDDO !ljj2
               ENDDO !lp
            ENDDO !mjj1
         ENDDO !mjj1

               ENDDO !ltil2
            ENDDO !mem
         ENDDO !lh
        ENDDO !mtil1
       ENDDO !ltil1

  555  CONTINUE
       call cpu_time(t0)
       t_gggint = t_gggint + t0-t1

#if (defined(CPP_MPI) && !defined(CPP_T90))
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      deallocate( ulovulo_sph, ulovulo_non)
      deallocate( ulovu_sph, uvulo_sph, ulovd_sph, dvulo_sph )
      deallocate( ulovu_non, uvulo_non, ulovd_non, dvulo_non )
      deallocate( plo, dplo, plo_b, dplo_b )
      deallocate( x )


      END SUBROUTINE wann_uHu_tlo
      END MODULE m_wann_uHu_tlo
