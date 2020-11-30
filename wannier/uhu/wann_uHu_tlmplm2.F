c*******************************************c
c   Muffin tin contribution to uHu-dmi      c
c  <u_{k+b1} | H_{k}^{mt} | u_{\theta+b2}>  c
c*******************************************c
c   Routine to set up T(lm,lmp) for all     c
c   pairs (b1,b2) and every atom n          c
c                                           c
c   Includes spherical and non-sph.         c
c   contributions to Hamiltonian:           c
c   H^{mt} ~ H^{sph} + H^{non-sph}          c
c                                           c
c   possible SOC contribution is dealt      c
c   with in a separate routine              c
c*******************************************c
c               J.-P. Hanke, Feb. 2016      c
c*******************************************c
c  TODO: local orbital part: epar(jspin)    c
c*******************************************c
      MODULE m_wann_uHu_tlmplm2
      CONTAINS
      SUBROUTINE wann_uHu_tlmplm2(
     >                  memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd,
     >                  ntype,dx,rmsh,jri,lmax,ntypsy,natd,
     >                  lnonsph,lmd,lmplmd,clnu,mlh,nmem,llh,nlh,neq,
     >                  irank,mlotot,mlolotot,
     >                  vr,nlod,llod,loplod,ello,llo,nlo,lo1l,l_dulo,
     >                  ulo_der,f,g,flo,f_b,g_b,flo_b,
     >                  kdiff,kdiff2,nntot,nntot2,bmat,bbmat,vr0,epar,
     >                  invsat,jspin,jspin_b,
     >                  l_skip_sph,l_skip_non,l_skip_loc,
     <                  tuu,tud,tdu,tdd,
     >                  tuulo,tulou,tdulo,tulod,tuloulo)

      USE m_intgr, ONLY : intgr3
      USE m_radflo
      USE m_radfun
      USE m_tlo
      USE m_sphbes
      USE m_ylm
      USE m_gaunt, ONLY: gaunt1
      USE m_dujdr
      USE m_wann_uHu_radintsra5
      USE m_wann_uHu_tlo
      USE m_constants

      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: memd,nlhd,ntypsd,ntypd,jmtd,lmaxd,jspd
      INTEGER, INTENT (IN) :: lmd,lmplmd,ntype,irank
      INTEGER, INTENT (IN) :: nlod,llod,loplod,natd,mlotot,mlolotot
      INTEGER, INTENT (IN) :: nntot,nntot2,jspin,jspin_b
      LOGICAL, INTENT (IN) :: l_skip_sph,l_skip_non,l_skip_loc
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (OUT):: 
     >       tdd(0:lmd,0:lmd,ntypd,nntot*nntot2),
     >       tdu(0:lmd,0:lmd,ntypd,nntot*nntot2),
     >       tud(0:lmd,0:lmd,ntypd,nntot*nntot2),
     >       tuu(0:lmd,0:lmd,ntypd,nntot*nntot2)
      COMPLEX, INTENT (OUT):: 
     >       tdulo(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tuulo(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tulou(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tulod(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tuloulo(nlod,-llod:llod,nlod,-llod:llod,ntypd,nntot*nntot2)
      COMPLEX, INTENT (IN) :: clnu(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: lo1l(0:llod,ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: mlh(memd,0:nlhd,ntypsd),nlh(ntypsd)
      INTEGER, INTENT (IN) :: nmem(0:nlhd,ntypsd),llh(0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: jri(ntypd),lmax(ntypd),ntypsy(natd)
      INTEGER, INTENT (IN) :: lnonsph(ntypd),ulo_der(nlod,ntypd)
      INTEGER, INTENT (IN) :: invsat(natd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: vr(jmtd,0:nlhd,ntypd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd)
      REAL,    INTENT (IN) :: epar(0:lmaxd,ntypd,2)
      REAL,    INTENT (IN) :: vr0(jmtd,ntypd)
      REAL,    INTENT (IN) :: f  (ntypd,jmtd,2,0:lmaxd),
     >                        g  (ntypd,jmtd,2,0:lmaxd),
     >                        f_b(ntypd,jmtd,2,0:lmaxd),
     >                        g_b(ntypd,jmtd,2,0:lmaxd),
     >                        flo  (ntypd,jmtd,2,nlod),
     >                        flo_b(ntypd,jmtd,2,nlod)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN) :: kdiff (3,nntot)
      REAL,    INTENT (IN) :: kdiff2(3,nntot2)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd)
C     ..
C     .. Local Scalars ..
      COMPLEX cil,cci,cil1,fac1,fac2,fac3,cilp,cciljj1,ciljj2
      INTEGER i,ii,l,l2,lamda,lh,lm,lmin,lmp,lp
      INTEGER lp1,lpl,m,mem,mems,mp,mu,n,nh,noded,nodeu,nrec,nsym,na,ne
      INTEGER ljj1,ljj2,mjj1,mjj2,lo,lop
      INTEGER ltil1,ltil2,mtil1,mtil2
      INTEGER mlo, mlolo, iu
      INTEGER :: lmini,lmini2,lmini3
      INTEGER :: lmaxi,lmaxi2,lmaxi3
      INTEGER :: lmx,lmx2,lmx3
      INTEGER :: lzb,ll,llp,lljj1,lljj2,lmjj1,lmjj2
      INTEGER :: ikpt_b,ikpt_b2
      INTEGER :: lwn,indexx
      REAL gs,rk1,rk2,sign
      REAL gc1,gc2,gc3,invsfct,e,temp_val1,temp_val2
      REAL t1,t0,t_sphbes,t_ylm,t_radint,t_gggint,t_sphint,t_lo
      REAL tt1,tt2,tt3
      LOGICAL :: l_nns
C     ..
C     .. Local Arrays ..
      COMPLEX, ALLOCATABLE :: yl1(:),yl2(:)
      REAL, ALLOCATABLE :: dvd_non(:,:,:,:,:),dvu_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: uvd_non(:,:,:,:,:),uvu_non(:,:,:,:,:)
      REAL, ALLOCATABLE :: dvd_sph(:,:,:,:,:),dvu_sph(:,:,:,:,:)
      REAL, ALLOCATABLE :: uvd_sph(:,:,:,:,:),uvu_sph(:,:,:,:,:)
      REAL, ALLOCATABLE :: p(:,:,:,:),p_b(:,:,:,:),dp_b(:,:,:,:)
      REAL, ALLOCATABLE :: q(:,:,:,:),q_b(:,:,:,:),dq_b(:,:,:,:)
      REAL, ALLOCATABLE :: dp(:,:,:,:),dq(:,:,:,:)
      REAL, ALLOCATABLE :: x(:),vv(:),vm(:)
      REAL, ALLOCATABLE :: jlpp(:),jj1(:,:),jj2(:,:)
      REAL :: bpt(3),bpt2(3),bkrot(3)

#if (defined(CPP_MPI) && !defined(CPP_T90))
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,max,mod
C     ..

      t_sphbes = 0.
      t_ylm    = 0.
      t_radint = 0.
      t_sphint = 0.
      t_gggint = 0.
      t_lo = 0.
      tt1 = 0.
      tt2 = 0.
      tt3 = 0.

      cci = cmplx(0.,-1.)
      tulou = cmplx(0.,0.); tulod = cmplx(0.,0.)
      tdulo = cmplx(0.,0.); tuulo = cmplx(0.,0.); tuloulo = cmplx(0.,0.)
      tuu = cmplx(0.,0.); tdu = cmplx(0.,0.)
      tud = cmplx(0.,0.); tdd = cmplx(0.,0.)

      allocate( yl1((lmaxd+1)**2), yl2((lmaxd+1)**2) )
      allocate( dvd_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( dvu_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uvd_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uvu_non(1:nlhd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( dvd_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( dvu_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uvd_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uvu_sph(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( p(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( p_b(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( q(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( q_b(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( dp(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( dp_b(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( dq(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( dq_b(jmtd,2,0:lmaxd,0:lmaxd) )
      allocate( x(jmtd),vv(jmtd) )
      allocate( vm(jmtd) )
      allocate( jlpp(0:lmaxd), jj1(0:lmaxd,jmtd), jj2(0:lmaxd,jmtd) )

c***************************************c
c  begin 1st neighbor loop:   <k+b1|    c
c***************************************c
      do ikpt_b=1,nntot
       bpt = kdiff(:,ikpt_b)
       rk1 = sqrt(dot_product(bpt,matmul(bbmat,bpt)))
!       rk1= sqrt(dotirp(bpt,bpt,bbmat))
       !if(irank.eq.0) write(*,*)'rk1',rk1

c***************************************c
c  begin 2nd neighbor loop:   |k+b2>    c
c***************************************c
       do ikpt_b2=1,nntot2
        indexx = ikpt_b2 + (ikpt_b-1)*nntot2
        bpt2 = kdiff2(:,ikpt_b2)
       rk2 = sqrt(dot_product(bpt2,matmul(bbmat,bpt2)))
!        rk2= sqrt(dotirp(bpt2,bpt2,bbmat))
        if(rk2.ne.0.0) stop 'dmi but rk2.ne.0'

        na = 1
        mlo = 1 ; mlolo = 1

      ! loop over atoms
      DO 210 n = 1,ntype
       lwn = lmax(n)
       nsym = ntypsy(na)
       nh = nlh(nsym)

       !if(indexx.eq.1 .and. irank.eq.0) then
       ! write(*,*)'atom n=',n
       ! write(*,*)'spins=',jspin,jspin_b
       ! do l=0,lwn
       !  write(*,*)'l=',l,epar(l,n,jspin),epar(l,n,jspin_b)
       ! enddo
       !endif

       l_nns = l_skip_non .OR. (lnonsph(n).LT.0) .OR.
     >         ( (invsat(na).ne.0).and.(invsat(na).ne.1) )
       if(invsat(na).eq.0) invsfct = 1.0
       if(invsat(na).eq.1) invsfct = 2.0
       if(indexx.eq.1 .and. irank.eq.0) write(*,*)'invsfct=',invsfct

       ! set up spherical Bessel functions
       ! jj1(l,b1*r) and jj2(l,b2*r)
       call cpu_time(t0)
       do i=1,jri(n)
          gs = rk1*rmsh(i,n)
          call sphbes(lwn,gs,jlpp)
          jj1(:,i) = jlpp(:)

          gs = rk2*rmsh(i,n)
          call sphbes(lwn,gs,jlpp)
          jj2(:,i) = jlpp(:)

c          if(rk1.eq.0.0) then
c           if(jj1(0,i).ne.1.0) write(*,*)'jj1.ne.1',i,rmsh(i,n),jj1(0,i)
c           if(ANY(jj1(1:lwn,i).ne.0.0)) write(*,*)'jj1.ne.0',i,rmsh(i,n)
c          endif
          if(jj2(0,i).ne.1.0) stop 'dmi but jj2.ne.1'
          if(ANY(jj2(1:lwn,i).ne.0.0)) stop 'dmi but jj2.ne.0'
         
       enddo
       call cpu_time(t1)
       t_sphbes=t_sphbes + t1-t0

       ! set up spherical harmonics
       ! yl1(lm) for b1 and yl2(lm) for b2
       yl1 = cmplx(0.,0.)
       bkrot=MATMUL(bpt,bmat)
       call ylm4(lwn,bkrot,yl1)

       yl2 = cmplx(0.,0.)
       bkrot=MATMUL(bpt2,bmat)
       if(ANY(bkrot.ne.0.0)) stop 'dmi but bkrot.ne.0'
       call ylm4(lwn,bkrot,yl2)
       !if(indexx.eq.1 .and. irank.eq.0) write(*,*)'yl2',yl2(1)

       call cpu_time(t0)
       t_ylm=t_ylm + t0-t1

       ! set up products of radial functions and sph. Bessel
       DO ljj1 = 0, lwn
        DO l = 0, lwn
         DO i=1,jri(n)
          p(i,:,l,ljj1) = f(n,i,:,l)*jj1(ljj1,i)
          q(i,:,l,ljj1) = g(n,i,:,l)*jj1(ljj1,i)
          p_b(i,:,l,ljj1) = f_b(n,i,:,l)*jj2(ljj1,i)
          q_b(i,:,l,ljj1) = g_b(n,i,:,l)*jj2(ljj1,i)
         ENDDO
         CALL dujdr(jmtd,jri(n),rmsh(1,n),dx(n),f(n,:,:,l),
     >              jj1,rk1,ljj1,lmaxd,dp(:,:,l,ljj1))
         CALL dujdr(jmtd,jri(n),rmsh(1,n),dx(n),g(n,:,:,l),
     >              jj1,rk1,ljj1,lmaxd,dq(:,:,l,ljj1))
         CALL dujdr(jmtd,jri(n),rmsh(1,n),dx(n),f_b(n,:,:,l),
     >              jj2,rk2,ljj1,lmaxd,dp_b(:,:,l,ljj1))
         CALL dujdr(jmtd,jri(n),rmsh(1,n),dx(n),g_b(n,:,:,l),
     >              jj2,rk2,ljj1,lmaxd,dq_b(:,:,l,ljj1))
        ENDDO
       ENDDO

       if(.false. .and. irank.eq.0 .and. n.eq.1) then
        open(876,file='out_radff')
        open(877,file='out_radgg')
        open(878,file='out_raddf')
        open(879,file='out_raddg')
        open(880,file='out_radjj')
        do i=1,jri(n)
         do l=0,lwn
          write(876,*)rmsh(i,n),f(n,i,1,l),f(n,i,2,l)
          write(877,*)rmsh(i,n),g(n,i,1,l),g(n,i,2,l)
          write(878,*)rmsh(i,n),dp(i,1,l,1),dp(i,2,l,1)
          write(879,*)rmsh(i,n),dq(i,1,l,1),dq(i,2,l,1)
          write(880,*)rmsh(i,n),jj1(l,i),jj2(l,i)
         enddo
        enddo
        close(876)
        close(877)
        close(878)
        close(879)
        close(880)
        stop
       endif

c****************************************************c
c   compute radial integrals                         c
c   <u(l')jj1(ljj1) | v(lamda,nu) | u(l)jj2(ljj2)>   c
c   where v(lamda,nu) ~ V^{sph} + V^{non-sph}        c
c****************************************************c
       ! spherical part
       IF(.not.l_skip_sph) THEN
       uvu_sph = 0.
       uvd_sph = 0.
       dvu_sph = 0.
       dvd_sph = 0.
       DO ljj1 = 0, lwn
        DO lp = 0, lwn
         lmini = abs(lp-ljj1)
         lmaxi = lp+ljj1
         lmx = min(lmaxi,lwn)
         DO lh=lmini,lmx
          if(mod(lmaxi+lh,2).eq.1) cycle
          DO ljj2 = 0,0!lwn
             lmini2= abs(lh-ljj2)
             lmaxi2= lh+ljj2
             lmx2 = min(lmaxi2,lwn)
             DO l = lmini2, lmx2
                if(mod(lmaxi2+l,2).eq.1) cycle
                e=(epar(lp,n,jspin)+epar(l,n,jspin_b))/2.0!epar(lh,n,jspin)
                lzb = lh

                !if(irank.eq.0) write(*,*)lh,l,lp,ljj2,ljj1

                !if(irank.eq.0) write(*,*)'uvu',epar(lp,n,jspin)
                call wann_uHu_radintsra5(jmtd,jri(n),rmsh(1,n),dx(n),
     >            e,vr0(1,n),
     >            p(:,:,lp,ljj1),p_b(:,:,l,ljj2),
     >            dp(:,:,lp,ljj1),dp_b(:,:,l,ljj2),
     >            lmaxd,lzb,uvu_sph(lh,l,lp,ljj2,ljj1),irank)

                !if(irank.eq.0) write(*,*)'uvd',1.0
                call wann_uHu_radintsra5(jmtd,jri(n),rmsh(1,n),dx(n),
     >            e,vr0(1,n),
     >            p(:,:,lp,ljj1),q_b(:,:,l,ljj2),
     >            dp(:,:,lp,ljj1),dq_b(:,:,l,ljj2),
     >            lmaxd,lzb,uvd_sph(lh,l,lp,ljj2,ljj1),irank)

                !if(irank.eq.0) write(*,*)'dvu',0.0
                call wann_uHu_radintsra5(jmtd,jri(n),rmsh(1,n),dx(n),
     >            e,vr0(1,n),
     >            q(:,:,lp,ljj1),p_b(:,:,l,ljj2),
     >            dq(:,:,lp,ljj1),dp_b(:,:,l,ljj2),
     >            lmaxd,lzb,dvu_sph(lh,l,lp,ljj2,ljj1),irank)

c                ! TODO: check influence of symmetrization on OM and DMI
c                uvd_sph(lh,l,lp,ljj2,ljj1) = (uvd_sph(lh,l,lp,ljj2,ljj1)
c     >                                  +dvu_sph(lh,l,lp,ljj2,ljj1))/2.0
c                dvu_sph(lh,l,lp,ljj2,ljj1) = uvd_sph(lh,l,lp,ljj2,ljj1)

                !if(irank.eq.0) write(*,*)'dvd',epar(lp,n,jspin)
                call wann_uHu_radintsra5(jmtd,jri(n),rmsh(1,n),dx(n),
     >            e,vr0(1,n),
     >            q(:,:,lp,ljj1),q_b(:,:,l,ljj2),
     >            dq(:,:,lp,ljj1),dq_b(:,:,l,ljj2),
     >            lmaxd,lzb,dvd_sph(lh,l,lp,ljj2,ljj1),irank)
             ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       ENDIF


       ! nonspherical part
       IF(.not.l_nns) then
       uvu_non = 0.
       uvd_non = 0.
       dvu_non = 0.
       dvd_non = 0.
       DO ljj1 = 0, lwn
        DO ljj2 = 0, 0!lwn
         DO lp = 0, lwn
c          DO l = 0, lp
          DO l = 0, lwn
               DO lh = 1, nh
                     DO i = 1,jri(n)
                     x(i)=(p(i,1,lp,ljj1)*p_b(i,1,l,ljj2)
     >                    +p(i,2,lp,ljj1)*p_b(i,2,l,ljj2))
     >                    * vr(i,lh,n)
                     ENDDO
                     CALL intgr3(x(1:jri(n)),rmsh(1,n),dx(n),jri(n),
     >                           uvu_non(lh,l,lp,ljj2,ljj1))
                     DO i = 1,jri(n)
                     x(i)=(q(i,1,lp,ljj1)*p_b(i,1,l,ljj2)
     >                    +q(i,2,lp,ljj1)*p_b(i,2,l,ljj2))
     >                   * vr(i,lh,n)
                     ENDDO
                     CALL intgr3(x(1:jri(n)),rmsh(1,n),dx(n),jri(n),
     >                           dvu_non(lh,l,lp,ljj2,ljj1))
                     DO i = 1,jri(n)
                     x(i)=(p(i,1,lp,ljj1)*q_b(i,1,l,ljj2)
     >                    +p(i,2,lp,ljj1)*q_b(i,2,l,ljj2))
     >                   * vr(i,lh,n)
                     ENDDO
                     CALL intgr3(x(1:jri(n)),rmsh(1,n),dx(n),jri(n),
     >                           uvd_non(lh,l,lp,ljj2,ljj1))
                     DO i = 1,jri(n)
                     x(i)=(q(i,1,lp,ljj1)*q_b(i,1,l,ljj2)
     >                    +q(i,2,lp,ljj1)*q_b(i,2,l,ljj2))
     >                   * vr(i,lh,n)
                     ENDDO
                     CALL intgr3(x(1:jri(n)),rmsh(1,n),dx(n),jri(n),
     >                           dvd_non(lh,l,lp,ljj2,ljj1))
               ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       ! symmetry can be used if f=f_b and g=g_b
       ! i.e. if jspin = jspin_b
c       DO lp=0,lwn
c        DO l=lp+1,lwn
c            uvu_non(:,l,lp,:,:) = uvu_non(:,lp,l,:,:)
c            dvd_non(:,l,lp,:,:) = dvd_non(:,lp,l,:,:)
c            uvd_non(:,l,lp,:,:) = dvu_non(:,lp,l,:,:)
c            dvu_non(:,l,lp,:,:) = uvd_non(:,lp,l,:,:)
c        ENDDO
c       ENDDO
       ENDIF

       call cpu_time(t0)
       t_radint = t_radint + t0-t1 

       tuu(:,:,n,indexx) = cmplx(0.,0.)
       tdu(:,:,n,indexx) = cmplx(0.,0.)
       tud(:,:,n,indexx) = cmplx(0.,0.)
       tdd(:,:,n,indexx) = cmplx(0.,0.)

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
       DO ljj1=0,lwn
        lljj1 = ljj1 * (ljj1 + 1)
        cciljj1 = cci**ljj1
        DO mjj1=-ljj1,ljj1
         lmjj1 = lljj1 + mjj1
         fac1 = cciljj1*conjg(yl1(lmjj1 + 1))
         if(fac1.eq.0.0) cycle

         DO lp=0,lwn
          llp = lp*(lp+1)
          lmini = abs(ljj1-lp)
          lmaxi = ljj1+lp
          lmx = min(lmaxi,lwn)
          DO mp=-lp,lp
           lmp = llp+mp
           mtil1=mjj1+mp
           IF(abs(mtil1).gt.lmx) cycle

           DO ltil1=lmini,lmx
            ! Gaunt conditions (G1):
            if((mod(lmaxi+ltil1,2).eq.1).or.
     >         (abs(mtil1).gt.ltil1)) cycle

             gc1 = gaunt1(ltil1, ljj1, lp, 
     >                    mtil1, mjj1, mp, lmaxd)

             cil1 =  (ImagUnit**lp) * fac1 * gc1

         DO ljj2=0,0!lwn
          lljj2 = ljj2 * (ljj2 + 1)
          lmini2 = abs(ltil1-ljj2)
          lmaxi2 = ltil1+ljj2
          lmx2 = min(lmaxi2,lwn)
          DO mjj2=-ljj2,ljj2
           lmjj2 = lljj2 + mjj2
           m=mjj2+mtil1
           IF(abs(m).gt.lmx2) cycle

           DO l=lmini2,lmx2
             ! Gaunt conditions (G2):
             if((mod(lmaxi2+l,2).eq.1).or.
     >          (abs(m).gt.l)) cycle
            ll = l*(l+1)
            lm = ll+m

             gc2 = gaunt1( l, ljj2, ltil1,
     >                     m, mjj2, mtil1, lmaxd)

             ! set up prefactor
             cil =  cil1* ( ImagUnit ** (ljj2 - l) ) 
     +              * gc2 * conjg( yl2(lmjj2 + 1) )

             if(cil.eq.0.0) cycle

             ! additional factor from symmetrization (below)
c             cil = 0.5 * cil

             ! compute T-coefficients
             ! using symmetrized uvu,uvd,dvu,dvd integrals
             tuu(lm, lmp, n, indexx)
     >         = tuu(lm, lmp, n, indexx)
     >         + cil * ( uvu_sph(ltil1, l, lp, ljj2, ljj1) )

             tdd(lm, lmp, n, indexx)
     >         = tdd(lm, lmp, n, indexx)
     >         + cil * ( dvd_sph(ltil1, l, lp, ljj2, ljj1) )

             tud(lm, lmp, n, indexx)
     >         = tud(lm, lmp, n, indexx)
     >         + cil * ( uvd_sph(ltil1, l, lp, ljj2, ljj1) )

             tdu(lm, lmp, n, indexx)
     >         = tdu(lm, lmp, n, indexx)
     >         + cil * ( dvu_sph(ltil1, l, lp, ljj2, ljj1) )
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
       DO ltil2 = 0, lwn         
        DO mtil2 = -ltil2, ltil2
         DO lh = 1, nh
            lamda = llh(lh,nsym)
            lmini = abs(ltil2 - lamda)
            lmaxi = ltil2 + lamda
            lmx = min(lmaxi,lwn)
            mems = nmem(lh,nsym)
            DO mem = 1, mems
               mu = mlh(mem,lh,nsym)
               mtil1 = mtil2+mu
               IF(abs(mtil1).GT.lmx) CYCLE

               DO ltil1 = lmini, lmx

               ! Gaunt conditions (G1):
               if((mod(lmaxi+ltil1,2).eq.1).or.
     >            (abs(mtil1).gt.ltil1)) cycle

                  gc1 = gaunt1(ltil1, lamda, ltil2, 
     >                         mtil1,    mu, mtil2, lmaxd)
                  fac1 = conjg(clnu(mem,lh,nsym)) * gc1 * invsfct
                  if(fac1.eq.0.0) CYCLE

         DO lp = 0, lnonsph(n)
            llp = lp * (lp + 1)
            cilp = ImagUnit**lp
            DO mp = -lp, lp
               lmp = llp + mp
               mjj1 = mtil1 - mp
               IF(abs(mjj1).GT.lwn) CYCLE

               DO ljj1 = 0, lwn
                  lmini2 = abs(lp-ljj1)
                  lmaxi2 = lp+ljj1

                  ! Gaunt conditions (G2):
                  if((lmini2.gt.ltil1).or.(lmaxi2.lt.ltil1).or.
     >               (mod(lmaxi2+ltil1,2).eq.1).or.
     >               (abs(mjj1).gt.ljj1)) cycle

                  lljj1 = ljj1 * (ljj1 + 1)
                  lmjj1 = lljj1 + mjj1

                  gc2 = gaunt1(ltil1,  ljj1,  lp,
     >                         mtil1,  mjj1,  mp, lmaxd)
                  fac2 = fac1 * cilp * (cci**ljj1) * gc2
     >                 * conjg( yl1(lmjj1 + 1) )
                  IF(fac2.eq.0.0) CYCLE

         DO ljj2 = 0, 0!lwn
            lljj2 = ljj2 * (ljj2 + 1)
            ciljj2 = ImagUnit**ljj2
            lmini3 = abs(ltil2 - ljj2)
            lmaxi3 = ltil2 + ljj2
            lmx3 = min(lmaxi3,lnonsph(n))
            DO mjj2 = -ljj2,ljj2
               m = mjj2+mtil2
               IF(abs(m).GT.lmx3) CYCLE
               lmjj2 = lljj2 + mjj2
               fac3 = fac2 * ciljj2 * conjg(yl2(lmjj2 + 1))
               if(fac3.eq.0.0) CYCLE

            DO l = lmini3,lmx3

                  ! Gaunt conditions (G3):
                  if((mod(lmaxi3+l,2).eq.1).or.
     >               (abs(m).gt.l)) cycle

                  ll = l * (l + 1)
                  lm = ll + m

                  gc3 = gaunt1(l,  ljj2,  ltil2,
     >                         m,  mjj2,  mtil2, lmaxd)

                  ! set up prefactor
                   cil = fac3 * (cci**l) * gc3

                  ! compute T-coefficients
                  tuu(lm, lmp, n, indexx)
     >              = tuu(lm, lmp, n, indexx)
     >              + cil * uvu_non(lh, l, lp, ljj2, ljj1)

                  tdd(lm, lmp, n, indexx)
     >              = tdd(lm, lmp, n, indexx)
     >              + cil * dvd_non(lh, l, lp, ljj2, ljj1)

                  tud(lm, lmp, n, indexx)
     >              = tud(lm, lmp, n, indexx)
     >              + cil * uvd_non(lh, l, lp, ljj2, ljj1)

                  tdu(lm, lmp, n, indexx)
     >              = tdu(lm, lmp, n, indexx)
     >              + cil * dvu_non(lh, l, lp, ljj2, ljj1)
               ENDDO !l
            ENDDO !mjj2
         ENDDO !ljj2
               ENDDO !lp
            ENDDO !mjj1
         ENDDO !mjj1
               ENDDO !ltil1
            ENDDO !mem
         ENDDO !lh
        ENDDO !mtil2
       ENDDO !ltil2

  555  CONTINUE
       call cpu_time(t0)
       t_gggint = t_gggint + t0-t1

c************** LOCAL ORBITAL CONTRIBUTION ***************c
c   determine contributions to T-matrix due to loc. orb.  c
c*********************************************************c
       IF((nlo(n).GE.1).AND.(.NOT.l_skip_loc)) THEN
        !write(*,*)'TODO: include loc. orb. for DMI'
          call wann_uHu_tlo(
     >             memd,nlhd,ntypsd,jmtd,lmaxd,dx(n),rmsh(1,n),jri(n),
     >             lmax(n),ntypsy,natd,lnonsph(n),lmd,lmplmd,clnu,
     >             mlh,nmem,llh,nlh,irank,vr(1,0,n),nlod,llod,
     >             loplod,ello(1,n),
     >             llo(1,n),nlo(n),lo1l(0,n),l_dulo(1,n),ulo_der(1,n),
     >             f(n,1:,1:,0:),g(n,1:,1:,0:),flo(n,1:,1:,1:),
     >             f_b(n,1:,1:,0:),g_b(n,1:,1:,0:),flo_b(n,1:,1:,1:),
     >             vr0(1,n),epar(0,n,jspin),l_skip_sph,
     >             l_nns,rk1,rk2,invsfct,yl1,yl2,jj1,jj2,
     >             p,dp,p_b,dp_b,q,dq,q_b,dq_b,
     >             tuulo(0,1,-llod,n,indexx),
     >             tdulo(0,1,-llod,n,indexx),
     >             tulou(0,1,-llod,n,indexx),
     >             tulod(0,1,-llod,n,indexx),
     >             tuloulo(1,-llod,1,-llod,n,indexx))
          call cpu_time(t1)
          t_lo = t_lo + t1-t0
       ENDIF

       na = na + neq(n)
  210 ENDDO !loop over atoms

c        write(*,*)ikpt_b,ikpt_b2
c        write(*,*)'t_radint=',t_radint-tt1
c        write(*,*)'t_sphint=',t_sphint-tt2
c        write(*,*)'t_gggint=',t_gggint-tt3
c        tt1 = t_radint
c        tt2 = t_sphint
c        tt3 = t_gggint

       enddo !ikpt_b2
      enddo !ikpt_b

#if (defined(CPP_MPI) && !defined(CPP_T90))
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      deallocate( yl1, yl2 )
      deallocate( dvd_sph, dvu_sph, uvd_sph, uvu_sph ) 
      deallocate( dvd_non, dvu_non, uvd_non, uvu_non ) 
      deallocate( p, dp, p_b, dp_b, q, dq, q_b, dq_b )
      deallocate( x,vv )
      deallocate( vm )
      deallocate( jlpp, jj1, jj2 )

c      if(irank.eq.0) then
c       write(*,*)'t_sphbes=',t_sphbes
c       write(*,*)'t_ylm   =',t_ylm
c       write(*,*)'t_radint=',t_radint
c       write(*,*)'t_sphint=',t_sphint
c       write(*,*)'t_gggint=',t_gggint
c       write(*,*)'t_locorb=',t_lo
c      endif

      END SUBROUTINE wann_uHu_tlmplm2
      END MODULE m_wann_uHu_tlmplm2
