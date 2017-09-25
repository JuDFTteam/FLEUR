c*****************************************c
c    Spin-orbit coupling cont. to uHu     c
c  < u_{k+b1} | V_{k}^{soc} | u_{k+b2} >  c
c*****************************************c
c    Routine to set up T(lm,lmp) for all  c
c    pairs (b1,b2) and every atom n       c
c*****************************************c
c                J.-P. Hanke, Dec. 2015   c
c*****************************************c
      MODULE m_wann_uHu_soc
      CONTAINS
      SUBROUTINE wann_uHu_soc(
     >                  input,atoms,
     >                  ntypd,jmtd,lmaxd,jspd,ntype,dx,rmsh,jri,lmax,
     >                  natd,lmd,lmplmd,neq,irank,nlod,llod,loplod,
     >                  ello,llo,nlo,lo1l,l_dulo,ulo_der,f,g,flo,f_b,
     >                  g_b,flo_b,kdiff,kdiff2,nntot,nntot2,
     >                  bmat,bbmat,vr,epar,jspin,
     >                  jspin_b,jspins,l_spav,theta,phi,alph,beta,
     >                  l_noco,l_skip_loc,
     <                  tuu,tud,tdu,tdd,tuulo,tulou,tdulo,tulod,tuloulo)

      USE m_intgr, ONLY : intgr0
      USE m_sphbes
      USE m_ylm
      USE m_gaunt, ONLY: gaunt1
      USE m_sointg
      USE m_constants
      USE m_wann_uHu_soc_tlo
      USE m_types
      USE m_anglso

      IMPLICIT NONE

      REAL :: sgml
      EXTERNAL sgml
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,max,mod

      TYPE(t_input),INTENT(IN)   :: input
      TYPE(t_atoms),INTENT(IN)   :: atoms

C     .. Scalar Arguments ..
      LOGICAL, INTENT (IN) :: l_noco,l_skip_loc
      INTEGER, INTENT (IN) :: ntypd,jmtd,lmaxd,jspd,jspins
      INTEGER, INTENT (IN) :: lmd,lmplmd,ntype,irank,jspin,jspin_b
      INTEGER, INTENT (IN) :: nlod,llod,loplod,natd
      INTEGER, INTENT (IN) :: nntot,nntot2
      REAL,    INTENT (IN) :: theta,phi
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (INOUT):: 
     >       tdd(0:lmd,0:lmd,ntypd,nntot*nntot2),
     >       tdu(0:lmd,0:lmd,ntypd,nntot*nntot2),
     >       tud(0:lmd,0:lmd,ntypd,nntot*nntot2),
     >       tuu(0:lmd,0:lmd,ntypd,nntot*nntot2)
      COMPLEX, INTENT (INOUT)::
     >       tdulo(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tuulo(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tulou(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tulod(0:lmd,nlod,-llod:llod,ntypd,nntot*nntot2),
     >       tuloulo(nlod,-llod:llod,nlod,-llod:llod,ntypd,nntot*nntot2)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: lo1l(0:llod,ntypd),neq(ntypd)
      INTEGER, INTENT (IN) :: jri(ntypd),lmax(ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: alph(ntypd),beta(ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,max(jspd,2))
      REAL,    INTENT (IN) :: epar(0:lmaxd,ntypd,max(jspd,2))
      REAL,    INTENT (IN) :: vr(jmtd,ntypd,jspd)
      REAL,    INTENT (IN) :: f  (ntypd,jmtd,2,0:lmaxd),
     >                        g  (ntypd,jmtd,2,0:lmaxd),
     >                        f_b(ntypd,jmtd,2,0:lmaxd),
     >                        g_b(ntypd,jmtd,2,0:lmaxd),
     >                        flo  (ntypd,jmtd,2,nlod),
     >                        flo_b(ntypd,jmtd,2,nlod)
      REAL,    INTENT (IN) :: bmat(3,3),bbmat(3,3)
      REAL,    INTENT (IN) :: kdiff (3,nntot)
      REAL,    INTENT (IN) :: kdiff2(3,nntot2)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd),l_spav
C     ..
C     .. Local Scalars ..
      REAL gs,rk1,rk2,socscreen_fac
      REAL gc1,gc2,c,e,r0
      REAL t0,t1,t11,t22,t33,t44,t55,t66
      COMPLEX cil,ci,cci,cciljj1,ciljj2,fac1,fac2,fac3,dump1,dump2
      INTEGER i,l,l2,lh,lm,lmp,lp
      INTEGER lp1,lpl,m,mp,mu,n,nh,noded,nodeu,nrec,na,ne
      INTEGER ljj1,ljj2,mjj1,mjj2
      INTEGER ltil1,mtil1,mtil2
      INTEGER mlo, mlolo, iu, indexx
      INTEGER lmini,lmini2, lmaxi,lmaxi2,lmx,lmx2
      INTEGER ll,llp,lljj1,lljj2,lmjj1,lmjj2
      INTEGER lltil1,lmtil1,lmtil2,sp1,sp2
      INTEGER mmin,mmax,lwn, ikpt_b,ikpt_b2
      LOGICAL l_socscreen
C     ..
C     .. Local Arrays ..
      REAL, ALLOCATABLE :: dvd(:,:,:,:,:),dvu(:,:,:,:,:)
      REAL, ALLOCATABLE :: uvd(:,:,:,:,:),uvu(:,:,:,:,:)
      REAL, ALLOCATABLE :: p(:,:,:),p_b(:,:,:),q(:,:,:),q_b(:,:,:)
      REAL, ALLOCATABLE :: x(:)
      REAL, ALLOCATABLE :: jlpp(:),jj1(:,:),jj2(:,:)
      REAL, ALLOCATABLE :: v0(:),vso(:,:)
      REAL :: bpt(3),bpt2(3),bkrot(3)
      COMPLEX, ALLOCATABLE :: yl1(:),yl2(:)
      COMPLEX, ALLOCATABLE :: angso(:,:,:)
      INTEGER :: ispjsp(2)
      DATA ispjsp/1,-1/

#if (defined(CPP_MPI) && !defined(CPP_T90))
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif

      t11 = 0.; t22 = 0.; t33 = 0.; t44 = 0.; t55 = 0.; t66 = 0.
      ci  = cmplx(0.0, 1.0)
      cci = cmplx(0.0,-1.0)
      c = c_light(1.0)
      sp1 = ispjsp(jspin)
      sp2 = ispjsp(jspin_b)

      l_socscreen = .false. ; socscreen_fac = 1.0
      INQUIRE(file='socscreen',exist=l_socscreen)
      IF(l_socscreen) THEN
       OPEN(757,file='socscreen')
       READ(757,*)socscreen_fac
       CLOSE(757)
       IF(irank.eq.0) write(*,*)'socscreen with factor',socscreen_fac
      ENDIF

      allocate( dvd(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( dvu(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uvd(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uvu(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( p(jmtd,0:lmaxd,0:lmaxd), p_b(jmtd,0:lmaxd,0:lmaxd) )
      allocate( q(jmtd,0:lmaxd,0:lmaxd), q_b(jmtd,0:lmaxd,0:lmaxd) )
      allocate( x(jmtd) )
      allocate( jlpp(0:lmaxd), jj1(0:lmaxd,jmtd), jj2(0:lmaxd,jmtd) )
      allocate( v0(jmtd), vso(jmtd,2) )
      allocate( angso(0:lmd,0:lmd,ntypd) )
      allocate( yl1((lmaxd+1)**2), yl2((lmaxd+1)**2) )

      call cpu_time(t0)
      angso = cmplx(0.,0.)
c*******************************************c
c  set up < l m sp1 | sigma L | l mp sp2 >  c
c*******************************************c
       DO n=1,ntype
        DO l=1,lmaxd              ! TODO: exclude s-states?
         ll=l*(l+1)
         DO m=-l,l
          lm = ll+m
          DO mp=-l,l
           lmp = ll+mp
           if(l_noco) then
            angso(lm,lmp,n) = anglso(beta(n),alph(n),l,m,sp1,l,mp,sp2)
           else
            angso(lm,lmp,n) = anglso(theta,phi,l,m,sp1,l,mp,sp2)
           endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
       call cpu_time(t1)
       t11 = t11 + t1-t0

c***************************************c
c  begin 1st neighbor loop:   <k+b1|    c
c***************************************c
      do ikpt_b=1,nntot
       bpt = kdiff(:,ikpt_b)
       rk1 = sqrt(dot_product(bpt,matmul(bbmat,bpt)))
!       rk1= sqrt(dotirp(bpt,bpt,bbmat))

c***************************************c
c  begin 2nd neighbor loop:   |k+b2>    c
c***************************************c
       do ikpt_b2=1,nntot2
        indexx = ikpt_b2 + (ikpt_b-1)*nntot2

        bpt2 = kdiff2(:,ikpt_b2)
        rk2 = sqrt(dot_product(bpt2,matmul(bbmat,bpt2)))
!        rk2= sqrt(dotirp(bpt2,bpt2,bbmat))

      ! loop over atoms
      DO 210 n = 1,ntype
       lwn = lmax(n)
       r0 = rmsh(1,n)

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
       enddo
       call cpu_time(t1)
       t22 = t22 + t1-t0

       ! set up spherical harmonics
       ! yl1(lm) for b1 and yl2(lm) for b2
       yl1 = cmplx(0.,0.)
       bkrot=MATMUL(bpt,bmat)
!       call cotra3(bpt,bkrot,bmat)
       call ylm4(lwn,bkrot,yl1)

       yl2 = cmplx(0.,0.)
       bkrot=MATMUL(bpt2,bmat)
!       call cotra3(bpt2,bkrot,bmat)
       call ylm4(lwn,bkrot,yl2)
       call cpu_time(t0)
       t33 = t33 + t0-t1

       ! set up products of radial functions and sph. Bessel
       DO ljj1 = 0, lwn
        DO l = 0, lwn
         DO i=1,jri(n)
          p(i,l,ljj1) = f(n,i,1,l)*jj1(ljj1,i)
          q(i,l,ljj1) = g(n,i,1,l)*jj1(ljj1,i)
          p_b(i,l,ljj1) = f_b(n,i,1,l)*jj2(ljj1,i)
          q_b(i,l,ljj1) = g_b(n,i,1,l)*jj2(ljj1,i)
         ENDDO
        ENDDO
       ENDDO

c****************************************************c
c   compute radial integrals                         c
c     <u(l')jj1(ljj1) | vso(lh) | u(l)jj2(ljj2)>     c
c****************************************************c
       uvu = 0.0; uvd = 0.0; dvu = 0.0; dvd = 0.0
       DO lp = 0, lwn
        DO l = 0, lwn

        ! construct vso
        v0 = 0.
        IF (jspins.EQ.1) THEN
           v0(1:jri(n)) = vr(1:jri(n),n,1)
           e = ( epar(l,n,1) + epar(lp,n,1) )/2.
        ELSE
           DO i = 1,jri(n)
             v0(i) = (vr(i,n,1)+vr(i,n,jspins))/2.
           END DO
           e = ( epar(l,n,1)+epar(l,n,jspins) 
     >          +epar(lp,n,1)+epar(lp,n,jspins) )/4.
        END IF

        CALL sointg(n,e,vr(:,n,:),v0,atoms,input,vso)  

        ! spin-averaging
        if(l_spav)then
            DO i= 1,jmtd
             vso(i,1)= (vso(i,1)+vso(i,2))/2.
             vso(i,2)= vso(i,1)
           ENDDO        
        endif

        if(l_socscreen) vso(:,:) = vso(:,:)*socscreen_fac

        ! calculate radial integrals
         DO ljj1 = 0, lwn
              lmini = abs(lp-ljj1)
              lmaxi = lp+ljj1
              lmx = min(lmaxi,lwn)
          DO ljj2 = 0, lwn
           DO lh = lmini, lmx
              lmini2= abs(lh-ljj2)
              lmaxi2= lh+ljj2
              if((lmini2.gt.l) .or. (lmaxi2.lt.l).or.
     >           (mod(lmaxi+lh,2).eq.1).or.
     >           (mod(lmaxi2+l,2).eq.1) ) CYCLE

               DO i=1,jri(n)
                x(i) = p(i,lp,ljj1)*p_b(i,l,ljj2)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(n)),r0,dx(n),jri(n),
     >                     uvu(lh,l,lp,ljj2,ljj1))

               DO i=1,jri(n)
                x(i) = q(i,lp,ljj1)*p_b(i,l,ljj2)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(n)),r0,dx(n),jri(n),
     >                     dvu(lh,l,lp,ljj2,ljj1))

               DO i=1,jri(n)
                x(i) = p(i,lp,ljj1)*q_b(i,l,ljj2)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(n)),r0,dx(n),jri(n),
     >                     uvd(lh,l,lp,ljj2,ljj1))

               DO i=1,jri(n)
                x(i) = q(i,lp,ljj1)*q_b(i,l,ljj2)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(n)),r0,dx(n),jri(n),
     >                     dvd(lh,l,lp,ljj2,ljj1))
           ENDDO
          ENDDO
         ENDDO
        ENDDO!l
       ENDDO!lp
       call cpu_time(t1)
       t44 = t44 +t1-t0

c***************** SOC CONTRIBUTION **********************c
c   compute product of the two Gaunt coefficients         c
c   with the radial integrals (+prefactors)               c
c*********************************************************c
c   We deal with two Gaunt coefficients:                  c
c  G1 = G( (ltil1, mtil1), (ljj1, mjj1), (lp, mp)       ) c
c  G2 = G( (l , m)       , (ljj2, mjj2), (ltil1, mtil2) ) c
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
        DO mjj1 = -ljj1,ljj1
         lmjj1 = lljj1 + mjj1
         fac1 = cciljj1 * conjg(yl1(lmjj1+1))
         IF(fac1.eq.0.0) CYCLE

         DO lp=0,lwn
          llp = lp*(lp+1)
          lmini = abs(lp-ljj1)
          lmaxi = lp+ljj1
          lmx = min(lmaxi,lwn)
          DO mp=-lp,lp
           lmp = llp+mp
           mtil1=mjj1+mp
           IF(abs(mtil1).GT.lmx) CYCLE

           DO ltil1=lmini,lmx
            ! Gaunt conditions (G1):
            if((mod(lmaxi+ltil1,2).eq.1).or.
     >         (abs(mtil1).gt.ltil1)) cycle

             lltil1 = ltil1*(ltil1+1)
             lmtil1 = lltil1+mtil1
             mmin = max(-ltil1,mtil1-1)
             mmax = min( ltil1,mtil1+1)

             gc1 = gaunt1(ltil1, ljj1, lp, 
     >                    mtil1, mjj1, mp, lmaxd)
             fac2 = fac1 * gc1 * (ci**lp)

        DO ljj2=0,lwn
         lljj2 = ljj2 * (ljj2 + 1)
         ciljj2 = ci**ljj2
         lmini2 = abs(ltil1-ljj2)
         lmaxi2 = ltil1+ljj2
         lmx2 = min(lmaxi2,lwn)
         DO mjj2 = -ljj2,ljj2
          lmjj2 = lljj2 + mjj2
          fac3 = fac2 * ciljj2 * conjg(yl2(lmjj2+1))
          IF(fac3.eq.0.0) CYCLE

          DO mtil2=mmin,mmax     ! mtil1-1 <= mtil2 <= mtil1+1
           lmtil2 = lltil1+mtil2
           m=mjj2+mtil2
           IF(abs(m).GT.lmx2) cycle

           DO l=lmini2,lmx2
             ! Gaunt conditions (G2):
             if((mod(lmaxi2+l,2).eq.1).or.
     >          (abs(m).gt.l)) cycle
             ll = l*(l+1)
             lm = ll+m

             gc2 = gaunt1( l, ljj2, ltil1,
     >                     m, mjj2, mtil2, lmaxd)

             ! set up prefactor
             cil = fac3 * (cci**l) * gc2 * conjg(angso(lmtil1,lmtil2,n))

             ! compute T-coefficients
             tuu(lm, lmp, n, indexx)
     >         = tuu(lm, lmp, n, indexx)
     >         + cil * uvu(ltil1, l, lp, ljj2, ljj1)

             tdd(lm, lmp, n, indexx)
     >         = tdd(lm, lmp, n, indexx)
     >         + cil * dvd(ltil1, l, lp, ljj2, ljj1)

             tud(lm, lmp, n, indexx)
     >         = tud(lm, lmp, n, indexx)
     >         + cil * uvd(ltil1, l, lp, ljj2, ljj1)

             tdu(lm, lmp, n, indexx)
     >         = tdu(lm, lmp, n, indexx)
     >         + cil * dvu(ltil1, l, lp, ljj2, ljj1)
           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
        ENDDO!mtil2
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1
       call cpu_time(t0)
       t55 = t55 + t0-t1

c************** LOCAL ORBITAL CONTRIBUTION ***************c
c   determine contributions to T-matrix due to loc. orb.  c
c*********************************************************c
       IF((nlo(n).GE.1).AND.(.NOT.l_skip_loc)) THEN

          call wann_uHu_soc_tlo(
     >           input,atoms,
     >           ntypd,jmtd,lmaxd,jspd,ntype,dx,rmsh,jri,lmax,
     >           natd,lmd,irank,nlod,llod,loplod,ello,llo,nlo,
     >           lo1l,l_dulo,ulo_der,flo,flo_b,kdiff,kdiff2,
     >           nntot,nntot2,
     >           vr,epar,jspin,jspin_b,jspins,l_spav,theta,phi,
     >           yl1,yl2,jj1,jj2,p,p_b,q,q_b,n,angso(0,0,n),
     >           l_socscreen,socscreen_fac,
     >           tuulo(0,1,-llod,n,indexx),
     >           tdulo(0,1,-llod,n,indexx),
     >           tulou(0,1,-llod,n,indexx),
     >           tulod(0,1,-llod,n,indexx),
     >           tuloulo(1,-llod,1,-llod,n,indexx))

          call cpu_time(t1)
          t66 = t66 + t1-t0
       ENDIF

  210 ENDDO !loop over atoms

       enddo !ikpt_b2
      enddo !ikpt_b

#if (defined(CPP_MPI) && !defined(CPP_T90))
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      deallocate( dvd, dvu, uvd, uvu )
      deallocate( p, p_b, q, q_b )
      deallocate( x )
      deallocate( jlpp, jj1, jj2 )
      deallocate( yl1, yl2 )
      deallocate( v0, vso, angso )

      if(irank.eq.0) then
       write(*,*)'t_angso =',t11
       write(*,*)'t_sphbes=',t22
       write(*,*)'t_ylm   =',t33
       write(*,*)'t_radint=',t44
       write(*,*)'t_gaunt =',t55
       write(*,*)'t_locorb=',t66
      endif

      END SUBROUTINE wann_uHu_soc
      END MODULE m_wann_uHu_soc
