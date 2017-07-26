c*****************************************c
c    Spin-orbit coupling cont. to uHu     c
c  < u_{k+b1} | V_{k}^{soc} | u_{k+b2} >  c
c*****************************************c
c    Routine to set up T(lm,lmp) for all  c
c    pairs (b1,b2) and every atom n       c
c    due to local orbitals                c
c*****************************************c
c                J.-P. Hanke, Dec. 2015   c
c*****************************************c
      MODULE m_wann_uHu_soc_tlo
      CONTAINS
      SUBROUTINE wann_uHu_soc_tlo(
     >                  ntypd,jmtd,lmaxd,jspd,
     >                  ntype,dx,rmsh,jri,lmax,natd,
     >                  lmd,irank,nlod,llod,
     >                  loplod,ello,llo,nlo,lo1l,l_dulo,ulo_der,
     >                  flo,flo_b,kdiff,kdiff2,nntot,nntot2,
     >                  vr,epar,jspin,jspin_b,jspins,l_spav,theta,phi,
     >                  yl1,yl2,jj1,jj2,p,p_b,q,q_b,ntyp,angso,
     >                  l_socscreen,socscreen_fac,
     <                  tuulo,tdulo,tulou,tulod,tuloulo)

      USE m_intgr, ONLY : intgr0
      USE m_sphbes
      USE m_dotir
      USE m_ylm
      USE m_cotra
      USE m_gaunt, ONLY: gaunt1
      USE m_sointg
      USE m_constants, ONLY : c_light
      

      IMPLICIT NONE
C     .. Intrinsic Functions ..
      INTRINSIC abs,cmplx,max,mod
C     ..
C     .. Scalar Arguments ..
      INTEGER, INTENT (IN) :: ntypd,jmtd,lmaxd,jspd,jspins
      INTEGER, INTENT (IN) :: lmd,ntype,irank,jspin,jspin_b
      INTEGER, INTENT (IN) :: nlod,llod,loplod,natd
      INTEGER, INTENT (IN) :: nntot,nntot2,ntyp
      REAL,    INTENT (IN) :: theta,phi,socscreen_fac
      LOGICAL, INTENT (IN) :: l_socscreen
C     ..
C     .. Array Arguments ..
      COMPLEX, INTENT (INOUT)::
     >        tdulo(0:lmd,nlod,-llod:llod),
     >        tuulo(0:lmd,nlod,-llod:llod),
     >        tulou(0:lmd,nlod,-llod:llod),
     >        tulod(0:lmd,nlod,-llod:llod),
     >        tuloulo(nlod,-llod:llod,nlod,-llod:llod)
      COMPLEX, INTENT (IN) :: angso(0:lmaxd,0:lmaxd)
      COMPLEX, INTENT (IN) :: yl1((lmaxd+1)**2),yl2((lmaxd+1)**2)
      INTEGER, INTENT (IN) :: llo(nlod,ntypd),nlo(ntypd)
      INTEGER, INTENT (IN) :: lo1l(0:llod,ntypd)
      INTEGER, INTENT (IN) :: jri(ntypd),lmax(ntypd)
      INTEGER, INTENT (IN) :: ulo_der(nlod,ntypd)
      REAL,    INTENT (IN) :: dx(ntypd),rmsh(jmtd,ntypd)
      REAL,    INTENT (IN) :: ello(nlod,ntypd,max(jspd,2))
      REAL,    INTENT (IN) :: epar(0:lmaxd,ntypd,max(jspd,2))
      REAL,    INTENT (IN) :: vr(jmtd,ntypd,jspd)
      REAL,    INTENT (IN) :: flo  (ntypd,jmtd,2,nlod),
     >                        flo_b(ntypd,jmtd,2,nlod)
      REAL,    INTENT (IN) :: p(jmtd,0:lmaxd,0:lmaxd),
     >                        q(jmtd,0:lmaxd,0:lmaxd),
     >                        p_b(jmtd,0:lmaxd,0:lmaxd),
     >                        q_b(jmtd,0:lmaxd,0:lmaxd)
      REAL,    INTENT (IN) :: kdiff (3,nntot)
      REAL,    INTENT (IN) :: kdiff2(3,nntot2)
      REAL,    INTENT (IN) :: jj1(0:lmaxd,jmtd),jj2(0:lmaxd,jmtd)
      LOGICAL, INTENT (IN) :: l_dulo(nlod,ntypd),l_spav
C     ..
C     .. Local Scalars ..
      REAL gc1,gc2,c,e,r0
      COMPLEX cil,ci
      INTEGER i,l,l2,lh,lm,lmp,lp,lo,lop
      INTEGER lp1,lpl,m,mp,mu
      INTEGER ljj1,ljj2,mjj1,mjj2
      INTEGER ltil1,mtil1,mtil2
      INTEGER lmini,lmini2, lmaxi,lmaxi2
      INTEGER ll,llp,lljj1,lljj2,lmjj1,lmjj2
      INTEGER lltil1,lmtil1,lmtil2,sp1,sp2
      INTEGER mmin,mmax,lwn, ikpt_b,ikpt_b2
C     ..
C     .. Local Arrays ..
      REAL, ALLOCATABLE :: dvulo(:,:,:,:,:),ulovd(:,:,:,:,:)
      REAL, ALLOCATABLE :: uvulo(:,:,:,:,:),ulovu(:,:,:,:,:)
      REAL, ALLOCATABLE :: ulovulo(:,:,:,:,:)
      REAL, ALLOCATABLE :: plo(:,:,:),plo_b(:,:,:)
      REAL, ALLOCATABLE :: x(:)
      REAL, ALLOCATABLE :: v0(:),vso(:,:)
      INTEGER :: ispjsp(2)
      DATA ispjsp/1,-1/

#if (defined(CPP_MPI) && !defined(CPP_T90))
      INCLUDE 'mpif.h'
      INTEGER ierr(3)
#endif

      ci = cmplx(0.0,1.0)
      c = c_light(1.0)
      sp1 = ispjsp(jspin)
      sp2 = ispjsp(jspin_b)

      allocate( ulovulo(0:lmaxd,0:lmaxd,0:lmaxd,nlod,nlod) )
      allocate( ulovu(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( ulovd(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( uvulo(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( dvulo(0:lmaxd,0:lmaxd,0:lmaxd,0:lmaxd,nlod) )
      allocate( x(jmtd) )
      allocate( v0(jmtd), vso(jmtd,2) )
      allocate( plo(jmtd,0:lmaxd,nlod), plo_b(jmtd,0:lmaxd,nlod) )

      lwn = lmax(ntyp)
      r0 = rmsh(1,ntyp)

       DO lo = 1, nlo(ntyp)
        DO ljj1 = 0, lwn
         DO i=1,jri(ntyp)
          plo(i,ljj1,lo) = flo(ntyp,i,1,lo)*jj1(ljj1,i)
          plo_b(i,ljj1,lo) = flo_b(ntyp,i,1,lo)*jj2(ljj1,i)
         ENDDO
        ENDDO
       ENDDO

c****************************************************c
c   compute radial integrals                         c
c     <u(l')jj1(ljj1) | vso(lh) | u(l)jj2(ljj2)>     c
c****************************************************c

! lo-lo
       DO lop = 1, nlo(ntyp)
        lp = llo(lop,ntyp)
        DO lo = 1, nlo(ntyp)
         l = llo(lop,ntyp)

        ! construct vso
        v0 = 0.
        IF (jspins.EQ.1) THEN
           v0(1:jri(ntyp)) = vr(1:jri(ntyp),ntyp,1)
           e = ( ello(lo,ntyp,1) + ello(lop,ntyp,1) )/2.
        ELSE
           DO i = 1,jri(ntyp)
             v0(i) = (vr(i,ntyp,1)+vr(i,ntyp,jspins))/2.
           END DO
           e = ( ello(lo,ntyp,1) + ello(lo,ntyp,jspins)
     >          +ello(lop,ntyp,1) + ello(lop,ntyp,jspins) )/4.
        END IF

        CALL sointg(
     >         e,vr(:,ntyp,:),v0,r0,dx(ntyp),jri(ntyp),jmtd,c,jspins,
     <         vso)  

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
          DO ljj2 = 0, lwn
           DO lh=0,lwn
              lmini = abs(lp-ljj1)
              lmini2= abs(lh-ljj2)
              lmaxi = lp+ljj1
              lmaxi2= lh+ljj2
              if((lmini.gt.lh) .or. (lmaxi.lt.lh).or.
     >           (lmini2.gt.l) .or. (lmaxi2.lt.l).or.
     >           (mod(lp+lh+ljj1,2).ne.0).or.
     >           (mod( l+lh+ljj2,2).ne.0) ) then

               ulovulo(lh,ljj2,ljj1,lo,lop) = 0.

              else

               DO i=1,jri(ntyp)
                x(i) = plo(i,ljj1,lop)*plo_b(i,ljj2,lo)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(ntyp)),r0,dx(ntyp),jri(ntyp),
     >                     ulovulo(lh,ljj2,ljj1,lo,lop))

              endif
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO!lh

! apw-lo
       DO lop = 1, nlo(ntyp)
        lp = llo(lop,ntyp)
        DO l = 0, lwn

        ! construct vso
        v0 = 0.
        IF (jspins.EQ.1) THEN
           v0(1:jri(ntyp)) = vr(1:jri(ntyp),ntyp,1)
           e = ( epar(l,ntyp,1) + ello(lop,ntyp,1) )/2.
        ELSE
           DO i = 1,jri(ntyp)
             v0(i) = (vr(i,ntyp,1)+vr(i,ntyp,jspins))/2.
           END DO
           e = ( epar(l,ntyp,1) + epar(l,ntyp,jspins)
     >          +ello(lop,ntyp,1) + ello(lop,ntyp,jspins) )/4.
        END IF

        CALL sointg(
     >         e,vr(:,ntyp,:),v0,r0,dx(ntyp),jri(ntyp),jmtd,c,jspins,
     <         vso)  

        ! spin-averaging
        if(l_spav)then
            DO i= 1,jmtd
             vso(i,1)= (vso(i,1)+vso(i,2))/2.
             vso(i,2)= vso(i,1)
           ENDDO        
        endif

        ! calculate radial integrals
         DO ljj1 = 0, lwn
          DO ljj2 = 0, lwn
           DO lh=0,lwn
              lmini = abs(lp-ljj1)
              lmini2= abs(lh-ljj2)
              lmaxi = lp+ljj1
              lmaxi2= lh+ljj2
              if((lmini.gt.lh) .or. (lmaxi.lt.lh).or.
     >           (lmini2.gt.l) .or. (lmaxi2.lt.l).or.
     >           (mod(lp+lh+ljj1,2).ne.0).or.
     >           (mod( l+lh+ljj2,2).ne.0) ) then

               ulovu(lh,l,ljj2,ljj1,lop) = 0.
               ulovd(lh,l,ljj2,ljj1,lop) = 0.

              else

               DO i=1,jri(ntyp)
                x(i) = plo(i,ljj1,lop)*p_b(i,l,ljj2)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(ntyp)),r0,dx(ntyp),jri(ntyp),
     >                     ulovu(lh,l,ljj2,ljj1,lop))
               DO i=1,jri(ntyp)
                x(i) = plo(i,ljj1,lop)*q_b(i,l,ljj2)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(ntyp)),r0,dx(ntyp),jri(ntyp),
     >                     ulovd(lh,l,ljj2,ljj1,lop))

              endif

              lmini = abs(l-ljj1)
              lmini2= abs(lh-ljj2)
              lmaxi = l+ljj1
              lmaxi2= lh+ljj2
              if((lmini.gt.lh) .or. (lmaxi.lt.lh).or.
     >           (lmini2.gt.lp) .or. (lmaxi2.lt.lp).or.
     >           (mod( l+lh+ljj1,2).ne.0).or.
     >           (mod(lp+lh+ljj2,2).ne.0) ) then

               uvulo(lh,l,ljj2,ljj1,lop) = 0.
               dvulo(lh,l,ljj2,ljj1,lop) = 0.

              else

               DO i=1,jri(ntyp)
                x(i) = p(i,l,ljj1)*plo_b(i,ljj2,lop)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(ntyp)),r0,dx(ntyp),jri(ntyp),
     >                     uvulo(lh,l,ljj2,ljj1,lop))
               DO i=1,jri(ntyp)
                x(i) = q(i,l,ljj1)*plo_b(i,ljj2,lop)*vso(i,jspin)
               ENDDO
               call intgr0(x(1:jri(ntyp)),r0,dx(ntyp),jri(ntyp),
     >                     dvulo(lh,l,ljj2,ljj1,lop))
             
              endif
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO!lh

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

! lo-lo
       DO lop=1,nlo(ntyp)
        lp = llo(lop,ntyp)
        DO mp=-lp,lp
         DO ltil1=0,lwn
          lltil1 = ltil1*(ltil1+1)
          DO mtil1=-ltil1,ltil1
           lmtil1 = lltil1+mtil1
           mmin = max(-ltil1,mtil1-1)
           mmax = min( ltil1,mtil1+1)
           mjj1=mtil1-mp
           DO ljj1=0,lwn
            lljj1 = ljj1 * (ljj1 + 1)
            lmjj1 = lljj1 + mjj1
            lmini = abs(lp-ljj1)
            lmaxi = lp+ljj1
            ! Gaunt conditions (G1):
            if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >         (mod(lp+ltil1+ljj1,2).ne.0).or.
     >         (abs(mjj1).gt.ljj1)) cycle

             gc1 = gaunt1(ltil1, ljj1, lp, 
     >                    mtil1, mjj1, mp, lmaxd)

        DO lo=1,nlo(ntyp)
         l = llo(lo,ntyp)
         DO m=-l,l
          DO mtil2=mmin,mmax     ! mtil1-1 <= mtil2 <= mtil1+1
           lmtil2 = lltil1+mtil2
           mjj2=m-mtil2
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
     >                     m, mjj2, mtil2, lmaxd)

             ! set up prefactor
             cil =   ( ci ** (lp - l + ljj1 + ljj2) )
     +              *( (-1.0) ** ljj1 )
     +              * gc1 * gc2
     +              * conjg( yl1(lmjj1 + 1) )
     +              * conjg( yl2(lmjj2 + 1) )
     +              * conjg( angso(lmtil1,lmtil2) )

             ! compute T-coefficients
             tuloulo(lo, m, lop, mp) = tuloulo(lo, m, lop, mp)
     >         + cil * ulovulo(ltil1, ljj2, ljj1, lo, lop)

           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
        ENDDO!mtil2
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1

! apw-lo
       DO lop=1,nlo(ntyp)
        lp = llo(lop,ntyp)
        DO mp=-lp,lp
         DO ltil1=0,lwn
          lltil1 = ltil1*(ltil1+1)
          DO mtil1=-ltil1,ltil1
           lmtil1 = lltil1+mtil1
           mmin = max(-ltil1,mtil1-1)
           mmax = min( ltil1,mtil1+1)
           mjj1=mtil1-mp
           DO ljj1=0,lwn
            lljj1 = ljj1 * (ljj1 + 1)
            lmjj1 = lljj1 + mjj1
            lmini = abs(lp-ljj1)
            lmaxi = lp+ljj1
            ! Gaunt conditions (G1):
            if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >         (mod(lp+ltil1+ljj1,2).ne.0).or.
     >         (abs(mjj1).gt.ljj1)) cycle

             gc1 = gaunt1(ltil1, ljj1, lp, 
     >                    mtil1, mjj1, mp, lmaxd)

        DO mtil2=mmin,mmax     ! mtil1-1 <= mtil2 <= mtil1+1
         lmtil2 = lltil1+mtil2
         DO l=0,lwn
          ll = l*(l+1)
          DO m=-l,l
           lm = ll+m
           mjj2=m-mtil2
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
     >                     m, mjj2, mtil2, lmaxd)

             ! set up prefactor
             cil =   ( ci ** (lp - l + ljj1 + ljj2) )
     +              *( (-1.0) ** ljj1 )
     +              * gc1 * gc2
     +              * conjg( yl1(lmjj1 + 1) )
     +              * conjg( yl2(lmjj2 + 1) )
     +              * conjg( angso(lmtil1,lmtil2) )

             ! compute T-coefficients
             tulou(lm, lop, mp) = tulou(lm, lop, mp)
     >         + cil * ulovu(ltil1, l, ljj2, ljj1, lop)

             tulod(lm, lop, mp) = tulod(lm, lop, mp)
     >         + cil * ulovd(ltil1, l, ljj2, ljj1, lop)

           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
        ENDDO!mtil2
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1


       DO l=0,lwn
        ll = l*(l+1)
        DO m=-l,l
         lm = ll+m
         DO ltil1=0,lwn
          lltil1 = ltil1*(ltil1+1)
          DO mtil1=-ltil1,ltil1
           lmtil1 = lltil1+mtil1
           mmin = max(-ltil1,mtil1-1)
           mmax = min( ltil1,mtil1+1)
           mjj1=mtil1-m
           DO ljj1=0,lwn
            lljj1 = ljj1 * (ljj1 + 1)
            lmjj1 = lljj1 + mjj1
            lmini = abs(l-ljj1)
            lmaxi = l+ljj1
            ! Gaunt conditions (G1):
            if((lmini.gt.ltil1).or.(lmaxi.lt.ltil1).or.
     >         (mod(l+ltil1+ljj1,2).ne.0).or.
     >         (abs(mjj1).gt.ljj1)) cycle

             gc1 = gaunt1(ltil1, ljj1, l, 
     >                    mtil1, mjj1, m, lmaxd)

        DO mtil2=mmin,mmax     ! mtil1-1 <= mtil2 <= mtil1+1
         lmtil2 = lltil1+mtil2
         DO lop=1,nlo(ntyp)
          lp = llo(lop,ntyp)
          DO mp=-lp,lp
           mjj2=mp-mtil2
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
     >                     mp, mjj2, mtil2, lmaxd)

             ! set up prefactor
             cil =   ( ci ** (l - lp + ljj1 + ljj2) )
     +              *( (-1.0) ** ljj1 )
     +              * gc1 * gc2
     +              * conjg( yl1(lmjj1 + 1) )
     +              * conjg( yl2(lmjj2 + 1) )
     +              * conjg( angso(lmtil1,lmtil2) )

             ! compute T-coefficients
             tuulo(lm, lop, mp) = tuulo(lm, lop, mp)
     >         + cil * uvulo(ltil1, l, ljj2, ljj1, lop)

             tdulo(lm, lop, mp) = tdulo(lm, lop, mp)
     >         + cil * dvulo(ltil1, l, ljj2, ljj1, lop)

           ENDDO!ljj2
          ENDDO!m
         ENDDO!l
        ENDDO!mtil2
           ENDDO!ljj1
          ENDDO!mp
         ENDDO!lp
        ENDDO!mtil1
       ENDDO!ltil1

#if (defined(CPP_MPI) && !defined(CPP_T90))
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

      deallocate( ulovulo )
      deallocate( uvulo, dvulo, ulovu, ulovd )
      deallocate( plo,plo_b )
      deallocate( x )
      deallocate( v0, vso )

      END SUBROUTINE wann_uHu_soc_tlo
      END MODULE m_wann_uHu_soc_tlo
