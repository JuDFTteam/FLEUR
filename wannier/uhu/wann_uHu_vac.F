c***************************************c
c   Vacuum contribution to uHu matrix   c
c   in FLAPW film mode                  c
c***************************************c
c   < u_{k+b1} | H_{k} | u_{k+b2} >     c
c***************************************c
c              J.-P. Hanke, Dec. 2015   c
c***************************************c
      MODULE m_wann_uHu_vac
      USE m_juDFT

      CONTAINS
      SUBROUTINE wann_uHu_vac(
     >     chi,l_noco,l_soc,zrfs,jspins,nlotot,qss,nbnd,z1,nmzxyd,
     >     nmzd,n2d,nv2d,k1d,k2d,k3d,n3d,nvac,ig,rgphs,nmzxy,
     >     nmz,delz,ig2,nq2,kv2,area,bmat,bbmat,evac,evac_b,
     >     bkpt,bkpt_b,vzxy,vz,nslibd,nslibd_b,jspin,jspin_b,
     >     ico,k1,k2,k3,k1_b,k2_b,k3_b,jspd,nvd,nbasfcn,neigd,
     >     zMat,zMat_b,nv,nv_b,omtil,gb,gb2,sign2,uHu)
#include "cpp_double.h"
      USE m_types
      use m_constants, only : pimach
      use m_intgr, only : intgz0
      use m_d2fdz2cmplx

      implicit none

      TYPE(t_zmat), INTENT(IN) :: zMat, zMat_b

c     .. scalar Arguments..
      logical, intent (in) :: l_noco,l_soc,zrfs
      integer, intent (in) :: nlotot,jspin_b,n2d,jspins,ico
      integer, intent (in) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n3d,nbnd
      integer, intent (in) :: nmzxy,nmz,nslibd,nslibd_b,nvac
      integer, intent (in) :: jspin,jspd,nvd,nq2
      integer, intent (in) :: nbasfcn,neigd
      real,    intent (in) :: delz,z1,omtil,area
      complex, intent (in) :: chi

c     ..array arguments..
      real,    intent (in) :: bkpt(3),bkpt_b(3),evac(2),evac_b(2)
      real,    intent (in) :: qss(3)
      real,    intent (in) :: vz(nmzd,2,4)
      real,    intent (in) :: bbmat(3,3),bmat(3,3)
      real,    intent (in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig2(n3d),nv(jspd),nv_b(jspd)
      integer, intent (in) :: kv2(2,n2d)
      integer, intent (in) :: gb(3),gb2(3)
      integer, intent (in) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      integer, intent (in) :: k1_b(nvd,jspd),k2_b(nvd,jspd)
      integer, intent (in) :: k3_b(nvd,jspd)
      complex, intent (in) :: vzxy(nmzxyd,n2d-1,2)
      complex, intent (inout) :: uHu(nbnd,nbnd)

c     ..local arrays..
      complex, allocatable :: ac(:,:),bc(:,:),ac_b(:,:),bc_b(:,:)
      complex, allocatable :: vv(:),fac(:),fac_b(:),fac1(:)
      complex, allocatable :: du(:),due(:),du_b(:),due_b(:)
      complex, allocatable :: vxy(:,:),vxy_help(:)
      real,    allocatable :: dt(:),dte(:)
      real,    allocatable :: t(:),te(:),tei(:)
      real,    allocatable :: u(:,:),ue(:,:),v(:)
      real,    allocatable :: dt_b(:),dte_b(:)
      real,    allocatable :: t_b(:),te_b(:),tei_b(:)
      real,    allocatable :: u_b(:,:),ue_b(:,:)
      real,    allocatable :: zmsh(:)

c     ..local scalars..
      logical tail
      real wronk,arg,zks,zks_b,tpi,vz0(2),vz0_b(2),scale,evacp,ev,const
      real xx(nmz),xximag(nmz),phase,phase2,rk,rk_b
      real :: qss1,qss2,xv,yv,s(3)
      integer i,m,l,j,k,n,nv2,nv2_b,ivac,n2,n2_b,sign,ik,sign2
      integer :: lp,np0,np1,addnoco,addnoco2,jspin2,jspin2_b,jspin2H
      complex :: av,bv,ic,c_1
c      complex :: tuu,tud,tdu,tdd
      complex, allocatable :: tuu(:,:),tud(:,:),tdu(:,:),tdd(:,:)
      complex, allocatable :: mat(:,:)
      integer, allocatable :: kvac1(:),kvac2(:),map2(:)
      integer, allocatable :: kvac1_b(:),kvac2_b(:),map2_b(:)
      integer symvac,symvacvac,igvm2,igvm2i
      integer :: j1,j2,i1,i2,ind3,ind2
      character(len=20) :: fname
c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt

c..   external functions
      external vacudz,vacuz

      allocate ( ac(nv2d,nslibd),bc(nv2d,nslibd),
     +           ac_b(nv2d,nslibd_b),bc_b(nv2d,nslibd_b),
     +           dt(nv2d),dte(nv2d),t(nv2d),te(nv2d),
     +           tei(nv2d),u(nmzd,nv2d),ue(nmzd,nv2d),
     +           dt_b(nv2d),dte_b(nv2d),t_b(nv2d),te_b(nv2d),
     +           tei_b(nv2d),u_b(nmzd,nv2d),ue_b(nmzd,nv2d),
     +           v(3),kvac1(nv2d),kvac2(nv2d),map2(nvd),
     +           kvac1_b(nv2d),kvac2_b(nv2d),map2_b(nvd),
     +           vv(nmzd),zmsh(nmzd),vxy(nmzxyd,n2d-1),vxy_help(n2d-1),
     +           fac1(nmzd), fac(nmzd), fac_b(nmzd),
     +           du(nmzd), due(nmzd), du_b(nmzd), due_b(nmzd),
     +           tuu(nv2d,nv2d),tud(nv2d,nv2d),tdu(nv2d,nv2d),
     +           tdd(nv2d,nv2d),mat(nslibd_b,nv2d) )

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)
      tail = .true.
      np0 = nmzxy + 1
      np1 = nmz + 1

      jspin2 = jspin
      jspin2_b = jspin_b
      jspin2H = ico
      if(l_soc.and.jspins.eq.1) then
         jspin2 = 1
         jspin2_b = 1
         jspin2H = 1
      endif

c.. determining the indexing array (in-plane stars)
c.. for the k-point

      wronk = 2.0
      const = 1.0 / ( sqrt(omtil)*wronk )

      do ivac = 1,2
         vz0(ivac) = vz(nmz,ivac,jspin2)
         vz0_b(ivac) = vz(nmz,ivac,jspin2_b)
c         write(*,*)'vz0',ivac,vz0(ivac),vz0_b(ivac)
      enddo

      n2 = 0 ; n2_b = 0

      addnoco=0
      addnoco2=0
      if(l_noco.and.jspin.eq.2)then
        addnoco=nv(1)+nlotot
      endif
      if(l_noco.and.jspin_b.eq.2)then
        addnoco2=nv_b(1)+nlotot
      endif

      do 40 k = 1,nv(jspin)
         do 30 j = 1,n2
            if ( k1(k,jspin).eq.kvac1(j) .and.
     +          k2(k,jspin).eq.kvac2(j) ) then
                map2(k) = j
                goto 40
             endif 
 30      continue
         n2 = n2 + 1
         
         IF(n2>nv2d) then
            write(*,*)n2,nv2d,'jspin',jspin
         endif
 
         IF (n2>nv2d)  CALL juDFT_error("wannier uHu vac",calledby
     +        ="wann_uHu_vac")

         kvac1(n2) = k1(k,jspin)
         kvac2(n2) = k2(k,jspin)
         map2(k) = n2
 40   continue


c.. and for the b-point
 
      do 41 k = 1,nv_b(jspin_b)
         do 31 j = 1,n2_b
            if ( k1_b(k,jspin_b).eq.kvac1_b(j) .and.
     +          k2_b(k,jspin_b).eq.kvac2_b(j) ) then
                map2_b(k) = j
                goto 41
             endif
 31      continue
         n2_b = n2_b + 1
     
         IF(n2_b>nv2d) then
            write(*,*)n2_b,nv2d,'jspin_b',jspin_b
         endif

         IF (n2_b>nv2d)  CALL juDFT_error("wannier uHu vac",calledby
     +        ="wann_uHu_vac")

         kvac1_b(n2_b) = k1_b(k,jspin_b)
         kvac2_b(n2_b) = k2_b(k,jspin_b)
         map2_b(k) = n2_b
 41   continue




c...cycle by the vacua
      do 140 ivac = 1,2!nvac

      sign = 3. - 2.*ivac
      evacp = evac(ivac)

      vxy(:,:) = vzxy(:,:,ivac)
      IF(l_noco) THEN ! symmetrization
       IF (nvac.EQ.1 .and. ivac.EQ.2 .AND.(.NOT.zrfs) ) THEN
        DO i=1,nmzxy
         DO igvm2 = 2,nq2
          igvm2i = ig2(ig(-kv2(1,igvm2),-kv2(2,igvm2),0))
          vxy_help(igvm2-1) = vxy(i,igvm2i-1)
         ENDDO
         DO igvm2 = 2,nq2
          vxy(i,igvm2-1) = vxy_help(igvm2-1)
         ENDDO
        ENDDO
       ENDIF
      ENDIF

      tuu = cmplx(0.,0.)
      tud = cmplx(0.,0.)
      tdu = cmplx(0.,0.)
      tdd = cmplx(0.,0.)

      nv2 = n2 ; nv2_b = n2_b

c.. the body of the routine

      qss1=0.0
      qss2=0.0
      if(l_noco.and.jspin.eq.1)then
        qss1=-qss(1)/2.0
        qss2=-qss(2)/2.0
      elseif(l_noco.and.jspin.eq.2)then
        qss1=qss(1)/2.0
        qss2=qss(2)/2.0
      endif

      do ik = 1,nv2
         v(1) = bkpt(1) + kvac1(ik) + qss1
         v(2) = bkpt(2) + kvac2(ik) + qss2
         v(3) = 0.
         ev = evacp - 0.5*dot_product(v,matmul(bbmat,v))
!         ev = evacp - 0.5*dotirp(v,v,bbmat)
         call vacuz(ev,vz(1,ivac,jspin2),vz0(ivac),nmz,delz,t(ik),
     +        dt(ik),u(1,ik))
         call vacudz(ev,vz(1,ivac,jspin2),vz0(ivac),nmz,delz,te(ik),
     +        dte(ik),tei(ik),ue(1,ik),dt(ik),
     +        u(1,ik))
         scale = wronk/ (te(ik)*dt(ik)-dte(ik)*t(ik))
         te(ik) = scale*te(ik)
         dte(ik) = scale*dte(ik)
         tei(ik) = scale*tei(ik)
         do j = 1,nmz
            ue(j,ik) = scale*ue(j,ik)
         enddo
      enddo
c-----> construct a and b coefficients for the k-point
       symvacvac=1
c       if (nvac==1) symvacvac=2
       do symvac=1,symvacvac
         do i = 1,nv2d
            do n = 1,nslibd
               ac(i,n) = cmplx(0.0,0.0)
               bc(i,n) = cmplx(0.0,0.0)
            enddo   
            do n = 1,nslibd_b
               ac_b(i,n) = cmplx(0.0,0.0)
               bc_b(i,n) = cmplx(0.0,0.0)
            enddo   
         enddo   

        if (symvac==2) sign=-1.0   

      do k = 1,nv(jspin)
         l = map2(k)
         zks = k3(k,jspin)*bmat(3,3)*sign
         arg = zks*z1
         c_1 = cmplx(cos(arg),sin(arg)) * const
         av = -c_1 * cmplx( dte(l),zks*te(l) )
         bv =  c_1 * cmplx(  dt(l),zks* t(l) )
c-----> loop over basis functions

         IF(zMat%l_real) THEN
            do n = 1,nslibd
               ac(l,n) = ac(l,n) + zMat%z_r(k+addnoco,n)*av
               bc(l,n) = bc(l,n) + zMat%z_r(k+addnoco,n)*bv
            enddo
         ELSE
            do n = 1,nslibd
               ac(l,n) = ac(l,n) + zMat%z_c(k+addnoco,n)*av
               bc(l,n) = bc(l,n) + zMat%z_c(k+addnoco,n)*bv
            enddo
         END IF
      enddo


c...now for the k+b point

      evacp = evac_b(ivac)
      do ik = 1,nv2_b
         v(1) = bkpt_b(1) + kvac1_b(ik) + qss1
         v(2) = bkpt_b(2) + kvac2_b(ik) + qss2
         v(3) = 0.
         ev = evacp - 0.5*dot_product(v,matmul(bbmat,v))
!         ev = evacp - 0.5*dotirp(v,v,bbmat)
         call vacuz(ev,vz(1,ivac,jspin2_b),vz0_b(ivac),nmz,delz,
     +        t_b(ik),dt_b(ik),u_b(1,ik))
         call vacudz(ev,vz(1,ivac,jspin2_b),vz0_b(ivac),nmz,delz,
     +        te_b(ik),dte_b(ik),tei_b(ik),ue_b(1,ik),dt_b(ik),
     +        u_b(1,ik))
         scale = wronk/ (te_b(ik)*dt_b(ik)-dte_b(ik)*t_b(ik))
         te_b(ik) = scale*te_b(ik)
         dte_b(ik) = scale*dte_b(ik)
         tei_b(ik) = scale*tei_b(ik)
         do j = 1,nmz
            ue_b(j,ik) = scale*ue_b(j,ik)
         enddo
      enddo
c-----> construct a and b coefficients for the k+b point

      do k = 1,nv_b(jspin_b)
         l = map2_b(k)
         zks = k3_b(k,jspin_b)*bmat(3,3)*sign
         arg = zks*z1
         c_1 = cmplx(cos(arg),sin(arg)) * const
         av = -c_1 * cmplx( dte_b(l),zks*te_b(l) )
         bv =  c_1 * cmplx( dt_b(l),zks*t_b(l) )
c-----> loop over basis functions

         IF(zMat_b%l_real) THEN
            do n = 1,nslibd_b
               ac_b(l,n) = ac_b(l,n) + zMat_b%z_r(k+addnoco2,n)*av
               bc_b(l,n) = bc_b(l,n) + zMat_b%z_r(k+addnoco2,n)*bv
            enddo
         ELSE
            do n = 1,nslibd_b
               ac_b(l,n) = ac_b(l,n) + zMat_b%z_c(k+addnoco2,n)*av
               bc_b(l,n) = bc_b(l,n) + zMat_b%z_c(k+addnoco2,n)*bv
            enddo
         END IF
      enddo


      ! set up z-mesh and plane-wave factors
      zks  =  gb(3) *bmat(3,3)   ! TODO: different sign before
      zks_b= -gb2(3)*bmat(3,3)   ! TODO: different sign before
      do i=1,nmz
       zmsh(i)= ( z1 + (i-1)*delz )*sign
       fac(i)  = cmplx( cos(zmsh(i)*zks), sin(zmsh(i)*zks) )
       fac_b(i)= cmplx( cos(zmsh(i)*zks_b),sin(zmsh(i)*zks_b))
       fac1(i) = fac(i) * fac_b(i)
      enddo


      ! calculate uHu matrix elements
      do l = 1,nv2

       j1 = kvac1(l) - gb(1)
       j2 = kvac2(l) - gb(2)

       do lp = 1,nv2_b
        i1 = j1 - kvac1_b(lp) + gb2(1)
        i2 = j2 - kvac2_b(lp) + gb2(2)

c        i1 = -i1
c        i2 = -i2

        ind3 = ig(sign2*i1,sign2*i2,0)
        IF (ind3.EQ.0) CYCLE

        phase = rgphs(i1,i2,0)   ! TODO: sign also here?
        phase2= rgphs(sign2*i1,sign2*i2,0)
        if(phase.ne.phase2) then
           stop 'rgphs in wann_uHu_vac'
        endif

        ind2 = ig2(ind3)
        IF (ind2.EQ.0) CALL juDFT_error("error in map2 for 2-d stars",
     >                                calledby="wann_uHu_vac")
        ind2 = ind2 - 1
        IF (ind2.NE.0) THEN ! warping components G=/=0
                            ! need integral with vxy and plane-waves
           ! tuu
           DO i=1,nmzxy
              xx(np0-i) = u(i,l)*u_b(i,lp)
     >                  * real( fac1(i)*vxy(i,ind2) )
              xximag(np0-i) = u(i,l)*u_b(i,lp)
     >                      * aimag( fac1(i)*vxy(i,ind2) )
           ENDDO
           CALL intgz0(xx,delz,nmzxy,xv,tail)
           CALL intgz0(xximag,delz,nmzxy,yv,tail)
           tuu(lp,l) = phase*cmplx(xv,yv)

           ! tud
           DO i=1,nmzxy
              xx(np0-i) = u(i,l)*ue_b(i,lp)
     >                  * real( fac1(i)*vxy(i,ind2) )
              xximag(np0-i) = u(i,l)*ue_b(i,lp)
     >                      * aimag( fac1(i)*vxy(i,ind2) )
           ENDDO
           CALL intgz0(xx,delz,nmzxy,xv,tail)
           CALL intgz0(xximag,delz,nmzxy,yv,tail)
           tud(lp,l) = phase*cmplx(xv,yv)

           ! tdu
           DO i=1,nmzxy
              xx(np0-i) = ue(i,l)*u_b(i,lp) 
     >                  * real( fac1(i)*vxy(i,ind2) )
              xximag(np0-i) = ue(i,l)*u_b(i,lp)
     >                      * aimag( fac1(i)*vxy(i,ind2) )
           ENDDO
           CALL intgz0(xx,delz,nmzxy,xv,tail)
           CALL intgz0(xximag,delz,nmzxy,yv,tail)
           tdu(lp,l) = phase*cmplx(xv,yv)

           ! tdd
           DO i=1,nmzxy
              xx(np0-i) = ue(i,l)*ue_b(i,lp)
     >                  * real( fac1(i)*vxy(i,ind2) )
              xximag(np0-i) = ue(i,l)*ue_b(i,lp)
     >                      * aimag( fac1(i)*vxy(i,ind2) )
           ENDDO
           CALL intgz0(xx,delz,nmzxy,xv,tail)
           CALL intgz0(xximag,delz,nmzxy,yv,tail)
           tdd(lp,l) = phase*cmplx(xv,yv)

        ELSE ! non-warping components G==0
             ! need integral with H(z) = -1/2 d2/dz2 + vz + (G+k)^2

           IF ( (ico.EQ.1) .OR. (ico.EQ.2) ) THEN ! spin-diagonal

              ! determine second derivative of (u*fac) etc.
              CALL d2fdz2cmplx(nmzd,nmz,zmsh,delz,
     >                         u(:,l),fac,du)
              CALL d2fdz2cmplx(nmzd,nmz,zmsh,delz,
     >                         ue(:,l),fac,due)
              CALL d2fdz2cmplx(nmzd,nmz,zmsh,delz,
     >                         u_b(:,lp),fac_b,du_b)
              CALL d2fdz2cmplx(nmzd,nmz,zmsh,delz,
     >                         ue_b(:,lp),fac_b,due_b)

              ! determine |G+k+b1|^2 and |G'+k+b2|^2
              s(1) = bkpt(1) + kvac1(l)
              s(2) = bkpt(2) + kvac2(l)
              s(3) = 0.0
              rk = dot_product(s,matmul(bbmat,s))
!              rk = dotirp(s,s,bbmat)

              s(1) = bkpt_b(1) + kvac1_b(lp)
              s(2) = bkpt_b(2) + kvac2_b(lp)
              s(3) = 0.0 
              rk_b = dot_product(s,matmul(bbmat,s))
!              rk_b = dotirp(s,s,bbmat)

c              ! no kinetic energy if jspin =/= jspin_b
c              if(jspin.ne.jspin_b) then
c               rk = 0.0
c               rk_b = 0.0
c               du=cmplx(0.0,0.0)
c               due=cmplx(0.0,0.0)
c               du_b=cmplx(0.0,0.0)
c               due_b=cmplx(0.0,0.0)
c              endif

c              rk_b=0.0
c              du_b=cmplx(0.,0.)
c              due_b=cmplx(0.,0.)

              ! tuu
              DO i=1,nmz
               xx(np1-i) = u(i,l)*u_b(i,lp) 
     >           * real(fac1(i)) * ( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*real(  fac(i)*u(i,l)*du_b(i) 
     >                       + du(i)*u_b(i,lp)*fac_b(i) )

               xximag(np1-i) = u(i,l)*u_b(i,lp) 
     >           * aimag(fac1(i)) *( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*aimag(  fac(i)*u(i,l)*du_b(i) 
     >                        + du(i)*u_b(i,lp)*fac_b(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tuu(lp,l) = cmplx(xv,yv)

              ! tud
              DO i=1,nmz
               xx(np1-i) = u(i,l)*ue_b(i,lp) 
     >           * real(fac1(i)) * ( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*real(  fac(i)*u(i,l)*due_b(i) 
     >                       + du(i)*ue_b(i,lp)*fac_b(i) )

               xximag(np1-i) = u(i,l)*ue_b(i,lp) 
     >           * aimag(fac1(i)) *( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*aimag(  fac(i)*u(i,l)*due_b(i) 
     >                        + du(i)*ue_b(i,lp)*fac_b(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tud(lp,l) = cmplx(xv,yv)

              ! tdu
              DO i=1,nmz
               xx(np1-i) = ue(i,l)*u_b(i,lp) 
     >           * real(fac1(i)) * ( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*real(  fac(i)*ue(i,l)*du_b(i) 
     >                       + due(i)*u_b(i,lp)*fac_b(i) )

               xximag(np1-i) = ue(i,l)*u_b(i,lp) 
     >           * aimag(fac1(i)) *( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*aimag(  fac(i)*ue(i,l)*du_b(i) 
     >                        + due(i)*u_b(i,lp)*fac_b(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tdu(lp,l) = cmplx(xv,yv)

              ! tdd
              DO i=1,nmz
               xx(np1-i) = ue(i,l)*ue_b(i,lp) 
     >           * real(fac1(i)) * ( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*real(  fac(i)*ue(i,l)*due_b(i) 
     >                       + due(i)*ue_b(i,lp)*fac_b(i) )

               xximag(np1-i) = ue(i,l)*ue_b(i,lp) 
     >           * aimag(fac1(i)) *( vz(i,ivac,jspin2H)+0.25*(rk+rk_b) )
     >           - 0.25*aimag(  fac(i)*ue(i,l)*due_b(i) 
     >                        + due(i)*ue_b(i,lp)*fac_b(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tdd(lp,l) = cmplx(xv,yv)

           ELSE ! spin-off-diagonal

              DO i=1,nmz
                 vv(i) = cmplx( vz(i,ivac,3), vz(i,ivac,4) )
              ENDDO
              IF ( ico.EQ.4 ) vv = conjg(vv)

              ! tuu
              DO i=1,nmz
                 xx(np1-i) = u(i,l)*u_b(i,lp)
     >                     * real( fac1(i)*vv(i) )
                 xximag(np1-i) = u(i,l)*u_b(i,lp)
     >                         * aimag( fac1(i)*vv(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tuu(lp,l) = cmplx(xv,yv)

              ! tud
              DO i=1,nmz
                 xx(np1-i) = u(i,l)*ue_b(i,lp)
     >                     * real( fac1(i)*vv(i) )
                 xximag(np1-i) = u(i,l)*ue_b(i,lp)
     >                         * aimag( fac1(i)*vv(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tud(lp,l) = cmplx(xv,yv)

              ! tdu
              DO i=1,nmz
                 xx(np1-i) = ue(i,l)*u_b(i,lp)
     >                     * real( fac1(i)*vv(i) )
                 xximag(np1-i) = ue(i,l)*u_b(i,lp)
     >                         * aimag( fac1(i)*vv(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tdu(lp,l) = cmplx(xv,yv)

              ! tdd
              DO i=1,nmz
                 xx(np1-i) = ue(i,l)*ue_b(i,lp)
     >                     * real( fac1(i)*vv(i) )
                 xximag(np1-i) = ue(i,l)*ue_b(i,lp)
     >                         * aimag( fac1(i)*vv(i) )
              ENDDO
              CALL intgz0(xx,delz,nmz,xv,tail)
              CALL intgz0(xximag,delz,nmz,yv,tail)
              tdd(lp,l) = cmplx(xv,yv)

           ENDIF !(ico.eq.1).or.(ico.eq.2)

        ENDIF !(ind2.ne.0)


        ! contraction of integrals with a,b coefficients
        ! yields contribution to uHu matrix        
c        do i = 1,nslibd
c           do j = 1,nslibd_b
c              uHu(i,j) = uHu(i,j) + area*
c     *            ( ac(l,i)*conjg(ac_b(lp,j))*tuu +
c     +              ac(l,i)*conjg(bc_b(lp,j))*tud +
c     *              bc(l,i)*conjg(ac_b(lp,j))*tdu +
c     +              bc(l,i)*conjg(bc_b(lp,j))*tdd )
c           enddo
c        enddo


       enddo ! lp
      enddo ! l

      if(.false.) then
      call CPP_BLAS_cgemm('C','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    ac_b(1,1),nv2d,tuu(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','T',nslibd,nslibd_b,nv2d,chi,
     >                    ac(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)

      call CPP_BLAS_cgemm('C','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    bc_b(1,1),nv2d,tud(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','T',nslibd,nslibd_b,nv2d,chi,
     >                    ac(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)

      call CPP_BLAS_cgemm('C','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    ac_b(1,1),nv2d,tdu(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','T',nslibd,nslibd_b,nv2d,chi,
     >                    bc(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)

      call CPP_BLAS_cgemm('C','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    bc_b(1,1),nv2d,tdd(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','T',nslibd,nslibd_b,nv2d,chi,
     >                    bc(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)
      endif

      if(.true.)then ! standard so far
      call CPP_BLAS_cgemm('T','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    ac_b(1,1),nv2d,tuu(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','C',nslibd,nslibd_b,nv2d,chi,
     >                    ac(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)

      call CPP_BLAS_cgemm('T','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    bc_b(1,1),nv2d,tud(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','C',nslibd,nslibd_b,nv2d,chi,
     >                    ac(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)

      call CPP_BLAS_cgemm('T','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    ac_b(1,1),nv2d,tdu(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','C',nslibd,nslibd_b,nv2d,chi,
     >                    bc(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)

      call CPP_BLAS_cgemm('T','N',nslibd_b,nv2d,nv2d,cmplx(area),
     >                    bc_b(1,1),nv2d,tdd(1,1),nv2d,cmplx(0.0),
     >                    mat(1,1),nslibd_b)
      call CPP_BLAS_cgemm('T','C',nslibd,nslibd_b,nv2d,chi,
     >                    bc(1,1),nv2d,mat(1,1),nslibd_b,cmplx(1.0),
     >                    uHu,nbnd)
      endif

      enddo !symvac
c... cycle by the vacua finishes
 140  enddo      

      deallocate ( ac,bc,dt,dte,t,te,tei,u,ue,
     +             v,kvac1,kvac2,map2 )
      deallocate ( ac_b,bc_b,dt_b,dte_b,t_b,te_b,tei_b,u_b,ue_b,
     +             kvac1_b,kvac2_b,map2_b )
      deallocate ( vxy, vxy_help, vv, zmsh, fac1, fac, fac_b )
      deallocate ( du,due,du_b,due_b )
      deallocate ( tuu,tud,tdu,tdd,mat )

c      call fleur_end("stop")

      END SUBROUTINE wann_uHu_vac
      END MODULE m_wann_uHu_vac
