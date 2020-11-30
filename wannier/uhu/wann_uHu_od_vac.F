c***************************************c
c   Vacuum contribution to uHu matrix   c
c   in one-dimensional FLAPW mode       c
c***************************************c
c   < u_{k+b1} | H_{k} | u_{k+b2} >     c
c***************************************c
c              J.-P. Hanke, Dec. 2015   c
c***************************************c
      MODULE m_wann_uHu_od_vac
      CONTAINS
      SUBROUTINE wann_uHu_od_vac(
     >     oneD,vacuum,stars,cell,
     >     chi,l_noco,l_soc,jspins,nlotot,
     >     nbnd,z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >     ig,nmzxy,nmz,delz,ig2,
     >     bbmat,evac,evac_b,bkpt,bkpt_b,odi,vxy,vz,
     >     nslibd,nslibd_b,jspin,jspin_b,ico,
     >     k1,k2,k3,k1_b,k2_b,k3_b,
     >     jspd,nvd,area,nbasfcn,neigd,
     >     zMat,zMat_b,nv,nv_b,sk2,phi2,omtil,gb,gb2,qss,sign2,
     <     uHu)

      use m_constants
      use m_types
      use m_od_abvac
      use m_cylbes
      use m_dcylbs
      use m_intgr, only : intgz0
      use m_d2fdz2

      implicit none


      TYPE(t_oneD),INTENT(IN)        :: oneD
      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_stars),INTENT(IN)       :: stars
      TYPE(t_cell),INTENT(IN)        :: cell

      TYPE(t_mat), INTENT(IN) :: zMat, zMat_b

c     .. scalar Arguments..
      logical, intent (in) :: l_noco,l_soc
      integer, intent (in) :: jspins,nlotot,ico
      integer, intent (in) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n3d
      integer, intent (in) :: nmzxy,nmz,n2d,nbnd
      integer, intent (in) :: nslibd,nslibd_b
      integer, intent (in) :: jspin,jspin_b,jspd,nvd
      integer, intent (in) :: nbasfcn,neigd
      real,    intent (in) :: delz,z1,evac,area,omtil,evac_b
      complex, intent (in) :: chi
      type (od_inp), intent (in) :: odi


c     ..array arguments..
      real,    intent (in) :: bkpt(:),bkpt_b(:),qss(:) !bkpt(3),bkpt_b(3),qss(3)
      real,    intent (in) :: sk2(:),phi2(:) !sk2(n2d),phi2(n2d)
      integer, intent (in) :: ig(-k1d:,-k2d:,-k3d:) !ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig2(:),nv(:),nv_b(:) !ig2(n3d),nv(jspd),nv_b(jspd)
      integer, intent (in) :: gb(3),gb2(3)
      real,    intent (in) :: vz(nmzd,2,4),bbmat(3,3)
      integer, intent (in) :: k1(:,:),k2(:,:),k3(:,:) !k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      integer, intent (in) :: k1_b(:,:),k2_b(:,:),k3_b(:,:) !k1_b(nvd,jspd),k2_b(nvd,jspd),k3_b(nvd,jspd)
      complex, intent (inout) :: uHu(:,:) !uHu(nbnd,nbnd)

      complex, intent (in):: vxy(nmzxyd,odi%n2d-1,2)


c     ..basis wavefunctions in the vacuum
      real,    allocatable :: udz(:,:)
      real,    allocatable :: uz(:,:)
      real,    allocatable :: dudz(:,:)
      real,    allocatable :: duz(:,:)
      real,    allocatable :: u(:,:,:)
      real,    allocatable :: ud(:,:,:)
      real,    allocatable :: ddnv(:,:)

      real,    allocatable :: udz_b(:,:)
      real,    allocatable :: uz_b(:,:)
      real,    allocatable :: dudz_b(:,:)
      real,    allocatable :: duz_b(:,:)
      real,    allocatable :: u_b(:,:,:)
      real,    allocatable :: ud_b(:,:,:)
      real,    allocatable :: ddnv_b(:,:)

c     ..local scalars..
      logical tail
      real wronk,wronk1,arg,zks,tpi,zz,rk,rk_b
      real zks0,zks0_b,arg_b,fac1,xv,yv,s(3)
      integer i,m,l,j,k,irec3,irec2,n,nv2,mp,ispin,np0,ind1,ind3
      integer nv2_b,np1,lp,addnoco,addnoco2,i3,j3
      integer sign,ivac,m1,m2,jspin2,jspin2_b,jspin2H,sign2
      complex tuu,tud,tdu,tdd
      complex avac,bvac,ic,phasfc
      complex, allocatable :: acof(:,:,:),bcof(:,:,:)
      integer, allocatable :: kvac3(:),map1(:)
      complex, allocatable :: acof_b(:,:,:),bcof_b(:,:,:)
      integer, allocatable :: kvac3_b(:),map1_b(:)
      real, allocatable :: bess(:),dbss(:)
      real, allocatable :: gbess(:,:),gbess_b(:,:),besss(:)
      real, allocatable :: zmsh(:),xx(:),xximag(:),v1(:)
      real, allocatable :: gdu(:),gdud(:),gdu_b(:),gdud_b(:)
      REAL qssbti(3,2)

      addnoco=0
      addnoco2=0
      if(l_noco.and.(jspin.eq.2))then
         addnoco  = nv(1)   + nlotot
      endif
      if(l_noco.and.(jspin_b.eq.2))then
         addnoco2 = nv_b(1) + nlotot
      endif

      allocate ( udz(nv2d,-odi%mb:odi%mb),uz(nv2d,-odi%mb:odi%mb),
     +           dudz(nv2d,-odi%mb:odi%mb),
     +           duz(nv2d,-odi%mb:odi%mb),u(nmzd,nv2d,-odi%mb:odi%mb),
     +           ud(nmzd,nv2d,-odi%mb:odi%mb),ddnv(nv2d,-odi%mb:odi%mb),
     +           udz_b(nv2d,-odi%mb:odi%mb),uz_b(nv2d,-odi%mb:odi%mb),
     +           dudz_b(nv2d,-odi%mb:odi%mb),
     +           duz_b(nv2d,-odi%mb:odi%mb),
     +           u_b(nmzd,nv2d,-odi%mb:odi%mb),
     +           ud_b(nmzd,nv2d,-odi%mb:odi%mb),
     +           ddnv_b(nv2d,-odi%mb:odi%mb),
     +           bess(-odi%mb:odi%mb),dbss(-odi%mb:odi%mb),
     +           besss(-odi%M:odi%M),gbess(-odi%M:odi%M,nmzd),
     +           gbess_b(-odi%M:odi%M,nmzd),
     +           acof(nv2d,-odi%mb:odi%mb,nslibd),
     +           bcof(nv2d,-odi%mb:odi%mb,nslibd),
     +           acof_b(nv2d,-odi%mb:odi%mb,nslibd_b),
     +           bcof_b(nv2d,-odi%mb:odi%mb,nslibd_b),
     +           kvac3(nv2d),map1(nvd),
     +           kvac3_b(nv2d),map1_b(nvd),
     +           zmsh(nmz),xx(nmz),xximag(nmz),v1(nmzd),
     +           gdu(nmzd),gdud(nmzd),gdu_b(nmzd),gdud_b(nmzd) )

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)

      tail = .true.
      np0 = nmzxy + 1
      np1 = nmz + 1
      ivac = 1

      jspin2 = jspin
      jspin2_b = jspin_b
      jspin2H = ico
      if(l_soc.and.jspins.eq.1) then
         jspin2 = 1
         jspin2_b = 1
         jspin2H = 1
      endif

      sign = 1.0
      IF (ico.EQ.4) sign=-1.0  ! for c.c. of vz


      acof(:,:,:) = cmplx(0.,0.) ; bcof(:,:,:) = cmplx(0.,0.)
      acof_b(:,:,:) = cmplx(0.,0.) ; bcof_b(:,:,:) = cmplx(0.,0.)

      nv2 = 0 ; nv2_b = 0

      do 20 k = 1,nv(jspin)
         do 10 j = 1,nv2
            if (k3(k,jspin).eq.kvac3(j)) then
               map1(k) = j
               goto 20
            endif
 10      continue
         nv2 = nv2 + 1
         if (nv2.gt.nv2d) stop 'nv2d'
         kvac3(nv2) = k3(k,jspin)
         map1(k) = nv2
 20   continue

      do 21 k = 1,nv_b(jspin_b)
         do 11 j = 1,nv2_b
            if (k3_b(k,jspin_b).eq.kvac3_b(j)) then
               map1_b(k) = j
               goto 21
            endif
 11      continue
         nv2_b = nv2_b + 1
         if (nv2_b.gt.nv2d) stop 'nv2d'
         kvac3_b(nv2_b) = k3_b(k,jspin_b)
         map1_b(k) = nv2_b
 21   continue

      wronk = 2.0

c...for the k-point

      qssbti(1,1) = - qss(1)/2.
      qssbti(2,1) = - qss(2)/2.
      qssbti(1,2) = + qss(1)/2.
      qssbti(2,2) = + qss(2)/2.
      qssbti(3,1) = - qss(3)/2.
      qssbti(3,2) = + qss(3)/2.
      DO ispin = 1,1 ! jspins
      CALL od_abvac(
     >      cell,vacuum,stars,oneD,
     >      qssbti(3,jspin),odi%n2d,
     >      wronk,evac,bkpt,odi%M,odi%mb,
     >      vz(1,ivac,jspin2),kvac3,nv2,
     <      uz(:,-odi%mb:),duz(:,-odi%mb:),
     <      u(:,:,-odi%mb:),udz(1:,-odi%mb:),
     <      dudz(:,-odi%mb:),ddnv(:,-odi%mb:),
     <      ud(:,:,-odi%mb:))
      ENDDO

      do k = 1,nv(jspin)
         l = map1(k)
         irec3 = ig(k1(k,jspin),k2(k,jspin),k3(k,jspin))
         if (irec3.ne.0) then
            irec2 = ig2(irec3)
            zks = sk2(irec2)*z1
            arg = phi2(irec2)
            call cylbes(odi%mb,zks,bess)
            call dcylbs(odi%mb,zks,bess,dbss)
            do m = -odi%mb,odi%mb
               wronk1 = uz(l,m)*dudz(l,m) -
     -              udz(l,m)*duz(l,m)
               avac = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *              cmplx(dudz(l,m)*bess(m) -
     +              udz(l,m)*sk2(irec2)*dbss(m),0.0)/
     /              ((wronk1)*sqrt(omtil))
               bvac = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *              cmplx(-duz(l,m)*bess(m) +
     -              uz(l,m)*sk2(irec2)*dbss(m),0.0)/
     /              ((wronk1)*sqrt(omtil))
               IF(zMat%l_real) THEN
                  do n = 1,nslibd
                      acof(l,m,n) = acof(l,m,n) +
     +                   zMat%data_r(k+addnoco,n)*avac
c     +                    conjg(zMat%z_r(k,n))*avac
                      bcof(l,m,n) = bcof(l,m,n) +
     +                   zMat%data_r(k+addnoco,n)*bvac
c     +                    conjg(zMat%z_r(k,n))*bvac
                  enddo
               ELSE
                  do n = 1,nslibd
                      acof(l,m,n) = acof(l,m,n) +
     +                   zMat%data_c(k+addnoco,n)*avac
c     +                    conjg(zMat%z_c(k,n))*avac
                      bcof(l,m,n) = bcof(l,m,n) +
     +                   zMat%data_c(k+addnoco,n)*bvac
c     +                    conjg(zMat%z_c(k,n))*bvac
                  enddo
               END IF
            enddo      ! -mb:mb
         endif
      enddo

c...for the b-point

      DO ispin = 1,1 ! jspins
      call od_abvac(
     >      cell,vacuum,stars,oneD,
     >      qssbti(3,jspin_b),odi%n2d,
     >      wronk,evac_b,bkpt_b,odi%M,odi%mb,
     >      vz(1,ivac,jspin2_b),kvac3_b,nv2_b,
     <      uz_b(:,-odi%mb:),duz_b(:,-odi%mb:),u_b(:,:,-odi%mb:),
     <      udz_b(:,-odi%mb:),
     <      dudz_b(:,-odi%mb:),ddnv_b(:,-odi%mb:),
     <      ud_b(:,:,-odi%mb:))
      ENDDO

      do k = 1,nv_b(jspin_b)
         l = map1_b(k)
         irec3 = ig(k1_b(k,jspin_b),k2_b(k,jspin_b),k3_b(k,jspin_b))
         if (irec3.ne.0) then
            irec2 = ig2(irec3)
            zks = sk2(irec2)*z1
            arg = phi2(irec2)
            call cylbes(odi%mb,zks,bess)
            call dcylbs(odi%mb,zks,bess,dbss)
            do m = -odi%mb,odi%mb
               wronk1 = uz_b(l,m)*dudz_b(l,m) -
     -              udz_b(l,m)*duz_b(l,m)
               avac = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *              cmplx(dudz_b(l,m)*bess(m) -
     +              udz_b(l,m)*sk2(irec2)*dbss(m),0.0)/
     /              ((wronk1)*sqrt(omtil))
               bvac = exp(-cmplx(0.0,m*arg))*(ic**m)*
     &              cmplx(-duz_b(l,m)*bess(m) +
     -              uz_b(l,m)*sk2(irec2)*dbss(m),0.0)/
     /              ((wronk1)*sqrt(omtil))

               IF(zMat_b%l_real) THEN
                  do n = 1,nslibd_b
                      acof_b(l,m,n) = acof_b(l,m,n) +
     +                   zMat_b%data_r(k+addnoco2,n)*avac
c     +                    conjg(zMat_b%z_r(k,n))*avac
                      bcof_b(l,m,n) = bcof_b(l,m,n) +
     +                   zMat_b%data_r(k+addnoco2,n)*bvac
c     +                    conjg(zMat_b%z_r(k,n))*bvac
                  enddo
               ELSE
                  do n = 1,nslibd_b
                      acof_b(l,m,n) = acof_b(l,m,n) +
     +                   zMat_b%data_c(k+addnoco2,n)*avac
c     +                    conjg(zMat_b%z_c(k,n))*avac
                      bcof_b(l,m,n) = bcof_b(l,m,n) +
     +                   zMat_b%data_c(k+addnoco2,n)*bvac
c     +                    conjg(zMat_b%z_c(k,n))*bvac
                  enddo
               END IF
            enddo      ! -mb:mb
         endif
      enddo          ! k = 1,nv

c  now actually computing the uHu matrix

      irec3 = ig(gb(1),gb(2),gb(3))
      if (irec3.eq.0) stop 'Gb is not in the list of Gs'
      irec2 = ig2(irec3)
      zks0 = sk2(irec2)
      arg = phi2(irec2)

      irec3 = ig(gb2(1),gb2(2),gb2(3))
      if (irec3.eq.0) stop 'Gb2 is not in the list of Gs'
      irec2 = ig2(irec3)
      zks0_b = sk2(irec2)
      arg_b = phi2(irec2)

      gbess(:,:) = 0.
      gbess_b(:,:) = 0.
      do i = 1,nmz
         ! set up mesh
         zmsh(i) = z1 + (i-1)*delz
         zz = sqrt(zmsh(i))

         ! transformation u --> v = sqrt(z)*u
         ! thus we can use simplified 'pseudopotential'
         ! and skip r in integration dx dy = r dr dphi
         u(i,:,:) = zz*u(i,:,:)
         ud(i,:,:) = zz*ud(i,:,:)
         u_b(i,:,:) = zz*u_b(i,:,:)
         ud_b(i,:,:) = zz*ud_b(i,:,:)

         ! cyl. Bessel at G(k+b1)
         zks = zks0*zmsh(i)
         besss(:) = 0.
         call cylbes(odi%M,zks,besss)
         do m = -odi%M,odi%M
            gbess(m,i) = besss(m)
         enddo

         ! cyl. Bessel at G(k+b2)
         zks = zks0_b*zmsh(i)
         besss(:) = 0.
         call cylbes(odi%M,zks,besss)
         do m = -odi%M,odi%M
            gbess_b(m,i) = besss(m)
         enddo
      enddo


      ! calculate uHu matrix elements
      do l = 1,nv2

         j3 = kvac3(l) - gb(3)

      do lp = 1,nv2_b

         i3 = j3 - kvac3_b(lp) + gb2(3)

         ind3 = ig(0,0,sign2*i3)
         IF (ind3.EQ.0) CYCLE

         DO m = -odi%mb, odi%mb
          DO mp = -odi%mb, odi%mb
           IF ((mp.EQ.m) .OR. ((iabs(m ).LE.odi%m_cyl) .AND.
     >                         (iabs(mp).LE.odi%m_cyl))) THEN

            DO m1 = -odi%mb, odi%mb
             DO m2 = -odi%mb, odi%mb

              ind1 = odi%ig(sign2*i3,m-m1-mp+m2)  ! TODO: sign also for m-part?
              if(sign2.ne.1) stop 'sign2, check m-part, wann_uHu_od_vac'

              IF(ind1.EQ.0) CYCLE

              ind1 = ind1 - 1
              IF(ind1.NE.0) THEN ! warping components Gz=/=0, m=/=0
                                 ! need integral with vxy and cyl. Bessel
                 ! tuu
                 DO i=1,nmzxy
                    fac1 = u  (i,l ,m )*gbess  (m1,i)
     >                    *u_b(i,lp,mp)*gbess_b(m2,i)
                    xx    (np0-i) = fac1*real (vxy(i,ind1,ivac))
                    xximag(np0-i) = fac1*aimag(vxy(i,ind1,ivac))
                 ENDDO
                 call intgz0(xx,delz,nmzxy,xv,tail)
                 call intgz0(xximag,delz,nmzxy,yv,tail)
                 tuu = cmplx(xv,-yv)

                 ! tud
                 DO i=1,nmzxy
                    fac1 = u   (i,l ,m )*gbess  (m1,i)
     >                    *ud_b(i,lp,mp)*gbess_b(m2,i)
                    xx    (np0-i) = fac1*real (vxy(i,ind1,ivac))
                    xximag(np0-i) = fac1*aimag(vxy(i,ind1,ivac))
                 ENDDO
                 call intgz0(xx,delz,nmzxy,xv,tail)
                 call intgz0(xximag,delz,nmzxy,yv,tail)
                 tud = cmplx(xv,-yv)

                 ! tdu
                 DO i=1,nmzxy
                    fac1 = ud (i,l ,m )*gbess  (m1,i)
     >                    *u_b(i,lp,mp)*gbess_b(m2,i)
                    xx    (np0-i) = fac1*real (vxy(i,ind1,ivac))
                    xximag(np0-i) = fac1*aimag(vxy(i,ind1,ivac))
                 ENDDO
                 call intgz0(xx,delz,nmzxy,xv,tail)
                 call intgz0(xximag,delz,nmzxy,yv,tail)
                 tdu = cmplx(xv,-yv)

                 ! tdd
                 DO i=1,nmzxy
                    fac1 = ud  (i,l ,m )*gbess  (m1,i)
     >                    *ud_b(i,lp,mp)*gbess_b(m2,i)
                    xx    (np0-i) = fac1*real (vxy(i,ind1,ivac))
                    xximag(np0-i) = fac1*aimag(vxy(i,ind1,ivac))
                 ENDDO
                 call intgz0(xx,delz,nmzxy,xv,tail)
                 call intgz0(xximag,delz,nmzxy,yv,tail)
                 tdd = cmplx(xv,-yv)

              ELSE ! non-warping components Gz==0, m==0
                   ! need integral with XXXXXX

                 IF ((ico.EQ.1) .OR. (ico.EQ.2)) THEN ! spin-diagonal

                    ! determine second derivative of (u*gbess) etc.
                    CALL d2fdz2(nmzd,nmz,zmsh,delz,
     >                          u(:,l,m),gbess(m1,:),gdu)
                    CALL d2fdz2(nmzd,nmz,zmsh,delz,
     >                          ud(:,l,m),gbess(m1,:),gdud)
                    CALL d2fdz2(nmzd,nmz,zmsh,delz,
     >                          u_b(:,lp,mp),gbess_b(m2,:),gdu_b)
                    CALL d2fdz2(nmzd,nmz,zmsh,delz,
     >                          ud_b(:,lp,mp),gbess_b(m2,:),gdud_b)

                    ! determine |G+k+b1|^2 and |G'+k+b2|^2
                    s(1) = 0.0
                    s(2) = 0.0
                    s(3) = bkpt(3) + kvac3(l) + qssbti(3,jspin)!-gb(3)
                    rk = dot_product(s,matmul(bbmat,s))
!                    rk   = dotirp(s,s,bbmat)

                    s(1) = 0.0
                    s(2) = 0.0
                    s(3) = bkpt_b(3) + kvac3_b(lp) + qssbti(3,jspin_b)!-gb2(3)
                    rk_b = dot_product(s,matmul(bbmat,s))
!                    rk_b = dotirp(s,s,bbmat)

                    ! construct symmetrized 'pseudopotential'
                    ! TODO: valid to simply symmetrize?
                    DO i=1,nmzd
                       v1(i) = vz(i,ivac,jspin2H)
     >                       + (m*m+mp*mp)/(4.*zmsh(i)*zmsh(i))
     >                       - 1./(8.*zmsh(i)*zmsh(i))
                    ENDDO

                    ! tuu
                    DO i=1,nmz
                       fac1 = u  (i,l ,m )*gbess  (m1,i)
     >                       *u_b(i,lp,mp)*gbess_b(m2,i)
                       xx(np1-i) = fac1*( 0.25*(rk+rk_b) + v1(i) )
     >                   -0.25*( gdu  (i)*u_b(i,lp,mp)*gbess_b(m2,i)
     >                          +gdu_b(i)*u  (i,l ,m )*gbess  (m1,i) )
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    tuu = cmplx(xv,0.0)

                    ! tud
                    DO i=1,nmz
                       fac1 = u   (i,l ,m )*gbess  (m1,i)
     >                       *ud_b(i,lp,mp)*gbess_b(m2,i)
                       xx(np1-i) = fac1*( 0.25*(rk+rk_b) + v1(i) )
     >                   -0.25*( gdu   (i)*ud_b(i,lp,mp)*gbess_b(m2,i)
     >                          +gdud_b(i)*u   (i,l ,m )*gbess  (m1,i) )
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    tud = cmplx(xv,0.0)

                    ! tdu
                    DO i=1,nmz
                       fac1 = ud (i,l ,m )*gbess  (m1,i)
     >                       *u_b(i,lp,mp)*gbess_b(m2,i)
                       xx(np1-i) = fac1*( 0.25*(rk+rk_b) + v1(i) )
     >                   -0.25*( gdud (i)*u_b(i,lp,mp)*gbess_b(m2,i)
     >                          +gdu_b(i)*ud (i,l ,m )*gbess  (m1,i) )
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    tdu = cmplx(xv,0.0)

                    ! tdd
                    DO i=1,nmz
                       fac1 = ud  (i,l ,m )*gbess  (m1,i)
     >                       *ud_b(i,lp,mp)*gbess_b(m2,i)
                       xx(np1-i) = fac1*( 0.25*(rk+rk_b) + v1(i) )
     >                   -0.25*( gdud  (i)*ud_b(i,lp,mp)*gbess_b(m2,i)
     >                          +gdud_b(i)*ud  (i,l ,m )*gbess  (m1,i) )
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    tdd = cmplx(xv,0.0)

                 ELSE ! spin-off-diagonal

                    ! tuu
                    DO i=1,nmz
                       fac1 = u  (i,l ,m )*gbess  (m1,i)
     >                       *u_b(i,lp,mp)*gbess_b(m2,i)
                       xx    (np1-i) = fac1*vz(i,ivac,3)
                       xximag(np1-i) = fac1*vz(i,ivac,4)*sign
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    call intgz0(xximag,delz,nmz,yv,tail)
                    tuu = cmplx(xv,-yv)

                    ! tud
                    DO i=1,nmz
                       fac1 = u   (i,l ,m )*gbess  (m1,i)
     >                       *ud_b(i,lp,mp)*gbess_b(m2,i)
                       xx    (np1-i) = fac1*vz(i,ivac,3)
                       xximag(np1-i) = fac1*vz(i,ivac,4)*sign
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    call intgz0(xximag,delz,nmz,yv,tail)
                    tud = cmplx(xv,-yv)

                    ! tdu
                    DO i=1,nmz
                       fac1 = ud (i,l ,m )*gbess  (m1,i)
     >                       *u_b(i,lp,mp)*gbess_b(m2,i)
                       xx    (np1-i) = fac1*vz(i,ivac,3)
                       xximag(np1-i) = fac1*vz(i,ivac,4)*sign
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    call intgz0(xximag,delz,nmz,yv,tail)
                    tdu = cmplx(xv,-yv)

                    ! tdd
                    DO i=1,nmz
                       fac1 = ud  (i,l ,m )*gbess  (m1,i)
     >                       *ud_b(i,lp,mp)*gbess_b(m2,i)
                       xx    (np1-i) = fac1*vz(i,ivac,3)
                       xximag(np1-i) = fac1*vz(i,ivac,4)*sign
                    ENDDO
                    call intgz0(xx,delz,nmz,xv,tail)
                    call intgz0(xximag,delz,nmz,yv,tail)
                    tdd = cmplx(xv,-yv)

                 ENDIF ! ((ico.EQ.1) .OR. (ico.EQ.2))

              ENDIF ! (ind1.NE.0)

              ! determine phase factor
              phasfc = chi*exp(-cmplx(0.0,m2*arg_b-m1*arg))

              ! contraction of integrals with a,b coefficients
              ! yields contribution to uHu matrix
              do i = 1,nslibd
                do j = 1,nslibd_b
                   uHu(i,j) = uHu(i,j) + phasfc*area*(
     *             acof(l,m,i)*conjg(acof_b(lp,mp,j))*tuu +
     +             acof(l,m,i)*conjg(bcof_b(lp,mp,j))*tud +
     +             bcof(l,m,i)*conjg(acof_b(lp,mp,j))*tdu +
     +             bcof(l,m,i)*conjg(bcof_b(lp,mp,j))*tdd )
                enddo
              enddo

             ENDDO ! m2
            ENDDO  ! m1
           ENDIF ! noncyl. contributions
          ENDDO ! mp
         ENDDO  ! m

      enddo  ! lp
      enddo  ! l

      deallocate ( udz,uz,dudz,duz,u,ud,ddnv,bess,dbss,acof,bcof )
      deallocate ( udz_b,uz_b,dudz_b,duz_b,u_b,ud_b,ddnv_b,gbess,besss )
      deallocate ( acof_b,bcof_b,gbess_b )
      deallocate ( zmsh,xx,xximag,v1 )
      deallocate ( gdu, gdud, gdu_b, gdud_b )

      END SUBROUTINE wann_uHu_od_vac
      END MODULE m_wann_uHu_od_vac
