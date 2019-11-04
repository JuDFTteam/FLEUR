!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_mmkb_od_vac
c**************************************************************
c      Determines the overlap matrix Mmn(k,k+b) in the vacuum
c      for the wannier functions.
c      For more details see routine wannier.F and wann_mmk0_od_vac.F
c
c      Y. Mokrousov, F. Freimuth
c***************************************************************
      CONTAINS
      SUBROUTINE wann_mmkb_od_vac(
     >     oneD,vacuum,stars,cell,
     >     vacchi,l_noco,nlotot,
     >     nbnd,z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >     ig,nmzxy,nmz,delz,ig2,n2d_1,
     >     bbmat,evac,evac_b,bkpt,bkpt_b,MM,vM,vz,vz_b,odi,
     >     nslibd,nslibd_b,jspin,jspin_b,k1,k2,k3,k1_b,k2_b,k3_b,
     >     jspd,nvd,area,nbasfcn,neigd,
     >     zMat,zMat_b,nv,nv_b,sk2,phi2,omtil,gb,qss,
     >     l_q, sign_q,
     <     mmn)
      use m_constants
      use m_types
      use m_od_abvac
      use m_cylbes
      use m_dcylbs
      use m_intgr, only : intgz0

      implicit none


      TYPE(t_oneD),INTENT(IN)        :: oneD
      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_stars),INTENT(IN)       :: stars
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_mat), INTENT(IN)        :: zMat, zMat_b

c     .. scalar Arguments..
      logical, intent (in) :: l_noco
      integer, intent (in) :: nlotot
      integer, intent (in) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n3d
      integer, intent (in) :: nmzxy,nmz,MM,n2d,vM,nbnd
      integer, intent (in) :: nslibd,nslibd_b
      integer, intent (in) :: n2d_1,jspin,jspin_b,jspd,nvd
      integer, intent (in) :: nbasfcn,neigd
      real,    intent (in) :: delz,z1,evac,area,omtil,evac_b
      complex, intent (in) :: vacchi
      type (od_inp), intent (in) :: odi

      logical, intent (in) :: l_q
      integer, intent (in) :: sign_q

c     ..array arguments..
      real,    intent (in) :: bkpt(:),bkpt_b(:),qss(:) !bkpt(3),bkpt_b(3),qss(3)
      real,    intent (in) :: sk2(:),phi2(:) !sk2(n2d),phi2(n2d)
      integer, intent (in) :: ig(-k1d:,-k2d:,-k3d:) !ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig2(:),nv(:),nv_b(:) !ig2(n3d),nv(jspd),nv_b(jspd)
      integer, intent (in) :: gb(:) !gb(3)
      real,    intent (in) :: vz(:),vz_b(:),bbmat(:,:) !vz(nmzd),bbmat(3,3)
      integer, intent (in) :: k1(:,:),k2(:,:),k3(:,:) !k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      integer, intent (in) :: k1_b(:,:),k2_b(:,:),k3_b(:,:) !k1_b(nvd,jspd),k2_b(nvd,jspd),k3_b(nvd,jspd)
      complex, intent (inout) :: mmn(:,:) !mmn(nbnd,nbnd)

      !logical, intent (in) :: l_q   ! dealing with q points?

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
      real wronk,wronk1,arg,zks,tpi,uuo,udo,duo,ddo,zz,xx(nmz),zks0
      integer i,m,l,j,k,irec3,irec2,n,nv2,mp,ispin
      integer nv2_b,np1,lprime,addnoco,addnoco2
      complex avac,bvac,ic,phasfc
      complex, allocatable :: acof(:,:,:),bcof(:,:,:)
      integer, allocatable :: kvac3(:),map1(:)
      complex, allocatable :: acof_b(:,:,:),bcof_b(:,:,:)
      integer, allocatable :: kvac3_b(:),map1_b(:)
      real, allocatable :: bess(:),dbss(:)
      real, allocatable :: gbess(:,:),besss(:)
      REAL qssbti(3,2)

      addnoco=0
      addnoco2=0
      if(l_noco.and.(jspin.eq.2))then
         addnoco  = nv(1)   + nlotot
      endif
      if(l_noco.and.(jspin_b.eq.2))then
         addnoco2 = nv_b(1) + nlotot
      endif

      allocate ( udz(nv2d,-vM:vM),uz(nv2d,-vM:vM),dudz(nv2d,-vM:vM),
     +           duz(nv2d,-vM:vM),u(nmzd,nv2d,-vM:vM),
     +           ud(nmzd,nv2d,-vM:vM),ddnv(nv2d,-vM:vM),
     +           udz_b(nv2d,-vM:vM),uz_b(nv2d,-vM:vM),
     +           dudz_b(nv2d,-vM:vM),
     +           duz_b(nv2d,-vM:vM),u_b(nmzd,nv2d,-vM:vM),
     +           ud_b(nmzd,nv2d,-vM:vM),ddnv_b(nv2d,-vM:vM),
     +           bess(-odi%mb:odi%mb),dbss(-odi%mb:odi%mb),
     +           besss(-odi%M:odi%M),gbess(-odi%M:odi%M,nmzd),
     +           acof(nv2d,-odi%mb:odi%mb,nslibd),
     +           bcof(nv2d,-odi%mb:odi%mb,nslibd),
     +           acof_b(nv2d,-odi%mb:odi%mb,nslibd_b),
     +           bcof_b(nv2d,-odi%mb:odi%mb,nslibd_b),
     +           kvac3(nv2d),map1(nvd),
     +           kvac3_b(nv2d),map1_b(nvd) )

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)

      tail = .true.
      np1 = nmz + 1

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
     >      vz,kvac3,nv2,
     <      uz(:,-vM:),duz(:,-vM:),u(:,:,-vM:),udz(:,-vM:),
     <      dudz(:,-vM:),ddnv(:,-vM:),ud(:,:,-vM:))
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
               IF (zMat%l_real) THEN
                  do n = 1,nslibd
                      acof(l,m,n) = acof(l,m,n) +
     +                   zMat%data_r(k+addnoco,n)*avac
c     +                    conjg(zMat%data_r(k,n))*avac
                      bcof(l,m,n) = bcof(l,m,n) +
     +                   zMat%data_r(k+addnoco,n)*bvac
c     +                    conjg(zMat%data_r(k,n))*bvac
                  enddo
               ELSE
                  do n = 1,nslibd
                      acof(l,m,n) = acof(l,m,n) +
     +                   zMat%data_c(k+addnoco,n)*avac
c     +                    conjg(zMat%data_c(k,n))*avac
                      bcof(l,m,n) = bcof(l,m,n) +
     +                   zMat%data_c(k+addnoco,n)*bvac
c     +                    conjg(zMat%data_c(k,n))*bvac
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
     >      vz_b,kvac3_b,nv2_b,
     <      uz_b(1:,-vM:),duz_b(1:,-vM:),u_b(1:,1:,-vM:),udz_b(1:,-vM:),
     <      dudz_b(1:,-vM:),ddnv_b(1:,-vM:),ud_b(1:,1:,-vM:))
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
               IF (zMat_b%l_real) THEN
                  do n = 1,nslibd_b
                      acof_b(l,m,n) = acof_b(l,m,n) +
     +                   zMat_b%data_r(k+addnoco2,n)*avac
c     +                    conjg(zMat_b%data_r(k,n))*avac
                      bcof_b(l,m,n) = bcof_b(l,m,n) +
     +                   zMat_b%data_r(k+addnoco2,n)*bvac
c     +                    conjg(zMat_b%data_r(k,n))*bvac
                  enddo
               ELSE
                  do n = 1,nslibd_b
                      acof_b(l,m,n) = acof_b(l,m,n) +
     +                   zMat_b%data_c(k+addnoco2,n)*avac
c     +                    conjg(zMat_b%data_c(k,n))*avac
                      bcof_b(l,m,n) = bcof_b(l,m,n) +
     +                   zMat_b%data_c(k+addnoco2,n)*bvac
c     +                    conjg(zMat_b%data_c(k,n))*bvac
                  enddo
               END IF
            enddo      ! -mb:mb
         endif
      enddo          ! k = 1,nv

c  now actually computing the Mmn matrix

      irec3 = ig(gb(1),gb(2),gb(3))
      if (irec3.eq.0) stop 'Gb is not in the list of Gs'
      irec2 = ig2(irec3)
      zks0 = sk2(irec2)
      arg = phi2(irec2)

      !write(*,*)'mmkb_od_vac zks0',zks0

      gbess(:,:) = 0.
      do i = 1,nmz
         zz = z1+(i-1)*delz
         zks = zks0*zz
         besss(:) = 0.
         call cylbes(odi%M,zks,besss)
         do m = -odi%M,odi%M
            gbess(m,i) = besss(m)
c           gbess(i) = 1.
         enddo
      enddo

      do l = 1,nv2
      do lprime = 1,nv2_b

      if (kvac3(l).eq.(kvac3_b(lprime)-gb(3))) then

       do m = -odi%mb,odi%mb
        do mp = -odi%mb,odi%mb

         phasfc = exp(-cmplx(0.0,(mp-m)*arg))
         if(l_q) phasfc =  (sign_q)**(mp-m) * phasfc

         do i = 1,nmz
           zz = z1+(i-1)*delz
           xx(np1-i) = zz*u(i,l,m)*u_b(i,lprime,mp)*gbess(mp-m,i)
         enddo
         call intgz0(xx,delz,nmz,uuo,tail)

         do i = 1,nmz
           zz = z1+(i-1)*delz
           xx(np1-i) = zz*u(i,l,m)*ud_b(i,lprime,mp)*gbess(mp-m,i)
         enddo
         call intgz0(xx,delz,nmz,udo,tail)

         do i = 1,nmz
           zz = z1+(i-1)*delz
           xx(np1-i) = zz*ud(i,l,m)*u_b(i,lprime,mp)*gbess(mp-m,i)
         enddo
         call intgz0(xx,delz,nmz,duo,tail)

         do i = 1,nmz
           zz = z1+(i-1)*delz
           xx(np1-i) = zz*ud(i,l,m)*ud_b(i,lprime,mp)*gbess(mp-m,i)
         enddo
         call intgz0(xx,delz,nmz,ddo,tail)

         do i = 1,nslibd
           do j = 1,nslibd_b
            mmn(i,j) = mmn(i,j) + phasfc*area*(
     *          acof(l,m,i)*conjg(acof_b(lprime,mp,j))*uuo +
     +          acof(l,m,i)*conjg(bcof_b(lprime,mp,j))*udo +
     +          bcof(l,m,i)*conjg(acof_b(lprime,mp,j))*duo +
     +          bcof(l,m,i)*conjg(bcof_b(lprime,mp,j))*ddo )*vacchi
           enddo
         enddo

        enddo !mp
       enddo  !m

      endif  ! kvac(k)=kvac(k+b)

      enddo  ! lprime
      enddo  ! l

      deallocate ( udz,uz,dudz,duz,u,ud,ddnv,bess,dbss,acof,bcof )
      deallocate ( udz_b,uz_b,dudz_b,duz_b,u_b,ud_b,ddnv_b,gbess,besss )
      deallocate ( acof_b,bcof_b )

      END SUBROUTINE wann_mmkb_od_vac
      END MODULE m_wann_mmkb_od_vac
