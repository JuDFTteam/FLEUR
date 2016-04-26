      MODULE m_wann_mmk0_od_vac
      use m_juDFT
c**************************************************************
c      Determines the overlap matrix Mmn(k) in the vacuum
c      for the wannier functions.
c      For more details see routine wannier.F 
c
c      Y. Mokrousov, F. Freimuth 
c*************************************************************** 

      CONTAINS

      SUBROUTINE wann_mmk0_od_vac(
     >     l_noco,nlotot,
     >     z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >     ig,nmzxy,nmz,delz,ig2,n2d_1,
     >     bbmat,evac,bkpt,MM,vM,vz,odi,
     >     nslibd,jspin,k1,k2,k3,jspd,nvd,area,
     >     nbasfcn,neigd,z,nv,sk2,phi2,omtil,qss,
     <     mmn)

      use m_constants, only : pimach
      use m_od_types, only : od_inp
      use m_od_abvac
      use m_cylbes
      use m_dcylbs

      implicit none
c     .. scalar Arguments..
      logical, intent (in) :: l_noco
      integer, intent (in) :: nlotot
      integer, intent (in) :: nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n3d
      integer, intent (in) :: nmzxy,nmz,MM,n2d,vM,nslibd
      integer, intent (in) :: n2d_1,jspin,jspd,nvd
      integer, intent (in) :: nbasfcn,neigd
      real,    intent (in) :: delz,z1,evac,area,omtil
      type (od_inp), intent (in) :: odi

c     ..array arguments..
      real,    intent (in) :: bkpt(3),qss(3)
      real,    intent (in) :: sk2(:),phi2(:) !sk2(n2d),phi2(n2d)
      integer, intent (in) :: ig(-k1d:,-k2d:,-k3d:) !ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig2(:),nv(:) !ig2(n3d),nv(jspd)
      real,    intent (in) :: vz(:),bbmat(3,3)!vz(nmzd),bbmat(3,3)
      integer, intent (in) :: k1(:,:),k2(:,:),k3(:,:)!k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      complex, intent (inout) :: mmn(:,:)!mmn(nslibd,nslibd)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex, intent (in):: z(:,:) !z(nbasfcn,neigd)
#else
      real,    intent (in):: z(:,:) !z(nbasfcn,neigd)
#endif

c     ..basis wavefunctions in the vacuum
      real,    allocatable :: udz(:,:)
      real,    allocatable :: uz(:,:)
      real,    allocatable :: dudz(:,:)
      real,    allocatable :: duz(:,:)
      real,    allocatable :: u(:,:,:)
      real,    allocatable :: ud(:,:,:)
      real,    allocatable :: ddnv(:,:)

c     ..local scalars..
      real :: wronk,wronk1,arg,zks,tpi
      integer :: i,m,l,j,k,irec3,irec2,n,nv2,ispin,addnoco
      complex :: avac,bvac,ic
      complex, allocatable :: acof(:,:,:),bcof(:,:,:)
      integer, allocatable :: kvac3(:),map1(:)
      real, allocatable :: bess(:),dbss(:)
      real :: qssbti(3,2)

c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt

      allocate ( udz(nv2d,-vM:vM),uz(nv2d,-vM:vM),
     +           dudz(nv2d,-vM:vM),
     +           duz(nv2d,-vM:vM),u(nmzd,nv2d,-vM:vM),
     +           ud(nmzd,nv2d,-vM:vM),
     +           ddnv(nv2d,-vM:vM),
     +           bess(-odi%mb:odi%mb),dbss(-odi%mb:odi%mb),
     +           acof(nv2d,-odi%mb:odi%mb,nslibd),
     +           bcof(nv2d,-odi%mb:odi%mb,nslibd),
     +           kvac3(nv2d),map1(nvd) )

      acof=cmplx(0.0,0.0)
      bcof=cmplx(0.0,0.0)

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)

      nv2 = 0

      do 20 k = 1,nv(jspin)
         do 10 j = 1,nv2
            if (k3(k,jspin).eq.kvac3(j)) then
               map1(k) = j
               goto 20
            endif
 10      continuE
         nv2 = nv2 + 1
         IF (nv2>nv2d)  CALL juDFT_error("nv2d",calledby
     +        ="wann_mmk0_od_vac")
         kvac3(nv2) = k3(k,jspin)
         map1(k) = nv2
 20   continue

      wronk = 2.0

      qssbti(1,1) = - qss(1)/2.     ! noco (ss) case not
      qssbti(2,1) = - qss(2)/2.     ! implemented, just
      qssbti(1,2) = + qss(1)/2.     ! for compatibility
      qssbti(2,2) = + qss(2)/2.
      qssbti(3,1) = - qss(3)/2.
      qssbti(3,2) = + qss(3)/2.
      DO ispin = 1,1 ! jspins
      CALL od_abvac(
     >      z1,nmzxyd,nmzd,nv2d,k1d,k2d,k3d,n2d,n3d,
     >      ig,odi%ig,tpi,qssbti(3,jspin),
     >      nmzxy,nmz,delz,ig2,odi%n2d,
     >      bbmat,wronk,evac,bkpt,odi%M,odi%mb,
     >      vz,kvac3(1),nv2,
     <      uz(1,-vM),duz(1,-vM),u(1,1,-vM),udz(1,-vM),
     <      dudz(1,-vM),ddnv(1,-vM),ud(1,1,-vM))
      ENDDO

      addnoco=0
      if(l_noco.and.(jspin.eq.2))then
         addnoco= nv(1)+nlotot
      endif
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
               do n = 1,nslibd
                   acof(l,m,n) = acof(l,m,n) +
     +                z(k+addnoco,n)*avac
c     +                 conjg(z(k,n))*avac
                   bcof(l,m,n) = bcof(l,m,n) +
     +                z(k+addnoco,n)*bvac
c     +                 conjg(z(k,n))*bvac
               enddo
            enddo      ! -mb:mb
         endif 
      enddo          ! k = 1,nv

c  now actually computing the Mmn matrix

      do l = 1,nv2
       do m = -odi%mb,odi%mb
         do i = 1,nslibd
           do j = 1,nslibd 
            mmn(i,j) = mmn(i,j) + 
     +                area*(acof(l,m,i)*conjg(acof(l,m,j))
     +              + ddnv(l,m)*bcof(l,m,i)*conjg(bcof(l,m,j)))  
           enddo 
         enddo
       enddo
      enddo

      deallocate ( udz,uz,dudz,duz,u,ud,ddnv,bess,dbss,acof,bcof )

      END SUBROUTINE wann_mmk0_od_vac
      END MODULE m_wann_mmk0_od_vac
