c********************************************************
c     calculate a-, and b-coefficients of 1d-vacuum
c     Y. Mokrousov from Frank's routine, January 2007
c********************************************************
      module m_wann_1dvacabcof
      contains
      subroutine wann_1dvacabcof(
     >    DIMENSION,oneD,vacuum,stars,cell,
     >    nv2d,nslibd,nmzd,nmz,omtil,vz,
     >    nv,bkpt,z1,odi,ods,
     >    nvd,k1,k2,k3,evac,
     >    bbmat,delz,bmat,nbasfcn,neigd,zMat,
     >    n2d,n3d,ig,nmzxy,nmzxyd,ig2,sk2,
     >    phi2,k1d,k2d,k3d,
     <    ac,bc,u,ue,addnoco,l_ss,qss,jspin)

      use m_types
      use m_od_abvac
      use m_constants
      use m_cylbes
      use m_dcylbs
      implicit none

      TYPE(t_dimension),INTENT(IN)   :: DIMENSION
      TYPE(t_oneD),INTENT(IN)        :: oneD
      TYPE(t_vacuum),INTENT(IN)      :: vacuum
      TYPE(t_stars),INTENT(IN)       :: stars
      TYPE(t_cell),INTENT(IN)        :: cell
      TYPE(t_zMat),INTENT(IN)        :: zMat

      integer,intent(in)::nv2d,n2d,n3d
      integer,intent(in)::nslibd
      integer,intent(in)::k1d,k2d,k3d
      integer,intent(in)::nmzd,nmzxyd
      integer,intent(in)::nmz,nmzxy
      integer,intent(in)::nbasfcn,neigd
      real,intent(in)::omtil
      real,intent(in)::vz(nmzd,2)
      real,intent(in)::evac(2)
      real,intent(in)::bbmat(3,3)
      real,intent(in)::delz
      real,intent(in)::bmat(3,3)
      real,intent(in)::z1
      integer,intent(in)::nv
      integer,intent(in)::nvd
      integer,intent(in)::k1(nvd)
      integer,intent(in)::k2(nvd)
      integer,intent(in)::k3(nvd)
      integer,intent(in)::ig(-k1d:k1d,-k2d:k2d,-k3d:k3d),ig2(n3d)
      integer,intent(in)::addnoco,jspin
      real,intent(in) :: sk2(n2d),phi2(n2d)
      real,intent(in)::bkpt(3),qss(3)
      logical,intent(in)::l_ss

      type (od_inp), intent (in) :: odi
      type (od_sym), intent (in) :: ods

      complex,intent(out)::ac(nv2d,-odi%mb:odi%mb,nslibd)
      complex,intent(out)::bc(nv2d,-odi%mb:odi%mb,nslibd)
      real,intent(out)::u(nmzd,nv2d,-odi%mb:odi%mb)
      real,intent(out)::ue(nmzd,nv2d,-odi%mb:odi%mb)

      real wronk,const,wronk_1
      complex av,bv,ic
      real,    allocatable :: dt(:,:),dte(:,:)
      real,    allocatable :: t(:,:),te(:,:),tei(:,:)
      integer n2,k,nv2,ik,n,l,ispin
      real vz0,evacp,ev,zks,arg,tpi
      integer kvac3(nv2d),map1(nvd),i,j,irec3,irec2,m,nvac
      real bess(-odi%mb:odi%mb),dbss(-odi%mb:odi%mb)
      real qssbti(3,2)
         ic = cmplx(0.,1.)
         tpi = 2.*pimach()
         wronk = 2.0
         const = 1.0 / ( sqrt(omtil)*wronk )
         allocate (dt(nv2d,-odi%mb:odi%mb),dte(nv2d,-odi%mb:odi%mb),
     &             t(nv2d,-odi%mb:odi%mb),te(nv2d,-odi%mb:odi%mb),
     &             tei(nv2d,-odi%mb:odi%mb))

         nvac = 1
         vz0 = vz(nmz,nvac)

         n2 = 0
         do 35 k = 1,nv
            do 45 j = 1,n2
               if (k3(k).eq.kvac3(j)) then
                  map1(k) = j
                  goto 35
               end if
 45         continue
            n2 = n2 + 1
            if (n2.gt.nv2d) stop 'wann_plot:vac'
            kvac3(n2) =  k3(k)
            map1(k) = n2
 35      continue
         nv2 = n2

         ac(:,:,:) = cmplx(0.,0.) ; bc (:,:,:) = cmplx(0.,0.)

         evacp = evac(1)
         if(l_ss) then
            qssbti(1,1)=-qss(1)/2.
            qssbti(2,1)=-qss(2)/2.
            qssbti(3,1)=-qss(3)/2.
            qssbti(1,2)=+qss(1)/2.
            qssbti(2,2)=+qss(2)/2.
            qssbti(3,2)=+qss(3)/2.
         else
            qssbti(:,:)=0.0
         endif
         ispin=1
         call od_abvac(
     >           cell,vacuum,DIMENSION,stars,oneD,
     >           qssbti(3,jspin),odi%n2d,wronk,evacp,bkpt,
     >           odi%M,odi%mb,vz(1,nvac),kvac3(1),nv2,
     >           t(1,-odi%mb),dt(1,-odi%mb),u(1,1,-odi%mb),
     <           te(1,-odi%mb),dte(1,-odi%mb),tei(1,-odi%mb),
     <           ue(1,1,-odi%mb))

         do k = 1,nv
            l = map1(k)
            irec3 = ig(k1(k),k2(k),k3(k))
            if (irec3.ne.0) then
               irec2 = ig2(irec3)
               zks = sk2(irec2)*z1
               arg = phi2(irec2)
               call cylbes(odi%mb,zks,bess)
               call dcylbs(odi%mb,zks,bess,dbss)
               do m = -odi%mb,odi%mb
                  wronk_1 = t(l,m)*dte(l,m) - te(l,m)*dt(l,m)
                  av = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *                 cmplx(dte(l,m)*bess(m) -
     +                 te(l,m)*sk2(irec2)*dbss(m),0.0)/
     /                 ((wronk_1)*sqrt(omtil))
                  bv = exp(-cmplx(0.0,m*arg))*(ic**m)*
     *                 cmplx(-dt(l,m)*bess(m) +
     -                 t(l,m)*sk2(irec2)*dbss(m),0.0)/
     /                 ((wronk_1)*sqrt(omtil))
                  IF (zMat%l_real) THEN
                     do n = 1,nslibd
                        ac(l,m,n) = ac(l,m,n) + zMat%z_r(k+addnoco,n)*av
                        bc(l,m,n) = bc(l,m,n) + zMat%z_r(k+addnoco,n)*bv
                     end do
                  ELSE
                     do n = 1,nslibd
                        ac(l,m,n) = ac(l,m,n) + zMat%z_c(k+addnoco,n)*av
                        bc(l,m,n) = bc(l,m,n) + zMat%z_c(k+addnoco,n)*bv
                     end do
                  END IF
               end do      ! -mb:mb
            end if
         end do         ! k = 1,nv

         deallocate (dt,dte,t,te,tei)

      end subroutine
      end module m_wann_1dvacabcof



