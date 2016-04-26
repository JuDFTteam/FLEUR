      MODULE m_wann_mmk0_vac
      use m_juDFT
c**************************************************************
c      Determines the overlap matrix Mmn(k) in the vacuum
c      in the film case for the wannier functions.
c      For more details see routine wannier.F 
c
c      Y.Mokrousov, F. Freimuth
c*************************************************************** 
      CONTAINS
      SUBROUTINE wann_mmk0_vac(
     >     l_noco,nlotot,qss,
     >     z1,nmzd,nv2d,k1d,k2d,k3d,n3d,nvac,
     >     ig,nmz,delz,ig2,area,bmat,
     >     bbmat,evac,bkpt,vz,
     >     nslibd,jspin,k1,k2,k3,jspd,nvd,
     >     nbasfcn,neigd,z,nv,omtil,
     <     mmn)

      use m_constants, only : pimach
      use m_dotir
      implicit none

c     .. scalar Arguments..
      logical, intent (in) :: l_noco
      integer, intent (in) :: nlotot
      real,    intent (in) :: qss(:) !qss(3)
      integer, intent (in) :: nmzd,nv2d,k1d,k2d,k3d,n3d
      integer, intent (in) :: nmz,nslibd,nvac
      integer, intent (in) :: jspin,jspd,nvd
      integer, intent (in) :: nbasfcn,neigd
      real,    intent (in) :: delz,z1,omtil,area

c     ..array arguments..
      real,    intent (in) :: bkpt(:) !bkpt(3)
      real,    intent (in) :: evac(:) !evac(2)
      integer, intent (in) :: ig(-k1d:,-k2d:,-k3d:) !ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig2(:) !ig2(n3d)
      integer, intent (in) :: nv(:) !nv(jspd)
      real,    intent (in) :: vz(:,:) !vz(nmzd,2)
      real,    intent (in) :: bbmat(3,3),bmat(3,3)
      integer, intent (in) :: k1(:,:) !k1(nvd,jspd)
      integer, intent (in) :: k2(:,:) !k2(nvd,jspd)
      integer, intent (in) :: k3(:,:) !k3(nvd,jspd)
      complex, intent (inout) :: mmn(:,:) !mmn(nslibd,nslibd)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex, intent (in):: z(:,:) !z(nbasfcn,neigd)
#else
      real,    intent (in):: z(:,:) !z(nbasfcn,neigd)
#endif

c     ..basis wavefunctions in the vacuum
      complex, allocatable :: ac(:,:),bc(:,:)
      real,    allocatable :: dt(:),dte(:)
      real,    allocatable :: t(:),te(:),tei(:)
      real,    allocatable :: u(:,:),ue(:,:),v(:)

c     ..local scalars..
      real wronk,arg,zks,tpi,vz0(2),scale,evacp,ev,const
      real :: qss1,qss2
      integer i,m,l,j,k,n,nv2,ivac,n2,sign,ik,symvac,addnoco
      integer symvacvac
      complex av,bv,ic,c_1
      integer, allocatable :: kvac1(:),kvac2(:),map2(:)
c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt

c..   external functions
      external vacudz,vacuz
c      real dotirp
c      external dotirp !module now

      allocate ( ac(nv2d,neigd),bc(nv2d,neigd),dt(nv2d),
     +           dte(nv2d),t(nv2d),te(nv2d),tei(nv2d),
     +           u(nmzd,nv2d),ue(nmzd,nv2d),
     +           v(3),kvac1(nv2d),kvac2(nv2d),map2(nvd) )

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)

c.. determining the indexing array (in-plane stars)

      wronk = 2.0
      const = 1.0 / ( sqrt(omtil)*wronk )

      do ivac = 1,2
         vz0(ivac) = vz(nmz,ivac)
      enddo

      addnoco=0
      if(l_noco .and. jspin.eq.2)then
        addnoco=nv(1)+nlotot
      endif

      n2 = 0
      do 40 k = 1,nv(jspin)
         do 30 j = 1,n2
            if ( k1(k,jspin).eq.kvac1(j) .and.
     +          k2(k,jspin).eq.kvac2(j) ) then
                map2(k) = j
                goto 40
             endif 
 30      continue
         n2 = n2 + 1
    
         IF (n2>nv2d)  CALL juDFT_error("wannier Mmn vac",calledby
     +        ="wann_mmk0_vac")

         kvac1(n2) = k1(k,jspin)
         kvac2(n2) = k2(k,jspin)
         map2(k) = n2
 40   continue

c...cycle by the vacua
      do 140 ivac = 1,nvac


       sign = 3. - 2.*ivac
       evacp = evac(ivac)

       nv2 = n2

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
         ev = evacp - 0.5*dotirp(v,v,bbmat)
         call vacuz(ev,vz(1,ivac),vz0(ivac),nmz,delz,t(ik),dt(ik),
     +        u(1,ik))
         call vacudz(ev,vz(1,ivac),vz0(ivac),nmz,delz,te(ik),
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
c-----> construct a and b coefficients

       symvacvac=1
       if (nvac==1) symvacvac=2
       do symvac=1,symvacvac
        do 60 n = 1,nslibd
            do 50 i = 1,nv2d
               ac(i,n) = cmplx(0.0,0.0)
               bc(i,n) = cmplx(0.0,0.0)
 50         continue
 60     continue

        if (symvac==2) sign=-1.0   
          
        do k = 1,nv(jspin)
          l = map2(k)
          zks = k3(k,jspin)*bmat(3,3)*sign
          arg = zks*z1
          c_1 = cmplx(cos(arg),sin(arg)) * const
          av = -c_1 * cmplx( dte(l),zks*te(l) )
          bv =  c_1 * cmplx(  dt(l),zks* t(l) )
c-----> loop over basis functions
          do n = 1,nslibd
            ac(l,n) = ac(l,n) + z(k+addnoco,n)*av
            bc(l,n) = bc(l,n) + z(k+addnoco,n)*bv
          enddo
        enddo

        do l = 1,nv2 
         do i = 1,nslibd
            do j = 1,nslibd 
               mmn(i,j) = mmn(i,j) +
     +               area*(ac(l,i)*conjg( ac(l,j))
     +              + tei(l)*bc(l,i)*conjg( bc(l,j)))
            enddo 
         enddo
        enddo
       enddo !symvac 

c... cycle by the vacua finishes
 140  enddo      

      deallocate ( ac,bc,dt,dte,t,te,tei,u,ue,
     +             v,kvac1,kvac2,map2 )

      end subroutine wann_mmk0_vac
      end module m_wann_mmk0_vac
