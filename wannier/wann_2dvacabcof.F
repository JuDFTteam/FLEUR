!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_2dvacabcof
      use m_juDFT
c********************************************************
c     calculate a-, and b-coefficients of 2d-vacuum
c     Frank Freimuth, November 2006
c********************************************************
      CONTAINS
      SUBROUTINE wann_2dvacabcof(
     >nv2d,nslibd,nvac,nmzd,nmz,omtil,vz,nv,bkpt,z1,
     >nvd,k1,k2,k3,evac,bbmat,delz,bmat,nbasfcn,neigd,z,
     <ac,bc,u,ue)
      
      implicit none
      integer,intent(in)::nv2d
      integer,intent(in)::nslibd
      integer,intent(in)::nvac
      integer,intent(in)::nmzd
      integer,intent(in)::nmz
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
      real,intent(in)::bkpt(3)
      complex,intent(out)::ac(nv2d,nslibd,2)
      complex,intent(out)::bc(nv2d,nslibd,2)
      real,intent(out)::u(nmzd,nv2d,nvac)
      real,intent(out)::ue(nmzd,nv2d,nvac)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      COMPLEX, INTENT (IN) :: z(nbasfcn,neigd)
#else
      REAL, INTENT (IN) :: z(nbasfcn,neigd)
#endif
      real wronk,const
      complex c_1,av,bv
      real,    allocatable :: dt(:),dte(:)
      real,    allocatable :: t(:),te(:),tei(:)
      integer ivac,n2,k,nv2,ik,jvac,symvac,symvacvac,n,l
      real vz0(2),evacp,sign,v(3),ev,scale,zks,arg
      integer kvac1(nv2d),kvac2(nv2d),map2(nvd),i,j

         wronk = 2.0
         const = 1.0 / ( sqrt(omtil)*wronk )
         allocate (dt(nv2d),dte(nv2d),t(nv2d),te(nv2d),tei(nv2d))

         do ivac = 1,2
            vz0(ivac) = vz(nmz,ivac)
         enddo


         n2 = 0 

         do 40 k = 1,nv
            do 30 j = 1,n2
               if ( k1(k).eq.kvac1(j) .and.
     +          k2(k).eq.kvac2(j) ) then
                map2(k) = j
                goto 40
               endif 
 30         continue
            n2 = n2 + 1
            IF (n2>nv2d)  CALL juDFT_error("wann_plot: vac",calledby
     +           ="wann_2dvacabcof")
            kvac1(n2) = k1(k)
            kvac2(n2) = k2(k)
            map2(k) = n2
 40      continue
         nv2=n2
         do ivac=1,nvac  
            evacp=evac(ivac)
            sign=3-2*ivac
            do ik = 1,nv2
               v(1) = bkpt(1) + kvac1(ik)
               v(2) = bkpt(2) + kvac2(ik)
               v(3) = 0.
               ev = evacp - 0.5*dot_product(matmul(v,bbmat),v)
               call vacuz(ev,vz(1,ivac),vz0(ivac),nmz,delz,t(ik),
     +        dt(ik),
     +        u(1,ik,ivac))
               call vacudz(ev,vz(1,ivac),vz0(ivac),nmz,delz,te(ik),
     +        dte(ik),tei(ik),ue(1,ik,ivac),dt(ik),
     +        u(1,ik,ivac))
               scale = wronk/ (te(ik)*dt(ik)-
     -                dte(ik)*t(ik))
               te(ik) = scale*te(ik)
               dte(ik) = scale*dte(ik)
               tei(ik) = scale*tei(ik)
               do j = 1,nmz
                  ue(j,ik,ivac) = scale*ue(j,ik,ivac)
               enddo
            enddo

c            do l=1,nv2
c               do j=1,nmz
c                  if (abs(ue(j,l,ivac)).gt.10)then
c                     print*,"l=",l
c                     print*,"j=",j
c                     print*,"ue(j,l,ivac)=",ue(j,l,ivac)
c                  endif   
c               enddo   
c            enddo   
      
            jvac=ivac
            symvacvac=1
            if (nvac==1) symvacvac=2
            do symvac=1,symvacvac
               if(symvac==2) then
                  sign=-1.0
                  jvac=2
               endif   

               do i = 1,nv2d
                  do n = 1,nslibd
                     ac(i,n,jvac) = cmplx(0.0,0.0)
                     bc(i,n,jvac) = cmplx(0.0,0.0)
                  enddo   
               enddo   

               do k = 1,nv
                  l = map2(k)
                  zks = k3(k)*bmat(3,3)*sign
                  arg = zks*z1
                  c_1 = cmplx(cos(arg),sin(arg)) * const
                  av = -c_1 * cmplx( dte(l),zks*te(l) )
                  bv =  c_1 * cmplx(  dt(l),zks* t(l) )
c-----> loop over basis functions
                  do n = 1,nslibd
                     ac(l,n,jvac) = ac(l,n,jvac) + z(k,n)*av
                     bc(l,n,jvac) = bc(l,n,jvac) + z(k,n)*bv
                  enddo
           
               enddo
            enddo !symvac    
         enddo               !loop over ivac
      deallocate (dt,dte,t,te,tei)
      END SUBROUTINE
      END MODULE m_wann_2dvacabcof
