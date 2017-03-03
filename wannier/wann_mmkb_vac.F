!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_mmkb_vac
      use m_juDFT
c**************************************************************
c      Determines the overlap matrix Mmn(k) in the vacuum
c      in the film case for the wannier functions.
c      For more details see routine wannier.F 
c
c      Y. Mokrousov, F. Freimuth
c*************************************************************** 
      CONTAINS
      SUBROUTINE wann_mmkb_vac(
     >     l_noco,nlotot,qss,
     >     nbnd,z1,nmzd,nv2d,k1d,k2d,k3d,n3d,nvac,
     >     ig,nmz,delz,ig2,area,bmat,
     >     bbmat,evac,bkpt,bkpt_b,vz,
     >     nslibd,nslibd_b,jspin,
     >     k1,k2,k3,k1_b,k2_b,k3_b,jspd,nvd,
     >     nbasfcn,neigd,z,z_b,nv,nv_b,omtil,gb,
     <     mmn)

      use m_constants, only : pimach
      use m_intgr, only : intgz0
      implicit none
c     .. scalar Arguments..
      logical, intent (in) :: l_noco
      integer, intent (in) :: nlotot
      real,    intent (in) :: qss(3)
      integer, intent (in) :: nmzd,nv2d,k1d,k2d,k3d,n3d,nbnd
      integer, intent (in) :: nmz,nslibd,nslibd_b,nvac
      integer, intent (in) :: jspin,jspd,nvd
      integer, intent (in) :: nbasfcn,neigd
      integer, intent (in) :: gb(3)
      real,    intent (in) :: delz,z1,omtil,area

c     ..array arguments..
      real,    intent (in) :: bkpt(3),bkpt_b(3),evac(2)
      integer, intent (in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent (in) :: ig2(n3d),nv(jspd),nv_b(jspd)
      real,    intent (in) :: vz(nmzd,2),bbmat(3,3),bmat(3,3)
      integer, intent (in) :: k1(nvd,jspd),k2(nvd,jspd),k3(nvd,jspd)
      integer,intent(in)::k1_b(nvd,jspd),k2_b(nvd,jspd),k3_b(nvd,jspd)
      complex, intent (inout) :: mmn(nbnd,nbnd)

#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
      complex, intent (in):: z(nbasfcn,neigd),z_b(nbasfcn,neigd)
#else
      real,    intent (in):: z(nbasfcn,neigd),z_b(nbasfcn,neigd)
#endif

c     ..basis wavefunctions in the vacuum
      complex, allocatable :: ac(:,:),bc(:,:),ac_b(:,:),bc_b(:,:)
      real,    allocatable :: dt(:),dte(:)
      real,    allocatable :: t(:),te(:),tei(:)
      real,    allocatable :: u(:,:),ue(:,:),v(:)
      real,    allocatable :: dt_b(:),dte_b(:)
      real,    allocatable :: t_b(:),te_b(:),tei_b(:)
      real,    allocatable :: u_b(:,:),ue_b(:,:)

c     ..local scalars..
      logical tail
      real wronk,arg,zks,tpi,vz0(2),scale,evacp,ev,const
      real uu,ud,du,dd,xx(nmz),xximag(nmz)
      real :: uuimag,udimag,duimag,ddimag,qss1,qss2
      integer i,m,l,j,k,n,nv2,nv2_b,ivac,n2,n2_b,sign,ik
      integer :: lprime,np1,addnoco,addnoco2
      complex :: av,bv,ic,c_1
      integer, allocatable :: kvac1(:),kvac2(:),map2(:)
      integer, allocatable :: kvac1_b(:),kvac2_b(:),map2_b(:)
      integer symvac,symvacvac
c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt

c..   external functions
      external vacudz,vacuz
c      real dotirp
c      external dotirp  ! module now

      allocate ( ac(nv2d,nslibd),bc(nv2d,nslibd),
     +           ac_b(nv2d,nslibd_b),bc_b(nv2d,nslibd_b),
     +           dt(nv2d),dte(nv2d),t(nv2d),te(nv2d),
     +           tei(nv2d),u(nmzd,nv2d),ue(nmzd,nv2d),
     +           dt_b(nv2d),dte_b(nv2d),t_b(nv2d),te_b(nv2d),
     +           tei_b(nv2d),u_b(nmzd,nv2d),ue_b(nmzd,nv2d),
     +           v(3),kvac1(nv2d),kvac2(nv2d),map2(nvd),
     +           kvac1_b(nv2d),kvac2_b(nv2d),map2_b(nvd) )

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)

      tail = .true.
      np1 = nmz + 1

c.. determining the indexing array (in-plane stars)
c.. for the k-point

      wronk = 2.0
      const = 1.0 / ( sqrt(omtil)*wronk )

      do ivac = 1,2
         vz0(ivac) = vz(nmz,ivac)
      enddo

      n2 = 0 ; n2_b = 0

      addnoco=0
      addnoco2=0
      if(l_noco.and.jspin.eq.2)then
        addnoco=nv(1)+nlotot
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
        
         IF (n2>nv2d)  CALL juDFT_error("wannier Mmn vac",calledby
     +        ="wann_mmkb_vac")

         kvac1(n2) = k1(k,jspin)
         kvac2(n2) = k2(k,jspin)
         map2(k) = n2
 40   continue

c.. and for the b-point
 
      do 41 k = 1,nv_b(jspin)
         do 31 j = 1,n2_b
            if ( k1_b(k,jspin).eq.kvac1_b(j) .and.
     +          k2_b(k,jspin).eq.kvac2_b(j) ) then
                map2_b(k) = j
                goto 41
             endif
 31      continue
         n2_b = n2_b + 1
     
         IF (n2_b>nv2d)  CALL juDFT_error("wannier Mmn vac",calledby
     +        ="wann_mmkb_vac")

         kvac1_b(n2_b) = k1_b(k,jspin)
         kvac2_b(n2_b) = k2_b(k,jspin)
         map2_b(k) = n2_b
 41   continue

c...cycle by the vacua
      do 140 ivac = 1,nvac



      sign = 3. - 2.*ivac
      evacp = evac(ivac)

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
c-----> construct a and b coefficients for the k-point
       symvacvac=1
       if (nvac==1) symvacvac=2
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
         do n = 1,nslibd
            ac(l,n) = ac(l,n) + z(k+addnoco,n)*av
            bc(l,n) = bc(l,n) + z(k+addnoco,n)*bv
         enddo
      enddo

c...now for the k+b point

      do ik = 1,nv2_b
         v(1) = bkpt_b(1) + kvac1_b(ik) + qss1
         v(2) = bkpt_b(2) + kvac2_b(ik) + qss2
         v(3) = 0.
         ev = evacp - 0.5*dot_product(v,matmul(bbmat,v))
         call vacuz(ev,vz(1,ivac),vz0(ivac),nmz,delz,t_b(ik),dt_b(ik),
     +        u_b(1,ik))
         call vacudz(ev,vz(1,ivac),vz0(ivac),nmz,delz,te_b(ik),
     +        dte_b(ik),tei_b(ik),ue_b(1,ik),dt_b(ik),
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

      do k = 1,nv_b(jspin)
         l = map2_b(k)
         zks = k3_b(k,jspin)*bmat(3,3)*sign
         arg = zks*z1
         c_1 = cmplx(cos(arg),sin(arg)) * const
         av = -c_1 * cmplx( dte_b(l),zks*te_b(l) )
         bv =  c_1 * cmplx( dt_b(l),zks*t_b(l) )
c-----> loop over basis functions
         do n = 1,nslibd_b
            ac_b(l,n) = ac_b(l,n) + z_b(k+addnoco2,n)*av
            bc_b(l,n) = bc_b(l,n) + z_b(k+addnoco2,n)*bv
         enddo
      enddo


      do l = 1,nv2
      do lprime = 1,nv2_b
      if (kvac1(l).eq.(kvac1_b(lprime)-gb(1))
     & .and. kvac2(l).eq.(kvac2_b(lprime)-gb(2)))then
         zks = gb(3)*bmat(3,3)*sign

         do i = 1,nmz
             xx(np1-i) = u(i,l)*u_b(i,lprime)*
     *          cos((z1+(i-1)*delz)*zks)
            xximag(np1-i) = u(i,l)*u_b(i,lprime)*
     *          sin((z1+(i-1)*delz)*zks)
         enddo   
         call intgz0(xx,delz,nmz,uu,tail)
         call intgz0(xximag,delz,nmz,uuimag,tail)

         do i = 1,nmz
            xx(np1-i) = u(i,l)*ue_b(i,lprime)*
     *   cos((z1+(i-1)*delz)*zks)
            xximag(np1-i) = u(i,l)*ue_b(i,lprime)*
     *          sin((z1+(i-1)*delz)*zks)
         enddo   
         call intgz0(xx,delz,nmz,ud,tail)
         call intgz0(xximag,delz,nmz,udimag,tail)

         do i = 1,nmz
            xx(np1-i) = ue(i,l)*u_b(i,lprime)*
     *   cos((z1+(i-1)*delz)*zks)
            xximag(np1-i) = ue(i,l)*u_b(i,lprime)*
     *          sin((z1+(i-1)*delz)*zks)
         enddo   
         call intgz0(xx,delz,nmz,du,tail)
         call intgz0(xximag,delz,nmz,duimag,tail)
         do i = 1,nmz
            xx(np1-i) = ue(i,l)*ue_b(i,lprime)*
     *   cos((z1+(i-1)*delz)*zks)
            xximag(np1-i) = ue(i,l)*ue_b(i,lprime)*
     *          sin((z1+(i-1)*delz)*zks)
         enddo   
         call intgz0(xx,delz,nmz,dd,tail)
         call intgz0(xximag,delz,nmz,ddimag,tail)

         do i = 1,nslibd
            do j = 1,nslibd_b
               mmn(i,j) = mmn(i,j) + area*(
     *  ac(l,i)*conjg(ac_b(lprime,j))*cmplx(uu,uuimag) +
     +  ac(l,i)*conjg(bc_b(lprime,j))*cmplx(ud,udimag) +
     *  bc(l,i)*conjg(ac_b(lprime,j))*cmplx(du,duimag) +
     +  bc(l,i)*conjg(bc_b(lprime,j))*cmplx(dd,ddimag) )
            enddo
         enddo





      endif ! kvac1=kvac1_b and kvac2=kvac2_b
      enddo ! lprime
      enddo ! l

      enddo !symvac
c... cycle by the vacua finishes
 140  enddo      

      deallocate ( ac,bc,dt,dte,t,te,tei,u,ue,
     +             v,kvac1,kvac2,map2 )
      deallocate ( ac_b,bc_b,dt_b,dte_b,t_b,te_b,tei_b,u_b,ue_b,
     +             kvac1_b,kvac2_b,map2_b )

      END SUBROUTINE wann_mmkb_vac
      END MODULE m_wann_mmkb_vac
