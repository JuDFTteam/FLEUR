!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_plot_vac
      use m_juDFT
c**************************************************************
c   Calculates the lattice periodic part of the Bloch function
c   in the vacuum. Used for plotting.  
c   FF, Sep. '06    
c        
c*************************************************************** 
      CONTAINS
      SUBROUTINE wann_plot_vac(point,
     >     z1,nmzd,nv2d,n3d,nvac,
     >     nmz,delz,bmat,bbmat,
     >     evac,bkpt,vz,
     >     jspin,
     >     k1,k2,k3,nvd,
     >     nbasfcn,neigd,nv,omtil,nslibd,ac,bc,u,ue,
     <     value)

      use m_constants, only : pimach
      implicit none

c     .. scalar Arguments..
      integer, intent (in) :: nmzd,nv2d,n3d,nslibd
      integer, intent (in) :: nmz,nvac
      integer, intent (in) :: jspin,nvd
      integer, intent (in) :: nbasfcn,neigd
      real,    intent (in) :: delz,z1,omtil,point(3)

c     ..array arguments..
      real,    intent (in) :: bkpt(3),evac(2)
      integer, intent (in) :: nv
      real,    intent (in) :: vz(nmzd,2),bmat(3,3),bbmat(3,3)
      integer, intent (in) :: k1(nvd),k2(nvd),k3(nvd)
      complex, intent (out) :: value

c     ..basis wavefunctions in the vacuum
      complex, intent(in) :: ac(nv2d),bc(nv2d)
      real,    intent(in) :: u(nmzd,nv2d),ue(nmzd,nv2d)


c     ..local scalars..
      real wronk,arg,zks,tpi,vz0(2),scale,evacp,ev,const
      real uu,ud,du,dd,xx(nmz)
      integer i,m,l,j,k,n,nv2,nv2_b,ivac,n2,n2_b,sign,ik
      integer lprime,np1
      complex av,bv,ic,c_1
      integer, allocatable :: kvac1(:),kvac2(:),map2(:)
      complex value1
c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt


      allocate (kvac1(nv2d),kvac2(nv2d),map2(nvd))

      tpi = 2 * pimach() ; ic = cmplx(0.,1.)

      np1 = nmz + 1

c.. determining the indexing array (in-plane stars)
c.. for the k-point

      wronk = 2.0
      const = 1.0 / ( sqrt(omtil)*wronk )

      n2 = 0 

      do 40 k = 1,nv
         do 30 j = 1,n2
            if ( k1(k).eq.kvac1(j) .and.
     +          k2(k).eq.kvac2(j) ) then
                map2(k) = j
                goto 40
             endif 
 30      continue
         n2 = n2 + 1
         if (n2>nv2d)  CALL juDFT_error("wann_plot_vac: map",calledby
     +        ="wann_plot_vac")
         kvac1(n2) = k1(k)
         kvac2(n2) = k2(k)
         map2(k) = n2
 40   continue

 
      nv2 = n2 

c.. the body of the routine



      value=cmplx(0.0,0.0)
      value1=cmplx(0.0,0.0)
c      print*,"difference=",(abs(point(3))-z1)/delz
         i=(abs(point(3))-z1)/delz +1
  
 
   
         if (i.gt.nmz) then
            i=nmz
            print*,"i.gt.nmz in wann_plot_vac"
         endif

      do l = 1,nv2  !calculation for i


         arg=(kvac1(l)*bmat(1,1)+kvac2(l)*bmat(2,1))*point(1)+
     +    (kvac1(l)*bmat(1,2)+kvac2(l)*bmat(2,2))*point(2)
         c_1=cmplx(cos(arg),sin(arg))
            value =value+ (u(i,l)*ac(l)+ue(i,l)*bc(l))*c_1
c            print*,"value=",value
       if (real(value).gt.10.or.real(value).lt.-10)then
          print*,"value=",value
          print*,"i=",i
          print*,"u(i,l)=",u(i,l)
          print*,"ac(l)=",ac(l)
          print*,"bc(l)=",bc(l)
          print*,"ue(i,l)=",ue(i,l)
      endif
      enddo ! l

      
      i=i+1
      do l = 1,nv2

         arg=(kvac1(l)*bmat(1,1)+kvac2(l)*bmat(2,1))*point(1)+
     +    (kvac1(l)*bmat(1,2)+kvac2(l)*bmat(2,2))*point(2)
         c_1=cmplx(cos(arg),sin(arg))
            value1 =value1+ (u(i,l)*ac(l)+ue(i,l)*bc(l))*c_1

      enddo ! l

      value=(value1-value)*((abs(point(3))-z1)/delz+2-i)+value
    
      deallocate (kvac1,kvac2,map2 )

      END SUBROUTINE wann_plot_vac
      END MODULE m_wann_plot_vac
