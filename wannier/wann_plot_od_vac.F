      module m_wann_plot_od_vac
      contains
      subroutine wann_plot_od_vac(point,
     >     z1,nmzd,nv2d,odi,nmz,delz,bmat,
     >     bkpt,nvd,nv,omtil,k3,ac,bc,u,ue,
     <     value)
c**************************************************************
c   Calculates the lattice periodic part of the Bloch function
c   in the vacuum. Used for plotting.  
c   FF, Sep. '06  ---->  1D version YM January 07  
c*************************************************************** 
      use m_constants, only : pimach
      use m_angle 
      use m_types
      implicit none
c     .. scalar Arguments..
      integer, intent (in) :: nmzd,nv2d,nmz,nvd,nv
      real,    intent (in) :: delz,z1,omtil,point(3)

c     ..array arguments..
      integer, intent (in) :: k3(nvd)
      real,    intent (in) :: bkpt(3),bmat(3,3)
      complex, intent (out) :: value
      type (od_inp), intent (in) :: odi

c     ..basis wavefunctions in the vacuum
      complex, intent(in) :: ac(nv2d,-odi%mb:odi%mb)
      complex, intent(in) :: bc(nv2d,-odi%mb:odi%mb)
      real,    intent(in) :: u(nmzd,nv2d,-odi%mb:odi%mb)
      real,    intent(in) :: ue(nmzd,nv2d,-odi%mb:odi%mb)


c     ..local scalars..
      real arg,zks,tpi,r,phi
      integer i,m,l,j,k,n,nv2,n2
      integer np1
      complex c_1
      integer, allocatable :: kvac3(:),map1(:)
      complex value1
c     ..intrinsic functions..
      intrinsic aimag,cmplx,conjg,real,sqrt

      allocate (kvac3(nv2d),map1(nvd))

      tpi = 2 * pimach() 

      np1 = nmz + 1

c.. determining the indexing array (in-plane stars)
c.. for the k-point

      n2 = 0
      do 35 k = 1,nv
         do 45 j = 1,n2
            if (k3(k).eq.kvac3(j)) then
               map1(k) = j
               goto 35
            end if
 45      continue
         n2 = n2 + 1
         if (n2.gt.nv2d) stop 'wann_plot:vac'
         kvac3(n2) =  k3(k)
         map1(k) = n2
 35   continue
      nv2 = n2

c.. the body of the routine

      value = cmplx(0.0,0.0)
      value1 = cmplx(0.0,0.0)
      r = sqrt(point(1)**2+point(2)**2)
      phi = angle(point(1),point(2))

c      print*,"difference=",(abs(point(3))-z1)/delz

         i=(r-z1)/delz + 1

         if (i.gt.nmz) then
            i=nmz
            print*,"i.gt.nmz in wann_plot_od_vac"
         endif

      do l = 1,nv2  !calculation for i
        do m = -odi%mb,odi%mb
          arg=kvac3(l)*bmat(3,3)*point(3) + m*phi
          c_1=cmplx(cos(arg),sin(arg))
          value = value + (u(i,l,m)*ac(l,m)+ue(i,l,m)*bc(l,m))*c_1
c         print*,"value=",value
          if (real(value).gt.10.or.real(value).lt.-10)then
          print*,"value=",value
          print*,"i=",i
          print*,"u(i,l)=",u(i,l,m)
          print*,"ac(l)=",ac(l,m)
          print*,"bc(l)=",bc(l,m)
          print*,"ue(i,l)=",ue(i,l,m)
          endif
        enddo !m
      enddo ! l

      i=i+1

      do l = 1,nv2  !calculation for i
        do m = -odi%mb,odi%mb
          arg=(kvac3(l)*bmat(3,3))*point(3) + m*phi
          c_1=cmplx(cos(arg),sin(arg))
          value1 = value1 + (u(i,l,m)*ac(l,m)+ue(i,l,m)*bc(l,m))*c_1
        enddo !m
      enddo ! l

c     value=(value1-value)*(r/delz+2-i) + value
    
      deallocate ( kvac3,map1 )

      end subroutine wann_plot_od_vac
      end module m_wann_plot_od_vac
