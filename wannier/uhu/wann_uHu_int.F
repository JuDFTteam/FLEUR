c*****************************************c
c   Interstitial contribution to uHu      c
c  < u_{k+b1} | H_{k}^{int} | u_{k+b2} >  c
c*****************************************c
c  k1_b,  k2_b,  k3_b : G-vectors at k+b1 c
c  k1_b2, k2_b2, k3_b2: G-vectors at k+b2 c
c                                         c
c  gb : reciprocal latt. vector taking    c
c       k+b1 back to first Brillouin zone c
c  gb2: reciprocal latt. vector taking    c
c       k+b2 back to first Brillouin zone c
c                                         c
c  z_b : eigenvectors at k+b1             c
c  z_b2: eigenvectors at k+b2             c
c*****************************************c
c                J.-P. Hanke, Dec. 2015   c
c*****************************************c
      module m_wann_uHu_int
      contains
      subroutine wann_uHu_int(
     >               chi,nvd,k1d,k2d,k3d,
     >               n3d,nv_b,nv_b2,nbnd,neigd,
     >               nslibd_b,nslibd_b2,nbasfcn,
     >               addnoco,addnoco2,
     >               k1_b ,k2_b ,k3_b, gb,
     >               k1_b2,k2_b2,k3_b2,gb2,
     >               bkpt,bbmat,vpw,zMat_b,zMat_b2,
     >               rgphs,ustep,ig,l_kin,sign,uHu)
#include "cpp_double.h"

      USE m_types

      implicit none

      TYPE(t_zmat), INTENT(IN) :: zMat_b, zMat_b2

c     ..arguments..
      logical, intent(in) :: l_kin
      integer, intent(in) :: gb(3),gb2(3)
      integer, intent(in) :: k1d,k2d,k3d
      integer, intent(in) :: nvd,n3d,nv_b,nv_b2,nbnd,neigd
      integer, intent(in) :: nslibd_b,nslibd_b2,nbasfcn
      integer, intent(in) :: addnoco,addnoco2,sign
      integer, intent(in) :: k1_b(nvd) ,k2_b(nvd) ,k3_b(nvd)
      integer, intent(in) :: k1_b2(nvd),k2_b2(nvd),k3_b2(nvd)
      integer, intent(in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      complex, intent(in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      real,    intent(in) :: bkpt(3),bbmat(3,3)
      complex, intent(in) :: ustep(n3d),vpw(n3d)
      complex, intent(in) :: chi

      complex, intent(inout) :: uHu(nbnd,nbnd)

c     ..local variables..
      complex,allocatable :: vstep_c(:,:)
      complex,allocatable :: mat_c(:,:)
      complex :: th_c
      real,allocatable :: vstep_r(:,:)
      real,allocatable :: mat_r(:,:)
      real,allocatable :: uHu_tmp(:,:)
      real :: th_r

      real,allocatable :: rk_b(:),rk_b2(:)
      real :: phase,phase2,ekk,s(3)
      integer :: i,j,j1,j2,j3,i1,i2,i3,in,ind

      allocate( rk_b (nv_b)  )
      allocate( rk_b2(nv_b2) )

      IF(zMat_b%l_real) THEN
         allocate( uHu_tmp(nslibd_b,nslibd_b2) )
         allocate( vstep_r(nv_b2,nv_b) )
         allocate( mat_r(nv_b,nslibd_b2) )
         vstep_r(:,:) = 0.0
      ELSE
         allocate( vstep_c(nv_b2,nv_b) )
         allocate( mat_c(nv_b,nslibd_b2) )
         vstep_c(:,:) = CMPLX(0.0,0.0)
      END IF

      ! set up |k+G-G(k+b1)|^2
      do i=1,nv_b
       s(1) = bkpt(1) + k1_b(i) - gb(1)
       s(2) = bkpt(2) + k2_b(i) - gb(2)
       s(3) = bkpt(3) + k3_b(i) - gb(3)
       rk_b(i) = dot_product(s,matmul(bbmat,s))
!       rk_b(i) = dotirp(s,s,bbmat)
      enddo

      ! set up |k+G'-G(k+b2)|^2
      do i=1,nv_b2
       s(1) = bkpt(1) + k1_b2(i) - gb2(1)
       s(2) = bkpt(2) + k2_b2(i) - gb2(2)
       s(3) = bkpt(3) + k3_b2(i) - gb2(3)
       rk_b2(i) = dot_product(s,matmul(bbmat,s))
!       rk_b2(i) = dotirp(s,s,bbmat)
      enddo

      ! construct vstep(g,g') ~ V(g-g') 
      !                       + Theta(g-g')*[rk_b+rk_b2]
      do i=1,nv_b
       j1 = -k1_b(i) + gb(1) - gb2(1)
       j2 = -k2_b(i) + gb(2) - gb2(2)
       j3 = -k3_b(i) + gb(3) - gb2(3)
       do j=1,nv_b2
        i1 = j1 + k1_b2(j)
        i2 = j2 + k2_b2(j)
        i3 = j3 + k3_b2(j)
        in = ig(sign*i1,sign*i2,sign*i3)
        if(in.eq.0) cycle
        phase = rgphs(i1,i2,i3)    ! TODO: sign also here?
        phase2= rgphs(sign*i1,sign*i2,sign*i3)
        if(phase.ne.phase2) then
           stop 'rgphs in wann_uHu_int'
        endif

        ekk = rk_b(i) + rk_b2(j)

        IF(zMat_b%l_real) THEN
           th_r = phase*real(vpw(in))
           if(l_kin) th_r = th_r + phase*0.25*ekk*real(ustep(in))
           vstep_r(j,i) = th_r
        ELSE
           th_c = phase*conjg(vpw(in))
           if(l_kin) th_c = th_c + phase*0.25*ekk*conjg(ustep(in))
           vstep_c(j,i) = th_c
        END IF
       enddo
      enddo

      ! complex conjugate of (z(k+b1,g))^* vstep(g,g') z(k+b2,g')
      IF(zMat_b%l_real) THEN
         call CPP_BLAS_sgemm('T','N',nv_b,nslibd_b2,nv_b2,real(1.0),
     >               vstep_r,nv_b2,zMat_b2%z_r(1+addnoco2,1),nbasfcn,
     >               real(0.0),mat_r,nv_b)
         call CPP_BLAS_sgemm('T','N',nslibd_b,nslibd_b2,nv_b,
     >              real(1.0),zMat_b%z_r(1+addnoco,1),nbasfcn,
     >              mat_r,nv_b,real(0.0),uHu_tmp,nslibd_b)
         uHu(1:nslibd_b,1:nslibd_b2) = uHu(1:nslibd_b,1:nslibd_b2)
     >                            + uHu_tmp(1:nslibd_b,1:nslibd_b2)*chi
      ELSE
         call CPP_BLAS_cgemm('T','N',nv_b,nslibd_b2,nv_b2,cmplx(1.0),
     >               vstep_c,nv_b2,zMat_b2%z_c(1+addnoco2,1),nbasfcn,
     >               cmplx(0.0),mat_c,nv_b)
         mat_c = conjg(mat_c)
         call CPP_BLAS_cgemm('T','N',nslibd_b,nslibd_b2,nv_b,
     >              chi,zMat_b%z_c(1+addnoco,1),nbasfcn,
     >              mat_c,nv_b,cmplx(1.0),uHu,nbnd)
      END IF

      deallocate( rk_b, rk_b2 )

      IF(zMat_b%l_real) THEN
         deallocate( uHu_tmp )
         deallocate( vstep_r, mat_r )
      ELSE
         deallocate( vstep_c, mat_c )
      END IF

      end subroutine
      end module m_wann_uHu_int
