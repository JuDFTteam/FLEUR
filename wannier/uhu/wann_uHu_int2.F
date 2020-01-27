      module m_wann_uHu_int2
      contains
      subroutine wann_uHu_int2(
     >               nvd,k1d,k2d,k3d,
     >               n3d,k1,k2,k3,
     >               nv,
     >               k1_b,k2_b,k3_b,
     >               nv_b,
     >               rgphs,ustep,ig,gb,uHu)
#include "cpp_double.h"
      implicit none
      integer, intent(in) :: nvd,n3d,k1(nvd),k2(nvd),k3(nvd)
      integer, intent(in) :: nv
      integer, intent(in) :: nv_b,k1_b(nvd),k2_b(nvd),k3_b(nvd)
      integer, intent(in) :: k1d,k2d,k3d
      real,    intent(in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      complex, intent(in) :: ustep(n3d)
      integer, intent(in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent(in) :: gb(3)
      complex, intent(inout) :: uHu(nvd,nvd)!(nv,nv_b)

#if (!defined(CPP_INVERSION)||defined(CPP_SOC))
      complex :: stepf(nv_b,nv)
      complex phasust
#else
      real :: stepf(nv_b,nv)
      real phasust
#endif
      integer i,j1,j2,j3,i1,i2,i3,j,in,m,n
      real phase

      stepf(:,:)=0.0
      do i =1,nv
       do j = 1,nv_b
c-->     determine index and phase factor
         i1 = k1_b(j) - k1(i) - gb(1)
         i2 = k2_b(j) - k2(i) - gb(2)
         i3 = k3_b(j) - k3(i) - gb(3)
         in = ig(i1,i2,i3)
         if (in.eq.0) cycle
         phase = rgphs(i1,i2,i3)
#if ( !defined(CPP_INVERSION) || defined(CPP_SOC) )
         phasust=conjg(phase*ustep(in))
#else
         phasust=phase*real(ustep(in))
#endif
         stepf(j,i)=phasust
       enddo
      enddo
      uHu(1:nv,1:nv_b) = uHu(1:nv,1:nv_b) + conjg(transpose(stepf))

      end subroutine
      end module m_wann_uHu_int2
