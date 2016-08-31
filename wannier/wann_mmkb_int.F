!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_mmkb_int
      contains
      subroutine wann_mmkb_int(
     >               addnoco,addnoco2,nvd,k1d,k2d,k3d,
     >               n3d,k1,k2,k3,
     >               nv,neigd,nbasfcn,z,nslibd,
     >               k1_b,k2_b,k3_b,
     >               nv_b,z_b,nslibd_b,
     >               nbnd,rgphs,ustep,ig,gb,
     <               mmnk)
#include "cpp_double.h"
      implicit none
      integer, intent(in) :: addnoco,addnoco2
      integer, intent(in) :: nvd,n3d,k1(nvd),k2(nvd),k3(nvd)
      integer, intent(in) :: nv,neigd,nbasfcn,nslibd,nslibd_b
      integer, intent(in) :: nv_b,k1_b(nvd),k2_b(nvd),k3_b(nvd)
      integer, intent(in) :: nbnd
      integer, intent(in) :: k1d,k2d,k3d
      real,    intent(in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      complex, intent(in) :: ustep(n3d)
      integer, intent(in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent(in) :: gb(3)
#if (!defined(CPP_INVERSION)||defined(CPP_SOC))
      complex, intent(in) :: z(nbasfcn,neigd)
      complex, intent(in) :: z_b(nbasfcn,neigd)
#else
      real,intent(in) :: z(nbasfcn,neigd)
      real,intent(in) :: z_b(nbasfcn,neigd)
#endif      
      complex,intent(inout) :: mmnk(nbnd,nbnd)

#if (!defined(CPP_INVERSION)||defined(CPP_SOC))
      complex,allocatable::stepf(:,:)
      complex phasust
      complex,allocatable::phasusbmat(:,:)
#else
      real,allocatable::stepf(:,:)
      real,allocatable::phasusbmat(:,:)
      real phasust
      real,allocatable::mmnk_tmp(:,:)
#endif
      integer i,j1,j2,j3,i1,i2,i3,j,in,m,n
      real phase

      allocate(stepf(nv_b,nv))
      allocate(phasusbmat(nv,nslibd_b))
#if (!defined(CPP_INVERSION)||defined(CPP_SOC))
#else
      allocate(mmnk_tmp(nslibd,nslibd_b))
#endif
      stepf(:,:)=0.0
      do i =1,nv
       j1 =-k1(i)  - gb(1)
       j2 =-k2(i)  - gb(2)
       j3 =-k3(i)  - gb(3)
       do j = 1,nv_b
c-->     determine index and phase factor
         i1 =j1 + k1_b(j) 
         i2 =j2 + k2_b(j) 
         i3 =j3 + k3_b(j) 
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
#if (!defined(CPP_INVERSION)||defined(CPP_SOC))
      call CPP_BLAS_cgemm('T','N',nv,nslibd_b,nv_b,cmplx(1.0),
     &            stepf,nv_b,z_b(1+addnoco2,1),
c     &            stepf,nv_b,z_b,
     &            nbasfcn,cmplx(0.0),
     &            phasusbmat,nv)
      phasusbmat=conjg(phasusbmat)
      call CPP_BLAS_cgemm('T','N',nslibd,nslibd_b,nv,cmplx(1.0),
     &            z(1+addnoco,1),nbasfcn,phasusbmat,nv,cmplx(1.0),
c     &            z,nbasfcn,phasusbmat,nv,cmplx(1.0),
     &            mmnk,nbnd)
#else
      call CPP_BLAS_sgemm('T','N',nv,nslibd_b,nv_b,real(1.0),
     &            stepf,nv_b,z_b(1+addnoco2,1),nbasfcn,real(0.0),
c     &            stepf,nv_b,z_b,nbasfcn,real(0.0),
     &            phasusbmat,nv)
      call CPP_BLAS_sgemm('T','N',nslibd,nslibd_b,nv,real(1.0),
     &            z(1+addnoco,1),nbasfcn,phasusbmat,nv,real(0.0),
c     &            z,nbasfcn,phasusbmat,nv,real(0.0),
     &            mmnk_tmp,nbnd)
      mmnk(1:nslibd,1:nslibd_b)=mmnk(1:nslibd,1:nslibd_b)+
     &            mmnk_tmp(1:nslibd,1:nslibd_b)
#endif


      deallocate(stepf)
      deallocate(phasusbmat)
#if (!defined(CPP_INVERSION)||defined(CPP_SOC))
#else
      deallocate(mmnk_tmp)
#endif
      end subroutine
      end module m_wann_mmkb_int
