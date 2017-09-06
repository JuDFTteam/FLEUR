!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_mmkb_int
      contains
      subroutine wann_mmkb_int(
     >               interchi,addnoco,addnoco2,nvd,k1d,k2d,k3d,
     >               n3d,k1,k2,k3,
     >               nv,neigd,nbasfcn,zMat,nslibd,
     >               k1_b,k2_b,k3_b,
     >               nv_b,zMat_b,nslibd_b,
     >               nbnd,rgphs,ustep,ig,gb,
     <               mmnk)

      USE m_types

#include "cpp_double.h"
      implicit none

      TYPE(t_zmat), INTENT(IN) :: zMat, zMat_b

      integer, intent(in) :: addnoco,addnoco2
      integer, intent(in) :: nvd,n3d,k1(nvd),k2(nvd),k3(nvd)
      integer, intent(in) :: nv,neigd,nbasfcn,nslibd,nslibd_b
      integer, intent(in) :: nv_b,k1_b(nvd),k2_b(nvd),k3_b(nvd)
      integer, intent(in) :: nbnd
      integer, intent(in) :: k1d,k2d,k3d
      complex, intent(in) :: rgphs(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      complex, intent(in) :: ustep(n3d),interchi
      integer, intent(in) :: ig(-k1d:k1d,-k2d:k2d,-k3d:k3d)
      integer, intent(in) :: gb(3)
      complex, intent(inout) :: mmnk(nbnd,nbnd)

      complex,allocatable::stepf_c(:,:)
      complex,allocatable::phasusbmat_c(:,:)
      real,allocatable::stepf_r(:,:)
      real,allocatable::phasusbmat_r(:,:)
      real,allocatable::mmnk_tmp(:,:)
      integer i,j1,j2,j3,i1,i2,i3,j,in,m,n
      real phase

      IF(zMat%l_real) THEN
         allocate(mmnk_tmp(nslibd,nslibd_b))
         allocate(phasusbmat_r(nv,nslibd_b))
         allocate(stepf_r(nv_b,nv))
         stepf_r = 0.0
      ELSE
         allocate(phasusbmat_c(nv,nslibd_b))
         allocate(stepf_c(nv_b,nv))
         stepf_c = CMPLX(0.0,0.0)
      END IF

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
         IF (zMat%l_real) THEN
            stepf_r(j,i) = phase*real(ustep(in))
         ELSE
            stepf_c(j,i) = conjg(phase*ustep(in))
         END IF
       enddo
      enddo
      IF(zMat%l_real) THEN
         call CPP_BLAS_sgemm('T','N',nv,nslibd_b,nv_b,real(1.0),
     &               stepf_r,nv_b,zMat_b%z_r(1+addnoco2,1),nbasfcn,
c     &               stepf_r,nv_b,zMat_b%z_r,nbasfcn,
     &               real(0.0),phasusbmat_r,nv)
         call CPP_BLAS_sgemm('T','N',nslibd,nslibd_b,nv,real(1.0),
     &               zMat%z_r(1+addnoco,1),nbasfcn,phasusbmat_r,nv,
c     &               zMat%z_r,nbasfcn,phasusbmat_r,nv,
     &               real(0.0),mmnk_tmp,nbnd)
         mmnk(1:nslibd,1:nslibd_b)=mmnk(1:nslibd,1:nslibd_b)+
     &               mmnk_tmp(1:nslibd,1:nslibd_b)*interchi
      ELSE
         call CPP_BLAS_cgemm('T','N',nv,nslibd_b,nv_b,cmplx(1.0),
     &               stepf_c,nv_b,zMat_b%z_c(1+addnoco2,1),
c     &               stepf_c,nv_b,zMat_b%z_c,
     &               nbasfcn,cmplx(0.0),
     &               phasusbmat_c,nv)
         phasusbmat_c=conjg(phasusbmat_c)            
         call CPP_BLAS_cgemm('T','N',nslibd,nslibd_b,nv,interchi, 
     &               zMat%z_c(1+addnoco,1),nbasfcn,phasusbmat_c,nv,
c     &               zMat%z_c,nbasfcn,phasusbmat_c,nv,
     &               cmplx(1.0),mmnk,nbnd)
      END IF

      IF(zMat%l_real) THEN
         deallocate(mmnk_tmp)
         deallocate(phasusbmat_r,stepf_r)
      ELSE
         deallocate(phasusbmat_c,stepf_c)
      END IF

      end subroutine
      end module m_wann_mmkb_int
