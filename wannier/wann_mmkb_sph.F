!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_mmkb_sph
      contains
      subroutine wann_mmkb_sph(
c***********************************************************************
c   computes the Mmn(k,b) matrix elements which are the overlaps
c   between the Bloch wavefunctions at the state m and n
c   at the k-point and its nearest neighbor b (bpt), in the spheres
c                                Y.Mokrousov 15.6.06
c***********************************************************************
c   the gaunt coefficients impose certain conditions (which can be
c   seen in the gaunt.F) on the l,lp,lpp angular momenta, which can be
c   used in order to speedup the procedure. 
c   The conditions are the following:
c   1. |l - l''| <= l' <= l + l''
c   2. m' = m + m''
c   3. l + l' + l''  - is even
c*********************************************************************** 
c   modified for use in conjunction with wann_ujugaunt, which
c   calculates certain matrix elements that have to be treated only
c   once per calculation
c   Frank Freimuth, October 2006
c***********************************************************************
c   use BLAS matrix-matrix multiplication in apw-apw part to be faster
c   J.-P. Hanke, Dec. 2015
c***********************************************************************
     >        nbnd,llod,nslibd,nslibd_b,nlod,natd,ntypd,lmd,jmtd,
     >        taual,nop,lmax,
     >        ntype,neq,nlo,llo,acof,bcof,ccof,bbpt,
     >        acof_b,bcof_b,ccof_b,gb,bkpt,ujug,ujdg,djug,djdg,ujulog,
     >        djulog,ulojug,ulojdg,ulojulog,kdiff,nntot,chi,
     =        mmn)
#include "cpp_double.h"
      use m_juDFT
      use m_constants, only : pimach
      use m_matmul   , only : matmul3,matmul3r
      use m_sphbes
      use m_ylm
      use m_intgr, only : intgr3
      use m_gaunt, only: gaunt1
    
      implicit none

c     .. scalar arguments ..
      integer, intent (in) :: llod,nlod,natd,ntypd,lmd,nbnd
      integer, intent (in) :: nntot
      integer, intent (in) :: ntype,nslibd,nslibd_b,nop,jmtd

c     .. array arguments ..
      integer, intent (in)  :: neq(:) !(ntypd)
      integer, intent (in)  :: lmax(:) !(ntypd)
      integer, intent (in)  :: nlo(:) !(ntypd)
      integer, intent (in)  :: llo(:,:) !(nlod,ntypd)
      complex, intent (in)  :: chi(ntypd)
      real,    intent (in)  :: bbpt(:) !(3)
      real,    intent (in)  :: taual(:,:) !(3,natd)
      real,    intent (in)  :: bkpt(:) !(3)
      integer, intent (in)  :: gb(:) !(3)
      complex, intent (in)  :: ccof(-llod:,:,:,:) !(-llod:llod,nslibd,nlod,natd)
      complex, intent (in)  :: acof(:,0:,:) !(nslibd,0:lmd,natd)
      complex, intent (in)  :: bcof(:,0:,:) !(nslibd,0:lmd,natd)
      complex, intent (in)  :: ccof_b(-llod:,:,:,:) !(-llod:llod,nslibd_b,nlod,natd)
      complex, intent (in)  :: acof_b(:,0:,:) !(nslibd_b,0:lmd,natd)
      complex, intent (in)  :: bcof_b(:,0:,:) !(nslibd_b,0:lmd,natd)

      complex, intent (in)  :: ujug(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (in)  :: ujdg(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (in)  :: djug(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (in)  :: djdg(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (in)  :: ujulog(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (in)  :: djulog(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (in)  :: ulojug(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (in)  :: ulojdg(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (in)  :: ulojulog(:,-llod:,:,-llod:,:,:) !(1:nlod,-llod:llod,1:nlod,-llod:llod,1:ntype,1:nntot)

      real, intent (in)     :: kdiff(:,:) !(3,nntot)
      complex,intent(inout) :: mmn(:,:)!nbnd,nbnd) !(nbnd,nbnd)


c     .. local scalars ..
      integer i,lm,nn,n,na,j,lmp,l,lp,m,mp,lwn,lo,lop
      integer ll,llp
      real rph,cph,tpi,th,t1nn,t2nn,t3nn,dummy
      complex ic
      integer nene,lmd2
      complex :: fac1,fac2,fac3,fac4
      complex :: fac5,fac6,fac7,fac8
      complex :: mat(0:lmd,nslibd_b)
C     ..
C     .. local arrays ..
      real bpt(3)

C     ..
C     .. intrinsic functions ..
      intrinsic conjg,cmplx,sqrt,cos,sin

      ic = cmplx(0.,1.)
      tpi = 2* pimach()
      lmd2 = lmd+1
      
      na = 0

      bpt(:) = bbpt(:) + gb(:) - bkpt(:)
      do nene=1,nntot
         if(all(abs(bpt(:)-kdiff(:,nene)).lt.1e-4)) exit
      enddo  
      IF(nene==nntot+1) CALL juDFT_error
     +     ("cannot find matching nearest neighbor k",calledby
     +     ="wann_mmkb_sph")

      do n=1,ntype
       lwn = lmax(n)
         do nn = 1,neq(n) ! cycle by the atoms within the atom type
         na = na + 1
c...set up phase factors ( e^{ib\tau} )

         t1nn =  tpi*taual(1,na)
         t2nn =  tpi*taual(2,na)
         t3nn =  tpi*taual(3,na)

         th = bpt(1)*t1nn + bpt(2)*t2nn + bpt(3)*t3nn
         rph = 2*tpi*cos(th)
         cph = 2*tpi*sin(th)

c...contributions from the apws
         call CPP_BLAS_cgemm(
     >          'T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >          ujug(0,0,n,nene),lmd2,acof_b(1,0,na),nslibd_b,
     >          cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm(
     >          'N','N',nslibd,nslibd_b,lmd2,chi(n),
     >          acof(1,0,na),nslibd,mat(0,1),lmd2,
     >          cmplx(1.0),mmn,nbnd)

         call CPP_BLAS_cgemm(
     >          'T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >          ujdg(0,0,n,nene),lmd2,bcof_b(1,0,na),nslibd_b,
     >          cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm(
     >          'N','N',nslibd,nslibd_b,lmd2,chi(n),
     >          acof(1,0,na),nslibd,mat(0,1),lmd2,
     >          cmplx(1.0),mmn,nbnd)

         call CPP_BLAS_cgemm(
     >          'T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >          djug(0,0,n,nene),lmd2,acof_b(1,0,na),nslibd_b,
     >          cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm(
     >          'N','N',nslibd,nslibd_b,lmd2,chi(n),
     >          bcof(1,0,na),nslibd,mat(0,1),lmd2,
     >          cmplx(1.0),mmn,nbnd)

         call CPP_BLAS_cgemm(
     >          'T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >          djdg(0,0,n,nene),lmd2,bcof_b(1,0,na),nslibd_b,
     >          cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm(
     >          'N','N',nslibd,nslibd_b,lmd2,chi(n),
     >          bcof(1,0,na),nslibd,mat(0,1),lmd2,
     >          cmplx(1.0),mmn,nbnd)



         if (nlo(n).ge.1) then 
         
         do lo = 1,nlo(n) 
          l = llo(lo,n)
          ll = l*(l+1)
          do m = -l,l
           lm = ll + m

           do lp = 0,lwn
            llp = lp*(lp+1)
            do mp = -lp,lp
             lmp = llp + mp

             fac1=cmplx(rph,cph)*ujulog(lmp,lo,m,n,nene)*chi(n)
             fac2=cmplx(rph,cph)*djulog(lmp,lo,m,n,nene)*chi(n)
             fac3=cmplx(rph,cph)*ulojdg(lmp,lo,m,n,nene)*chi(n)
             fac4=cmplx(rph,cph)*ulojug(lmp,lo,m,n,nene)*chi(n)
             
c             if(fac1.ne.0)write(*,*)'fac1',lmp,lm
c             if(fac2.ne.0)write(*,*)'fac2',lmp,lm
c             if(fac3.ne.0)write(*,*)'fac3',lmp,lm
c             if(fac4.ne.0)write(*,*)'fac4',lmp,lm

                do j = 1,nslibd_b
                 do i = 1,nslibd
                  mmn(i,j) = mmn(i,j) +
     +                 ccof(m,i,lo,na)*
     *              (conjg(acof_b(j,lmp,na))*fac1+
     +               conjg(bcof_b(j,lmp,na))*fac2 )+

     +                 conjg(ccof_b(m,j,lo,na))*
     *              ( bcof(i,lmp,na)*fac3 +
     +                acof(i,lmp,na)*fac4 )
                 enddo 
                enddo  


            enddo ! mp
           enddo  ! lp

          enddo ! lo
         enddo  ! m lo

c...contributions from lo*lo
          
         do lo = 1,nlo(n)
          l = llo(lo,n)
          ll = l*(l+1)
          do m = -l,l
           lm = ll + m

           do lop = 1,nlo(n)
            lp = llo(lop,n)
            llp = lp*(lp+1)
            do mp = -lp,lp
             lmp = llp + mp

             fac1=cmplx(rph,cph)*ulojulog(lop,mp,lo,m,n,nene)*chi(n)

c             if(fac1.ne.0)write(*,*)'fac1',lop,lo,mp,m

                do j = 1,nslibd_b
                 do i = 1,nslibd

                  mmn(i,j) = mmn(i,j) +
     +               ccof(m,i,lo,na)*conjg(ccof_b(mp,j,lop,na))*
     *                        fac1

                 enddo 
                enddo  
                
            enddo ! mp lop
           enddo  ! lop

          enddo ! m lo
         enddo  ! lo           
         
         endif ! local orbitals on this atom

      enddo  ! atoms in the type

      enddo  ! atom type   


      end subroutine wann_mmkb_sph
      end module m_wann_mmkb_sph
