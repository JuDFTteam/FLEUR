c****************************************c
c   Muffin tin contribution to uHu       c
c  < u_{k+b1} | H_{k}^{mt} | u_{k+b2} >  c
c****************************************c
c   Use T(lmp,lm) calculated before for  c
c   every pair (b1,b2) and all atoms n   c
c                                        c
c   acof  ,bcof  ,ccof  : coefficients   c
c                         at k+b1        c
c   acof_b,bcof_b,ccof_b: coefficients   c
c                         at k+b2        c
c                                        c
c   bkpt   : k-point                     c
c   bkpt_b : (k+b1)-point                c
c   bkpt_b2: (k+b2)-point                c
c****************************************c
c               J.-P. Hanke, Dec. 2015   c
c****************************************c
      MODULE m_wann_uHu_sph
      contains
      subroutine wann_uHu_sph(
     >        chi,nbnd,llod,nslibd,nslibd_b,nlod,natd,ntypd,lmd,jmtd,
     >        taual,nop,lmax,
     >        ntype,neq,nlo,llo,acof,bcof,ccof,bkpt_b2,
     >        acof_b,bcof_b,ccof_b,bkpt_b,bkpt,gb_b,gb_b2,
     >        tuu,tud,tdu,tdd,tuulo,tulou,tdulo,tulod,tuloulo,
     >        kdiff,kdiff2,nntot,nntot2,uHu)
#include "cpp_double.h"
      USE m_juDFT
      use m_constants, only : pimach
      use m_matmul   , only : matmul3,matmul3r

      implicit none

c     .. scalar arguments ..
      integer, intent (in) :: llod,nlod,natd,ntypd,lmd,nbnd
      integer, intent (in) :: nntot,nntot2
      integer, intent (in) :: ntype,nslibd,nslibd_b,nop,jmtd
      complex, intent (in) :: chi

c     .. array arguments ..
      integer, intent (in)  :: neq(:) !(ntypd)
      integer, intent (in)  :: lmax(:) !(ntypd)
      integer, intent (in)  :: nlo(:) !(ntypd)
      integer, intent (in)  :: llo(:,:) !(nlod,ntypd)
      real,    intent (in)  :: bkpt_b(:) !(3)
      real,    intent (in)  :: bkpt_b2(:) !(3)
      real,    intent (in)  :: taual(:,:) !(3,natd)
      real,    intent (in)  :: bkpt(:) !(3)
      integer, intent (in)  :: gb_b(:) !(3)
      integer, intent (in)  :: gb_b2(:) !(3)
      complex, intent (in)  :: ccof(-llod:,:,:,:) !(-llod:llod,nslibd,nlod,natd)
      complex, intent (in)  :: acof(:,0:,:) !(nslibd,0:lmd,natd)
      complex, intent (in)  :: bcof(:,0:,:) !(nslibd,0:lmd,natd)
      complex, intent (in)  :: ccof_b(-llod:,:,:,:) !(-llod:llod,nslibd_b,nlod,natd)
      complex, intent (in)  :: acof_b(:,0:,:) !(nslibd_b,0:lmd,natd)
      complex, intent (in)  :: bcof_b(:,0:,:) !(nslibd_b,0:lmd,natd)

      complex, intent (in)  :: tuu(0:,0:,:,:)
      complex, intent (in)  :: tud(0:,0:,:,:)
      complex, intent (in)  :: tdu(0:,0:,:,:)
      complex, intent (in)  :: tdd(0:,0:,:,:)
      complex, intent (in)  :: tuulo(0:,:,-llod:,:,:)
      complex, intent (in)  :: tulou(0:,:,-llod:,:,:)
      complex, intent (in)  :: tdulo(0:,:,-llod:,:,:)
      complex, intent (in)  :: tulod(0:,:,-llod:,:,:)
      complex, intent (in)  :: tuloulo(:,-llod:,:,-llod:,:,:)

      real, intent (in)     :: kdiff(:,:),kdiff2(:,:)
      complex,intent(inout) :: uHu(:,:)


c     .. local scalars ..
      integer i,lm,nn,n,na,j,lmp,l,lp,m,mp,lwn,lo,lop
      integer ll,llp,lmd2
      real rph,cph,tpi,sqpi16,th,t1nn,t2nn,t3nn
      integer nene,nene2,indexx
      complex :: fac1,fac2,fac3,fac4
      complex :: mat(0:lmd,nslibd_b)
C     ..
C     .. local arrays ..
      real bpt(3)
      real bpt2(3)

C     ..
C     .. intrinsic functions ..
      intrinsic conjg,cmplx,sqrt,cos,sin

      tpi = 2* pimach()
      sqpi16 = 4*tpi*tpi
      lmd2 = lmd+1
      
      na = 0

      ! find neighbor k+b1
      bpt(:) = bkpt_b(:) + gb_b(:) - bkpt(:)
      do nene=1,nntot
         if(all(abs(bpt(:)-kdiff(:,nene)).lt.1e-4)) exit
      enddo
      IF(nene==nntot+1) CALL juDFT_error
     +     ("cannot find matching nearest neighbor k+b1",calledby
     +     ="wann_uHu_sph")

      ! find neighbor k+b2
      bpt2(:) = bkpt_b2(:) + gb_b2(:) - bkpt(:)
      do nene2=1,nntot2
         if(all(abs(bpt2(:)-kdiff2(:,nene2)).lt.1e-4)) exit
      enddo  
      IF(nene2==nntot2+1) CALL juDFT_error
     +     ("cannot find matching nearest neighbor k+b2",calledby
     +     ="wann_uHu_sph")

      indexx=nene2+(nene-1)*nntot2

c      if(nene2.ne.1) stop 'nene2.ne.1'
c      if(indexx.ne.nene) stop 'nene.ne.indexx'
c      if(ANY(bpt2.ne.0.0)) stop 'bpt2.ne.0'

      do n=1,ntype
       lwn = lmax(n)
         do nn = 1,neq(n) ! cycle by the atoms within the atom type
         na = na + 1
c...set up phase factors ( 4pi e^{ib2\tau} 4pi e^{-ib1\tau} )

         t1nn =  tpi*taual(1,na)
         t2nn =  tpi*taual(2,na)
         t3nn =  tpi*taual(3,na)

         th = (bpt2(1)-bpt(1))*t1nn
     >       +(bpt2(2)-bpt(2))*t2nn
     >       +(bpt2(3)-bpt(3))*t3nn
         rph = sqpi16*cos(th)
         cph = sqpi16*sin(th)


c...apw-apw
         call CPP_BLAS_cgemm('T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >                    tuu(0,0,n,indexx),lmd2,
     >                    acof_b(1,0,na),nslibd_b,
     >                    cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm('N','N',nslibd,nslibd_b,lmd2,chi,
     >                    acof(1,0,na),nslibd,mat(0,1),lmd2,
     >                    cmplx(1.0),uHu,nbnd)

         call CPP_BLAS_cgemm('T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >                    tud(0,0,n,indexx),lmd2,
     >                    bcof_b(1,0,na),nslibd_b,
     >                    cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm('N','N',nslibd,nslibd_b,lmd2,chi,
     >                    acof(1,0,na),nslibd,mat(0,1),lmd2,
     >                    cmplx(1.0),uHu,nbnd)

         call CPP_BLAS_cgemm('T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >                    tdu(0,0,n,indexx),lmd2,
     >                    acof_b(1,0,na),nslibd_b,
     >                    cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm('N','N',nslibd,nslibd_b,lmd2,chi,
     >                    bcof(1,0,na),nslibd,mat(0,1),lmd2,
     >                    cmplx(1.0),uHu,nbnd)

         call CPP_BLAS_cgemm('T','C',lmd2,nslibd_b,lmd2,cmplx(rph,cph),
     >                    tdd(0,0,n,indexx),lmd2,
     >                    bcof_b(1,0,na),nslibd_b,
     >                    cmplx(0.0),mat(0,1),lmd2)
         call CPP_BLAS_cgemm('N','N',nslibd,nslibd_b,lmd2,chi,
     >                    bcof(1,0,na),nslibd,mat(0,1),lmd2,
     >                    cmplx(1.0),uHu,nbnd)

         if(nlo(n).ge.1) then

c...apw-lo
         do lo = 1,nlo(n)
          l = llo(lo,n)
          do m = -l, l
           
           do lp = 0, lwn
            llp = lp*(lp+1)
            do mp = -lp, lp
             lmp = llp + mp

             fac1=cmplx(rph,cph)*tulou(lmp,lo,m,n,indexx)*chi
             fac2=cmplx(rph,cph)*tulod(lmp,lo,m,n,indexx)*chi
             fac3=cmplx(rph,cph)*tuulo(lmp,lo,m,n,indexx)*chi
             fac4=cmplx(rph,cph)*tdulo(lmp,lo,m,n,indexx)*chi

             do i = 1,nslibd
               do j = 1,nslibd_b
                 uHu(i,j) = uHu(i,j)
     >             + ccof(m,i,lo,na)* fac1 *conjg(acof_b(j,lmp,na))
     >             + ccof(m,i,lo,na)* fac2 *conjg(bcof_b(j,lmp,na))
     >             + acof(i,lmp,na) * fac3 *conjg(ccof_b(m,j,lo,na))
     >             + bcof(i,lmp,na) * fac4 *conjg(ccof_b(m,j,lo,na))
               enddo 
             enddo  

            enddo !mp
           enddo  !lp

          enddo ! m
         enddo  ! lo

c...lo-lo
         do lo = 1,nlo(n)
          l = llo(lo,n)
          do m = -l, l
           
           do lop = 1,nlo(n)
            lp = llo(lop,n)
            do mp = -lp, lp

             fac1=cmplx(rph,cph)*tuloulo(lop,mp,lo,m,n,indexx)*chi

             do i = 1,nslibd
               do j = 1,nslibd_b
                 uHu(i,j) = uHu(i,j)
     >             + ccof(m,i,lo,na)*fac1*conjg(ccof_b(mp,j,lop,na))
               enddo 
             enddo  

            enddo !mp
           enddo  !lop

          enddo ! m
         enddo  ! lo

         endif !(nlo(n).ge.1)

      enddo  ! atoms in the type

      enddo  ! atom type 

      end subroutine wann_uHu_sph
      end module m_wann_uHu_sph
