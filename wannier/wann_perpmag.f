!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


      module m_wann_perpmag
      contains
      SUBROUTINE wann_perpmag(
     >               nbnd,neq,mlh,nlhd,nlh,ntypsy,llh,lmax,
     >               nmem,ntype,ntypd,bbmat,bmat,
     >               nlod,llod,nlo,llo,flo,f,g,jri,rmsh,dx,jmtd,
     >               lmaxd,lmd,clnu,
     >               ujug,ujdg,djug,djdg,ujulog,djulog,
     >               ulojug,ulojdg,ulojulog,
     >               acof,bcof,ccof,
     <               perpmag)
c*************************************************************************
c    wann_perpmag calculates integrals of radial wave functions with
c    the exchange field and multiplies them with an angular factor.
c
c    Frank Freimuth, 2010
c************************************************************************* 

      use m_constants, only : pimach
      use m_matmul   , only : matmul3,matmul3r
      use m_sphbes
      use m_ylm
      use m_types, only : od_inp, od_sym
      use m_intgr, only : intgr3
      use m_gaunt, only: gaunt1

      IMPLICIT NONE
      integer, intent (in) :: nbnd
      integer, intent (in) :: neq(:) 
      INTEGER, INTENT (IN) :: mlh(:,0:,:)!(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN) :: nlhd
      integer, intent (in) :: nlh(:)
      integer, intent (in) :: ntypsy(:)
      integer, intent (in) :: nmem(:,:)
      integer, intent (in) :: ntype,ntypd
      integer, intent (in) :: llh(:,:)
      INTEGER, INTENT (IN) :: lmaxd,jmtd,lmd
      real,    intent (in) :: bbmat(:,:)!(3,3)
      real,    intent (in) :: bmat(:,:)!(3,3)
      integer, intent (in) :: lmax(:)!(ntypd)
      integer, intent (in) :: nlod,llod
      integer, intent (in) :: jri(:)!(ntypd)
      integer, intent (in) :: nlo(:)!(ntypd)
      integer, intent (in) :: llo(:,:)!(nlod,ntypd)    
      real,    intent (in) :: f(:,:,:,0:,:)!(ntypd,jmtd,2,0:lmaxd,2)
      real,    intent (in) :: g(:,:,:,0:,:)!(ntypd,jmtd,2,0:lmaxd,2)
      real,    intent (in) :: flo(:,:,:,:,:)!(ntypd,jmtd,2,nlod,2)
      real,    intent (in) :: rmsh(:,:)!(jmtd,ntypd)
      real,    intent (in) :: dx(:)!(ntypd)
      complex, intent (in) :: clnu(:,0:,:)!(memd,0:nlhd,ntypsd)
      complex, intent (in) :: ujug(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)
      complex, intent (in) :: ujdg(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)
      complex, intent (in) :: djug(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)
      complex, intent (in) :: djdg(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)

      complex, intent (in) :: ujulog(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (in) :: djulog(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (in) :: ulojug(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (in) :: ulojdg(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (in) :: ulojulog(:,-llod:,:,-llod:,:) !(1:nlod,-llod:llod,1:nlod,-llod:llod,1:ntype)

      complex, intent (in) :: acof(:,0:,:,:) !acof(noccbd,0:lmd,natd,jspd)
      complex, intent (in) :: bcof(:,0:,:,:) !bcof(noccbd,0:lmd,natd,jspd)
      complex, intent (in) :: ccof(-llod:,:,:,:,:)!ccof(-llod:llod,noccbd,nlod,natd,jspd)

      complex, intent (inout):: perpmag(:,:)

      integer :: nat,n,nn,l,m,lm,lp,mp,lpmp,i,j,lo,lop

      nat=0
      do n=1,ntype
       do nn=1,neq(n)
        nat=nat+1
        do l=0,lmax(n)
         do m=-l,l  
          lm=l*(l+1)+m  
          do lp=0,lmax(n)
           do mp=-lp,lp
            lpmp=lp*(lp+1)+mp              
            do i=1,nbnd
             do j=1,nbnd
              perpmag(j,i)=perpmag(j,i)+
     + conjg(acof(j,lpmp,nat,1))*acof(i,lm,nat,2)*ujug(lpmp,lm,n)+  
     + conjg(acof(j,lpmp,nat,1))*bcof(i,lm,nat,2)*ujdg(lpmp,lm,n)+
     + conjg(bcof(j,lpmp,nat,1))*acof(i,lm,nat,2)*djug(lpmp,lm,n)+
     + conjg(bcof(j,lpmp,nat,1))*bcof(i,lm,nat,2)*djdg(lpmp,lm,n)
             enddo !j
            enddo !i
           enddo !mp
          enddo !lp
         enddo !m
        enddo !l

        if(nlo(n).ge.1)then
         do lo=1,nlo(n)
          l=llo(lo,n)
          do m=-l,l
           
           do lop=1,nlo(n)
            lp=llo(lop,n)
            do mp=-lp,lp

             do i=1,nbnd
              do j=1,nbnd
               perpmag(j,i)=perpmag(j,i)+
     &            conjg(ccof(mp,j,lop,nat,1))*ccof(m,i,lo,nat,2)*
     &                    ulojulog(lop,mp,lo,m,n)
              enddo !j
             enddo !i

            enddo !mp
           enddo !lop   

          enddo !m
         enddo !lo 
        endif !lo-lo

        if(nlo(n).ge.1)then
         do lo=1,nlo(n)
          l=llo(lo,n)
          do m=-l,l
           
           do lp=0,lmax(n)
            do mp=-lp,lp
            lpmp=lp*(lp+1)+mp   

             do i=1,nbnd
              do j=1,nbnd
               perpmag(j,i)=perpmag(j,i)+
     &            conjg(acof(j,lpmp,nat,1))*ccof(m,i,lo,nat,2)*
     &                    ujulog(lpmp,lo,m,n)
               perpmag(j,i)=perpmag(j,i)+
     &            conjg(bcof(j,lpmp,nat,1))*ccof(m,i,lo,nat,2)*
     &                    djulog(lpmp,lo,m,n)

               perpmag(j,i)=perpmag(j,i)+
     &            acof(i,lpmp,nat,2)*conjg(ccof(m,j,lo,nat,1))*
     &                    ulojug(lpmp,lo,m,n)
               perpmag(j,i)=perpmag(j,i)+
     &            bcof(i,lpmp,nat,2)*conjg(ccof(m,j,lo,nat,1))*
     &                    ulojdg(lpmp,lo,m,n)

              enddo !j
             enddo !i

            enddo !mp
           enddo !lp   

          enddo !m
         enddo !lo 
        endif !lo-apw

  
       enddo !nn 
      enddo !n

      end subroutine wann_perpmag

      end module m_wann_perpmag
