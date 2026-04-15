!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------


      module m_wann_perpmag
            use m_juDFT
      contains
      SUBROUTINE wann_perpmag(
     >               nbnd,atoms,mlh,nlhd,nlh,ntypsy,llh,
     >               nmem,bbmat,bmat,
     >               flo,f,g,clnu,
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

      use m_types
      use m_constants, only : pimach
      use m_sphbes
      use m_ylm
      use m_intgr, only : intgr3
      use m_gaunt, only: gaunt1

      IMPLICIT NONE

      TYPE(t_atoms), INTENT(IN) :: atoms

      integer, intent (in) :: nbnd
      INTEGER, INTENT (IN) :: mlh(:,0:,:)
      INTEGER, INTENT (IN) :: nlhd
      integer, intent (in) :: nlh(:)
      integer, intent (in) :: ntypsy(:)
      integer, intent (in) :: nmem(:,:)
      integer, intent (in) :: llh(:,:)
      real,    intent (in) :: bbmat(:,:)
      real,    intent (in) :: bmat(:,:)
      real,    intent (in) :: f(:,:,:,0:,:)
      real,    intent (in) :: g(:,:,:,0:,:)
      real,    intent (in) :: flo(:,:,:,:,:)
      complex, intent (in) :: clnu(:,0:,:)
      complex, intent (in) :: ujug(0:,0:,1:)
      complex, intent (in) :: ujdg(0:,0:,1:)
      complex, intent (in) :: djug(0:,0:,1:)
      complex, intent (in) :: djdg(0:,0:,1:)

      complex, intent (in) :: ujulog(0:,:,-atoms%llod:,:)
      complex, intent (in) :: djulog(0:,:,-atoms%llod:,:)
      complex, intent (in) :: ulojug(0:,:,-atoms%llod:,:)
      complex, intent (in) :: ulojdg(0:,:,-atoms%llod:,:)
      complex, intent (in) ::
     >    ulojulog(:,-atoms%llod:,:,-atoms%llod:,:)

      complex, intent (in) :: acof(:,0:,:,:)
      complex, intent (in) :: bcof(:,0:,:,:)
      complex, intent (in) ::
     >    ccof(-atoms%llod:,:,:,:,:)

      complex, intent (inout):: perpmag(:,:)

      integer :: nat,n,nn,l,m,lm,lp,mp,lpmp,i,j,lo,lop

      call timestart("wann_perpmag")
      nat=0
      do n=1,atoms%ntype
       do nn=1,atoms%neq(n)
        nat=nat+1
        do l=0,atoms%lmax(n)
         do m=-l,l  
          lm=l*(l+1)+m  
          do lp=0,atoms%lmax(n)
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

        if(atoms%nlo(n).ge.1)then
         do lo=1,atoms%nlo(n)
          l=atoms%llo(lo,n)
          do m=-l,l
           
           do lop=1,atoms%nlo(n)
            lp=atoms%llo(lop,n)
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

        if(atoms%nlo(n).ge.1)then
         do lo=1,atoms%nlo(n)
          l=atoms%llo(lo,n)
          do m=-l,l
           
           do lp=0,atoms%lmax(n)
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

      call timestop("wann_perpmag")
      end subroutine wann_perpmag

      end module m_wann_perpmag
