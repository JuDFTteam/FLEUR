!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_gabeffgagaunt
      contains
      SUBROUTINE wann_gabeffgagaunt(
     >               vTot,
     >               l_perpmagatlres,perpmagl,
     >               neq,mlh,nlhd,nlh,ntypsy,llh,lmax,
     >               nmem,ntype,ntypd,bbmat,bmat,
     >               nlod,llod,nlo,llo,flo,f,g,jri,rmsh,dx,jmtd,
     >               lmaxd,lmd,clnu,
     <               ujug,ujdg,djug,djdg,ujulog,djulog,
     <               ulojug,ulojdg,ulojulog)
c*************************************************************************
c    wann_gabeffgagaunt calculates integrals of radial wave functions with
c    the exchange field and multiplies them with an angular factor.
c
c    Frank Freimuth, 2010
c************************************************************************* 
      USE m_types
      use m_constants
      use m_matmul   , only : matmul3,matmul3r
      use m_sphbes
      use m_ylm
      use m_intgr, only : intgr3
      use m_gaunt, only: gaunt1

      IMPLICIT NONE
      TYPE(t_potden),    INTENT(IN) :: vTot
      logical, intent (in)  :: l_perpmagatlres
      integer, intent (in)  :: perpmagl
      integer, intent (in)  :: neq(:) 
      INTEGER, INTENT (IN)  :: mlh(:,0:,:)!(memd,0:nlhd,ntypsd)
      INTEGER, INTENT (IN)  :: nlhd
      integer, intent (in)  :: nlh(:)
      integer, intent (in)  :: ntypsy(:)
      integer, intent (in)  :: nmem(0:,:)
      integer, intent (in)  :: ntype,ntypd
      integer, intent (in)  :: llh(0:,:)
      INTEGER, INTENT (IN)  :: lmaxd,jmtd,lmd
      real,    intent (in)  :: bbmat(:,:)!(3,3)
      real,    intent (in)  :: bmat(:,:)!(3,3)
      integer, intent (in)  :: lmax(:)!(ntypd)
      integer, intent (in)  :: nlod,llod
      integer, intent (in)  :: jri(:)!(ntypd)
      integer, intent (in)  :: nlo(:)!(ntypd)
      integer, intent (in)  :: llo(:,:)!(nlod,ntypd)    
      real,    intent (in)  :: f(:,:,:,0:,:)!(ntypd,jmtd,2,0:lmaxd,2)
      real,    intent (in)  :: g(:,:,:,0:,:)!(ntypd,jmtd,2,0:lmaxd,2)
      real,    intent (in)  :: flo(:,:,:,:,:)!(ntypd,jmtd,2,nlod,2)
      real,    intent (in)  :: rmsh(:,:)!(jmtd,ntypd)
      real,    intent (in)  :: dx(:)!(ntypd)
      complex, intent (in)  :: clnu(:,0:,:)!(memd,0:nlhd,ntypsd)

      complex, intent (out) :: ujug(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)
      complex, intent (out) :: ujdg(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)
      complex, intent (out) :: djug(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)
      complex, intent (out) :: djdg(0:,0:,1:)!(0:lmd,0:lmd,1:ntype)

      complex, intent (out) :: ujulog(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (out) :: djulog(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (out) :: ulojug(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (out) :: ulojdg(0:,:,-llod:,:) !(0:lmd,nlod,-llod:llod,1:ntype)
      complex, intent (out) :: ulojulog(:,-llod:,:,-llod:,:) !(1:nlod,-llod:llod,1:nlod,-llod:llod,1:ntype)

      real, allocatable :: djd(:,:,:),ujd(:,:,:),uju(:,:,:)
      real, allocatable :: dju(:,:,:),loju(:,:,:),lojd(:,:,:)
      real, allocatable :: ujulo(:,:,:),djulo(:,:,:),ulojulo(:,:,:)
      integer :: i,lwn,n,lpp,lop,lo,l,lp,lmini,lmaxi,m,mp,llpp,mpp
      integer :: lmpp,lminp,lmaxp,lm,lpmp,na,ns,lh,j,jmem
      real    :: rk,gs,jj(0:lmaxd,jmtd),x(jmtd)
      complex :: factor,ic
      logical :: l_doit
      real,allocatable :: beff(:,:,:)

      ic = cmplx(0.,1.)

      allocate( beff(jmtd,0:nlhd,ntype) )

      if(.false.)then !from file
        open(777,file='beff',form='unformatted')
        read(777)beff
        close(777)  
      else !from vTot
        beff=(vTot%mt(:,:,:,1)-vTot%mt(:,:,:,2))/2.0
! In vgen/vgen_finalize.F90 there is the line
! vTot%mt(:atoms%jri(n),0,n,js)  = atoms%rmsh(:atoms%jri(n),n)*vTot%mt(:atoms%jri(n),0,n,js)/sfp_const
! The following line removes this factor:
        do n=1,ntype
           beff(:jri(n),0,n)=beff(:jri(n),0,n)*sfp_const/rmsh(:jri(n),n)
        enddo   
      endif   

      if(l_perpmagatlres)then
       na=1  
       do n=1,ntype
        ns=ntypsy(na)  
        do lh = 0,nlh(ns)
         lpp=llh(lh,ns)
         if(lpp.ne.perpmagl)then
           beff(:,lh,n)=0.0 
         endif
        enddo !lh 
        na=na+neq(n)
       enddo !n
      endif !l_perpmagatlres  

      allocate ( djd(0:lmaxd,0:lmaxd,0:nlhd),
     +           dju(0:lmaxd,0:lmaxd,0:nlhd),
     +           ujd(0:lmaxd,0:lmaxd,0:nlhd),
     +           uju(0:lmaxd,0:lmaxd,0:nlhd), 

     +           ujulo(nlod,0:lmaxd,0:nlhd),
     +           djulo(nlod,0:lmaxd,0:nlhd),
     +           loju(nlod,0:lmaxd,0:nlhd),
     +           lojd(nlod,0:lmaxd,0:nlhd),
     +           ulojulo(nlod,nlod,0:nlhd) )


      na=1
      do n=1,ntype
       ns=ntypsy(na)  
       lwn = lmax(n)
       do lh = 0,nlh(ns)
        lpp=llh(lh,ns)  
c***************************************************************************
c...the local orbitals overlaps
c***************************************************************************
        if (nlo(n).GE.1) then
         do lop = 1,nlo(n)
          do lo = 1,nlo(n)
             l = llo(lo,n)
             lp = llo(lop,n)
             lmini = abs(lp - l)
             lmaxi = lp + l
c..the gaunt conditions
             if ((mod(l+lp+lpp,2).eq.1) .or. (lpp.LT.lmini) .or.
     +             (lpp.gt.lmaxi)) then
               ulojulo(lo,lop,lh) = 0. 
             else 
              do i = 1,jri(n)
                x(i) = ( flo(n,i,1,lo,1)*flo(n,i,1,lop,2)+
     +                   flo(n,i,2,lo,1)*flo(n,i,2,lop,2) )*beff(i,lh,n)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),ulojulo(lo,lop,lh))
             endif
          enddo
         enddo
        endif ! local orbitals 
c**************************************************************************
c...overlaps of the apws only
c**************************************************************************
        do lp = 0,lwn
         do l = 0,lwn
          lmini = abs(lp-l)
          lmaxi = lp + l
c..gaunt conditions
          if ((mod(l+lp+lpp,2).eq.1) .or. (lpp.LT.lmini) .or.
     +             (lpp.gt.lmaxi)) then
           uju(l,lp,lh) = 0.
           dju(l,lp,lh) = 0.
           ujd(l,lp,lh) = 0.
           djd(l,lp,lh) = 0.
          else
           do i = 1,jri(n)
                x(i) = ( f(n,i,1,l,1)*f(n,i,1,lp,2)+
     +                   f(n,i,2,l,1)*f(n,i,2,lp,2) )*beff(i,lh,n)
           enddo      
           call intgr3(x,rmsh(1:,n),dx(n),jri(n),uju(l,lp,lh))

           do i = 1,jri(n)
                x(i) = ( g(n,i,1,l,1)*f(n,i,1,lp,2)+
     +                   g(n,i,2,l,1)*f(n,i,2,lp,2) )*beff(i,lh,n)
           enddo      
           call intgr3(x,rmsh(1:,n),dx(n),jri(n),dju(l,lp,lh))

           do i = 1,jri(n)
                x(i) = ( f(n,i,1,l,1)*g(n,i,1,lp,2)+
     +                   f(n,i,2,l,1)*g(n,i,2,lp,2) )*beff(i,lh,n)
           enddo      
           call intgr3(x,rmsh(1:,n),dx(n),jri(n),ujd(l,lp,lh))

           do i = 1,jri(n)
                x(i) = ( g(n,i,1,l,1)*g(n,i,1,lp,2)+
     +                   g(n,i,2,l,1)*g(n,i,2,lp,2) )*beff(i,lh,n)
           enddo     
           call intgr3(x,rmsh(1:,n),dx(n),jri(n),djd(l,lp,lh))
          endif
         enddo ! l

c********************************************************************
c...overlaps of the lo's with the apws 
c********************************************************************
         if (nlo(n).GE.1) then
          do lo = 1,nlo(n)
             l = llo(lo,n)
             lmini = abs(lp-l)
             lmaxi = lp + l
c..gaunt conditions
             if ((mod(l+lp+lpp,2).eq.1) .OR. (lpp.lt.lmini) .or.
     +             (lpp.gt.lmaxi)) then
               ujulo(lo,lp,lh) = 0.
               djulo(lo,lp,lh) = 0.
                loju(lo,lp,lh) = 0.
                lojd(lo,lp,lh) = 0.
             else
              do i = 1,jri(n)
               x(i) = ( flo(n,i,1,lo,1)*f(n,i,1,lp,2)+
     +                  flo(n,i,2,lo,1)*f(n,i,2,lp,2) )*beff(i,lh,n)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),loju(lo,lp,lh))

              do i = 1,jri(n)
               x(i) = ( flo(n,i,1,lo,1)*g(n,i,1,lp,2)+
     +                  flo(n,i,2,lo,1)*g(n,i,2,lp,2) )*beff(i,lh,n)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),lojd(lo,lp,lh))

              do i = 1,jri(n)
               x(i) = ( flo(n,i,1,lo,2)*f(n,i,1,lp,1)+
     +                  flo(n,i,2,lo,2)*f(n,i,2,lp,1) )*beff(i,lh,n)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),ujulo(lo,lp,lh))

              do i = 1,jri(n)
               x(i) = ( flo(n,i,1,lo,2)*g(n,i,1,lp,1)+
     +                  flo(n,i,2,lo,2)*g(n,i,2,lp,1) )*beff(i,lh,n)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),djulo(lo,lp,lh))

             endif
          enddo !lo  
         endif  ! local orbitals  
        enddo !lp
       enddo !lh
c********************************************************************
c       multiply with gaunt-coefficient (apw-apw)
c********************************************************************
       ujug(:,:,n)=cmplx(0.0,0.0) 
       ujdg(:,:,n)=cmplx(0.0,0.0)
       djug(:,:,n)=cmplx(0.0,0.0) 
       djdg(:,:,n)=cmplx(0.0,0.0)
       do l = 0,lwn
          do m = -l,l
           lm=l*(l+1)+m  
           do lp = 0,lwn
            do mp = -lp,lp
             lpmp=lp*(lp+1)+mp  
             do lh = 0,nlh(ns)
               lpp=llh(lh,ns) 
               llpp = lpp*(lpp+1)
               mpp = mp - m
               lmpp = llpp + mpp 
               lmini = abs(l-lpp)
               lmaxi = l+lpp
               jmem=0
               do j=1,nmem(lh,ns)
                  if(mlh(j,lh,ns)==mpp)then
                     jmem=j
                     exit
                  endif
               enddo
               l_doit=           (lmini.le.lp)
               l_doit=l_doit.and.(lp.le.lmaxi)
               l_doit=l_doit.and.(mod(l+lp+lpp,2).eq.0)
               l_doit=l_doit.and.(abs(mpp).LE.lpp)
               l_doit=l_doit.and.(jmem.ne.0)
               if ( l_doit )then  
                  factor=ic**(l-lp)*
     *               gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)*
     *               clnu(jmem,lh,ns)     
                  ujug(lpmp,lm,n)=ujug(lpmp,lm,n)+
     +               factor*uju(lp,l,lh)
                  ujdg(lpmp,lm,n)=ujdg(lpmp,lm,n)+
     +               factor*ujd(lp,l,lh)
                  djug(lpmp,lm,n)=djug(lpmp,lm,n)+
     +               factor*dju(lp,l,lh)
                  djdg(lpmp,lm,n)=djdg(lpmp,lm,n)+
     +               factor*djd(lp,l,lh)
              
               endif
             enddo  ! lh
            enddo ! mp
           enddo  ! lp
          enddo  ! m
       enddo   ! l
c******************************************************************
c       multiply with the gaunt-coefficient (apw-lo)
c******************************************************************
       if(nlo(n).ge.1) then
        ulojug(:,:,:,n)=cmplx(0.0,0.0) 
        ulojdg(:,:,:,n)=cmplx(0.0,0.0)
        ujulog(:,:,:,n)=cmplx(0.0,0.0)
        djulog(:,:,:,n)=cmplx(0.0,0.0)

        do lo = 1,nlo(n)
          l = llo(lo,n)
          do m = -l,l
           lm=l*(l+1)+m
  
           do lp=0,lwn
            do mp = -lp,lp
             lpmp=lp*(lp+1)+mp
  
             do lh = 0,nlh(ns)
               lpp=llh(lh,ns) 
               llpp = lpp*(lpp+1)
               mpp = mp - m  !<mp|m+mpp>
               lmpp = llpp + mpp 
               lmini = abs(l-lpp)
               lmaxi = l+lpp
               jmem=0
               do j=1,nmem(lh,ns)
                  if(mlh(j,lh,ns)==mpp)then
                     jmem=j
                     exit
                  endif
               enddo
               l_doit=           (lmini.le.lp)
               l_doit=l_doit.and.(lp.le.lmaxi)
               l_doit=l_doit.and.(mod(l+lp+lpp,2).eq.0)
               l_doit=l_doit.and.(abs(mpp).LE.lpp)
               l_doit=l_doit.and.(jmem.ne.0)
               if ( l_doit )then  
                  factor=ic**(l-lp)*
     *               gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)*
     *               clnu(jmem,lh,ns)     

                  ujulog(lpmp,lo,m,n)=
     =               ujulog(lpmp,lo,m,n)+
     *                  ujulo(lo,lp,lpp)*factor
              
                  djulog(lpmp,lo,m,n)=
     =               djulog(lpmp,lo,m,n)+
     *                  djulo(lo,lp,lpp)*factor
               endif

               mpp = m - mp  !<m|mp+mpp>
               lmpp = llpp + mpp 
               lmini = abs(lp-lpp)
               lmaxi = lp+lpp
               jmem=0
               do j=1,nmem(lh,ns)
                  if(mlh(j,lh,ns)==mpp)then
                     jmem=j
                     exit
                  endif
               enddo
               l_doit=           (lmini.le.l)
               l_doit=l_doit.and.(l.le.lmaxi)
               l_doit=l_doit.and.(mod(l+lp+lpp,2).eq.0)
               l_doit=l_doit.and.(abs(mpp).LE.lpp)
               l_doit=l_doit.and.(jmem.ne.0)
               if ( l_doit )then  
                  factor=ic**(lp-l)*
     *               gaunt1(l,lpp,lp,m,mpp,mp,lmaxd)*
     *               clnu(jmem,lh,ns)     

                  ulojug(lpmp,lo,m,n)=
     =               ulojug(lpmp,lo,m,n)+
     *                  loju(lo,lp,lpp)*factor
              
                  ulojdg(lpmp,lo,m,n)=
     =               ulojdg(lpmp,lo,m,n)+
     *                  lojd(lo,lp,lpp)*factor
               endif

             enddo  ! lh
            enddo ! mp
           enddo  ! lop
          enddo  ! m
        enddo   ! lo
       endif ! local orbitals on this atom


c*************************************************************
c         multiply with the gaunt-coefficient (lo-lo)
c*************************************************************
       if(nlo(n).ge.1) then

        ulojulog(:,:,:,:,n)=cmplx(0.0,0.0) 
        do lo = 1,nlo(n)
          l = llo(lo,n)
          do m = -l,l
           lm=l*(l+1)+m
  
           do lop = 1,nlo(n)
            lp = llo(lop,n)
            do mp = -lp,lp
             lpmp=lp*(lp+1)+mp
  
             do lh = 0,nlh(ns)
               lpp=llh(lh,ns) 
               llpp = lpp*(lpp+1)
               mpp = mp - m !<mp|m+mpp>
               lmpp = llpp + mpp 
               lmini = abs(l-lpp)
               lmaxi = l+lpp
               jmem=0
               do j=1,nmem(lh,ns)
                  if(mlh(j,lh,ns)==mpp)then
                     jmem=j
                     exit
                  endif
               enddo
               l_doit=           (lmini.le.lp)
               l_doit=l_doit.and.(lp.le.lmaxi)
               l_doit=l_doit.and.(mod(l+lp+lpp,2).eq.0)
               l_doit=l_doit.and.(abs(mpp).LE.lpp)
               l_doit=l_doit.and.(jmem.ne.0)
               if ( l_doit )then  
                  factor=ic**(l-lp)*
     *               gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)*
     *               clnu(jmem,lh,ns)     

                  ulojulog(lop,mp,lo,m,n)=
     =               ulojulog(lop,mp,lo,m,n)+
     *                  ulojulo(lop,lo,lpp)*factor
              
               endif
             enddo  ! lh
            enddo ! mp
           enddo  ! lop
          enddo  ! m
        enddo   ! lo
       endif ! local orbitals on this atom


       na=na+neq(n)  
      enddo !ntype
      deallocate(djd,ujd,uju,dju)
      
      deallocate(ujulo,djulo,ulojulo,loju,lojd)

      end subroutine wann_gabeffgagaunt

      end module m_wann_gabeffgagaunt
