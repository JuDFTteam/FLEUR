!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_ujugaunt
      contains     
      SUBROUTINE wann_ujugaunt(
c***********************************************************************
c    wann_ujugaunt calculates integrals of radial wave functions with
c    bessel functions and multiplies them with an angular factor.
c    Calculating them only once gives some speed-up of wann_mmkb_sph.
c    Frank Freimuth, October 2006
c*********************************************************************** 
     >                          llod,nntot,kdiff,lmax,
     >                          ntype,ntypd,bbmat,bmat,
     >                          nlod,nlo,llo,flo,f,g,jri,rmsh,dx,jmtd,
     >                          lmaxd,lmd,
     <                          ujug,ujdg,djug,djdg,ujulog,djulog,
     <                          ulojug,ulojdg,ulojulog)
      use m_constants, only : pimach
      use m_matmul   , only : matmul3,matmul3r
      use m_sphbes
      use m_ylm
      use m_types, only : od_inp, od_sym
      use m_intgr, only : intgr3
      use m_gaunt, only: gaunt1
  
      IMPLICIT NONE
      integer, intent (in)  :: llod
      INTEGER, INTENT (IN)  :: nntot,ntype,ntypd
      INTEGER, INTENT (IN)  :: lmaxd,jmtd,lmd
      real,    intent (in)  :: kdiff(:,:) !(3,nntot)
      real,    intent (in)  :: bbmat(:,:) !(3,3)
      real,    intent (in)  :: bmat(:,:) !(3,3)
      integer, intent (in)  :: lmax(:) !(ntypd)
      integer, intent (in)  :: nlod
      integer, intent (in)  :: jri(:) !(ntypd)
      integer, intent (in)  :: nlo(:) !(ntypd)
      integer, intent (in)  :: llo(:,:) !(nlod,ntypd)    
      real,    intent (in)  :: f(:,:,:,0:) !(ntypd,jmtd,2,0:lmaxd)
      real,    intent (in)  :: g(:,:,:,0:) !(ntypd,jmtd,2,0:lmaxd)
      real,    intent (in)  :: flo(:,:,:,:) !(ntypd,jmtd,2,nlod)
      real,    intent (in)  :: rmsh(:,:) !(jmtd,ntypd)
      real,    intent (in)  :: dx(:) !(ntypd)

      complex, intent (out) :: ujug(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (out) :: ujdg(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (out) :: djug(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (out) :: djdg(0:,0:,:,:) !(0:lmd,0:lmd,1:ntype,1:nntot)
      complex, intent (out) :: ujulog(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (out) :: djulog(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (out) :: ulojug(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (out) :: ulojdg(0:,:,-llod:,:,:) !(0:lmd,nlod,-llod:llod,1:ntype,1:nntot)
      complex, intent (out) :: ulojulog(:,-llod:,:,-llod:,:,:) !(1:nlod,-llod:llod,1:nlod,-llod:llod,1:ntype,1:nntot)

      real, allocatable :: djd(:,:,:),ujd(:,:,:),uju(:,:,:)
      real, allocatable :: ujulo(:,:,:),djulo(:,:,:),ulojulo(:,:,:)
      integer           :: ikpt_b,i,lwn,n,lpp,lop,lo,l,lp
      integer           :: lmini,lmaxi,m,mp,llpp,mpp
      integer           :: lmpp,lminp,lmaxp,lm,lpmp
      real              :: rk,bpt(3),gs,jlpp(0:lmaxd)
      real              :: jj(0:lmaxd,jmtd),x(jmtd)
      real              :: bkrot(3)
      complex           :: ylmpp((lmaxd+1)**2),factor,ic

      ic = cmplx(0.,1.)

      allocate( djd(0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( ujd(0:lmaxd,0:lmaxd,0:lmaxd) )
      allocate( uju(0:lmaxd,0:lmaxd,0:lmaxd) ) 

      allocate( ujulo(nlod,0:lmaxd,0:lmaxd) )
      allocate( djulo(nlod,0:lmaxd,0:lmaxd) )

      allocate( ulojulo(nlod,nlod,0:lmaxd) )

      ujug = cmplx(0.0,0.0)
      ujdg = cmplx(0.0,0.0)
      djug = cmplx(0.0,0.0)
      djdg = cmplx(0.0,0.0)

      ujulog = cmplx(0.0,0.0)
      djulog = cmplx(0.0,0.0)
      ulojug = cmplx(0.0,0.0)
      ulojdg = cmplx(0.0,0.0)

      ulojulog = cmplx(0.0,0.0)

      do ikpt_b=1,nntot
        bpt(:)=kdiff(:,ikpt_b)
        rk = sqrt(dot_product(bpt,matmul(bbmat,bpt)))
        do n=1,ntype
         lwn = lmax(n)
c...generate the j_lpp(br) on the radial grid   
         do i = 1,jri(n)
           gs = rk*rmsh(i,n)
           call sphbes(lwn,gs,jlpp)
           jj(:,i) = jlpp(:)
         enddo
         do lpp = 0,lwn   ! lpp is the ang. momentum of the bessel function
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
               ulojulo(lo,lop,lpp) = 0. 
             else 
              do i = 1,jri(n)
                x(i) = ( flo(n,i,1,lo)*flo(n,i,1,lop)+
     +                   flo(n,i,2,lo)*flo(n,i,2,lop) )*jj(lpp,i)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),ulojulo(lo,lop,lpp))
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
             uju(l,lp,lpp) = 0.
             ujd(l,lp,lpp) = 0.
             djd(l,lp,lpp) = 0.
            else
             do i = 1,jri(n)
                x(i) = ( f(n,i,1,l)*f(n,i,1,lp)+
     +                   f(n,i,2,l)*f(n,i,2,lp) )*jj(lpp,i)
             enddo      
             call intgr3(x,rmsh(1:,n),dx(n),jri(n),uju(l,lp,lpp))

             do i = 1,jri(n)
                x(i) = ( f(n,i,1,l)*g(n,i,1,lp)+
     +                   f(n,i,2,l)*g(n,i,2,lp) )*jj(lpp,i)
             enddo      
             call intgr3(x,rmsh(1:,n),dx(n),jri(n),ujd(l,lp,lpp))

             do i = 1,jri(n)
                x(i) = ( g(n,i,1,l)*g(n,i,1,lp)+
     +                   g(n,i,2,l)*g(n,i,2,lp) )*jj(lpp,i)
             enddo     
             call intgr3(x,rmsh(1:,n),dx(n),jri(n),djd(l,lp,lpp))
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
               ujulo(lo,lp,lpp) = 0.
               djulo(lo,lp,lpp) = 0.
             else
              do i = 1,jri(n)
               x(i) = ( flo(n,i,1,lo)*f(n,i,1,lp)+
     +                  flo(n,i,2,lo)*f(n,i,2,lp) )*jj(lpp,i)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),ujulo(lo,lp,lpp))
              do i = 1,jri(n)
               x(i) = ( flo(n,i,1,lo)*g(n,i,1,lp)+
     +                  flo(n,i,2,lo)*g(n,i,2,lp) )*jj(lpp,i)
              enddo 
              call intgr3(x,rmsh(1:,n),dx(n),jri(n),djulo(lo,lp,lpp))
             endif
            enddo !lo  
           endif  ! local orbitals  
          enddo !lp
         enddo !lpp
c********************************************************************
c       multiply with gaunt-coefficient (apw-apw)
c********************************************************************
         bkrot=matmul(bpt,bmat)
         !call cotra3(bpt,bkrot,bmat)
         call ylm4(lwn,bkrot,ylmpp)
         do l = 0,lwn
          do m = -l,l
           lm=l*(l+1)+m  
           do lp = 0,lwn
            do mp = -lp,lp
             lpmp=lp*(lp+1)+mp  
             do lpp = 0,lwn
               llpp = lpp*(lpp+1)
               mpp = mp - m
               lmpp = llpp + mpp 
               lmini = abs(l-lpp)
               lmaxi = l+lpp
               if ((lmini.le.lp).and.(lp.le.lmaxi).and.
     &            (mod(l+lp+lpp,2).eq.0).and.(abs(mpp).LE.lpp))then  
                  factor=conjg(ylmpp(lmpp+1))*(ic**(l+lpp-lp))*
     *                     gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)
     
                ujug(lpmp,lm,n,ikpt_b)=ujug(lpmp,lm,n,ikpt_b)+
     +               factor*uju(l,lp,lpp)
                ujdg(lpmp,lm,n,ikpt_b)=ujdg(lpmp,lm,n,ikpt_b)+
     +               factor*ujd(l,lp,lpp)
                djug(lpmp,lm,n,ikpt_b)=djug(lpmp,lm,n,ikpt_b)+
     +               factor*ujd(lp,l,lpp)
                djdg(lpmp,lm,n,ikpt_b)=djdg(lpmp,lm,n,ikpt_b)+
     +               factor*djd(l,lp,lpp)
              
               endif
              enddo  ! lpp
            enddo ! mp
           enddo  ! lp
          enddo  ! m
         enddo   ! l
c******************************************************************
c       multiply with the gaunt-coefficient (apw-lo)
c******************************************************************
         if (nlo(n).ge.1) then 
         do lo = 1,nlo(n) 
          l = llo(lo,n)
          do m = -l,l
           lm=l*(l+1)+m  
           do lp = 0,lwn
            do mp = -lp,lp
              lpmp=lp*(lp+1)+mp 
              do lpp = 0,lwn
               llpp = lpp*(lpp+1)
               lmini = abs(l-lpp)
               lmaxi = l+lpp
               lminp = abs(lp-lpp)
               lmaxp = lp+lpp
               if ((lmini.le.lp).and.(lp.le.lmaxi).and.
     &            (mod(l+lp+lpp,2).eq.0).and.(abs(mp-m).LE.lpp)) then
                mpp = mp - m
                lmpp = llpp + mpp
                factor=conjg(ylmpp(lmpp+1))*(ic**(l+lpp-lp))*
     *                   gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)

                ujulog(lpmp,lo,m,n,ikpt_b)=ujulog(lpmp,lo,m,n,ikpt_b)+
     +               factor*ujulo(lo,lp,lpp)
                djulog(lpmp,lo,m,n,ikpt_b)=djulog(lpmp,lo,m,n,ikpt_b)+
     +               factor*djulo(lo,lp,lpp)
               endif

               if ((lminp.le.l).and.(l.le.lmaxp).and.
     &            (mod(l+lp+lpp,2).eq.0).and.(abs(m-mp).LE.lpp)) then
                mpp = m - mp
                lmpp = llpp + mpp
                factor=conjg(ylmpp(lmpp+1))*(ic**(lp+lpp-l))*
     *                   gaunt1(l,lpp,lp,m,mpp,mp,lmaxd)

                ulojug(lpmp,lo,m,n,ikpt_b)=ulojug(lpmp,lo,m,n,ikpt_b)+
     +               factor*ujulo(lo,lp,lpp)         
                ulojdg(lpmp,lo,m,n,ikpt_b)=ulojdg(lpmp,lo,m,n,ikpt_b)+
     +               factor*djulo(lo,lp,lpp)        
               endif
              enddo  ! lpp
            enddo ! mp
           enddo  ! lp
          enddo ! m lo
         enddo  ! lo
c*************************************************************
c         multiply with the gaunt-coefficient (lo-lo)
c*************************************************************
         do lo = 1,nlo(n)
          l = llo(lo,n)
          do m = -l,l
           lm=l*(l+1)+m  
           do lop = 1,nlo(n)
            lp = llo(lop,n)
            do mp = -lp,lp
              lpmp=lp*(lp+1)+mp 
              do lpp = 0,lwn
               llpp = lpp*(lpp+1)
               mpp = mp - m 
               lmpp = llpp + mpp
               lmini = abs(l-lpp)
               lmaxi = l+lpp
               if ((lmini.le.lp).and.(lp.le.lmaxi).and.
     &            (mod(l+lp+lpp,2).eq.0).and.(abs(mpp).LE.lpp))then  
                     ulojulog(lop,mp,lo,m,n,ikpt_b)=
     =                 ulojulog(lop,mp,lo,m,n,ikpt_b)+
     +                 conjg(ylmpp(lmpp+1))*(ic**(l+lpp-lp))*
     *                 gaunt1(lp,lpp,l,mp,mpp,m,lmaxd)*
     *                     ulojulo(lo,lop,lpp)
               endif
              enddo  ! lpp
            enddo ! mp lop
           enddo  ! lop
          enddo ! m lo
         enddo  ! lo           
         endif ! local orbitals on this atom
        enddo !ntype
      enddo !ikpt_b   
      deallocate(djd,ujd,uju,ujulo,djulo,ulojulo)

      end subroutine wann_ujugaunt
      end module m_wann_ujugaunt
