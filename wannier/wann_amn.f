!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_amn
      use m_juDFT
      contains
      subroutine wann_amn(
     >               chi,nslibd,nwfsd,ntypd,nlod,llod,llo,nlo,
     >               lmaxd,jmtd,lmd,neq,natd,ikpt,nbnd,
     >               rmsh,rmt,jri,dx,lmax,
     >               us,dus,uds,duds,flo,
     >               ff,gg,acof,bcof,ccof,l_nocosoc,jspin,
     &               l_amn2_in,
     <               amn,
     >               bkpt) 
c*********************************************************************
c...here: twf - trial wf                      Y. Mokrousov June'06
c..                                                  and August'06

c...For the moment calculates only the overlaps inside the spheres.

c...As an input file the routine reads the file 'proj', where
c...it reads the following parameters:
c...nwfs: total number of wfs
c...For each atom (not atom type) the following values
c...are specified:
c..   lwf, mrwf - specify the angular part of the twf (see tutorial)
c..   rwf       - specifies the radial index of the twf  
c..   alpha,beta,gamma - the euler angles of the twf 
c..   zona      - the diffusivity of the twf
c..   regio     - 2*Rmt for instance, will calculate the radial 
c..               integrals of the twf and the wavefunction in the
c..               sphere around the atom center with the radius of
c..               2*Rmt. Outside the MT the integrals with the   
c..               radial part of the planewave (in terms of bessels)
c..               will be calculated without seeing other MTs, vacuum
c..               boundaries whatsoever. But for the moment works
c..               only for 1*Rmt (so regio=1.0)
c*********************************************************************
c     l_amn2_in=.true. => Variant of standard wann_amn:
c     Projections are defined in file pro2.1/pro2.2 
c     Includes the possibility to define displacements 
c     of the MT-localized
c     trial orbitals onto atoms in neighboring unit cells. Otherwise
c     equal to standard case. May
c     be useful for covalent systems: The proj 
c     file specifies the localized
c     trial orbital for one partner of the covalent bond 
c     while proj2 does 
c     the same for the second. In many cases the second partner will be
c     located in a neighboring unit cell, however.
c     Frank Freimuth, March 2007
c************************************************************************* 
c     Inclusion of Noco&&soc; tidied up version
c     Frank Freimuth
c*********************************************************************
      use m_intgr, only : intgr3
      use m_constants,only:pimach
      use m_dwigner
      use m_wann_tlmw
      use m_wann_rad_twf
      use m_eulerrot
      implicit none

      integer, intent (in)    :: nwfsd,nslibd,ntypd,nlod,llod
      integer, intent (in)    :: jmtd,lmaxd,lmd,natd,ikpt,nbnd
      integer, intent (in)    :: jri(ntypd),nlo(ntypd),llo(nlod,ntypd)
      integer, intent (in)    :: neq(ntypd),lmax(ntypd),jspin
      real,    intent (in)    :: rmsh(jmtd,ntypd),rmt(ntypd),dx(ntypd)
      real,    intent (in)    :: us(0:lmaxd,ntypd),dus(0:lmaxd,ntypd)
      real,    intent (in)    :: uds(0:lmaxd,ntypd),duds(0:lmaxd,ntypd)
      real,    intent (in)    :: ff(ntypd,jmtd,2,0:lmaxd)
      real,    intent (in)    :: gg(ntypd,jmtd,2,0:lmaxd)
      real,    intent (in)    :: flo(ntypd,jmtd,2,nlod)
      complex, intent (in)    :: acof(nslibd,0:lmd,natd)
      complex, intent (in)    :: bcof(nslibd,0:lmd,natd)
      complex, intent (in)    :: ccof(-llod:llod,nslibd,nlod,natd)
      logical, intent (in)    :: l_nocosoc
      logical, intent (inout) :: l_amn2_in
      complex, intent (inout) :: amn(nbnd,nwfsd)
      real,intent(in),optional:: bkpt(3)
      complex, intent (in)    :: chi
c...local
      integer          :: nwf,nwfs,nat,j,ntyp,ne,l,m,lm,iatom,i,mp,lo
      integer          :: ind(nwfsd),ntp(natd),banddummy
      integer          :: lwf(nwfsd),mrwf(nwfsd),rwf(nwfsd),spi(nwfsd)
      real             :: alpha(nwfsd),beta(nwfsd),gamma(nwfsd)
      real             :: jwf(nwfsd),jmwf(nwfsd),pi
      real,allocatable :: posshifts(:,:)
      real,allocatable :: weights(:)
      real             :: zona(nwfsd),regio(nwfsd),amx(3,3,nwfsd)
      real             :: rr,vl,vld,r0,vlo
      real             :: rads(nwfsd,0:3,jmtd,2),vlpr(jmtd),vlprd(jmtd)
      complex          :: tlmwf(0:3,-3:3,nwfsd),tlmwft(0:3,-3:3,nwfsd)
      logical          :: l_oldproj,l_amn2,l_file
      complex          :: d_wgn(-3:3,-3:3,1:3),wign(-3:3,-3:3,3,nwfsd)
      complex          :: ci,factor
      real             :: arg
      real             :: mrott(3,3),bmatt(3,3),imx(3,3)
      character(len=2) :: spin12(0:2)
      data spin12/'  ', '.1', '.2'/
      character(len=6) :: filename

      l_amn2=l_amn2_in
      if(.not.l_amn2_in)then
         inquire(file='pro2.1',exist=l_amn2_in)
      endif
      pi=pimach()
      ci=(0.0,1.0)
c..generates an array giving the atom type for each atom
      ntp(:) = 0
      iatom = 0
      do ntyp = 1,ntypd
         do nat = 1,neq(ntyp)
            iatom = iatom + 1
            ntp(iatom) = ntyp
         enddo
      enddo 

      if(l_amn2)then
       do j=jspin,1,-1
         inquire(file=trim('pro2'//spin12(j)),exist=l_file)
         if(l_file)then
            filename='pro2'//spin12(j)
            exit
         endif
       enddo
      else   
c..reading the proj.1 / proj.2 / proj file
       do j=jspin,0,-1
         inquire(file=trim('proj'//spin12(j)),exist=l_file)
         if(l_file)then
            filename='proj'//spin12(j)
            exit
         endif
       enddo
      endif

      if(l_file)then
        open (203,file=trim(filename),status='old')
        rewind (203)
      else
         CALL juDFT_error("no proj/proj.1/proj.2",calledby="wann_amn")
      endif  

      if(l_amn2)then
       allocate(posshifts(3,nwfsd))
       allocate(weights(nwfsd))
      endif    

      if(l_nocosoc)then
       read (203,*)nwfs,banddummy,l_oldproj
       if(.not.l_oldproj)then
        do nwf=1,nwfs
         read(203,*)ind(nwf),lwf(nwf),jwf(nwf),jmwf(nwf),rwf(nwf)
         read(203,*)alpha(nwf),beta(nwf),gamma(nwf),zona(nwf),regio(nwf)
         if(l_amn2) read(203,*) posshifts(1:3,nwf),weights(nwf)
        enddo !nwf
       else
        do nwf = 1,nwfs
         read (203,*) 
     &            ind(nwf),lwf(nwf),mrwf(nwf),rwf(nwf),spi(nwf)
         read (203,*) 
     &            alpha(nwf),beta(nwf),gamma(nwf),zona(nwf),regio(nwf)
         if(l_amn2) read(203,*) posshifts(1:3,nwf),weights(nwf)
        enddo !nwf
       endif !oldproj
      else
       read (203,*) nwfs
       do nwf = 1,nwfs
         read (203,*) 
     &            ind(nwf),lwf(nwf),mrwf(nwf),rwf(nwf)
         read (203,*) 
     &            alpha(nwf),beta(nwf),gamma(nwf),zona(nwf),regio(nwf)
         if(l_amn2) read(203,*) posshifts(1:3,nwf),weights(nwf)
       enddo !nwf
      endif !l_nocosoc
      rewind (203)
      close (203)

      if (ikpt.eq.1) then
      write (6,*) 'Number of trial WFs:',nwfs
      write (6,*)
      do nwf = 1,nwfs
        write (6,*) 'Twfs N:',nwf,' Atom N:',ind(nwf)
        write (6,*) 'l=',lwf(nwf),' mr=',mrwf(nwf),' r=',rwf(nwf)
        write (6,*) 'zona=',zona(nwf),' region=',regio(nwf),'*Rmt'
        write (6,*) 'alpha=',alpha(nwf),
     &          ' beta=',beta(nwf),' gamma=',gamma(nwf)
        write (6,*)
      enddo 
      endif

c..generating the radial twf function

      rads(:,:,:,:) = 0.

      call wann_rad_twf(
     >         nwfs,jmtd,natd,ind,rwf,zona,regio,
     >         us,dus,uds,duds,ff,gg,lmaxd,ikpt,
     >         ntypd,ntp,jri,rmsh,dx,
     >         nlod,flo,llo,nlo, 
     <         rads)

      open (100,file='rads')
      do i = 1,jmtd
c       write (100,'(i3,2x,4f10.6)') i,rads(1,0:3,i,1)
        write (100,'(f10.6,2x,4f10.6)') rmsh(i,1),rads(1,0:3,i,1)
      enddo
      close(100)

c..now generate the coefficients in the expansion in lattice 
c..harmonics of the angular part of the twf
 
      tlmwft(:,:,:) = cmplx(0.,0.)
      tlmwf(:,:,:)  = cmplx(0.,0.)
      
      if(l_nocosoc.and..not.l_oldproj)then
         call soc_tlmw(nwfs,lwf,jwf,jmwf,jspin,tlmwf)
      else
         call wann_tlmw(
     >          nwfs,lwf,mrwf,
     >          l_nocosoc,spi,jspin,
     <          tlmwf)  
      endif

      call eulerrot(nwfs,alpha,beta,gamma,amx)

      imx(:,:) = 0. ; imx(1,1) = 1. ; imx(2,2) = 1. ; imx(3,3) = 1.

c..performing the wigner rotation
c..These rotations are specified in the proj file in terms of the 
c..euler angles. The wave functions inside the muffin-tins are represen-
c..-ted in terms of the global-frame-lattice-harmonics, therefore, 
c..the wigner transformation has to be applied to tlmwf, resulting
c..in the tlmwft matrix, which is suitable for further use  

c..given the euler angles, the following procedure has to be 
c..performed in order to match the two local frames:
c..1. Rotation around the z-axis by alpha
c..2. Rotation around the x-axis by beta
c..3. Rotation around the z-axis by gamma again.
      call d_wigner(nwfs,amx,imx,3,wign)

c..now we transform the tlmwf coefficients

      do nwf = 1,nwfs
         tlmwft(0,:,nwf) = tlmwf(0,:,nwf)
         do l = 1,3
            do m = -l,l
               do mp = -l,l
                  tlmwft(l,m,nwf) = tlmwft(l,m,nwf) + 
     +                     wign(mp,m,l,nwf)*tlmwf(l,mp,nwf)
               enddo 
            enddo
         enddo         
      enddo

c..calculating the amn matrix

      vlpr(:) = 0. ; vlprd(:) = 0.

c...sum by wfs, each of them is localized at a certain mt
      do nwf = 1,nwfs
         if(l_amn2)then
           arg=-bkpt(1)*posshifts(1,nwf)
           arg=arg-bkpt(2)*posshifts(2,nwf)
           arg=arg-bkpt(3)*posshifts(3,nwf)
           arg=2*pi*arg
           factor=cmplx(cos(arg),sin(arg))*weights(nwf)
         else
           factor=1.0
         endif
         factor = factor*chi

         nat = ind(nwf)
         ntyp = ntp(nat)
c...sum by bands
         do ne = 1,nslibd
c...sum by l,m
            do l = 0,min(lmax(ntyp),3)
               do j = 1,jri(ntyp)
                  vlpr(j) = ff(ntyp,j,1,l)*rads(nwf,l,j,1)+
     +                      ff(ntyp,j,2,l)*rads(nwf,l,j,2) 
                  vlprd(j) = gg(ntyp,j,1,l)*rads(nwf,l,j,1) +
     +                       gg(ntyp,j,2,l)*rads(nwf,l,j,2)
                  if (rwf(nwf).gt.0)then
                     vlpr(j) = vlpr(j)*rmsh(j,ntyp)
                     vlprd(j) = vlprd(j)*rmsh(j,ntyp)
                  endif

               enddo

c..these integrations are not necessary if rads is the lin.comb. 
c..of the u_l and \dot{u}_l, but are necessary for other ways of
c..constructing the radial part, therefore, we do it anyway
 
               r0 = rmsh(1,ntyp)
               call intgr3(vlpr,rmsh(1,ntyp),dx(ntyp),jri(ntyp),vl)
               call intgr3(vlprd,rmsh(1,ntyp),dx(ntyp),jri(ntyp),vld)
               do m = -l,l
                  lm = l*(l+1) + m
                  amn(ne,nwf) = amn(ne,nwf) +
     +               tlmwft(l,m,nwf)*conjg(( acof(ne,lm,nat)*vl + 
     +                        bcof(ne,lm,nat)*vld )*(ci)**l)*factor
     
               enddo
            enddo

c..local orbitals
            if (nlo(ntyp).ge.1) then
               do lo = 1,nlo(ntyp)
                  l = llo(lo,ntyp)
                  do j = 1,jri(ntyp)
                     vlpr(j) = flo(ntyp,j,1,lo)*rads(nwf,l,j,1) +
     +                         flo(ntyp,j,2,lo)*rads(nwf,l,j,2)
                     if (rwf(nwf).gt.0)then
                     vlpr(j) = vlpr(j)*rmsh(j,ntyp)
                     endif
                  enddo


                  r0 = rmsh(1,ntyp)
                  call intgr3(vlpr,rmsh(1,ntyp),dx(ntyp),jri(ntyp),vlo)
                  do m = -l,l
                     amn(ne,nwf) = amn(ne,nwf) +
     +          tlmwft(l,m,nwf)*conjg((ccof(m,ne,lo,nat)*vlo)*(ci)**l)*
     *               factor
                  enddo
               enddo
            endif 

         enddo
      enddo

      return

      end subroutine wann_amn

c*********************************************************************
c     Calculate the expansion of the angular part of the trial orbital
c     in terms of spherical harmonics.
c     Version for Soc: The trial orbital is a spinor.
c     Frank Freimuth, June 2007
c*********************************************************************
      subroutine soc_tlmw(nwfs,lwf,jwf,jmwf,jspin,tlmwf)
      use m_clebsch
      implicit none

      integer, intent (in)  :: nwfs
      integer, intent (in)  :: lwf(nwfs)
      real,    intent (in)  :: jwf(nwfs),jmwf(nwfs)
      integer, intent (in)  :: jspin
      complex, intent (out) :: tlmwf(0:3,-3:3,nwfs)

      integer nwf,l,m
      complex tlm(0:3,-3:3,1:7)
      real j,jm


      do nwf=1,nwfs
         l=lwf(nwf)
         IF(l<0) CALL juDFT_error("not yet implemented",calledby
     +        ="wann_amn")
         j=jwf(nwf)
         jm=jmwf(nwf)
         IF(j<0) CALL juDFT_error("jwf",calledby ="wann_amn")
         IF( ABS(jm) - j  >1e-10) CALL juDFT_error("jmwf",calledby
     +        ="wann_amn")
         if( abs( l + 0.5 -j ).gt. 1e-10 .and.
     &       abs( l - 0.5 -j ).gt. 1e-10 )
     &     CALL juDFT_error ('regula trianguli violata',
     &                     calledby="wann_amn")
         tlmwf(0:3,-3:3,nwf)=cmplx(0.0,0.0)
         do m=-l,l
            tlmwf(l,m,nwf)=clebsch(real(l),0.5,real(m),1.5-jspin,j,jm)
         enddo
      enddo

      end subroutine soc_tlmw

      end module m_wann_amn
