!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_dipole
      use m_juDFT
      contains
      subroutine wann_dipole(
     >               jspins,omtil,natd,pos,
     >               amat,
     >               ntype,neq,zatom)
c***************************************
c     Calculate electronic polarization.
c     Frank Freimuth
c***************************************
      use m_wann_readcenters
      implicit none
      integer,intent(in)           :: jspins
      real,intent(in)              :: omtil
      integer,intent(in)           :: natd
      real,intent(in)              :: pos(3,natd)
      real,intent(in)              :: amat(3,3)
      integer,intent(in)           :: ntype
      integer,intent(in)           :: neq(ntype)
      real,intent(in)              :: zatom(ntype)

      integer                      :: nwfs,i
      integer                      :: num_atoms,num_wann,num_wann2
      integer                      :: j,k,ind,jj,jspin
      real,allocatable             :: ioncharge(:)
      character(len=2),allocatable :: namat(:)
      character(len=2)             :: symbol
      real,allocatable             :: wann_centers1(:,:)
      real,allocatable             :: wann_centers2(:,:)
      integer,allocatable          :: wann_of_at(:)
      real                         :: polarization,charge
      integer                      :: num_symbols
      real                         :: ioni_polari(3)
      real                         :: shifted_polari(3)
      real                         :: size_polari
      real                         :: smallest_polari
      real                         :: final_polari(3)
      real                         :: electroni_polari(3,2)
      real                         :: elemchargmu,bohrtocm
      character*2                  :: namat2(0:103)
      character(len=2)             :: spin12(0:2)
      data spin12/'  ', '.1', '.2'/
      character(len=6)             :: filename
      logical                      :: l_file

      DATA namat2/'va',' H','He','Li','Be',
     +     ' B',' C',' N',' O',' F','Ne',
     +     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     +     ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se',
     +     'Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd',
     +     'Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce',
     +     'Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     +     'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb',
     +     'Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',' U','Np','Pu',
     +     'Am','Cm','Bk','Cf','Es','Fm','Md','No','Lw'/

      elemchargmu=1.60217646e-13
      bohrtocm=0.529177e-8

      open(666,file='polarization_out')

      num_atoms=0
      do j=1,ntype
         do k=1,neq(j)
            num_atoms=num_atoms+1
         enddo
      enddo

      allocate( namat(num_atoms) )
      ind=0
      do j=1,ntype
         do k=1,neq(j)
            ind=ind+1
            namat(ind)=namat2(nint(zatom(j)))
         enddo
      enddo

      do j=1,num_atoms
         print*,namat(j)," pos3=",pos(3,j)
      enddo

      allocate( ioncharge(num_atoms) )
      ioncharge=0.0
      open(400,file='IONS')
      read(400,*)num_symbols
      do j=1,num_symbols
         read(400,fmt=333)symbol,charge
         write(*,fmt=333)symbol,charge
         do k=1,num_atoms
            if(namat(k)==symbol)then
               ioncharge(k)=charge
            endif
         enddo
      enddo
      close(400)
 333  format(a2,1x,f10.6)

      open(300,file='ioncharge')
      do j=1,num_atoms
         write(300,*)ioncharge(j)
      enddo
      close(300)

      open(300,file='ioncharge')
      do j=1,num_atoms
         read(300,*)ioncharge(j)
      enddo
      close(300)

      ioni_polari=0.0
      do j=1,num_atoms
         ioni_polari(:)=
     &        ioni_polari(:)+ioncharge(j)*pos(:,j)
      enddo

      ioni_polari=ioni_polari/omtil*elemchargmu/((bohrtocm)**2)

      print*,"ioni_polari=",ioni_polari,"uC/cm2"
      write(666,*) "ioni_polari=",ioni_polari,"uC/cm2"

c*************************************************
c..reading the proj.1 / proj.2 / proj file
c*************************************************
      do j=1,0,-1
         inquire(file=trim('proj'//spin12(j)),exist=l_file)
         if(l_file)then
            filename='proj'//spin12(j)
            exit
         endif
      enddo
      if(l_file)then
        open (203,file=trim(filename),status='old')
        rewind (203)
      else
         CALL juDFT_error("no proj/proj.1/proj.2",calledby ="wann_dipole")
      endif  
      read(203,*)num_wann
      close(203)
      print*,"number of wannier functions= ",num_wann
      allocate(wann_centers1(3,num_wann))

      if(jspins.eq.2)then
        do j=2,0,-1
          inquire(file=trim('proj'//spin12(j)),exist=l_file)
          if(l_file)then
            filename='proj'//spin12(j)
            exit
          endif
        enddo
        if(l_file)then
          open (203,file=trim(filename),status='old')
          rewind (203)
        else
           CALL juDFT_error("no proj/proj.1/proj.2",calledby
     +          ="wann_dipole")
        endif  
        read(203,*)num_wann2
        close(203)
        print*,"number of wannier functions spin 2= ",num_wann2
        allocate(wann_centers2(3,num_wann2))
      endif

      electroni_polari=0.0
      call wann_readcenters(
     >         num_wann,1,
     <         wann_centers1(:,:))

      do k=1,num_wann
         electroni_polari(:,1)=
     &     electroni_polari(:,1)-wann_centers1(:,k)
      enddo

      if(jspins.eq.2)then
         call wann_readcenters(
     >         num_wann2,2,
     <         wann_centers2(:,:))
        do k=1,num_wann2
           electroni_polari(:,2)=
     &        electroni_polari(:,2)-wann_centers2(:,k)
        enddo
      endif   


      if(jspins.eq.1) electroni_polari = electroni_polari*2.0
      electroni_polari=
     &  electroni_polari/omtil*elemchargmu/((bohrtocm)**2)

      do j=1,jspins
         print*,"spin ",j," electronic_polarization=",
     &                electroni_polari(:,j)
         write(666,*)"spin ",j," electronic_polarization=",
     &                electroni_polari(:,j)
      enddo
      
      ioni_polari(:)=ioni_polari(:)+electroni_polari(:,1)
      if(jspins.eq.2)ioni_polari(:)=ioni_polari(:)
     &              +electroni_polari(:,2)

      print*,"total polarization=",ioni_polari(:)
      write(666,*)"total polarization=",ioni_polari(:)
!     Check if the polarization may be reduced by adding 
!     contributions due to electrons shifted by primitive
!     lattice translations.
      final_polari(:)=ioni_polari(:)
      size_polari=sqrt( (final_polari(1))**2 +
     +                  (final_polari(2))**2 +
     +                  (final_polari(3))**2 )
      smallest_polari=size_polari
      do i=-5,5
         do j=-5,5
            do k=-5,5
               shifted_polari(:)=ioni_polari(:)+
     +                  elemchargmu/((bohrtocm)**2)/omtil *
     *                  (amat(:,1)*k+amat(:,2)*j+amat(:,3)*i)
               size_polari=sqrt( (shifted_polari(1))**2 +
     +                           (shifted_polari(2))**2 +
     +                           (shifted_polari(3))**2 )
               if(size_polari.lt.smallest_polari)then
                  final_polari=shifted_polari
                  smallest_polari=size_polari
               endif
            enddo
         enddo
      enddo
      print*,"final polarization after shifting=",final_polari(:)
      write(666,*)"final polarization after shifting=",final_polari(:)
      close(666)



      end subroutine wann_dipole
      end module m_wann_dipole
