!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_dipole3
      contains
      subroutine wann_dipole3(
     >               jspins_in,
     >               cell,atoms,num_wann,
     >               l_nocosoc)
c***************************************
c     Calculate electronic polarization.
c     Frank Freimuth
c***************************************
      use m_wann_dipole_electronic
      use m_wann_dipole_takehome
      use m_wann_dipole_ionic
      USE m_types

      implicit none
      integer,intent(in)           :: jspins_in
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_atoms),INTENT(IN)     :: atoms
      integer,intent(in)           :: num_wann(2)
      logical,intent(in)           :: l_nocosoc

      integer                      :: jspins
      integer                      :: nwfs,i
      integer                      :: num_atoms
      integer                      :: j,k,ind,jj,jspin
      real,allocatable             :: ioncharge(:)
      character(len=2),allocatable :: namat(:)
      character(len=2)             :: symbol
      real,allocatable             :: wann_centers(:,:)
      integer,allocatable          :: wann_of_at(:)
      real                         :: polarization,charge
      integer                      :: num_symbols
      real                         :: ionic_moment(3)
      real                         :: shifted_polari(3)
      real                         :: size_polari
      real                         :: smallest_polari
      real                         :: final_polarization(3)
      real                         :: final_moment(3)
      real                         :: electronic_polari(3,2)
      real                         :: electronic_moment(3,2)
      real,parameter               :: elemchargmu=1.60217646e-13
      real,parameter               :: bohrtocm=0.529177e-8
      character*2                  :: namat2(0:103)
      character(len=2)             :: spin12(0:2)

      character(len=6)             :: filename
      logical                      :: l_file
      real                         :: pos_inv(3,atoms%nat)
      real                         :: taual_inv(3,atoms%nat)
      real                         :: coordinate
      integer                      :: nwf
      integer,allocatable          :: waind(:)
      integer                      :: yesorno

      data spin12/'  ', '.1', '.2'/
      call timestart("wann_dipole3")
      jspins=jspins_in
      yesorno=1
      if(l_nocosoc)then
         jspins=1
         yesorno=0
      endif
      open(666,file='polarization_out')

c-----calculate ionic contribution
      ionic_moment=0.0
      call wann_dipole_ionic(
     >         atoms,cell,.false.,
     >         .true.,
     <         ionic_moment)

c-----calculate electronic contribution
      electronic_moment=0.0
      call wann_dipole_electronic(
     >         atoms,cell,
     >         jspins,.false.,num_wann,
     <         electronic_moment)

c-----sum up terms
      write(666,*)"*****************************"
      write(666,*)" Sum of terms                "
      write(666,*)"*****************************"
      write(*,*)  "*****************************"
      write(*,*)  " Sum of  terms               "
      write(*,*)  "*****************************"
      final_moment(:)    = electronic_moment(:,1) +
     +            electronic_moment(:,jspins)*yesorno +
     +                     ionic_moment(:)
      write(*,  fmt=555)final_moment(:)
      write(666,fmt=555)final_moment(:)
      final_polarization = final_moment /
     /              cell%omtil*elemchargmu/((bohrtocm)**2)
      write(*,  fmt=777)final_polarization(:)
      write(666,fmt=777)final_polarization(:)

c-----Check if the polarization may be reduced by adding 
c-----contributions due to electrons shifted by primitive
c-----lattice translations. This is the solution to the 
c-----problem of Wannier functions being pulled away from 
c-----their host atom to a mirror host atom during the
c-----Wannierization process.
      call wann_dipole_takehome(
     >         jspins,.false.,
     >         cell,
     <         electronic_moment)

c-----sum up terms
      write(666,*)"*****************************"
      write(666,*)" Sum of terms                "
      write(666,*)"*****************************"
      write(*,*)  "*****************************"
      write(*,*)  " Sum of  terms               "
      write(*,*)  "*****************************"
      final_moment(:)    = electronic_moment(:,1) +
     +            electronic_moment(:,jspins)*yesorno +
     +                     ionic_moment(:)
      write(*,  fmt=555)final_moment(:)
      write(666,fmt=555)final_moment(:)
      final_polarization = final_moment /
     /              cell%omtil*elemchargmu/((bohrtocm)**2)
      write(*,  fmt=777)final_polarization(:)
      write(666,fmt=777)final_polarization(:)

      close(666)

 555  format("final moment = (",3f12.6,")a.u.")
 777  format("final polarization = (",3f12.6,")uC/cm**2")

      call timestop("wann_dipole3")
      end subroutine wann_dipole3
      end module m_wann_dipole3
