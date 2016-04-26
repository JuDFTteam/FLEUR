      module m_wann_dipole_takehome
      contains
      subroutine wann_dipole_takehome(
     >               jspins,l_absolute,
     >               amat,bmat,omtil,
     <               electronic_moment)
c***********************************************************
c-----Check if the polarization may be reduced by adding 
c-----contributions due to electrons shifted by primitive
c-----lattice translations. This is the solution to the 
c-----problem of Wannier functions being pulled away from 
c-----their host atom to a mirror host atom during the
c-----Wannierization process.
c                            Frank Freimuth
c***********************************************************
      use m_constants, only: pimach
      implicit none

      integer, intent(in) :: jspins
      logical, intent(in) :: l_absolute
      real, intent(in)    :: amat(3,3),bmat(3,3),omtil
      real, intent(inout) :: electronic_moment(3,2)

      real,parameter      :: elemchargmu=1.60217646e-13
      real,parameter      :: bohrtocm=0.529177e-8
      real                :: electronic_polarization(3,2)
      integer             :: jspin,k
      real                :: int_electronic_moment(3,2)
      real                :: tpi

      tpi=2.0*pimach()

      write(666,*)"*********************************"
      write(666,*)" Try to minimize electronic term "
      write(666,*)"*********************************"
      write(*,*)  "*********************************"
      write(*,*)  " Try to minimize electronic term "
      write(*,*)  "*********************************"

      do jspin=1,jspins
        int_electronic_moment(:,jspin)=
     =      matmul( bmat,electronic_moment(:,jspin) )/tpi
        do k=1,3
           int_electronic_moment(k,jspin)=
     =       int_electronic_moment(k,jspin) -
     -       nint( int_electronic_moment(k,jspin) )  
        enddo
        electronic_moment(:,jspin) =
     =         matmul(amat,int_electronic_moment(:,jspin))
        write(*,  fmt=555)jspin,electronic_moment(:,jspin)
        write(666,fmt=555)jspin,electronic_moment(:,jspin)
        if(.not.l_absolute)then
           electronic_polarization=
     &            electronic_moment/omtil*elemchargmu/((bohrtocm)**2)
           write(*,  fmt=777)jspin,electronic_polarization(:,jspin)
           write(666,fmt=777)jspin,electronic_polarization(:,jspin)
        endif
      enddo  

 555  format("spin ",i1.1," electronic moment = (",3f12.6,") a.u.")
 777  format("spin ",i1.1,
     &          " electronic polarization = (",3f12.6,") uC/cm**2")

      end subroutine wann_dipole_takehome
      end module m_wann_dipole_takehome
      
