      module m_wann_dipole_electronic
      use m_juDFT
      contains
      subroutine wann_dipole_electronic(
     >               natd,pos,omtil,
     >               jspins,l_absolute,num_wann,
     <               electronic_moment)
c*********************************************
c     Calculate the dipole moment due to 
c     Wannier-electrons.
c
c     This subroutine can operate in two modes:
c
c     i)  l_absolute=.true.  => Simply evaluate
c         the center of mass moment of all 
c         Wannier functions.
c 
c     ii) l_absolute=.false. => Sum up the dipoles
c         that arise from the displacements of the
c         Wannier functions' centers from their 
c         host atoms. 
c
c     Frank Freimuth
c*********************************************
      use m_wann_readcenters
      implicit none
      integer, intent(in)  :: natd
      real,    intent(in)  :: pos(3,natd)
      real,    intent(in)  :: omtil
      integer, intent(in)  :: jspins
      logical, intent(in)  :: l_absolute
      integer, intent(in)  :: num_wann(2)
      real,    intent(out) :: electronic_moment(3,2)

      real                 :: electronic_polarization(3,2)
      real,parameter       :: elemchargmu=1.60217646e-13
      real,parameter       :: bohrtocm=0.529177e-8
      integer              :: jspin,j
      logical              :: l_file
      character(len=2)     :: spin12(0:2)
      character(len=6)     :: filename
      integer              :: num_wann_dum
      integer, allocatable :: waind(:)
      integer              :: nwf,k
      real, allocatable    :: wann_centers(:,:)

      data spin12/'  ', '.1', '.2'/

      write(666,*)"*****************************"
      write(666,*)" Electronic terms            "
      write(666,*)"*****************************"
      write(*,*)  "*****************************"
      write(*,*)  " Electronic terms            "
      write(*,*)  "*****************************"

      electronic_moment=0.0
      
      do jspin=1,jspins
        if(.not.l_absolute)then 
c--------reading the proj.1 / proj.2 / proj file
         do j=jspin,0,-1
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
     +           ="wann_dipole_electronic")
         endif  
         read(203,*)num_wann_dum
         if(num_wann_dum/=num_wann(jspin))  CALL juDFT_error
     +        ("num_wann_dum",calledby ="wann_dipole_electronic")
         allocate( waind(num_wann(jspin)) )
         do nwf=1,num_wann(jspin)
          read(203,*)waind(nwf)
          read(203,*)
         enddo
         close(203)
        endif 
        print*,"number of wannier functions= ",num_wann(jspin)
        write(6,*)"number of wannier functions",num_wann(jspin)
        allocate( wann_centers(3,num_wann(jspin)) )
        call wann_readcenters(
     >         num_wann(jspin),jspin,
     <         wann_centers(:,:))
        if(l_absolute)then
         do k=1,num_wann(jspin)  
           electronic_moment(:,jspin)=
     &     electronic_moment(:,jspin)-wann_centers(:,k)
         enddo  
        else
         do k=1,num_wann(jspin)
           electronic_moment(:,jspin)=
     &     electronic_moment(:,jspin)+
     +        pos(:,waind(k))-
     -        wann_centers(:,k)
         enddo
         deallocate(waind)
        endif 
        deallocate(wann_centers)
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
      end subroutine wann_dipole_electronic
      end module m_wann_dipole_electronic
