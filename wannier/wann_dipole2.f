c********************************************************
c     Calculate the electric polarization.
c     Frank Freimuth, October 2006
c********************************************************
      module m_wann_dipole2
      use m_juDFT
      contains
      subroutine wann_dipole2(
     >               jspins,pos,omtil,natd,
     >               l_nocosoc)
      use m_wann_readcenters
      implicit none
      integer,intent(in) :: jspins
      real,intent(in)    :: omtil
      integer,intent(in) :: natd
      real,intent(in)    :: pos(3,natd)
      logical,intent(in) :: l_nocosoc

      integer            :: jspin,jspins2
      logical            :: l_file
      integer,allocatable:: ind(:)
      integer            :: nwf,nwfs
      real,allocatable   :: distance(:,:)
      real               :: moment(3,jspins)
      real               :: moment2(3)
      real               :: elemchargmu,bohrtocm
      integer            :: shifted,ishift,ind1,ind2
      real,allocatable   :: centers(:,:)
      real,allocatable   :: spreads(:)

      open(204,file='dipole_out')
      write(204,*)"*****************************************"
      write(204,*)"Calculation of the electric polarization."
      write(204,*)"*****************************************"
      elemchargmu=1.60217646e-13
      bohrtocm=0.529177e-8

      jspins2=jspins
      if(l_nocosoc)jspins2=1
      do jspin=1,jspins2
         inquire(file='proj',exist=l_file)
         IF(.NOT.l_file) CALL juDFT_error("proj not found",calledby
     +        ="wann_dipole2")
         open(203,file='proj',form='formatted')
         read(203,*)nwfs
         print*,"nwfs=",nwfs
         allocate(ind(nwfs))
         allocate(distance(3,nwfs))
         allocate(centers(3,nwfs))
         allocate(spreads(nwfs))
         do nwf=1,nwfs
          read(203,*)ind(nwf)
          read(203,*)
         enddo
         close(203)
         write(204,*)"spin=",jspin
         call wann_readcenters(nwfs,jspin,centers,spreads)
         write(204,*)"*****centers:******"
         do nwf=1,nwfs
           write(204,fmt=1000)nwf,centers(:,nwf),spreads(nwf)
         enddo
         moment(:,jspin)=0.0
         do nwf=1,nwfs
            distance(:,nwf)=pos(:,ind(nwf))-centers(:,nwf)
            moment(:,jspin)=moment(:,jspin)+distance(:,nwf)
         enddo
         if((jspins.eq.1).and..not.l_nocosoc) 
     &                    moment(:,jspin)=moment(:,jspin)*2
         moment(:,jspin)=moment(:,jspin)/omtil
         moment(:,jspin)=moment(:,jspin)*elemchargmu/((bohrtocm)**2)
         write(204,*)"polarization due to shift of wannierorbitals"
         write(204,*)moment(:,jspin),"uC/cm2"
         deallocate(ind,distance,centers,spreads)
      enddo!jspin
1000  format(2x,'WF centre and spread', 
     &       i5,2x,'(',f10.6,',',f10.6,',',f10.6,' )',f15.8)
      inquire(file='dipole',exist=l_file)
      if(.not.l_file)then
            write(204,*)"no file dipole found"
      else
         moment2(:)=0.0
         open(224,file='dipole',form='formatted')
         read(224,*)shifted
         do ishift=1,shifted
               read(224,*)ind1,ind2 !electron taken from 1 and moved to 2
               moment2(:)=moment2(:)+pos(:,ind1)-pos(:,ind2)
         enddo
         close(224)
         moment2(:)=moment2(:)/omtil
         moment2(:)=moment2(:)*elemchargmu/(bohrtocm)**2
         write(204,*)"polarization due to electrons 
     &              being catched by other atoms"
         write(204,*)moment2(:),"uC/cm2"

      endif
      write(204,*)"total moment:"
      if(jspins2.eq.2)moment(:,1)=moment(:,1)+moment(:,2)
      if(l_file)then
         moment(:,1)=moment(:,1)+moment2(:)
      endif
      write(204,*)moment(:,1),"uC/cm2"
      close(204)
      end subroutine wann_dipole2
      end module m_wann_dipole2
