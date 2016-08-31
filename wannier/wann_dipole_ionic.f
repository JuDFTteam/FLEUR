!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_dipole_ionic
      contains
      subroutine wann_dipole_ionic(
     >               natd,pos,omtil,
     >               amat,taual,ntype,
     >               neq,zatom,l_absolute,
     >               l_invsubtract,
     <               ionic_moment)
      use m_wann_ioncharge_gen
      implicit none
      real,intent(in)              :: omtil
      integer,intent(in)           :: natd
      real,intent(in)              :: pos(3,natd)
      real,intent(in)              :: amat(3,3)
      real,intent(in)              :: taual(3,natd)
      integer,intent(in)           :: ntype
      integer,intent(in)           :: neq(ntype)
      real,intent(in)              :: zatom(ntype)
      logical,intent(in)           :: l_absolute
      logical,intent(in)           :: l_invsubtract
      real,intent(out)             :: ionic_moment(3)

      integer                      :: i
      integer                      :: num_atoms
      integer                      :: j,k,ind,jj,jspin
      real,allocatable             :: ioncharge(:)
      real                         :: polarization
      real                         :: ionic_polarization(3)
      real,parameter               :: elemchargmu=1.60217646e-13
      real,parameter               :: bohrtocm=0.529177e-8
      character(len=6)             :: filename
      logical                      :: l_file
      real                         :: pos_inv(3,natd)
      real                         :: taual_inv(3,natd)
      real                         :: coordinate
      real                         :: ionchargesum

      write(666,*)"*****************************"
      write(666,*)" Ionic terms                 "
      write(666,*)"*****************************"
      write(*,*)  "*****************************"
      write(*,*)  " Ionic terms                 "
      write(*,*)  "*****************************"

c-----count atoms
      num_atoms=0
      do j=1,ntype
         do k=1,neq(j)
            num_atoms=num_atoms+1
         enddo
      enddo

c-----get charges of ions
      allocate( ioncharge(num_atoms) )
      call wann_ioncharge_gen(
     >         num_atoms,ntype,natd,
     >         neq,zatom,pos,
     <         ioncharge)
      ionchargesum=sum(ioncharge(1:num_atoms))
      write(666,*)"ionchargesum=",ionchargesum

c-----Try to find the atomic positions of the inversion symmetric counterpart.
c-----The ferroelectric system is assumed to be connected to an inversion 
c-----symmetric system by an adiabatic path. The polarization of the inversion
c-----symmetric system is zero and consequently ambiguities in the evaluation of
c-----the polarization can be removed by subtracting the calculated polarizations
c-----of the two systems. We do not care here, whether the calculated inversion
c-----symmetric positions are physically reasonable: For example two atoms might
c-----end up at the same location. But this does not matter.
      pos_inv=0.0
      if(l_invsubtract)then
        do j=1,num_atoms
          do k=1,3
               coordinate=taual(k,j)
               if(abs(coordinate).lt.0.125)then
                  taual_inv(k,j)=0.0
               elseif(abs(coordinate).lt.0.375)then
                  taual_inv(k,j)=0.25*coordinate/abs(coordinate)
               elseif(abs(coordinate).lt.0.625)then
                  taual_inv(k,j)=0.5*coordinate/abs(coordinate)
               elseif(abs(coordinate).lt.0.875)then
                  taual_inv(k,j)=0.75*coordinate/abs(coordinate)
               else
                  taual_inv(k,j)=coordinate/abs(coordinate)
               endif   
          enddo
        enddo
        do j=1,num_atoms
           pos_inv(:,j)=matmul(amat,taual_inv(:,j))
        enddo
        do j=1,num_atoms
          print*,"old position:"
          print*,pos(:,j)
          print*,"new position:"
          print*,pos_inv(:,j)
          print*,"************************"
        enddo
      endif

c-----compute the ionic moment
      ionic_moment=0.0
      do j=1,num_atoms
         ionic_moment(:)=
     &        ionic_moment(:)+ioncharge(j)*(pos(:,j)-pos_inv(:,j))
      enddo

      ionic_polarization =
     =   ionic_moment/omtil*elemchargmu/((bohrtocm)**2)

      write(*,  fmt=555)ionic_moment(:)
      write(666,fmt=555)ionic_moment(:)
      if(.not.l_absolute)then
        write(*,  fmt=777)ionic_polarization(:)
        write(666,fmt=777)ionic_polarization(:)
      endif

 555  format("ionic moment = (",3f12.6,")a.u.")
 777  format("ionic polarization = (",3f12.6,")uC/cm**2")

      end subroutine wann_dipole_ionic
      end module m_wann_dipole_ionic










