!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_ioncharge_gen
      use m_juDFT
      contains
      subroutine wann_ioncharge_gen(
     >               num_atoms,atoms,
     <               ioncharge)
c********************************************
c     Utility routine used to set up or read 
c     the file 'ioncharge'.
c     Frank Freimuth
c********************************************
      USE m_types
      implicit none
      integer, intent(in)          :: num_atoms
      TYPE(t_atoms), INTENT(IN)    :: atoms
      real,    intent(out)         :: ioncharge(num_atoms)

      character*2                  :: namat2(0:103)
      character(len=2),allocatable :: namat(:)
      integer                      :: ind,k,j
      logical                      :: l_file
      integer                      :: num_symbols
      real                         :: charge
      character(len=2)             :: symbol

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

      call timestart("wann_ioncharge_gen")

      allocate( namat(num_atoms) )
      ind=0
      do j=1,atoms%ntype
         do k=1,atoms%neq(j)
            ind=ind+1
            namat(ind)=namat2(nint(atoms%zatom(j)))
         enddo
      enddo

      inquire(file='ioncharge',exist=l_file)
      if(.not.l_file)then
       ioncharge=0.0
       open(400,file='IONS',status='old')
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
 333   format(a2,1x,f10.6)
       open(300,file='ioncharge')
       do j=1,num_atoms
          write(300,*)ioncharge(j)
       enddo
       close(300)
      endif 

      open(300,file='ioncharge')
      do j=1,num_atoms
         read(300,*)ioncharge(j)
      enddo
      close(300)

      write(666,*)"ionic charges:"
      do j=1,num_atoms
         write(666,fmt=111)namat(j),ioncharge(j),atoms%pos(:,j)
      enddo
 111  format(a2," ",f6.3," ",3f8.3)

      call timestop("wann_ioncharge_gen")
      end subroutine wann_ioncharge_gen
      end module m_wann_ioncharge_gen
