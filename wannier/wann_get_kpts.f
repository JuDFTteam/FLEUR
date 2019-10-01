!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_get_kpts
      use m_juDFT
      USE m_types
      contains
      subroutine wann_get_kpts(input,kpts,
     >               l_bzsym,film,l_onedimens,l_readkpts,
     <               nkpts,kpoints)
c********************************************************
c     Read in the k-points from kpts/w90kpts file.
c
c     Frank Freimuth
c********************************************************
      implicit none
      TYPE(t_input), INTENT(IN) :: input
      TYPE(t_kpts), INTENT(IN)  :: kpts
      logical,intent(in)  :: l_bzsym,film
      logical,intent(in)  :: l_onedimens,l_readkpts
      integer,intent(out) :: nkpts
      real,intent(inout),allocatable  :: kpoints(:,:)

      real             :: scale
      integer          :: at,j
      integer          :: iter,len,num_wann,num_bands,nn,i
      logical          :: l_file

      if(l_bzsym)then
         inquire(file='w90kpts',exist=l_file)
         IF(.NOT.l_file) CALL juDFT_error("where is w90kpts?",calledby
     +        ="wann_get_kpts")
         open(987,file='w90kpts',status='old',form='formatted')
         read(987,*)nkpts, scale
         write(6,*)"wann_get_kpts: nkpts=",nkpts
         if(l_readkpts)then
            IF(SIZE(kpoints,1)/=3) CALL juDFT_error("wann_get_kpts: 1"
     +           ,calledby ="wann_get_kpts")
            IF(SIZE(kpoints,2)/=nkpts)CALL juDFT_error("wann_get_kpts:2"
     +           ,calledby ="wann_get_kpts")
            do iter=1,nkpts
               read(987,*)kpoints(:,iter)
            enddo
         endif
         close(987)
      else
             nkpts = kpts%nkpt
            write(6,*)"wann_get_kpts: nkpts=",nkpts
            if(l_readkpts)then
               do iter=1,nkpts
                  kpoints(:,iter) = kpts%bk(:,iter)
               enddo
            endif
      endif

      IF (l_readkpts) THEN
         do iter=1,nkpts
            write(6,*)kpoints(:,iter)
         enddo
      END IF

      end subroutine wann_get_kpts
      end module m_wann_get_kpts
