!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_read_umatrix
      use m_juDFT
c************************************************************
c     Read in the transformation matrix between Bloch
c     and Wannier functions which is calculated by 
c     wannier90.
c     Frank Freimuth, October 2006
c************************************************************
      CONTAINS
      SUBROUTINE wann_read_umatrix2(
     >               fullnkpts,num_wann,num_bands,
     >               um_format,jspin,wan90version,
     <               have_disentangled,
     <               lwindow,ndimwin,
     <               u_matrix_opt,u_matrix_tmp,m_matrix)
      implicit none
      integer,intent(in)  :: fullnkpts,jspin
      integer,intent(in)  :: num_wann,num_bands
      logical,intent(in)  :: um_format
      integer,intent(in)  :: wan90version
      logical,intent(out) :: have_disentangled
      logical,intent(out) :: lwindow(num_bands,fullnkpts)
      integer,intent(out) :: ndimwin(fullnkpts)
      complex,intent(out) :: u_matrix_opt(num_bands,num_wann,fullnkpts)
      complex,intent(out) :: u_matrix_tmp(num_wann,num_wann,fullnkpts)
      complex,intent(out),optional :: m_matrix(:,:,:,:)

      integer,allocatable :: exclude_bands(:)
      integer             :: num_kpts,chk_unit,mp_grid(3),l
      integer             :: i,j,k,nkp,ntmp,num_nnmax,nntot
      logical             :: l_chk,l_umdat
      character(len=3)    :: spin12(2)
      character(len=33)   :: header
      real                :: tmp_latt(3,3)
      real,allocatable    :: tmp_kpt_latt(:,:)
      character(len=20)   :: checkpoint
      real                :: tmp_omi
      real                :: omega_invariant
      character(len=2)    :: spinspin12(2)
      data spin12/'WF1' , 'WF2'/
      data spinspin12/'.1','.2'/

c**************************************************************
c     read in chk
c**************************************************************

      write(6,*)"read in chk"
      l_chk=.false.
      inquire (file=spin12(jspin)//'.chk',exist=l_chk)
      IF(.NOT.l_chk)  CALL juDFT_error("file chk not found",calledby
     +     ="wann_read_umatrix")

      chk_unit=152
      open(chk_unit,file=spin12(jspin)//'.chk',status='old',
     &     form='unformatted')


      if(wan90version.eq.1)then !wannier90 version 1.1
        ! Read comment line
        read(chk_unit) header
        write(6,*)header
        ! Real and reciprocal lattice (units: Angstroem)
        read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Real lattice
        write(6,*)tmp_latt
        read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Reciprocal lattice
        write(6,*)tmp_latt
 
        read(chk_unit) num_kpts  ! K-points
        write(6,*)"num_kpts=",num_kpts
        if (num_kpts.ne.fullnkpts) CALL 
     &     juDFT_error("num_kpts.ne.fullnkpts",calledby
     +     ="wann_read_umatrix")
      elseif((wan90version.eq.2).or.(wan90version.eq.3)) then 
        !wannier90 version 1.2 or wannier90 version 2.0

        ! Read comment line
        read(chk_unit) header
        write(6,*)header

        read(chk_unit) ntmp 

        read(chk_unit) ntmp 
        allocate(exclude_bands(ntmp))
        read(chk_unit) (exclude_bands(i),i=1,ntmp) ! Excluded bands

        ! Real and reciprocal lattice (units: Angstroem)
        read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Real lattice
        write(6,*)tmp_latt
        read(chk_unit) ((tmp_latt(i,j),i=1,3),j=1,3)  ! Reciprocal lattice
        write(6,*)tmp_latt
 
        read(chk_unit) num_kpts                ! K-points
        write(6,*)"num_kpts=",num_kpts
        if (num_kpts.ne.fullnkpts) stop "num_kpts.ne.fullnkpts"
        read(chk_unit) (mp_grid(i),i=1,3)         ! M-P grid
      else
          CALL juDFT_error("unknown wan90version",calledby
     +     ="wann_read_umatrix")
      endif

      allocate(tmp_kpt_latt(3,fullnkpts))  
      read(chk_unit) ((tmp_kpt_latt(i,nkp),i=1,3),nkp=1,fullnkpts)
c      do nkp=1,fullnkpts
c       do i=1,3
c          if (abs(tmp_kpt_latt(i,nkp)-kpoints(i,nkp)).gt.1.0e-6) 
!     &   CALL juDFT_error("mismatch of k-point lattice",calledby
!     ="wann_read_umatrix")
c       enddo
c      enddo
      deallocate(tmp_kpt_latt)

      read(chk_unit) ntmp                ! nntot
c      if (ntmp.ne.nntot)  CALL juDFT_error("mismatch in nntot",calledby="wann_read_umatrix")
      read(chk_unit) ntmp                ! num_wann
      IF (ntmp/=num_wann)  CALL juDFT_error("mismatch in num_wann"
     +     ,calledby ="wann_read_umatrix")

      read(chk_unit) checkpoint             ! checkpoint
      checkpoint=adjustl(trim(checkpoint))
      write(6,*)checkpoint





      read(chk_unit) have_disentangled      
! whether a disentanglement has been performed

      if (have_disentangled) then
         write(6,*)"You used disentangling"
         write(6,*)"Reading in disentangling information"
         read(chk_unit) omega_invariant     ! omega invariant
         write(6,*)omega_invariant


       ! U matrix opt
         read(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,num_kpts)
         read(chk_unit) (ndimwin(nkp),nkp=1,num_kpts)
         read(chk_unit)(((u_matrix_opt(i,j,nkp),i=1,num_bands)
     &  ,j=1,num_wann),nkp=1,num_kpts)

      endif



      if(.not.wan90version.eq.0)then

         read(chk_unit) (((u_matrix_tmp(i,j,k),i=1,num_wann),
     &                       j=1,num_wann),k=1,num_kpts)

      endif

      if(present(m_matrix))then
         write(*,*)"nntot=",nntot
         write(*,*)"num_kpts=",num_kpts
         write(*,*)"num_wann=",num_wann
               read(chk_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),
     &                   j=1,num_wann),k=1,nntot),l=1,num_kpts)
      endif   

      close(chk_unit)

      if(wan90version.eq.0)then
c****************************************************************
c                     read in _um.dat (old version of wannier90)
c****************************************************************
      l_umdat=.false.
      inquire (file=spin12(jspin)//'_um.dat',exist=l_umdat)
      IF(.NOT.l_umdat) CALL juDFT_error("where is your um_dat?",calledby
     +     ="wann_read_umatrix")


      open(111,file=spin12(jspin)//'_um.dat',form='unformatted')
      write(6,*)"read in um_data"
      read(111)header
      read(111)tmp_omi
      read(111) ntmp,num_kpts,num_nnmax
      IF(ntmp/=num_wann)CALL judft_error("mismatch in num_wann",calledby
     +     ="wann_read_umatrix")

        read(111)(((u_matrix_tmp(i,j,k),i=1,num_wann),j=1,num_wann),
     &    k=1,num_kpts)
        close(111)
      endif !wannier90 version

c$$$      if (um_format) then
c$$$          write(6,*)"write um_data to formatted file WF1.umn/WF2.umn"
c$$$          open(222,
c$$$     &       file=spin12(jspin)//'.umn',form='formatted')
c$$$c          write(222,*)header
c$$$c          write(222,*)tmp_omi
c$$$c          write(222,*)num_wann,num_kpts,num_nnmax
c$$$          write(222,*)"transformation between Bloch and Wannier"
c$$$          write(222,*)num
c$$$          do k=1,num_kpts
c$$$             do i=1,num_wann
c$$$                do j=1,num_wann
c$$$                  write(222,'(i3,3x,i3,3x,i3,3x,f20.16,f20.16)')
c$$$     &                                 i,j,k,u_matrix_tmp(i,j,k)
c$$$                enddo
c$$$             enddo
c$$$          enddo
c$$$          close(222)
c$$$      endif !um_format       
      END SUBROUTINE wann_read_umatrix2

      SUBROUTINE wann_read_umatrix(
     >               fullnkpts,num_wann,num_bands,
     >               um_format,jspin,wan90version,
     <               have_disentangled,
     <               lwindow,ndimwin,u_matrix)
      implicit none
      integer,intent(in)  :: fullnkpts,jspin
      integer,intent(in)  :: num_wann,num_bands
      logical,intent(in)  :: um_format
      integer,intent(in)  :: wan90version
      logical,intent(out) :: have_disentangled
      logical,intent(out) :: lwindow(num_bands,fullnkpts)
      integer,intent(out) :: ndimwin(fullnkpts)
      complex,intent(out) :: u_matrix(num_bands,num_wann,fullnkpts)

      complex             :: u_matrix_opt(num_bands,num_wann,fullnkpts)
      complex             :: u_matrix_tmp(num_wann,num_wann,fullnkpts)
      integer             :: i,j,k
      character(len=3)    :: spin12(2)
      data spin12/'WF1' , 'WF2'/

      call wann_read_umatrix2(
     >         fullnkpts,num_wann,num_bands,
     >         um_format,jspin,wan90version,
     <         have_disentangled,
     <         lwindow,ndimwin,u_matrix_opt,
     <         u_matrix_tmp)
      
      u_matrix(:,:,:)=0.0
      if (have_disentangled) then 
c**************************************************************
c     calculate u_matrix for the entangled case
c*************************************************************
         do k=1,fullnkpts
            do i=1,num_wann
               do j=1,num_wann
           u_matrix(:,j,k)= u_matrix(:,j,k)+
     &   u_matrix_opt(:,i,k)*u_matrix_tmp(i,j,k)
               enddo
            enddo
         enddo
      else   
         IF(num_bands/=num_wann) CALL juDFT_error("num_bands/=num_wann"
     +        ,calledby ="wann_read_umatrix")
        u_matrix(:,:,:)=u_matrix_tmp(:,:,:)
      endif   

      if (um_format) then
          write(6,*)"write um_data to formatted file WF1.umn/WF2.umn"
          open(222,
     &       file=spin12(jspin)//'.umn',form='formatted')
c          write(222,*)header
c          write(222,*)tmp_omi
c          write(222,*)num_wann,num_kpts,num_nnmax
          write(222,*)"transformation between Bloch and Wannier"
          write(222,*)num_bands,fullnkpts,num_wann
          do k=1,fullnkpts
             do i=1,num_wann
                do j=1,num_bands
                  write(222,'(i5,1x,i5,1x,i5,1x,f18.12,1x,f18.12)')
     &                                 j,i,k,u_matrix(j,i,k)
                enddo
             enddo
          enddo
          close(222)
      endif !um_format       


      END SUBROUTINE wann_read_umatrix
      END MODULE m_wann_read_umatrix
