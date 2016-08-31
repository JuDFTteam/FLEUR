!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_fft5
      contains 
      subroutine wann_fft5(
     >          l_conjugate,
     >          dim1,dim2,
     >          filename,title,
     >          rvecnum,rvec,kpoints,
     >          jspin,nkpts,film,
     >          l_soc,band_min,band_max,
     >          neigd,wan90version)

c*************************************************
c     Transform 5-dimensional matrices from 
c     Bloch representation to Wannier representation.
c     
c     Frank Freimuth, February 2011
c*************************************************

      use m_constants, only:pimach
      use m_wann_read_umatrix

      implicit none
      logical, intent(in) :: l_conjugate
      integer, intent(in) :: dim1,dim2

      character,intent(in) :: filename*(*)
      character,intent(in) :: title*(*)

      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)

      integer, intent(in) :: jspin
      integer, intent(in) :: nkpts
      logical, intent(in) :: film

      logical, intent(in) :: l_soc
      integer,intent(in)  :: band_min(2),band_max(2),neigd
      integer, intent(in) :: wan90version

      integer             :: ikpt
      integer             :: kpts
      logical             :: l_file
      integer             :: num_wann,num_kpts,num_nnmax
      integer             :: kspin,kkspin
      integer             :: wann_shift,num_wann2
      integer             :: i,j,k,m,info,r1,r2,r3,dummy1
      integer             :: dummy2,dummy3,dummy4,dummy5,dummy6
      integer             :: hopmin,hopmax,counter,m1,m2
      integer             :: num_bands2
      integer,allocatable :: iwork(:)
      real,allocatable    :: energy(:,:),ei(:)
      real,allocatable    :: eigw(:,:),rwork(:)
      complex,allocatable :: work(:),vec(:,:)
      complex,allocatable :: u_matrix(:,:,:)
      complex,allocatable :: hwann(:,:,:,:,:)
      complex,allocatable :: hreal(:,:,:,:,:)
      complex,allocatable :: hsomtx(:,:,:,:,:)
      complex,allocatable :: hsomtx2(:,:,:,:,:)
      complex             :: fac,eulav,eulav1
      real                :: tmp_omi,rdotk,tpi,minenerg,maxenerg
      real, allocatable   :: minieni(:),maxieni(:)
      character           :: jobz,uplo
      integer             :: kpt,band,lee,lwork,lrwork,liwork,n,lda
      complex             :: value(4)
      logical             :: um_format
      logical             :: repro_eig
      logical             :: l_chk,l_proj
      logical             :: have_disentangled
      integer,allocatable :: ndimwin(:)
      logical,allocatable :: lwindow(:,:)
      integer             :: chk_unit,nkp,ntmp,ierr
      character(len=33)   :: header
      character(len=20)   :: checkpoint
      real                :: tmp_latt(3,3), tmp_kpt_latt(3,nkpts)
      real                :: omega_invariant
      complex,allocatable :: u_matrix_opt(:,:,:)
      integer             :: num_bands
      logical             :: l_umdat
      real,allocatable    :: eigval2(:,:)
      real,allocatable    :: eigval_opt(:,:)
      real                :: scale,a,b
      character(len=2)    :: spinspin12(0:2)
      character(len=3)    :: spin12(2)
      character(len=6)    :: filename2
      integer             :: jp,mp,kk,ii,jj,dir,dir2,rvecind
      integer             :: spin1,spin2
      complex,parameter   :: ci=(0.0,1.0)

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      tpi=2*pimach()

      write(6,*)"nkpts=",nkpts
c*****************************************************
c     get num_bands and num_wann from the proj file
c*****************************************************
      do j=jspin,0,-1
          inquire(file=trim('proj'//spinspin12(j)),exist=l_file)
          if(l_file)then
            filename2='proj'//spinspin12(j)
            exit
          endif
      enddo
      if(l_file)then
          open (203,file=trim(filename2),status='old')
          rewind (203)
      else
          stop 'no proj/proj.1/proj.2' 
      endif
      read (203,*) num_wann,num_bands
      close (203)
      write(6,*)'According to proj there are ',num_bands,' bands'
      write(6,*)"and ",num_wann," wannier functions."

c****************************************************************
c        read in chk
c****************************************************************
      num_kpts=nkpts
      allocate( u_matrix_opt(num_bands,num_wann,nkpts) )
      allocate( u_matrix(num_wann,num_wann,nkpts) )
      allocate( lwindow(num_bands,nkpts) )
      allocate( ndimwin(nkpts) )

      call wann_read_umatrix2(
     >       nkpts,num_wann,num_bands,
     >       um_format,jspin,wan90version,
     <       have_disentangled,
     <       lwindow(:,:),
     <       ndimwin(:),
     <       u_matrix_opt(:,:,:),
     <       u_matrix(:,:,:))
      num_bands2=num_bands

c****************************************************
c        Read the file "filename".
c**************************************************** 
      allocate( hsomtx(dim1,dim2,num_bands,num_bands,nkpts) )
      open(304,file=filename,form='formatted')
      read(304,*)
      read(304,*)
      if(l_conjugate)then
       do nkp=1,num_kpts
        do i=1,num_bands
         do j=1,num_bands
          do dir2=1,dim2  
           do dir=1,dim1  
            read(304,*)dummy1,dummy2,dummy3,dummy4,dummy5,a,b
            hsomtx(dir,dir2,j,i,nkp)=cmplx(a,-b)
           enddo !dir
          enddo !dir2
         enddo !j
        enddo !i
       enddo !nkp
      else
       do nkp=1,num_kpts
        do i=1,num_bands
         do j=1,num_bands
          do dir2=1,dim2  
           do dir=1,dim1  
            read(304,*)dummy1,dummy2,dummy3,dummy4,dummy5,a,b
            hsomtx(dir,dir2,j,i,nkp)=cmplx(a,b)
           enddo !dir
          enddo !dir2
         enddo !j
        enddo !i
       enddo !nkp
      endif
      close(304)

c****************************************************************
c     fft5: transformation to Bloch-like    
c****************************************************************
      allocate( hsomtx2(dim1,dim2,num_wann,num_wann,nkpts) )
      write(6,*)"fft5: transformation to Bloch-like"

      if(have_disentangled) then       
       hsomtx2=0.0  
       do nkp=1,num_kpts
        do j=1,num_wann
         do jp=1,num_wann  
          do m=1,ndimwin(nkp)
           do mp=1,ndimwin(nkp) 
            do dir2=1,dim2
             do dir=1,dim1
              hsomtx2(dir,dir2,jp,j,nkp)=hsomtx2(dir,dir2,jp,j,nkp)+ 
     &            conjg(u_matrix_opt(mp,jp,nkp))*
     &                  hsomtx(dir,dir2,mp,m,nkp)*
     &                  u_matrix_opt(m,j,nkp)
             enddo !dir
            enddo !dir2
           enddo !mp   
          enddo !m
         enddo !jp 
        enddo !j
       enddo !nkp
      else
       hsomtx2 = hsomtx
      end if !have_disentangled

      allocate(hwann(dim1,dim2,num_wann,num_wann,num_kpts))
      hwann=cmplx(0.0,0.0)
      wann_shift=0
      do k=1,num_kpts
       print*,"k=",k
       do m=1,num_wann
        do mp=1,num_wann
         do i=1,num_wann
          do j=1,num_wann
           do dir2=1,dim2
            do dir=1,dim1
             hwann(dir,dir2,mp,m,k)=hwann(dir,dir2,mp,m,k)+
     *        conjg(u_matrix(j,mp,k))*
     *                hsomtx2(dir,dir2,j,i,k)*
     *              u_matrix(i,m,k)
            enddo
           enddo
          enddo !j
         enddo !i     
        enddo !mp
       enddo !m
      enddo !k

c************************************************************
c        Calculate matrix elements in real space.
c***********************************************************      
      write(6,*)"transform to rs"
      allocate(hreal(dim1,dim2,num_wann,num_wann,rvecnum))
      hreal=cmplx(0.0,0.0)
      do rvecind=1,rvecnum
       do k=1,nkpts
        rdotk=tpi*(  kpoints(1,k)*rvec(1,rvecind)+
     +                 kpoints(2,k)*rvec(2,rvecind)+
     +                 kpoints(3,k)*rvec(3,rvecind)  )
        fac=cmplx(cos(rdotk),-sin(rdotk))

        do m2=1,num_wann
         do m1=1,num_wann
          do dir2=1,dim2
           do dir=1,dim1

               hreal(dir,dir2,m1,m2,rvecind)=
     &         hreal(dir,dir2,m1,m2,rvecind)+
     &            fac*hwann(dir,dir2,m1,m2,k)
           enddo !dir
          enddo !dir2
         enddo !m1
        enddo !m2

       enddo !k
      enddo !rvecind
      hreal=hreal/cmplx(real(nkpts),0.0)

      open(321,file=title,form='formatted')
      do rvecind=1,rvecnum
       r3=rvec(3,rvecind)
       r2=rvec(2,rvecind)
       r1=rvec(1,rvecind)

        do j=1,num_wann
         do i=1,num_wann
          do dir2=1,dim2
           do dir=1,dim1
            write(321,'(i3,1x,i3,1x,i3,1x,i3,1x,i3,
     &            1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &          r1,r2,r3,i,j,dir,dir2,
     &          hreal(dir,dir2,i,j,rvecind) 
           enddo !dir
          enddo !dir2
         enddo!i
        enddo !j

      enddo !rvecnum 
      close(321)

c      deallocate(lwindow,u_matrix_opt,ndimwin)
c      deallocate(u_matrix,hwann,hreal)

      end subroutine wann_fft5
      end module m_wann_fft5
