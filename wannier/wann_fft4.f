!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      module m_wann_fft4
      contains 
      subroutine wann_fft4(
     >          inputfilename,outputfilename,l_conjugate,
     >          rvecnum,rvec,kpoints,
     >          jspins_in,nkpts,l_bzsym,film,l_onedimens,
     >          l_soc,band_min,band_max,neigd,
     >          l_socmmn0,wan90version)

c*************************************************
c     Transform 4-dimensional matrices from 
c     Bloch representation to Wannier representation.
c     
c     Frank Freimuth, February 2011
c*************************************************

      use m_constants, only:pimach
      use m_wann_read_umatrix

      implicit none
      character,intent(in):: inputfilename*(*)
      character,intent(in):: outputfilename*(*)	
      logical,intent(in)  :: l_conjugate	
      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)
      integer, intent(in) :: jspins_in
      integer, intent(in) :: nkpts
      logical,intent (in) :: l_bzsym,l_soc
      logical,intent(in)  :: film
      integer,intent(in)  :: band_min(2),band_max(2),neigd
      logical, intent(in) :: l_socmmn0
      integer, intent(in) :: wan90version

      integer             :: ikpt,jspins
      integer             :: kpts
      logical             :: l_file
c      real                :: kpoints(3,nkpts)
      integer             :: num_wann,num_kpts,num_nnmax,jspin
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
      complex,allocatable :: u_matrix(:,:,:,:)
      complex,allocatable :: hwann(:,:,:,:)
      complex,allocatable :: hreal(:,:,:,:)
      complex,allocatable :: hsomtx(:,:,:,:)
      complex,allocatable :: hsomtx2(:,:,:,:)
      complex             :: fac,eulav,eulav1
      real                :: tmp_omi,rdotk,tpi,minenerg,maxenerg
      real, allocatable   :: minieni(:),maxieni(:)
      character           :: jobz,uplo
      integer             :: kpt,band,lee,lwork,lrwork,liwork,n,lda
      complex             :: value(4)
      logical             :: um_format
      logical             :: repro_eig
      logical             :: l_chk,l_proj,l_onedimens
      logical             :: have_disentangled
      integer,allocatable :: ndimwin(:,:)
      logical,allocatable :: lwindow(:,:,:)
      integer             :: chk_unit,nkp,ntmp,ierr
      character(len=33)   :: header
      character(len=20)   :: checkpoint
      real                :: tmp_latt(3,3), tmp_kpt_latt(3,nkpts)
      real                :: omega_invariant
      complex,allocatable :: u_matrix_opt(:,:,:,:)
      integer             :: num_bands
      logical             :: l_umdat
      real,allocatable    :: eigval2(:,:)
      real,allocatable    :: eigval_opt(:,:)
      real                :: scale,a,b
      character(len=2)    :: spinspin12(0:2)
      character(len=3)    :: spin12(2)
      character(len=6)    :: filename
      integer             :: jp,mp,kk,ii,jj,dir,rvecind
      integer             :: spin1,spin2

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      tpi=2*pimach()

      jspins=jspins_in
      if(l_soc)jspins=1

      write(6,*)"nkpts=",nkpts
c*****************************************************
c     get num_bands and num_wann from the proj file
c*****************************************************
      do j=1,0,-1
          inquire(file=trim('proj'//spinspin12(j)),exist=l_file)
          if(l_file)then
            filename='proj'//spinspin12(j)
            exit
          endif
      enddo
      if(l_file)then
          open (203,file=trim(filename),status='old')
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
      allocate( u_matrix_opt(num_bands,num_wann,nkpts,2) )
      allocate( u_matrix(num_wann,num_wann,nkpts,2) )
      allocate( lwindow(num_bands,nkpts,2) )
      allocate( ndimwin(nkpts,2) )

      do jspin=1,jspins  !spin loop
         call wann_read_umatrix2(
     >       nkpts,num_wann,num_bands,
     >       um_format,jspin,wan90version,
     <       have_disentangled,
     <       lwindow(:,:,jspin),
     <       ndimwin(:,jspin),
     <       u_matrix_opt(:,:,:,jspin),
     <       u_matrix(:,:,:,jspin))
         num_bands2=num_bands
      enddo !jspin   
      if(jspins.eq.1)then
         lwindow(:,:,2)        = lwindow(:,:,1)
         ndimwin(:,2)          = ndimwin(:,1)
         u_matrix_opt(:,:,:,2) = u_matrix_opt(:,:,:,1)
         u_matrix(:,:,:,2)     = u_matrix(:,:,:,1)
      endif

c****************************************************
c        Read the file "WF1.socspicom".
c**************************************************** 
      allocate( hsomtx(num_bands2,num_bands2,3,nkpts) )
      open(304,file=inputfilename,form='formatted')
      read(304,*)
      read(304,*)
      if(l_conjugate)then	
       do nkp=1,num_kpts
        do i=1,num_bands2
         do j=1,num_bands2
          do dir=1,3  
            read(304,*)dummy1,dummy2,dummy3,dummy4,a,b
            hsomtx(j,i,dir,nkp)=cmplx(a,-b)
          enddo !dir
         enddo !j
        enddo !i
       enddo !nkp
      else
       do nkp=1,num_kpts
        do i=1,num_bands2
         do j=1,num_bands2
          do dir=1,3
            read(304,*)dummy1,dummy2,dummy3,dummy4,a,b
		    hsomtx(j,i,dir,nkp)=cmplx(a,b)
          enddo !dir
         enddo !j
        enddo !i
       enddo !nkp
      endif		
      close(304)

c****************************************************************
c        Calculate matrix elements of SOC in the basis of
c        rotated Bloch functions.
c****************************************************************
      allocate( hsomtx2(num_wann,num_wann,3,nkpts) )
      write(6,*)"calculate matrix elements of SOC commutator
     &between wannier orbitals"

      if(have_disentangled) then       
       hsomtx2=0.0  
       do nkp=1,num_kpts
        print*,"nkp=",nkp  
        do dir=1,3  
           
        do j=1,num_wann
         do jp=1,num_wann  
            do m=1,ndimwin(nkp,1)
             do mp=1,ndimwin(nkp,1)  
              hsomtx2(jp,j,dir,nkp)=hsomtx2(jp,j,dir,nkp)+ 
     &            conjg(u_matrix_opt(mp,jp,nkp,1))*
     &                  hsomtx(mp,m,dir,nkp)*
     &                  u_matrix_opt(m,j,nkp,1)
             enddo !mp   
            enddo !m
         enddo !jp 
        enddo !j
        enddo !dir
       enddo !nkp
      else
       hsomtx2 = hsomtx
      end if !have_disentangled

      allocate(hwann(num_wann,num_wann,3,num_kpts))
      hwann=cmplx(0.0,0.0)
      wann_shift=0
      do k=1,num_kpts
         print*,"k=",k
       do dir=1,3  
       do m=1,num_wann
        do mp=1,num_wann
           do i=1,num_wann
            do j=1,num_wann
             hwann(mp,m,dir,k)=hwann(mp,m,dir,k)+
     *        conjg(u_matrix(j,mp,k,1))*
     *                hsomtx2(j,i,dir,k)*
     *              u_matrix(i,m,k,1)
            enddo !j
           enddo !i     
        enddo !mp
       enddo !m
       enddo
      enddo !k

c************************************************************
c        Calculate matrix elements in real space.
c***********************************************************      
      write(6,*)"calculate SOC-mat in rs"
      allocate(hreal(num_wann,num_wann,3,rvecnum))
      hreal=cmplx(0.0,0.0)
      do rvecind=1,rvecnum
       do k=1,nkpts
        rdotk=tpi*(  kpoints(1,k)*rvec(1,rvecind)+
     +                 kpoints(2,k)*rvec(2,rvecind)+
     +                 kpoints(3,k)*rvec(3,rvecind)  )
        fac=cmplx(cos(rdotk),-sin(rdotk))
        do dir=1,3
         do m2=1,num_wann
          do m1=1,num_wann
               hreal(m1,m2,dir,rvecind)=
     &         hreal(m1,m2,dir,rvecind)+
     &            fac*hwann(m1,m2,dir,k)
          enddo !m1
         enddo !m2
        enddo !dir
       enddo !k
      enddo !rvecind
      hreal=hreal/cmplx(real(nkpts),0.0)

      open(321,file=outputfilename,form='formatted')
      do rvecind=1,rvecnum
       r3=rvec(3,rvecind)
       r2=rvec(2,rvecind)
       r1=rvec(1,rvecind)

       do j=1,num_wann
        do i=1,num_wann
         do dir=1,3
            write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &            1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &          r1,r2,r3,i,j,dir,
     &          hreal(i,j,dir,rvecind) 
         enddo !dir
        enddo!i
       enddo !j

      enddo !rvecnum 
      close(321)

      deallocate(lwindow,u_matrix_opt,ndimwin)
      deallocate(u_matrix,hwann,hreal)

      end subroutine wann_fft4
      end module m_wann_fft4
