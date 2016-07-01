!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

c*************************************************
c     Calculate the matrix elements of the 
c     Pauli*Momentum matrix in real space from the
c     files WF1.chk (and WF1_um.dat) produced
c     by wannier90.
c     FF, January 2009
c*************************************************
      module m_wann_nabla_pauli_rs
      use m_juDFT
      contains 
      subroutine wann_nabla_pauli_rs(
     >               rvecnum,rvec,kpoints,
     >               jspins_in,nkpts,l_bzsym,film,l_onedimens,
     >               l_soc,band_min,band_max,neigd,
     >               l_socmmn0,wan90version)

      use m_constants, only:pimach
      use m_wann_read_umatrix
c$$$      use m_wann_wigner_seitz

      implicit none
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
      integer             :: dummy2,dummy3,dummy4
      integer             :: hopmin,hopmax,counter,m1,m2
      integer             :: num_bands2
      integer,allocatable :: iwork(:)
      real,allocatable    :: energy(:,:),ei(:)
      real,allocatable    :: eigw(:,:),rwork(:)
      complex,allocatable :: work(:),vec(:,:)
      complex,allocatable :: u_matrix(:,:,:),hwann(:,:,:,:,:)
      complex,allocatable :: hreal(:,:,:,:,:)
      complex,allocatable :: hrealsoc(:,:,:,:,:,:,:)
      complex,allocatable :: hwannsoc(:,:,:,:,:)
      complex,allocatable :: paulimat(:,:,:,:,:)
      complex,allocatable :: paulimat2(:,:,:,:,:)
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
      character(len=6)    :: filename
      integer             :: jp,mp,kk,kkk
      complex,parameter   :: ci=(0.0,1.0)
      integer             :: hopmin_z,hopmax_z
      integer             :: hopmin_y,hopmax_y
      integer             :: hopmin_x,hopmax_x
      integer             :: ii,spin,dir
      integer             :: rvecind,num(3),int_dummy
c      integer,allocatable :: rvec(:,:)

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      tpi=2*pimach()

      jspins=jspins_in
      if(l_soc)jspins=1

      write(6,*)"nkpts=",nkpts

c$$$c***************************************************
c$$$c     read in the kpoints from w90kpts or kpts
c$$$c***************************************************      
c$$$      if (l_bzsym) then
c$$$         l_file=.false.
c$$$         inquire(file='w90kpts',exist=l_file)
c$$$         if (.not.l_file) stop 'where is w90kpts?'
c$$$         open(177,file='w90kpts',form='formatted')
c$$$         read(177,*)kpts,scale
c$$$      else
c$$$         l_file=.false.
c$$$         inquire(file='kpts',exist=l_file)
c$$$         if(.not.l_file) stop 'where is kpts?'
c$$$         open(177,file='kpts',form='formatted')
c$$$         read(177,*)kpts,scale
c$$$      endif   
c$$$      if(kpts.ne.nkpts) stop 'mismatch in number of kpoints'
c$$$
c$$$      do ikpt = 1,nkpts 
c$$$         read(177,*)kpoints(:,ikpt)
c$$$      enddo 
c$$$
c$$$      if(film.and..not.l_onedimens)then
c$$$         kpoints(3,:)=0.0
c$$$      endif   
c$$$      kpoints=kpoints/scale
c$$$      close(177)

      do jspin=1,jspins  !spin loop
c*****************************************************
c     get num_bands and num_wann from the proj file
c*****************************************************
         do j=jspin,0,-1
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
            CALL juDFT_error("no proj/proj.1/proj.2",calledby
     +           ="wann_nabla_pauli_rs")
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
     >            nkpts,num_wann,num_bands,
     >            um_format,jspin,wan90version,
     <            have_disentangled,
     <            lwindow,ndimwin,
     <            u_matrix_opt,u_matrix)

         num_bands2=num_bands
         if(l_soc.and.l_socmmn0)then
          num_bands2=neigd
         endif
 
c****************************************************
c        Read the file "nabl.updown".
c****************************************************
         allocate( paulimat(3,2,num_bands2,num_bands2,nkpts) )
         open(304,file='nabl.updown',form='formatted')
         read(304,*)
         read(304,*)
         do nkp=1,num_kpts
          do i=1,num_bands2
           do j=1,num_bands2
            do k=1,3  
             read(304,*)dummy1,dummy2,dummy3,dummy4,a,b
             paulimat(k,1,i,j,nkp)=cmplx(a,b)
            enddo !k 
           enddo !j
          enddo !i
          do k=1,3
          paulimat(k,2,:,:,nkp)=ci*paulimat(k,1,:,:,nkp)
          paulimat(k,1,:,:,nkp)=paulimat(k,1,:,:,nkp)+
     &        transpose(conjg( paulimat(k,1,:,:,nkp) ))
          paulimat(k,2,:,:,nkp)=paulimat(k,2,:,:,nkp)+
     &        transpose(conjg( paulimat(k,2,:,:,nkp) ))
          enddo !k
         enddo !nkp
         close(304)


c****************************************************************
c        Calculate matrix elements of Pauli in the basis of
c        rotated Bloch functions.
c****************************************************************
         allocate( paulimat2(3,2,num_wann,num_wann,nkpts) )
         write(6,*)"calculate matrix elements of momentum operator
     &   between wannier orbitals"

         if(have_disentangled) then       
          paulimat2=0.0  
          do nkp=1,num_kpts
           do j=1,num_wann
            do jp=1,num_wann  
             do m=1,ndimwin(nkp)
              do mp=1,ndimwin(nkp)  
               do k=1,2
               do kk=1,3   
                paulimat2(kk,k,jp,j,nkp)=paulimat2(kk,k,jp,j,nkp)+ 
     &                 conjg(u_matrix_opt(mp,jp,nkp))*
     &                        paulimat(kk,k,mp,m,nkp)*
     &                       u_matrix_opt(m,j,nkp)
               enddo !kk 
               enddo !k 
              enddo !mp  
             enddo !m
            enddo !jp 
           enddo !j
          enddo !nkp
         else
          paulimat2(:,:,:,:,:)=paulimat(:,:,:,:,:)
         end if !have_disentangled

         allocate(hwann(3,2,num_wann,num_wann,num_kpts))
         hwann=cmplx(0.0,0.0)
         wann_shift=0
         if(l_socmmn0)then
            wann_shift=band_min(jspin)-1
         endif
         do k=1,num_kpts
          do m=1,num_wann
           do mp=1,num_wann
            do i=1,num_wann
             do j=1,num_wann
              do kk=1,2
               do kkk=1,3  
                hwann(kkk,kk,mp,m,k)=hwann(kkk,kk,mp,m,k)+
     *              conjg(u_matrix(j,mp,k))*
     *                paulimat2(kkk,kk,j,i,k)*
     *                   u_matrix(i,m,k)
               enddo !kkk 
              enddo !kk
             enddo !j
            enddo !i     
           enddo !mp
          enddo !m
         enddo !k

c************************************************************
c        Calculate matrix elements in real space.
c************************************************************      
         write(6,*)"calculate pauli-mat in rs"

c$$$       if(.false.)then !specify r-mesh by its boundaries
c$$$         hopmin_z=-5;hopmax_z=5
c$$$         hopmin_x=0;hopmax_x=0
c$$$         hopmin_y=0;hopmax_y=0
c$$$         rvecnum=(hopmax_z-hopmin_z+1)
c$$$         if(.not.l_onedimens.and.film)then
c$$$           hopmin_x=-5;hopmax_x=5
c$$$           hopmin_y=-5;hopmax_y=5
c$$$           hopmin_z=0;     hopmax_z=0
c$$$         else
c$$$           hopmin_x=-5;hopmax_x=5
c$$$           hopmin_y=-5;hopmax_y=5
c$$$         endif
c$$$         rvecnum=        (hopmax_z-hopmin_z+1)
c$$$         rvecnum=rvecnum*(hopmax_y-hopmin_y+1)
c$$$         rvecnum=rvecnum*(hopmax_x-hopmin_x+1)
c$$$
c$$$         allocate(rvec(3,rvecnum))
c$$$         rvecind=0
c$$$         do r3=hopmin_z,hopmax_z
c$$$          do r2=hopmin_y,hopmax_y
c$$$           do r1=hopmin_x,hopmax_x
c$$$            rvecind=rvecind+1
c$$$            if(rvecind.gt.rvecnum)stop 'wann_hopping:1'
c$$$            rvec(1,rvecind)=r1
c$$$            rvec(2,rvecind)=r2
c$$$            rvec(3,rvecind)=r3
c$$$           enddo !r1
c$$$          enddo !r2
c$$$         enddo !r3
c$$$       else !determine optimal r-mesh
c$$$         call wann_wigner_seitz(
c$$$     >       .true.,num,amat,0,
c$$$     <       rvecnum,rvec)
c$$$         allocate(rvec(3,rvecnum))
c$$$         call wann_wigner_seitz(
c$$$     >       .false.,num,amat,rvecnum,
c$$$     <       int_dummy,rvec)
c$$$
c$$$         open(333,file='wig_vectors',recl=1000)
c$$$         do ii=1,rvecnum
c$$$            write(333,*)ii,rvec(1,ii),rvec(2,ii),rvec(3,ii)
c$$$         enddo
c$$$         close(333)
c$$$
c$$$       endif

         allocate(hreal(3,2,num_wann,num_wann,rvecnum))
         hreal=cmplx(0.0,0.0)

         do rvecind=1,rvecnum
          do k=1,nkpts  
           rdotk=tpi*( kpoints(1,k)*rvec(1,rvecind)+
     &                 kpoints(2,k)*rvec(2,rvecind)+
     &                 kpoints(3,k)*rvec(3,rvecind) )
           fac=cmplx(cos(rdotk),-sin(rdotk))
           do m2=1,num_wann
            do m1=1,num_wann
             do spin=1,2  
              do dir=1,3  
               hreal(dir,spin,m1,m2,rvecind)=
     &                   hreal(dir,spin,m1,m2,rvecind)+
     &                   fac*hwann(dir,spin,m1,m2,k)
              enddo !dir 
             enddo !spin
            enddo !m1  
           enddo !m2  
          enddo !k
         enddo !rvecind
         hreal=hreal/cmplx(real(nkpts),0.0)

         open(321,file='rspaulisvx'//spinspin12(jspin),form='formatted')
          do rvecind=1,rvecnum
           r3=rvec(3,rvecind)
           r2=rvec(2,rvecind)
           r1=rvec(1,rvecind)
           do j=1,num_wann
            do i=1,num_wann
             do kk=1,3   
              write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &           1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &          r1,r2,r3,i,j,kk,hreal(kk,1,i,j,rvecind) 
             enddo !kk 
            enddo !i
           enddo !j   
          enddo !rvecnum          
         close(321)

         open(321,file='rspaulisvy'//spinspin12(jspin),form='formatted')
          do rvecind=1,rvecnum
           r3=rvec(3,rvecind)
           r2=rvec(2,rvecind)
           r1=rvec(1,rvecind)
           do j=1,num_wann
            do i=1,num_wann
             do kk=1,3   
              write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &           1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &          r1,r2,r3,i,j,kk,hreal(kk,2,i,j,rvecind) 
             enddo !kk 
            enddo !i
           enddo !j   
          enddo !rvecnum 
         close(321)


         deallocate(lwindow,u_matrix_opt,ndimwin)
         deallocate(u_matrix,hwann,hreal)
      enddo !jspin

      end subroutine wann_nabla_pauli_rs
      end module m_wann_nabla_pauli_rs
