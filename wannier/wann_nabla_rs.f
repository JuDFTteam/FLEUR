      module m_wann_nabla_rs
      use m_juDFT
      contains 
      subroutine wann_nabla_rs(
     >          rvecnum,rvec,kpoints,
     >          jspins_in,nkpts,l_bzsym,film,l_onedimens,
     >          l_soc,band_min,band_max,neigd,
     >          l_socmmn0,wan90version)
c*************************************************
c     Calculate the matrix elements of nabla
c     in real space from the
c     files WF1.chk (and WF1_um.dat) produced
c     by wannier90.
c     FF, December 2009
c*************************************************
      use m_constants, only:pimach
      use m_wann_read_umatrix

      implicit none
      integer, intent(in) :: rvecnum
      integer, intent(in) :: rvec(:,:)
      real,    intent(in) :: kpoints(:,:)
      integer, intent(in) :: jspins_in
      integer, intent(in) :: nkpts
      logical, intent(in) :: l_bzsym
      logical, intent(in) :: film
      logical, intent(in) :: l_onedimens

      logical, intent(in) :: l_soc
      integer, intent(in) :: band_min(2),band_max(2),neigd

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
      complex,allocatable :: u_matrix(:,:,:),hwann(:,:,:,:)
      complex,allocatable :: hreal(:,:,:,:)
      complex,allocatable :: hrealsoc(:,:,:,:,:,:,:)
      complex,allocatable :: hwannsoc(:,:,:,:,:)
      complex,allocatable :: nablamat(:,:,:,:,:)
      complex,allocatable :: nablamat2(:,:,:,:)
      complex             :: fac,eulav,eulav1
      real                :: tmp_omi,rdotk,tpi,minenerg,maxenerg
      real, allocatable   :: minieni(:),maxieni(:)
      character           :: jobz,uplo
      integer             :: kpt,band,lee,lwork,lrwork,liwork,n,lda
      complex             :: value(4)
      logical             :: um_format,l_she
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
      character(len=6)    :: filename
      integer             :: jp,mp,kk
      integer             :: rvecind

      data spinspin12/'  ','.1' , '.2'/
      data spin12/'WF1','WF2'/

      inquire(file='she',exist=l_she)

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
     +           ="wann_nabla_rs")
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
 
c************************************************
c        Read the files "WF1.nabl" and "WF2.nabl"
c************************************************
         allocate( nablamat(3,num_bands2,num_bands2,nkpts,2) )
         open(304,file=spin12(jspin)//'.nabl',form='formatted')
         read(304,*)
         read(304,*)
         do nkp=1,num_kpts
          do i=1,num_bands2
           do j=1,num_bands2
            do k=1,3  
             read(304,*)dummy1,dummy2,dummy3,dummy4,a,b
             nablamat(k,j,i,nkp,1)=cmplx(a,-b)
            enddo !k 
           enddo !j
          enddo !i
         enddo !nkp
         close(304)

         if(l_soc)then
          open(304,file=spin12(2)//'.nabl',form='formatted')
          read(304,*)
          read(304,*)
          do nkp=1,num_kpts
           do i=1,num_bands2
            do j=1,num_bands2
             do k=1,3  
              read(304,*)dummy1,dummy2,dummy3,dummy4,a,b
              nablamat(k,j,i,nkp,2)=cmplx(a,-b)
             enddo !k 
            enddo !j
           enddo !i
          enddo !nkp
          close(304)
         endif

c**************************************
c        Enforce Hermiticity
c**************************************
         do nkp=1,num_kpts
          do i=1,num_bands2
           do j=1,i
            do k=1,3
             nablamat(k,i,j,nkp,1)=( nablamat(k,i,j,nkp,1)+
     &                         conjg(nablamat(k,j,i,nkp,1)) )/2.0
             nablamat(k,j,i,nkp,1)=conjg(nablamat(k,i,j,nkp,1))
             if(l_soc)then
               nablamat(k,i,j,nkp,2)=( nablamat(k,i,j,nkp,2)+
     &                         conjg(nablamat(k,j,i,nkp,2)) )/2.0
               nablamat(k,j,i,nkp,2)=conjg(nablamat(k,i,j,nkp,2))
             endif
            enddo !k
           enddo !j
          enddo !i
         enddo !nkp

c**************************************
c        Nabla-sum
c**************************************
         if( l_soc )then
          if(l_she)then
           do nkp=1,num_kpts
            do i=1,num_bands2
             do j=1,num_bands2
              do k=1,3  
               nablamat(k,j,i,nkp,1)=nablamat(k,j,i,nkp,1)
     &                            -nablamat(k,j,i,nkp,2)
              enddo !k 
             enddo !j
            enddo !i
           enddo !nkp
          else   
           do nkp=1,num_kpts
            do i=1,num_bands2
             do j=1,num_bands2
              do k=1,3  
               nablamat(k,j,i,nkp,1)=nablamat(k,j,i,nkp,1)
     &                            +nablamat(k,j,i,nkp,2)
              enddo !k 
             enddo !j
            enddo !i
           enddo !nkp
          endif
         endif !jspins.eq.2

c****************************************************************
c        Calculate matrix elements of Nabla in the basis of
c        rotated Bloch functions.
c****************************************************************
         allocate( nablamat2(3,num_wann,num_wann,nkpts) )
         write(6,*)"calculate matrix elements of momentum operator
     &   between wannier orbitals"

         if(have_disentangled) then       
          nablamat2=0.0  
          do nkp=1,num_kpts
           do j=1,num_wann
            do jp=1,num_wann  
             do m=1,ndimwin(nkp)
              do mp=1,ndimwin(nkp)  
               do k=1,3  
                nablamat2(k,jp,j,nkp)=nablamat2(k,jp,j,nkp)+ 
     &                 conjg(u_matrix_opt(mp,jp,nkp))*
     &                        nablamat(k,mp,m,nkp,1)*
     &                       u_matrix_opt(m,j,nkp)
               enddo !k 
              enddo !mp  
             enddo !m
            enddo !jp 
           enddo !j
          enddo !nkp
         else
          nablamat2(:,:,:,:)=nablamat(:,:,:,:,1)
         end if                    !have_disentangled

         allocate(hwann(3,num_wann,num_wann,num_kpts))
         hwann=cmplx(0.0,0.0)
         wann_shift=0
         if(l_socmmn0)then
            wann_shift=band_min(jspin)-1
         endif
         do k=1,num_kpts
          do m=1,num_wann
           do i=1,num_wann
            do j=1,num_wann
             do mp=1,num_wann
              do kk=1,3  
                hwann(kk,mp,m,k)=hwann(kk,mp,m,k)+
     *              conjg(u_matrix(j,mp,k))*
     *                nablamat2(kk,j,i,k)*
     *                    u_matrix(i,m,k)
              enddo !kk
             enddo !mp   
            enddo !j
           enddo !i
          enddo !m
         enddo !k

c************************************************************
c        Calculate matrix elements in real space.
c***********************************************************      
         write(6,*)"calculate nabla-mat in rs"
c$$$         hopmin=-5
c$$$         hopmax=5
c$$$         allocate(hreal(3,num_wann,num_wann,hopmin:hopmax,
c$$$     &         hopmin:hopmax,hopmin:hopmax))
c$$$         hreal=cmplx(0.0,0.0)
c$$$         do r3=hopmin,hopmax
c$$$          do r2=hopmin,hopmax
c$$$           do r1=hopmin,hopmax
c$$$            do k=1,nkpts  
c$$$              rdotk=tpi*(kpoints(1,k)*r1+kpoints(2,k)*r2+
c$$$     &                                   kpoints(3,k)*r3)
c$$$              fac=cmplx(cos(rdotk),-sin(rdotk))/nkpts
c$$$              do i=1,num_wann
c$$$               do j=1,num_wann
c$$$                do kk=1,3
c$$$                 hreal(kk,j,i,r1,r2,r3)=
c$$$     &                   hreal(kk,j,i,r1,r2,r3)+
c$$$     &                   fac*hwann(kk,j,i,k)
c$$$                enddo !kk
c$$$               enddo !j
c$$$              enddo !i
c$$$            enddo !k
c$$$           enddo !r1
c$$$          enddo !r2
c$$$         enddo !r3


         allocate(hreal(3,num_wann,num_wann,rvecnum))
         hreal=cmplx(0.0,0.0)

         do rvecind=1,rvecnum
            do k=1,nkpts  
              rdotk=tpi*(kpoints(1,k)*rvec(1,rvecind)+
     &                   kpoints(2,k)*rvec(2,rvecind)+
     &                   kpoints(3,k)*rvec(3,rvecind) )
              fac=cmplx(cos(rdotk),-sin(rdotk))
              do i=1,num_wann
               do j=1,num_wann
                do kk=1,3
                 hreal(kk,j,i,rvecind)=
     &                   hreal(kk,j,i,rvecind)+
     &                   fac*hwann(kk,j,i,k)
                enddo !kk
               enddo !j
              enddo !i
            enddo !k
         enddo !rvecind
         hreal=hreal/cmplx(real(nkpts),0.0)

         if(l_she)then
          open(321,file='rs_sv'//spinspin12(jspin),form='formatted')
         else   
          open(321,file='rsnabla'//spinspin12(jspin),form='formatted')
         endif

         do rvecind=1,rvecnum
            r3=rvec(3,rvecind)
            r2=rvec(2,rvecind)
            r1=rvec(1,rvecind)             
            do j=1,num_wann
             do i=1,num_wann
              do kk=1,3  
               write(321,'(i3,1x,i3,1x,i3,1x,i3,
     &                 1x,i3,1x,i3,1x,f20.8,1x,f20.8)')
     &                 r1,r2,r3,i,j,kk,hreal(kk,i,j,rvecind) 
              enddo !kk  
             enddo !j
            enddo !i
         enddo !rvecind  
         close(321)

         deallocate(lwindow,u_matrix_opt,ndimwin)
         deallocate(u_matrix,hwann,hreal)
         deallocate(nablamat,nablamat2)
      enddo !jspin

      end subroutine wann_nabla_rs
      end module m_wann_nabla_rs
