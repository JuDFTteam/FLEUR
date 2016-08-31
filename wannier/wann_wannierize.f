!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_wann_wannierize
      use m_juDFT
c****************************************************
c     Call wannier90 subroutines needed for
c     wannierization from within Fleur.
c
c     Frank Freimuth
c**************************************************** 
      CONTAINS
      SUBROUTINE wann_wannierize(
     >          film,l_bzsym,jspins,
     >          natd,pos,
     >          amat,bmat,ntype,neq,zatom)

      use m_wann_read_umatrix
      implicit none
      logical,intent(in)  :: film
      logical,intent(in)  :: l_bzsym
      integer,intent(in)  :: jspins
      integer,intent(in)  :: natd
      real,intent(in)     :: pos(3,natd)
      real,intent(in)     :: amat(3,3)
      real,intent(in)     :: bmat(3,3)
      integer,intent(in)  :: ntype
      integer,intent(in)  :: neq(ntype)
      real,intent(in)     :: zatom(ntype)

      integer             :: num_wann
      real,allocatable    :: centers(:,:)      
      integer             :: nkpts,dim,i,j,at
      character(len=50)   :: seedname
      integer             :: jspin
      character(len=3)    :: spin12(2)
      integer             :: nntot
      data spin12/'WF1' , 'WF2'/
      integer             :: num(3)
      integer             :: num_kpts,num_wann2
      real                :: real_lattice(3,3)
      real                :: recip_lattice(3,3)
      integer             :: num_bands
      integer             :: num_atoms
      real                :: atoms_cart(3,natd)
      character(len=2)    :: atom_symbols(natd)
      logical             :: gamma_only
      complex,allocatable :: M_matrix(:,:,:,:)
      complex,allocatable :: A_matrix(:,:,:)
      real,allocatable    :: eigenvalues(:,:)
      complex,allocatable :: U_matrix(:,:,:)
      complex,allocatable :: U_matrix_opt(:,:,:)
      logical,allocatable :: lwindow(:,:)
      integer,allocatable :: ndimwin(:)
      real,allocatable    :: wann_spreads(:)
      real                :: spread(3),maxi,mini
      logical             :: l_file,l_bkpts
      integer             :: iter
      real                :: increm,compare
      real,allocatable    :: kpoints(:,:)
      real,parameter      :: bohr=0.5291772108
      character(len=2)    :: namat(0:103)
      real                :: realp,imagp
      real                :: scale
      integer             :: ikpt,ikpt_b,nwf,nwf2,i2,ikpt2
      logical             :: have_disentangled

      ! Taken from wannier90-1.2/src/wannier_lib.F90
      interface
       subroutine wannier_run(seed__name,mp_grid_loc,num_kpts_loc,
     +    real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_loc,
     +    num_wann_loc,nntot_loc,num_atoms_loc,atom_symbols_loc,
     +    atoms_cart_loc,gamma_only_loc,M_matrix_loc,A_matrix_loc,
     +    eigenvalues_loc,
     +    U_matrix_loc,U_matrix_opt_loc,lwindow_loc,wann_centres_loc,
     +    wann_spreads_loc,spread_loc)
         implicit none
         integer, parameter :: dp = selected_real_kind(15,300)
         character(len=*), intent(in) :: seed__name
         integer, dimension(3), intent(in) :: mp_grid_loc
         integer, intent(in) :: num_kpts_loc
         real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
         real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
         real(kind=dp), dimension(3,num_kpts_loc), intent(in) ::
     +         kpt_latt_loc
         integer, intent(in) :: num_bands_loc
         integer, intent(in) :: num_wann_loc
         integer, intent(in) :: nntot_loc
         integer, intent(in) :: num_atoms_loc
         character(len=*), dimension(num_atoms_loc), intent(in) ::
     +         atom_symbols_loc
         real(kind=dp), dimension(3,num_atoms_loc), intent(in) ::
     +         atoms_cart_loc
         logical, intent(in) :: gamma_only_loc
         complex(kind=dp), dimension(num_bands_loc,num_bands_loc,
     +         nntot_loc,num_kpts_loc), intent(in) :: M_matrix_loc
         complex(kind=dp),
     +      dimension(num_bands_loc,num_wann_loc,num_kpts_loc),
     +      intent(in) :: A_matrix_loc
         real(kind=dp), dimension(num_bands_loc,num_kpts_loc),
     +      intent(in) :: eigenvalues_loc
         complex(kind=dp),
     +      dimension(num_wann_loc,num_wann_loc,num_kpts_loc),
     +      intent(out) :: U_matrix_loc
         complex(kind=dp),
     +      dimension(num_bands_loc,num_wann_loc,num_kpts_loc),
     +      optional, intent(out) :: U_matrix_opt_loc
         logical, dimension(num_bands_loc,num_kpts_loc),
     +      optional, intent(out) :: lwindow_loc
         real(kind=dp), dimension(3,num_wann_loc),
     +      optional, intent(out) :: wann_centres_loc
         real(kind=dp), dimension(num_wann_loc), optional,
     +      intent(out) :: wann_spreads_loc
         real(kind=dp), dimension(3), optional,
     +      intent(out) :: spread_loc
       end subroutine wannier_run
      end interface

      DATA namat/'va',' h','he','li','be',' b',' c',' n',' o',' f','ne',
     +     'na','mg','al','si',' p',' s','cl','ar',' k','ca','sc','ti',
     +     ' v','cr','mn','fe','co','ni','cu','zn','ga','ge','as','se',
     +     'br','kr','rb','sr',' y','zr','nb','mo','tc','ru','rh','pd',
     +     'ag','cd','in','sn','sb','te',' j','xe','cs','ba','la','ce',
     +     'pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     +     'lu','hf','ta',' w','re','os','ir','pt','au','hg','tl','pb',
     +     'bi','po','at','rn','fr','ra','ac','th','pa',' u','np','pu',
     +     'am','cm','bk','cf','es','fm','md','no','lw'/


      gamma_only=.false.
      atoms_cart=pos*bohr
c**********************************************************
c     read the bkpts file
c**********************************************************
      l_bkpts = .false.
      inquire (file='bkpts',exist=l_bkpts)
      IF (.NOT.l_bkpts)  CALL juDFT_error("need bkpts for matrixmmn"
     +     ,calledby ="wann_wannierize")
      open (202,file='bkpts',form='formatted',status='old')
      rewind (202)
      read (202,'(i4)') nntot
      close(202)

c**********************************************************
c     information on atoms
c**********************************************************
      num_atoms=0
      do i=1,ntype
         at=nint(zatom(i))
         do j=1,neq(i)
            num_atoms=num_atoms+1
            atom_symbols(num_atoms)=namat(at)
         enddo !j
      enddo !i

c**********************************************************
c     read in kpoints from kpts/w90kpts file
c**********************************************************
      if(l_bzsym)then
         l_file=.false.
         inquire(file='w90kpts',exist=l_file)
         IF(.NOT.l_file) CALL juDFT_error("where is w90kpts?",calledby
     +        ="wann_wannierize")
         open(987,file='w90kpts',status='old',form='formatted')
         read(987,*)nkpts,scale
         print*,"nkpts=",nkpts
         allocate(kpoints(3,nkpts))
         do iter=1,nkpts
           read(987,*)kpoints(:,iter)
         enddo
         close(987)
         do iter=1,nkpts
           print*,kpoints(:,iter)
         enddo
         kpoints=kpoints/scale
      else
         l_file=.false.
         inquire(file='kpts',exist=l_file)
         IF(.NOT.l_file) CALL juDFT_error("where is kpts?",calledby
     +        ="wann_wannierize")
         open(987,file='kpts',status='old',form='formatted')
         read(987,*)nkpts,scale
         allocate(kpoints(3,nkpts))
         do iter=1,nkpts
            read(987,*)kpoints(:,iter)
         enddo   
         close(987)
         if(film) kpoints(3,:)=0.0
         kpoints=kpoints/scale
         do iter=1,nkpts
            print*,kpoints(:,iter)
         enddo
      endif
      num_kpts=nkpts
      allocate(ndimwin(num_kpts))
c*********************************************************
c           find out the structure of k-point set
c*********************************************************
      do dim=1,3
         maxi=maxval(kpoints(dim,:))
         mini=minval(kpoints(dim,:))
         if(mini==maxi)then
            num(dim)=1
         else   
            increm=maxi-mini
            do iter=1,nkpts
               compare=maxi-kpoints(dim,iter)
               if(abs(compare).lt.1e-6)cycle
               if(compare.lt.increm) then
                  increm=compare
               endif   
            enddo
            num(dim)=(maxi-mini)/increm+1.01
         endif   
      enddo
      print*,"num(:)=",num(:)
      IF(num(1)*num(2)*num(3)/=nkpts)  CALL juDFT_error
     +     ("mysterious kpoints",calledby ="wann_wannierize")

c********************************************************
c        proj file provides num_wann and num_bands
c********************************************************
      l_file=.false.
      inquire(file='proj',exist=l_file)
      IF(.NOT.l_file)  CALL juDFT_error("where is proj?",calledby
     +     ="wann_wannierize")
      open(712,file='proj',form='formatted',status='old')
      read(712,*)num_wann2,num_bands
      close(712)
      num_wann=num_wann2
      print*,"num_wann=",num_wann
      print*,"num_bands=",num_bands

      real_lattice  = transpose(amat)*bohr
      recip_lattice = bmat/bohr

      allocate( M_matrix(num_bands,num_bands,nntot,num_kpts) )
      allocate( A_matrix(num_bands,num_wann,num_kpts) )
      allocate( eigenvalues(num_bands,num_kpts) )
      allocate( U_matrix(num_wann,num_wann,num_kpts) )
      allocate( U_matrix_opt(num_bands,num_wann,num_kpts) )
      allocate( lwindow(num_bands,num_kpts) )
      do jspin=1,jspins
         seedname=spin12(jspin)
c******** read mmn-matrix
         open (305,file=spin12(jspin)//'.mmn',
     &             form='formatted',status='old')
         read (305,*)
         read (305,'(3i5)') 
         do ikpt = 1,num_kpts
          do ikpt_b = 1,nntot
           read (305,'(2i5,3x,3i4)') 
           do i = 1,num_bands
            do j = 1,num_bands
             read (305,*)realp,imagp
             m_matrix(j,i,ikpt_b,ikpt)=cmplx(realp,imagp)
            enddo !j
           enddo !i 
          enddo !ikpt_b
         enddo !ikpt
         close(305)
c******** read amn-matrix
         open (303,file=spin12(jspin)//'.amn',
     &             form='formatted',status='old')
         read (303,*) 
         read (303,'(3i5)') 
         do ikpt = 1,num_kpts
          do nwf = 1,num_wann
           do i = 1,num_bands
c            print*,"ikpt=",ikpt,"nwf=",nwf,"i=",i  
            read (303,'(3i5,3x,2f18.12)') i2,nwf2,ikpt2,realp,imagp
               a_matrix(i,nwf,ikpt)=cmplx(realp,imagp)
c            print*,"i2=",i2,"nwf2=",nwf2,"ikpt2=2",ikpt2
           enddo
          enddo
         enddo
         close (303)
c********* read eigenvalues
         open(303,file=spin12(jspin)//'.eig',
     &            form='formatted',status='old')
         do ikpt=1,num_kpts
            do i=1,num_bands
               read(303,*)nwf2,ikpt2,eigenvalues(i,ikpt)
            enddo
         enddo
         allocate( centers(3,num_wann) )
         allocate( wann_spreads(num_wann) )
         call wannier_run(
     >       seedname,num,num_kpts, 
     >       real_lattice,recip_lattice,kpoints,num_bands, 
     >       num_wann,nntot,num_atoms,atom_symbols, 
     >       atoms_cart,gamma_only,M_matrix,A_matrix,eigenvalues, 
     >       U_matrix,U_matrix_opt,lwindow,
     <       centers(:,:), 
     <       wann_spreads,spread)
         deallocate( centers ) 
         deallocate( wann_spreads )
c********read the u_matrix and write it to a formatted file
         call wann_read_umatrix(
     >       num_kpts,num_wann,num_bands,
     >       .true.,jspin,1,
     <       have_disentangled,
     <       lwindow,ndimwin,u_matrix_opt)

      enddo !jspin

      END SUBROUTINE wann_wannierize
      END MODULE m_wann_wannierize

