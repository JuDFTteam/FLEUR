      MODULE m_wann_postproc_setup
      CONTAINS
      SUBROUTINE wann_postproc_setup(
     >           natd,nkpts,kpoints,amat,bmat,
     >           num,num_bands,ntype,neq,
     >           zatom,taual,namat)
c******************************************************
c     Call the Wannier90 setup routine.
c     Purpose: get the file of nearest-neighbor kpoints
c     Frank Freimuth
c******************************************************
      IMPLICIT NONE

      integer, intent(in)         :: natd
      integer, intent(in)         :: nkpts
      real,    intent(in)         :: kpoints(3,nkpts)
      real,    intent(in)         :: amat(3,3)
      real,    intent(in)         :: bmat(3,3)
      integer, intent(in)         :: num(3)
      integer, intent(in)         :: num_bands
      integer, intent(in)         :: ntype
      integer, intent(in)         :: neq(ntype)
      real,    intent(in)         :: zatom(ntype)
      real,    intent(in)         :: taual(:,:)
      character(len=2),intent(in) :: namat(0:103)

      integer             :: i,j,at
      character(len=50)   :: seedname
      integer             :: num_atoms
      character(len=2)    :: atom_symbols(natd)
      logical             :: gamma_only
      logical             :: spinors
      integer,parameter   :: num_nnmax=12
      integer             :: nntot
      integer             :: nnlist(nkpts,num_nnmax)
      integer             :: nncell(3,nkpts,num_nnmax)
      integer             :: num_bands2
      integer             :: num_wann2
      real                :: proj_site(3,num_bands)
      integer             :: proj_l(num_bands)
      integer             :: proj_m(num_bands)
      integer             :: proj_radial(num_bands)
      real                :: proj_z(3,num_bands)
      real                :: proj_x(3,num_bands)
      real                :: proj_zona(num_bands)
      integer             :: exclude_bands(num_bands)
      integer             :: nn,ikpt
      real                :: pos(3,natd)

      ! Taken from wannier90-1.2/src/wannier_lib.F90
      interface
        subroutine wannier_setup(seed__name, mp_grid_loc, num_kpts_loc,
     +    real_lattice_loc, recip_lattice_loc, kpt_latt_loc,
     +    num_bands_tot, num_atoms_loc, atom_symbols_loc,
     +    atoms_cart_loc, gamma_only_loc, spinors_loc, nntot_loc,
     +    nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc,
     +    proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc,
     +    proj_z_loc, proj_x_loc, proj_zona_loc, exclude_bands_loc)
         implicit none
         integer, parameter :: dp = selected_real_kind(15,300)
         integer, parameter :: num_nnmax=12
         character(len=*), intent(in) :: seed__name
         integer, dimension(3), intent(in) :: mp_grid_loc
         integer, intent(in) :: num_kpts_loc
         real(kind=dp), dimension(3,3), intent(in) :: real_lattice_loc
         real(kind=dp), dimension(3,3), intent(in) :: recip_lattice_loc
         real(kind=dp), dimension(3,num_kpts_loc), intent(in) ::
     +      kpt_latt_loc
         integer, intent(in) :: num_bands_tot
         integer, intent(in) :: num_atoms_loc
         character(len=*), dimension(num_atoms_loc), intent(in) ::
     +      atom_symbols_loc
         real(kind=dp), dimension(3,num_atoms_loc), intent(in) ::
     +      atoms_cart_loc
         logical, intent(in) :: gamma_only_loc
         logical, intent(in) :: spinors_loc
         integer, intent(out) :: nntot_loc
         integer, dimension(num_kpts_loc,num_nnmax), intent(out) ::
     +      nnlist_loc
         integer,dimension(3,num_kpts_loc,num_nnmax), intent(out) ::
     +      nncell_loc
         integer, intent(out) :: num_bands_loc
         integer, intent(out) :: num_wann_loc
         real(kind=dp), dimension(3,num_bands_tot), intent(out) ::
     +      proj_site_loc
         integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
         integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
         integer, dimension(num_bands_tot), intent(out) ::
     +      proj_radial_loc
         real(kind=dp), dimension(3,num_bands_tot), intent(out) ::
     +      proj_z_loc
         real(kind=dp), dimension(3,num_bands_tot), intent(out) ::
     +      proj_x_loc
         real(kind=dp), dimension(num_bands_tot), intent(out) ::
     +      proj_zona_loc
         integer, dimension(num_bands_tot), intent(out) ::
     +      exclude_bands_loc
        end subroutine wannier_setup
      end interface

      do j=1,natd
         pos(:,j)=matmul(amat(:,:),taual(:,j))
      enddo

      seedname='WF1'
      gamma_only=.false.
      spinors=.false.

      num_atoms=0
      do i=1,ntype
         at=nint(zatom(i))
         do j=1,neq(i)
            num_atoms=num_atoms+1
            atom_symbols(num_atoms)=namat(at)
         enddo !j
      enddo !i

c**********************************************************
c     Call Wannier90 routine for preparation.
c**********************************************************
      call wannier_setup(
     >     seedname,num,
     >     nkpts,
     >     transpose(amat),bmat,
     >     kpoints,num_bands, 
     >     num_atoms,atom_symbols,pos, 
     >     gamma_only,spinors, 
     >     nntot,nnlist,nncell,num_bands2,
     >     num_wann2,
     >     proj_site,proj_l,proj_m,
     >     proj_radial,proj_z, 
     >     proj_x,proj_zona,exclude_bands)

c******************************************************
c           write bkpts
c******************************************************
      open(202,file='bkpts',form='formatted')
      write (202,'(i4)') nntot
      do ikpt=1,nkpts
        do nn=1,nntot
          write (202,'(2i6,3x,3i4)') 
c     &     ikpt,bpt(nn,ikpt),(gb(i,nn,ikpt),i=1,3)
     &     ikpt,nnlist(ikpt,nn),(nncell(i,ikpt,nn),i=1,3)
        enddo
      enddo
      close(202)

      END SUBROUTINE wann_postproc_setup
      END MODULE m_wann_postproc_setup
