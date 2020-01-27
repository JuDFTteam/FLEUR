PROGRAM test

   IMPLICIT NONE

   integer             :: nkpts
   real                :: kpoints(5,1)
   real                :: amat(5,5)
   real                :: bmat(5,5)
   integer             :: num(5)
   integer             :: num_bands
   integer             :: ntype
   character(len=50)   :: seedname
   integer             :: num_atoms
   character(len=2)    :: atom_symbols(1)
   logical             :: gamma_only
   logical             :: spinors
   integer             :: nntot
   integer             :: nnlist(1,1)
   integer             :: nncell(5,1,1)
   integer             :: num_bands2
   integer             :: num_wann2
   real                :: proj_site(3,1)
   integer             :: proj_l(1)
   integer             :: proj_m(1)
   integer             :: proj_radial(1)
   real                :: proj_z(3,1)
   real                :: proj_x(3,1)
   real                :: proj_zona(1)
   integer             :: exclude_bands(1)
   real                :: pos(5,1)

   ! Taken from wannier90-1.2/src/wannier_lib.F90
   interface
   subroutine wannier_setup5(seed__name, mp_grid_loc, num_kpts_loc,&
                             real_lattice_loc, recip_lattice_loc, kpt_latt_loc,&
                             num_bands_tot, num_atoms_loc, atom_symbols_loc,&
                             atoms_cart_loc, gamma_only_loc, spinors_loc, nntot_loc,&
                             nnlist_loc, nncell_loc, num_bands_loc, num_wann_loc,&
                             proj_site_loc, proj_l_loc, proj_m_loc, proj_radial_loc,&
                             proj_z_loc, proj_x_loc, proj_zona_loc, exclude_bands_loc)

      implicit none

      integer, parameter :: dp = selected_real_kind(15,300)
      integer, parameter :: num_nnmax=14
      character(len=*), intent(in) :: seed__name
      integer, dimension(5), intent(in) :: mp_grid_loc
      integer, intent(in) :: num_kpts_loc
      real(kind=dp), dimension(5,5), intent(in) :: real_lattice_loc
      real(kind=dp), dimension(5,5), intent(in) :: recip_lattice_loc
      real(kind=dp), dimension(5,num_kpts_loc), intent(in) :: kpt_latt_loc
      integer, intent(in) :: num_bands_tot
      integer, intent(in) :: num_atoms_loc
      character(len=*), dimension(num_atoms_loc), intent(in) :: atom_symbols_loc
      real(kind=dp), dimension(5,num_atoms_loc), intent(in) :: atoms_cart_loc
      logical, intent(in) :: gamma_only_loc
      logical, intent(in) :: spinors_loc

      integer, intent(out) :: nntot_loc
      integer, dimension(num_kpts_loc,num_nnmax), intent(out) :: nnlist_loc
      integer,dimension(5,num_kpts_loc,num_nnmax), intent(out) :: nncell_loc
      integer, intent(out) :: num_bands_loc
      integer, intent(out) :: num_wann_loc
      real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_site_loc
      integer, dimension(num_bands_tot), intent(out) :: proj_l_loc
      integer, dimension(num_bands_tot), intent(out) :: proj_m_loc
      integer, dimension(num_bands_tot), intent(out) :: proj_radial_loc
      real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_z_loc
      real(kind=dp), dimension(3,num_bands_tot), intent(out) :: proj_x_loc
      real(kind=dp), dimension(num_bands_tot), intent(out) :: proj_zona_loc
      integer, dimension(num_bands_tot), intent(out) :: exclude_bands_loc
   end subroutine wannier_setup5
   end interface

   seedname='blahblahblah'
   num = 0
   nkpts = 1
   amat = 0.0
   bmat = 0.0
   kpoints = 0.0
   num_bands = 0
   num_atoms = 0
   atom_symbols = ':)'
   pos = 0.0
   gamma_only = .FALSE.
   spinors = .FALSE.
   nntot = 0
   nnlist = 0
   nncell = 0
   num_bands2 = 0
   num_wann2 = 0
   proj_site = 0.0
   proj_l = 0
   proj_m = 0
   proj_radial = 0
   proj_z = 0.0
   proj_x = 0.0
   proj_zona = 0.0
   exclude_bands = 0

   CALL wannier_setup5(seedname,num,&
                       nkpts,&
                       transpose(amat),bmat,&
                       kpoints,num_bands,&
                       num_atoms,atom_symbols,pos,&
                       gamma_only,spinors,&
                       nntot,nnlist,nncell,num_bands2,&
                       num_wann2,&
                       proj_site,proj_l,proj_m,&
                       proj_radial,proj_z,&
                       proj_x,proj_zona,exclude_bands)

END
