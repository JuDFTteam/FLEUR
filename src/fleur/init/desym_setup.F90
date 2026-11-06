!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_desym_setup
   use m_juDFT
   use m_constants
   use m_make_sym
   use m_make_sphhar
   use m_make_stars
   use m_lapwdim
   use m_setupMPI
   use m_types_atoms
   use m_types_cell
   use m_types_desym_map
   use m_types_enpara
   use m_types_fleurinput
   use m_types_forcetheo
   use m_types_gfinp
   use m_types_kpts
   use m_types_lapw
   use m_types_mpi
   use m_types_noco
   use m_types_nococonv
   use m_types_sphhar
   use m_types_stars
   use m_types_sym
   implicit none
   private
   public :: desym_setup

contains

   subroutine desym_setup(fmpi, fi, enpara, nococonv, forcetheo, &
                          fmpi_nosym, fi_nosym, stars_nosym, sphhar_nosym, enpara_nosym, nococonv_nosym, map)
      !! Builds in memory the setup of fi without symmetry: identity only, one type per atom, full k-mesh in bkf order.
      !! Atom indices are kept; map relates both setups.
      type(t_mpi),        intent(in)  :: fmpi
      type(t_fleurinput), intent(in)  :: fi
      type(t_enpara),     intent(in)  :: enpara
      type(t_nococonv),   intent(in)  :: nococonv
      class(t_forcetheo), intent(in)  :: forcetheo
      type(t_mpi),        intent(out) :: fmpi_nosym
      type(t_fleurinput), intent(out) :: fi_nosym
      type(t_stars),      intent(out) :: stars_nosym
      type(t_sphhar),     intent(out) :: sphhar_nosym
      type(t_enpara),     intent(out) :: enpara_nosym
      type(t_nococonv),   intent(out) :: nococonv_nosym
      type(t_desym_map),  intent(out) :: map

      type(t_gfinp) :: gfinp_empty
      integer, allocatable :: parent_type(:)
      integer :: nbasfcn, nvd_save, nv2d_save, nbasfcn_save

      call timestart("desym_setup")

      if (fi%noco%l_noco .or. fi%noco%l_soc) call juDFT_error("Desymmetrization is not implemented for noco/SOC.", calledby="desym_setup")
      if (fi%input%film) call juDFT_error("Desymmetrization is not implemented for films yet.", calledby="desym_setup")
      if (fi%atoms%n_u + fi%atoms%n_hia > 0) call juDFT_error("Desymmetrization is not implemented for LDA+U yet.", calledby="desym_setup")

      fi_nosym = fi
      fi_nosym%gfinp = gfinp_empty

      ! The modified input types are built on rank 0 and broadcast, as in fleur_init
      if (fmpi%irank == 0) then
         call fi%atoms%split_types(fi_nosym%atoms, map%parent_type, map%ldau_parent)
         parent_type = map%parent_type

         call build_sym_nosym(fi%cell, fi%input%film, fi_nosym%sym)
         call make_sym(fi_nosym%sym, fi_nosym%cell, fi_nosym%atoms, fi_nosym%noco, fi_nosym%input, fi_nosym%gfinp)

         map%type_op = fi%sym%ngopr

         call split_noco_types(parent_type, fi_nosym%noco)
         if (allocated(fi%enparaXML%qn_el)) then
            deallocate(fi_nosym%enparaXML%qn_el)
            allocate(fi_nosym%enparaXML%qn_el(0:3, size(parent_type), size(fi%enparaXML%qn_el,3)))
            fi_nosym%enparaXML%qn_el   = fi%enparaXML%qn_el(:,parent_type,:)
            fi_nosym%enparaXML%qn_ello = fi%enparaXML%qn_ello(:,parent_type,:)
         end if

         if (allocated(fi_nosym%vacuum%mrot2)) deallocate(fi_nosym%vacuum%mrot2)
         if (allocated(fi_nosym%vacuum%tau2)) deallocate(fi_nosym%vacuum%tau2)
         call fi_nosym%vacuum%init(fi_nosym%sym)

         call build_kpts(fi, fi_nosym, map)
      end if

      call fi_nosym%sym%mpi_bc(fmpi%mpi_comm, 0)
      call fi_nosym%atoms%mpi_bc(fmpi%mpi_comm, 0)
      call fi_nosym%noco%mpi_bc(fmpi%mpi_comm, 0)
      call fi_nosym%vacuum%mpi_bc(fmpi%mpi_comm, 0)
      call fi_nosym%kpts%mpi_bc(fmpi%mpi_comm, 0)
      call fi_nosym%enparaXML%mpi_bc(fmpi%mpi_comm, 0)
      call map%mpi_bc(fmpi%mpi_comm, 0)
      parent_type = map%parent_type

      fi_nosym%input%l_real = fi%input%l_real .and. fi_nosym%sym%invs

      enpara_nosym = enpara
      ! Arrays starting at l=0 are allocated first, so that the assignment keeps their lower bound
      deallocate(enpara_nosym%el0, enpara_nosym%el1, enpara_nosym%lchange, enpara_nosym%qn_el)
      allocate(enpara_nosym%el0(0:fi%atoms%lmaxd, size(parent_type), size(enpara%el0,3)))
      allocate(enpara_nosym%el1(0:fi%atoms%lmaxd, size(parent_type), size(enpara%el1,3)))
      allocate(enpara_nosym%lchange(0:fi%atoms%lmaxd, size(parent_type), size(enpara%lchange,3)))
      allocate(enpara_nosym%qn_el(0:3, size(parent_type), size(enpara%qn_el,3)))
      enpara_nosym%el0     = enpara%el0(:,parent_type,:)
      enpara_nosym%el1     = enpara%el1(:,parent_type,:)
      enpara_nosym%ello0   = enpara%ello0(:,parent_type,:)
      enpara_nosym%ello1   = enpara%ello1(:,parent_type,:)
      enpara_nosym%skiplo  = enpara%skiplo(parent_type,:)
      enpara_nosym%lchange = enpara%lchange(:,parent_type,:)
      enpara_nosym%llochg  = enpara%llochg(:,parent_type,:)
      enpara_nosym%qn_el   = enpara%qn_el(:,parent_type,:)
      enpara_nosym%qn_ello = enpara%qn_ello(:,parent_type,:)
      if (allocated(enpara%vr)) enpara_nosym%vr = enpara%vr(:,parent_type,:)
      if (allocated(enpara%vr_core)) enpara_nosym%vr_core = enpara%vr_core(:,parent_type,:)

      nococonv_nosym = nococonv
      if (allocated(nococonv%alph))     nococonv_nosym%alph     = nococonv%alph(parent_type)
      if (allocated(nococonv%beta))     nococonv_nosym%beta     = nococonv%beta(parent_type)
      if (allocated(nococonv%alphRlx))  nococonv_nosym%alphRlx  = nococonv%alphRlx(parent_type)
      if (allocated(nococonv%betaRlx))  nococonv_nosym%betaRlx  = nococonv%betaRlx(parent_type)
      if (allocated(nococonv%alphPrev)) nococonv_nosym%alphPrev = nococonv%alphPrev(parent_type)
      if (allocated(nococonv%betaPrev)) nococonv_nosym%betaPrev = nococonv%betaPrev(parent_type)
      if (allocated(nococonv%b_con))    nococonv_nosym%b_con    = nococonv%b_con(:,parent_type)

      ! No storeStructureIfNew/transform_by_moving_atoms: the setup must not touch cdn.hdf
      call make_sphhar(.false., fi_nosym%atoms, sphhar_nosym, fi_nosym%sym, fi_nosym%cell)
      call make_stars(stars_nosym, fi_nosym%sym, fi_nosym%atoms, fi_nosym%vacuum, sphhar_nosym, fi_nosym%input, fi_nosym%cell, fi_nosym%noco, fmpi)

      nvd_save = lapw_dim_nvd
      nv2d_save = lapw_dim_nv2d
      nbasfcn_save = lapw_dim_nbasfcn
      ! Needed for DFPT: k+q with the full k-mesh and the irreducible q-points can need more plane waves than with the IBZ k-points
      call lapw_dim(fi_nosym%kpts, fi_nosym%cell, fi_nosym%input, fi_nosym%noco, nococonv_nosym, forcetheo, fi_nosym%atoms, nbasfcn, fi_nosym%dfpt)
      if (lapw_dim_nvd < nvd_save .or. lapw_dim_nv2d < nv2d_save .or. lapw_dim_nbasfcn < nbasfcn_save) &
         call juDFT_error("LAPW basis dimensions of the desymmetrized setup are smaller than in the symmetric setup.", calledby="desym_setup")

      fmpi_nosym = fmpi
      if (allocated(fmpi_nosym%k_list)) deallocate(fmpi_nosym%k_list)
      if (allocated(fmpi_nosym%ev_list)) deallocate(fmpi_nosym%ev_list)
      if (allocated(fmpi_nosym%coulomb_owner)) deallocate(fmpi_nosym%coulomb_owner)
      call setupMPI(fi_nosym%kpts%nkpt, fi_nosym%input%neig, nbasfcn, fmpi_nosym, fi_nosym%input%l_real, fi_nosym%noco%l_noco)

      if (fmpi%irank == 0) then
         write(oUnit, '(a)') 'Desymmetrized setup:'
         write(oUnit, '(a,i0,a,i0)') '  symmetry operations: ', fi%sym%nop, ' -> ', fi_nosym%sym%nop
         write(oUnit, '(a,i0,a,i0)') '  atom types:          ', fi%atoms%ntype, ' -> ', fi_nosym%atoms%ntype
         write(oUnit, '(a,i0,a,i0)') '  k-points:            ', fi%kpts%nkpt, ' -> ', fi_nosym%kpts%nkpt
         write(oUnit, '(a,i0)')      '  stars:               ', stars_nosym%ng3
      end if

      call timestop("desym_setup")

   end subroutine desym_setup

   subroutine build_sym_nosym(cell, film, sym_nosym)
      !! Sets up sym_nosym with the identity as the only operation.
      type(t_cell), intent(in)  :: cell
      logical,      intent(in)  :: film
      type(t_sym),  intent(out) :: sym_nosym

      integer :: i

      sym_nosym%nop  = 1
      sym_nosym%nop2 = 1
      allocate(sym_nosym%mrot(3,3,1), sym_nosym%tau(3,1))
      sym_nosym%mrot = 0
      do i = 1, 3
         sym_nosym%mrot(i,i,1) = 1
      end do
      sym_nosym%tau = 0.0
      call sym_nosym%init(cell, film)

   end subroutine build_sym_nosym

   subroutine split_noco_types(parent_type, noco)
      !! Re-indexes the per-type arrays of noco from the old to the new atom types.
      integer,      intent(in)    :: parent_type(:)
      type(t_noco), intent(inout) :: noco

      if (allocated(noco%mix_RelaxWeightOffD)) noco%mix_RelaxWeightOffD = noco%mix_RelaxWeightOffD(parent_type)
      if (allocated(noco%l_constrained))   noco%l_constrained   = noco%l_constrained(parent_type)
      if (allocated(noco%l_fixedMoment))   noco%l_fixedMoment   = noco%l_fixedMoment(parent_type)
      if (allocated(noco%l_unrestrictMT))  noco%l_unrestrictMT  = noco%l_unrestrictMT(parent_type)
      if (allocated(noco%l_spinoffd_ldau)) noco%l_spinoffd_ldau = noco%l_spinoffd_ldau(parent_type)
      if (allocated(noco%l_alignMT))       noco%l_alignMT       = noco%l_alignMT(parent_type)
      if (allocated(noco%alph_inp))        noco%alph_inp        = noco%alph_inp(parent_type)
      if (allocated(noco%beta_inp))        noco%beta_inp        = noco%beta_inp(parent_type)
      if (allocated(noco%socscale))        noco%socscale        = noco%socscale(parent_type)

   end subroutine split_noco_types

   subroutine build_kpts(fi, fi_nosym, map)
      !! Sets the k-points of fi_nosym to the full k-mesh of fi in bkf order with equal weights.
      type(t_fleurinput), intent(in)    :: fi
      type(t_fleurinput), intent(inout) :: fi_nosym
      type(t_desym_map),  intent(inout) :: map

      integer :: nkptf

      if (fi%kpts%nkptf == 0) call juDFT_error("Full k-mesh of the symmetric setup is missing.", calledby="desym_setup")
      nkptf = fi%kpts%nkptf

      associate(kpts_nosym => fi_nosym%kpts)
         kpts_nosym%nkpt  = nkptf
         kpts_nosym%bk    = fi%kpts%bkf(:,:nkptf)
         kpts_nosym%wtkpt = spread(1.0 / real(nkptf), 1, nkptf)
         kpts_nosym%nkptf = 0
         kpts_nosym%ntet  = 0
         if (allocated(kpts_nosym%bkf))       deallocate(kpts_nosym%bkf)
         if (allocated(kpts_nosym%bkp))       deallocate(kpts_nosym%bkp)
         if (allocated(kpts_nosym%bksym))     deallocate(kpts_nosym%bksym)
         if (allocated(kpts_nosym%EIBZ))      deallocate(kpts_nosym%EIBZ)
         if (allocated(kpts_nosym%ntetra))    deallocate(kpts_nosym%ntetra)
         if (allocated(kpts_nosym%voltet))    deallocate(kpts_nosym%voltet)
         if (allocated(kpts_nosym%tetraList)) deallocate(kpts_nosym%tetraList)
         kpts_nosym%l_set_eibz = .false.
         call kpts_nosym%init(fi_nosym%sym, fi_nosym%input%film, .false., .false.)
         call kpts_nosym%initTetra(fi_nosym%input, fi_nosym%cell, fi_nosym%sym, fi_nosym%noco%l_soc .or. fi_nosym%noco%l_ss)
         if (kpts_nosym%nkptf /= nkptf) call juDFT_error("k-mesh of the desymmetrized setup is not closed.", calledby="desym_setup")
      end associate

      map%kpt_parent = fi%kpts%bkp(:nkptf)
      map%kpt_op     = fi%kpts%bksym(:nkptf)

   end subroutine build_kpts

end module m_desym_setup
