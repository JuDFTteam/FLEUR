!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_metagga
   PUBLIC  :: calc_EnergyDen
   PRIVATE :: calc_EnergyDen_auxillary_weights, &
              calc_kinEnergyDen_pw, &
              calc_kinEnergyDen_mt

   type t_RS_potden
      REAL, ALLOCATABLE :: is(:,:), mt(:,:)
   end type t_RS_potden

CONTAINS
   SUBROUTINE calc_kinEnergyDen_pw(EnergyDen_rs, vTot_rs, den_rs, kinEnergyDen_RS)
      USE m_juDFT_stop
      !use m_cdngen
      IMPLICIT NONE
      REAL, INTENT(in)                 :: den_RS(:,:), EnergyDen_RS(:,:), vTot_RS(:,:)
      REAL, INTENT(inout), allocatable :: kinEnergyDen_RS(:,:)
#ifdef CPP_LIBXC
      REAL, PARAMETER                  :: eps = 1e-15

      kinEnergyDen_RS = EnergyDen_RS - vTot_RS * den_RS
#else
      CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
   END SUBROUTINE calc_kinEnergyDen_pw

   SUBROUTINE calc_kinEnergyDen_mt(EnergyDen_RS, vTot_rs, vTot0_rs, core_den_rs, val_den_rs, &
                                   atm_idx, nsp, kinEnergyDen_RS)
      USE m_juDFT_stop
      USE m_juDFT_string
      implicit none
      REAL, INTENT(in)                 :: EnergyDen_RS(:,:), vTot_rs(:,:), vTot0_rs(:,:), core_den_rs(:,:), val_den_rs(:,:)
      INTEGER, intent(in)              :: atm_idx, nsp
      REAL, INTENT(inout)              :: kinEnergyDen_RS(:,:)

#ifdef CPP_LIBXC
      kinEnergyDen_RS = EnergyDen_RS - (vTot0_rs * core_den_rs + vTot_rs * val_den_rs)
#else
      CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
   END SUBROUTINE calc_kinEnergyDen_mt


   SUBROUTINE calc_EnergyDen(eig_id, mpi, kpts, noco, input, banddos, cell, atoms, enpara, stars, &
         vacuum, DIMENSION, sphhar, sym, vTot, oneD, results, EnergyDen)
      ! calculates the energy density
      ! EnergyDen = \sum_i n_i(r) \varepsilon_i
      ! where n_i(r) is the one-particle density
      ! and \varepsilon_i are the eigenenergies

      USE m_types_setup
      USE m_types_potden
      USE m_types_kpts
      USE m_types_mpi
      USE m_types_enpara
      USE m_types_misc
      USE m_types_regionCharges
      USE m_types_dos
      USE m_types_cdnval
      USE m_cdnval

      IMPLICIT NONE

      INTEGER,           INTENT(in)           :: eig_id
      TYPE(t_mpi),       INTENT(in)           :: mpi
      TYPE(t_kpts),      INTENT(in)           :: kpts
      TYPE(t_noco),      INTENT(in)           :: noco
      TYPE(t_input),     INTENT(in)           :: input
      TYPE(t_banddos),   INTENT(in)           :: banddos
      TYPE(t_cell),      INTENT(in)           :: cell
      TYPE(t_atoms),     INTENT(in)           :: atoms
      TYPE(t_enpara),    INTENT(in)           :: enpara
      TYPE(t_stars),     INTENT(in)           :: stars
      TYPE(t_vacuum),    INTENT(in)           :: vacuum
      TYPE(t_dimension), INTENT(in)           :: DIMENSION
      TYPE(t_sphhar),    INTENT(in)           :: sphhar
      TYPE(t_sym),       INTENT(in)           :: sym
      TYPE(t_potden),    INTENT(in)           :: vTot
      TYPE(t_oneD),      INTENT(in)           :: oneD
      TYPE(t_results),   INTENT(in)           :: results
      TYPE(t_potden),    INTENT(inout)        :: EnergyDen

      ! local
      INTEGER                         :: jspin

      TYPE(t_regionCharges)           :: regCharges
      TYPE(t_dos)                     :: dos
      TYPE(t_moments)                 :: moments
      TYPE(t_results)                 :: tmp_results
      TYPE(t_cdnvalJob)               :: cdnvalJob
      TYPE(t_potden)                  :: aux_den, real_den


      CALL regCharges%init(input, atoms)
      CALL dos%init(input,        atoms, DIMENSION, kpts, vacuum)
      CALL moments%init(input,    atoms)
      tmp_results = results

      DO jspin = 1,input%jspins
         CALL cdnvalJob%init(mpi,input,kpts,noco,results,jspin)


         ! replace brillouin weights with auxillary weights
         CALL calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, cdnvalJob%weights)

         CALL cdnval(eig_id, mpi, kpts, jspin, noco, input, banddos, cell, atoms, &
            enpara, stars, vacuum, DIMENSION, sphhar, sym, vTot, oneD, cdnvalJob, &
            EnergyDen, regCharges, dos, tmp_results, moments)
      ENDDO

   END SUBROUTINE calc_EnergyDen

   SUBROUTINE calc_EnergyDen_auxillary_weights(eig_id, kpts, jspin, f_ik)
      USE m_types_kpts
      USE m_eig66_io
      IMPLICIT NONE
      ! calculates new (auxillary-)weights as
      ! f_iks = w_iks * E_iks
      !, where  f_iks are the new (auxillary-)weights
      ! w_iks are the weights used in brillouin zone integration
      ! E_iks are the eigen energies

      INTEGER,      INTENT(in)        :: eig_id
      INTEGER,      INTENT(in)        :: jspin
      TYPE(t_kpts), INTENT(in)        :: kpts
      REAL,         INTENT(inout)     :: f_ik(:,:) ! f_ik(band_idx, kpt_idx)

      ! local vars
      REAL                       :: e_i(SIZE(f_ik,dim=1))
      INTEGER                    :: ikpt

      DO ikpt = 1,kpts%nkpt
         CALL read_eig(eig_id,ikpt,jspin, eig=e_i)
         f_ik(:,ikpt) = f_ik(:,ikpt) * e_i
      ENDDO
   END SUBROUTINE calc_EnergyDen_auxillary_weights

   subroutine set_zPrime(dim_idx, zMat, kpt, lapw, cell, zPrime)
      USE m_types
      USE m_constants
      implicit none
      INTEGER, intent(in)      :: dim_idx
      TYPE (t_mat), intent(in) :: zMat
      REAL, intent(in)         :: kpt(3) 
      TYPE(t_lapw), intent(in) :: lapw
      TYPE(t_cell), intent(in) :: cell
      TYPE (t_mat)             :: zPrime

      REAL                     :: k_plus_g(3), fac
      INTEGER                  :: basis_idx

      call zPrime%free()
      call zPrime%init(zMat)

      do basis_idx = 1,size(lapw%gvec,dim=2)
         k_plus_g = kpt + lapw%gvec(:,basis_idx,1)
         k_plus_g = internal_to_rez(cell, k_plus_g)

         fac = k_plus_g(dim_idx)
         if(zPrime%l_real) then
            zPrime%data_r(basis_idx,:) =            fac * zMat%data_r(basis_idx,:) 
         else
            zPrime%data_c(basis_idx,:) = ImagUnit * fac * zMat%data_c(basis_idx,:) 
         endif
      enddo
   end subroutine set_zPrime

   function internal_to_rez(cell, vec) result(res)
      use m_types
      implicit none
      type(t_cell), intent(in) :: cell
      real, intent(in)      :: vec(3)
      real                  :: res(3)

      res = matmul(transpose(cell%bmat), vec)
   end function internal_to_rez

   subroutine set_kinED(mpi,   sphhar, atoms, sym, core_den, val_den, xcpot, &
                        input, noco,   stars, cell,     den,     EnergyDen, vTot)
      use m_types
      implicit none
      TYPE(t_mpi),INTENT(IN)       :: mpi
      TYPE(t_sphhar),INTENT(IN)    :: sphhar
      TYPE(t_atoms),INTENT(IN)     :: atoms
      TYPE(t_sym), INTENT(IN)      :: sym
      TYPE(t_potden),INTENT(IN)    :: core_den, val_den
      CLASS(t_xcpot),INTENT(INOUT) :: xcpot
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_potden),INTENT(IN)    :: den, EnergyDen, vTot

      call set_kinED_is(xcpot, input, noco, stars, sym, cell, den, EnergyDen, vTot)
      call set_kinED_mt(mpi,   sphhar,    atoms, sym, core_den, val_den, &
                           xcpot, EnergyDen, input, vTot)
   end subroutine set_kinED

   subroutine set_kinED_is(xcpot, input, noco, stars, sym, cell, den, EnergyDen, vTot)
      use m_types
      use m_pw_tofrom_grid
      implicit none
      CLASS(t_xcpot),INTENT(INOUT) :: xcpot
      TYPE(t_input),INTENT(IN)     :: input
      TYPE(t_noco),INTENT(IN)      :: noco
      TYPE(t_stars),INTENT(IN)     :: stars
      TYPE(t_sym), INTENT(IN)      :: sym
      TYPE(t_cell),INTENT(IN)      :: cell
      TYPE(t_potden),INTENT(IN)    :: den, EnergyDen, vTot

      !local arrays
      REAL, ALLOCATABLE            :: den_rs(:,:), ED_rs(:,:), vTot_rs(:,:)
      TYPE(t_gradients)            :: tmp_grad
      
      CALL init_pw_grid(xcpot,stars,sym,cell)

      CALL pw_to_grid(xcpot, input%jspins, noco%l_noco, stars, &
                      cell,  EnergyDen%pw, tmp_grad,    ED_rs)
      CALL pw_to_grid(xcpot, input%jspins, noco%l_noco, stars, &
                      cell,  vTot%pw,      tmp_grad,    vTot_rs)
      CALL pw_to_grid(xcpot, input%jspins, noco%l_noco, stars, &
                      cell,  den%pw,       tmp_grad,    den_rs)

      CALL finish_pw_grid()

      xcpot%kinED%is  = ED_RS - vTot_RS * den_RS
      xcpot%kinED%set = .True.
   end subroutine set_kinED_is

   subroutine set_kinED_mt(mpi,   sphhar,    atoms, sym, core_den, val_den, &
                           xcpot, EnergyDen, input, vTot)
      use m_types
      use m_mt_tofrom_grid
      implicit none
      TYPE(t_mpi),INTENT(IN)         :: mpi
      TYPE(t_sphhar),INTENT(IN)      :: sphhar
      TYPE(t_atoms),INTENT(IN)       :: atoms
      TYPE(t_sym), INTENT(IN)        :: sym
      TYPE(t_potden),INTENT(IN)      :: core_den, val_den, EnergyDen, vTot
      CLASS(t_xcpot),INTENT(INOUT)   :: xcpot
      TYPE(t_input),INTENT(IN)       :: input

      INTEGER                        :: jr, loc_n, n, n_start, n_stride, cnt
      REAL,ALLOCATABLE               :: vTot_mt(:,:,:), ED_rs(:,:), vTot_rs(:,:), vTot0_rs(:,:),&
                                        core_den_rs(:,:), val_den_rs(:,:)
      TYPE(t_gradients)              :: tmp_grad
      TYPE(t_sphhar)                 :: tmp_sphhar

#ifdef CPP_MPI
      n_start=mpi%irank+1
      n_stride=mpi%isize
#else
      n_start=1
      n_stride=1
#endif
      CALL init_mt_grid(input%jspins,atoms,sphhar,xcpot,sym)
      loc_n = 0
      allocate(ED_rs(atoms%nsp()*atoms%jmtd, input%jspins))
      allocate(vTot_rs, mold=ED_rs)
      allocate(vTot0_rs, mold=ED_rs)
      allocate(core_den_rs, mold=ED_rs)
      allocate(val_den_rs, mold=ED_rs)

      call xcpot%kinED%alloc_mt(atoms%nsp()*atoms%jmtd, input%jspins, &
                                n_start,                atoms%ntype,  n_stride)
      loc_n = 0
      do n = n_start,atoms%ntype,n_stride
         loc_n = loc_n + 1

         if(.not. allocated(vTot_mt)) then
            allocate(vTot_mt(lbound(vTot%mt, dim=1):ubound(vTot%mt, dim=1),&
                             lbound(vTot%mt, dim=2):ubound(vTot%mt, dim=2),&
                             lbound(vTot%mt, dim=4):ubound(vTot%mt, dim=4)))
         endif
         
         do jr=1,atoms%jri(n)
            vTot_mt(jr,0:,:) = vTot%mt(jr,0:,n,:) * atoms%rmsh(jr,n)**2
         enddo
         CALL mt_to_grid(xcpot, input%jspins, atoms, sphhar, EnergyDen%mt(:, 0:, n, :), &
                         n,     tmp_grad,     ED_rs)
         CALL mt_to_grid(xcpot, input%jspins, atoms, sphhar, vTot_mt(:,0:,:), &
                         n,     tmp_grad,     vTot_rs)
         
         tmp_sphhar%nlhd = sphhar%nlhd
         tmp_sphhar%nlh  = [(0, cnt=1,size(sphhar%nlh))]

         CALL mt_to_grid(xcpot, input%jspins, atoms, tmp_sphhar, vTot_mt(:,0:0,:), &
                         n,     tmp_grad,     vTot0_rs)
         CALL mt_to_grid(xcpot, input%jspins, atoms, sphhar, &
                         core_den%mt(:,0:,n,:), n, tmp_grad, core_den_rs)
         CALL mt_to_grid(xcpot, input%jspins, atoms, sphhar, &
                         val_den%mt(:,0:,n,:), n, tmp_grad, val_den_rs)
         
         xcpot%kinED%mt(:,:,loc_n) = ED_RS - (vTot0_rs * core_den_rs + vTot_rs * val_den_rs)
      enddo
      xcpot%kinED%set = .True.
      CALL finish_mt_grid()
   end subroutine set_kinED_mt
END MODULE m_metagga
