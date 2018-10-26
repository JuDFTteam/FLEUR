!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_metagga
   PUBLIC  :: calc_EnergyDen
   PRIVATE :: calc_EnergyDen_auxillary_weights

   type t_RS_potden
      REAL, ALLOCATABLE :: is(:,:), mt(:,:)
   end type t_RS_potden

CONTAINS
   SUBROUTINE calc_kinEnergyDen(EnergyDen_rs, vTot_rs, den_rs, kinEnergyDen_RS)
#ifdef CPP_LIBXC
      IMPLICIT NONE
      REAL, INTENT(in)                 :: den_RS(:,:), EnergyDen_RS(:,:), vTot_RS(:,:)
      REAL, INTENT(inout), allocatable :: kinEnergyDen_RS(:,:)
      REAL, PARAMETER                  :: eps = 1e-15

      !implicit allocation
      kinEnergyDen_RS = EnergyDen_RS - vTot_RS * den_RS

      if(any(kinEnergyDen_RS < eps)) then
         write (6,*) "         lowest kinetic energy density cutoff = ", minval(kinEnergyDen_RS)
         kinEnergyDen_RS = max(kinEnergyDen_RS, eps)
      endif
   
#else
      USE m_juDFT_stop
      CALL juDFT_error("MetaGGA require LibXC",hint="compile Fleur with LibXC (e.g. by giving '-external libxc' to ./configure")
#endif
   END SUBROUTINE calc_kinEnergyDen

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
      REAL                       :: w_i(SIZE(f_ik,dim=1)), e_i(SIZE(f_ik,dim=1))
      INTEGER                    :: ikpt

      DO ikpt = 1,kpts%nkpt
         CALL read_eig(eig_id,ikpt,jspin, eig=e_i, w_iks=w_i)

         f_ik(:,ikpt) = e_i * w_i
      ENDDO
   END SUBROUTINE calc_EnergyDen_auxillary_weights
   
END MODULE m_metagga
