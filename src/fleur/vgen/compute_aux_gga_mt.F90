!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!> Compute the auxiliary GGA XC spherical potential for MetaGGA calculations.
!!
!! When a MetaGGA functional is used, the radial basis functions (u_l, u̇_l)
!! must be generated from a GGA potential rather than the non-multiplicative
!! MetaGGA potential. This routine computes the spherical (l=0) component
!! of the auxiliary GGA XC potential for each atom type.
!!
!! Reference: Doumont et al., Phys. Rev. B 106, 235159 (2022), Section II.B
MODULE m_compute_aux_gga_mt

   USE m_juDFT
   IMPLICIT NONE

CONTAINS

   SUBROUTINE compute_aux_gga_mt(xcpot, atoms, sphhar, sym, input, noco, fmpi, den, auxGGA_vxc_sph)
      USE m_mt_tofrom_grid
      USE m_libxc_postprocess_gga
      USE m_types
#ifdef CPP_MPI
      USE mpi
#endif

      IMPLICIT NONE

      CLASS(t_xcpot), INTENT(IN)    :: xcpot
      TYPE(t_atoms), INTENT(IN)     :: atoms
      TYPE(t_sphhar), INTENT(IN)    :: sphhar
      TYPE(t_sym), INTENT(IN)       :: sym
      TYPE(t_input), INTENT(IN)     :: input
      TYPE(t_noco), INTENT(IN)      :: noco
      TYPE(t_mpi), INTENT(IN)       :: fmpi
      TYPE(t_potden), INTENT(IN)    :: den
      REAL, INTENT(OUT)             :: auxGGA_vxc_sph(:,:,:)  ! (jmtd, ntype, jspins)

      ! Local variables
      TYPE(t_gradients) :: grad
      REAL, ALLOCATABLE :: ch(:,:), v_xc(:,:), v_x(:,:)
      REAL, ALLOCATABLE :: aux_mt(:,:,:)
      INTEGER :: n, nsp, n_start, n_stride
#ifdef CPP_MPI
      INTEGER :: ierr
#endif

      CALL timestart("compute_aux_gga_mt")

      auxGGA_vxc_sph = 0.0

      ! Initialize MT grid with gradient support (needed for GGA)
      CALL init_mt_grid(input%jspins, atoms, sphhar, .TRUE., sym)

      nsp = atoms%nsp()

#ifdef CPP_MPI
      n_start = fmpi%irank + 1
      n_stride = fmpi%isize
#else
      n_start = 1
      n_stride = 1
#endif

      DO n = n_start, atoms%ntype, n_stride
         ALLOCATE(ch(nsp*atoms%jri(n), input%jspins))
         ALLOCATE(v_xc(nsp*atoms%jri(n), input%jspins))
         ALLOCATE(v_x(nsp*atoms%jri(n), input%jspins))
         ALLOCATE(aux_mt(atoms%jri(n), 0:sphhar%nlhd, input%jspins))

         ! Allocate gradient arrays
         CALL xcpot%alloc_gradients(SIZE(ch, 1), input%jspins, grad)

         ! Convert density to real-space grid and compute gradients
         CALL mt_to_grid(.TRUE., input%jspins, atoms, sym, sphhar, .TRUE., &
                         den%mt(:, 0:, n, :), n, noco, grad, ch)

         ! Apply density cutoffs (same as in main XC computation)
         CALL xcpot%apply_cutoffs(1.E-6, ch, grad)

         ! Evaluate auxiliary GGA functional
         CALL xcpot%get_vxc(input%jspins, ch, v_xc, v_x, grad, l_aux=.TRUE.)

         ! Apply GGA postprocessing: v_xc -= 2*div(vsigma * grad(rho))
         CALL libxc_postprocess_gga_mt(xcpot, atoms, sym, sphhar, noco, n, v_xc, grad, atom_num=n)

         ! Project auxiliary v_xc back to spherical harmonics representation
         aux_mt = 0.0
         CALL mt_from_grid(atoms, sym, sphhar, n, input%jspins, v_xc, aux_mt)

         ! Extract the l=0 (spherical) component
         auxGGA_vxc_sph(1:atoms%jri(n), n, :) = aux_mt(1:atoms%jri(n), 0, :)

         DEALLOCATE(ch, v_xc, v_x, aux_mt)
      ENDDO

      CALL finish_mt_grid()

#ifdef CPP_MPI
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, auxGGA_vxc_sph, SIZE(auxGGA_vxc_sph), &
                         MPI_DOUBLE_PRECISION, MPI_SUM, fmpi%mpi_comm, ierr)
#endif

      CALL timestop("compute_aux_gga_mt")

   END SUBROUTINE compute_aux_gga_mt

END MODULE m_compute_aux_gga_mt
