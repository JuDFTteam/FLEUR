!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
!>  The two Bloch-basis inputs Wannier90 needs from one spin channel: the projections
!>  A_mn(k) and the neighbour overlaps M_mn(k,b).
!>
!>  This is where the matrix-element layer is consumed to wannierise, and the reason it is
!>  its own module: the pass reads states and their augmentation coefficients for every k
!>  this rank owns, and none of that belongs in the routine that orchestrates the run.
!>
!>  The anchor is taken per k here and released by the CALLER, once, after every channel is
!>  done -- see matrix_element_release_anchor in m_wannierlib_main. It is a single module
!>  slot that each call re-points, so the release is not per k; keeping it with the caller
!>  also keeps it valid across the pair-overlap passes that run later.
MODULE m_wannierlib_build_amn_mmn
   USE m_juDFT
   USE m_types, ONLY: t_stars
   USE m_types_atoms
   USE m_types_cell
   USE m_types_vacuum
   USE m_types_input
   USE m_types_kpts
   USE m_types_lapw
   USE m_types_noco
   USE m_types_nococonv
   USE m_types_sym
   USE m_types_enpara
   USE m_types_mpi
   USE m_types_potden
   USE m_types_usdus
   USE m_types_mat
   USE m_types_radfun
   USE m_types_abc
   USE m_types_wannierlib
   USE m_types_wgauge_manifold, ONLY: t_wgauge_manifold
   USE m_types_wgauge_bmesh, ONLY: t_wgauge_bmesh
   USE m_types_spinor_layout, ONLY: radial_slot
   USE m_matrix_element_factory, ONLY: matrix_element_states
   USE m_melem_ujugaunt
   USE m_wannierlib_amn
   USE m_wannierlib_mmnkb
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: wannierlib_build_amn_mmn, &
             wannierlib_reduce_amn, wannierlib_gather_mmn
CONTAINS

   SUBROUTINE wannierlib_build_amn_mmn(this, manifold, bmesh, atoms, cell, input, kpts, sym, &
                                       noco, nococonv, stars, enpara, fmpi, vtot, eig_id, &
                                       radfun, usdus, distk, kdiff, nntot_w90, jspin, &
                                       l_wannierlib_spinors, amn, mmn, vacuum)
      TYPE(t_wannierlib_wannierize), INTENT(IN) :: this
      TYPE(t_wgauge_manifold), INTENT(IN) :: manifold
      TYPE(t_wgauge_bmesh),    INTENT(IN) :: bmesh
      TYPE(t_atoms),    INTENT(IN) :: atoms
      TYPE(t_cell),     INTENT(IN) :: cell
      TYPE(t_vacuum),   INTENT(IN) :: vacuum
      TYPE(t_input),    INTENT(IN) :: input
      TYPE(t_kpts),     INTENT(IN) :: kpts
      TYPE(t_sym),      INTENT(IN) :: sym
      TYPE(t_noco),     INTENT(IN) :: noco
      TYPE(t_nococonv), INTENT(IN) :: nococonv
      TYPE(t_stars),    INTENT(IN) :: stars
      TYPE(t_enpara),   INTENT(IN) :: enpara
      TYPE(t_mpi),      INTENT(IN) :: fmpi
      TYPE(t_potden),   INTENT(IN) :: vtot
      INTEGER,          INTENT(IN) :: eig_id
      TYPE(t_radfun),   INTENT(IN) :: radfun(:)
      TYPE(t_usdus),    INTENT(IN) :: usdus
      INTEGER, INTENT(IN) :: distk(:)          !> (nkptf) owning rank of each k
      REAL,    INTENT(IN) :: kdiff(:, :)
      INTEGER, INTENT(IN) :: nntot_w90
      INTEGER, INTENT(IN) :: jspin             !> the channel being wannierised
      LOGICAL, INTENT(IN) :: l_wannierlib_spinors
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: amn(:, :, :)     !> (num_bands, num_wann, nkptf)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: mmn(:, :, :, :)  !> (nb, nb, nntot, nk_loc)

      INTEGER :: ikpt, ib, ierr, irec, ik_local, nk_local
      INTEGER :: jspin_comp, jspin_rad, jspin_rad_done
      INTEGER, ALLOCATABLE :: ev_list(:)
      COMPLEX, ALLOCATABLE :: ujug(:, :, :, :, :, :)
      TYPE(t_lapw) :: lapw
      TYPE(t_mat),    POINTER :: zmat_p(:)    !> into the factory cache
      TYPE(t_abc),    POINTER :: abc_p(:, :)  !> likewise

      ! calculate the  matrices for all k-points
      ALLOCATE (amn(this%num_bands, this%num_wann, kpts%nkptf), stat=ierr, source=cmplx(0.0, 0.0))
      IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating amn buffer', calledby='wannierlib_build_amn_mmn')

      nk_local = COUNT(distk == fmpi%irank)
      ALLOCATE(mmn(this%num_bands, this%num_bands, nntot_w90, nk_local), stat=ierr, source=cmplx(0.0, 0.0))
      IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating local mmn buffer', calledby='wannierlib_build_amn_mmn')

      jspin_rad_done = -1
      DO jspin_comp = MERGE(1, jspin, l_wannierlib_spinors), MERGE(2, jspin, l_wannierlib_spinors)
         ! jspin_comp = the eig record (spinor up/down); jspin_rad = the radial index.
         jspin_rad = radial_slot(radfun, jspin_comp)
         !> These depend on the radial set and on nothing else in this loop. With a single
         !> potential both spinor components read the same set, so the second pass would
         !> rebuild an identical array -- and it is the largest thing built here, lm by lm
         !> by radial pair by type by neighbour.
         IF (jspin_rad /= jspin_rad_done) THEN
            CALL melem_ujugaunt(atoms, cell, nntot_w90, kdiff, radfun, radfun, jspin_rad, jspin_rad, .FALSE., 1, ujug)
            jspin_rad_done = jspin_rad
         END IF

         ev_list = [(ib, ib = this%min_band, this%max_band)]
         ik_local = 0
         DO ikpt = 1, kpts%nkptf
            IF (distk(ikpt) /= fmpi%irank) CYCLE   ! each rank computes only its k-slice -> parallel eigenvector I/O
            !> This k stays put while its neighbours come and go: mmnkb asks the factory
            !> for one per b, and without the anchor this one would be the oldest by the
            !> third of them and be overwritten under the pointers held here.
            CALL lapw%init(input, noco, nococonv, kpts, atoms, sym, ikpt, cell)
            CALL matrix_element_states(eig_id, ikpt, input, atoms, sym, cell, noco, &
                                       nococonv, enpara, lapw, vtot, fmpi, zmat_p, abc_p, &
                                       ev_list=ev_list, &
                                       l_both_spinors=(noco%l_soc .AND. .NOT. noco%l_noco), &
                                       kpts=kpts, l_anchor=.TRUE.)
            !> Non-collinearly the whole spinor is one record; otherwise each channel is
            !> its own and the spin block is reached by row offset further down.
            irec = MERGE(1, jspin_comp, noco%l_noco)

            CALL wannierlib_amn(this, atoms, kpts, ikpt, usdus, radfun, abc_p(jspin_comp, :), l_wannierlib_spinors, jspin_comp, jspin_rad, amn(:, :, ikpt))

            ik_local = ik_local + 1
            CALL wannierlib_mmnkb(manifold, bmesh, ikpt, kpts, &
                                  ujug, atoms, cell, input, sym, noco, nococonv, &
                                  abc_p(jspin_comp, :), jspin_comp, jspin_rad, eig_id, stars, lapw, &
                                  zmat_p(irec), mmn, ik_local, enpara, vtot, fmpi, vacuum, radfun)
         END DO

      END DO
      IF (ALLOCATED(ujug)) DEALLOCATE (ujug)
   END SUBROUTINE wannierlib_build_amn_mmn


   ! Complete the distributed amn: each rank filled only its distk k-slice (zeros elsewhere),
   ! so an MPI_ALLREDUCE(SUM) reassembles the full amn on every rank (needed by w90_set_u_opt).
   SUBROUTINE wannierlib_reduce_amn(fmpi, amn)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      COMPLEX, INTENT(INOUT) :: amn(:, :, :)
#ifdef CPP_MPI
      INTEGER :: ierr
      IF (fmpi%isize > 1) &
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, amn, SIZE(amn), MPI_DOUBLE_COMPLEX, MPI_SUM, fmpi%mpi_comm, ierr)
#endif
   END SUBROUTINE wannierlib_reduce_amn

   ! Reassemble the full-mesh mmn so it can be written as one file. amn is allocated over
   ! every k and left at zero outside each rank's slice, so summing completes it; mmn is
   ! stored compactly as (nb, nb, nntot, nk_local) indexed by POSITION within the slice, so
   ! the slices have to be placed at their global k before being summed. The local order is
   ! ascending global k (wannierlib_build_amn_mmn fills it in that order), and every k is
   ! owned by exactly one rank, so the sum is a copy. This is the one full-mesh buffer the
   ! distributed post-processing otherwise avoids, so it is allocated only when asked for.
   SUBROUTINE wannierlib_gather_mmn(fmpi, distk, nkptf, mmn_loc, mmn_full)
#ifdef CPP_MPI
      use mpi
#endif
      TYPE(t_mpi), INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: distk(:)             ! (nkptf) owning rank of each k
      INTEGER, INTENT(IN) :: nkptf
      COMPLEX, INTENT(IN) :: mmn_loc(:, :, :, :)  ! (nb, nb, nntot, nk_local)
      COMPLEX, ALLOCATABLE, INTENT(OUT) :: mmn_full(:, :, :, :)
      INTEGER :: ikpt, ik_local, ierr
#ifdef CPP_MPI
      INTEGER :: mpi_err
#endif
      ALLOCATE (mmn_full(SIZE(mmn_loc, 1), SIZE(mmn_loc, 2), SIZE(mmn_loc, 3), nkptf), &
                stat=ierr, source=CMPLX(0.0, 0.0))
      IF (ierr /= 0) CALL juDFT_error('wannierlib failed allocating the full mmn buffer', &
                                      calledby='wannierlib_gather_mmn')
      ik_local = 0
      DO ikpt = 1, nkptf
         IF (distk(ikpt) /= fmpi%irank) CYCLE
         ik_local = ik_local + 1
         mmn_full(:, :, :, ikpt) = mmn_loc(:, :, :, ik_local)
      END DO
#ifdef CPP_MPI
      IF (fmpi%isize > 1) &
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, mmn_full, SIZE(mmn_full), MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, fmpi%mpi_comm, mpi_err)
#endif
   END SUBROUTINE wannierlib_gather_mmn
END MODULE m_wannierlib_build_amn_mmn
