!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_add_vnonlocal
   USE m_judft
   USE m_types
   use m_types_mpimat
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
!     This module is the driver routine for the calculation of the Hartree    c
!     Fock exchange term by using the mixed basis set.                        c
!                                                                             c
!     hsfock                                                                  c
!         |                                                                   c
!         |- symm.F:                                                          c
!         |  calculates the irreducible representation                        c
!         |                                                                   c
!         |- wavefproducts.F:                 s      s*                       c
!         |  computes the repsentation of phi    phi       in the mixed basis c
!         |                                  n,k    n',k+q                    c
!         |                                                                   c
!         |- exchange.F:                                                      c
!         |  calculates valence-valence part of the exchange matrix (mat_ex), c
!         |                                                                   c
!         |- exchange_core.F                                                  c
!         |  calculate valence-core contribution                              c
!                                                                             c
!     variables:                                                              c
!         fi%kpts%nkptf   :=   number of kpoints                              c
!         fi%kpts%nkpt   :=   number of irreducible kpoints                   c
!         nbands  :=   number of bands for which the exchange matrix (mat_ex) c
!                      in the space of the wavefunctions is calculated        c
!         te_hfex :=   hf exchange contribution to the total energy           c
!         mnobd   :=   maximum number of occupied bands                       c
!         parent  :=   parent(ikpt) points to the symmetry equivalent point   c
!                      under the little group of kpoint nk                    c
!         symop   :=   symop(ikpt) points to the symmetry operation, which    c
!                      maps parent(ikpt) on ikpt                              c
!                                                                             c
!                                                                             c
!                                               M.Betzinger (09/07)           c
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
CONTAINS
   SUBROUTINE add_vnonlocal(nk, lapw, fi, hybdat, jsp,&
                            xcpot, fmpi, nococonv, hmat)
      USE m_constants
      USE m_symm_hf, ONLY: symm_hf
      USE m_intgrf, ONLY: intgrf, intgrf_init
      USE m_exchange_valence_hf
      USE m_exchange_core
      USE m_symmetrizeh
      USE m_wrapper
      USE m_hsefunctional, ONLY: exchange_vccvHSE, exchange_ccccHSE
      USE m_io_hybrid
      use m_glob_tofrom_loc
      IMPLICIT NONE

      type(t_fleurinput), intent(in) :: fi
      CLASS(t_xcpot), INTENT(IN)     :: xcpot
      TYPE(t_hybdat), INTENT(INOUT)  :: hybdat
      TYPE(t_lapw), INTENT(IN)       :: lapw
      type(t_mpi), intent(in)        :: fmpi
      type(t_nococonv),intent(in)    :: nococonv
      TYPE(t_mat), INTENT(INOUT)     :: hmat

      INTEGER, INTENT(IN)    :: jsp
      INTEGER, INTENT(IN)    :: nk

      ! local scalars
      INTEGER                   :: iband, nbasfcn, i, i0, j, ierr, iband_loc, pe_iband
      INTEGER                   :: tempI, tempJ
      integer, allocatable      :: list(:)
      REAL                      :: a_ex
      class(t_mat), allocatable :: tmp, z
      COMPLEX                   :: exch(fi%input%neig)

      call timestart("add_vnonlocal")
      call timestart("apply v_x")
      ! initialize weighting factor for HF exchange part
      a_ex = xcpot%get_exchange_weight()

      nbasfcn = lapw%nv(jsp) + fi%atoms%nlotot
      
      IF (hmat%l_real) THEN
         DO i = fmpi%n_rank+1,hmat%matsize1,fmpi%n_size
            i0=(i-1)/fmpi%n_size+1
            DO  j = 1,MIN(i,hmat%matsize1)
               hmat%data_r(j,i0) = hmat%data_r(j, i0) - a_ex * hybdat%v_x(nk, jsp)%data_r(j, i0)
            enddo
         enddo
      else
         DO i = fmpi%n_rank+1,hmat%matsize1,fmpi%n_size
            i0=(i-1)/fmpi%n_size+1
            DO  j = 1,MIN(i,hmat%matsize1)
               hmat%data_c(j,i0) = hmat%data_c(j, i0) - a_ex * CONJG(hybdat%v_x(nk, jsp)%data_c(j, i0))
            enddo
         enddo
      endif
      call timestop("apply v_x")

      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::z, tmp)
      ELSE
         ALLOCATE (t_mpimat::z, tmp)
      END IF

      CALL z%init(hmat%l_real, nbasfcn, hybdat%nbands(nk, jsp), fmpi%sub_comm, .false.)
      list = [(i, i= fmpi%n_rank+1,hybdat%nbands(nk,jsp), fmpi%n_size )]
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv,  fi%input, nk, jsp, z, list=list)

#ifdef CPP_MPI
      call timestart("post add_vnonl read_z barrier")
      call MPI_Barrier(fmpi%mpi_comm, ierr)
      call timestop("post add_vnonl read_z barrier")
#endif

      ! calculate exchange contribution of current k-point nk to total energy (te_hfex)
      ! in the case of a spin-unpolarized calculation the factor 2 is added in eigen.F90

      exch = 0
      select type(vx =>hybdat%v_x(nk, jsp))
      class is (t_mat)
         if(nbasfcn /= vx%matsize2) call juDFT_error("these dimension should match. is this a spin issue?")
      class is (t_mpimat)
         if(nbasfcn /= vx%global_size2) call juDFT_error("these dimension should match. is this a spin issue?")
      end select
      !z%matsize1 = MIN(z%matsize1, hybdat%v_x(nk, jsp)%matsize2)
      call tmp%init(hmat%l_real, nbasfcn, hybdat%nbands(nk, jsp), fmpi%sub_comm, .false.)
      IF (hybdat%v_x(nk, jsp)%l_real) then
         CALL hybdat%v_x(nk, jsp)%multiply(z, tmp)
      else
         CALL hybdat%v_x(nk, jsp)%multiply(z, tmp, transA="T")
      endif
      ! WRITE (oUnit, '(A)') "          K-points,   iband,    exch - div (eV), div (eV),  exch (eV)"
      DO iband = 1, hybdat%nbands(nk, jsp)
         call glob_to_loc(fmpi, iband, pe_iband, iband_loc)
         if(pe_iband == fmpi%n_rank) then
            IF (z%l_real) THEN
               exch(iband) = dot_product(z%data_r(:z%matsize1, iband_loc), tmp%data_r(:, iband_loc))
            ELSE
               exch(iband) = dot_product(z%data_c(:z%matsize1, iband_loc), tmp%data_c(:, iband_loc))
            END IF
         endif
#ifdef CPP_MPI
         call MPI_Bcast(exch(iband), 1, MPI_DOUBLE_COMPLEX, pe_iband, fmpi%sub_comm, ierr)
#endif
         IF (iband <= hybdat%nobd(nk,jsp) .AND. allocated(hybdat%results%w_iks)) THEN ! Issue #769
            hybdat%results%te_hfex%valence = hybdat%results%te_hfex%valence - real(a_ex*hybdat%results%w_iks(iband, nk, jsp)*exch(iband))
         END IF
         ! IF (hybdat%l_calhf) THEN
         !    WRITE (oUnit, '(      ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,4X,3F15.5)') &
         !       fi%kpts%bkf(:, nk), iband, (REAL(exch(iband)) - hybdat%div_vv(iband, nk, jsp))*(-hartree_to_ev_const), &
         !       hybdat%div_vv(iband, nk, jsp)*(-hartree_to_ev_const), REAL(exch(iband))*(-hartree_to_ev_const)
         ! END IF
      END DO
      call timestop("add_vnonlocal")
   END SUBROUTINE add_vnonlocal
END MODULE m_add_vnonlocal
