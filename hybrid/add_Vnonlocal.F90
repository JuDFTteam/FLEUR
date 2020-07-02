!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_add_vnonlocal
   USE m_judft
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
!         fi%kpts%nkptf   :=   number of kpoints                                      c
!         fi%kpts%nkpt   :=   number of irreducible kpoints                          c
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
   SUBROUTINE add_vnonlocal(nk, lapw, fi, hybdat, jsp, results,&
                            xcpot, fmpi, nococonv, hmat)

      USE m_types
      USE m_constants
      USE m_symm_hf, ONLY: symm_hf
      USE m_intgrf, ONLY: intgrf, intgrf_init
      USE m_exchange_valence_hf
      USE m_exchange_core
      USE m_symmetrizeh
      USE m_wrapper
      USE m_hsefunctional, ONLY: exchange_vccvHSE, exchange_ccccHSE
      USE m_io_hybinp
      use m_judft

      IMPLICIT NONE

      type(t_fleurinput), intent(in) :: fi
      TYPE(t_results), INTENT(INOUT) :: results
      CLASS(t_xcpot), INTENT(IN)     :: xcpot
      TYPE(t_hybdat), INTENT(INOUT)  :: hybdat
      TYPE(t_lapw), INTENT(IN)       :: lapw
      type(t_mpi), intent(in)        :: fmpi
      type(t_nococonv),intent(in)    :: nococonv
      TYPE(t_mat), INTENT(INOUT)     :: hmat

      INTEGER, INTENT(IN)    :: jsp
      INTEGER, INTENT(IN)    :: nk

      ! local scalars
      INTEGER                 :: n, nn, iband, nbasfcn, i, i0, j
      REAL                    :: a_ex
      TYPE(t_mat)             :: tmp, v_x, z
      COMPLEX                 :: exch(fi%input%neig, fi%input%neig)

      call timestart("add_vnonlocal")

      ! initialize weighting factor for HF exchange part
      a_ex = xcpot%get_exchange_weight()

      nbasfcn = MERGE(lapw%nv(1) + lapw%nv(2) + 2*fi%atoms%nlotot, lapw%nv(1) + fi%atoms%nlotot, fi%noco%l_noco)
      CALL v_x%init(hmat%l_real, nbasfcn, nbasfcn)
      CALL read_v_x(v_x, fi%kpts%nkpt*(jsp - 1) + nk)

      
      ! add non-local x-potential to the hamiltonian hmat

      ! write (*,*) "shape(hmat%data_r)", shape(hmat%data_r)
      ! write (*,*) "shape(v_x%data_r)", shape(v_x%data_r)
      ! if(any( shape(hmat%data_r) /= shape(v_x%data_r) ) ) then 
      !    if(all(shape(v_x%data_r) /= 0)) then
      !       call judft_error("shapes don't agree")
      !    endif
      ! endif

      if(v_x%matsize1 > 0) then 
         call v_x%u2l()
      endif
      
      IF (hmat%l_real) THEN
         DO i = fmpi%n_rank+1,v_x%matsize1,fmpi%n_size
            i0=(i-1)/fmpi%n_size+1
            DO  j = 1,MIN(i,v_x%matsize1) 
               hmat%data_r(j,i0) = hmat%data_r(j, i0) - a_ex * v_x%data_r(j, i)
            enddo
         enddo
      else         
         DO i = fmpi%n_rank+1,v_x%matsize1,fmpi%n_size
            i0=(i-1)/fmpi%n_size+1
            DO  j = 1,MIN(i,v_x%matsize1) 
               hmat%data_c(j,i0) = hmat%data_c(j, i0) - a_ex * v_x%data_c(j, i)
            enddo
         enddo
      endif

      CALL z%init(hmat%l_real, nbasfcn, fi%input%neig)

      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv,  fi%input, nk, jsp, z)

      ! calculate exchange contribution of current k-point nk to total energy (te_hfex)
      ! in the case of a spin-unpolarized calculation the factor 2 is added in eigen.F90
      IF (.NOT. v_x%l_real) v_x%data_c = conjg(v_x%data_c)
      exch = 0
      z%matsize1 = MIN(z%matsize1, v_x%matsize2)

      CALL v_x%multiply(z, tmp)

      DO iband = 1, hybdat%nbands(nk)
         IF (z%l_real) THEN
            exch(iband, iband) = dot_product(z%data_r(:z%matsize1, iband), tmp%data_r(:, iband))
         ELSE
            exch(iband, iband) = dot_product(z%data_c(:z%matsize1, iband), tmp%data_c(:, iband))
         END IF
         IF (iband <= hybdat%nobd(nk,jsp)) THEN
            results%te_hfex%valence = results%te_hfex%valence - a_ex*results%w_iks(iband, nk, jsp)*exch(iband, iband)
         END IF
         IF (hybdat%l_calhf) THEN
            WRITE (oUnit, '(      ''  ('',F5.3,'','',F5.3,'','',F5.3,'')'',I4,4X,3F15.5)') &
               fi%kpts%bkf(:, nk), iband, (REAL(exch(iband, iband)) - hybdat%div_vv(iband, nk, jsp))*(-hartree_to_ev_const), &
               hybdat%div_vv(iband, nk, jsp)*(-hartree_to_ev_const), REAL(exch(iband, iband))*(-hartree_to_ev_const)
         END IF
      END DO
      call timestop("add_vnonlocal")
   END SUBROUTINE add_vnonlocal

END MODULE m_add_vnonlocal
