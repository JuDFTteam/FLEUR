!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_hsfock

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
!         parent  :=   parent(ikpt) points to the symmetry equivalent point   c
!                      under the little group of kpoint nk                    c
!         symop   :=   symop(ikpt) points to the symmetry operation, which    c
!                      maps parent(ikpt) on ikpt                              c
!                                                                             c
!                                                                             c
!                                               M.Betzinger (09/07)           c
! c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c

CONTAINS

   SUBROUTINE hsfock(fi, k_pack, mpdata, lapw, jsp, hybdat, &
                     eig_irr, nococonv, stars, &
                     results, xcpot, fmpi)

      use m_ex_to_vx
      USE m_judft
      USE m_types
      USE m_intgrf
      USE m_wrapper
      USE m_io_hybinp
      USE m_hsefunctional
      USE m_symm_hf
      USE m_exchange_valence_hf
      USE m_exchange_core
      USE m_symmetrizeh
      use m_work_package
      USE m_eig66_data
      use m_eig66_mpi
      IMPLICIT NONE

      type(t_fleurinput), intent(in)    :: fi
      type(t_k_package), intent(in)     :: k_pack
      TYPE(t_xcpot_inbuild), INTENT(IN)    :: xcpot
      TYPE(t_mpi), INTENT(IN)    :: fmpi
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_lapw), INTENT(IN)    :: lapw
      type(t_stars), intent(in)   :: stars
      TYPE(t_mpdata), intent(inout)  :: mpdata
      TYPE(t_hybdat), INTENT(INOUT) :: hybdat
      TYPE(t_results), INTENT(INOUT) :: results

      ! scalars
      INTEGER, INTENT(IN)    :: jsp

      ! arrays
      REAL, INTENT(IN)    :: eig_irr(:, :)

      ! local scalars
      INTEGER                 ::  l, itype
      INTEGER                 ::  iband, nk
      INTEGER                 ::  ikpt, ikpt0
      INTEGER                 ::  nbasfcn
      INTEGER                 ::  nsymop
      INTEGER                 ::  ncstd
      INTEGER                 ::  ok
      REAL                    ::  a_ex

      ! local arrays
      INTEGER                 ::  nsest(hybdat%nbands(k_pack%nk )), indx_sest(hybdat%nbands(k_pack%nk ), hybdat%nbands(k_pack%nk ))
      INTEGER                 ::  rrot(3, 3, fi%sym%nsym)
      INTEGER                 ::  psym(fi%sym%nsym) ! Note: psym is only filled up to index nsymop

      INTEGER, ALLOCATABLE     ::  parent(:)
      INTEGER, ALLOCATABLE     ::  n_q(:)

      complex                  :: c_phase_k(hybdat%nbands(k_pack%nk ))
      REAL                     ::  wl_iks(fi%input%neig, fi%kpts%nkptf)
      TYPE(t_mat)              :: ex, z_k
      TYPE(t_data_MPI), POINTER, ASYNCHRONOUS :: d

      CALL timestart("total time hsfock")
      nk = k_pack%nk 
      ! initialize weighting factor for HF exchange part
      a_ex = xcpot%get_exchange_weight()
      ncstd = sum([((hybdat%nindxc(l, itype)*(2*l + 1)*fi%atoms%neq(itype), l=0, hybdat%lmaxc(itype)), itype=1, fi%atoms%ntype)])
      IF(nk == 1 .and. fmpi%irank == 0) WRITE(*, *) 'calculate new HF matrix'
      IF(nk == 1 .and. jsp == 1 .and. fi%input%imix > 10) CALL system('rm -f broyd*')
      ! calculate all symmetrie operations, which yield k invariant

      allocate(parent(fi%kpts%nkptf), stat=ok)
      IF(ok /= 0) call judft_error('mhsfock: failure allocation parent')
      parent = 0

      nbasfcn = lapw%hyb_num_bas_fun(fi)
      call z_k%init(fi%sym%invs, nbasfcn, fi%input%neig)
      call read_z(fi%atoms, fi%cell, hybdat, fi%kpts, fi%sym, fi%noco, nococonv,  fi%input, nk, jsp, z_k, &
                   c_phase=c_phase_k)    
#ifdef(CPP_MPI) 
      call timestart("Post read_z Barrier: hsfock")
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call timestop("Post read_z Barrier: hsfock")
#endif

      CALL symm_hf_init(fi, nk, nsymop, rrot, psym)

      CALL symm_hf(fi, nk, hybdat, eig_irr, mpdata, lapw, nococonv, z_k, c_phase_k, jsp, &
                   rrot, nsymop, psym, n_q, parent, nsest, indx_sest)

      ! remove weights(wtkpt) in w_iks
      DO ikpt = 1, fi%kpts%nkptf
         DO iband = 1, fi%input%neig
            ikpt0 = fi%kpts%bkp(ikpt)
            wl_iks(iband, ikpt) = results%w_iks(iband, ikpt0, jsp)/(fi%kpts%wtkpt(ikpt0)*fi%kpts%nkptf)
         END DO
      END DO

      ! calculate contribution from valence electrons to the
      ! HF exchange
      ex%l_real = fi%sym%invs
      CALL exchange_valence_hf(k_pack, fi, z_k, c_phase_k, mpdata, jsp, hybdat, lapw, eig_irr, results, &
                               n_q, wl_iks, xcpot, nococonv, stars, nsest, indx_sest, fmpi, ex)

      if(.not. allocated(hybdat%v_x)) allocate(hybdat%v_x(fi%kpts%nkpt, fi%input%jspins))

      if(k_pack%submpi%root()) then
         ! calculate contribution from the core states to the HF exchange
         CALL timestart("core exchange calculation")
         IF(xcpot%is_name("hse") .OR. xcpot%is_name("vhse")) THEN
            call judft_error('HSE not implemented in hsfock')
         ELSE
            CALL exchange_vccv1(nk, fi%input, fi%atoms, fi%cell, fi%kpts, fi%sym, fi%noco, nococonv, fi%oneD, &
                              mpdata, fi%hybinp, hybdat, jsp, &
                              lapw, nsymop, nsest, indx_sest, fmpi, a_ex, results, ex)

            CALL exchange_cccc(nk, fi%atoms, hybdat, ncstd, fi%sym, fi%kpts, a_ex, results)
         END IF

         CALL timestop("core exchange calculation")

         call ex_to_vx(fi, nk, jsp, nsymop, psym, hybdat, lapw, z_k, ex, hybdat%v_x(nk, jsp))
         call hybdat%v_x(nk, jsp)%u2l()
      endif
      call hybdat%v_x(nk,jsp)%bcast(0, k_pack%submpi%comm)
      hybdat%l_addhf = .True.
      CALL timestop("total time hsfock")
   END SUBROUTINE hsfock
END MODULE m_hsfock
