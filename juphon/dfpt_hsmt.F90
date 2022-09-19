!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_hsmt
   USE m_juDFT

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_hsmt(atoms, sym, enpara, iSpin, iDir, iDtype, input, fmpi, &
                      & noco, nococonv, cell, lapw, lapwq, usdus, td, tdV1, hmat, smat)
                      
      !> Setup of the MT part of the Hamiltonian and the overlap perturbation matrices
      !! Adapted from hsmt()
      !!
      !! There are two parts to this. For each atom, the part from the perturbed
      !! potential is calculated via
      !! 1. The non-spherical part in hsmt_nonsph()
      !! 2. The LO part in hsmt_lo() [with no smat passed]
      !!
      !! Additionally, ONLY for the perturbed atom, we need the unperturbed Hamiltonian
      !! and overlap with a prefactor of i(G'-G-q). This is done by first evaluating them
      !! and then passing the prefactor in a postprocess routine.
      !!
      !! The necessary noco logic is already implemented here similar to the base case
      !! in hsmt().
      !!
      !! DFPT-specific variables:
      !! - td, tdV1: Local matrix elements for the unperturbed Hamiltonian and
      !! the perturbed potential respectively.
      !! - lapwq: Set of LAPW basis vectors shifted by q.
      !! - iDir: Displacement direction.
      !! - iDtype: Type of the displaced atom.

      USE m_types
      USE m_types_mpimat
      USE m_hsmt_nonsph
      USE m_hsmt_sph
      USE m_hsmt_lo
      USE m_hsmt_distspins
      USE m_hsmt_fjgj
      USE m_hsmt_spinor
      USE m_hsmt_offdiag
      USE m_matrix_pref

      IMPLICIT NONE

      TYPE(t_mpi),      INTENT(IN)    :: fmpi
      TYPE(t_input),    INTENT(IN)    :: input
      TYPE(t_noco),     INTENT(IN)    :: noco
      TYPE(t_nococonv), INTENT(IN)    :: nococonv
      TYPE(t_sym),      INTENT(IN)    :: sym
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_lapw),     INTENT(IN)    :: lapw, lapwq
      TYPE(t_tlmplm),   INTENT(IN)    :: td, tdV1
      TYPE(t_usdus),    INTENT(IN)    :: usdus
      CLASS(t_mat),     INTENT(INOUT) :: smat(:,:),hmat(:,:)

      INTEGER, INTENT(IN) :: iSpin, iDir, iDtype

      TYPE(t_fjgj) :: fjgj

      INTEGER :: ilSpinPr, ilSpin
      INTEGER :: igSpinPr, igSpin, n
      COMPLEX :: chi(2,2),chi_one

      CLASS(t_mat), ALLOCATABLE :: smat_tmp, hmat_tmp, s1mat_tmp(:,:), h1mat_tmp(:,:)

      IF (noco%l_noco.AND..NOT.noco%l_ss) THEN
         IF (fmpi%n_size==1) THEN
            ALLOCATE(t_mat::hmat_tmp)
            ALLOCATE(t_mat::smat_tmp)
         ELSE
            ALLOCATE(t_mpimat::hmat_tmp)
            ALLOCATE(t_mpimat::smat_tmp)
         END IF
         CALL smat_tmp%init(hmat(1,1))
         CALL hmat_tmp%init(hmat(1,1))
         !$acc enter data copyin(smat_tmp,hmat_tmp)create(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
      END IF

      h1mat_tmp = hmat
      s1mat_tmp = smat

      DO ilSpinPr = MERGE(1, 1, noco%l_noco), MERGE(2, 1, noco%l_noco)
         DO ilSpin = MERGE(1, 1, noco%l_noco), MERGE(2, 1, noco%l_noco)
            CALL h1mat_tmp(ilSpinPr, ilSpin)%reset(CMPLX(0.0,0.0))
            CALL s1mat_tmp(ilSpinPr, ilSpin)%reset(CMPLX(0.0,0.0))
         END DO
      END DO

      CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,iSpin,noco)
      !$acc data copyin(fjgj) create(fjgj%fj,fjgj%gj)
      igSpinPr = 1; igSpin = 1; chi_one = 1.0 ! Defaults in non-noco case
      DO n = 1, atoms%ntype
         DO ilSpinPr = MERGE(1,iSpin,noco%l_noco), MERGE(2,iSpin,noco%l_noco)
            CALL timestart("fjgj coefficients")
            CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,n,ilSpinPr)
            !$acc update device(fjgj%fj,fjgj%gj)
            CALL timestop("fjgj coefficients")
            DO ilSpin = ilSpinPr, MERGE(2,iSpin,noco%l_noco)
               IF (.NOT.noco%l_noco) THEN
                  IF (n.EQ.iDtype) THEN
                     CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,s1mat_tmp(1,1),h1mat_tmp(1,1),.TRUE.,.TRUE.,lapwq)
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,h1mat_tmp(1,1),.FALSE.,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,h1mat_tmp(1,1),.FALSE.,.TRUE.,s1mat_tmp(1,1),lapwq)
                  END IF
                  CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat(1,1),.FALSE.,lapwq)
                  CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,.TRUE.,lapwq=lapwq)
                  !CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,smat(1,1))
               ELSE
                  IF (ilSpinPr==ilSpin) THEN !local spin-diagonal contribution
                     CALL hsmt_spinor(ilSpinPr,n,nococonv,chi)
                     IF (n.EQ.iDtype) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq)
                        CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat_tmp,hmat_tmp,.TRUE.,.TRUE.,lapwq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,smat_tmp,lapwq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                        CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                        CALL timestop("hsmt_distspins")
                     END IF
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,lapwq=lapwq)
                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat)
                     CALL timestop("hsmt_distspins")
                  ELSE IF (noco%l_unrestrictMT(n)) THEN
                     !2,1
                     CALL hsmt_spinor(3,n,nococonv,chi)
                     IF (n.EQ.iDtype) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,lapwq=lapwq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                        CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                        CALL timestop("hsmt_distspins")
                     END IF
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.FALSE.,.TRUE.,lapwq=lapwq)

                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat)
                     CALL timestop("hsmt_distspins")

                     !1,2
                     CALL hsmt_spinor(4,n,nococonv,chi)
                     IF (n.EQ.iDtype) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,lapwq=lapwq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                        CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                        CALL timestop("hsmt_distspins")
                     END IF
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.FALSE.,.TRUE.,lapwq=lapwq)

                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat)
                     CALL timestop("hsmt_distspins")
                  END IF
               END IF
            END DO
         END DO
      END DO
      !$acc end data
      ! TODO: Does this need some ACC magic?
      DO igSpinPr=MERGE(1,1,noco%l_noco),MERGE(2,1,noco%l_noco)
         DO igSpin=MERGE(1,1,noco%l_noco),MERGE(2,1,noco%l_noco)
            CALL matrix_pref(fmpi, cell%bmat, lapwq%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapwq%bkpt, lapw%bkpt, &
                           & lapwq%nv(igSpinPr), lapw%nv(igSpin), iDir, &
                           & h1mat_tmp(igSpinPr,igSpin), s1mat_tmp(igSpinPr,igSpin), hmat(igSpinPr,igSpin), smat(igSpinPr,igSpin))
         END DO
      END DO
      IF (noco%l_noco) THEN
         !$acc exit data delete(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
         !$acc exit data delete(smat_tmp,hmat_tmp)
      END IF
      RETURN
   END SUBROUTINE dfpt_hsmt
END MODULE m_dfpt_hsmt
