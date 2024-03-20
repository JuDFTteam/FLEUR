!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_hsmt
   USE m_juDFT

   IMPLICIT NONE

CONTAINS
   SUBROUTINE dfpt_hsmt(atoms, sym, juphon, enpara, iSpin, iDir, iDtype, input, fmpi, &
                      & noco, nococonv, cell, lapw, lapwq, usdus, td, tdV1, hmat, smat, nk, killcont)

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
      TYPE(t_juphon),   INTENT(IN)    :: juphon
      TYPE(t_cell),     INTENT(IN)    :: cell
      TYPE(t_atoms),    INTENT(IN)    :: atoms
      TYPE(t_enpara),   INTENT(IN)    :: enpara
      TYPE(t_lapw),     INTENT(IN)    :: lapw, lapwq
      TYPE(t_tlmplm),   INTENT(IN)    :: td, tdV1
      TYPE(t_usdus),    INTENT(IN)    :: usdus
      CLASS(t_mat),     INTENT(INOUT) :: smat(:,:),hmat(:,:)

      INTEGER, INTENT(IN) :: iSpin, iDir, iDtype, nk, killcont(3)

      TYPE(t_fjgj) :: fjgj, fjgjq

      INTEGER :: ilSpinPr, ilSpin, nspins, i, j
      INTEGER :: igSpinPr, igSpin, n
      COMPLEX :: chi(2,2),chi_one

      CLASS(t_mat), ALLOCATABLE :: smat_tmp, hmat_tmp, s1mat_tmp(:,:), h1mat_tmp(:,:)

      !TODO: All of the openACC is most certainly scuffed for DFPT. Fix it someday.
      !      But wait until it is right and proper in the main code!
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

      nspins = MERGE(2, 1, noco%l_noco)
      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::s1mat_tmp(nspins, nspins), h1mat_tmp(nspins, nspins))
      ELSE
         ALLOCATE (t_mpimat::s1mat_tmp(nspins, nspins), h1mat_tmp(nspins, nspins))
      END IF

      DO i = 1, nspins
         DO j = 1, nspins
            CALL s1mat_tmp(i, j)%init(.FALSE., lapwq%nv(i) + atoms%nlotot, lapw%nv(j) + atoms%nlotot, fmpi%sub_comm, .false.)
            CALL h1mat_tmp(i, j)%init(s1mat_tmp(i, j))
         END DO
      END DO

      CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,iSpin,noco)
      CALL fjgjq%alloc(MAXVAL(lapwq%nv),atoms%lmaxd,iSpin,noco)
      !!$acc data copyin(fjgj) create(fjgj%fj,fjgj%gj)
      !!$acc data copyin(fjgjq) create(fjgjq%fj,fjgjq%gj)
      igSpinPr = 1; igSpin = 1; chi_one = 1.0 ! Defaults in non-noco case
      DO n = 1, atoms%ntype
         DO ilSpinPr = MERGE(1,iSpin,noco%l_noco), MERGE(2,iSpin,noco%l_noco)
            CALL timestart("fjgj coefficients")
            CALL fjgjq%calculate(input,atoms,cell,lapwq,noco,usdus,n,ilSpinPr)
            !$acc update device(fjgjq%fj,fjgjq%gj)
            CALL timestop("fjgj coefficients")
            DO ilSpin = ilSpinPr, MERGE(2,iSpin,noco%l_noco)
               CALL timestart("fjgjq coefficients")
               CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,n,ilSpin)
               CALL timestop("fjgjq coefficients")

               IF (.NOT.noco%l_noco) THEN
                  IF (n.EQ.iDtype .AND. juphon%l_phonon) THEN
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,h1mat_tmp(1,1),.TRUE.,lapwq,fjgjq)
                     CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,s1mat_tmp(1,1),h1mat_tmp(1,1),.TRUE.,.TRUE.,lapwq,fjgjq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,h1mat_tmp(1,1),.FALSE.,.TRUE.,.TRUE.,s1mat_tmp(1,1),lapwq,fjgjq)
                  END IF
                  IF (killcont(1)/=0) THEN
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat(1,1),.FALSE.,lapwq,fjgjq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                     !CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,smat(1,1))
                  END IF
               ELSE
                  ! TODO: Everything from here onwards  most certainly has the wrong spin logic.
                  IF (ilSpinPr==ilSpin) THEN !local spin-diagonal contribution
                     CALL hsmt_spinor(ilSpinPr,n,nococonv,chi)
                     IF (n.EQ.iDtype .AND. juphon%l_phonon) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat_tmp,hmat_tmp,.TRUE.,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.TRUE.,smat_tmp,lapwq,fjgjq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                        CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                        CALL timestop("hsmt_distspins")
                     END IF
                     IF (killcont(1)/=0) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,smat)
                        CALL hsmt_distspins(chi,hmat_tmp,hmat)
                        CALL timestop("hsmt_distspins")
                     END IF
                  ELSE IF (noco%l_unrestrictMT(n)) THEN
                     !2,1
                     CALL hsmt_spinor(3,n,nococonv,chi)
                     IF (n.EQ.iDtype .AND. juphon%l_phonon) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                        CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                        CALL timestop("hsmt_distspins")
                     END IF
                     IF (killcont(1)/=0) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.FALSE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,smat)
                        CALL hsmt_distspins(chi,hmat_tmp,hmat)
                        CALL timestop("hsmt_distspins")
                     END IF

                     !1,2
                     CALL hsmt_spinor(4,n,nococonv,chi)
                     IF (n.EQ.iDtype .AND. juphon%l_phonon) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                        CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                        CALL timestop("hsmt_distspins")
                     END IF
                     IF (killcont(1)/=0) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.FALSE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                        CALL timestart("hsmt_distspins")
                        CALL hsmt_distspins(chi,smat_tmp,smat)
                        CALL hsmt_distspins(chi,hmat_tmp,hmat)
                        CALL timestop("hsmt_distspins")
                     END IF
                  END IF
               END IF
            END DO
         END DO
      END DO
      !!$acc end data

      ! TODO: Does this need some ACC magic?
      IF (juphon%l_phonon) THEN
         DO igSpinPr=MERGE(1,1,noco%l_noco),MERGE(2,1,noco%l_noco)
            DO igSpin=MERGE(1,1,noco%l_noco),MERGE(2,1,noco%l_noco)
               CALL matrix_pref(fmpi, atoms, cell%bmat, lapwq%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapwq, lapw, &
                              & nk, lapwq%nv(igSpinPr), lapw%nv(igSpin), iDtype, iDir, &
                              & h1mat_tmp(igSpinPr,igSpin), s1mat_tmp(igSpinPr,igSpin), hmat(igSpinPr,igSpin), smat(igSpinPr,igSpin),killcont(2:3))
               CALL h1mat_tmp(igSpinPr,igSpin)%free()
               CALL s1mat_tmp(igSpinPr,igSpin)%free()
            END DO
         END DO
      END IF
      IF (noco%l_noco) THEN
         !$acc exit data delete(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
         !$acc exit data delete(smat_tmp,hmat_tmp)
      END IF
      RETURN
   END SUBROUTINE dfpt_hsmt

   SUBROUTINE dfpt_dynmat_hsmt(atoms, sym, enpara, iSpin, iDir_row, iDtype_row, iDir_col, iDtype_col, input, fmpi, &
                      & noco, nococonv, cell, lapw, lapwq, usdus, td, tdV1,&
                      hmat1, smat1, hmat1q, smat1q, hmat2, smat2, nk, killcont, vmat2)

      USE m_types
      USE m_types_mpimat
      USE m_hsmt_nonsph
      USE m_hsmt_sph
      USE m_hsmt_lo
      USE m_hsmt_distspins
      USE m_hsmt_fjgj
      USE m_hsmt_spinor
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
      CLASS(t_mat),     INTENT(INOUT) :: hmat1(:,:),smat1(:,:), hmat1q(:,:),smat1q(:,:), hmat2(:,:),smat2(:,:)
      
      CLASS(t_mat), OPTIONAL, INTENT(INOUT) :: vmat2(:,:)

      INTEGER, INTENT(IN) :: iSpin, iDir_row, iDtype_row, iDir_col, iDtype_col, nk, killcont(7)

      TYPE(t_fjgj) :: fjgj, fjgjq

      INTEGER :: ilSpinPr, ilSpin, nspins, i, j
      INTEGER :: igSpinPr, igSpin
      COMPLEX :: chi(2,2),chi_one

      CLASS(t_mat), ALLOCATABLE :: smat_tmp, hmat_tmp, s1mat_tmp(:,:), h1mat_tmp(:,:)
      CLASS(t_mat), ALLOCATABLE :: s1qmat_tmp(:,:), h1qmat_tmp(:,:), s2mat_tmp(:,:), h2mat_tmp(:,:)

      IF (noco%l_noco.AND..NOT.noco%l_ss) THEN
         IF (fmpi%n_size==1) THEN
            ALLOCATE(t_mat::hmat_tmp)
            ALLOCATE(t_mat::smat_tmp)
         ELSE
            ALLOCATE(t_mpimat::hmat_tmp)
            ALLOCATE(t_mpimat::smat_tmp)
         END IF
         CALL smat_tmp%init(hmat1(1,1))
         CALL hmat_tmp%init(hmat1(1,1))
         !$acc enter data copyin(smat_tmp,hmat_tmp)create(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
      END IF

      nspins = MERGE(2, 1, noco%l_noco)
      IF (fmpi%n_size == 1) THEN
         ALLOCATE (t_mat::s1mat_tmp(nspins, nspins), h1mat_tmp(nspins, nspins))
         ALLOCATE (t_mat::s1qmat_tmp(nspins, nspins), h1qmat_tmp(nspins, nspins))
         ALLOCATE (t_mat::s2mat_tmp(nspins, nspins), h2mat_tmp(nspins, nspins))
      ELSE
         ALLOCATE (t_mpimat::s1mat_tmp(nspins, nspins), h1mat_tmp(nspins, nspins))
         ALLOCATE (t_mpimat::s1qmat_tmp(nspins, nspins), h1qmat_tmp(nspins, nspins))
         ALLOCATE (t_mpimat::s2mat_tmp(nspins, nspins), h2mat_tmp(nspins, nspins))
      END IF

      DO i = 1, nspins
         DO j = 1, nspins
            CALL s1mat_tmp(i, j)%init(.FALSE., lapw%nv(i) + atoms%nlotot, lapw%nv(j) + atoms%nlotot, fmpi%sub_comm, .false.)
            CALL h1mat_tmp(i, j)%init(s1mat_tmp(i, j))
            CALL s1qmat_tmp(i, j)%init(.FALSE., lapwq%nv(i) + atoms%nlotot, lapw%nv(j) + atoms%nlotot, fmpi%sub_comm, .false.)
            CALL h1qmat_tmp(i, j)%init(s1qmat_tmp(i, j))
            IF (.NOT.PRESENT(vmat2)) THEN
               CALL s2mat_tmp(i, j)%init(.FALSE., lapw%nv(i) + atoms%nlotot, lapw%nv(j) + atoms%nlotot, fmpi%sub_comm, .false.) 
            ELSE
               CALL s2mat_tmp(i, j)%init(.FALSE., lapwq%nv(i) + atoms%nlotot, lapw%nv(j) + atoms%nlotot, fmpi%sub_comm, .false.)
            END IF
            CALL h2mat_tmp(i, j)%init(s2mat_tmp(i, j))
         END DO
      END DO

      CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,iSpin,noco)
      CALL fjgjq%alloc(MAXVAL(lapwq%nv),atoms%lmaxd,iSpin,noco)
      !$acc data copyin(fjgj) create(fjgj%fj,fjgj%gj)
      !$acc data copyin(fjgjq) create(fjgjq%fj,fjgjq%gj)
      igSpinPr = 1; igSpin = 1; chi_one = 1.0 ! Defaults in non-noco case
      DO ilSpinPr = MERGE(1,iSpin,noco%l_noco), MERGE(2,iSpin,noco%l_noco)
         CALL timestart("fjgj coefficients")
         CALL fjgjq%calculate(input,atoms,cell,lapwq,noco,usdus,iDtype_col,ilSpinPr)
         !$acc update device(fjgjq%fj,fjgjq%gj)
         CALL timestop("fjgj coefficients")
         DO ilSpin = ilSpinPr, MERGE(2,iSpin,noco%l_noco)
            CALL timestart("fjgjq coefficients")
            CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,iDtype_col,ilSpin)
            CALL timestop("fjgjq coefficients")
            IF (.NOT.noco%l_noco) THEN
               CALL hsmt_sph(iDtype_col,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(iDtype_col,ilSpinPr),usdus,fjgj,s1qmat_tmp(1,1),h1qmat_tmp(1,1),.TRUE.,.TRUE.,lapwq,fjgjq)
               CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,h1qmat_tmp(1,1),.FALSE.,lapwq,fjgjq)
               CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,iDtype_col,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,h1qmat_tmp(1,1),.FALSE.,.TRUE.,.TRUE.,s1qmat_tmp(1,1),lapwq,fjgjq)

               CALL hsmt_sph(iDtype_col,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(iDtype_col,ilSpinPr),usdus,fjgj,s1mat_tmp(1,1),h1mat_tmp(1,1),.TRUE.,.TRUE.,lapw,fjgj)
               CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,h1mat_tmp(1,1),.FALSE.,lapw,fjgj)
               CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,iDtype_col,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,h1mat_tmp(1,1),.FALSE.,.TRUE.,.TRUE.,s1mat_tmp(1,1),lapw,fjgj)
               IF (killcont(1)/=0) THEN
                  IF (.NOT.PRESENT(vmat2)) THEN
                     CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,h2mat_tmp(1,1),.FALSE.,lapw,fjgj)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,iDtype_col,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,h2mat_tmp(1,1),.FALSE.,.TRUE.,.FALSE.,lapwq=lapw,fjgjq=fjgj)
                  ELSE
                     CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,h2mat_tmp(1,1),.FALSE.,lapwq,fjgjq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,iDtype_col,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,h2mat_tmp(1,1),.FALSE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                  END IF
               END IF
            ELSE
               RETURN
               ! NOCO_DFPT
               ! TODO: I did not even try to do the right logic here yet.
               IF (ilSpinPr==ilSpin) THEN !local spin-diagonal contribution
                  CALL hsmt_spinor(ilSpinPr,iDtype_col,nococonv,chi)
                  CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                  CALL hsmt_sph(iDtype_col,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(iDtype_col,ilSpinPr),usdus,fjgj,smat_tmp,hmat_tmp,.TRUE.,.TRUE.,lapwq,fjgjq)
                  CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,iDtype_col,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.TRUE.,smat_tmp,lapwq,fjgjq)
                  CALL timestart("hsmt_distspins")
                  CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                  CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                  CALL timestop("hsmt_distspins")
                  IF (killcont(1)/=0) THEN
                     CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,iDtype_col,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat1)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat1)
                     CALL timestop("hsmt_distspins")
                  END IF
               ELSE IF (noco%l_unrestrictMT(iDtype_col)) THEN
                  !2,1
                  CALL hsmt_spinor(3,iDtype_col,nococonv,chi)
                  CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                  CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,iDtype_col,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                  CALL timestart("hsmt_distspins")
                  CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                  CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                  CALL timestop("hsmt_distspins")
                  IF (killcont(1)/=0) THEN
                     CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,iDtype_col,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.FALSE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat1)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat1)
                     CALL timestop("hsmt_distspins")
                  END IF

                  !1,2
                  CALL hsmt_spinor(4,iDtype_col,nococonv,chi)
                  CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                  CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,iDtype_col,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                  CALL timestart("hsmt_distspins")
                  CALL hsmt_distspins(chi,smat_tmp,s1mat_tmp)
                  CALL hsmt_distspins(chi,hmat_tmp,h1mat_tmp)
                  CALL timestop("hsmt_distspins")
                  IF (killcont(1)/=0) THEN
                     CALL hsmt_nonsph(iDtype_col,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.TRUE.,lapwq,fjgjq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,iDtype_col,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.FALSE.,.TRUE.,.FALSE.,lapwq=lapwq,fjgjq=fjgjq)
                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat1)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat1)
                     CALL timestop("hsmt_distspins")
                  END IF
               END IF
            END IF
         END DO
      END DO
      !$acc end data
      !$acc end data

      ! TODO: Does this need some ACC magic?
      DO igSpinPr=MERGE(1,1,noco%l_noco),MERGE(2,1,noco%l_noco)
         DO igSpin=MERGE(1,1,noco%l_noco),MERGE(2,1,noco%l_noco)
            CALL matrix_pref(fmpi, atoms, cell%bmat, lapwq%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapwq, lapw, &
                           & nk, lapwq%nv(igSpinPr), lapw%nv(igSpin), iDtype_col, iDir_col, &
                           & h1qmat_tmp(igSpinPr,igSpin), s1qmat_tmp(igSpinPr,igSpin), hmat1q(igSpinPr,igSpin), smat1q(igSpinPr,igSpin),killcont(2:3))
            CALL h1qmat_tmp(igSpinPr,igSpin)%free()
            CALL s1qmat_tmp(igSpinPr,igSpin)%free()
            CALL matrix_pref(fmpi, atoms, cell%bmat, lapw%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapw, lapw, &
                           & nk, lapw%nv(igSpinPr), lapw%nv(igSpin), iDtype_col, iDir_col, &
                           & h1mat_tmp(igSpinPr,igSpin), s1mat_tmp(igSpinPr,igSpin), hmat1(igSpinPr,igSpin), smat1(igSpinPr,igSpin),killcont(4:5))
            CALL h1mat_tmp(igSpinPr,igSpin)%free()
            CALL s1mat_tmp(igSpinPr,igSpin)%free()
            IF (.NOT.PRESENT(vmat2)) THEN
               CALL matrix_pref(fmpi, atoms, cell%bmat, lapw%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapw, lapw, &
                              & nk, lapw%nv(igSpinPr), lapw%nv(igSpin), iDtype_col, iDir_col, &
                              & h2mat_tmp(igSpinPr,igSpin), s2mat_tmp(igSpinPr,igSpin), hmat2(igSpinPr,igSpin), smat2(igSpinPr,igSpin),[1,0])
            ELSE
               CALL matrix_pref(fmpi, atoms, cell%bmat, lapwq%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapwq, lapw, &
                              & nk, lapwq%nv(igSpinPr), lapw%nv(igSpin), iDtype_col, iDir_col, &       
                              & h2mat_tmp(igSpinPr,igSpin), s2mat_tmp(igSpinPr,igSpin), vmat2(igSpinPr,igSpin), smat1q(igSpinPr,igSpin),[1,0])
            END IF
            CALL h2mat_tmp(igSpinPr,igSpin)%free()
            CALL s2mat_tmp(igSpinPr,igSpin)%free()
            IF (iDtype_row==iDtype_col) THEN
               CALL matrix_pref(fmpi, atoms, cell%bmat, lapw%gvec(:, :, igSpinPr), lapw%gvec(:,:,igSpin), lapw, lapw, &
                              & nk, lapw%nv(igSpinPr), lapw%nv(igSpin), iDtype_row, iDir_row, &
                              & hmat1(igSpinPr,igSpin), smat1(igSpinPr,igSpin), hmat2(igSpinPr,igSpin), smat2(igSpinPr,igSpin),killcont(6:7))
            END IF
         END DO
      END DO
      IF (noco%l_noco) THEN
         !$acc exit data delete(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
         !$acc exit data delete(smat_tmp,hmat_tmp)
      END IF
      RETURN
   END SUBROUTINE dfpt_dynmat_hsmt
END MODULE m_dfpt_hsmt
