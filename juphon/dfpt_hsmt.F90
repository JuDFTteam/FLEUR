!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_hsmt
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  !> Setup of MT-part of the Hamiltonian and the overlap matrix
  !!
  !! Here the MT-components are added to the matrices.
  !! 1. The spherical part in hsmt_sph()
  !! 2. The non-spherical part in hsmt_nonsph()
  !! 3. The lo-part in hsmt_lo()
  !!
  !! - In the case of a noco-calculation (but not spin-spiral), first a temporary matrix is set-up
  !! for each atom in its local spin-frame and this matrix is the rotated into the global frame and added to the full matrix
  !! - In the spin-spiral case, a loop over the global spin is performed and the four parts of the matrix are calculated one-by-one
  !! @todo
  !! The off-diagonal contribution in first-variation soc and constraint calculations is still missing
  ! DFPT: Handle H/S elements with shifted k+q on the lhs i(k+G-k'-G'-q) prefactor through sph,
  !       V1 elements fully with nonsph. LO?

  SUBROUTINE dfpt_hsmt(atoms,sym,enpara,&
       isp,iDir,iDtype,input,fmpi,noco,nococonv,cell,lapw,lapwq,usdus,td,smat,hmat,tdV1)
    USE m_types
    USE m_types_mpimat
    USE m_hsmt_nonsph
    USE m_hsmt_sph
    USE m_hsmt_lo
    USE m_hsmt_distspins
    USE m_hsmt_fjgj
    USE m_hsmt_spinor
    USE m_hsmt_soc_offdiag
    USE m_hsmt_mtNocoPot_offdiag
    USE m_hsmt_offdiag
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: fmpi
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
    TYPE(t_nococonv),INTENT(IN)   :: nococonv
    TYPE(t_sym),INTENT(IN)        :: sym
    TYPE(t_cell),INTENT(IN)       :: cell
    TYPE(t_atoms),INTENT(IN)      :: atoms
    TYPE(t_enpara),INTENT(IN)     :: enpara
    TYPE(t_lapw),INTENT(IN)       :: lapw,lapwq
    TYPE(t_tlmplm),INTENT(IN)     :: td, tdV1
    TYPE(t_usdus),INTENT(IN)      :: usdus
    CLASS(t_mat),INTENT(INOUT)    :: smat(:,:),hmat(:,:)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp, iDir, iDtype  !This is the global spin in a collinear calculation

    !locals
    TYPE(t_fjgj)::fjgj
    INTEGER :: ilSpinPr,ilSpin !local spin in atom
    INTEGER :: igSpinPr,igSpin,n
    COMPLEX :: chi(2,2),chi_one

    CLASS(t_mat),ALLOCATABLE::smat_tmp,hmat_tmp

    !
    IF (noco%l_noco.AND..NOT.noco%l_ss) THEN
       IF (fmpi%n_size==1) THEN
          ALLOCATE(t_mat::hmat_tmp)
          ALLOCATE(t_mat::smat_tmp)
       ELSE
          ALLOCATE(t_mpimat::hmat_tmp)
          ALLOCATE(t_mpimat::smat_tmp)
       ENDIF
       CALL smat_tmp%init(hmat(1,1))
       CALL hmat_tmp%init(hmat(1,1))
       !$acc enter data copyin(smat_tmp,hmat_tmp)create(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
    ENDIF

    CALL fjgj%alloc(MAXVAL(lapw%nv),atoms%lmaxd,isp,noco)
    !$acc data copyin(fjgj) create(fjgj%fj,fjgj%gj)
    igSpinPr=1;igSpin=1;chi_one=1.0 !Defaults in non-noco case
    DO n=1,atoms%ntype
       DO ilSpinPr=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
          CALL timestart("fjgj coefficients")
          CALL fjgj%calculate(input,atoms,cell,lapw,noco,usdus,n,ilSpinPr)
          !$acc update device(fjgj%fj,fjgj%gj)
          CALL timestop("fjgj coefficients")
          DO ilSpin=ilSpinPr,MERGE(2,isp,noco%l_noco)
               IF (.NOT.noco%l_noco) THEN
                  IF (n.EQ.iDtype) THEN
                     CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat(1,1),hmat(1,1),.FALSE.,.TRUE.,cell%bmat,iDir,lapwq)
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat(1,1),.FALSE.,.TRUE.,iDir,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,iDir,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,.TRUE.,.TRUE.,smat(1,1),lapwq)
                  END IF
                  CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat(1,1),.FALSE.,.FALSE.,iDir,lapwq)
                  CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,iDir,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,.FALSE.,.TRUE.,lapwq=lapwq)
                  !CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,smat(1,1))
               ELSE
                  IF (ilSpinPr==ilSpin) THEN !local spin-diagonal contribution
                     IF (n.EQ.iDtype) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,.TRUE.,iDir,lapwq)
                        CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat_tmp,hmat_tmp,.TRUE.,.TRUE.,cell%bmat,iDir,lapwq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,iDir,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.TRUE.,smat_tmp,lapwq)
                     END IF
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.FALSE.,.FALSE.,iDir,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,iDir,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.FALSE.,.FALSE.,.TRUE.,lapwq=lapwq)
                     CALL hsmt_spinor(ilSpinPr,n,nococonv,chi)
                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat)
                     CALL timestop("hsmt_distspins")
                  ELSE IF (noco%l_unrestrictMT(n)) THEN
                     !2,1
                     IF (n.EQ.iDtype) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,.TRUE.,iDir,lapwq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,iDir,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.TRUE.,lapwq=lapwq)
                     END IF
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,2,1,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.FALSE.,.FALSE.,iDir,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,iDir,chi_one,2,1,igSpinPr,igSpin,hmat_tmp,.FALSE.,.FALSE.,.TRUE.,lapwq=lapwq)
                     CALL hsmt_spinor(3,n,nococonv,chi)
                     CALL timestart("hsmt_distspins")
                     CALL hsmt_distspins(chi,smat_tmp,smat)
                     CALL hsmt_distspins(chi,hmat_tmp,hmat)
                     CALL timestop("hsmt_distspins")

                     !1,2
                     IF (n.EQ.iDtype) THEN
                        CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.,.TRUE.,iDir,lapwq)
                        CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,iDir,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.TRUE.,.TRUE.,.TRUE.,lapwq=lapwq)
                     END IF
                     CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,2,chi_one,noco,nococonv,cell,lapw,tdV1,fjgj,hmat_tmp,.FALSE.,.FALSE.,iDir,lapwq)
                     CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,tdV1,fjgj,n,iDir,chi_one,1,2,igSpinPr,igSpin,hmat_tmp,.FALSE.,.FALSE.,.TRUE.,lapwq=lapwq)
                     CALL hsmt_spinor(3,n,nococonv,chi)
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
      IF (noco%l_noco) THEN
         !$acc exit data delete(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
         !$acc exit data delete(smat_tmp,hmat_tmp)
      END IF
      RETURN
   END SUBROUTINE dfpt_hsmt
END MODULE m_dfpt_hsmt
