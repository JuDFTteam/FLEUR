!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_hsmt
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

  SUBROUTINE hsmt(atoms,sym,enpara,&
       isp,input,fmpi,noco,nococonv,cell,lapw,usdus,td,smat,hmat)
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
    TYPE(t_lapw),INTENT(IN)       :: lapw
    TYPE(t_tlmplm),INTENT(IN)     :: td
    TYPE(t_usdus),INTENT(IN)      :: usdus
    CLASS(t_mat),INTENT(INOUT)    :: smat(:,:),hmat(:,:)
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: isp  !This is the global spin in a collinear calculation

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
              !This is for collinear calculations: the (1,1) element of the matrices is all
              !that is needed and allocated

              CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,CONJG(chi_one),lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat(1,1),hmat(1,1),.FALSE.)
              CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat(1,1),.FALSE.)
              CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(1,1),.FALSE.,smat(1,1))
            ELSEIF(noco%l_noco.AND..NOT.noco%l_ss) THEN
              !The NOCO but non-spinspiral setup follows:
              !The Matrix-elements are first calculated in the local frame of the atom and
              !stored in tmp-variables. Then these are distributed (rotated) into the 2x2
              !global spin-matrices.
              IF (ilSpinPr==ilSpin) THEN !local spin-diagonal contribution
                !initialize the non-LO part of hmat_tmp matrix with zeros
                CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpinPr,1,1,chi_one,noco,nococonv,cell,lapw,td,fjgj,hmat_tmp,.TRUE.)
                !initialize the smat_tmp matrix with zeros
                CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,1,1,CONJG(chi_one),lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat_tmp,hmat_tmp,.TRUE.)
                !initialize the LO part of the matrices with zeros
                CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,n,chi_one,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat_tmp,.TRUE.,smat_tmp)
                CALL hsmt_spinor(ilSpinPr,n,nococonv,chi)
                CALL timestart("hsmt_distspins")
                CALL hsmt_distspins(chi,smat_tmp,smat)
                CALL hsmt_distspins(chi,hmat_tmp,hmat)
                CALL timestop("hsmt_distspins")
              ELSE !Add off-diagonal contributions to Hamiltonian if needed
                IF (noco%l_unrestrictMT(n).OR.noco%l_spinoffd_ldau(n)) THEN
                  CALL hsmt_mtNocoPot_offdiag(n,input,fmpi,sym,atoms,noco,nococonv,cell,lapw,usdus,td,fjgj,igSpinPr,igSpin,hmat_tmp,hmat)
                ENDIF
                IF (noco%l_constrained(n)) CALL hsmt_offdiag(n,atoms,fmpi,nococonv,lapw,td,usdus,fjgj,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat)
                IF (noco%l_soc) CALL hsmt_soc_offdiag(n,atoms,cell,fmpi,nococonv,lapw,sym,usdus,td,fjgj,hmat)
              ENDIF
            ELSE
              !In the spin-spiral case the loop over the interstitial=global spin has to
              !be performed explicitely
              CALL hsmt_spinor(ilSpinPr,n,nococonv,chi)
              DO igSpinPr=1,2
                DO igSpin=1,2
                  IF (ilSpinPr==ilSpin) THEN !local diagonal spin
                    CALL hsmt_sph(n,atoms,fmpi,ilSpinPr,input,nococonv,igSpinPr,igSpin,CONJG(chi(igSpinPr,igSpin)),&
                    lapw,enpara%el0,td%e_shift(n,ilSpinPr),usdus,fjgj,smat(igSpinPr,igSpin),hmat(igSpinPr,igSpin),.FALSE.)
                    CALL hsmt_nonsph(n,fmpi,sym,atoms,ilSpinPr,ilSpin,igSpinPr,igSpin,chi(igSpinPr,igSpin),noco,nococonv,cell,&
                    lapw,td,fjgj,hmat(igSpinPr,igSpin),.FALSE.)
                    CALL hsmt_lo(input,atoms,sym,cell,fmpi,noco,nococonv,lapw,usdus,td,fjgj,&
                    n,chi(igSpinPr,igSpin),ilSpinPr,ilSpin,igSpinPr,igSpin,hmat(igSpinPr,igSpin),.FALSE.,smat(igSpinPr,igSpin))
                  ELSE
                    IF (any(noco%l_unrestrictMT).OR.noco%l_spinoffd_ldau(n)) call hsmt_mtNocoPot_offdiag(n,input,fmpi,sym,atoms,noco,nococonv,cell,lapw,usdus,td,fjgj,igSpinPr,igSpin,hmat_tmp,hmat)
                    IF (any(noco%l_constrained)) CALL hsmt_offdiag(n,atoms,fmpi,nococonv,lapw,td,usdus,fjgj,ilSpinPr,ilSpin,igSpinPr,igSpin,hmat)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      END DO
      !$acc end data
      if (noco%l_noco.AND..NOT.noco%l_ss) then
         !$acc exit data delete(smat_tmp%data_c,smat_tmp%data_r,hmat_tmp%data_c,hmat_tmp%data_r)
         !$acc exit data delete(smat_tmp,hmat_tmp)
      endif
      RETURN
    END SUBROUTINE hsmt
  END MODULE m_hsmt
