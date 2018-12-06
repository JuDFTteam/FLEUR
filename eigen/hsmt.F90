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
       ispin,input,mpi,noco,cell,lapw,usdus,td,smat,hmat)
    USE m_types
    USE m_hsmt_nonsph
    USE m_hsmt_sph
    use m_hsmt_lo
    USE m_hsmt_distspins
    USE m_hsmt_fjgj
    USE m_hsmt_spinor
    USE m_hsmt_soc_offdiag
    USE m_hsmt_mtNocoPot_offdiag
    USE m_hsmt_offdiag
    IMPLICIT NONE
    TYPE(t_mpi),INTENT(IN)        :: mpi
    TYPE(t_input),INTENT(IN)      :: input
    TYPE(t_noco),INTENT(IN)       :: noco
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
    INTEGER, INTENT (IN) :: ispin  
    
    !locals
#ifdef CPP_GPU
    REAL, ALLOCATABLE,MANAGED    :: fj(:,:,:,:),gj(:,:,:,:)
#else
    REAL, ALLOCATABLE    :: fj(:,:,:,:),gj(:,:,:,:)
#endif

    INTEGER :: iintsp,jintsp,n
    COMPLEX :: chi(2,2),chi_one

    TYPE(t_mat)::smat_tmp,hmat_tmp

    !
    IF (noco%l_noco.AND..NOT.noco%l_ss) THEN
       CALL smat_tmp%alloc(smat(1,1)%l_real,smat(1,1)%matsize1,smat(1,1)%matsize2)
       CALL hmat_tmp%alloc(smat(1,1)%l_real,smat(1,1)%matsize1,smat(1,1)%matsize2)
    ENDIF
    
    ALLOCATE(fj(MAXVAL(lapw%nv),0:atoms%lmaxd,input%jspins,MERGE(2,1,noco%l_noco)))
    ALLOCATE(gj(MAXVAL(lapw%nv),0:atoms%lmaxd,input%jspins,MERGE(2,1,noco%l_noco)))

    iintsp=1;jintsp=1;chi_one=1.0 !Defaults in non-noco case
       DO n=1,atoms%ntype
          CALL timestart("fjgj coefficients")
          CALL hsmt_fjgj(input,atoms,cell,lapw,noco,usdus,n,ispin,fj,gj)
          CALL timestop("fjgj coefficients")
          IF (.NOT.noco%l_noco) THEN
             !This is for collinear calculations: the (1,1) element of the matrices is all
             !that is needed and allocated
             CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,1,1,chi_one,lapw,enpara%el0,&
                           td%e_shift(n,ispin),usdus,fj(:,0:,ispin,:),gj(:,0:,ispin,:),smat(1,1),hmat(1,1))
             CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,1,1,chi_one,noco,cell,lapw,td,&
                              fj(:,0:,ispin,:),gj(:,0:,ispin,:),hmat(1,1))
             CALL hsmt_lo(input,atoms,sym,cell,mpi,noco,lapw,usdus,td,fj(:,0:,ispin,:),gj(:,0:,ispin,:),&
                          n,chi_one,ispin,iintsp,jintsp,hmat(1,1),smat(1,1))
          ELSEIF(noco%l_noco.AND..NOT.noco%l_ss) THEN
             !The NOCO but non-spinspiral setup follows:
             !The Matrix-elements are first calculated in the local frame of the atom and
             !stored in tmp-variables. Then these are distributed (rotated) into the 2x2
             !global spin-matrices.
             CALL hmat_tmp%clear();CALL smat_tmp%clear()
             CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,1,1,chi_one,lapw,enpara%el0,td%e_shift(n,ispin),&
                  usdus,fj(:,0:,ispin,:),gj(:,0:,ispin,:),smat_tmp,hmat_tmp)
             CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,1,1,chi_one,noco,cell,lapw,td,&
                  fj(:,0:,ispin,:),gj(:,0:,ispin,:),hmat_tmp)
             CALL hsmt_lo(input,atoms,sym,cell,mpi,noco,lapw,usdus,td,fj(:,0:,ispin,:),gj(:,0:,ispin,:),&
                  n,chi_one,ispin,iintsp,jintsp,hmat_tmp,smat_tmp)
             CALL hsmt_spinor(ispin,n,noco,chi)
             CALL hsmt_distspins(chi,smat_tmp,smat)
             CALL hsmt_distspins(chi,hmat_tmp,hmat)
             !Add off-diagonal contributions to Hamiltonian if needed
             IF (ispin==1.AND.noco%l_mtNocoPot) THEN
                CALL hsmt_mtNocoPot_offdiag()
             ENDIF
             IF (ispin==1.and.noco%l_soc) &
                  CALL hsmt_soc_offdiag(n,atoms,mpi,noco,lapw,usdus,td,fj(:,0:,:,iintsp),gj(:,0:,:,iintsp),hmat)
             IF (noco%l_constr) &
                  CALL hsmt_offdiag(n,atoms,mpi,ispin,noco,lapw,td,usdus,fj(:,0:,ispin,:),gj(:,0:,ispin,:),hmat)
          ELSE
             !In the spin-spiral case the loop over the interstitial=global spin has to
             !be performed explicitely
             CALL hsmt_spinor(ispin,n,noco,chi)
             DO iintsp=1,2
                DO jintsp=1,2
                   CALL hsmt_sph(n,atoms,mpi,ispin,input,noco,iintsp,jintsp,chi(iintsp,jintsp),&
                                 lapw,enpara%el0,td%e_shift(n,ispin),usdus,fj(:,0:,ispin,:),gj(:,0:,ispin,:),&
                                 smat(iintsp,jintsp),hmat(iintsp,jintsp))
                   CALL hsmt_nonsph(n,mpi,sym,atoms,ispin,iintsp,jintsp,chi(iintsp,jintsp),noco,cell,&
                                    lapw,td,fj(:,0:,ispin,:),gj(:,0:,ispin,:),hmat(iintsp,jintsp))
                   CALL hsmt_lo(input,atoms,sym,cell,mpi,noco,lapw,usdus,td,fj(:,0:,ispin,:),gj(:,0:,ispin,:),&
                          n,chi(iintsp,jintsp),ispin,iintsp,jintsp,hmat(iintsp,jintsp),smat(iintsp,jintsp))
                ENDDO
             ENDDO
          ENDIF
          
    END DO

    
    RETURN
  END SUBROUTINE hsmt
END MODULE m_hsmt
