!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_hssetup
CONTAINS
  SUBROUTINE eigen_hssetup(isp,mpi,DIMENSION,oned,hybrid,enpara,input,vacuum,noco,jij,sym,&
       stars,cell,kpts,sphhar,atoms,ud,td,v,bkpt,lapw,smat_final,hmat_final)
    USE m_hs_int
    USE m_hsvac
    USE m_od_hsvac
    USE m_hsmt
    USE m_types
    IMPLICIT NONE
    INTEGER,INTENT(IN)           :: isp
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_oneD),INTENT(IN)      :: oneD
    TYPE(t_hybrid),INTENT(INOUT) :: hybrid
    TYPE(t_enpara),INTENT(INOUT) :: enpara
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_jij),INTENT(IN)       :: jij
    TYPE(t_sym),INTENT(IN)       :: sym  
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_kpts),INTENT(IN)      :: kpts
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_usdus),INTENT(IN)     :: ud
    TYPE(t_tlmplm),INTENT(IN)    :: td
    TYPE(t_lapw),INTENT(IN)      :: lapw
    TYPE(t_potden),INTENT(IN)    :: v
    TYPE(t_lapwmat),INTENT(INOUT)  :: smat_final,hmat_final
    REAL,INTENT(IN)              :: bkpt(3)
    

    
    TYPE(t_lapwmat),ALLOCATABLE :: smat(:,:),hmat(:,:)
    INTEGER :: i,j,ispin
    
    !Matrices for Hamiltonian and Overlapp
    !In noco case we need 4-matrices for each spin channel
    IF (noco%l_noco) THEN
       ALLOCATE(smat(2,2),hmat(2,2))
       DO i=1,2
          DO j=1,2
             CALL smat(i,j)%ALLOC(.FALSE.,lapw%nv(i)+atoms%nlotot,lapw%num_local_cols(j))
             CALL hmat(i,j)%ALLOC(.FALSE.,lapw%nv(i)+atoms%nlotot,lapw%num_local_cols(j))
          ENDDO
       ENDDO
    ELSE
       ALLOCATE(smat(1,1),hmat(1,1))
       CALL smat(1,1)%alloc(smat_final%l_real,lapw%nv(1)+atoms%nlotot,lapw%num_local_cols(1))
       CALL hmat(1,1)%alloc(smat_final%l_real,lapw%nv(1)+atoms%nlotot,lapw%num_local_cols(1))
    ENDIF
    CALL timestart("Interstitial part")
    !Generate interstitial part of Hamiltonian
    CALL hs_int(input,noco,stars,lapw,mpi,cell,isp,bkpt,v%pw,smat,hmat)
    CALL timestop("Interstitial part")
    CALL timestart("MT part")
      !MT-part of Hamiltonian. In case of noco, we need an loop over the local spin of the atoms
    DO ispin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
       CALL hsmt(atoms,sphhar,sym,enpara,ispin,input,mpi,noco,cell,lapw,ud,td,smat,hmat)
    ENDDO
    CALL timestop("MT part")
   
    !Vacuum contributions
    IF (input%film) THEN
       CALL timestart("Vacuum part")
       CALL hsvac(vacuum,stars,DIMENSION, atoms,mpi,isp,input,v,enpara%evac0,cell,&
            lapw,sym, noco,jij,hmat,smat)
       CALL timestop("Vacuum part")
    ENDIF
    !Now copy the data into final matrix
    IF (mpi%n_size==1) THEN
       IF (.NOT.noco%l_noco) THEN
          smat_final%matsize1=smat(1,1)%matsize1
          smat_final%matsize2=smat(1,1)%matsize2
          hmat_final%matsize1=smat(1,1)%matsize1
          hmat_final%matsize2=smat(1,1)%matsize2
          IF(smat_final%l_real) THEN
             CALL move_ALLOC(hmat(1,1)%data_r,hmat_final%data_r)
             CALL move_ALLOC(smat(1,1)%data_r,smat_final%data_r)
          ELSE
             CALL move_ALLOC(hmat(1,1)%data_c,hmat_final%data_c)
             CALL move_ALLOC(smat(1,1)%data_c,smat_final%data_c)
          ENDIF
       ELSE
          CALL smat_final%alloc(.FALSE.,lapw%nv_tot+2*atoms%nlotot,lapw%nv_tot+2*atoms%nlotot)
          CALL hmat_final%alloc(.FALSE.,lapw%nv_tot+2*atoms%nlotot,lapw%nv_tot+2*atoms%nlotot)
          !up-up
          smat_final%data_c(:lapw%nv(1)+atoms%nlotot,:lapw%nv(1)+atoms%nlotot)=smat(1,1)%data_c
          hmat_final%data_c(:lapw%nv(1)+atoms%nlotot,:lapw%nv(1)+atoms%nlotot)=hmat(1,1)%data_c
          !down-down
          smat_final%data_c(lapw%nv(1)+atoms%nlotot+1:,lapw%nv(1)+atoms%nlotot+1:)=smat(2,2)%data_c
          hmat_final%data_c(lapw%nv(1)+atoms%nlotot+1:,lapw%nv(1)+atoms%nlotot+1:)=hmat(2,2)%data_c
          !off-diag
          DO i=1,smat(1,2)%matsize1 !First map U-part of smat&hmat(2,1) into smat(1,2)
             DO j=i+1,smat(1,2)%matsize2
                smat(2,1)%data_c(j,i)=CONJG(smat(1,2)%data_c(i,j))
                hmat(2,1)%data_c(j,i)=CONJG(hmat(1,2)%data_c(i,j))
             ENDDO
          ENDDO
          smat_final%data_c(:lapw%nv(1)+atoms%nlotot,lapw%nv(1)+atoms%nlotot+1:)=smat(2,1)%data_c
          hmat_final%data_c(:lapw%nv(1)+atoms%nlotot,lapw%nv(1)+atoms%nlotot+1:)=hmat(2,1)%data_c
       ENDIF
    ELSE
       !CALL eigen_redist_matrix(mpi,lapw,atoms,noco,smat_final,smat)
       !CALL eigen_redist_matrix(mpi,lapw,atoms,noco,hmat_final,hmat)
    ENDIF
  END SUBROUTINE eigen_hssetup
END MODULE m_eigen_hssetup
          
       
