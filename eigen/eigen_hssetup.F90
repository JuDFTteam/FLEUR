!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_hssetup
CONTAINS
  SUBROUTINE eigen_hssetup(isp,mpi,DIMENSION,hybrid,enpara,input,vacuum,noco,jij,sym,&
       stars,cell,sphhar,atoms,ud,td,v,lapw,l_real,smat_final,hmat_final)
    USE m_hs_int
    USE m_hsvac
    USE m_od_hsvac
    USE m_hsmt
    USE m_types
    USE m_types_mpimat
    USE m_eigen_redist_matrix
    IMPLICIT NONE
    INTEGER,INTENT(IN)           :: isp
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_hybrid),INTENT(IN)    :: hybrid
    TYPE(t_enpara),INTENT(IN)    :: enpara
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_vacuum),INTENT(IN)    :: vacuum
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_jij),INTENT(IN)       :: jij
    TYPE(t_sym),INTENT(IN)       :: sym  
    TYPE(t_stars),INTENT(IN)     :: stars
    TYPE(t_cell),INTENT(IN)      :: cell
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_usdus),INTENT(IN)     :: ud
    TYPE(t_tlmplm),INTENT(IN)    :: td
    TYPE(t_lapw),INTENT(IN)      :: lapw
    TYPE(t_potden),INTENT(IN)    :: v
    CLASS(t_mat),ALLOCATABLE,INTENT(INOUT)   :: smat_final,hmat_final
    LOGICAL,INTENT(IN)           :: l_real
    

    
    CLASS(t_mat),ALLOCATABLE :: smat(:,:),hmat(:,:)
    INTEGER :: i,j,ispin,nspins
    
    !Matrices for Hamiltonian and Overlapp
    !In noco case we need 4-matrices for each spin channel
    nspins=MERGE(2,1,noco%l_noco)
    IF (mpi%n_size==1) THEN
       ALLOCATE(t_mat::smat(nspins,nspins),hmat(nspins,nspins))
    ELSE
       ALLOCATE(t_mpimat::smat(nspins,nspins),hmat(nspins,nspins))
    ENDIF
    DO i=1,nspins
       DO j=1,nspins
          CALL smat(i,j)%init(l_real,lapw%nv(i)+atoms%nlotot,lapw%nv(i)+atoms%nlotot,mpi%sub_comm,.false.)
          CALL hmat(i,j)%init(l_real,lapw%nv(i)+atoms%nlotot,lapw%nv(i)+atoms%nlotot,mpi%sub_comm,.false.)
       ENDDO
    ENDDO

    
    CALL timestart("Interstitial part")
    !Generate interstitial part of Hamiltonian
    CALL hs_int(input,noco,stars,lapw,mpi,cell,isp,v%pw,smat,hmat)
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
       ALLOCATE(t_mat::smat_final,hmat_final)
    ELSE
       ALLOCATE(t_mpimat::smat_final,hmat_final)
    ENDIF
    
    CALL eigen_redist_matrix(mpi,lapw,atoms,smat,smat_final)
    CALL eigen_redist_matrix(mpi,lapw,atoms,hmat,hmat_final)
    
  END SUBROUTINE eigen_hssetup
END MODULE m_eigen_hssetup

#if 1==2
IF (mpi%n_size==1) THEN
       IF (.NOT.noco%l_noco) THEN
          smat_final%matsize1=smat(1,1)%matsize1
          smat_final%matsize2=smat(1,1)%matsize2
          hmat_final%matsize1=smat(1,1)%matsize1
          hmat_final%matsize2=smat(1,1)%matsize2
          IF(l_real) THEN
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
#endif
       
