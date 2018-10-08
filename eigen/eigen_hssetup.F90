!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_eigen_hssetup
CONTAINS
  !> The setup of the Hamiltonian and Overlap matrices are performed here
  !!
  !! The following steps are executed:
  !! 1. The matrices are a allocated (in the noco-case these are 2x2-arrays of matrices)
  !! 2. The Interstitial contribution is calculated (in hs_int())
  !! 3. The MT-part is calculated (in hsmt() )
  !! 4. The vacuum part is added (in hsvac())
  !! 5. The matrices are copied to the final matrix, in the noco-case the full matrix is constructed from the 4-parts.
  
  SUBROUTINE eigen_hssetup(isp,mpi,DIMENSION,hybrid,enpara,input,vacuum,noco,sym,&
       stars,cell,sphhar,atoms,ud,td,v,lapw,l_real,smat_final,hmat_final)
    USE m_types
    USE m_types_mpimat
    USE m_types_gpumat
    USE m_hs_int
    USE m_hsvac
    USE m_od_hsvac
    USE m_hsmt
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
       IF (judft_was_argument("-gpu")) THEN
          ALLOCATE(t_gpumat::smat(nspins,nspins),hmat(nspins,nspins))
       ELSE
          ALLOCATE(t_mat::smat(nspins,nspins),hmat(nspins,nspins))
       ENDIF
    ELSE
       ALLOCATE(t_mpimat::smat(nspins,nspins),hmat(nspins,nspins))
    ENDIF
    DO i=1,nspins
       DO j=1,nspins
          CALL smat(i,j)%init(l_real,lapw%nv(i)+atoms%nlotot,lapw%nv(j)+atoms%nlotot,mpi%sub_comm,.false.)
          CALL hmat(i,j)%init(smat(i,j))
       ENDDO
    ENDDO

    
    CALL timestart("Interstitial part")
    !Generate interstitial part of Hamiltonian
    CALL hs_int(input,noco,stars,lapw,mpi,cell,isp,v%pw_w,smat,hmat)
    CALL timestop("Interstitial part")
    CALL timestart("MT part")
      !MT-part of Hamiltonian. In case of noco, we need an loop over the local spin of the atoms
    DO ispin=MERGE(1,isp,noco%l_noco),MERGE(2,isp,noco%l_noco)
       CALL hsmt(atoms,sym,enpara,ispin,input,mpi,noco,cell,lapw,ud,td,smat,hmat)
    ENDDO
    CALL timestop("MT part")
   
    !Vacuum contributions
    IF (input%film) THEN
       CALL timestart("Vacuum part")
       CALL hsvac(vacuum,stars,DIMENSION, atoms,mpi,isp,input,v,enpara%evac,cell,&
            lapw,sym, noco,hmat,smat)
       CALL timestop("Vacuum part")
    ENDIF
    !Now copy the data into final matrix
    ! Collect the four noco parts into a single matrix
    ! In collinear case only a copy is done
    ! In the parallel case also a redistribution happens
    ALLOCATE(smat_final,source=smat(1,1))
    ALLOCATE(hmat_final,source=smat(1,1))
    CALL eigen_redist_matrix(mpi,lapw,atoms,smat,smat_final)
    CALL eigen_redist_matrix(mpi,lapw,atoms,hmat,hmat_final,smat_final)
    
  END SUBROUTINE eigen_hssetup
END MODULE m_eigen_hssetup
       
