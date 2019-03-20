!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mt_setup

CONTAINS
  SUBROUTINE mt_setup(atoms,sym,sphhar,input,noco,enpara,gOnsite,inden,vTot,mpi,results,DIMENSION,td,ud)
    USE m_types
    USE m_usetup
    USE m_tlmplm_cholesky
    USE m_tlmplm_store
    USE m_spnorb
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_mpi),INTENT(IN)       :: mpi
    TYPE(t_dimension),INTENT(IN) :: DIMENSION
    TYPE(t_enpara),INTENT(INOUT) :: enpara
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_sym),INTENT(IN)       :: sym  
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_greensf),INTENT(IN)   :: gOnsite
    TYPE(t_potden),INTENT(IN)    :: inDen
    TYPE(t_potden),INTENT(INOUT) :: vTot
    TYPE(t_tlmplm),INTENT(INOUT) :: td
    TYPE(t_usdus),INTENT(INOUT)  :: ud

    INTEGER:: jsp

    IF ((atoms%n_u.GT.0)) THEN
       CALL u_setup(sym,atoms,sphhar,input,enpara%el0(0:,:,:),inDen,vTot,mpi,results)
    END IF
    !IF((atoms%n_hia.GT.0)) THEN
    !   CALL hia_setup(gOnsite,sym,atoms,sphhar,input,enpara%el0(0:,:,:),inDen,vTot,mpi,results)
    !END IF

    CALL timestart("tlmplm")
    CALL td%init(DIMENSION%lmplmd,DIMENSION%lmd,atoms%ntype,atoms%lmaxd,atoms%llod,SUM(atoms%nlo),&
         DOT_PRODUCT(atoms%nlo,atoms%nlo+1)/2,MERGE(4,input%jspins,noco%l_mtNocoPot),&
         (noco%l_noco.AND.noco%l_soc.AND..NOT.noco%l_ss).OR.noco%l_constr)!l_offdiag

    DO jsp=1,MERGE(4,input%jspins,noco%l_mtNocoPot)
       !CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,mpi,vTot%mt(:,0,1,jsp),input,vTot%mmpMat, td,ud)
       CALL tlmplm_cholesky(sphhar,atoms,noco,enpara,jsp,jsp,mpi,vTot,input,td,ud)
       IF (input%l_f) CALL write_tlmplm(td,vTot%mmpMat,atoms%n_u>0,jsp,jsp,input%jspins)
    END DO
    CALL timestop("tlmplm")

    !Setup of soc parameters for first-variation SOC
    IF (noco%l_soc.AND.noco%l_noco.AND..NOT.noco%l_ss) THEN
       CALL spnorb(atoms,noco,input,mpi,enpara,vTot%mt,ud,td%rsoc,.FALSE.)
    END IF

  END SUBROUTINE mt_setup
END MODULE m_mt_setup
