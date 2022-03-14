!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mt_setup

CONTAINS
  SUBROUTINE mt_setup(atoms,sym,sphhar,input,noco,nococonv,enpara,hub1inp,hub1data,inden,vTot,vx,fmpi,results,td,ud,alpha_hybrid)
    USE m_types
    USE m_tlmplm_cholesky
    USE m_spnorb
    IMPLICIT NONE
    TYPE(t_results),INTENT(INOUT):: results
    TYPE(t_mpi),INTENT(IN)       :: fmpi

    TYPE(t_enpara),INTENT(IN) :: enpara
    TYPE(t_input),INTENT(IN)     :: input
    TYPE(t_noco),INTENT(IN)      :: noco
    TYPE(t_nococonv),INTENT(IN)  :: nococonv
    TYPE(t_sym),INTENT(IN)       :: sym
    TYPE(t_sphhar),INTENT(IN)    :: sphhar
    TYPE(t_atoms),INTENT(IN)     :: atoms
    TYPE(t_potden),INTENT(IN)    :: inDen
    TYPE(t_potden),INTENT(IN)    :: vTot,vx
    TYPE(t_tlmplm),INTENT(INOUT) :: td
    TYPE(t_usdus),INTENT(INOUT)  :: ud
    TYPE(t_hub1inp),INTENT(IN)   :: hub1inp
    TYPE(t_hub1data),INTENT(INOUT)::hub1data
    REAl,INTENT(IN)               :: alpha_hybrid

    INTEGER:: jsp


    CALL timestart("tlmplm")
    CALL td%init(atoms,input%jspins,(noco%l_noco.AND.noco%l_soc.AND..NOT.noco%l_ss).OR.any(noco%l_constrained))!l_offdiag

    DO jsp=1,MERGE(4,input%jspins,any(noco%l_unrestrictMT).OR.any(noco%l_spinoffd_ldau))
       !CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,fmpi,vTot%mt(:,0,1,jsp),input,vTot%mmpMat, td,ud)
       CALL tlmplm_cholesky(sphhar,atoms,sym,noco,nococonv,enpara,jsp,fmpi,vTot,vx,inDen,input,hub1inp,hub1data,td,ud,alpha_hybrid)
    END DO
    CALL timestop("tlmplm")

    !Setup of soc parameters for first-variation SOC
    IF (noco%l_soc.AND.noco%l_noco.AND..NOT.noco%l_ss) THEN
       CALL spnorb(atoms,noco,nococonv,input,fmpi,enpara,vTot%mt,ud,td%rsoc,.FALSE.,hub1inp,hub1data)
    END IF

  END SUBROUTINE mt_setup
END MODULE m_mt_setup
