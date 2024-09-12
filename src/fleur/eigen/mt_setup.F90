!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mt_setup

CONTAINS
  SUBROUTINE mt_setup(atoms,sym,sphhar,input,noco,nococonv,enpara,hub1inp,hub1data,inden,vTot,vx,fmpi,td,ud,alpha_hybrid,l_dfptmod)
    USE m_types
    USE m_tlmplm_cholesky
    USE m_spnorb
    IMPLICIT NONE

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
    LOGICAL,INTENT(IN),OPTIONAL :: l_dfptmod

    INTEGER:: jsp


   
  !l_offdiag

    
       !CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,fmpi,vTot%mt(:,0,1,jsp),input,vTot%mmpMat, td,ud)
       IF (PRESENT(l_dfptmod)) THEN
          CALL tlmplm_cholesky(sphhar,atoms,sym,noco,nococonv,enpara,jsp,fmpi,vTot,vx,inDen,input,hub1inp,hub1data,td,ud,alpha_hybrid,l_dfptmod)
       ELSE
         CALL tlmplm_cholesky(sphhar,atoms,sym,noco,nococonv,enpara,jsp,fmpi,vTot,vx,inDen,input,hub1inp,hub1data,td,ud,alpha_hybrid)
      END IF
   
    CALL timestop("tlmplm")

    

  END SUBROUTINE mt_setup
END MODULE m_mt_setup
