!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mt_setup

CONTAINS
  SUBROUTINE mt_setup(atoms,sym,sphhar,input,noco,enpara,v,mpi,results,DIMENSION,td,ud)
    USE m_usetup
    USE m_tlmplm_cholesky
    USE m_tlmplm_store
    USE m_types
    USE m_socinit
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
    TYPE(t_potden),INTENT(IN)    :: v
    TYPE(t_tlmplm),INTENT(INOUT) :: td
    TYPE(t_usdus),INTENT(INOUT)  :: ud

    COMPLEX,ALLOCATABLE:: vs_mmp(:,:,:,:)
    INTEGER:: jsp
    INTEGER, PARAMETER :: lmaxb=3


    IF ((atoms%n_u.GT.0)) THEN
       ALLOCATE( vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins) )
       CALL u_setup(sym,atoms,lmaxb,sphhar,input, enpara%el0(0:,:,:),v%mt,mpi, vs_mmp,results)
    ELSE
       ALLOCATE( vs_mmp(-lmaxb:-lmaxb,-lmaxb:-lmaxb,1,2) )
    ENDIF

    CALL timestart("tlmplm")
    CALL td%init(DIMENSION%lmplmd,DIMENSION%lmd,atoms%ntype,atoms%lmaxd,atoms%llod,SUM(atoms%nlo),&
         DOT_PRODUCT(atoms%nlo,atoms%nlo+1)/2,input%jspins,&
         (noco%l_noco.AND.noco%l_soc.AND..NOT.noco%l_ss).OR.noco%l_constr)!l_offdiag

    DO jsp=1,input%jspins
       !CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,mpi,v%mt(:,0,1,jsp),input,vs_mmp, td,ud)
       CALL tlmplm_cholesky(sphhar,atoms,noco,enpara, jsp,jsp,mpi,v%mt,input,vs_mmp, td,ud)
       IF (input%l_f) CALL write_tlmplm(td,vs_mmp,atoms%n_u>0,1,jsp,input%jspins)
    END DO
    CALL timestop("tlmplm")

    !Setup of soc parameters for first-variation SOC
    IF (noco%l_soc.AND.noco%l_noco.AND..NOT.noco%l_ss) &
         CALL socinit(mpi,atoms,sphhar,enpara,input,v%mt,noco,ud,td%rsoc)
  
    

  END SUBROUTINE mt_setup
END MODULE m_mt_setup
