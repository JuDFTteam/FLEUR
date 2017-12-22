!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_mt_setup

CONTAINS
  SUBROUTINE mt_setup(jsp,atoms,sym,sphhar,input,noco,enpara,v,mpi,results,DIMENSION,td,ud)
    USE m_usetup
    USE m_tlmplm_cholesky
    USE m_tlmplm_store
     USE m_types
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
     INTEGER,INTENT(IN)           :: jsp
     
    COMPLEX,ALLOCATABLE:: vs_mmp(:,:,:,:)
    INTEGER:: isp
    INTEGER, PARAMETER :: lmaxb=3

    
    IF ((atoms%n_u.GT.0)) THEN
       ALLOCATE( vs_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u,input%jspins) )
       CALL u_setup(sym,atoms,lmaxb,sphhar,input, enpara%el0(0:,:,:),v%mt,mpi, vs_mmp,results)
    ELSE
       ALLOCATE( vs_mmp(-lmaxb:-lmaxb,-lmaxb:-lmaxb,1,2) )
    ENDIF
    
    CALL timestart("tlmplm")
    CALL td%init(DIMENSION%lmplmd,DIMENSION%lmd,atoms%ntype,atoms%lmaxd,atoms%llod,SUM(atoms%nlo),dot_PRODUCT(atoms%nlo,atoms%nlo+1)/2,MERGE(2,1,noco%l_noco))

    !CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,mpi,v%mt(:,0,1,jsp),input,vs_mmp, td,ud)
    CALL tlmplm_cholesky(sphhar,atoms,DIMENSION,enpara, jsp,1,mpi,v%mt(:,0,1,jsp),input, td,ud)
    IF (input%l_f) CALL write_tlmplm(td,vs_mmp,atoms%n_u>0,1,jsp,input%jspins)
    CALL timestop("tlmplm")
    
    !---> pk non-collinear
    !--->       call tlmplm again for the second spin direction in
    !--->       each MT, because the t-matrices are needed for both
    !--->       spins at once in hsmt
    IF (noco%l_noco) THEN
       isp = 2
       CALL timestart("tlmplm")
       !CALL tlmplm(sphhar,atoms,DIMENSION,enpara,isp,isp,mpi, v%mt(1,0,1,isp),lh0,input, td,ud)
       IF (input%l_f) CALL write_tlmplm(td,vs_mmp,atoms%n_u>0,2,2,input%jspins)
       CALL timestop("tlmplm")
    ENDIF
    
  END SUBROUTINE mt_setup
END MODULE m_mt_setup
