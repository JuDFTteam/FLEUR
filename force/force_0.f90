MODULE m_force0
CONTAINS
  SUBROUTINE force_0(results)
! ************************************************************
! Presetting force components
! Also presets lm quantities for Pulay terms
! al la Yu et al equa A12,A17,A20
! ************************************************************
!
    USE m_types
      IMPLICIT NONE 
      TYPE(t_results),INTENT(INOUT)   :: results
!     ..
!     preset force components to zero and save force
!
      results%force_old(:,:)=results%force(:,:,1)
      !sum second spin if present
      if (size(results%force,3)>1) results%force_old(:,:)=results%force_old(:,:)+ results%force(:,:,2)
      results%force=0.0
      END SUBROUTINE force_0
      END MODULE m_force0
