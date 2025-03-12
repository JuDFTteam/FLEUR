!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_cdnval

USE m_juDFT
#ifdef CPP_MPI
use mpi
#endif

CONTAINS

SUBROUTINE cdnval_mt()

   TYPE(t_input),INTENT(IN)             :: input
   TYPE(t_usdus),INTENT(IN)             :: usdus
   TYPE(t_lapw),INTENT(IN)              :: lapw
    
   TYPE(t_noco),INTENT(IN)              :: noco
   TYPE(t_nococonv),INTENT(IN)          :: nococonv
   TYPE(t_sym),INTENT(IN)               :: sym
   TYPE(t_cell),INTENT(IN)              :: cell
   TYPE(t_atoms),INTENT(IN)             :: atoms
   TYPE(t_mat),INTENT(IN)               :: zMat
   TYPE(t_force),OPTIONAL,INTENT(INOUT) :: force

   noccbd,eig,sphhar,we,ne

   CALL timestart("cdnval_mt")

      DO ispin = lbound(denmatrix,1),ubound(denmatrix,1)
         CALL abcof(input,atoms,sym,cell,lapw,noccbd,usdus,noco,nococonv,ispin,&
                    eigVecCoeffs%abcof(:,0:,0,:,ispin),eigVecCoeffs%abcof(:,0:,1,:,ispin),&
                    eigVecCoeffs%ccof(-atoms%llod:,:,:,:,ispin),zMat,eig,force)
        
        
      END DO ! end loop over ispin
      IF (noco%l_mperp) then
        call timestart("denCoeffsOffdiag%calcCoefficients")
      
        call timestop("denCoeffsOffdiag%calcCoefficients")
      endif

     
   CALL timestop("cdnval_mt")

END SUBROUTINE cdnval_mt

END MODULE m_cdnval
