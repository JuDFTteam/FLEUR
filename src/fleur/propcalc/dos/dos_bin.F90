!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dosbin

#ifdef CPP_MPI
   USE mpi 
#endif

CONTAINS
   SUBROUTINE dos_bin(jspins, wtkpt, e, eig, qal, g, energyShift)
      !! This subroutine generates the idos, the ldos, the partial
      !! ldos in the spheres and the z-dependent dos for the vacuum.

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: jspins
      REAL,    INTENT(IN)  :: wtkpt(:),e(:)
      REAL,    INTENT(IN)  :: eig(:,:,:),qal(:,:,:)
      REAL,    INTENT(OUT) :: g(:,:)

      REAL,    OPTIONAL, INTENT(IN) :: energyShift

      INTEGER :: nl, k, j, i, js
      REAL    :: de, wk, emin, shift

      de = abs(e(2)-e(1))
      g = 0.0
      shift = 0.0
      IF(PRESENT(energyShift)) shift = energyShift
      emin = minval(e)
      ! put weights in the right bins
      DO js=1, size(g,2)
         DO k = 1 , size(qal,2)
            wk = wtkpt(k)/de
            DO j = 1 , size(eig,1)
               i = NINT((eig(j,k,min(2,js))-shift-emin)/de) + 1   !no third energy in noco case
               IF ( (i.LE.size(g,1)) .AND. (i.GE.1) ) THEN
                  g(i,js) = g(i,js) + wk*qal(j,k,js)* 2.0/jspins
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE dos_bin

   SUBROUTINE dos_bin_transport(fmpi,jspins,wtkpt,egrid,eig,eigq,matrix_element,g,energyShift)
      !! This subroutine generates the idos,
      !! Almost identical to upper subroutine
      !! However, we have have e1 and e2 
      !! Evaluation of delta-distribution is $\delta(\epsilon_{\nu , k }  - shift - \epsilon_{\nu' , k+q })$
      USE m_types_mpi

      IMPLICIT NONE 

      TYPE(t_mpi),INTENT(IN) :: fmpi
      INTEGER, INTENT(IN)  :: jspins
      REAL,    INTENT(IN)  :: wtkpt(:),egrid(:)
      REAL, INTENT(IN)     :: eig(:,:,:) , eigq(:,:,:) , matrix_element(:,:,:,:) !(nu',nu,kpts,jspin)
      REAL, INTENT(OUT)    :: g(:,:) 
      REAL, OPTIONAL, INTENT(IN) :: energyShift
      !LOGICAL , OPTIONAL, INTENT(IN) :: spinDeg


      REAL :: de, wk , emin, shift  , degen
      INTEGER :: jsp, nk , nk_i , nu , iNupr , i 
      !LOGICAL :: l_deg


      !degen = 2.0
      !l_deg = .TRUE. 
      
      !IF (PRESENT(spinDeg)) l_deg = spinDeg
      !IF(.NOT.l_deg) degen = 1.0 
      

      de = abs(egrid(2)-egrid(1))
      g = 0.0
      shift = 0.0 
      IF (PRESENT(energyShift)) shift=energyShift
      emin = MINVAL(egrid)
      DO jsp = 1, SIZE(matrix_element,4)
         DO nk_i = 1 , size(matrix_element,3)
            nk = fmpi%k_list(nk_i)
            wk = wtkpt(nk) / de 
               DO nu = 1 , SIZE(matrix_element,2)
                  DO iNupr = 1 , size(matrix_element,1)
                     ! make sure if we hit emin we are at the first bin
                     i = NINT((eig(nu,nk,jsp) -eigq(iNupr,nk,jsp) -shift - emin)/de) + 1
                     IF ( (i.LE.size(g,1)) .AND. (i.GE.1) ) THEN 
                        g(i,jsp) = g(i,jsp) + wk*matrix_element(iNupr,nu,nk_i,jsp)* 2.0/jspins
                     END IF 
                  END DO ! iNupr
               END DO ! nu 
         END DO  ! k 
      END DO  ! jsp 

   END SUBROUTINE dos_bin_transport

   SUBROUTINE dos_bin_double(fmpi,jspins,wtkpt,egrid,eig,eigq,matrix_element,smearing,g,energyShift)
      !! This subroutine calculates the idos 
      !! Here we have to evalulation we have to take, one for k and one for k'
      !! \delta(eig - shift) \delta(eigq - shift)

      USE m_constants
      USE m_types_mpi
      USE m_smooth
      
      IMPLICIT NONE 

      TYPE(t_mpi) , INTENT(IN) :: fmpi
      INTEGER, INTENT(IN) :: jspins 
      REAL,    INTENT(IN)  :: wtkpt(:),egrid(:)
      REAL, INTENT(IN)     :: eig(:,:,:) , eigq(:,:,:) , matrix_element(:,:,:,:) !(nu',nu,kpts,jspin)
      REAL , INTENT(IN)    :: smearing 
      REAL, INTENT(OUT)    :: g(:,:)
      REAL, OPTIONAL, INTENT(IN) :: energyShift


      REAL :: de, wk , emin, shift  
      REAL, ALLOCATABLE :: gridk(:) , gridq(:)
      INTEGER :: jsp, nk, nk_i , nu , iNupr , i , j  

      ALLOCATE(gridk(size(egrid)))
      ALLOCATE(gridq(size(egrid)))

      


      de = abs(egrid(2)-egrid(1))
      g = 0.0
      shift = 0.0

      IF (PRESENT(energyShift)) shift=energyShift
      emin = MINVAL(egrid)
      DO jsp = 1, SIZE(matrix_element,4)
         DO nk_i = 1 , size(matrix_element,3)
            nk = fmpi%k_list(nk_i)
            wk = wtkpt(nk) / de
               DO nu = 1 , SIZE(matrix_element,2)
                  ! delta(eps_nu(k) - E_F)
                  gridk = 0.0
                  i = NINT((eig(nu,nk,jsp) -shift - emin)/de) + 1
                  IF ( (i.GT.size(gridk,1)) .OR. (i.LT.1) ) CYCLE
                  gridk(i) = 1.0 / de
                  CALL smooth(egrid,gridk,smearing,size(egrid))
                  DO iNupr = 1 , size(matrix_element,1)
                     ! delta(eps_nu'(k+q) - E_F)
                     gridq = 0.0
                     j = NINT((eigq(iNupr,nk,jsp) -shift - emin)/de) + 1
                     IF ( (j.GT.size(gridq,1)) .OR. (j.LT.1) ) CYCLE
                     gridq(j) = 1.0 / de
                     CALL smooth(egrid,gridq,smearing,size(egrid))
                     g(:,jsp) = g(:,jsp) + wk*de*matrix_element(iNupr,nu,nk_i,jsp)*2.0/jspins &
                                          * gridk(:)*gridq(:)
                  END DO ! iNupr
               END DO ! nu
         END DO  ! k
      END DO  ! jsp

   END SUBROUTINE dos_bin_double
END MODULE m_dosbin
