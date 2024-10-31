MODULE m_dosbin
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
      DO js=1, size(qal,3)
         DO k = 1 , size(qal,2)
            wk = wtkpt(k)/de
            DO j = 1 , size(eig,1)
               i = NINT((eig(j,k,js)-shift-emin)/de) + 1
               IF ( (i.LE.size(g,1)) .AND. (i.GE.1) ) THEN
                  g(i,js) = g(i,js) + wk*qal(j,k,js)* 2.0/jspins
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE dos_bin

   SUBROUTINE dos_bin_transport(jspins,wtkpt,egrid,eig,eigq,matrix_element,g,energyShift)
      !! This subroutine generates the idos,
      !! Almost identical to upper subroutine
      !! However, we have have e1 and e2 
      !! Evaluation of delta-distribution is $\delta(\epsilon_{\nu , k }  - shift - \epsilon_{\nu' , k+q })$
      IMPLICIT NONE 

      INTEGER, INTENT(IN)  :: jspins
      REAL,    INTENT(IN)  :: wtkpt(:),egrid(:)
      REAL, INTENT(IN)     :: eig(:,:,:) , eigq(:,:,:) , matrix_element(:,:,:,:) !(nu',nu,kpts,jspin)
      REAL, INTENT(OUT)    :: g(:,:) 
      REAL, OPTIONAL, INTENT(IN) :: energyShift
      !LOGICAL , OPTIONAL, INTENT(IN) :: spinDeg


      REAL :: de, wk , emin, shift  , degen
      INTEGER :: jsp, nk , nu , iNupr , i 
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
         DO nk = 1 , size(matrix_element,3)
            wk = wtkpt(nk) / de 
               DO nu = 1 , SIZE(matrix_element,2)
                  DO iNupr = 1 , size(matrix_element,1)
                     ! make sure if we hit emin we are at the first bin
                     i = NINT((eig(nu,nk,jsp) -eigq(iNupr,nk,jsp) -shift - emin)/de) + 1
                     IF ( (i.LE.size(g,1)) .AND. (i.GE.1) ) THEN 
                        g(i,jsp) = g(i,jsp) + wk*matrix_element(iNupr,nu,nk,jsp)* 2.0/jspins
                     END IF 
                  END DO ! iNupr
               END DO ! nu 
         END DO  ! k 
      END DO  ! jsp 

   END SUBROUTINE dos_bin_transport

   SUBROUTINE dos_bin_double(jspins,wtkpt,egrid,eig,eigq,matrix_element,smearing,g,energyShift)
      !! This subroutine calculates the idos 
      !! Here we have to evalulation we have to take, one for k and one for k'
      !! We already do gaussian smearing in this routine

      !! WARNING: still needs to be tested
      USE m_constants
      
      IMPLICIT NONE 

      INTEGER, INTENT(IN) :: jspins 
      REAL,    INTENT(IN)  :: wtkpt(:),egrid(:)
      REAL, INTENT(IN)     :: eig(:,:,:) , eigq(:,:,:) , matrix_element(:,:,:,:) !(nu',nu,kpts,jspin)
      REAL , INTENT(IN)    :: smearing 
      REAL, INTENT(OUT)    :: g(:,:)
      REAL, OPTIONAL, INTENT(IN) :: energyShift


      REAL :: de, wk , emin, shift  
      REAL, ALLOCATABLE :: gauss(:) , gaussq(:)
      INTEGER :: jsp, nk , nu , iNupr , i , j , grid2 
      LOGICAL :: l_true

      ALLOCATE(gauss(size(egrid)))
      ALLOCATE(gaussq(size(egrid)))
      gauss = 0.0
      gaussq = 0.0



      de = abs(egrid(2)-egrid(1))
      g = 0.0
      shift = 0.0 

      IF (PRESENT(energyShift)) shift=energyShift
      emin = MINVAL(egrid)
      DO jsp = 1, SIZE(matrix_element,4)
         DO nk = 1 , size(matrix_element,3)
            wk = wtkpt(nk) / de 
               DO nu = 1 , SIZE(matrix_element,2)
                  gauss(:) =  1/SQRT(tpi_const*smearing) *EXP( -(eGrid(:)-eig(nu,nk,jsp) + energyShift )**2/(2*smearing)) 
                  DO iNupr = 1 , size(matrix_element,1)
                     gaussq(:) = 1/SQRT(tpi_const*smearing) *EXP( -(eGrid(:)-eigq(iNupr,nk,jsp) + energyShift )**2/(2*smearing)) 
                     g(:,jsp) = g(:,jsp) + wk*matrix_element(iNupr,nu,nk,jsp) * gauss(:) * gaussq(:) * 2.0/jspins 
                  END DO ! iNupr
               END DO ! nu 
         END DO  ! k 
      END DO  ! jsp 

   END SUBROUTINE dos_bin_double
END MODULE m_dosbin
