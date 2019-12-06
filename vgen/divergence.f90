!--------------------------------------------------------------------------------  
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_divergence
   USE m_types
   PRIVATE
   PUBLIC :: mt_div, pw_div, divergence, divergence2, mt_grad, pw_grad, divpotgrad, divpotgrad2

CONTAINS
   SUBROUTINE mt_div(jspins,n,atoms,sphhar,noco,sym,xcB,div)
   !-----------------------------------------------------------------------------!
   !By use of the cartesian components of a field, its radial/angular derivati-  !
   !ves in the muffin tin at each spherical grid point and the corresponding an- !
   !gles:                                                                        !
   !                                                                             !
   !Make the divergence of said field in real space and store it as a source     !
   !density, again represented by mt-coefficients in a potden.                   !
   !                                                                             !
   !Code by A. Neukirchen, September 2019                                        !
   !-----------------------------------------------------------------------------! 
   USE m_mt_tofrom_grid

   IMPLICIT NONE
   
   INTEGER, INTENT(IN)                         :: jspins, n
   TYPE(t_atoms), INTENT(IN)                   :: atoms
   TYPE(t_sphhar), INTENT(IN)                  :: sphhar
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_noco), INTENT(IN)                    :: noco
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: xcB
   TYPE(t_potden), INTENT(INOUT)               :: div

   TYPE(t_gradients)                           :: gradx, grady, gradz

   REAL, ALLOCATABLE :: div_temp(:, :)
   REAL, ALLOCATABLE :: thet(:), phi(:)
   REAL :: r,th,ph
   INTEGER :: jr, k, nsp, kt

   nsp = atoms%nsp()

   ALLOCATE (gradx%gr(3,atoms%jri(n)*nsp,jspins),grady%gr(3,atoms%jri(n)*nsp,jspins),gradz%gr(3,atoms%jri(n)*nsp,jspins))
   ALLOCATE (div_temp(atoms%jri(n)*nsp,jspins))

   CALL init_mt_grid(jspins, atoms, sphhar, .TRUE., sym)

   CALL mt_to_grid(.TRUE., jspins, atoms, sphhar, xcB(1)%mt(:,0:,n,:), n, noco,gradx)
   CALL mt_to_grid(.TRUE., jspins, atoms, sphhar, xcB(2)%mt(:,0:,n,:), n, noco,grady)
   CALL mt_to_grid(.TRUE., jspins, atoms, sphhar, xcB(3)%mt(:,0:,n,:), n, noco,gradz)

   kt = 0
   DO jr = 1, atoms%jri(n)
      r=atoms%rmsh(jr, n)
      DO k = 1, nsp
         th = thet(k)
         ph = phi(k)
         div_temp(kt+nsp,1) = (SIN(th)*COS(ph)*gradx%gr(1,kt+nsp,jspins) + SIN(th)*SIN(ph)*grady%gr(1,kt+nsp,jspins) + COS(th)*gradz%gr(1,kt+nsp,jspins))&
                             +(COS(th)*COS(ph)*gradx%gr(2,kt+nsp,jspins) + COS(th)*SIN(ph)*grady%gr(2,kt+nsp,jspins) - SIN(th)*gradz%gr(2,kt+nsp,jspins))/r&
                             -(SIN(ph)*gradx%gr(3,kt+nsp,jspins)         + COS(ph)*grady%gr(3,kt+nsp,jspins))/(r*SIN(th))
      ENDDO ! k
   
      kt = kt+nsp
   ENDDO ! jr
    
   CALL mt_from_grid(atoms, sphhar, n, jspins, div_temp, div%mt(:,0:,n,:))
   
   CALL finish_mt_grid

   END SUBROUTINE mt_div

   SUBROUTINE pw_div(ifftxc3,jspins,stars,cell,noco,sym,xcB,div)
   !-----------------------------------------------------------------------------!
   !By use of the cartesian components of a field, its radial/angular derivati-  !
   !ves in the muffin tin at each spherical grid point and the corresponding an- !
   !gles:                                                                        !
   !                                                                             !
   !Make the divergence of said field in real space and store it as a source     !
   !density, again represented by mt-coefficients in a potden.                   !
   !                                                                             !
   !Code by A. Neukirchen, September 2019                                        !
   !-----------------------------------------------------------------------------! 
   USE m_pw_tofrom_grid

   IMPLICIT NONE
   
   INTEGER, INTENT(IN)                         :: jspins, ifftxc3
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_noco), INTENT(IN)                    :: noco
   TYPE(t_stars),INTENT(IN)                    :: stars
   TYPE(t_cell),INTENT(IN)                     :: cell
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: xcB
   TYPE(t_potden), INTENT(INOUT)               :: div

   TYPE(t_gradients)                           :: gradx, grady, gradz

   REAL, ALLOCATABLE :: div_temp(:, :)
   INTEGER :: i, nsp

   nsp = 3*ifftxc3

   ALLOCATE (gradx%gr(3,nsp,jspins),grady%gr(3,nsp,jspins),gradz%gr(3,nsp,jspins))
   ALLOCATE (div_temp(nsp,jspins))

   CALL init_pw_grid(.TRUE.,stars,sym,cell)

   CALL pw_to_grid(.TRUE.,jspins,noco%l_noco,stars,cell,xcB(1)%pw,gradx)
   CALL pw_to_grid(.TRUE.,jspins,noco%l_noco,stars,cell,xcB(2)%pw,grady)
   CALL pw_to_grid(.TRUE.,jspins,noco%l_noco,stars,cell,xcB(3)%pw,gradz)

   DO i = 1, nsp
         div_temp(i,1)=gradx%gr(1,i,1)+grady%gr(2,i,1)+gradz%gr(3,i,1)
   ENDDO ! i
    
   CALL pw_from_grid(.TRUE.,stars,.TRUE.,div_temp,div%pw,div%pw_w)

   CALL finish_pw_grid()
   
   END SUBROUTINE pw_div

   SUBROUTINE divergence(jspins,n,ifftxc3,atoms,sphhar,sym,stars,cell,vacuum,noco,xcB,div)
   USE m_types
   IMPLICIT NONE
   
   INTEGER, INTENT(IN)                         :: jspins, n, ifftxc3
   TYPE(t_atoms), INTENT(IN)                   :: atoms
   TYPE(t_sphhar), INTENT(IN)                  :: sphhar
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_noco), INTENT(IN)                    :: noco
   TYPE(t_stars),INTENT(IN)                    :: stars
   TYPE(t_cell),INTENT(IN)                     :: cell
   TYPE(t_vacuum),INTENT(IN)                   :: vacuum
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: xcB
   TYPE(t_potden), INTENT(OUT)                 :: div

   CALL div%init(stars,atoms,sphhar,vacuum,noco,jspins,1001)

   CALL mt_div(jspins,n,atoms,sphhar,noco,sym,xcB,div)
   CALL pw_div(ifftxc3,jspins,stars,cell,noco,sym,xcB,div)
      
   END SUBROUTINE divergence

   SUBROUTINE divergence2(stars,atoms,sphhar,vacuum,sym,cell,noco,bxc,div)
      USE m_lh_tofrom_lm
      USE m_gradYlm

      !--------------------------------------------------------------------------
      ! Use the interstitial/vacuum divergence subroutine and an external MT-gra-
      ! dient routine from juPhon to assemble the divergence of a field into a
      ! t_potden variable. The MT-gradient is first calculated in sperical har-
      ! monics coefficients.                                         
      !--------------------------------------------------------------------------   

      IMPLICIT NONE

      TYPE(t_stars),                INTENT(IN)    :: stars
      TYPE(t_atoms),                INTENT(IN)    :: atoms
      TYPE(t_sphhar),               INTENT(IN)    :: sphhar
      TYPE(t_vacuum),               INTENT(IN)    :: vacuum
      TYPE(t_sym),                  INTENT(IN)    :: sym
      TYPE(t_cell),                 INTENT(IN)    :: cell
      TYPE(t_noco),                 INTENT(IN)    :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(INOUT) :: bxc
      TYPE(t_potden),               INTENT(INOUT) :: div

      INTEGER :: i,iType,indmax, lh
      COMPLEX, ALLOCATABLE :: flm(:,:,:),grsflm1(:,:,:,:),grsflm2(:,:,:,:),grsflm3(:,:,:,:),divflm(:,:,:) ! (iR,lm,n[,x,i])

      indmax=(atoms%lmaxd+1)**2
      
      ALLOCATE(flm(atoms%jmtd,indmax,atoms%ntype))
      ALLOCATE(divflm(atoms%jmtd,indmax,atoms%ntype))

      DO i=1,3
         DO iType=1, atoms%ntype
            CALL lh_to_lm(atoms, sphhar, iType, bxc(i)%mt(:,:,iType,1), flm(:,:,iType))
         END DO
         IF (i==1) THEN
            CALL gradYlm(atoms,flm,grsflm1)
         ELSE IF (i==2) THEN
            CALL gradYlm(atoms,flm,grsflm2)
         ELSE
            CALL gradYlm(atoms,flm,grsflm3)
         END IF
      END DO

      DEALLOCATE(flm)

      CALL divYlm(grsflm1(:,:,:,:),grsflm2(:,:,:,:),grsflm3(:,:,:,:), divflm)

      DO iType=1, atoms%ntype
         CALL lh_from_lm(atoms, sphhar, iType, divflm(:,1:indmax,iType), div%mt(:,0:,iType,1))
      END DO

      DEALLOCATE(divflm,grsflm1,grsflm2,grsflm3)

      CALL pw_div(stars,sym,cell,noco,bxc,div)

      
   END SUBROUTINE divergence2


   SUBROUTINE mt_grad(n,atoms,sphhar,sym,den,gradphi)
      !-----------------------------------------------------------------------------
      !By use of the cartesian components of a field, its radial/angular derivati-  
      !ves in the muffin tin at each spherical grid point and the corresponding an- 
      !gles:                                                                        
      !                                                                             
      !Make the divergence of said field in real space and store it as a source     
      !density, again represented by mt-coefficients in a potden.                   
      !                                                                             
      !Code by A. Neukirchen, September 2019                                        
      !----------------------------------------------------------------------------- 
      USE m_mt_tofrom_grid
      USE m_constants

      IMPLICIT NONE
   
      INTEGER, INTENT(IN)                         :: n
      TYPE(t_atoms), INTENT(IN)                   :: atoms
      TYPE(t_sphhar), INTENT(IN)                  :: sphhar
      TYPE(t_sym), INTENT(IN)                     :: sym
      TYPE(t_potden), INTENT(IN)                  :: den
      TYPE(t_potden), dimension(3), INTENT(INOUT) :: gradphi

      TYPE(t_potden)                              :: denloc
      TYPE(t_gradients)                           :: grad

      REAL, ALLOCATABLE :: thet(:), phi(:), grad_temp(:, :, :)
      REAL :: r, th, ph, eps
      INTEGER :: i, jr, k, nsp, kt, lh, lhmax

      nsp = atoms%nsp()
      lhmax=sphhar%nlh(atoms%ntypsy(SUM(atoms%neq(:n - 1)) + 1))
      eps=1.e-10

      ALLOCATE (grad%gr(3,atoms%jri(n)*nsp,1))
      ALLOCATE (grad_temp(atoms%jri(n)*nsp,1,3))
      ALLOCATE (thet(atoms%nsp()),phi(atoms%nsp()))

      denloc=den

      DO lh=0, lhmax
         denloc%mt(:,lh,n,1) = denloc%mt(:,lh,n,1)*atoms%rmsh(:, n)**2
      END DO ! lh

      CALL init_mt_grid(1, atoms, sphhar, .TRUE., sym, thet, phi)

      CALL mt_to_grid(.TRUE., 1, atoms, sphhar, denloc%mt(:,0:,n,:), n, grad)

      kt = 0
      DO jr = 1, atoms%jri(n)
         r=atoms%rmsh(jr, n)
         DO k = 1, nsp
            th = thet(k)
            ph = phi(k)
            grad_temp(kt+k,1,1) = (SIN(th)*COS(ph)*grad%gr(1,kt+k,1) + COS(th)*COS(ph)*grad%gr(2,kt+k,1)/r - SIN(ph)*grad%gr(3,kt+k,1)/(r*SIN(th)))/(4.0*pi_const)
            grad_temp(kt+k,1,2) = (SIN(th)*SIN(ph)*grad%gr(1,kt+k,1) + COS(th)*SIN(ph)*grad%gr(2,kt+k,1)/r + COS(ph)*grad%gr(3,kt+k,1)/(r*SIN(th)))/(4.0*pi_const)
            grad_temp(kt+k,1,3) = (        COS(th)*grad%gr(1,kt+k,1) -         SIN(th)*grad%gr(2,kt+k,1)/r                                        )/(4.0*pi_const)
         ENDDO ! k
         kt = kt+nsp
      ENDDO ! jr

      DO i=1,3 
         CALL mt_from_grid(atoms, sphhar, n, 1, grad_temp(:,:,i), gradphi(i)%mt(:,0:,n,:))
         DO lh=0, lhmax
            gradphi(i)%mt(:,lh,n,1) = gradphi(i)%mt(:,lh,n,1)*atoms%rmsh(:, n)**2
            !IF ((sphhar%llh(lh,1)/=0).AND.(sphhar%llh(lh,1)/=2)) THEN
            !   gradphi(i)%mt(:,lh,n,1) = 0.0
            !END IF
         END DO ! lh
      END DO ! i

      CALL finish_mt_grid
   
   END SUBROUTINE mt_grad

   SUBROUTINE pw_grad(stars,cell,noco,sym,den,gradphi)
      !-----------------------------------------------------------------------------
      !By use of the cartesian components of a field and its cartesian derivatives  
      !in the interstitial/vacuum region at each grid point:                        
      !                                                                             
      !Make the divergence of said field in real space and store it as a source     
      !density, again represented by pw-coefficients in a potden.                   
      !                                                                             
      !Code by A. Neukirchen, September 2019                                        
      !----------------------------------------------------------------------------- 
      USE m_constants
      USE m_pw_tofrom_grid

      IMPLICIT NONE
   
      TYPE(t_sym), INTENT(IN)                     :: sym
      TYPE(t_noco), INTENT(IN)                    :: noco
      TYPE(t_stars),INTENT(IN)                    :: stars
      TYPE(t_cell),INTENT(IN)                     :: cell
      TYPE(t_potden), dimension(3), INTENT(INOUT) :: gradphi
      TYPE(t_potden), INTENT(IN)                  :: den

      TYPE(t_gradients)                           :: grad

      REAL, ALLOCATABLE :: grad_temp(:, :, :)
      INTEGER :: i,ifftxc3

      ifftxc3=stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft

      ALLOCATE (grad_temp(ifftxc3,1,3))

      CALL init_pw_grid(.TRUE.,stars,sym,cell)

      CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,den%pw,grad)

      DO i = 1, ifftxc3
            grad_temp(i,1,:)=grad%gr(:,i,1)/(4.0*pi_const)
      ENDDO ! i

      DO i=1,3 
         CALL pw_from_grid(.TRUE.,stars,.TRUE.,grad_temp(:,:,i),gradphi(i)%pw,gradphi(i)%pw_w)
      END DO

      CALL finish_pw_grid()
   
   END SUBROUTINE pw_grad

   SUBROUTINE divpotgrad(stars,atoms,sphhar,vacuum,sym,cell,noco,pot,grad)

      USE m_types
      USE m_constants
      USE m_mt_tofrom_grid
      IMPLICIT NONE
      !-----------------------------------------------------------------------------
      !Use the two gradient subroutines above to now put the complete gradient      
      !of a potential into a t_potden variable.                                     
      !-----------------------------------------------------------------------------  
      TYPE(t_stars),INTENT(IN)                    :: stars 
      TYPE(t_atoms), INTENT(IN)                   :: atoms
      TYPE(t_sphhar), INTENT(IN)                  :: sphhar
      TYPE(t_vacuum),INTENT(IN)                   :: vacuum
      TYPE(t_sym), INTENT(IN)                     :: sym
      TYPE(t_cell),INTENT(IN)                     :: cell
      TYPE(t_noco), INTENT(IN)                    :: noco
      TYPE(t_potden), INTENT(IN)                  :: pot
      TYPE(t_potden), dimension(3), INTENT(INOUT) :: grad

      INTEGER :: n

      DO n=1,atoms%ntype
         CALL mt_grad(n,atoms,sphhar,sym,pot,grad)
      END DO
      CALL pw_grad(stars,cell,noco,sym,pot,grad)

   END SUBROUTINE divpotgrad

   SUBROUTINE divpotgrad2(stars,atoms,sphhar,vacuum,sym,cell,noco,pot,grad)

      USE m_types
      USE m_lh_tofrom_lm
      USE m_gradYlm
      USE m_constants

      !--------------------------------------------------------------------------
      ! Use the interstitial/vacuum gradient subroutine and an external MT-gra-
      ! dient routine from juPhon to assemble the gradient of a potenital into a
      ! t_potden variable. The MT-gradient is first calculated in sperical har-
      ! monics coefficients.                                         
      !--------------------------------------------------------------------------   

      IMPLICIT NONE
 
      TYPE(t_stars),INTENT(IN)                    :: stars 
      TYPE(t_atoms), INTENT(IN)                   :: atoms
      TYPE(t_sphhar), INTENT(IN)                  :: sphhar
      TYPE(t_vacuum),INTENT(IN)                   :: vacuum
      TYPE(t_sym), INTENT(IN)                     :: sym
      TYPE(t_cell),INTENT(IN)                     :: cell
      TYPE(t_noco), INTENT(IN)                    :: noco
      TYPE(t_potden), INTENT(IN)                  :: pot
      TYPE(t_potden), dimension(3), INTENT(INOUT) :: grad

      TYPE(t_potden)                              :: denloc
      INTEGER :: i,iType,indmax,lh,lhmax
      COMPLEX, ALLOCATABLE :: flm(:,:,:),grsflm(:,:,:,:) ! (iR,lm,n[,x,i])

      indmax=(atoms%lmaxd+1)**2

      ALLOCATE(flm(atoms%jmtd,indmax,atoms%ntype))

      denloc=pot

      DO iType=1,atoms%ntype
         lhmax=sphhar%nlh(atoms%ntypsy(SUM(atoms%neq(:iType - 1)) + 1))
         DO lh=0, lhmax
            denloc%mt(:,lh,iType,1) = denloc%mt(:,lh,iType,1)*atoms%rmsh(:, iType)**2
         END DO ! lh
         CALL lh_to_lm(atoms, sphhar, iType, denloc%mt(:,:,iType,1), flm(:,:,iType))
      END DO

      CALL gradYlm(atoms,flm,grsflm)

      DEALLOCATE(flm)
    
      DO i=1,3
         DO iType=1,atoms%ntype
            CALL lh_from_lm(atoms, sphhar, iType, grsflm(:,1:indmax,iType,i)/(4.0*pi_const), grad(i)%mt(:,0:,iType,1))
         END DO
      END DO

      DEALLOCATE(grsflm)

      CALL pw_grad(stars,cell,noco,sym,pot,grad)

   END SUBROUTINE divpotgrad2
END MODULE m_divergence
