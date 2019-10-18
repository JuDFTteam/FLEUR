!--------------------------------------------------------------------------------  
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_divergence
   USE m_types
   PRIVATE
   PUBLIC :: mt_div, pw_div, divergence, mt_grad, pw_grad, divpotgrad

CONTAINS
   SUBROUTINE mt_div(n,atoms,sphhar,sym,bxc,div)
      USE m_mt_tofrom_grid

      !--------------------------------------------------------------------------
      !By use of the cartesian components of a field, its radial/angular derivati-  
      !ves in the muffin tin at each spherical grid point and the corresponding an- 
      !gles:                                                                        
      !                                                                             
      !Make the divergence of said field in real space and store it as a source     
      !density, again represented by mt-coefficients in a potden.                   
      !                                                                             
      !Code by A. Neukirchen, September 2019                                        
      !-------------------------------------------------------------------------- 

      IMPLICIT NONE
   
      INTEGER,                      INTENT(IN)    :: n
      TYPE(t_atoms),                INTENT(IN)    :: atoms
      TYPE(t_sphhar),               INTENT(IN)    :: sphhar
      TYPE(t_sym),                  INTENT(IN)    :: sym
      TYPE(t_potden), DIMENSION(3), INTENT(IN)    :: bxc
      TYPE(t_potden),               INTENT(INOUT) :: div
   
      TYPE(t_gradients)                           :: gradx, grady, gradz

      REAL, ALLOCATABLE :: thet(:), phi(:), div_temp(:, :) 
      REAL, ALLOCATABLE :: xcBx_mt(:,:,:), xcBy_mt(:,:,:), xcBz_mt(:,:,:)
      REAL :: r, th, ph
      INTEGER :: jr, k, nsp, kt, i

      nsp = atoms%nsp()

      ALLOCATE (gradx%gr(3,atoms%jri(n)*nsp,1),grady%gr(3,atoms%jri(n)*nsp,1), &
                gradz%gr(3,atoms%jri(n)*nsp,1))
      ALLOCATE (div_temp(atoms%jri(n)*nsp,1))
      ALLOCATE (thet(atoms%nsp()),phi(atoms%nsp()))

      ALLOCATE (xcBx_mt(lbound(bxc(1)%mt, dim=1):ubound(bxc(1)%mt, dim=1), &
                        lbound(bxc(1)%mt, dim=2):ubound(bxc(1)%mt, dim=2), &
                        lbound(bxc(1)%mt, dim=4):ubound(bxc(1)%mt, dim=4)))

      ALLOCATE (xcBy_mt(lbound(bxc(2)%mt, dim=1):ubound(bxc(2)%mt, dim=1), &
                        lbound(bxc(2)%mt, dim=2):ubound(bxc(2)%mt, dim=2), &
                        lbound(bxc(2)%mt, dim=4):ubound(bxc(2)%mt, dim=4)))

      ALLOCATE (xcBz_mt(lbound(bxc(3)%mt, dim=1):ubound(bxc(3)%mt, dim=1), &
                        lbound(bxc(3)%mt, dim=2):ubound(bxc(3)%mt, dim=2), &
                        lbound(bxc(3)%mt, dim=4):ubound(bxc(3)%mt, dim=4)))

      CALL init_mt_grid(1, atoms, sphhar, .TRUE., sym, thet, phi)

      !DO jr=1,atoms%jri(n)
      !   xcBx_mt(jr,0:,:) = xcB(1)%mt(jr,0:,n,:) * atoms%rmsh(jr,n)**2
      !   xcBy_mt(jr,0:,:) = xcB(2)%mt(jr,0:,n,:) * atoms%rmsh(jr,n)**2
      !   xcBz_mt(jr,0:,:) = xcB(3)%mt(jr,0:,n,:) * atoms%rmsh(jr,n)**2
      !ENDDO

      CALL mt_to_grid(.TRUE., 1, atoms, sphhar, bxc(1)%mt(:,0:,n,:), n, gradx)
      CALL mt_to_grid(.TRUE., 1, atoms, sphhar, bxc(2)%mt(:,0:,n,:), n, grady)
      CALL mt_to_grid(.TRUE., 1, atoms, sphhar, bxc(3)%mt(:,0:,n,:), n, gradz)

      kt = 0
      DO jr = 1, atoms%jri(n)
         r =atoms%rmsh(jr, n)
         DO k = 1, nsp
            th = thet(k)
            ph = phi(k)
            div_temp(kt+k,1) = (SIN(th)*COS(ph)*gradx%gr(1,kt+k,1) + SIN(th)*SIN(ph)*grady%gr(1,kt+k,1) + COS(th)*gradz%gr(1,kt+k,1))*r**2&
                               +(COS(th)*COS(ph)*gradx%gr(2,kt+k,1) + COS(th)*SIN(ph)*grady%gr(2,kt+k,1) - SIN(th)*gradz%gr(2,kt+k,1))*r&
                               -(SIN(ph)*gradx%gr(3,kt+k,1)         - COS(ph)*grady%gr(3,kt+k,1))*r/SIN(th)
         ENDDO ! k
         kt = kt+nsp
      ENDDO ! jr
    
      CALL mt_from_grid(atoms, sphhar, n, 1, div_temp, div%mt(:,0:,n,:))
   
      CALL finish_mt_grid
   
   END SUBROUTINE mt_div

   SUBROUTINE pw_div(stars,sym,cell,noco,bxc,div)
      USE m_pw_tofrom_grid
      !--------------------------------------------------------------------------
      !By use of the cartesian components of a field and its cartesian derivatives  
      !in the interstitial/vacuum region at each grid point:                        
      !                                                                             
      !Make the divergence of said field in real space and store it as a source     
      !density, again represented by pw-coefficients in a potden.                   
      !                                                                             
      !Code by A. Neukirchen, September 2019                                        
      !--------------------------------------------------------------------------


      IMPLICIT NONE

      TYPE(t_stars),                INTENT(IN)    :: stars   
      TYPE(t_sym),                  INTENT(IN)    :: sym
      TYPE(t_cell),                 INTENT(IN)    :: cell
      TYPE(t_noco),                 INTENT(IN)    :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(IN)    :: bxc
      TYPE(t_potden),               INTENT(INOUT) :: div

      TYPE(t_gradients)                           :: gradx,grady,gradz

      REAL, ALLOCATABLE :: div_temp(:, :)
      INTEGER :: i,ifftxc3

      ifftxc3=stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft

      ALLOCATE (div_temp(ifftxc3,1))

      CALL init_pw_grid(.TRUE.,stars,sym,cell)

      CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,bxc(1)%pw,gradx)
      CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,bxc(2)%pw,grady)
      CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,bxc(3)%pw,gradz)

      DO i = 1, ifftxc3
         div_temp(i,1)=gradx%gr(1,i,1)+grady%gr(2,i,1)+gradz%gr(3,i,1)
      ENDDO ! i
       
      CALL pw_from_grid(.TRUE.,stars,.TRUE.,div_temp,div%pw,div%pw_w)

      CALL finish_pw_grid()
   
   END SUBROUTINE pw_div

   SUBROUTINE divergence(stars,atoms,sphhar,vacuum,sym,cell,noco,bxc,div)

      !--------------------------------------------------------------------------
      ! Use the two divergence subroutines above to now put the complete diver-  
      ! gence of a field into a t_potden variable.                                         
      !--------------------------------------------------------------------------   

      IMPLICIT NONE

      TYPE(t_stars),                INTENT(IN)    :: stars
      TYPE(t_atoms),                INTENT(IN)    :: atoms
      TYPE(t_sphhar),               INTENT(IN)    :: sphhar
      TYPE(t_vacuum),               INTENT(IN)    :: vacuum
      TYPE(t_sym),                  INTENT(IN)    :: sym
      TYPE(t_cell),                 INTENT(IN)    :: cell
      TYPE(t_noco),                 INTENT(IN)    :: noco
      TYPE(t_potden), DIMENSION(3), INTENT(IN)    :: bxc
      TYPE(t_potden),               INTENT(INOUT) :: div

      INTEGER                                     :: n

      DO n=1,atoms%ntype
         CALL mt_div(n,atoms,sphhar,sym,bxc,div)
      END DO

      CALL pw_div(stars,sym,cell,noco,bxc,div)
      
   END SUBROUTINE divergence

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
   USE m_constants
   USE m_mt_tofrom_grid

   IMPLICIT NONE
   
   INTEGER, INTENT(IN)                         :: n
   TYPE(t_atoms), INTENT(IN)                   :: atoms
   TYPE(t_sphhar), INTENT(IN)                  :: sphhar
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: gradphi
   TYPE(t_potden), INTENT(IN)                  :: den
   TYPE(t_gradients)                           :: grad

   REAL, ALLOCATABLE :: grad_temp(:, :, :)
   REAL, ALLOCATABLE :: thet(:), phi(:)
   REAL :: r,th,ph
   INTEGER :: i, jr, k, nsp, kt

   nsp = atoms%nsp()

   ALLOCATE (grad%gr(3,atoms%jri(n)*nsp,1))
   ALLOCATE (grad_temp(atoms%jri(n)*nsp,1,3))
   ALLOCATE (thet(atoms%nsp()),phi(atoms%nsp()))

   CALL init_mt_grid(1, atoms, sphhar, .TRUE., sym, thet, phi)

   CALL mt_to_grid(.TRUE., 1, atoms, sphhar, den%mt(:,0:,n,:), n, grad)

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
   END DO
   
   CALL finish_mt_grid
   
   END SUBROUTINE mt_grad

   SUBROUTINE pw_grad(ifftxc3,jspins,stars,cell,noco,sym,den,gradphi)
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
   
   INTEGER, INTENT(IN)                         :: jspins, ifftxc3
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_noco), INTENT(IN)                    :: noco
   TYPE(t_stars),INTENT(IN)                    :: stars
   TYPE(t_cell),INTENT(IN)                     :: cell
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: gradphi
   TYPE(t_potden), INTENT(IN)                  :: den
   TYPE(t_gradients)                           :: grad

   REAL, ALLOCATABLE :: grad_temp(:, :, :)
   INTEGER :: i

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


   SUBROUTINE divpotgrad(jspins,n,ifftxc3,atoms,sphhar,sym,stars,cell,vacuum,noco,pot,grad)

   USE m_types
   USE m_constants
   USE m_mt_tofrom_grid
   IMPLICIT NONE
   !-----------------------------------------------------------------------------
   !Use the two gradient subroutines above to now put the complete gradient      
   !of a potential into a t_potden variable.                                     
   !-----------------------------------------------------------------------------   
   INTEGER, INTENT(IN)                         :: jspins, n, ifftxc3
   TYPE(t_atoms), INTENT(IN)                   :: atoms
   TYPE(t_sphhar), INTENT(IN)                  :: sphhar
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_noco), INTENT(IN)                    :: noco
   TYPE(t_stars),INTENT(IN)                    :: stars
   TYPE(t_cell),INTENT(IN)                     :: cell
   TYPE(t_vacuum),INTENT(IN)                   :: vacuum
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: grad
   TYPE(t_potden), INTENT(IN)                  :: pot

   TYPE(t_gradients)                           :: dummygrad   
   TYPE(t_potden)                              :: den   
   REAL, ALLOCATABLE                           :: ch(:, :)
   REAL                                        :: r
   REAL, ALLOCATABLE                           :: r2Array(:)
   INTEGER :: jr, k, nsp, kt

   CALL den%init(stars,atoms,sphhar,vacuum,noco,1,POTDEN_TYPE_DEN)
   ALLOCATE(den%pw_w,mold=den%pw)

   nsp = atoms%nsp()
   ALLOCATE(r2Array(nsp*atoms%jri(n)))
   ALLOCATE(ch(nsp*atoms%jri(n),1))

   CALL init_mt_grid(1, atoms, sphhar, .TRUE., sym)
   CALL mt_to_grid(.TRUE., 1, atoms, sphhar, pot%mt(:,0:,n,:), n, dummygrad, ch)

   kt = 0
   DO jr = 1, atoms%jri(n)
      r=atoms%rmsh(jr, n)
      DO k = 1, nsp
         r2Array(kt+nsp) = r*r
      ENDDO ! k
      kt = kt+nsp
   ENDDO ! jr

   ch(:,1)=ch(:,1)*r2Array

   CALL mt_from_grid(atoms, sphhar, n, 1, ch, den%mt(:,0:,n,:))

   CALL finish_mt_grid

   den%pw    = pot%pw
   den%pw_w  = pot%pw_w
   den%vacz  = pot%vacz
   den%vacxy = pot%vacxy

   CALL mt_grad(n,atoms,sphhar,sym,pot,grad)
   CALL pw_grad(ifftxc3,jspins,stars,cell,noco,sym,pot,grad)

   END SUBROUTINE divpotgrad

END MODULE m_divergence
