!--------------------------------------------------------------------------------  
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_divergence
   USE m_types
   PRIVATE
   PUBLIC :: mt_div, pw_div, divergence

CONTAINS
   SUBROUTINE mt_div(jspins,n,atoms,sphhar,sym,xcB,div)
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

   CALL mt_div(jspins,n,atoms,sphhar,sym,xcB,div)
   CALL pw_div(ifftxc3,jspins,stars,cell,noco,sym,xcB,div)
      
   END SUBROUTINE divergence

END MODULE m_divergence
