!--------------------------------------------------------------------------------  
! Copyright (c) 2019 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_divergence
   USE m_types
   PRIVATE
   PUBLIC :: mt_div, pw_div, divergence!, divpotgrad

CONTAINS
   SUBROUTINE mt_div(n,atoms,sphhar,sym,xcB,div)
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
   
   INTEGER, INTENT(IN)                         :: n
   TYPE(t_atoms), INTENT(IN)                   :: atoms
   TYPE(t_sphhar), INTENT(IN)                  :: sphhar
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_potden), dimension(3), INTENT(IN)    :: xcB
   TYPE(t_potden), INTENT(INOUT)               :: div

   TYPE(t_gradients)                           :: gradx, grady, gradz

   REAL, ALLOCATABLE :: div_temp(:, :)
   REAL, ALLOCATABLE :: thet(:), phi(:)
   REAL :: r,th,ph
   INTEGER :: jr, k, nsp, kt

   nsp = atoms%nsp()

   ALLOCATE (gradx%gr(3,atoms%jri(n)*nsp,1),grady%gr(3,atoms%jri(n)*nsp,1),gradz%gr(3,atoms%jri(n)*nsp,1))
   ALLOCATE (div_temp(atoms%jri(n)*nsp,1))
   ALLOCATE (thet(atoms%nsp()),phi(atoms%nsp()))

   CALL init_mt_grid(1, atoms, sphhar, .TRUE., sym, thet, phi)

   CALL mt_to_grid(.TRUE., 1, atoms, sphhar, xcB(1)%mt(:,0:,n,:), n, gradx)
   CALL mt_to_grid(.TRUE., 1, atoms, sphhar, xcB(2)%mt(:,0:,n,:), n, grady)
   CALL mt_to_grid(.TRUE., 1, atoms, sphhar, xcB(3)%mt(:,0:,n,:), n, gradz)

   kt = 0
   DO jr = 1, atoms%jri(n)
      r=atoms%rmsh(jr, n)
      DO k = 1, nsp
         th = thet(k)
         ph = phi(k)
         div_temp(kt+nsp,1) = (SIN(th)*COS(ph)*gradx%gr(1,kt+nsp,1) + SIN(th)*SIN(ph)*grady%gr(1,kt+nsp,1) + COS(th)*gradz%gr(1,kt+nsp,1))&
                             +(COS(th)*COS(ph)*gradx%gr(2,kt+nsp,1) + COS(th)*SIN(ph)*grady%gr(2,kt+nsp,1) - SIN(th)*gradz%gr(2,kt+nsp,1))/r&
                             -(SIN(ph)*gradx%gr(3,kt+nsp,1)         + COS(ph)*grady%gr(3,kt+nsp,1))/(r*SIN(th))
      ENDDO ! k
   
      kt = kt+nsp
   ENDDO ! jr
    
   CALL mt_from_grid(atoms, sphhar, n, 1, div_temp, div%mt(:,0:,n,:))
   
   CALL finish_mt_grid
   
   END SUBROUTINE mt_div

   SUBROUTINE pw_div(ifftxc3,jspins,stars,cell,noco,sym,xcB,div)
   !-----------------------------------------------------------------------------!
   !By use of the cartesian components of a field and its carthesian derivatives !
   !in the interstitial/vacuum region at each grid point:                        !
   !                                                                             !
   !Make the divergence of said field in real space and store it as a source     !
   !density, again represented by pw-coefficients in a potden.                   !
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
   TYPE(t_potden), dimension(3), INTENT(IN)    :: xcB
   TYPE(t_potden), INTENT(INOUT)               :: div

   TYPE(t_gradients)                           :: gradx, grady, gradz

   REAL, ALLOCATABLE :: div_temp(:, :)
   INTEGER :: i, nsp

   nsp = ifftxc3
!   CALL xcpot%alloc_gradients(ifftxc3,1,gradx)
!   CALL xcpot%alloc_gradients(ifftxc3,2,gradx)
!   CALL xcpot%alloc_gradients(ifftxc3,3,gradx)
!   ALLOCATE (gradx%gr(3,nsp,1))
!   ALLOCATE (grady%gr(3,nsp,1))
!   ALLOCATE (gradz%gr(3,nsp,1))
!   ALLOCATE (gradx%agrt(nsp),grady%agrt(nsp),gradz%agrt(nsp))
!   ALLOCATE (gradx%agru(nsp),grady%agru(nsp),gradz%agru(nsp))
!   ALLOCATE (gradx%agrd(nsp),grady%agrd(nsp),gradz%agrd(nsp))
!   ALLOCATE (gradx%gggrt(nsp),grady%gggrt(nsp),gradz%gggrt(nsp))
!   ALLOCATE (gradx%gggru(nsp),grady%gggru(nsp),gradz%gggru(nsp))
!   ALLOCATE (gradx%gggrd(nsp),grady%gggrd(nsp),gradz%gggrd(nsp))
!   ALLOCATE (gradx%gzgr(nsp),grady%gzgr(nsp),gradz%gzgr(nsp))
!   ALLOCATE (gradx%g2rt(nsp),grady%g2rt(nsp),gradz%g2rt(nsp))
!   ALLOCATE (gradx%g2ru(nsp),grady%g2ru(nsp),gradz%g2ru(nsp))
!   ALLOCATE (gradx%g2rd(nsp),grady%g2rd(nsp),gradz%g2rd(nsp))

   ALLOCATE (div_temp(ifftxc3,1))

   CALL init_pw_grid(.TRUE.,stars,sym,cell)

   CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,xcB(1)%pw,gradx)
   CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,xcB(2)%pw,grady)
   CALL pw_to_grid(.TRUE.,1,.FALSE.,stars,cell,xcB(3)%pw,gradz)

!   DEALLOCATE(gradx%agrt,gradx%agru,gradx%agrd,gradx%gggrt,gradx%gggru,gradx%gggrd,gradx%gzgr,gradx%g2rt,gradx%g2ru,gradx%g2rd)
!   DEALLOCATE(grady%agrt,grady%agru,grady%agrd,grady%gggrt,grady%gggru,grady%gggrd,grady%gzgr,grady%g2rt,grady%g2ru,grady%g2rd)
!   DEALLOCATE(gradz%agrt,gradz%agru,gradz%agrd,gradz%gggrt,gradz%gggru,gradz%gggrd,gradz%gzgr,gradz%g2rt,gradz%g2ru,gradz%g2rd)

   DO i = 1, ifftxc3
         div_temp(i,1)=gradx%gr(1,i,1)+grady%gr(2,i,1)+gradz%gr(3,i,1)
   ENDDO ! i
    
   CALL pw_from_grid(.TRUE.,stars,.TRUE.,div_temp,div%pw,div%pw_w)

   CALL finish_pw_grid()
   
   END SUBROUTINE pw_div

   SUBROUTINE divergence(jspins,n,ifftxc3,atoms,sphhar,sym,stars,cell,vacuum,noco,xcB,div)

   USE m_types
   USE m_constants
   IMPLICIT NONE
   !-----------------------------------------------------------------------------!
   !Use the two divergence subroutines above to now put the complete divergence  !
   !of a field into a t_potden variable.                                         !
   !-----------------------------------------------------------------------------!   
   INTEGER, INTENT(IN)                         :: jspins, n, ifftxc3
   TYPE(t_atoms), INTENT(IN)                   :: atoms
   TYPE(t_sphhar), INTENT(IN)                  :: sphhar
   TYPE(t_sym), INTENT(IN)                     :: sym
   TYPE(t_noco), INTENT(IN)                    :: noco
   TYPE(t_stars),INTENT(IN)                    :: stars
   TYPE(t_cell),INTENT(IN)                     :: cell
   TYPE(t_vacuum),INTENT(IN)                   :: vacuum
   TYPE(t_potden), dimension(3), INTENT(IN)    :: xcB
   TYPE(t_potden), INTENT(INOUT)                 :: div

   CALL mt_div(n,atoms,sphhar,sym,xcB,div)
   CALL pw_div(ifftxc3,jspins,stars,cell,noco,sym,xcB,div)
      
   END SUBROUTINE divergence

!   SUBROUTINE divpotgrad(jspins,n,ifftxc3,atoms,sphhar,sym,stars,cell,vacuum,noco,xcB,div)
! 
!   END SUBROUTINE divpotgrad

END MODULE m_divergence
