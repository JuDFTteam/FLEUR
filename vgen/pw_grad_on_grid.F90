!--------------------------------------------------------------------------------  
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_pw_grad_on_grid
   USE m_types
   PRIVATE
   REAL,PARAMETER:: d_15=1.e-15

CONTAINS
   SUBROUTINE pw_do_div(cell,stars,n,atoms,sphhar,sym,xcB,div)
   !-----------------------------------------------------------------------------!
   !By use of the cartesian components of a field and its cartesian derivatives  !
   !in the interstitial space and vacuum:                                        !
   !                                                                             !
   !Make the divergence of said field in real space and store it as a source     !
   !density, again represented by pw-coefficients in a potden.                   !
   !                                                                             !
   !Code by A. Neukirchen, September 2019                                        !
   !-----------------------------------------------------------------------------! 
   USE m_grdrsis
   USE m_pw_tofrom_grid

   IMPLICIT NONE
   
   CLASS(t_xcpot),INTENT(IN)     :: xcpot
   INTEGER,INTENT(IN)            :: jspins
   LOGICAL,INTENT(IN)            :: l_noco
   TYPE(t_stars),INTENT(IN)      :: stars
   TYPE(t_cell),INTENT(IN)       :: cell
   TYPE(t_potden), dimension(3), INTENT(INOUT) :: xcB
   TYPE(t_potden), INTENT(INOUT)               :: div

   TYPE(t_gradients)                           :: gradx, grady, gradz

   REAL :: r,th,ph
   INTEGER :: jr, k, nsp, kt

   ALLOCATE (gradx%gr(3,stars%kxc1_fft*stars%kxc2_fft*stars%kxc3_fft,jspins),grady%gr(3,atoms%jri(n)*nsp,jspins),gradz%gr(3,atoms%jri(n)*nsp,jspins))
   ALLOCATE (div_temp(atoms%jri(n)*nsp,jspins))

   CALL grdrsis(,cell,stars%kxc1_fft,stars%kxc2_fft,stars%kxc3_fft,gradx)
   CALL mt_do_grad(xcpot, jspins, n, atoms, sphhar, sym, xcB(2)%mt(:,0:,n,:), grady)
   CALL mt_do_grad(xcpot, jspins, n, atoms, sphhar, sym, xcB(3)%mt(:,0:,n,:), gradz)

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
   
   CALL pw_from_grid(xcpot,stars,l_pw_w,v_in,v_out_pw,v_out_pw_w) 
   CALL mt_from_grid(atoms, sphhar, n, jspins, div_temp, div%mt(:,0:,n,:))
   
   DEALLOCATE (ylh, wt, rx, thet, phi, ylht, ylhtt, ylhf, ylhff, ylhtf)
   
   END SUBROUTINE pw_do_div
END MODULE m_pw_grad_on_grid
