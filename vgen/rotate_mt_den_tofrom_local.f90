!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_rotate_mt_den_tofrom_local
  USE m_juDFT
  USE m_types
  USE m_constants
  use m_mt_tofrom_grid
  IMPLICIT NONE
CONTAINS
  SUBROUTINE rotate_mt_den_to_local(atoms,sphhar,sym,den)
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_sphhar),INTENT(IN) :: sphhar
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_potden),INTENT(INOUT) :: den
    
    
    TYPE(t_xcpot_inbuild)    :: xcpot !local xcpot that is LDA to indicate we do not need gradients
    TYPE(t_gradients)        :: grad

    INTEGER :: n,nsp,imesh,i,b_ind
    REAL    :: rho_11,rho_22,rho_21r,rho_21i,mx,my,mz,magmom
    REAL    :: rhotot,rho_up,rho_down,theta,phi
    REAL,ALLOCATABLE :: ch(:,:)
    REAL    :: eps=1E-10
    nsp=atoms%nsp()
    ALLOCATE(ch(nsp*atoms%jmtd,4),&
             den%theta_mt(nsp*atoms%jmtd,atoms%ntype),&
             den%phi_mt(nsp*atoms%jmtd,atoms%ntype))
    CALL xcpot%init("vwn",.FALSE.,1)

    CALL init_mt_grid(4,atoms,sphhar,xcpot%needs_grad(),sym)
    DO n=1,atoms%ntype
       CALL mt_to_grid(xcpot%needs_grad(),4,atoms,sphhar,den%mt(:,0:,n,:),n,grad,ch)
       DO imesh = 1,nsp*atoms%jri(n)
    
          rho_11  = ch(imesh,1)
          rho_22  = ch(imesh,2)
          rho_21r = ch(imesh,3)
          rho_21i = ch(imesh,4)
          mx      =  2*rho_21r
          my      = -2*rho_21i
          mz      = (rho_11-rho_22)
          magmom  = SQRT(mx**2 + my**2 + mz**2)
          rhotot  = rho_11 + rho_22
          rho_up  = (rhotot + magmom)/2
          rho_down= (rhotot - magmom)/2

          IF (ABS(mz) .LE. eps) THEN
             theta = pi_const/2
          ELSEIF (mz .GE. 0.0) THEN
             theta = ATAN(SQRT(mx**2 + my**2)/mz)
          ELSE
             theta = ATAN(SQRT(mx**2 + my**2)/mz) + pi_const
          ENDIF

          IF (ABS(mx) .LE. eps) THEN
             IF (ABS(my) .LE. eps) THEN
                phi = 0.0
             ELSEIF (my .GE. 0.0) THEN
                phi = pi_const/2
             ELSE
                phi = -pi_const/2
             ENDIF
          ELSEIF (mx .GE. 0.0) THEN
             phi = ATAN(my/mx)
          ELSE
             IF (my .GE. 0.0) THEN
                phi = ATAN(my/mx) + pi_const
             ELSE
                phi = ATAN(my/mx) - pi_const
             ENDIF
          ENDIF
          
          ch(imesh,1) = rho_up
          ch(imesh,2) = rho_down
          den%theta_mt(imesh,n) = theta
          den%phi_mt(imesh,n) = phi
       ENDDO
       den%mt(:,0:,n,:)=0.0
       CALL mt_from_grid(atoms,sphhar,n,2,ch,den%mt(:,0:,n,:))
       DO i=1,atoms%jri(n)
          den%mt(i,:,n,:)=den%mt(i,:,n,:)*atoms%rmsh(i,n)**2
       ENDDO
    END DO
    CALL finish_mt_grid()
  END SUBROUTINE rotate_mt_den_to_local

  SUBROUTINE rotate_mt_den_from_local(atoms,sphhar,sym,den,vtot)
    TYPE(t_atoms),INTENT(IN)  :: atoms
    TYPE(t_sphhar),INTENT(IN) :: sphhar
    TYPE(t_sym),INTENT(IN)    :: sym
    TYPE(t_potden),INTENT(IN) :: den
    TYPE(t_potden),INTENT(INOUT) :: vtot
    
    TYPE(t_xcpot_inbuild)     :: xcpot !local xcpot that is LDA to indicate we do not need gradients
    TYPE(t_gradients) :: grad
    
    INTEGER :: n,nsp,imesh,i
    REAL    :: vup,vdown,veff,beff
    REAL    :: theta,phi
    REAL,ALLOCATABLE :: ch(:,:)
    
    nsp=atoms%nsp()
    ALLOCATE(ch(nsp*atoms%jmtd,4))
    CALL xcpot%init("vwn",.FALSE.,1)

    CALL init_mt_grid(4,atoms,sphhar,xcpot%needs_grad(),sym)
    DO n=1,atoms%ntype
       CALL mt_to_grid(xcpot%needs_grad(),4,atoms,sphhar,vtot%mt(:,0:,n,:),n,grad,ch)
       DO imesh = 1,nsp*atoms%jri(n)
          vup   = ch(imesh,1)
          vdown = ch(imesh,2)
          theta = den%theta_mt(imesh,n)
          phi   = den%phi_mt(imesh,n)
          veff  = (vup + vdown)/2.0
          beff  = (vup - vdown)/2.0
          ch(imesh,1) = veff + beff*COS(theta)
          ch(imesh,2) = veff - beff*COS(theta)
          ch(imesh,3) = beff*SIN(theta)*COS(phi)
          ch(imesh,4) = beff*SIN(theta)*SIN(phi)
       ENDDO

       vtot%mt(:,0:,n,:)=0.0

       CALL mt_from_grid(atoms,sphhar,n,4,ch,vtot%mt(:,0:,n,:))

       DO i=1,atoms%jri(n)
          vtot%mt(i,:,n,:)=vtot%mt(i,:,n,:)*atoms%rmsh(i,n)**2
       ENDDO
    END DO
    CALL finish_mt_grid()
  END SUBROUTINE rotate_mt_den_from_local

END MODULE m_rotate_mt_den_tofrom_local
