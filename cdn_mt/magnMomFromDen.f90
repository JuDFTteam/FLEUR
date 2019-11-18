!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and avhttps://gcc.gnu.org/onlinedocs/gfortran/SQRT.htmlailable as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!------------------------------------------------------------------------------
! This routine calculates the magnetic moments and the corresponding directions 
! (angles) according to the Atoms in the system. 
! 
!
! Robin Hilgers, Nov '19

MODULE m_magnMomFromDen
CONTAINS
SUBROUTINE magnMomFromDen(input,atoms,noco,den,moments)
   USE m_constants
   USE m_types
   USE m_constants
   USE m_intgr
   USE m_juDFT
   IMPLICIT NONE


   TYPE(t_input), INTENT(IN)     ::  input
   TYPE(t_atoms), INTENT(INOUT)  ::  atoms
   TYPE(t_noco), INTENT(IN)      ::  noco
   TYPE(t_potden),INTENT(IN)     ::  den
   REAL, INTENT(OUT)             ::  moments(atoms%ntype,3)

   INTEGER                       ::  jsp,i,j
   REAL                          ::  phi,theta, mx,my,mz
   REAL                          :: eps=1E-10

   REAL, ALLOCATABLE             ::  dummyResults(:,:)
  IF(.NOT.ALLOCATED(atoms%phi_mt_avg).AND..NOT.ALLOCATED(atoms%theta_mt_avg)) THEN
     ALLOCATE(atoms%phi_mt_avg(atoms%ntype),atoms%theta_mt_avg(atoms%ntype))
  ELSE IF (ALLOCATED(atoms%phi_mt_avg).AND.ALLOCATED(atoms%theta_mt_avg))THEN

  ELSE
     CALL judft_error("Allocation error", calledby="magnMomFromDen.f90")
  END IF

  ALLOCATE(dummyResults(SIZE(den%mt,3),SIZE(den%mt,4)))

  IF(noco%l_mtNocoPot) THEN
     jsp=4
  ELSE 
     jsp=input%jspins
  END IF
!!Loop over Spins and Atoms
   DO i=1, atoms%ntype 
      DO j=1, jsp
!!Integration over r
           CALL intgr3(den%mt(:,0,i,j), atoms%rmsh(:,i),atoms%dx(i),atoms%jri(i),dummyResults(i,j))
!!Considering Lattice harmonics integral (Only L=0 component does not vanish and has a factor of sqrt(4*Pi))
          dummyResults(i,j)=dummyResults(i,j)*sfp_const
      END DO
   END DO 
!!Assign results
   DO i=1 , atoms%ntype
   IF (noco%l_mtNocoPot) THEN
      moments(i,1:2)=2*dummyResults(i,3:4)
   END IF
      moments(i,3)=dummyResults(i,1)-dummyResults(i,2)
   END DO
DEALLOCATE(dummyResults)

!!Calculation of Angles
   DO i=1 , atoms%ntype
      mx=moments(i,1)
      my=moments(i,2)
      mz=moments(i,3)
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
     atoms%phi_mt_avg(i)=phi
     atoms%theta_mt_avg(i)=theta
   ENDDO




END SUBROUTINE magnMomFromDen
END MODULE m_magnMomFromDen
