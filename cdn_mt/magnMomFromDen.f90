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

IMPLICIT NONE

CONTAINS

SUBROUTINE magnMomFromDen(input,atoms,noco,den,moments)
   
   USE m_constants
   USE m_types
   USE m_intgr
   USE m_juDFT
   USE m_polangle



   TYPE(t_input), INTENT(IN)     ::  input
   TYPE(t_atoms), INTENT(INOUT)  ::  atoms
   TYPE(t_noco), INTENT(IN)      ::  noco
   TYPE(t_potden),INTENT(IN)     ::  den
   REAL, INTENT(OUT)             ::  moments(3,atoms%ntype)

   INTEGER                       ::  jsp,i,j
   REAL                          ::  mx,my,mz

   REAL, ALLOCATABLE             ::  dummyResults(:,:)

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
      moments(1:2,i)=2*dummyResults(i,3:4)
   END IF
      moments(i,3)=dummyResults(i,1)-dummyResults(i,2)
   END DO
DEALLOCATE(dummyResults)

!!Calculation of Angles
   DO i=1 , atoms%ntype
      mx=moments(1,i)
      my=moments(2,i)
      mz=moments(3,i)
      CALL pol_angle(mx,my,mz,atoms%theta_mt_avg(i),atoms%phi_mt_avg(i))
      IF(mx<0) atoms%theta_mt_avg(i)=-atoms%theta_mt_avg(i)
   ENDDO

END SUBROUTINE magnMomFromDen
END MODULE m_magnMomFromDen
