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
SUBROUTINE magnMomFromDen(input,atoms,noco,den)
   USE m_constants
   USE m_types
   USE m_constants
   USE m_intgr
   IMPLICIT NONE

!   TYPE(t_dimension), INTENT(IN) ::  dimension
   TYPE(t_input), INTENT(IN)     ::  input
   TYPE(t_atoms), INTENT(IN)     ::  atoms
   TYPE(t_noco), INTENT(IN)      ::  noco
   TYPE(t_potden),INTENT(IN)     ::  den

   INTEGER                       ::  jsp,lmax,i,l,j
   REAL, ALLOCATABLE             ::  dummyResults(:,:)
   REAL                          ::  moments(3)
   !!!Declare IntG1

  lmax=0

  ALLOCATE(dummyResults(SIZE(den%mt,3),SIZE(den%mt,4)))

  IF(noco%l_mtNocoPot) THEN
     jsp=4
  ELSE 
     jsp=input%jspins
  END IF

   DO i=1, atoms%ntype 
      DO j=1, jsp
!!!!Integration
           CALL intgr3(den%mt(:,0,i,j), atoms%rmsh(:,i),atoms%dx(i),atoms%jri(i),dummyResults(i,j))
!!!!Considering Lattice harmonics integral (Only L=0 component does not vanish and has a factor of sqrt(4*Pi))
          dummyResults(i,j)=dummyResults(i,j)*sfp_const
      END DO
   END DO 

!!Print Results
DO i=1 , atoms%ntype
   moments(3)=dummyResults(i,1)-dummyResults(i,2)
   moments(1:2)=2*dummyResults(i,3:4)
   write(*,*) "Magnetic Moment of Atom "
   write(*,*) i 
   write(*,*) " my=" 
   write(*,*) moments(2)!/ (4*pimach()*atoms%rmt**3) *3

   write(*,*) " mx=" 
   write(*,*) moments(1)!/(4*pimach()*atoms%rmt**3)*3

   write(*,*) " mz=" 
   write(*,*) moments(3)

END DO

!!!!Normalization?


 
 DEALLOCATE(dummyResults)


END SUBROUTINE magnMomFromDen
END MODULE m_magnMomFromDen
