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
SUBROUTINE magnMomFromDen(dimension,input,atoms,noco,den,moments)
   USE m_constants
   USE m_types
   USE m_intgr
   IMPLICIT NONE

   TYPE(t_dimension), INTENT(IN) ::  dimension
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
         !DO l=1, lmax
            CALL intgr3(den%mt(:,0,i,j), atoms%rmsh(1,i),atoms%dx(i),atoms%jri(i),dummyResults(i,j))
         !END DO
!!!!Considering Lattice harmonics integral (Only L=0 component does not vanish and has a factor of sqrt(4*Pi))
        dummyResults(:,:)=SQRT(4*pimach())*dummyResults(:,:) 
      END DO
   END DO 
!!Print Results
DO i=1 , atoms%ntype
   moments(1)=dummyResults(i,1)-dummyResults(i,2)
   moments(2:3)=dummyResults(i,3:4)
   write(*,*) "Magnetic Moment of Atom "
   write(*,*) i 
   write(*,*) " mx=" 
   write(*,*) moments(1) 
   write(*,*) " my=" 
   write(*,*) moments(2) 
   write(*,*) " mz=" 
   write(*,*) moments(3) 		
END DO

!!!!Normalization?


 
 DEALLOCATE(dummyResults)


END SUBROUTINE magnMomFromDen
END MODULE m_magnMomFromDen
