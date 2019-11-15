!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
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

   USE m_types
   USE m_intgr
   IMPLICIT NONE

   TYPE(t_dimension), INTENT(IN) :: dimension
   TYPE(t_input), INTENT(IN)     :: input
   TYPE(t_atoms), INTENT(IN)     :: atoms
   TYPE(t_noco), INTENT(IN)   :: noco
   TYPE(t_potden),INTENT(IN)     :: den

   INTEGER, ALLOCATABLE :: moments(:,:)
   INTEGER              :: jsp,lmax,i,l,j
   REAL, ALLOCATABLE    :: Dummy(:), DummyResults(:,:,:)
   REAL                 :: Result
   !!!Declare IntG1

  lmax=SIZE(den%mt,2)

  ALLOCATE(Dummy(SIZE(den%mt,1)))
  ALLOCATE(DummyResults(SIZE(den%mt,2),SIZE(den%mt,3),SIZE(den%mt,4)))

  IF(noco%l_mtNocoPot) THEN
     jsp=4
  ELSE 
     jsp=input%jspins
  END IF

   DO i=1, atoms%ntype 
      DO j=1, jsp
!!!!Integration
         DO l=1, lmax
            Dummy(:)=den%mt(:,l,i,j)
            CALL intgr3( Dummy(:), atoms%rmsh(1,i),atoms%dx(i),atoms%jri(i),Result)
            DummyResults(l,i,j)=Result
         END DO
!!!!SUMMATION considering Lattice harmonics 
      END DO
   END DO 

!!!!Normalization?




END SUBROUTINE magnMomFromDen
END MODULE m_magnMomFromDen
