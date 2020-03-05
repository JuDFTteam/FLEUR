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
! Modified for usability with the potential matrix; A. Neukirchen, Dec '19

MODULE m_magnMomFromDen

IMPLICIT NONE

CONTAINS
SUBROUTINE magnMomFromDen(input,atoms,noco,den,moments,theta_mt_avg,phi_mt_avg)
   USE m_constants
   USE m_types
   USE m_types_fleurinput
   USE m_intgr
   USE m_juDFT
   USE m_polangle
   USE m_constants

   TYPE(t_input), INTENT(IN)     ::  input
   TYPE(t_atoms), INTENT(IN)     ::  atoms
   TYPE(t_noco), INTENT(IN)      ::  noco
   TYPE(t_potden),INTENT(IN)     ::  den
   REAL, INTENT(OUT)             ::  moments(3,atoms%ntype)
   REAL,INTENT(OUT)              :: theta_mt_avg(atoms%ntype)
   REAL,INTENT(OUT)              :: phi_mt_avg(atoms%ntype)

   INTEGER                       ::  jsp,i,j,ir
   REAL                          ::  mx,my,mz

   TYPE(t_potden)                ::  denloc
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
           IF (den%potdenType<=1000) THEN
              CALL denloc%copyPotDen(den)
              DO ir=1, atoms%jri(i)
                 denloc%mt(ir,:,i,j)=den%mt(ir,:,i,j)*atoms%rmsh(ir,i)**2
              END DO
           ELSE
              CALL denloc%copyPotDen(den)
           END IF
           CALL intgr3(denloc%mt(:,0,i,j), atoms%rmsh(:,i),atoms%dx(i),atoms%jri(i),dummyResults(i,j))
!!Considering Lattice harmonics integral (Only L=0 component does not vanish and has a factor of sqrt(4*Pi))
          dummyResults(i,j)=dummyResults(i,j)*sfp_const
      END DO
   END DO
!!Assign results
   DO i=1 , atoms%ntype
   IF (noco%l_mtNocoPot) THEN
      moments(1:2,i)=2*dummyResults(i,3:4)
   END IF
      moments(3,i)=dummyResults(i,1)-dummyResults(i,2)
   END DO
   
   

DEALLOCATE(dummyResults)

IF (den%potdenType<=1000) THEN
   moments=moments/2
ELSE
   moments(2,:)=-moments(2,:)
END IF

!!Calculation of Angles
   DO i=1 , atoms%ntype
      mx=moments(1,i)
      my=moments(2,i)
      mz=moments(3,i)
      IF (den%potdenType>1000) THEN
         CALL pol_angle(mx,my,mz,theta_mt_avg(i),phi_mt_avg(i))
      END IF
   ENDDO

END SUBROUTINE magnMomFromDen
END MODULE m_magnMomFromDen
