!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_checkdopall

CONTAINS

SUBROUTINE checkDOPAll(input,dimension,sphhar,stars,atoms,sym,vacuum,oneD,&
                       cell,potden,ispin)

   USE m_sphpts
   USE m_checkdop
   USE m_types
   USE m_cylpts
   USE m_points
   USE m_juDFT

   IMPLICIT NONE

   TYPE(t_input),INTENT(IN)     :: input
   TYPE(t_dimension),INTENT(IN) :: dimension
   TYPE(t_sphhar),intent(in)    :: sphhar      
   TYPE(t_stars),INTENT(IN)     :: stars
   TYPE(t_atoms),INTENT(IN)     :: atoms
   TYPE(t_sym),INTENT(IN)       :: sym
   TYPE(t_vacuum),INTENT(IN)    :: vacuum
   TYPE(t_oneD),INTENT(IN)      :: oneD
   TYPE(t_cell),INTENT(IN)      :: cell
   TYPE(t_potden),INTENT(IN)    :: potden

   INTEGER, INTENT(IN)          :: ispin

   INTEGER                      :: npd, nat, n, ivac
   REAL                         :: signum

   REAL                         :: xp(3,dimension%nspd)

   CALL timestart("checkDOPAll")

   IF ((input%film).AND.(.NOT.oneD%odi%d1)) THEN
      !--->             vacuum boundaries
      npd = min(dimension%nspd,25)
      CALL points(xp,npd)
      DO ivac = 1,vacuum%nvac
         signum = 3.0 - 2.0*ivac
         xp(3,:npd) = signum*cell%z1/cell%amat(3,3)
         CALL checkdop(xp,npd,0,0,ivac,1,ispin,dimension,atoms,&
                       sphhar,stars,sym,vacuum,cell,oneD,potden)
      END DO
   ELSE IF (oneD%odi%d1) THEN
      npd = min(dimension%nspd,25)
      CALL cylpts(xp,npd,cell%z1)
      CALL checkdop(xp,npd,0,0,ivac,1,ispin,dimension,atoms,&
                    sphhar,stars,sym,vacuum,cell,oneD,potden)
   END IF

   !--->          m.t. boundaries
   nat = 1
   DO n = 1, atoms%ntype
      CALL sphpts(xp,dimension%nspd,atoms%rmt(n),atoms%pos(1,atoms%nat))
      CALL checkdop(xp,dimension%nspd,n,nat,0,-1,ispin,&
                    dimension,atoms,sphhar,stars,sym,vacuum,cell,oneD,potden)
      nat = nat + atoms%neq(n)
   END DO

   CALL timestop("checkDOPAll")

END SUBROUTINE checkDOPAll

END MODULE m_checkdopall
