!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

      MODULE m_inpnoco

!**********************************************************************
!     This subroutine reads the  Euler angles of the  magnetic field 
!     directions and other noncollinear parameters from the file nocoinp
!
!                                                 Philipp Kurz 98/01/28
!******** ABBREVIATIONS ***********************************************
!     alpha,beta:Euler angles of the local magnetic field direction of
!                each atom(-type). 
!**********************************************************************
      CONTAINS
      SUBROUTINE inpnoco(atoms,input,vacuum,jij,noco)

      USE m_constants, ONLY : tpi_const
      USE m_rwnoco
      USE m_types
      USE m_nocoInputCheck
      IMPLICIT NONE
      TYPE(t_atoms),INTENT(INOUT) ::atoms
      TYPE(t_input),INTENT(INOUT) ::input
      TYPE(t_vacuum),INTENT(IN)   ::vacuum
      TYPE(t_Jij),INTENT(INOUT)   ::Jij
      TYPE(t_noco),INTENT(INOUT)  ::noco

!     ..
!     .. Local Scalars ..
      INTEGER itype,iatom

      OPEN (24,file='nocoinp',form='formatted',status='old')

      WRITE (6,*)'This is a non-collinear calculation. The magnetic'
      WRITE (6,*)'moments of the atoms have a fixed direction.'
      WRITE (6,*)'The Euler-angles alpha and beta of this direction'
      WRITE (6,*)'are equal to the polar angles of the magnetic'
      WRITE (6,*)'moment vector phi and theta respectively.'
      WRITE (6,*)

      CALL rw_noco_read(atoms,jij,noco,input)

      CALL nocoInputCheck(atoms,input,vacuum,jij,noco)
         
      IF (.not.jij%l_j.and.noco%l_ss) THEN
!
!--->    the angle beta is relative to the spiral in a spin-spiral
!--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
!--->    that means that the moments are "in line" with the spin-spiral
!--->    (beta = qss * taual). note: this means that only atoms within
!--->    a plane perpendicular to qss can be equivalent!
         iatom = 1
         DO itype = 1,atoms%ntype
            noco%phi = tpi_const*dot_product(noco%qss,atoms%taual(:,iatom))
            noco%alph(itype) = noco%alph(itype) + noco%phi
            iatom = iatom + atoms%neq(itype)
         ENDDO
      ENDIF

      WRITE (6,*)
      CLOSE (24)

      END SUBROUTINE inpnoco
      END MODULE m_inpnoco
