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
      SUBROUTINE inpnoco(atoms,input,sym,vacuum,noco)

      USE m_constants
      USE m_rwnoco
      USE m_types_atoms
      USE m_types_input
      USE m_types_sym
      USE m_types_vacuum
      USE m_types_noco
      Use m_nocoInputCheck
      IMPLICIT NONE
      TYPE(t_atoms),INTENT(INOUT) ::atoms
      TYPE(t_input),INTENT(INOUT) ::input
      TYPE(t_sym),INTENT(IN)      ::sym
      TYPE(t_vacuum),INTENT(IN)   ::vacuum
      TYPE(t_noco),INTENT(INOUT)  ::noco

!     ..
!     .. Local Scalars ..
      INTEGER itype,iatom

      OPEN (24,file='nocoinp',form='formatted',status='old')

      WRITE (oUnit,*)'This is a non-collinear calculation. The magnetic'
      WRITE (oUnit,*)'moments of the atoms have a fixed direction.'
      WRITE (oUnit,*)'The Euler-angles alpha and beta of this direction'
      WRITE (oUnit,*)'are equal to the polar angles of the magnetic'
      WRITE (oUnit,*)'moment vector phi and theta respectively.'
      WRITE (oUnit,*)

      CALL rw_noco_read(atoms,noco,input)

      CALL nocoInputCheck(atoms,input,sym,vacuum,noco)

      IF (noco%l_ss) THEN
!
!--->    the angle beta is relative to the spiral in a spin-spiral
!--->    calculation, i.e. if beta = 0 for all atoms in the unit cell
!--->    that means that the moments are "in line" with the spin-spiral
!--->    (beta = qss * taual). note: this means that only atoms within
!--->    a plane perpendicular to qss can be equivalent!
         DO itype = 1,atoms%ntype
            iatom = atoms%firstAtom(itype)
            noco%phi_inp = tpi_const*dot_product(noco%qss_inp,atoms%taual(:,iatom))
            noco%alph_inp(itype) = noco%alph_inp(itype) + noco%phi_inp
         ENDDO
      ENDIF

      WRITE (oUnit,*)
      CLOSE (24)

      END SUBROUTINE inpnoco
      END MODULE m_inpnoco
