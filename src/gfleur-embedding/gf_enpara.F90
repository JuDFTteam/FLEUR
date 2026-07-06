!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_enpara
      USE m_constants, ONLY: oUnit
!----------------------------------------------------------------
!In the port to current FLEUR the energy parameters come from each
!layer's inp.xml (enpara%init_enpara + enpara%update). This module only
!keeps the optional "enpara_atoms" override: a text file with one line
!per element (zatom  s p d f quantum numbers + 4 lo quantum numbers)
!that is applied across all layers - useful to enforce consistent
!energy parameters in multi-layer setups. The quantum numbers are
!written into enpara%qn_el/qn_ello; the energies are then derived from
!the potential by the standard enpara%update.
!----------------------------------------------------------------
      use m_juDFT
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE gf_apply_enpara_atoms(atoms,jspins,enpara)
!*****************************************************************
!     DESC:applies the enpara_atoms override for one layer
!     Daniel Wortmann, Wed Aug 28 10:40:06 2002
!*****************************************************************
      USE m_gf_types
      IMPLICIT NONE
      !<--Arguments
      TYPE(t_atoms),INTENT(IN)     :: atoms
      INTEGER,INTENT(IN)           :: jspins
      TYPE(t_enpara),INTENT(INOUT) :: enpara
      !>
      LOGICAL              :: l_exist
      CHARACTER(LEN = 200) :: line
      INTEGER              :: l,n,i,zatom,lo
      INTEGER              :: ello_in(4),el_in(4)
      INTEGER              :: ello(4,112),el(4,112)

      INQUIRE(FILE ="enpara_atoms",EXIST = l_exist)
      IF (.NOT.l_exist) RETURN

      el = 0; ello = 0
      !<-- read all lines in enpara_atoms
      OPEN(40,FILE="enpara_atoms",STATUS ="old")
      DO
         READ(40,"(a)",END = 99) line
         IF (line(1:1) =="#") CYCLE
         READ(line,*,ERR = 999) zatom,el_in,ello_in
         el(:,zatom) = el_in
         ello(:,zatom) = ello_in
      ENDDO
   99 CLOSE(40)
      !>
      !<--assign enparas to atoms
      WRITE(oUnit,*)
      WRITE(oUnit,*) "Energy-parameter quantum numbers from enpara_atoms"
      WRITE(oUnit,*) "Atom   s  p  d  f  lo's"
      DO n=1,atoms%ntype
         IF (ALL(el(:,atoms%nz(n))==0)) CYCLE
         DO l=0,3
            enpara%qn_el(l,n,:jspins) = el(l+1,atoms%nz(n))
         ENDDO
         DO lo=1,atoms%nlo(n)
            enpara%qn_ello(lo,n,:jspins) =                              &
     &           ello(atoms%llo(lo,n)+1,atoms%nz(n))
         ENDDO
         WRITE (oUnit,"(i5,1x,4(i2,1x))") n,(enpara%qn_el(i,n,1),i=0,3)
         IF (atoms%nlo(n) >= 1) THEN
            WRITE (oUnit,"(7x,4(i2,1x))")                               &
     &           (enpara%qn_ello(lo,n,1),lo = 1,atoms%nlo(n))
         ENDIF
      ENDDO
      !>
      RETURN
  999 WRITE(*,*) "Error reading enpara_atoms:"
      WRITE(*,*) line
      CALL juDFT_error("Problem in enpara_atoms")
      END SUBROUTINE
      END
