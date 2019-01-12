!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_genNewNocoInp

CONTAINS

SUBROUTINE genNewNocoInp(input,atoms,noco,noco_new)

   USE m_juDFT
   USE m_types
   USE m_constants
   USE m_rwnoco

   IMPLICIT NONE

   TYPE(t_input),INTENT(IN)         :: input
   TYPE(t_atoms),INTENT(IN)         :: atoms
   TYPE(t_noco),INTENT(IN)          :: noco
   TYPE(t_noco),INTENT(INOUT)       :: noco_new

   INTEGER                          :: iAtom, iType
   REAL                             :: alphdiff

   IF (.NOT.noco%l_mperp) THEN
      CALL juDFT_error ("genNewNocoInp without noco%l_mperp" ,calledby ="genNewNocoInp")
   END IF
   iAtom = 1
   DO iType = 1, atoms%ntype
      IF (noco%l_ss) THEN
         alphdiff = 2.0*pi_const*(noco%qss(1)*atoms%taual(1,iAtom) + &
                                  noco%qss(2)*atoms%taual(2,iAtom) + &
                                  noco%qss(3)*atoms%taual(3,iAtom) )
         noco_new%alph(iType) = noco%alph(iType) - alphdiff
         DO WHILE (noco_new%alph(iType) > +pi_const)
            noco_new%alph(iType)= noco_new%alph(iType) - 2.0*pi_const
         END DO
         DO WHILE (noco_new%alph(iType) < -pi_const)
            noco_new%alph(iType)= noco_new%alph(iType) + 2.0*pi_const
         END DO
      ELSE
         noco_new%alph(iType) = noco%alph(iType)
      END IF
      iatom= iatom + atoms%neq(iType)
   END DO

   OPEN (24,file='nocoinp',form='formatted', status='unknown')
   REWIND (24)
   CALL rw_noco_write(atoms,noco_new, input)
   CLOSE (24)

END SUBROUTINE genNewNocoInp

END MODULE m_genNewNocoInp
