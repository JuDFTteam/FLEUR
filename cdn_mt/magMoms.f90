!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_magMoms

CONTAINS

SUBROUTINE magMoms(dimension,input,atoms,noco,vTot,chmom,qa21)

   USE m_types
   USE m_xmlOutput
   USE m_m_perp

   IMPLICIT NONE

   TYPE(t_dimension), INTENT(IN) :: dimension
   TYPE(t_input), INTENT(IN)     :: input
   TYPE(t_atoms), INTENT(IN)     :: atoms
   TYPE(t_noco), INTENT(INOUT)   :: noco
   TYPE(t_potden),INTENT(IN)     :: vTot

   REAL, INTENT(INOUT)           :: chmom(atoms%ntype,dimension%jspd)
   COMPLEX, INTENT(IN)           :: qa21(atoms%ntype)

   INTEGER                       :: iType, j, iRepAtom
   REAL                          :: smom
   CHARACTER(LEN=20)             :: attributes(4)

   WRITE (6,FMT=8020)
   WRITE (16,FMT=8020)

   CALL openXMLElement('magneticMomentsInMTSpheres',(/'units'/),(/'muBohr'/))
   iRepAtom = 1
   DO iType = 1, atoms%ntype
      smom = chmom(iType,1) - chmom(iType,input%jspins)
      WRITE (6,FMT=8030) iType,smom, (chmom(iType,j),j=1,input%jspins)
      WRITE (16,FMT=8030) iType,smom, (chmom(iType,j),j=1,input%jspins)
      attributes = ''
      WRITE(attributes(1),'(i0)') iType
      WRITE(attributes(2),'(f15.10)') smom
      WRITE(attributes(3),'(f15.10)') chmom(iType,1)
      WRITE(attributes(4),'(f15.10)') chmom(iType,2)
      CALL writeXMLElementFormPoly('magneticMoment',(/'atomType      ','moment        ','spinUpCharge  ',&
                                                      'spinDownCharge'/),&
                                   attributes,reshape((/8,6,12,14,6,15,15,15/),(/4,2/)))

      IF (noco%l_mperp) THEN
         !calculate the perpendicular part of the local moment
         !and relax the angle of the local moment or calculate
         !the constraint B-field.
         CALL m_perp(atoms,iType,iRepAtom,noco,vTot%mt(:,0,:,:),chmom,qa21)
      END IF
      iRepAtom= iRepAtom + atoms%neq(iType)
   END DO
   CALL closeXMLElement('magneticMomentsInMTSpheres')

   8020 FORMAT (/,/,2x,'-->  magnetic moments in the spheres:',/,2x,&
                'mm -->   type',t22,'moment',t33,'spin-up',t43,'spin-down')
   8030 FORMAT (2x,'--> mm',i8,2x,3f12.5)

END SUBROUTINE magMoms

END MODULE m_magMoms
