!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_orbMagMoms

CONTAINS

SUBROUTINE orbMagMoms(input,atoms,noco,nococonv,clmom)

   USE m_types
   USE m_constants
   USE m_xmlOutput

   IMPLICIT NONE

   TYPE(t_input), INTENT(IN) :: input
   TYPE(t_atoms), INTENT(IN)     :: atoms
   TYPE(t_noco), INTENT(IN)      :: noco
   TYPE(t_nococonv), INTENT(IN)      :: nococonv


   REAL, INTENT(INOUT)           :: clmom(3,atoms%ntype,input%jspins)

   INTEGER                       :: iType, j
   REAL                          :: thetai, phii, slmom, slxmom, slymom,szglobal,syglobal,sxglobal
   CHARACTER(LEN=20)             :: attributes(4)


   thetai = nococonv%theta
   phii   = nococonv%phi
   WRITE (oUnit,FMT=9020)
   CALL openXMLElement('orbitalMagneticMomentsInMTSpheres',(/'units'/),(/'muBohr'/))
   DO iType = 1, atoms%ntype
      ! magn. moment(-)
      slxmom = clmom(1,iType,1)+clmom(1,iType,2)
      slymom = clmom(2,iType,1)+clmom(2,iType,2)
      slmom =  clmom(3,iType,1)+clmom(3,iType,2)

      IF (any(noco%l_unrestrictMT)) THEN
         szglobal=clmom(3,iType,1)+clmom(3,iType,2)
         sxglobal=clmom(1,iType,1)+clmom(1,iType,2)
         syglobal=clmom(2,iType,1)+clmom(2,iType,2)
      END IF

      IF (noco%l_noco) THEN
         thetai = nococonv%beta(iType)
         phii   = nococonv%alph(iType)

         !Fix of sign of moment in first variation calculations. Perhaps it would be better to understand this :-(
         !slxmom=-1*slxmom
         !slymom=-1*slymom
         !slmom=-1*slmom
      END IF
      ! rotation: orbital moment || spin moment (extended to incude phi - hopefully)
      slmom   = cos(thetai)*slmom + sin(thetai)*(cos(phii)*slxmom + sin(phii)*slymom)
      clmom(3,iType,1) = cos(thetai)*clmom(3,iType,1) + &
                         sin(thetai)*(cos(phii)*clmom(1,iType,1) + &
                         sin(phii)*clmom(2,iType,1))
      clmom(3,iType,2) = cos(thetai)*clmom(3,iType,2) + &
                         sin(thetai)*(cos(phii)*clmom(1,iType,2) + &
                         sin(phii)*clmom(2,iType,2))

      WRITE (oUnit,FMT=8030) iType,slmom,(clmom(3,iType,j),j=1,2)
      attributes = ''
      WRITE(attributes(1),'(i0)') iType
      WRITE(attributes(2),'(f15.10)') slmom
      WRITE(attributes(3),'(f15.10)') clmom(3,iType,1)
      WRITE(attributes(4),'(f15.10)') clmom(3,iType,2)
      CALL writeXMLElementFormPoly('orbMagMoment',(/'atomType      ','moment        ','spinUpCharge  ',&
                                                    'spinDownCharge'/),&
                                   attributes,reshape((/8,6,12,14,6,15,15,15/),(/4,2/)))
      IF (any(noco%l_unrestrictMT)) THEN
         WRITE(oUnit,FMT=8032) iType, sxglobal,syglobal,szglobal,sqrt(sxglobal**2+syglobal**2+szglobal**2)
      END IF
   END DO
   CALL closeXMLElement('orbitalMagneticMomentsInMTSpheres')

   9020 FORMAT (/,/,10x,'orb. magnetic moments in the spheres:',/,10x,&
                'type',t22,'moment',t33,'spin-up',t43,'spin-down')
   8030 FORMAT (2x,'--> mm',i8,2x,3f12.5)
   IF (any(noco%l_unrestrictMT)) THEN
      8032 FORMAT (2x,'Atom: ',i8,' --> Orbital moment in global frame:', ' mx=',f9.5, ' my=',f9.5,' mz=',f9.5,' |m|=',f9.5)
   END IF
END SUBROUTINE orbMagMoms

END MODULE m_orbMagMoms
