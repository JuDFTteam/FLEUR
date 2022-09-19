!--------------------------------------------------------------------------------
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_magMoms
IMPLICIT NONE
CONTAINS

SUBROUTINE magMoms(input,atoms,noco,nococonv,vTot,moments,den)

   USE m_types
   USE m_constants
   USE m_xmlOutput
   USE m_m_perp
   USE m_intgr, ONLY : intgr3



   TYPE(t_input), INTENT(IN)       :: input
   TYPE(t_atoms), INTENT(IN)       :: atoms
   TYPE(t_noco), INTENT(IN)        :: noco
   TYPE(t_nococonv), INTENT(INOUT) :: nococonv
   TYPE(t_potden),INTENT(IN),OPTIONAL  :: vTot
   TYPE(t_moments),INTENT(IN),OPTIONAL :: moments
   TYPE(t_potden),INTENT(IN),OPTIONAL  :: den

   INTEGER                       :: iType, j
   REAL                          :: sval,stot,scor,smom,off1,off2,up,down
   CHARACTER(LEN=20)             :: attributes(4)
   COMPLEX                       :: off_diag


   if (present(moments)) THEN
     WRITE (oUnit,FMT=8000)
     DO iType = 1,atoms%ntype
        sval = moments%svdn(iType,1) - moments%svdn(iType,input%jspins)
        stot = moments%stdn(iType,1) - moments%stdn(iType,input%jspins)
        scor = stot - sval
        WRITE (oUnit,FMT=8010) iType,stot,sval,scor,moments%svdn(iType,1),moments%stdn(iType,1)
     END DO
   endif
   8000 FORMAT (/,/,10x,'spin density at the nucleus:',/,10x,'type',t25,&
                'total',t42,'valence',t65,'core',t90,&
                'majority valence and total density',/)
   8010 FORMAT (i13,2x,3e20.8,5x,2e20.8)

   write(ounit,*)
   write(oUnit,'(a,19x,a,19x,a,19x,a)') "Magnetic moments |","Global Frame","  | ","Local Frame"
   write(oUnit,*) "------------------------------------------------------------------------------------------------------------------------"
   write(oUnit,'(a,5x,a,2(" | ",5(a,5x)))') "Atom ","|m|   ","mx   ","my   ","mz   ","alpha","beta ","mx   ","my   ","mz   ","alpha","beta "
   
   !WRITE (oUnit,FMT=8020)

   CALL openXMLElement('magneticMomentsInMTSpheres',(/'units'/),(/'muBohr'/))
   DO iType = 1, atoms%ntype
      if (present(moments)) THEN
         up=moments%chmom(iType,1)
         down=moments%chmom(iType,input%jspins)
      elseif(present(den)) THEN
         CALL intgr3(den%mt(:,0,iType,1),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),up)
         CALL intgr3(den%mt(:,0,iType,input%jspins),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),down)
         up=up*sfp_const;down=down*sfp_const
      else
         return !no data found
      endif   
      smom = up-down
      !WRITE (oUnit,FMT=8030) iType,smom, up,down
      attributes = ''
      WRITE(attributes(1),'(i0)') iType
      WRITE(attributes(2),'(f15.10)') smom
      WRITE(attributes(3),'(f15.10)') up
      WRITE(attributes(4),'(f15.10)') down
      CALL writeXMLElementFormPoly('magneticMoment',(/'atomType      ','moment        ','spinUpCharge  ',&
                                                      'spinDownCharge'/),&
                                   attributes,reshape((/8,6,12,14,6,15,15,15/),(/4,2/)))

      off_diag=0.0
      if (present(moments).and.noco%l_mperp) THEN
          off_diag=moments%qa21(itype)
      elseif (present(den)) THEN
         if (size(den%mt,4)==4) THEN
            CALL intgr3(den%mt(:,0,iType,3),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),off1)
            CALL intgr3(den%mt(:,0,iType,4),atoms%rmsh(:,iType),atoms%dx(iType),atoms%jri(iType),off2)
            off_diag=cmplx(off1,off2)*sfp_const
         endif
      endif      
      !calculate the perpendicular part of the local moment
      !and relax the angle of the local moment or calculate
      !the constraint B-field.
      CALL m_perp(atoms,iType,noco,nococonv,(/up,down/),off_diag,vTot)
         
   END DO
   CALL closeXMLElement('magneticMomentsInMTSpheres')

   8020 FORMAT (/,/,2x,'-->  magnetic moments in the spheres:',/,2x,&
                'mm -->   type',t22,'moment',t33,'spin-up',t43,'spin-down')
   8030 FORMAT (2x,'--> mm',i8,2x,3f12.5)
   write(oUnit,*) "------------------------------------------------------------------------------------------------------------------------"
   write(oUnit,*)
END SUBROUTINE magMoms

END MODULE m_magMoms
