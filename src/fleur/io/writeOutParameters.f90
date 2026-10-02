!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_writeOutParameters

   USE m_xmlOutput
   USE m_types_lapw
   USE m_types_atoms
   USE m_types_banddos
   USE m_types_cell
   USE m_types_enpara
   USE m_types_hybinp
   USE m_types_input
   USE m_types_kpts
   USE m_types_mpi
   USE m_types_noco
   USE m_types_sliceplot
   USE m_types_sphhar
   USE m_types_stars
   USE m_types_sym
   USE m_types_vacuum
   USE m_types_xcpot
IMPLICIT NONE
   PRIVATE
   PUBLIC :: writeoutparameters

CONTAINS

SUBROUTINE writeOutParameters(fmpi,input,sym,stars,atoms,vacuum,kpts,&
                               hybinp,cell,banddos,sliceplot,xcpot,&
                              noco,enpara,sphhar)


   TYPE(t_mpi),       INTENT(IN) :: fmpi
   TYPE(t_input),     INTENT(IN) :: input
   TYPE(t_sym),       INTENT(IN) :: sym
   TYPE(t_stars),     INTENT(IN) :: stars
   TYPE(t_atoms),     INTENT(IN) :: atoms
   TYPE(t_vacuum),    INTENT(IN) :: vacuum
   TYPE(t_kpts),      INTENT(IN) :: kpts
    
   TYPE(t_hybinp),    INTENT(IN) :: hybinp
   TYPE(t_cell),      INTENT(IN) :: cell
   TYPE(t_banddos),   INTENT(IN) :: banddos
   TYPE(t_sliceplot), INTENT(IN) :: sliceplot
   CLASS(t_xcpot),    INTENT(IN) :: xcpot
   TYPE(t_noco),      INTENT(IN) :: noco

   TYPE(t_enpara),    INTENT(IN) :: enpara
   TYPE(t_sphhar),    INTENT(IN) :: sphhar

   INTEGER            :: i
   REAL               :: sumWeight
   CHARACTER(LEN=20)  :: attributes(7)

   CALL openXMLElementNoAttributes('numericalParameters')

   WRITE(attributes(1),'(i0)') atoms%nat
   WRITE(attributes(2),'(i0)') atoms%ntype
   WRITE(attributes(3),'(i0)') atoms%jmtd
   WRITE(attributes(4),'(i0)') atoms%n_u
   WRITE(attributes(5),'(i0)') atoms%n_hia
   WRITE(attributes(6),'(i0)') atoms%n_opc
   CALL writeXMLElementFormPoly('atomsInCell',(/'nat  ','ntype','jmtd ','n_u  ','n_hia', 'n_opc'/),&
                                attributes(:6),reshape((/3,6,6,6,6,6,8,8,8,8,8/),(/5,2/)))

   WRITE(attributes(1),'(i0)') lapw_dim_nvd
   WRITE(attributes(2),'(i0)') atoms%lmaxd
   WRITE(attributes(3),'(i0)') atoms%nlotot
   CALL writeXMLElementFormPoly('basis',(/'nvd   ','lmaxd ','nlotot'/),&
                                attributes(:3),reshape((/9,6,6,8,8,8/),(/3,2/)))

   WRITE(attributes(1),'(i0)') stars%ng3
   WRITE(attributes(2),'(i0)') stars%ng2
   CALL writeXMLElementFormPoly('density',(/'ng3','ng2'/),&
                                attributes(:2),reshape((/7,6,8,8/),(/2,2/)))

   WRITE(attributes(1),'(i0)') input%neig
   CALL writeXMLElementFormPoly('bands',(/'numbands'/),&
                                attributes(:1),reshape((/9,8/),(/1,2/)))

   WRITE(attributes(1),'(f0.8)') cell%vol
   WRITE(attributes(2),'(f0.8)') cell%volint
   IF(input%film) THEN
      WRITE(attributes(3),'(f0.8)') cell%omtil
      WRITE(attributes(4),'(f0.8)') cell%area
      WRITE(attributes(5),'(f0.8)') cell%z1
      CALL openXMLElementFormPoly('volumes',(/'unitCell    ', 'interstitial', 'omegaTilda  ', 'surfaceArea ', 'z1          '/),&
                                  attributes(:5),reshape((/8,12,10,11,2,10,10,10,10,10/),(/5,2/)))
   ELSE
      CALL openXMLElementFormPoly('volumes',(/'unitCell    ', 'interstitial'/),&
                                  attributes(:2),reshape((/8,12,10,10/),(/2,2/)))
   END IF
   DO i = 1, atoms%ntype
      WRITE(attributes(1),'(i0)') i
      WRITE(attributes(2),'(f0.8)') atoms%rmt(i)
      WRITE(attributes(3),'(f0.8)') atoms%volmts(i)
      CALL writeXMLElementFormPoly('mtVolume',(/'atomType','mtRadius','mtVolume'/),&
                                   attributes(:3),reshape((/8,8,8,5,10,10/),(/3,2/)))
   END DO
   CALL closeXMLElement('volumes')

   sumWeight = SUM(kpts%wtkpt(:kpts%nkpt))
   WRITE(attributes(1),'(f0.8)') sumWeight
   WRITE(attributes(2),'(i0)') kpts%nkpt
   CALL openXMLElementFormPoly('kPointList',(/ 'weightScale', 'count      '/),&
                               attributes(:2),reshape((/11,5,10,10/),(/2,2/)))
   DO i = 1, kpts%nkpt
      WRITE(attributes(1),'(f12.6)') kpts%wtkpt(i)
      WRITE(attributes(2),'(f12.6)') kpts%bk(1,i)
      WRITE(attributes(3),'(f12.6)') kpts%bk(2,i)
      WRITE(attributes(4),'(f12.6)') kpts%bk(3,i)
      CALL writeXMLElementForm('kPoint', (/'weight'/), attributes(:1),reshape((/6,12/),(/1,2/)),attributes(2:4))
   END DO
   CALL closeXMLElement('kPointList')

   CALL closeXMLElement('numericalParameters')

END SUBROUTINE writeOutParameters

END MODULE m_writeOutParameters
