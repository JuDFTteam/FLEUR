MODULE m_writeOutParameters

IMPLICIT NONE

CONTAINS

SUBROUTINE writeOutParameters(mpi,input,sym,stars,atoms,vacuum,obsolete,kpts,&
                              oneD,hybrid,cell,banddos,sliceplot,xcpot,&
                              noco,dimension,enpara,sphhar)

   USE m_types
   USE m_xmlOutput

   TYPE(t_mpi),       INTENT(IN) :: mpi
   TYPE(t_input),     INTENT(IN) :: input
   TYPE(t_sym),       INTENT(IN) :: sym
   TYPE(t_stars),     INTENT(IN) :: stars 
   TYPE(t_atoms),     INTENT(IN) :: atoms
   TYPE(t_vacuum),    INTENT(IN) :: vacuum
   TYPE(t_obsolete),  INTENT(IN) :: obsolete
   TYPE(t_kpts),      INTENT(IN) :: kpts
   TYPE(t_oneD),      INTENT(IN) :: oneD
   TYPE(t_hybrid),    INTENT(IN) :: hybrid
   TYPE(t_cell),      INTENT(IN) :: cell
   TYPE(t_banddos),   INTENT(IN) :: banddos
   TYPE(t_sliceplot), INTENT(IN) :: sliceplot
   CLASS(t_xcpot),    INTENT(IN) :: xcpot
   TYPE(t_noco),      INTENT(IN) :: noco
   TYPE(t_dimension), INTENT(IN) :: dimension
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
   CALL writeXMLElementFormPoly('atomsInCell',(/'nat  ','ntype','jmtd ','n_u  '/),&
                                attributes(:4),reshape((/3,6,6,6,8,8,8,8/),(/4,2/)))

   WRITE(attributes(1),'(i0)') dimension%nvd
   WRITE(attributes(2),'(i0)') atoms%lmaxd
   WRITE(attributes(3),'(i0)') atoms%nlotot
   CALL writeXMLElementFormPoly('basis',(/'nvd   ','lmaxd ','nlotot'/),&
                                attributes(:3),reshape((/9,6,6,8,8,8/),(/3,2/)))

   WRITE(attributes(1),'(i0)') stars%ng3
   WRITE(attributes(2),'(i0)') stars%ng2
   CALL writeXMLElementFormPoly('density',(/'ng3','ng2'/),&
                                attributes(:2),reshape((/7,6,8,8/),(/2,2/)))

   WRITE(attributes(1),'(i0)') dimension%neigd
   CALL writeXMLElementFormPoly('bands',(/'numbands'/),&
                                attributes(:1),reshape((/9,8/),(/1,2/)))

   sumWeight = SUM(kpts%wtkpt(:kpts%nkpt))
   WRITE(attributes(1),'(f0.8)') kpts%posScale
   WRITE(attributes(2),'(f0.8)') sumWeight
   WRITE(attributes(3),'(i0)') kpts%nkpt
   CALL openXMLElementFormPoly('kPointList',(/'posScale   ', 'weightScale', 'count      '/),&
                               attributes(:3),reshape((/8,11,5,10,10,5/),(/3,2/)))
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
