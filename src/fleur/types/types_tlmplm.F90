!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_tlmplm
   use m_types_rsoc
  IMPLICIT NONE
  PRIVATE
 

  TYPE t_tlmplm
     COMPLEX,ALLOCATABLE :: tuloulo_newer(:,:,:,:,:,:,:) !mp,m,lop,lo,ntype,ispin,jspin
     COMPLEX,ALLOCATABLE :: h_loc_LO(:,:,:,:,:)    !lm,lmp,ntype,ispin,jspin
     COMPLEX,ALLOCATABLE :: h_LO(:,:,:,:,:)    !lmp,m,lo+mlo,ispin,jspin
     COMPLEX,ALLOCATABLE :: h_LO2(:,:,:,:,:)    !lmp,m,lo+mlo,ispin,jspin
     COMPLEX,ALLOCATABLE :: h_loc(:,:,:,:,:)    !lm,lmp,ntype,ispin,jspin
     COMPLEX,ALLOCATABLE :: h_loc_nonsph(:,:,:,:,:)    !lm,lmp,ntype,ispin,jspin
     INTEGER,ALLOCATABLE :: h_loc2(:)
     INTEGER,ALLOCATABLE :: h_loc2_nonsph(:)

     REAL,ALLOCATABLE    :: e_shift(:,:)
     TYPE(t_rsoc)        :: rsoc
   CONTAINS
     PROCEDURE,PASS :: init => tlmplm_init
  END TYPE t_tlmplm
  PUBLIC t_tlmplm
CONTAINS
  SUBROUTINE tlmplm_init(td,atoms,jspins)
    USE m_judft
    USE m_types_atoms
    CLASS(t_tlmplm),INTENT(INOUT):: td
    TYPE(t_atoms)                :: atoms
    INTEGER,INTENT(in)           :: jspins
    INTEGER :: err(7),lmd
    err = 0
    lmd=atoms%lmaxd*(atoms%lmaxd+2)

    td%h_loc2=atoms%lmax*(atoms%lmax+2)+1
    td%h_loc2_nonsph=atoms%lnonsph*(atoms%lnonsph+2)+1
    IF (ALLOCATED(td%h_loc)) &
         DEALLOCATE(td%tuloulo_newer,td%h_loc,td%e_shift,td%h_loc_nonsph,td%h_loc_LO,td%h_lo,td%h_lo2)
    ALLOCATE(td%tuloulo_newer(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,atoms%nlod,atoms%nlod,atoms%ntype,jspins,jspins), stat=err(1));td%tuloulo_newer=0.0
    ALLOCATE(td%h_loc(0:2*lmd+1,0:2*lmd+1,atoms%ntype,jspins,jspins),stat=err(2));td%h_loc=0.0
    ALLOCATE(td%h_loc_nonsph(0:MAXVAL(td%h_loc2_nonsph)*2-1,0:MAXVAL(td%h_loc2_nonsph)*2-1,atoms%ntype,jspins,jspins),stat=err(3));td%h_loc_nonsph=0.0
    ALLOCATE(td%h_loc_lo(0:MAXVAL(td%h_loc2_nonsph)*2-1,0:MAXVAL(td%h_loc2_nonsph)*2-1,atoms%ntype,jspins,jspins),stat=err(4));td%h_loc_lo=0.0
    ALLOCATE(td%h_lo(0:MAXVAL(td%h_loc2_nonsph)*2-1,-atoms%llod:atoms%llod,SUM(atoms%nlo),jspins,jspins),stat=err(5));td%h_lo=0.0
    ALLOCATE(td%h_lo2(0:MAXVAL(td%h_loc2_nonsph)*2-1,-atoms%llod:atoms%llod,SUM(atoms%nlo),jspins,jspins),stat=err(6));td%h_lo2=0.0
    ALLOCATE(td%e_shift(atoms%ntype,jspins),stat=err(7))
    IF (ANY(err.NE.0)) THEN
       WRITE (*,*) 'an error occured during allocation of'
       WRITE (*,*) 'the tlmplm local matrix elements'
       WRITE (*,'(9i7)') err(:)
       CALL juDFT_error("eigen: Error during allocation of tlmplm",calledby ="types_tlmplm")
    ENDIF
  END SUBROUTINE tlmplm_init

END MODULE m_types_tlmplm
