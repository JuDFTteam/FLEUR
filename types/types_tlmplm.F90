!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_tlmplm
  IMPLICIT NONE
  PRIVATE
  TYPE t_rsoc
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: rsopp,rsoppd,rsopdp,rsopdpd     !(atoms%ntype,atoms%lmaxd,2,2)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: rsoplop,rsoplopd,rsopdplo,rsopplo!(atoms%ntype,atoms%nlod,2,2)
     REAL,ALLOCATABLE,DIMENSION(:,:,:,:,:) :: rsoploplop !(atoms%ntype,atoms%nlod,nlod,2,2)
     COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:,:,:)::soangl
  END TYPE t_rsoc

  TYPE t_tlmplm
     COMPLEX,ALLOCATABLE :: tdulo(:,:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuulo(:,:,:,:,:)
     COMPLEX,ALLOCATABLE :: ulotu(:,:,:,:,:)
     COMPLEX,ALLOCATABLE :: ulotd(:,:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuloulo(:,:,:,:,:)
     !(-llod:llod,-llod:llod,mlolotot,tspin)
     COMPLEX,ALLOCATABLE :: h_loc(:,:,:,:,:)    !lm,lmp,ntype,ispin,jspin
     INTEGER,ALLOCATABLE :: h_loc2(:)
     INTEGER,ALLOCATABLE :: h_loc2_nonsph(:)
     COMPLEX,ALLOCATABLE :: h_loc_nonsph(:,:,:,:,:)    !lm,lmp,ntype,ispin,jspin

     COMPLEX,ALLOCATABLE :: h_off(:,:,:,:,:)      !l,lp,ntype,ispin,jspin)
     REAL,ALLOCATABLE    :: e_shift(:,:)
     !COMPLEX,ALLOCATABLE :: h_loc_sp(:,:,:,:)   !l,lp,ntype,ispin,jspin
     !COMPLEX,ALLOCATABLE :: h_locLO(:,:,:,:,:)  !lm+mlo,mlo,ntype,ispin,jspin
     TYPE(t_rsoc)        :: rsoc
     ! For juPhon:
     INTEGER,ALLOCATABLE :: ind(:,:,:,:)
   CONTAINS
     PROCEDURE,PASS :: init => tlmplm_init
  END TYPE t_tlmplm
  PUBLIC t_tlmplm,t_rsoc
CONTAINS
  SUBROUTINE tlmplm_init(td,atoms,jspins,l_offdiag)
    USE m_judft
    USE m_types_atoms
    CLASS(t_tlmplm),INTENT(INOUT):: td
    TYPE(t_atoms)                :: atoms
    INTEGER,INTENT(in)           :: jspins
    LOGICAL,INTENT(IN)           :: l_offdiag
    INTEGER :: err(9),lmd,mlolotot
    err = 0
    mlolotot=DOT_PRODUCT(atoms%nlo,atoms%nlo+1)/2
    lmd=atoms%lmaxd*(atoms%lmaxd+2)
    !lmplmd=(lmd*(lmd+3))/2

    td%h_loc2=atoms%lmax*(atoms%lmax+2)+1
    td%h_loc2_nonsph=atoms%lnonsph*(atoms%lnonsph+2)+1
    IF (ALLOCATED(td%h_loc)) &
         DEALLOCATE(td%tdulo,td%tuulo,td%ulotd,td%ulotu,&
         td%tuloulo,td%h_loc,td%e_shift,td%h_off,td%h_loc_nonsph)
    !    ALLOCATE(td%tuu(0:lmplmd,ntype,jspins),stat=err)
    !    ALLOCATE(td%tud(0:lmplmd,ntype,jspins),stat=err)
    !    ALLOCATE(td%tdd(0:lmplmd,ntype,jspins),stat=err)
    !    ALLOCATE(td%tdu(0:lmplmd,ntype,jspins),stat=err)
    ALLOCATE(td%tdulo(0:lmd,-atoms%llod:atoms%llod,SUM(atoms%nlo),jspins,jspins),stat=err(1));td%tdulo=0.0
    ALLOCATE(td%tuulo(0:lmd,-atoms%llod:atoms%llod,SUM(atoms%nlo),jspins,jspins),stat=err(2));td%tuulo=0.0
    ALLOCATE(td%ulotd(0:lmd,-atoms%llod:atoms%llod,SUM(atoms%nlo),jspins,jspins),stat=err(8));td%ulotd=0.0
    ALLOCATE(td%ulotu(0:lmd,-atoms%llod:atoms%llod,SUM(atoms%nlo),jspins,jspins),stat=err(9));td%ulotu=0.0
    ALLOCATE(td%tuloulo(-atoms%llod:atoms%llod,-atoms%llod:atoms%llod,MAX(mlolotot,1),jspins,jspins), stat=err(3));td%tuloulo=0.0
    ALLOCATE(td%h_loc(0:2*lmd+1,0:2*lmd+1,atoms%ntype,jspins,jspins),stat=err(5));td%h_loc=0.0
    ALLOCATE(td%h_loc_nonsph(0:MAXVAL(td%h_loc2_nonsph)*2-1,0:MAXVAL(td%h_loc2_nonsph)*2-1,atoms%ntype,jspins,jspins),stat=err(6));td%h_loc_nonsph=0.0

    ALLOCATE(td%e_shift(atoms%ntype,jspins),stat=err(7))
    IF (l_offdiag) THEN
       ALLOCATE(td%h_off(0:2*atoms%lmaxd+1,0:2*atoms%lmaxd+1,atoms%ntype,2,2),stat=err(4))
    ELSE
       ALLOCATE(td%h_off(1,1,1,1,1),stat=err(4))
    END IF
    td%h_off=0.0
    IF (ANY(err.NE.0)) THEN
       WRITE (*,*) 'an error occured during allocation of'
       WRITE (*,*) 'the tlmplm local matrix elements'
       WRITE (*,'(9i7)') err(:)
       CALL juDFT_error("eigen: Error during allocation of tlmplm",calledby ="types_tlmplm")
    ENDIF
  END SUBROUTINE tlmplm_init

END MODULE m_types_tlmplm
