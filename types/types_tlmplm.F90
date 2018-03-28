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
     COMPLEX,ALLOCATABLE,DIMENSION(:,:,:,:,:,:,:)::soangl
  END TYPE t_rsoc

    TYPE t_tlmplm
     COMPLEX,ALLOCATABLE :: tdd(:,:,:)
     COMPLEX,ALLOCATABLE :: tdu(:,:,:)
     !(0:lmplmd,ntypd,tspin)
     COMPLEX,ALLOCATABLE :: tud(:,:,:)
     COMPLEX,ALLOCATABLE :: tuu(:,:,:)
     !(0:lmplmd,ntypd,tspin)
     INTEGER,ALLOCATABLE :: ind(:,:,:,:)
     !(0:lmd,0:lmd,ntypd,tspin)
     COMPLEX,ALLOCATABLE :: tdulo(:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuulo(:,:,:,:)
     !(0:lmd,-llod:llod,mlotot,tspin)
     COMPLEX,ALLOCATABLE :: tuloulo(:,:,:,:)
     !(-llod:llod,-llod:llod,mlolotot,tspin)
     COMPLEX,ALLOCATABLE :: h_loc(:,:,:,:)
     COMPLEX,ALLOCATABLE :: h_off(:,:,:,:)
     REAL,ALLOCATABLE    :: e_shift(:,:)
     TYPE(t_rsoc)        :: rsoc
   CONTAINS
     PROCEDURE,PASS :: init => tlmplm_init
  END TYPE t_tlmplm
  PUBLIC t_tlmplm,t_rsoc
CONTAINS
  SUBROUTINE tlmplm_init(td,lmplmd,lmd,ntype,lmaxd,llod,mlotot,mlolotot,jspins,l_offdiag)
    USE m_judft
    CLASS(t_tlmplm),INTENT(INOUT):: td
    INTEGER,INTENT(in)           :: lmplmd,lmd,ntype,lmaxd,llod,mlotot,mlolotot,jspins
    LOGICAL,INTENT(IN)           :: l_offdiag
    INTEGER :: err

    IF (ALLOCATED(td%tuu)) &
         DEALLOCATE(td%tuu,td%tud,td%tdd,td%tdu,td%tdulo,td%tuulo,&
         td%tuloulo,td%ind,td%h_loc,td%e_shift,td%h_off)
    ALLOCATE(td%tuu(0:lmplmd,ntype,jspins),stat=err)
    ALLOCATE(td%tud(0:lmplmd,ntype,jspins),stat=err)
    ALLOCATE(td%tdd(0:lmplmd,ntype,jspins),stat=err)
    ALLOCATE(td%tdu(0:lmplmd,ntype,jspins),stat=err)
    ALLOCATE(td%tdulo(0:lmd,-llod:llod,mlotot,jspins),stat=err)
    ALLOCATE(td%tuulo(0:lmd,-llod:llod,mlotot,jspins),stat=err)
    ALLOCATE(td%tuloulo(-llod:llod,-llod:llod,MAX(mlolotot,1),jspins), stat=err)
    ALLOCATE(td%ind(0:lmd,0:lmd,ntype,jspins),stat=err )
    ALLOCATE(td%h_loc(0:2*lmaxd*(lmaxd+2)+1,0:2*lmaxd*(lmaxd+2)+1,ntype,jspins))
    ALLOCATE(td%e_shift(ntype,jspins))
    IF (l_offdiag) THEN
       ALLOCATE(td%h_off(0:2*lmaxd+1,0:2*lmaxd+1,ntype,2))
    ELSE
       ALLOCATE(td%h_off(1,1,1,1))
    END IF
    IF (err.NE.0) THEN
       WRITE (*,*) 'an error occured during allocation of'
       WRITE (*,*) 'the tlmplm%tuu, tlmplm%tdd etc.: ',err,'  size: ',mlotot
       CALL juDFT_error("eigen: Error during allocation of tlmplm, tdd  etc.",calledby ="types_tlmplm")
    ENDIF
  END SUBROUTINE tlmplm_init
    
END MODULE m_types_tlmplm
