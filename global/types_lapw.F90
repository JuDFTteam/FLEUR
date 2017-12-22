!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_lapw
  TYPE t_lapw
     INTEGER :: nv(2),num_local_cols(2)
     INTEGER :: nv_tot
     INTEGER :: nmat
     INTEGER :: nlotot
     INTEGER,ALLOCATABLE:: k1(:,:)
     INTEGER,ALLOCATABLE:: k2(:,:)
     INTEGER,ALLOCATABLE:: k3(:,:)
     INTEGER,ALLOCATABLE:: gvec(:,:,:) !replaces k1,k2,k3
     INTEGER,ALLOCATABLE:: kp(:,:)
     REAL,ALLOCATABLE::rk(:,:)
     REAL,ALLOCATABLE::gk(:,:,:)
     REAL,ALLOCATABLE::vk(:,:,:)
     INTEGER,ALLOCATABLE::matind(:,:)
   CONTAINS
     PROCEDURE,PASS :: spinblock =>lapw_spinblock
  END TYPE t_lapw
CONTAINS
  SUBROUTINE lapw_spinblock(lapw,l_noco,spinblock,isp,ispin,jspin,ii,jj,vpw_spin)
    IMPLICIT NONE
    CLASS(t_lapw),INTENT(IN):: lapw
    LOGICAL,INTENT(In):: l_noco
    INTEGER,INTENT(IN):: spinblock,isp
    INTEGER,INTENT(OUT):: ispin,jspin,ii,jj,vpw_spin
    IF (l_noco) THEN
       SELECT CASE (spinblock)
       CASE (1) !UPUP
          ispin=1;jspin=1;jj=0;ii=0
          vpw_spin=1
       CASE(2) !DOWNDOWN
          ispin=2;jspin=2;jj=lapw%nv(1)+lapw%nlotot;ii=lapw%nv(1)+lapw%nlotot
          vpw_spin=2
       CASE(3) !DOWNUP
          ispin=2;jspin=1;jj=0;ii=lapw%nv(1)+lapw%nlotot
          vpw_spin=3
       END SELECT
    ELSE
       jj=0;ii=0
       vpw_spin=isp !input spin in collinear case
       ispin=isp
       jspin=isp      
    ENDIF
  END SUBROUTINE lapw_spinblock
END MODULE m_types_lapw
