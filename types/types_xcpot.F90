!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_types_xcpot
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: t_xcpot,t_gradients
       
  TYPE,ABSTRACT :: t_xcpot
     REAL :: gmaxxc
   CONTAINS
     PROCEDURE        :: is_gga=>xcpot_is_gga
     PROCEDURE        :: is_hybrid=>xcpot_is_hybrid
     PROCEDURE        :: get_exchange_weight=>xcpot_get_exchange_weight
     PROCEDURE        :: get_vxc=>xcpot_get_vxc
     PROCEDURE        :: get_exc=>xcpot_get_exc
     PROCEDURE        :: broadcast=>xcpot_broadcast
     PROCEDURE,NOPASS :: alloc_gradients=>xcpot_alloc_gradients
  END TYPE t_xcpot

  TYPE t_gradients
     !Naming convention:
     !t,u,d as last letter for total,up,down
     !agr for absolute value of gradient
     !g2r for laplacien of gradient
     !+??
     REAL,ALLOCATABLE :: agrt(:),agru(:),agrd(:) 
     REAL,ALLOCATABLE :: g2ru(:),g2rd(:),gggrt(:)
     REAL,ALLOCATABLE :: gggru(:),gzgr(:),g2rt(:)
     REAL,ALLOCATABLE :: gggrd(:),grgru(:),grgrd(:)
     !These are the contracted Gradients used in libxc
     REAL,ALLOCATABLE :: sigma(:,:)
  END TYPE t_gradients

CONTAINS
  
  LOGICAL FUNCTION xcpot_is_gga(xcpot)
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    xcpot_is_gga=.false.
  END FUNCTION xcpot_is_gga

  LOGICAL FUNCTION xcpot_is_hybrid(xcpot)
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    xcpot_is_hybrid=.FALSE.
  END FUNCTION xcpot_is_hybrid
  
  FUNCTION xcpot_get_exchange_weight(xcpot) RESULT(a_ex)
    USE m_judft
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(IN):: xcpot
    REAL:: a_ex
    a_ex=-1
  END FUNCTION xcpot_get_exchange_weight

  SUBROUTINE xcpot_get_vxc(xcpot,jspins,rh,vxc,vx,grad)
    CLASS(t_xcpot),INTENT(IN) :: xcpot
    INTEGER, INTENT (IN)     :: jspins
    !--> charge density
    REAL,INTENT (IN)         :: rh(:,:)
    !---> xc potential
    REAL, INTENT (OUT)       :: vxc (:,:),vx(:,:)
    TYPE(t_gradients),OPTIONAL,INTENT(IN)::grad
  END SUBROUTINE xcpot_get_vxc

  
  SUBROUTINE xcpot_get_exc(xcpot,jspins,rh,exc,grad)
    CLASS(t_xcpot),INTENT(IN) :: xcpot
    INTEGER, INTENT (IN)     :: jspins
    !--> charge density
    REAL,INTENT (IN)         :: rh(:,:)
    !---> xc energy density
    REAL, INTENT (OUT)       :: exc (:)
    TYPE(t_gradients),OPTIONAL,INTENT(IN)::grad
  END SUBROUTINE xcpot_get_exc

  SUBROUTINE xcpot_alloc_gradients(ngrid,jspins,grad)
    INTEGER, INTENT (IN)         :: jspins,ngrid
    TYPE(t_gradients),INTENT(OUT):: grad

    !For the in-build xc-pots
    ALLOCATE(grad%agrt(ngrid),grad%agru(ngrid),grad%agrd(ngrid))
    ALLOCATE(grad%g2ru(ngrid),grad%g2rd(ngrid),grad%gggrt(ngrid))
    ALLOCATE(grad%gggru(ngrid),grad%gzgr(ngrid),grad%g2rt(ngrid))
    ALLOCATE(grad%gggrd(ngrid),grad%grgru(ngrid),grad%grgrd(ngrid))
  END SUBROUTINE xcpot_alloc_gradients

  SUBROUTINE xcpot_broadcast(xcpot,mpi)
    USE m_types_mpi
    IMPLICIT NONE
    CLASS(t_xcpot),INTENT(INOUT) :: xcpot
    TYPE(t_mpi),INTENT(IN)       :: mpi
  END SUBROUTINE xcpot_broadcast
  
END MODULE m_types_xcpot
