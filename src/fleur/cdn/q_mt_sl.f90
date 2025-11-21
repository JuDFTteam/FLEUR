!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_qmtsl
   implicit none
CONTAINS
  !***********************************************************************
  ! Calculates the mt-spheres contribution to the layer charge for states
  !  {En} at the current k-point.
  !                                      Yury Koroteev 2003
  !                     from eparas.F  by  Philipp Kurz 99/04
  !
  !***********************************************************************
  !
  SUBROUTINE q_mt_sl(itype,jsp,ikpt,atoms,ev_list,ne,abc,radfun,slab)
    USE m_types_setup
    USE m_types_abc
    USE m_types_radfun
    USE m_types_slab
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)        :: atoms
    TYPE(t_abc),INTENT(IN)          :: abc
    TYPE(t_slab), INTENT(INOUT)     :: slab
    TYPE(t_radfun), INTENT(IN)      :: radfun
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN) :: jsp,itype
    INTEGER, INTENT (IN) :: ne,ikpt 

    INTEGER, INTENT (IN) :: ev_list(:)

    !     ..
    !     .. Local Scalars ..
    INTEGER:: i,l,natom,m,j,jj
    INTEGER:: lm,ll1,nl
    REAL :: sabd
    


    
    DO i = 1,ne  
      sabd = 0.0  !sum over l,m,natoms,j,jj
      DO l = 0,atoms%lmax(itype)
         ll1 = l* (l+1)
         DO m = -l,l
            lm = ll1 + m
            DO natom = 1, atoms%neq(itype)
            DO j = 1, radfun%n_r(l)
               DO jj = 1, radfun%n_r(l)
                  sabd = sabd + abc%cof(i, lm, j, natom)*CONJG(abc%cof(i, lm, jj, natom))*radfun%integral(j, jj, l, jsp, jsp)
               END DO
            END DO
            ENDDO
         enddo
      enddo
      !Map to slabs
       DO nl = 1,slab%nsl
          slab%qmtsl(nl,ev_list(i),ikpt,jsp) = sabd/atoms%neq(itype)*slab%nmtsl(itype,nl)
       ENDDO
    ENDDO
    
  END SUBROUTINE q_mt_sl
END MODULE m_qmtsl
