MODULE m_forcea21U
CONTAINS
  SUBROUTINE force_a21_U(nobd,atoms,lmaxb,i_u,itype,isp,we,ne,&
       usdus,v_mmp, acof,bcof,ccof,aveccof,bveccof,cveccof, a21)
    !
    !***********************************************************************
    ! This subroutine calculates the lda+U contribution to the HF forces, 
    ! similar to the A21 term, according to eqn. (22) of F. Tran et al.
    ! Comp.Phys.Comm. 179 (2008) 784-790
    !***********************************************************************
    !
    USE m_types
    IMPLICIT NONE

    TYPE(t_usdus),INTENT(IN)   :: usdus
    TYPE(t_atoms),INTENT(IN)   :: atoms
    !     ..
    !     .. Scalar Arguments ..
    INTEGER, INTENT (IN)    :: nobd   
    INTEGER, INTENT (IN)    :: itype,isp,ne,lmaxb
    INTEGER, INTENT (INOUT) :: i_u ! on input: index for the first U for atom type "itype or higher"
                                   ! on exit: index for the first U for atom type "itype+1 or higher"
    !     ..
    !     .. Array Arguments ..
    REAL,    INTENT (IN) :: we(nobd) 
    COMPLEX, INTENT (IN) :: v_mmp(-lmaxb:lmaxb,-lmaxb:lmaxb,atoms%n_u)
    COMPLEX, INTENT (IN) :: acof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (IN) :: bcof(:,0:,:)!(nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (IN) :: ccof(-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%nat)
    COMPLEX, INTENT (IN) :: aveccof(:,:,0:,:)!(3,nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (IN) :: bveccof(:,:,0:,:)!(3,nobd,0:dimension%lmd,atoms%nat)
    COMPLEX, INTENT (IN) :: cveccof(3,-atoms%llod:atoms%llod,nobd,atoms%nlod,atoms%nat)
    REAL, INTENT (INOUT) :: a21(3,atoms%nat)
    !     ..
    !     .. Local Scalars ..
    COMPLEX v_a,v_b,v_c,p1,p2,p3
    INTEGER lo,lop,l,lp ,mp,lm,lmp,iatom,ie,i,m
    !     ..
    !     ..
    !*************** ABBREVIATIONS *****************************************
    ! ccof       : coefficient of the local orbital function (u_lo*Y_lm)
    ! cveccof    : is defined equivalently to aveccof, but with the LO-fct.
    ! for information on nlo,llo,uulon,dulon, and uloulopn see
    ! comments in setlomap.
    !***********************************************************************

    IF (atoms%lda_u(i_u)%atomType.GT.itype) RETURN

    DO WHILE (atoms%lda_u(i_u)%atomType.EQ.itype)

       l = atoms%lda_u(i_u)%l

       !
       ! Add contribution for the regular LAPWs (like force_a21, but with
       ! the potential matrix, v_mmp, instead of the tuu, tdd ...)
       !
       DO m = -l,l
          lm = l* (l+1) + m
          DO mp = -l,l
             lmp = l* (l+1) + mp
             v_a = v_mmp(m,mp,i_u) 
             v_b = v_mmp(m,mp,i_u) * usdus%ddn(l,itype,isp) 
             DO iatom = sum(atoms%neq(:itype-1))+1,sum(atoms%neq(:itype))
                DO ie = 1,ne
                   DO i = 1,3
                      p1 = (CONJG(acof(ie,lm,iatom)) * v_a) * aveccof(i,ie,lmp,iatom)
                      p2 = (CONJG(bcof(ie,lm,iatom)) * v_b) * bveccof(i,ie,lmp,iatom) 
                      a21(i,iatom) = a21(i,iatom) + 2.0*AIMAG(p1 + p2) * we(ie)/atoms%neq(itype)
                   END DO
                END DO
             END DO
          END DO ! mp
       END DO   ! m

       !
       ! If there are also LOs on this atom, with the same l as
       ! the one of LDA+U, add another few terms
       !
       DO lo = 1,atoms%nlo(itype)
          IF (l == atoms%llo(lo,itype)) THEN
             DO m = -l,l
                lm = l* (l+1) + m
                DO mp = -l,l
                   lmp = l* (l+1) + mp
                   v_a = v_mmp(m,mp,i_u)
                   v_b = v_mmp(m,mp,i_u) * usdus%uulon(lo,itype,isp)
                   v_c = v_mmp(m,mp,i_u) * usdus%dulon(lo,itype,isp)
                   DO iatom =  sum(atoms%neq(:itype-1))+1,sum(atoms%neq(:itype))
                      DO ie = 1,ne
                         DO i = 1,3
                            p1 = v_a * (CONJG(ccof(m,ie,lo,iatom)) * cveccof(i,mp,ie,lo,iatom))
                            p2 = v_b * (CONJG(acof(ie,lm,iatom)) * cveccof(i,mp,ie,lo,iatom) + &
                                        CONJG(ccof(m,ie,lo,iatom)) * aveccof(i,ie,lmp,iatom))
                            p3 = v_c * (CONJG(bcof(ie,lm,iatom)) * cveccof(i,mp,ie,lo,iatom) + &
                                        CONJG(ccof(m,ie,lo,iatom)) * bveccof(i,ie,lmp,iatom))
                            a21(i,iatom) = a21(i,iatom) + 2.0*AIMAG(p1 + p2 + p3)*we(ie)/atoms%neq(itype)
                         END DO
                      END DO
                   END DO
                END DO
             END DO
          END IF   ! l == atoms%llo(lo,itype)
       END DO     ! lo = 1,atoms%nlo

       i_u = i_u + 1
    END DO

  END SUBROUTINE force_a21_U
END MODULE m_forcea21U
