MODULE m_gen_map

IMPLICIT NONE

CONTAINS

!
! build up field map(iatom,isym), which contains the number of the atom, on 
! which the atom iatom is mapped via the symmetry operation isym
! tvec is the translation, which maps R R_a + tau back in the unit cell
!
SUBROUTINE gen_map(atoms,sym,oneD,hybrid)
  USE m_types
  TYPE(t_atoms),INTENT(IN) :: atoms
  TYPE(t_sym),INTENT(IN)   :: sym
  TYPE(t_oneD),INTENT(IN)  :: oneD
  TYPE(t_hybrid),INTENT(INOUT)::hybrid
  ! private scalars
  INTEGER                           :: iatom,iatom0,itype,ieq,isym,iisym,ieq1
  INTEGER                           :: ratom,ok
  ! private arrays
  REAL                              :: rtaual(3)

  ALLOCATE( hybrid%map(atoms%nat,sym%nsym) , stat = ok )
  IF( ok .ne. 0 ) STOP 'gen_map: error during allocation of map'

  ALLOCATE( hybrid%tvec(3,atoms%nat,sym%nsym) , stat = ok )
  IF( ok .ne. 0 ) STOP 'gen_map: error during allocation of tvec'

  iatom  = 0
  iatom0 = 0
  DO itype = 1,atoms%ntype
    DO ieq = 1,atoms%neq(itype)
      iatom = iatom + 1
      DO isym = 1,sym%nsym

        IF( isym .le. sym%nop ) THEN
          iisym = isym
        ELSE
          iisym = isym - sym%nop
        END IF

        rtaual(:) = matmul(sym%mrot(:,:,iisym),atoms%taual(:,iatom)) + sym%tau(:,iisym)

        ratom = 0
        DO ieq1 = 1,atoms%neq(itype)
          IF( all(abs(modulo(rtaual-atoms%taual(:,iatom0 + ieq1)+10.0**-12,1.0)).lt. 10.0**-10) ) THEN
            ratom              = iatom0 + ieq1
            hybrid%map (  iatom,isym) = ratom
            hybrid%tvec(:,iatom,isym) = nint ( rtaual-atoms%taual(:,ratom) )
            CYCLE
          END IF
        END DO
        IF( ratom .eq. 0 ) STOP 'eigen_hf: ratom not found'

      END DO
    END DO
    iatom0 = iatom0 + atoms%neq(itype)
  END DO

END SUBROUTINE gen_map

END MODULE m_gen_map
