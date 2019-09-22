!--------------------------------------------------------------------------------
! Copyright (c) 2016 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_umix
  USE m_juDFT
  !
  ! mix the old and new density matrix for the lda+U method
  !                                                 gb.2001
  ! --------------------------------------------------------
  ! Extension to multiple U per atom type by G.M. 2017
CONTAINS
  SUBROUTINE u_mix(input,atoms,noco,n_mmp_in,n_mmp_out)

    USE m_types
    USE m_cdn_io
    USE m_nmat_rot
    USE m_xmlOutput

    ! ... Arguments

    IMPLICIT NONE
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_atoms),INTENT(IN)   :: atoms
    TYPE(t_noco),INTENT(IN)   :: noco
    COMPLEX, INTENT (INOUT)    :: n_mmp_out(-3:3,-3:3,atoms%n_u,MERGE(3,input%jspins,noco%l_mperp))
    COMPLEX, INTENT (INOUT)    :: n_mmp_in (-3:3,-3:3,atoms%n_u,MERGE(3,input%jspins,noco%l_mperp))
    !
    ! ... Locals ...
    INTEGER j,k,iofl,l,itype,ios,i_u,jsp
    REAL alpha,spinf,gam,del,sum1,sum2,sum3,mix_u, uParam, jParam
    REAL    zero(atoms%n_u)
    CHARACTER(LEN=20)   :: attributes(6)
    COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:)
    !
    ! check for possible rotation of n_mmp
    !
    !zero=0.0
    !CALL nmat_rot(zero,-atoms%lda_u%theta,-atoms%lda_u%phi,3,atoms%n_u,input%jspins,atoms%lda_u%l,n_mmp_out)

    ! Write out n_mmp_out to out.xml file

    CALL openXMLElementNoAttributes('ldaUDensityMatrix')
    DO jsp = 1, MERGE(3,input%jspins,noco%l_mperp)
       DO i_u = 1, atoms%n_u
          l = atoms%lda_u(i_u)%l
          itype = atoms%lda_u(i_u)%atomType
          uParam = atoms%lda_u(i_u)%u
          jParam = atoms%lda_u(i_u)%j
          attributes = ''
          WRITE(attributes(1),'(i0)') jsp
          WRITE(attributes(2),'(i0)') itype
          WRITE(attributes(3),'(i0)') i_u
          WRITE(attributes(4),'(i0)') l
          WRITE(attributes(5),'(f15.8)') uParam
          WRITE(attributes(6),'(f15.8)') jParam
          CALL writeXMLElementMatrixPoly('densityMatrixFor',&
                                         (/'spin    ','atomType','uIndex  ','l       ','U       ','J       '/),&
                                         attributes,n_mmp_out(-l:l,-l:l,i_u,jsp))
       END DO
    END DO
    CALL closeXMLElement('ldaUDensityMatrix')

    ! exit subroutine if density matrix does not exist
    IF(.NOT.ANY(n_mmp_in(:,:,:,:).NE.0.0)) THEN
       RETURN
    END IF

    IF (input%ldauLinMix) THEN

       ! mix here straight with given mixing factors

       ALLOCATE (n_mmp(-3:3,-3:3,MAX(1,atoms%n_u),MERGE(3,input%jspins,noco%l_mperp)))
       n_mmp = CMPLX(0.0,0.0)

       alpha = input%ldauMixParam
       spinf = input%ldauSpinf

       sum1 = 0.0
       IF (input%jspins.EQ.1) THEN
          DO i_u = 1, atoms%n_u
             DO j = -3,3
                DO k = -3,3
                   sum1 = sum1 + ABS(n_mmp_out(k,j,i_u,1) - n_mmp_in(k,j,i_u,1))
                   n_mmp(k,j,i_u,1) = alpha * n_mmp_out(k,j,i_u,1) + (1.0-alpha) * n_mmp_in(k,j,i_u,1)
                END DO
             END DO
          END DO
          WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
       ELSE
          sum2 = 0.0
          sum3 = 0.0
          gam = 0.5 * alpha * (1.0 + spinf)
          del = 0.5 * alpha * (1.0 - spinf)
          DO i_u = 1,atoms%n_u
             DO j = -3,3
                DO k = -3,3
                   sum1 = sum1 + ABS(n_mmp_out(k,j,i_u,1) - n_mmp_in(k,j,i_u,1))
                   sum2 = sum2 + ABS(n_mmp_out(k,j,i_u,2) - n_mmp_in(k,j,i_u,2))
                   IF(noco%l_mperp) sum3 = sum3 + ABS(n_mmp_out(k,j,i_u,3) - n_mmp_in(k,j,i_u,3))

                   n_mmp(k,j,i_u,1) =       gam * n_mmp_out(k,j,i_u,1) + &
                                      (1.0-gam) * n_mmp_in (k,j,i_u,1) + &
                                            del * n_mmp_out(k,j,i_u,2) - &
                                            del * n_mmp_in (k,j,i_u,2)

                   n_mmp(k,j,i_u,2) =       gam * n_mmp_out(k,j,i_u,2) + &
                                      (1.0-gam) * n_mmp_in (k,j,i_u,2) + &
                                            del * n_mmp_out(k,j,i_u,1) - &
                                            del * n_mmp_in (k,j,i_u,1)
                   IF(noco%l_mperp) THEN
                      n_mmp(k,j,i_u,3) =       alpha * n_mmp_out(k,j,i_u,3) + &
                                         (1.0-alpha) * n_mmp_in (k,j,i_u,3)
                   ENDIF
                END DO
             END DO
          END DO
          WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
          WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
          IF(noco%l_mperp) WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 3 =',sum3
       ENDIF
       n_mmp_in = n_mmp
       DEALLOCATE (n_mmp)
    ELSE ! input%ldauLinMix

       ! only calculate distance

       sum1 = 0.0
       DO i_u = 1, atoms%n_u
          DO j = -3,3
             DO k = -3,3
                sum1 = sum1 + ABS(n_mmp_out(k,j,i_u,1) - n_mmp_in(k,j,i_u,1))
             END DO
          END DO
       END DO
       IF (input%jspins.EQ.1) THEN
          WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
       ELSE
          sum2 = 0.0
          WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
          DO i_u = 1, atoms%n_u
             DO j = -3,3
                DO k = -3,3
                   sum2 = sum2 + ABS(n_mmp_out(k,j,i_u,2) - n_mmp_in(k,j,i_u,2))
                END DO
             END DO
          END DO
          !DO j=-3,3
          !   WRITE(6,'(14f12.6)') (n_mmp_in(k,j,1,2),k=-3,3)
          !END DO
          WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
          !DO j=-3,3
          !   WRITE(6,'(14f12.6)') (n_mmp_out(k,j,1,2),k=-3,3)
          !END DO
          IF(noco%l_mperp) THEN
            !Spin off-diagonal
            sum3 = 0.0
            DO i_u = 1, atoms%n_u
               DO j = -3,3
                  DO k = -3,3
                     sum3 = sum3 + ABS(n_mmp_out(k,j,i_u,3) - n_mmp_in(k,j,i_u,3))
                  END DO
               END DO
            END DO
            WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 3 =',sum3
          ENDIF
       END IF
    END IF ! input%ldauLinMix

  END SUBROUTINE u_mix
END MODULE m_umix
