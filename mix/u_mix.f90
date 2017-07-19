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
  SUBROUTINE u_mix(atoms,jspins,n_mmp_new)

    USE m_types
    USE m_nmat_rot
    USE m_xmlOutput

    ! ... Arguments

    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    INTEGER, INTENT (IN)       :: jspins 
    COMPLEX, INTENT (INOUT)    :: n_mmp_new(-3:3,-3:3,atoms%n_u,jspins)
    !
    ! ... Locals ...
    INTEGER j,k,iofl,l,itype,ios,i_u,jsp,lty(atoms%n_u)
    REAL alpha,spinf,gam,del,sum1,sum2,mix_u, uParam, jParam
    REAL    theta(atoms%n_u),phi(atoms%n_u),zero(atoms%n_u)
    LOGICAL n_exist
    CHARACTER(LEN=20)   :: attributes(6)
    COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:),n_mmp_old(:,:,:,:)
    !
    ! check for possible rotation of n_mmp
    !
    INQUIRE (file='n_mmp_rot',exist=n_exist)
    IF (n_exist) THEN
       OPEN (68,file='n_mmp_rot',status='old',form='formatted')
       DO i_u = 1, atoms%n_u
          l = atoms%lda_u(i_u)%l
          READ(68,*,iostat=ios) theta(i_u),phi(i_u)
          IF (ios == 0) THEN
             lty(i_u) = l
          ELSE
             IF (i_u == 1)  CALL juDFT_error("ERROR reading n_mmp_rot", calledby ="u_mix")
             theta(i_u) = theta(i_u-1) ; phi(i_u) = phi(i_u-1)
             lty(i_u) = lty(i_u-1)
          END IF
       END DO
       CLOSE (68)
       zero = 0.0
       CALL nmat_rot(zero,-theta,-phi,3,atoms%n_u,jspins,lty,n_mmp_new)
    END IF

    ! Write out n_mmp_new to out.xml file

    CALL openXMLElementNoAttributes('ldaUDensityMatrix')
    DO jsp = 1, jspins
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
                                         attributes,n_mmp_new(-l:l,-l:l,atoms%n_u,jsp))
       END DO
    END DO
    CALL closeXMLElement('ldaUDensityMatrix')

    !
    ! check for LDA+U and open density-matrix - file
    !
    INQUIRE (file='n_mmp_mat',exist=n_exist)
    OPEN (69,file='n_mmp_mat',status='unknown',form='formatted')


    IF (n_exist) THEN
       ALLOCATE (n_mmp_old(-3:3,-3:3,atoms%n_u,jspins))
       ALLOCATE (    n_mmp(-3:3,-3:3,atoms%n_u,jspins))
       READ (69,9000) n_mmp_old(:,:,:,:)

       READ (69,'(2(6x,f5.3))',IOSTAT=iofl) alpha,spinf
       IF ( iofl == 0 ) THEN
          !
          ! mix here straight with given mixing factors 
          !
          REWIND (69)
          sum1 = 0.0
          IF (jspins.EQ.1) THEN
             DO i_u = 1, atoms%n_u
                DO j = -3,3
                   DO k = -3,3
                      sum1 = sum1 + ABS(n_mmp_new(k,j,i_u,1) - n_mmp_old(k,j,i_u,1))
                      n_mmp(k,j,i_u,1) = alpha * n_mmp_new(k,j,i_u,1) + (1.0-alpha) * n_mmp_old(k,j,i_u,1)
                   END DO
                END DO
             END DO
             WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
          ELSE
             sum2 = 0.0
             gam = 0.5 * alpha * (1.0 + spinf)
             del = 0.5 * alpha * (1.0 - spinf)
             DO i_u = 1,atoms%n_u
                DO j = -3,3
                   DO k = -3,3
                      sum1 = sum1 + ABS(n_mmp_new(k,j,i_u,1) - n_mmp_old(k,j,i_u,1))
                      sum2 = sum2 + ABS(n_mmp_new(k,j,i_u,2) - n_mmp_old(k,j,i_u,2))

                      n_mmp(k,j,i_u,1) =       gam * n_mmp_new(k,j,i_u,1) + &
                                         (1.0-gam) * n_mmp_old(k,j,i_u,1) + &
                                               del * n_mmp_new(k,j,i_u,2) - &
                                               del * n_mmp_old(k,j,i_u,2)

                      n_mmp(k,j,i_u,2) =       gam * n_mmp_new(k,j,i_u,2) + &
                                         (1.0-gam) * n_mmp_old(k,j,i_u,2) + &
                                               del * n_mmp_new(k,j,i_u,1) - &
                                               del * n_mmp_old(k,j,i_u,1)
                   END DO
                END DO
             END DO
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
          ENDIF
          WRITE (69,9000) n_mmp
          WRITE (69,'(2(a6,f5.3))') 'alpha=',alpha,'spinf=',spinf

       ELSEIF (iofl > 0 ) THEN
          !
          ! read error ; stop
          !
          WRITE (6,*) 'ERROR READING mixing factors in n_mmp_mat'
          WRITE (6,'(2(a6,f5.3))') 'alpha=',alpha,'spinf=',spinf
          CALL juDFT_error("ERROR READING n_mmp_mat", calledby ="u_mix")
       ELSE
          !
          ! calculate distance and write new n_mmp to mix in broyden.F
          !
          sum1 = 0.0
          DO i_u = 1, atoms%n_u
             DO j = -3,3
                DO k = -3,3
                   sum1 = sum1 + ABS(n_mmp_new(k,j,i_u,1) - n_mmp_old(k,j,i_u,1))
                END DO
             END DO
          END DO
          IF (jspins.EQ.1) THEN
             WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
          ELSE
             sum2 = 0.0
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
             DO i_u = 1, atoms%n_u
                DO j = -3,3
                   DO k = -3,3
                      sum2 = sum2 + ABS(n_mmp_new(k,j,i_u,2) - n_mmp_old(k,j,i_u,2))
                   END DO
                END DO
             END DO
             DO j=-3,3
                WRITE(6,'(14f12.6)') (n_mmp_old(k,j,1,2),k=-3,3)
             END DO
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
             DO j=-3,3
                WRITE(6,'(14f12.6)') (n_mmp_new(k,j,1,2),k=-3,3)
             END DO
          END IF
          REWIND(69)
          WRITE (69,9000) n_mmp_old
          WRITE (69,9000) n_mmp_new
       END IF !  iofl == 0 

       DEALLOCATE ( n_mmp_old,n_mmp )
    ELSE
       !
       ! first time with lda+u; write new n_mmp  
       !
       WRITE (69,9000) n_mmp_new
       WRITE (69,'(2(a6,f5.3))') 'alpha=',0.05,'spinf=',1.0
    END IF

9000 FORMAT(7f20.13)

    CLOSE (69)
  END SUBROUTINE u_mix
END MODULE m_umix
