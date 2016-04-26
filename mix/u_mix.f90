MODULE m_umix
  USE m_juDFT
  !
  ! mix the old and new density matrix for the lda+U method
  !                                                 gb.2001
CONTAINS
  SUBROUTINE u_mix(atoms,&
       &                 jspins,n_mmp_new)

    USE m_types, ONLY : t_utype
    USE m_nmat_rot
    !
    ! ... Arguments
    USE m_types
    IMPLICIT NONE
    TYPE(t_atoms),INTENT(IN)   :: atoms
    INTEGER, INTENT (IN)       :: jspins 
    COMPLEX, INTENT (INOUT)    :: n_mmp_new(-3:3,-3:3,atoms%n_u,jspins)
    !
    ! ... Locals ...
    INTEGER n,j,k,iofl,l,itype,ios,lty(atoms%ntype)
    REAL alpha,spinf,gam,del,sum1,sum2,mix_u
    REAL    theta(atoms%ntype),phi(atoms%ntype),zero(atoms%ntype)
    LOGICAL n_exist
    COMPLEX,ALLOCATABLE :: n_mmp(:,:,:,:),n_mmp_old(:,:,:,:)
    !
    ! check for possible rotation of n_mmp
    !
    INQUIRE (file='n_mmp_rot',exist=n_exist)
    IF (n_exist) THEN
       OPEN (68,file='n_mmp_rot',status='old',form='formatted')
       n = 1
       DO itype = 1,atoms%ntype
          l = atoms%lda_u(itype)%l
          IF (l.GE.0) THEN
             READ(68,*,iostat=ios) theta(n),phi(n)
             IF (ios == 0) THEN
                lty(n) = l ; n = n + 1
             ELSE
                IF (n == 1)  CALL juDFT_error("ERROR reading n_mmp_rot"&
                     &                ,calledby ="u_mix")
                theta(n) = theta(n-1) ; phi(n) = phi(n-1)
                lty(n) = lty(n-1) ; n = n + 1
             ENDIF
          ENDIF
       ENDDO
       CLOSE (68)
       zero = 0.0
       CALL nmat_rot(zero,-theta,-phi,3,atoms%n_u,jspins,lty,n_mmp_new)
    ENDIF
    !
    ! check for LDA+U and open density-matrix - file
    !
    INQUIRE (file='n_mmp_mat',exist=n_exist)
    OPEN (69,file='n_mmp_mat',status='unknown',form='formatted')


    IF (n_exist) THEN
       ALLOCATE (  n_mmp_old(-3:3,-3:3,atoms%n_u,jspins) )
       ALLOCATE (      n_mmp(-3:3,-3:3,atoms%n_u,jspins) )
       READ (69,9000) n_mmp_old(:,:,:,:)

       READ (69,'(2(6x,f5.3))',IOSTAT=iofl) alpha,spinf
       IF ( iofl == 0 ) THEN
          !
          ! mix here straight with given mixing factors 
          !
          REWIND (69)
          sum1 = 0.0
          IF (jspins.EQ.1) THEN
             DO  n = 1,atoms%n_u
                DO j = -3,3
                   DO k = -3,3
                      sum1 = sum1 + ABS( n_mmp_new(k,j,n,1) - &
                           &                                  n_mmp_old(k,j,n,1) )
                      n_mmp(k,j,n,1) = alpha  * n_mmp_new(k,j,n,1) +&
                           &                          (1. - alpha) * n_mmp_old(k,j,n,1)
                   ENDDO
                ENDDO
             ENDDO
             WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
          ELSE
             sum2 = 0.0
             gam = 0.5 * alpha * ( 1.0 + spinf )
             del = 0.5 * alpha * ( 1.0 - spinf )
             DO  n = 1,atoms%n_u
                DO j = -3,3
                   DO k = -3,3
                      sum1 = sum1 + ABS( n_mmp_new(k,j,n,1) - &
                           &                                  n_mmp_old(k,j,n,1) )
                      sum2 = sum2 + ABS( n_mmp_new(k,j,n,2) - &
                           &                                  n_mmp_old(k,j,n,2) )
                      n_mmp(k,j,n,1) = gam  * n_mmp_new(k,j,n,1) +&
                           &                          (1. - gam) * n_mmp_old(k,j,n,1) +&
                           &                                del  * n_mmp_new(k,j,n,2) -&
                           &                                del  * n_mmp_old(k,j,n,2)
                      n_mmp(k,j,n,2) = gam  * n_mmp_new(k,j,n,2) +&
                           &                          (1. - gam) * n_mmp_old(k,j,n,2) +&
                           &                                del  * n_mmp_new(k,j,n,1) -&
                           &                                del  * n_mmp_old(k,j,n,1)
                   ENDDO
                ENDDO
             ENDDO
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
          CALL juDFT_error("ERROR READING n_mmp_mat",calledby ="u_mix")
       ELSE
          !
          ! calculate distance and write new n_mmp to mix in broyden.F
          !
          sum1 = 0.0
          DO  n = 1,atoms%n_u
             DO j = -3,3
                DO k = -3,3
                   sum1 = sum1 + ABS( n_mmp_new(k,j,n,1) -&
                        &                                n_mmp_old(k,j,n,1) )
                ENDDO
             ENDDO
          ENDDO
          IF (jspins.EQ.1) THEN
             WRITE (6,'(a16,f12.6)') 'n_mmp distance =',sum1
          ELSE
             sum2 = 0.0
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 1 =',sum1
             DO  n = 1,atoms%n_u
                DO j = -3,3
                   DO k = -3,3
                      sum2 = sum2 + ABS( n_mmp_new(k,j,n,2) -&
                           &                                  n_mmp_old(k,j,n,2) )
                   ENDDO
                ENDDO
             ENDDO
             DO j=-3,3
                WRITE(6,'(14f12.6)') (n_mmp_old(k,j,1,2),k=-3,3)
             ENDDO
             WRITE (6,'(a23,f12.6)') 'n_mmp distance spin 2 =',sum2
             DO j=-3,3
                WRITE(6,'(14f12.6)') (n_mmp_new(k,j,1,2),k=-3,3)
             ENDDO
          ENDIF
          REWIND(69)
          WRITE (69,9000) n_mmp_old
          WRITE (69,9000) n_mmp_new
       ENDIF !  iofl == 0 

       DEALLOCATE ( n_mmp_old,n_mmp )
    ELSE
       !
       ! first time with lda+u; write new n_mmp  
       !
       WRITE (69,9000) n_mmp_new
       WRITE (69,'(2(a6,f5.3))') 'alpha=',0.05,'spinf=',1.0
    ENDIF

9000 FORMAT(7f20.13)

    CLOSE (69)
  END SUBROUTINE u_mix
END MODULE m_umix
