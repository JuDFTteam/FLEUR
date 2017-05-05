      MODULE m_rwsymfile
      use m_juDFT
!----------------------------------------------------------------------!
!     writes spacegroup operations                                     ! 
!     and                                                              |  
!     rec. lattice vectors (for external k-point generator)            !
!----------------------------------------------------------------------!
      CONTAINS
      SUBROUTINE rw_symfile(
     >                      rw,symfh,symfn,nopd,bmat,
     X                      mrot,tau,nop,nop2,symor)

      IMPLICIT NONE

!===> Arguments
      CHARACTER(len=1), INTENT (IN)    :: rw
      CHARACTER(len=7), INTENT (IN)    :: symfn
      INTEGER,          INTENT (IN)    :: nopd,symfh
      REAL,             INTENT (IN)    :: bmat(3,3)
      INTEGER,          INTENT (INOUT) :: nop,nop2
      INTEGER,          INTENT (INOUT) :: mrot(3,3,nopd)
      REAL,             INTENT (INOUT) :: tau(3,nopd)
      LOGICAL,          INTENT (INOUT) :: symor

!===> Variables
      INTEGER i,j,n,ios,no2,no3,gen1,isrt(nop)
      REAL    t,d
      LOGICAL ex,op,l_exist
      CHARACTER(len=3) :: type
      CHARACTER(len=7) :: sym2fn

      sym2fn = 'sym.out'

      IF ( SCAN(rw,'wW') > 0 ) THEN

!===> write symfile

        OPEN (symfh, file=sym2fn, status='unknown', err=911, iostat=ios)
        WRITE (symfh,*) nop,nop2,symor,'    ! nop,nop2,symor '
        
        IF (nop == 2*nop2) THEN  ! film-calculation
          i = 1 ; j = nop2 + 1
          DO n = 1, nop
            IF (mrot(3,3,n) == 1) THEN
              isrt(n) = i ; i = i + 1
            ELSE
              isrt(n) = j ; j = j + 1
            ENDIF
          ENDDO
        ELSE
          DO n = 1, nop
            isrt(n) = n
          ENDDO
        ENDIF

        DO n = 1, nop
           WRITE (symfh,'(a1,i3)') '!', n
           WRITE (symfh,'(3i5,5x,f10.5)')
     &           ((mrot(i,j,isrt(n)),j=1,3),tau(i,isrt(n)),i=1,3)
        ENDDO

!        WRITE (symfh,*) '! reciprocal lattice vectors'
!        WRITE (symfh,'(3f25.15)') ((bmat(i,j),j=1,3),i=1,3)

      ELSEIF ( SCAN(rw,'rR') > 0 ) THEN

!===> read symfile
        INQUIRE(FILE=TRIM(ADJUSTL(symfn)),EXIST=l_exist)
        IF(.NOT.l_exist) THEN
           CALL juDFT_error("File "//TRIM(ADJUSTL(symfn))//
     +                      " is missing.",calledby="rw_symfile")
        END IF
        OPEN (symfh, file=trim(symfn),status='old',err=911,iostat=ios)
        READ (symfh,*) nop,nop2,symor
        IF (symfn.EQ.'sym.out') THEN
          gen1 = 0
        ELSEIF (trim(symfn).EQ.'sym') THEN
          gen1 = 1
        ELSE
           CALL juDFT_error("symfn should be sym or sym.out",calledby
     +          ="rw_symfile")
        ENDIF
        DO n = 1 + gen1, nop + gen1
          READ (symfh,*)
          READ (symfh,*) 
     &         ((mrot(i,j,n),j=1,3),tau(i,n),i=1,3)
        ENDDO
        IF (symor) THEN
          DO n=1,nop
            t= tau(1,n)**2 + tau(2,n)**2 + tau(3,n)**2
            IF (t > 1.e-8)  CALL juDFT_error("not symmorphic",calledby
     +           ="rw_symfile")
          ENDDO
        ELSE
          DO n=1,nop
            DO i = 1,3
             IF (ABS(tau(i,n)-0.33333) < 0.00001) THEN
               tau(i,n) = 1./3.
             ENDIF
             IF (ABS(tau(i,n)+0.33333) < 0.00001) THEN
               tau(i,n) = -1./3.
             ENDIF
             IF (ABS(tau(i,n)-0.66667) < 0.00001) THEN
               tau(i,n) = 2./3.
             ENDIF
             IF (ABS(tau(i,n)+0.66667) < 0.00001) THEN
               tau(i,n) = -2./3.
             ENDIF
             IF (ABS(tau(i,n)) > 0.00001) THEN
             IF (ABS(ABS(tau(i,n))-0.5) > 0.00001) THEN
                CALL juDFT_error("complex :: phases not fully tested!"
     +               ,calledby ="rw_symfile")
             ENDIF
             ENDIF
            ENDDO
          ENDDO
        ENDIF

        DO n = 1,nop
!
! Determine the kind of symmetry operation we have here
!
          d = det(mrot(:,:,n))
          t =  mrot(1,1,n) + mrot(2,2,n) + mrot(3,3,n)

          IF (d.EQ.-1) THEN
            type = 'm  '
            IF (t.EQ.-3) type = 'I  '
          ELSEIF (d.EQ.1) THEN
            IF (t.EQ.-1) type = 'c_2'
            IF (t.EQ. 0) type = 'c_3'
            IF (t.EQ. 1) type = 'c_4'
            IF (t.EQ. 2) type = 'c_6'
            IF (t.EQ. 3) type = 'E  '
          ELSE
             CALL juDFT_error("determinant =/= +/- 1",calledby
     +            ="rw_symfile")
          ENDIF
 
          WRITE (6,FMT=8020) n, type
 8020     FORMAT (/,1x,i3,' : ',a3)
          DO i = 1,3
             WRITE (6,FMT=8030) (mrot(i,j,n),j=1,3),tau(i,n)
          ENDDO
 8030     FORMAT (5x,3i3,3x,f4.1)
        ENDDO

      ELSE
         CALL juDFT_error("ERROR! rw_symfile #1",calledby="rw_symfile")
      ENDIF
      CLOSE (symfh)
      RETURN

! === errors
 911  CONTINUE
      WRITE(*,*) 'Error in inquire. IOS=',ios
 912  CONTINUE
      WRITE(*,*) 'Error in open. IOS=',ios
       CALL juDFT_error("i/o ERROR",calledby="rw_symfile")

      END SUBROUTINE rw_symfile

!--------------------------------------------------------------------
      INTEGER FUNCTION det(m)
        INTEGER m(3,3)
        det = m(1,1)*m(2,2)*m(3,3) + m(1,2)*m(2,3)*m(3,1) +
     +        m(2,1)*m(3,2)*m(1,3) - m(1,3)*m(2,2)*m(3,1) -
     +        m(2,3)*m(3,2)*m(1,1) - m(2,1)*m(1,2)*m(3,3)
      END FUNCTION det
!--------------------------------------------------------------------


      END MODULE m_rwsymfile
