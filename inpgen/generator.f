      MODULE m_generator
      use m_juDFT
!***********************************************************************
!     determines a set of generators for the group defined by the
!     nops matrices mrot. the set is not unique.
!***********************************************************************
      CONTAINS
      SUBROUTINE generator(
     >                     nops,mrot,tau,outfh,errfh)

      USE m_closure, ONLY : close_pt

      IMPLICIT NONE

!==> Arguments
      INTEGER, INTENT (IN) :: nops,outfh,errfh
      INTEGER, INTENT (IN) :: mrot(3,3,nops)
      REAL,    INTENT (IN) :: tau(3,nops)

!==> Locals
      INTEGER i,j,k,n,ngen,nops_sub,ncyl,ncl(1)
      INTEGER igenerator(nops)
      INTEGER mtable(nops,nops),nrow(nops)
      INTEGER trace(nops),mdet(nops),mtrd(nops)
      INTEGER nfactor(nops)

!--->    for trivial cases, print out and end

      IF ( nops .eq. 1 ) THEN
         WRITE (outfh,'(//," Space group can be generated using the",
     &                " identity only (no generators)")')
         RETURN
      ENDIF

!--->    generate multiplication table
      CALL close_pt(
     >              nops,mrot,
     <              mtable)

      IF ( nops .eq. 2 ) then
        ngen = 1
        igenerator(1) = nops
        GOTO 200
      ENDIF

!--->    determine the trace and determinant of each operation
      DO n=1,nops
         trace(n) = mrot(1,1,n) + mrot(2,2,n) + mrot(3,3,n)
         mdet(n) =
     &    mrot(1,1,n)*(mrot(2,2,n)*mrot(3,3,n)-mrot(3,2,n)*mrot(2,3,n))
     &   +mrot(1,2,n)*(mrot(3,1,n)*mrot(2,3,n)-mrot(2,1,n)*mrot(3,3,n))
     &   +mrot(1,3,n)*(mrot(2,1,n)*mrot(3,2,n)-mrot(3,1,n)*mrot(2,2,n))
         mtrd(n) = trace(n)*mdet(n)
      ENDDO

      ngen = 0
      nops_sub = nops
      nfactor(1:nops) = 1

!--->    check whether inversion exits (tr = -3); if so, a generator
      DO n=1,nops
         IF ( trace(n) == -3 ) THEN
            ngen = ngen + 1
            igenerator(ngen) = n
            nops_sub = nops_sub/2
            WHERE ( mdet == -1 ) nfactor = 0  ! get factor group
         ENDIF
      ENDDO

!--->    look for first 6, bar{6} operation
      IF ( mod(nops_sub,6) == 0 ) THEN
         DO n=1,nops
            IF ( (nfactor(n)==1) .and. (abs(trace(n))==2) ) THEN
               ngen = ngen + 1
               igenerator(ngen) = n
               nops_sub = nops_sub/6
               EXIT
            ENDIF
         ENDDO
      ENDIF

!--->    get first 3-fold axis not yet included
      IF ( mod(nops_sub,3) == 0 ) THEN
         DO n=1,nops
            IF ( (nfactor(n)==1) .and. (trace(n)==0) ) THEN
               ngen = ngen + 1
               igenerator(ngen) = n
               nops_sub = nops_sub/3
               EXIT
            ENDIF
         ENDDO
      ENDIF

!--->    get first 4-fold axis not yet included;
!--->    for cubic case, do not use, use 2 and m instead
      IF ( ( mod(nops_sub,4) == 0 ) .and. ( any(mtrd == 1) ) ) then
         IF ( nops < 24 ) THEN
            DO n=1,nops
               IF ( (nfactor(n)==1) .and. (mtrd(n)==1) ) THEN
                  ngen = ngen + 1
                  igenerator(ngen) = n
                  nops_sub = nops_sub/4
                  EXIT
               ENDIF
            ENDDO
         ELSE   ! cubic case: give the C_4 in terms of C_3 and m
            k = igenerator(ngen)  ! the 3-fold generator
            nrow(:) = mtable(:,k)
            ncyl = 0
            DO n=1,nops           ! loop over the 3-fold rotations
               IF ( (nfactor(n)==1) .and. (trace(n)==0) ) THEN
!-->                              ! find j such that n=mtable(j,k)
                   ncl = maxloc( nrow , MASK = nrow .eq. n )
                   if( mtrd(ncl(1)) == -1 ) then
                     ncyl = ncyl + 1
                     ngen = ngen + 1
                     igenerator(ngen) = ncl(1)
                     nops_sub = nops_sub/2
                     IF ( ncyl .ge. 2) EXIT
                   ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF

!--->    generate group determined by generators up to now
      nfactor = 0
      nfactor(1) = 1

      DO i=1,ngen
         k=1
         nrow = nfactor
         DO   ! multiply
            k = mtable( igenerator(i), k )
            IF ( k == 1 ) EXIT
            DO j=1,nops
               IF ( nfactor(j) == 1 ) THEN
                  nrow( mtable(k,j) ) = 1
               ENDIF
            ENDDO

         ENDDO
         nfactor = nrow
      ENDDO

!--->    at this point, only operations left to consider are 2, m
!--->    nfactor contains the elements of the factor group to this point

      DO
         IF ( nops_sub == 1 ) EXIT
!--->       check for first 2 or m operation not yet included
         nrow = nfactor
         DO n=1,nops
            IF ( (nfactor(n)==0) .and. (mtrd(n) == -1) ) THEN
               ngen = ngen + 1
               igenerator(ngen) = n
               nops_sub = nops_sub/2
               DO j=1,nops             ! generate larger group
                  if( nfactor(j) == 1 ) nrow( mtable(n,j) ) = 1
               ENDDO
               EXIT
            ENDIF
  100       CONTINUE
         ENDDO
         nfactor = nrow
      ENDDO

 200  CONTINUE

!--->    output

      IF (ngen == 1 ) THEN
         WRITE (outfh,'(//," Space group can be generated using",i2,
     &                " generator: ",10i4)') ngen,igenerator(1:ngen)
         WRITE (outfh,'(/,"   generator  (in lattice coordinates):")')
      ELSE
         WRITE (outfh,'(//," Space group can be generated using",i2,
     &                " generators:",10i4)') ngen,igenerator(1:ngen)
         WRITE (outfh,'(/,"   generators (in lattice coordinates):")')
      ENDIF
      WRITE (outfh,'(/,"&gen",i10)') ngen
      DO n=1,ngen
         WRITE (outfh,*)
         WRITE (outfh,'(3i5,5x,f10.5)') 
     &  ( ( mrot(i,j,igenerator(n)),j=1,3 ),tau(i,igenerator(n)),i=1,3 )
      ENDDO
      WRITE (outfh,'(/,"/ ! end generators",/)')

!--->    test to make sure the generators do generate the group
      nfactor = 0
      nfactor(1) = 1
      DO n=1,ngen
         k=1
         nrow = nfactor
         DO   ! multiply
            k = mtable( igenerator(n), k )
            if( k == 1 ) exit
            DO j=1,nops
               IF ( nfactor(j) == 1 ) THEN
                  i = mtable(k,j)
                  IF ( nrow(i) == 0 ) THEN
                     nrow( i ) = 1
                  ELSE
                     WRITE (errfh,'("generators: Error: multiple ops")')
                     CALL juDFT_error("Multiple ops.",calledby
     +                    ="generator")
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         nfactor = nrow
      ENDDO
      IF ( any(nfactor == 0) ) THEN
         WRITE (errfh,'("generators: Error: Group not covered")')
         CALL juDFT_error("Group not covered.",calledby ="generator")
      ENDIF

      END SUBROUTINE generator
      END MODULE m_generator
