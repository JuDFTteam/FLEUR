      MODULE m_closure
      use m_juDFT
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Contains 3 subroutines that more or less check the closure:
!     closure :    checks whether the space group operations close
!     close_pt:    checks that the point group of the bravais 
!                  lattice closes
!     check_close: additionally calculate the multiplication table,
!                  inverse operations and also determines the type 
!                  of every operation                    mw99,gs00
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      CONTAINS
      SUBROUTINE closure(
     >                   mops,mrot,tau,nops,index_op,
     <                   lclose)

      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: mops           ! number of operations of the bravais lattice
      INTEGER, INTENT (IN)  :: nops           ! number of operations in space group
      INTEGER, INTENT (IN)  :: mrot(3,3,mops) ! refer to the operations of the 
      REAL,    INTENT (IN)  :: tau(3,mops)    ! bravais lattice
      INTEGER, INTENT (IN)  :: index_op(nops) ! mapping function between space group 
                                              ! op's and those of the bravais lattice
      LOGICAL, INTENT (OUT) :: lclose

      REAL    ttau(3),eps7
      INTEGER i,ii,j,jj,k,kk,mp(3,3),map(nops)

      eps7 = 1.0e-7

!--->    loop over all operations
      DO jj=1,nops
         j = index_op(jj)

         map(1:nops) = 0

!--->    multiply {R_j|t_j}{R_i|t_i}
         DO ii=1,nops
            i = index_op(ii)
            mp = matmul( mrot(:,:,j) , mrot(:,:,i) )
            ttau = tau(:,j) + matmul( mrot(:,:,j) , tau(:,i) )
            ttau = ttau - anint( ttau - eps7 )

!--->    determine which operation this is
            DO kk=1,nops
              k = index_op(kk)
              IF ( all( mp(:,:) == mrot(:,:,k) ) .AND.         
     &             all( abs( ttau(:)-tau(:,k) ) < eps7 ) ) THEN
                 IF ( map(ii) .eq. 0 ) THEN
                    map(ii) = kk
                 ELSE
                    write(6,*)'ERROR Closure: Multiplying ', jj,' with '
     &                        ,kk, ' and with ',map(ii)
                    write(6,*) 'yields the same matrix'
                    lclose = .false.
                    RETURN
                 ENDIF
              ENDIF
            ENDDO

            IF (map(ii).eq.0) THEN
               write(6,*)'ERROR Closure:',ii,' times',jj,' leaves group'
               lclose = .false.
               RETURN
            ENDIF
         ENDDO
      ENDDO

      lclose = .true.

      END SUBROUTINE closure
!*********************************************************************

      SUBROUTINE close_pt(
     >                    nops,mrot,
     <                    mtable)

      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nops,mrot(3,3,nops)
      INTEGER, INTENT (OUT) :: mtable(nops,nops)   ! table(i,j) = {R_i|0}{R_j|0}

      INTEGER              :: i,j,k,mp(3,3),map(nops)

!---> loop over all operations
      DO j=1,nops

         map(1:nops) = 0

!--->    multiply {R_j|0}{R_i|0}
         DO i=1,nops
            mp = matmul( mrot(:,:,j) , mrot(:,:,i) )

!--->       determine which operation this is
            DO k = 1, nops
              IF ( all( mp(:,:)==mrot(:,:,k) ) ) THEN
                 IF ( map(i) .eq. 0 ) THEN
                    map(i) = k
                 ELSE
                    WRITE (6,'(" Symmetry error : multiple ops")')
                    CALL juDFT_error("close_pt: Multiple ops (Bravais)"
     +                   ,calledby ="closure")
                 ENDIF
              ENDIF
            ENDDO

            IF (map(i).eq.0) THEN
               WRITE (6,'(" Group not closed (Bravais lattice)")')
               WRITE (6,'(" operation j=",i2,"  map=",12i4,:/,
     &                  (21x,12i4))')  j, map(1:nops)
               CALL juDFT_error("close_pt:Not closed",calledby="closure"
     +              )
            ENDIF
         ENDDo
         mtable(j,1:nops) = map(1:nops)
      ENDDO

      END SUBROUTINE close_pt
!*********************************************************************

      SUBROUTINE check_close(
     >                       nops,mrot,tau,
     <                       multtab,inv_op,optype)

      IMPLICIT NONE

!===> Arguments
      INTEGER, INTENT (IN)  :: nops
      INTEGER, INTENT (IN)  :: mrot(3,3,nops)
      REAL,    INTENT (IN)  :: tau(3,nops)
      INTEGER, INTENT (OUT) :: inv_op(nops)
      INTEGER, INTENT (OUT) :: multtab(nops,nops)
      INTEGER, INTENT (OUT) :: optype(nops)

!===> Local Variables
      REAL    ttau(3)
      INTEGER i,j,n,k,mp(3,3),mdet,mtr

      REAL,    PARAMETER :: eps=1.0e-7
      INTEGER, PARAMETER :: cops(-1:3)=(/ 2, 3, 4, 6, 1 /)

      inv_op(1:nops) = 0

      multtab = 0

!--->    loop over all operations
      DO j=1,nops

!--->    multiply {R_j|t_j}{R_i|t_i}
         DO i=1,nops
            mp = matmul( mrot(:,:,j) , mrot(:,:,i) )
            ttau = tau(:,j) + matmul( mrot(:,:,j) , tau(:,i) )
            ttau = ttau - anint( ttau - eps )

!--->       determine which operation this is
            DO k=1,nops
              IF ( all( mp(:,:) == mrot(:,:,k) ) .and.
     &             all( abs( ttau(:)-tau(:,k) ) < eps ) ) THEN
                 IF ( multtab(j,i) .eq. 0 ) THEN
                    multtab(j,i) = k
                    IF (k .eq. 1) inv_op(j)=i
                 ELSE
                    WRITE(6,'(" Symmetry error: multiple ops")')
                    CALL juDFT_error("check_close: Multiple ops",
     +                   calledby ="closure")
                 ENDIF
              ENDIF
            ENDDO

            IF (multtab(j,i).eq.0) THEN
               WRITE (6,'(" Group not closed")')
               WRITE (6,'("  j , i =",2i4)') j,i
               CALL juDFT_error("check_close: Not closed",calledby
     +              ="closure")
            ENDIF
         ENDDO
      ENDDO

!--->      determine the type of each operation
      DO n = 1, nops
         mtr = mrot(1,1,n) + mrot(2,2,n) + mrot(3,3,n)
         mdet =
     &    mrot(1,1,n)*(mrot(2,2,n)*mrot(3,3,n)-mrot(3,2,n)*mrot(2,3,n))
     &   +mrot(1,2,n)*(mrot(3,1,n)*mrot(2,3,n)-mrot(2,1,n)*mrot(3,3,n))
     &   +mrot(1,3,n)*(mrot(2,1,n)*mrot(3,2,n)-mrot(3,1,n)*mrot(2,2,n))

         optype(n) = mdet*cops(mdet*mtr)

      ENDDO

      END SUBROUTINE check_close
      END MODULE m_closure
