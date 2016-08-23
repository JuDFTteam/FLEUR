      MODULE m_divi
      CONTAINS
      SUBROUTINE divi(
     >                nkpt,rltv,film,nop,nop2,
     <                n)

      IMPLICIT NONE

      INTEGER, INTENT (IN) :: nkpt      ! maximum number of k-points
      INTEGER, INTENT (IN) :: nop,nop2  ! number of symmetry operations (3D,2D)
      REAL,    INTENT (IN) :: rltv(3,3) ! reciprocal lattice axes
      LOGICAL, INTENT (IN) :: film      ! .true. -> 2D k-points
      INTEGER, INTENT(OUT) :: n(3)      ! # of k-points in each direction

      INTEGER ntes,i,ntot
      REAL    b(3)

      DO i = 1,3
        b(i) = sqrt( rltv(i,1)**2 + rltv(i,2)**2 + rltv(i,3)**2 )
      ENDDO

      IF (film) THEN
        ntot = nkpt * nop2
      ELSE
        ntot = nkpt * nop
      ENDIF

      IF (film) THEN
        ntes = ntot
  5     n(1) = nint( sqrt( ntes*b(1)/b(2) ) )
        n(2) = nint( n(1)*b(2)/b(1) )
        IF (n(1)*n(2).GT.ntot) THEN
          n(2) = nint(  sqrt( ntes*b(2)/b(1) ) )
          n(1) = nint( n(2)*b(1)/b(2) )
          IF (n(1)*n(2).GT.ntot) THEN
            ntes = ntes - 1
            GOTO 5
          ENDIF
        ENDIF
        n(3) = 0
        IF (nop2.GE.4) THEN
          n(1) = 2*NINT(n(1)/2.0+0.1)
          n(2) = 2*NINT(n(2)/2.0+0.1)
        ENDIF
        WRITE (*,*) n(1),n(2),n(1)*n(2)
      ELSE
        ntes = ntot
 10     n(1) = nint( (ntes*b(1)**2/(b(2)*b(3)))**(1./3.) ) 
        n(2) = nint( n(1)*b(2)/b(1) )
        n(3) = nint( n(1)*b(3)/b(1) )
        IF (n(1)*n(2)*n(3).GT.ntot) THEN
          n(2) = nint( (ntes*b(2)**2/(b(1)*b(3)))**(1./3.) )
          n(1) = nint( n(2)*b(1)/b(2) )
          n(3) = nint( n(2)*b(3)/b(2) )
          IF (n(1)*n(2)*n(3).GT.ntot) THEN
            n(3) = nint( (ntes*b(3)**2/(b(1)*b(2)))**(1./3.) )
            n(2) = nint( n(3)*b(2)/b(3) )
            n(1) = nint( n(3)*b(1)/b(3) )
            IF (n(1)*n(2)*n(3).GT.ntot) THEN
              ntes = ntes - 1
              GOTO 10
            ENDIF
          ENDIF
        ENDIF
        IF (nop.GE.8) THEN
          n(1) = MAX(2*NINT(n(1)/2.0-0.1),1)
          n(2) = MAX(2*NINT(n(2)/2.0-0.1),1)
          n(3) = MAX(2*NINT(n(3)/2.0-0.1),1)
        ENDIF
        WRITE (*,*) n(1),n(2),n(3),n(1)*n(2)*n(3)
      ENDIF

      RETURN
      END SUBROUTINE divi
      END MODULE m_divi
