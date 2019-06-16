MODULE m_bravaissymm
  use m_juDFT
  !********************************************************************
  !     determines the point group of the bravais lattice given the
  !     lattice vectors. the idea is to determine all the lattice
  !     vectors that have the same length as a_{1,2,3}, and then use
  !     these to determine the possible rotation matrices.
  !     these rotation matrices are in lattice coordinates.    mw 12-99
  !********************************************************************
CONTAINS
  SUBROUTINE bravais_symm(cell,nops,mrot)
    USE m_types_cell
    IMPLICIT NONE

    !==> Arguments
    TYPE(t_cell),INTENT(in) :: cell
    INTEGER, INTENT(OUT) :: nops, mrot(:,:,:)    ! point group operations

    !==> Locals
    REAL    amet(3,3),b1,b2,b3,d1,d2,d3,dmax,dt
    INTEGER i,k,k1,k2,k3,m1,m2,m3,n1,n2,n3
    INTEGER irot(3,3)
    
    INTEGER,PARAMETER::neig12=12! max. number of lattice vectors with same length
    ! (max occurs for close-packed fcc: 12)
    INTEGER lv1(3,neig12),lv2(3,neig12),lv3(3,neig12)

    REAL, PARAMETER :: eps=1.0e-9
 

    !---> distances for the lattice vectors
    d1 = cell%aamat(1,1)
    d2 = cell%aamat(2,2)
    d3 = cell%aamat(3,3)
    b1 = ( cell%bmat(1,1) )**2 + ( cell%bmat(1,2) )**2 + ( cell%bmat(1,3) )**2
    b2 = ( cell%bmat(2,1) )**2 + ( cell%bmat(2,2) )**2 + ( cell%bmat(2,3) )**2
    b3 = ( cell%bmat(3,1) )**2 + ( cell%bmat(3,2) )**2 + ( cell%bmat(3,3) )**2

    !---> determine the cutoffs along each direction a_i:
    dmax = max( d1,d2,d3)

    m1 = nint( dmax * b1 )
    m2 = nint( dmax * b2 )
    m3 = nint( dmax * b3 )

    !---->loop over all possible lattice vectors to find those with the
    !---->length, i.e., ones that could be rotations
    n1 = 1
    n2 = 1
    n3 = 1

    lv1(1:3,1) = (/ 1,0,0 /)
    lv2(1:3,1) = (/ 0,1,0 /)
    lv3(1:3,1) = (/ 0,0,1 /)

    DO k3=-m3,m3
       DO k2=-m2,m2
          DO k1=-m1,m1

             dt = distance2(k1,k2,k3)

             !---->    check if the same length
             IF ( abs( dt - d1 ) < eps ) THEN
                IF (.not.( k1==1 .and. k2==0 .and. k3==0 ) ) THEN
                   n1 = n1+1
                   IF(n1>neig12)  CALL juDFT_error("n1>neig12", calledby ="bravais_symm")
                   lv1(1,n1) = k1
                   lv1(2,n1) = k2
                   lv1(3,n1) = k3
                ENDIF
             ENDIF

             IF ( abs( dt - d2 ) < eps ) THEN
                IF (.not.( k1==0 .and. k2==1 .and. k3==0 ) ) THEN
                   n2 = n2+1
                   IF(n2>neig12)  CALL juDFT_error("n2>neig12",calledby="bravais_symm")
                   lv2(1,n2) = k1
                   lv2(2,n2) = k2
                   lv2(3,n2) = k3
                ENDIF
             ENDIF

             IF ( abs( dt - d3 ) < eps ) THEN
                IF (.not.( k1==0 .and. k2==0 .and. k3==1 ) ) THEN
                   n3 = n3+1
                   IF(n3>neig12)  CALL juDFT_error("n3>neig12",calledby="bravais_symm")
                   lv3(1,n3) = k1
                   lv3(2,n3) = k2
                   lv3(3,n3) = k3
                ENDIF
             ENDIF

          ENDDO
       ENDDO
    ENDDO

    !---> the possible rotation matrices are given by the matrix of
    !---> column vectors of lv_{1,2,3}
    nops = 0
    DO k3 = 1,n3
       DO k2 = 1,n2
          DO k1 = 1,n1

             !--->          check whether determinant is +/-1 (needs to be for rotation)
             IF ( abs(mdet(k1,k2,k3)) .NE. 1 ) CYCLE

             !--->          check whether this maintains lengths correctly
             !--->          if M is the metric, then must have R^T M R = M 
             irot = reshape( (/ lv1(:,k1),lv2(:,k2),lv3(:,k3) /) , (/ 3,3 /) )
             IF ( any( abs(matmul( transpose(irot), matmul(cell%aamat,irot) ) - cell%aamat) > eps ) ) CYCLE

             nops = nops + 1
             IF ( nops > SIZE(mrot,3) )  CALL juDFT_error("nop > size(mrot)", calledby="bravais_symm")
             mrot(:,:,nops) = irot

          ENDDO
       ENDDO
    ENDDO

    WRITE (6,'(//," Point group of the Bravais lattice has ",i2," operations")') nops

    RETURN

  CONTAINS   ! INTERNAL routines

    REAL FUNCTION distance2(l1,l2,l3)
      !*********************************************************************
      !     calculates the magnitude square for a vector (l1,l2,l3) given in
      !     lattice units
      !*********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: l1,l2,l3

      distance2 = l1*(l1*cell%aamat(1,1) + 2*l2*cell%aamat(2,1)) + l2*(l2*cell%aamat(2,2) + 2*l3*cell%aamat(3,2)) + l3*(l3*cell%aamat(3,3) + 2*l1*cell%aamat(1,3))

      RETURN
    END FUNCTION distance2

    INTEGER FUNCTION mdet(k1,k2,k3)
      !*********************************************************************
      !     determines the determinant for possible rotation matrix
      !     ( lv1(:,k1) ; lv2(:,k2) ; lv3(:,k3) )
      !*********************************************************************
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: k1,k2,k3

      mdet = lv1(1,k1)*( lv2(2,k2)*lv3(3,k3) - lv2(3,k2)*lv3(2,k3) ) + lv1(2,k1)*( lv2(3,k2)*lv3(1,k3) - lv2(1,k2)*lv3(3,k3) ) + lv1(3,k1)*( lv2(1,k2)*lv3(2,k3) - lv2(2,k2)*lv3(1,k3) )

      RETURN
    END FUNCTION mdet

  END SUBROUTINE bravais_symm
END MODULE m_bravaissymm
