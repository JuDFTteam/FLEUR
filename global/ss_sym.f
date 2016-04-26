      MODULE m_sssym
!-----------------------------------------------------------------------!
! tests the compatibility of the symmetry elements with the axis (q)    !
! defining a spin-sprial                                          gb`02 !
!-----------------------------------------------------------------------!
      CONTAINS
      SUBROUTINE ss_sym(
     >                  nop,mrot,qss,
     <                  error)
      
      USE m_constants, ONLY : pimach
      IMPLICIT NONE

      INTEGER, INTENT (IN)  :: nop, mrot(3,3,nop)
      REAL,    INTENT (IN)  :: qss(3)
      LOGICAL, iNTENT (OUT) :: error(nop)

      INTEGER iop,i,j
      REAL    qn,test
      REAL    q1(3),rrot(3,3)

!
! --> loop over symmetry elements
!
      error(:) = .false.
      qn = ( qss(1)**2 + qss(2)**2 + qss(3)**2 )

      IF (qn.LT.0.0000001) THEN
        WRITE(*,*) 'qss = 0 ; not a spin-spiral!'
        RETURN
      ELSE
        qn = 1.0 / qn
      ENDIF

      DO iop = 1, nop

        DO i=1,3
           DO j=1,3
              rrot(i,j)= REAL(mrot(i,j,iop))
           ENDDO
        ENDDO
!
! ----> rotate qss by symmetry element and form the dot-product
!       with unrotated vector (q1 . qss) 
!
!        CALL cotra3(qss,q1,rrot)
         q1=matmul(qss,rrot)  
!
! ----> if qss is unchanged, accept this symmetry element
!
        test = (qss(1)-q1(1))**2+(qss(2)-q1(2))**2+(qss(3)-q1(3))**2
        IF (abs(test).GT.0.0000001) THEN
          error(iop) = .true.
          WRITE (6,100) iop
        ENDIF
      ENDDO
 100  FORMAT ('Symmetry element no.',i3,' incompatible with axis qss')
      
      IF ( ANY(error(:)) ) THEN
        WRITE (6,*) 'symmetry incompatible with Spin Spiral Axis [qss]'
      ENDIF
      END SUBROUTINE ss_sym
      END MODULE m_sssym

