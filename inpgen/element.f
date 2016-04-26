      MODULE m_element
      CONTAINS
!  ---------------------------------------------------------------------
      INTEGER FUNCTION z_namat(element)
      USE m_constants,ONLY: namat_const
      IMPLICIT NONE

      CHARACTER(len=2)  element
      CHARACTER(len=2)  ele
      INTEGER           adiff,n

      adiff = IACHAR('a') - IACHAR('A')

      ele = ADJUSTL( element ) ! takes care of single letter elements

      IF ( LLT(ele(2:2),'a') .AND. ele(2:2).NE.' ' ) 
     &     ele(2:2) = achar( iachar(ele(2:2))+adiff )
      IF ( LGE(ele(1:1),'a') ) 
     &     ele(1:1) = achar( iachar(ele(1:1))-adiff )

      z_namat = -1
      DO n = 0, size(namat_const)
        IF ( ele == namat_const(n) ) THEN
          z_namat = n
          EXIT
        ENDIF
      ENDDO

      END FUNCTION z_namat
!  ---------------------------------------------------------------------
      END MODULE m_element
