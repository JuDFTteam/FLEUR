!--------------------------------------------------------------------------------
! Copyright (c) 2025 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_film_sym
  USE m_juDFT
  IMPLICIT NONE
CONTAINS
  SUBROUTINE film_sym(nop,mrot,error)
    INTEGER,INTENT(in)::nop,mrot(:,:,:)
    LOGICAL,INTENT(inout)::error(:)

    INTEGER:: n,nn

    outer:DO n=1,nop
       IF (ANY(mrot(1:2,3,n).NE.0).OR.ANY(mrot(3,1:2,n).NE.0)) error(n)=.TRUE.
       IF (mrot(3,3,n)==-1) THEN !check if 2d-sym without invs or z-refl. also exists
        DO nn=1,nop
          if (mrot(3,3,nn)==1.and.all(mrot(1:2,1:2,n)==mrot(1:2,1:2,nn))) cycle outer !z-refl OK
          if (mrot(3,3,nn)==1.and.all(-1*mrot(1:2,1:2,n)==mrot(1:2,1:2,nn))) cycle outer !invs ok
        ENDDO
        error(n)=.true.  
       ENDIF 
    ENDDO outer
  END SUBROUTINE film_sym
END MODULE m_film_sym
