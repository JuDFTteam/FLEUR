MODULE m_film_sym
  IMPLICIT NONE
  USE m_juDFT
CONTAINS
  SUBROUTINE film_sym(nop,mrot,error)
    INTEGER,INTENT(in)::nop,mrot(:,:,:)
    LOGICAL,INTENT(inout)::error(:)

    INTEGER:: n

    DO n=1,nop
       IF (ANY(mrot(1:2,3).NE.0).OR.ANY(mrot(3,1:2).NE.0)) error(n)=.TRUE.
    END DO
  END SUBROUTINE film_sym
END MODULE m_film_sym
