MODULE m_a_pulay
  USE m_juDFT
  PRIVATE
  REAL :: distance(3)
  REAL :: alpha(4)=(/0.0,0.0,0.0,0.0/)
  
  PUBLIC a_pulay
CONTAINS
  SUBROUTINE a_pulay(alpha_in,fm,sm)
    USE m_pulay
    USE m_types_mat
    USE m_types_mixvector
    IMPLICIT NONE
    REAL,INTENT(IN)                 :: alpha_in
    TYPE(t_mixvector),INTENT(IN)    :: fm(:)
    TYPE(t_mixvector),INTENT(INOUT) :: sm(:)


    TYPE(t_mixvector) :: mdf
    TYPE(t_mat)       :: a
    INTEGER           :: n,nn
    IF (alpha(1)==0.0) alpha=(/1.,2.,0.5,1./)*alpha_in

   
       
    SELECT CASE(SIZE(sm))
    CASE(1)
       mdf=fm(1)%apply_metric()
       sm(1)=sm(1)+alpha(1)*fm(1) !Simple mixing
    CASE(2)
       mdf=fm(2)%apply_metric()
       distance(1)=fm(2).dot.mdf
       sm(2)=sm(1)+alpha(2)*fm(1) !Simple mixing with double alpha
    CASE(3)
       mdf=fm(3)%apply_metric()
       distance(2)=fm(3).dot.mdf
       sm(3)=sm(1)+alpha(3)*fm(1) !Simple mixing with half alpha
    CASE(4)
       !Find best distance
       mdf=fm(4)%apply_metric()
       distance(3)=fm(4).dot.mdf
       !IF (distance(2)==MINVAL(distance)) alpha(4)=2*alpha(1)
       !IF (distance(3)==MINVAL(distance)) alpha(4)=0.5*alpha(1)
       !WRITE(*,'(a,3f9.4)'),"A-Pulay:",distance
       !WRITE(*,*) alpha(1)
       CALL pulay(alpha_in,fm,sm,0)
!!$       CALL a%alloc(.TRUE.,5,5)
!!$       DO n=1,4
!!$          mdf=fm(n)%apply_metric()
!!$          DO nn=1,n
!!$             a%data_r(n,nn)=mdf.dot.fm(nn)
!!$             a%data_r(nn,n)=a%data_r(n,nn)
!!$          ENDDO
!!$       ENDDO
!!$       a%data_r(:,5)=1.0
!!$       a%data_r(5,:)=1.0
!!$       a%data_r(5,5)=.0
!!$       PRINT *,"a:",a%data_r(1,1)
!!$       PRINT *,"a:",a%data_r(2,2)
!!$       PRINT *,"a:",a%data_r(3,3)
!!$       PRINT *,"a:",a%data_r(4,4)
!!$
!!$       CALL a%inverse()
!!$       PRINT *,"p:",a%data_r(:4,5)
!!$       mdf=0.0*mdf
!!$       DO n=1,4
!!$          mdf=mdf+a%data_r(n,5)*(sm(n)+alpha(n)*fm(n))
!!$       ENDDO
!!$       sm(4)=mdf
    END SELECT
    
  END SUBROUTINE a_pulay
END MODULE m_a_pulay
