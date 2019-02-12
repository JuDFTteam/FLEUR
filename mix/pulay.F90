MODULE m_pulay
  USE m_juDFT

CONTAINS
  SUBROUTINE pulay(alpha,fm,sm,simple_steps)
    USE m_types_mat
    USE m_types_mixvector
    USE m_mixing_history
    IMPLICIT NONE
    TYPE(t_mixvector),INTENT(IN)    :: fm(:)
    TYPE(t_mixvector),INTENT(INOUT) :: sm(:)
    INTEGER,INTENT(IN)              :: simple_steps
    REAL,INTENT(IN)                 :: alpha
    ! Locals
    INTEGER           :: h_len,n,nn
    REAL,ALLOCATABLE  :: b(:)
    TYPE(t_mat)       :: a
    TYPE(t_mixvector),ALLOCATABLE :: df(:),ds(:)
    TYPE(t_mixvector) ::mdf


    h_len=SIZE(fm)-1
    
    IF (h_len==0) THEN
       sm(h_len+1)=sm(h_len+1)+alpha*fm(h_len+1) !Simple mixing
       RETURN !No history present
    ENDIF
    IF (simple_steps>0) THEN
       IF (MOD(h_len,simple_steps).NE.0) THEN
          !Simple mixing step of periodic Pulay
          sm(h_len+1)=sm(h_len+1)+alpha*fm(h_len+1) 
          RETURN 
       ENDIF
    ENDIF
    CALL a%alloc(.TRUE.,h_len,h_len)

    ALLOCATE(df(h_len),ds(h_len),b(h_len))

    DO n=1,h_len
       df(n)=fm(n+1)-fm(n)
       ds(n)=sm(n+1)-sm(n)
       mdf=df(n).apply_metric()
       b(n)=mdf.dot.fm(h_len+1)
       DO nn=1,n
          a%data_r(n,nn)=mdf.dot.df(nn)
          a%data_r(nn,n)=a%data_r(n,nn)
       ENDDO
    ENDDO
    
    call a%inverse()

    b=MATMUL(a%data_r,b)
    sm(h_len+1)=sm(h_len+1)+alpha*fm(h_len+1) !Simple mixing
    DO n=1,h_len
       sm(h_len+1)=sm(h_len+1)-b(n)*(ds(n)+alpha*df(n))
    ENDDO
   
    IF (simple_steps>0) CALL mixing_history_limit(0)
  END SUBROUTINE pulay
END MODULE m_pulay
