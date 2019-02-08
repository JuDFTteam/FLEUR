MODULE m_pulay
  USE m_juDFT

CONTAINS
  SUBROUTINE pulay(alpha,fm,sm)
    USE m_types_mat
    USE m_types_mixvector
    IMPLICIT NONE
    REAL,INTENT(IN)                 :: alpha
    TYPE(t_mixvector),INTENT(IN)    :: fm(:)
    TYPE(t_mixvector),INTENT(INOUT) :: sm(:)

    ! Locals
    INTEGER           :: h_len,n,nn
    REAL,ALLOCATABLE  :: b(:)
    TYPE(t_mat)       :: a
    TYPE(t_mixvector),ALLOCATABLE :: df(:),ds(:)

    h_len=SIZE(fm)-1
    
    sm(h_len+1)=sm(h_len+1)+alpha*fm(h_len+1) !Simple mixing
    
    PRINT *,"len:",h_len
    DO n=1,h_len
       PRINT *,"D:",sm(n).dot.sm(n),fm(n).dot.fm(n)
    ENDDO
    IF (h_len==1) RETURN
    
    CALL a%alloc(.TRUE.,h_len,h_len)

    ALLOCATE(df(h_len),ds(h_len),b(h_len))

    DO n=1,h_len
       df(n)=fm(n+1)-fm(n)
       ds(n)=sm(n+1)-sm(n)
       b(n)=df(n).dot.fm(h_len+1)
       DO nn=1,n
          a%data_r(n,nn)=df(n).dot.df(nn)
          a%data_r(nn,n)=a%data_r(n,nn)
       ENDDO
    ENDDO
    
    call a%inverse()

    PRINT *,"b0:",b
    b=MATMUL(a%data_r,b)
    
    PRINT *,"b1:",b
    
    DO n=1,h_len
       sm(h_len+1)=sm(h_len+1)-b(n)*(ds(n)+alpha*df(n))
    ENDDO
    

 
  END SUBROUTINE pulay
END MODULE m_pulay
