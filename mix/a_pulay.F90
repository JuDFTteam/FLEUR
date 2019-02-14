MODULE m_a_pulay
  USE m_juDFT
  PRIVATE
  REAL :: distance(3)
  PUBLIC a_pulay
CONTAINS
  SUBROUTINE a_pulay(alpha,fm,sm)
    USE m_pulay
    USE m_types_mat
    USE m_types_mixvector
    USE m_mixing_history
    IMPLICIT NONE
    REAL,INTENT(IN)                 :: alpha
    TYPE(t_mixvector),INTENT(IN)    :: fm(:)
    TYPE(t_mixvector),INTENT(INOUT) :: sm(:)


    REAL    :: fac(2)=(/1.0,0.5/)
    INTEGER,parameter :: parallel_images=2
    INTEGER,PARAMETER :: maxiter=4,mode=2
    INTEGER :: hlen,image,i,ii,local_length
    INTEGER :: local_hist(maxiter*parallel_images)

    hlen=SIZE(sm) !Number of results in history
    image=MOD(hlen-1,parallel_images)+1
    PRINT *,"HI:",hlen,image
    IF (hlen<parallel_images) THEN
       !create inital parallel images by simple mixing of first result
       sm(hlen)=sm(1)+(alpha*fac(image))*fm(1)
       RETURN
    ENDIF

    SELECT CASE(mode)
    CASE(1)
       local_length=0
       DO i=image,hlen,parallel_images
          local_length=local_length+1
          local_hist(local_length)=i
       ENDDO
       
       IF (local_length>maxiter) THEN !recreate to enable cross-image mixing
          local_length=hlen+1-parallel_images
          local_hist(1)=image
          local_hist(2:local_length)=(/(i,i=parallel_images+1,hlen)/)
       END IF
       sm(hlen)=simple_pulay(alpha*fac(image),fm(local_hist(:local_length)),sm(local_hist(:local_length)))
       IF (hlen==parallel_images*(maxiter+1)) CALL mixing_history_limit(parallel_images)
    case(2)
       local_length=hlen-image
       local_hist(1:local_length)=(/(i,i=1,local_length)/)
       local_length=local_length+1
       local_hist(local_length)=hlen
       sm(hlen)=simple_pulay(alpha*fac(image),fm(local_hist(:local_length)),sm(local_hist(:local_length)))
       IF (hlen==parallel_images*(maxiter)) CALL mixing_history_limit(parallel_images)
    END SELECT

    PRINT *,"P:",local_hist(:local_length)
    
    
  END SUBROUTINE a_pulay

  FUNCTION simple_pulay(alpha,fm,sm)RESULT(sm_out)
    USE m_types_mixvector
    IMPLICIT NONE
    REAL,INTENT(IN)                 :: alpha
    TYPE(t_mixvector),INTENT(IN)    :: fm(:)
    TYPE(t_mixvector),INTENT(IN)    :: sm(:)
    TYPE(t_mixvector)               :: sm_out
    
    
    TYPE(t_mixvector)   :: mdf
    REAL,ALLOCATABLE    :: a(:,:),b(:),work(:)
    INTEGER,ALLOCATABLE :: ipiv(:)
    INTEGER             :: n,nn,info,lwork,hlen
    
    

    hlen=SIZE(sm)
    IF (hlen==1) THEN
       sm_out=sm(1)+alpha*fm(1)
       RETURN
    ENDIF

    ALLOCATE(a(hlen+1,hlen+1),b(hlen+1),ipiv(hlen+1),work((hlen+1)**2))
    DO n=1,hlen
       mdf=fm(n)%apply_metric()
       DO nn=1,n
          a(n,nn)=mdf.dot.fm(nn)
          a(nn,n)=a(n,nn)
       ENDDO
    ENDDO

    a(:,hlen+1)=1.0
    a(hlen+1,:)=1.0
    a(hlen+1,hlen+1)=0.0
    b=0.0;b(hlen+1)=1.
    CALL  DSYSV( 'Upper', hlen+1, 1, a, SIZE(a,1), ipiv, b, SIZE(b,1), work, SIZE(work), INFO )
    sm_out=0.0*mdf
    DO n=1,hlen
       sm_out=sm_out+b(n)*(sm(n)+alpha*fm(n))
    ENDDO
  
  END FUNCTION SIMPLE_PULAY

END MODULE m_a_pulay
