MODULE m_broyden
  USE m_juDFT
  !################################################################
  !     IMIX = 3 : BROYDEN'S FIRST METHOD
  !     IMIX = 5 : BROYDEN'S SECOND METHOD
  !     IMIX = 7 : GENERALIZED ANDERSEN METHOD
  !     sm   : input charge density of iteration m
  !            afterwards update rho(m+1)
  !     fm   : output minus input charge density of iteration m
  !     sm1  : input charge density of iteration m-1
  !     fm1   : output minus inputcharge density of iteration m-1
  !################################################################
CONTAINS
  SUBROUTINE broyden(alpha,fm,sm)
    USE m_types
    USE m_types_mixvector
    IMPLICIT NONE

    real,INTENT(IN)                 :: alpha
    TYPE(t_mixvector),INTENT(IN)    :: fm(:)
    TYPE(t_mixvector),INTENT(INOUT) :: sm(:)

    ! Locals
    INTEGER           :: n,it,hlen
    REAL              :: fmvm,vmnorm
    REAL,ALLOCATABLE  :: am(:),dfivi(:)
    TYPE(t_mixvector) :: fm1,sm1,ui,um,vi,vm
    TYPE(t_mixvector),allocatable :: u_store(:),v_store(:)

    hlen=size(fm)
    IF (hlen<2) THEN !Do a simple mixing step
       sm(hlen)=sm(hlen)+alpha*fm(hlen)
       RETURN
    ENDIF

    ALLOCATE(u_store(hlen-2),v_store(hlen-2))
    do it=1,hlen-2
       call u_store(it)%alloc()
       call v_store(it)%alloc()
    enddo
     
    CALL fm1%alloc()
    CALL sm1%alloc()
    CALL ui%alloc()
    CALL um%alloc()
    CALL vi%alloc()
    CALL vm%alloc()
    
    ALLOCATE (am(hlen-1),dfivi(hlen-1))
    dfivi = 0.0
    am  = 0.0
    DO n=2,hlen
       sm1 = sm(n) - sm(n-1)
       fm1 = fm(n) - fm(n-1)
       !     |vi> = w|vi> 
       !     loop to generate um : um = sm1 + alpha*fm1 - \sum <fm1|w|vi> ui
       um = alpha * fm1 + sm1
       
       DO it = n-2,1,-1
          ui=u_store(it)
          vi=v_store(it)
          
          am(it) = vi.dot.fm1
          ! calculate um(:) = -am(it)*ui(:) + um(:)
          um=um-am(it)*ui
          WRITE(6,FMT='(5x,"<vi|w|Fm> for it",i2,5x,f10.6)')it,am(it) 
       END DO

       ! calculate vm = alpha*wfm1 -\sum <fm1|w|vi> <fi1|w|vi><vi|
       ! convolute fm1 with the metrik and store in vm
       vm=fm1%apply_metric()
       DO it = n-2,1,-1
          vi=v_store(it)
          ! calculate vm(:) = -am(it)*dfivi*vi(:) + vm
          vm=vm-am(it)*dfivi(it)*vi
       END DO

       vmnorm=fm1.dot.vm
       ! if (vmnorm.lt.tol_10) stop

       ! calculate vm(:) = (1.0/vmnorm)*vm(:)
       vm=(1.0/vmnorm)*vm
     
       ! save dfivi(mit) for next iteration
       dfivi(n-1) = vmnorm
       IF (n<hlen) u_store(n-1)=um
       IF (n<hlen) v_store(n-1)=vm
    enddo
    ! update rho(m+1)
    ! calculate <fm|w|vm>
    fmvm = vm.dot.fm(hlen)
    ! calculate sm(:) = (1.0-fmvm)*um(:) + sm
    sm(hlen)=sm(hlen)+(1.0-fmvm)*um
 
  END SUBROUTINE broyden
END MODULE m_broyden
