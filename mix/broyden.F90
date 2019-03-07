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
  SUBROUTINE broyden(input,fm,sm)
    USE m_types
    USE m_broyd_io
    USE m_types_mixvector
    IMPLICIT NONE

    TYPE(t_input),INTENT(IN)        :: input
    TYPE(t_mixvector),INTENT(IN)    :: fm
    TYPE(t_mixvector),INTENT(INOUT) :: sm

    
    
    ! Local Scalars
    INTEGER         :: i,it,k,nit,iread,nmaph, mit
    REAL            :: bm,dfivi,fmvm,smnorm,vmnorm,alphan
    LOGICAL         :: l_exist
    REAL, PARAMETER :: one=1.0, zero=0.0

    ! Local Arrays
    REAL,ALLOCATABLE  :: am(:)
    TYPE(t_mixvector) :: fm1,sm1,ui,um,vi,vm

    
    dfivi = zero

    CALL fm1%alloc()
    CALL sm1%alloc()
    CALL ui%alloc()
    CALL um%alloc()
    CALL vi%alloc()
    CALL vm%alloc()
    
    ALLOCATE (am(input%maxiter+1))
    am  = 0.0
    l_exist = initBroydenHistory(input,hybrid,nmap) ! returns true if there already exists a Broyden history
    IF (l_exist) THEN
       ! load input charge density (sm1) and difference of 
       ! in and out charge densities (fm1) from previous iteration (m-1)

       CALL readLastIterInAndDiffDen(mit,alphan,sm1,fm1)
       IF (ABS(input%alpha-alphan) > 0.0001) THEN
          WRITE (6,*) 'mixing parameter has been changed; reset'
          WRITE (6,*) 'broyden algorithm or set alpha to',alphan
          CALL juDFT_error("mixing parameter (input) changed", calledby ="broyden")
       END IF

       ! generate F_m   - F_(m-1)  ... sm1
       !      and rho_m - rho_(m-1) .. fm1
       sm1 = sm - sm1
       fm1 = fm - fm1
    ELSE
       mit=1
    END IF

    ! save F_m and rho_m for next iteration
    nit = mit +1
    IF (nit > input%maxiter+1) nit = 1
    CALL writeLastIterInAndDiffDen(hybrid,nmap,nit,input%alpha,sm,fm)

    IF (.NOT.l_exist) THEN 
       !     update for rho for mit=1 is straight mixing
       !     sm = sm + alpha*fm
       sm=sm+input%alpha*fm
    ELSE
       !     |vi> = w|vi> 
       !     loop to generate um : um = sm1 + alpha*fm1 - \sum <fm1|w|vi> ui
       um = input%alpha * fm1 + sm1
       
       iread = MIN(mit-1,input%maxiter+1)
       DO it = 2, iread
          CALL readUVec(input,hybrid,nmap,it-mit,mit,ui)
          CALL readVVec(input,hybrid,nmap,it-mit,mit,dfivi,vi)

          am(it) = vi.dot.fm
          ! calculate um(:) = -am(it)*ui(:) + um(:)
          um=um-am(it)*ui
          WRITE(6,FMT='(5x,"<vi|w|Fm> for it",i2,5x,f10.6)')it,am(it) 
       END DO

       IF (input%imix.EQ.7) THEN
          !****************************************
          !      generalized anderson method
          !****************************************

          ! calculate vm = alpha*wfm1 -\sum <fm1|w|vi> <fi1|w|vi><vi|
          ! convolute fm1 with the metrik and store in vm
          vm=fm1%apply_metric()
          DO it = 2,iread
             CALL readVVec(input,hybrid,nmap,it-mit,mit,dfivi,vi)
             ! calculate vm(:) = -am(it)*dfivi*vi(:) + vm
             vm=vm-am(it)*dfivi*vi
          END DO

          vmnorm=fm1.dot.vm
          ! if (vmnorm.lt.tol_10) stop

          ! calculate vm(:) = (1.0/vmnorm)*vm(:)
          vm=(1.0/vmnorm)*vm
     
          ! save dfivi(mit) for next iteration
          dfivi = vmnorm
       ELSE
          CALL judft_error("Only generalized Anderson implemented")
       END IF

       ! write um,vm and dfivi on file broyd.?

       CALL writeUVec(input,hybrid,nmap,mit,um)
       CALL writeVVec(input,hybrid,nmap,mit,dfivi,vm)

       ! update rho(m+1)
       ! calculate <fm|w|vm>
       fmvm = vm.dot.fm
       ! calculate sm(:) = (1.0-fmvm)*um(:) + sm
       sm=sm+(one-fmvm)*um
    END IF
 
  END SUBROUTINE broyden
END MODULE m_broyden
