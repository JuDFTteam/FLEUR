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
  SUBROUTINE broyden(cell,stars,atoms,vacuum,sphhar,input,noco,oneD,sym,&
                     mmap,nmaph,mapmt,mapvac2,nmap,fm,mit,sm,lpot)

#include"cpp_double.h"

    USE m_metric
    USE m_types

    IMPLICIT NONE

    TYPE(t_oneD),INTENT(IN)    :: oneD
    TYPE(t_input),INTENT(IN)   :: input
    TYPE(t_vacuum),INTENT(IN)  :: vacuum
    TYPE(t_noco),INTENT(IN)    :: noco
    TYPE(t_sym),INTENT(IN)     :: sym
    TYPE(t_stars),INTENT(IN)   :: stars
    TYPE(t_cell),INTENT(IN)    :: cell
    TYPE(t_sphhar),INTENT(IN)  :: sphhar
    TYPE(t_atoms),INTENT(IN)   :: atoms

    ! Scalar Arguments
    INTEGER, INTENT (IN)        :: mmap,nmap
    INTEGER, INTENT (IN)        :: mapmt,mapvac2    
    INTEGER, INTENT (INOUT)     :: mit
    LOGICAL,OPTIONAL,INTENT(IN) :: lpot

    ! Array Arguments
    REAL,    INTENT (IN)    :: fm(nmap) 
    REAL,    INTENT (INOUT) :: sm(nmap)

    ! Local Scalars
    INTEGER         :: i,it,k,nit,npos,iread,nmaph
    REAL            :: bm,dfivi,fmvm,smnorm,vmnorm,alphan
    LOGICAL         :: l_pot
    REAL, PARAMETER :: one=1.0, zero=0.0

    ! Local Arrays
    REAL, ALLOCATABLE :: am(:)
    REAL, ALLOCATABLE :: fm1(:),sm1(:),ui(:),um(:),vi(:),vm(:)

    ! External Functions
    REAL CPP_BLAS_sdot
    EXTERNAL CPP_BLAS_sdot

    ! External Subroutines
    EXTERNAL CPP_BLAS_saxpy,CPP_BLAS_sscal

    dfivi = zero

    l_pot = .FALSE.
    IF (PRESENT(lpot)) l_pot = lpot

    ALLOCATE (fm1(mmap),sm1(mmap),ui(mmap),um(mmap),vi(mmap),vm(mmap))
    ALLOCATE (am(input%maxiter+1))

    fm1 = 0.0
    sm1 = 0.0
    ui  = 0.0
    um  = 0.0
    vi  = 0.0
    vm  = 0.0
    am  = 0.0

    IF (mit.NE.1) THEN
       ! load input charge density (sm1) and difference of 
       ! in and out charge densities (fm1) from previous iteration (m-1)

       REWIND 57
       READ (57) mit, alphan, (fm1(i),i=1,nmap), (sm1(i),i=1,nmap)
       IF (ABS(input%alpha-alphan) > 0.0001) THEN
          WRITE (6,*) 'mixing parameter has been changed; reset'
          WRITE (6,*) 'broyden algorithm or set alpha to',alphan
          CALL juDFT_error("mixing parameter (input) changed", calledby ="broyden")
       END IF

       ! generate F_m   - F_(m-1)  ... sm1
       !      and rho_m - rho_(m-1) .. fm1
       sm1(1:nmap) = sm(1:nmap) - sm1(1:nmap)
       fm1(1:nmap) = fm(1:nmap) - fm1(1:nmap)
    END IF

    ! save F_m and rho_m for next iteration
    REWIND 57
    nit = mit +1
    IF (nit > input%maxiter+1) nit = 1
    WRITE (57) nit,input%alpha,fm,sm

    IF (mit.EQ.1) THEN 
       !     update for rho for mit=1 is straight mixing
       !     sm = sm + alpha*fm
       CALL CPP_BLAS_saxpy(nmap,input%alpha,fm,1,sm,1)
    ELSE
       !     |vi> = w|vi> 
       !     loop to generate um : um = sm1 + alpha*fm1 - \sum <fm1|w|vi> ui
       DO k = 1, nmap
          um(k) = input%alpha * fm1(k) + sm1(k)
       END DO
       iread = MIN(mit-1,input%maxiter+1)
       DO it = 2, iread
          READ (59,rec=(it-1)*2-1) (ui(i),i=1,nmap)
          READ (59,rec=(it-1)*2)   (vi(i),i=1,nmap),dfivi
          am(it) = CPP_BLAS_sdot(nmap,vi,1,fm1,1)
          CALL CPP_BLAS_saxpy(nmap,-am(it),ui,1,um,1)
          WRITE(6,FMT='(5x,"<vi|w|Fm> for it",i2,5x,f10.6)')it,am(it) 
       END DO

       IF (input%imix.EQ.3) THEN
          !****************************************
          !        broyden's first method
          !****************************************

          ! convolute drho(m) with the metric: |fm1> = w|sm1>
          CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                      mmap,nmaph,mapmt,mapvac2,sm1,fm1,l_pot)

          ! calculate the norm of sm1 : <sm1|w|sm1>
          smnorm = CPP_BLAS_sdot(nmap,sm1,1,fm1,1)

          ! generate vm = alpha*sm1  - \sum <ui|w|sm1> vi
          vm(:) = input%alpha * fm1(:)
          DO it = 2,iread
             READ (59,rec=2*(it-1)-1) (ui(i),i=1,nmap)
             READ (59,rec=2*(it-1))   (vi(i),i=1,nmap), dfivi
             bm = CPP_BLAS_sdot(nmap,ui,1,fm1,1)
             CALL CPP_BLAS_saxpy(nmap,-bm,vi,1,vm,1)
             !write(6,FMT='(5x,"<ui|w|Fm> for it",i2,5x,f10.6)') it, bm 
          END DO

          ! complete evaluation of vm
          ! vmnorm = <um|w|sm1>-<sm1|w|sm1>
          vmnorm = CPP_BLAS_sdot(nmap,fm1,1,um,1) - smnorm
          ! if (vmnorm.lt.tol_10) stop

          CALL CPP_BLAS_sscal(nmap,one/vmnorm,vm,1)

       ELSE IF (input%imix.EQ.5) THEN
          !****************************************
          !      broyden's second method
          !****************************************

          ! multiply fm1 with metric matrix and store in vm:  w |fm1>  
          CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                      mmap,nmaph,mapmt,mapvac2,fm1,vm,l_pot)

          ! calculate the norm of fm1 and normalize vm it: vm = wfm1 / <fm1|w|fm1>
          vmnorm = one / CPP_BLAS_sdot(nmap,fm1,1,vm,1)
          CALL CPP_BLAS_sscal(nmap,vmnorm,vm,1)

       ELSE IF (input%imix.EQ.7) THEN
          !****************************************
          !      generalized anderson method
          !****************************************

          ! calculate vm = alpha*wfm1 -\sum <fm1|w|vi> <fi1|w|vi><vi|
          ! convolute fm1 with the metrik and store in vm
          CALL metric(cell,atoms,vacuum,sphhar,input,noco,stars,sym,oneD,&
                      mmap,nmaph,mapmt,mapvac2,fm1,vm,l_pot)

          DO it = 2,iread
             READ (59,rec=2*(it-1)-1) (ui(i),i=1,nmap)
             READ (59,rec=2*(it-1))   (vi(i),i=1,nmap), dfivi
             CALL CPP_BLAS_saxpy(nmap,-am(it)*dfivi,vi,1,vm,1)
          END DO

          vmnorm = CPP_BLAS_sdot(nmap,fm1,1,vm,1)
          ! if (vmnorm.lt.tol_10) stop

          CALL CPP_BLAS_sscal(nmap,one/vmnorm,vm,1)

          ! save dfivi(mit) for next iteration
          dfivi = vmnorm

       END IF

       ! write um,vm and dfivi on file broyd.?

       npos=mit-1
       IF (mit.GT.input%maxiter+1) npos = MOD(mit-2,input%maxiter)+1

       WRITE (59,rec=2*npos-1) (um(i),i=1,nmap)
       WRITE (59,rec=2*npos)   (vm(i),i=1,nmap), dfivi

       ! update rho(m+1)
       ! calculate <fm|w|vm>
       fmvm = CPP_BLAS_sdot(nmap,vm,1,fm,1)
       CALL CPP_BLAS_saxpy(nmap,one-fmvm,um,1,sm,1)
    END IF

    mit = mit + 1
    DEALLOCATE (fm1,sm1,ui,um,vi,vm,am)

  END SUBROUTINE broyden
END MODULE m_broyden
