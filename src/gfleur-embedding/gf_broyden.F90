!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_broyden
      USE m_constants, ONLY: oUnit
    use m_juDFT
    IMPLICIT NONE
!---------------------------------------------------------------        
!  modified version:                                                    
!  IO is encapulated here                                               
!                     Daniel Wortmann                                   
!---------------------------------------------------------------        
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
    SUBROUTINE gf_broyden(filename,l_potmix,l_surface,fmpi,                      &
        &     imix,maxiter,alpha,fm,stars,atoms,sphhar,cell,sym,jspins     &
        &     ,num_layers,sm)
        USE m_gf_metric
        USE m_gf_types
        IMPLICIT NONE
        !     .. Scalar Arguments ..
        LOGICAL,INTENT(IN)       :: l_potmix,l_surface
        INTEGER, INTENT (IN)     :: imix,maxiter
        REAL,    INTENT (IN)     :: alpha
        CHARACTER*(*),INTENT(IN) :: filename
        !     ..
        !   for metric
        TYPE(t_stars),INTENT(IN)   :: stars(:)
        TYPE(t_atoms),INTENT(IN)   :: atoms(:)
        TYPE(t_sphhar),INTENT(IN)  :: sphhar(:)
        TYPE(t_cell),INTENT(IN)    :: cell(:)
        TYPE(t_sym),INTENT(IN)     :: sym(:)
        TYPE(t_gfmpi),INTENT(IN)     :: fmpi
        INTEGER,INTENT(IN)         :: jspins,num_layers
                                                                        
        !     .. Array Arguments ..
        REAL   ,INTENT(IN)      :: fm(:)
        REAL,    INTENT (INOUT) :: sm(:)
        !     ..
        !     .. Local Scalars ..
        INTEGER :: i,it,k,nit,npos,iread,irecl,nmap,mit
        REAL    :: bm,dfivi,fmvm,one,smnorm,vmnorm,zero,alphan
        LOGICAL :: lexist
        !     ..
        !     .. Local Arrays ..
        REAL, ALLOCATABLE :: am(:)
        REAL, ALLOCATABLE :: fm1(:),sm1(:),ui(:),um(:),vi(:),vm(:)
        !     ..
        !     .. External Functions ..
        REAL :: ddot
        EXTERNAL ddot
        !     ..
        !     .. External Subroutines ..
        EXTERNAL daxpy,dscal
        !     ..
        !     .. Data statements ..
        DATA one/1.0E0/,zero/0.0E0/
        !
        nmap = SIZE(fm)
        dfivi = zero
        INQUIRE (FILE = filename ,EXIST = lexist)
        mit = 0
        IF (.NOT.lexist) mit = 1
        irecl = (2*nmap+1)*8
        OPEN (57,FILE =filename,FORM ='unformatted',STATUS ='unknown')
        OPEN (59,FILE =filename//'.'//CHAR(imix+48),ACCESS ='direct',     &
            &     RECL  = irecl,FORM ='unformatted',STATUS='unknown')
                                                                        
                                                                        
                                                                        
        ALLOCATE (fm1(nmap),sm1(nmap),ui(nmap),um(nmap),vi(nmap),vm(nmap))
        ALLOCATE ( am(maxiter+1) )
        !
        IF (mit/=1) THEN
            !
            !     load input charge density (sm1) and  difference of
            !     in and out charge densities (fm1) from previous iteration (m-1)
                                                                        
            REWIND 57
            READ (57) mit,alphan,(fm1(i),i=1,nmap),(sm1(i),i=1,nmap)
            IF ( abs(alpha-alphan) > 0.0001 ) THEN
                WRITE (oUnit,*) 'mixing parameter has been changed; reset'
                WRITE (oUnit,*) 'broyden algorithm or set alpha to',alphan
                CALL juDFT_error("broyden: mixing parameter (alpha) changed",calledby="gf_broyden.F90")
            ENDIF
            !
            !     loop to generate F_m   - F_(m-1)  ... sm1
            !                  and rho_m - rho_(m-1) .. fm1
            !
            DO k = 1,nmap
                sm1(k) = sm(k) - sm1(k)
                fm1(k) = fm(k) - fm1(k)
            END DO
        END IF
        !
        !     save F_m and rho_m for next iteration
        !
        REWIND 57
        nit = mit +1
        IF (nit > maxiter+1) nit = 1
        IF (fmpi%pe0) WRITE (57) nit,alpha,fm,sm
        !
        IF (mit==1) THEN
            !
            !     update for rho for mit = 1 is straight mixing
            !     sm                     = sm + alpha*fm
            !
            CALL daxpy(nmap,alpha,fm,1,sm,1)
        ELSE
            !
            !     |vi> = w|vi>
            !     loop to generate um : um = sm1 + alpha*fm1 - \sum <fm1|w|vi> ui
            !
            DO k = 1,nmap 
                um(k) = alpha*fm1(k) + sm1(k)
            END DO 
            iread = MIN(mit-1,maxiter+1)
            DO it = 2,iread 
                READ (59,REC = it-1) (ui(i),i=1,nmap),(vi(i),i=1,nmap),dfivi
                am(it)  = ddot(nmap,vi,1,fm1,1)
                CALL daxpy(nmap,-am(it),ui,1,um,1)
                WRITE(oUnit,FMT ='(5x,"<vi|w|Fm> for it",i2,5x,f10.6)')it,am(it)
            END DO 
            !
            !     ****************************************
            !     broyden's first method
            !     ****************************************
            !
            IF (imix==3) THEN
                !
                !     convolute drho(m) with the metric: |fm1> = w|sm1>
                !
                fm1 = gf_apply_metric(l_surface,l_potmix,fmpi,stars,atoms,cell,sphhar   &
                    &          ,sym,jspins,num_layers,sm1)
                !
                !     calculate the norm of sm1 : <sm1|w|sm1>
                smnorm = ddot(nmap,sm1,1,fm1,1)
                !
                !     loop to generate vm = alpha*sm1  - \sum <ui|w|sm1> vi
                !
                DO k = 1,nmap
                    vm(k) = alpha*fm1(k)
                END DO
                DO it = 2,iread
                    READ (59,REC=it-1)(ui(i),i=1,nmap),                      &
                        &                           (vi(i),i=1,nmap),dfivi
                    bm = ddot(nmap,ui,1,fm1,1)
                    CALL daxpy(nmap,-bm,vi,1,vm,1)
                    WRITE(oUnit,FMT='(5x,"<ui|w|Fm> for it",i2,5x,f10.6)') it, bm
                END DO
                !
                !     complete evaluation of vm
                !     vmnorm = <um|w|sm1>-<sm1|w|sm1>
                !
                vmnorm = ddot(nmap,fm1,1,um,1) - smnorm
                !     *   if (vmnorm<tol_10) stop
                CALL dscal(nmap,one/vmnorm,vm,1)
            !     write bm(it)
            !
            ELSE IF (imix==5) THEN
                !
                !     ****************************************
                !     broyden's second method
                !     ****************************************
                !
                !     --> multiply fm1 with metric matrix and store in vm:  w |fm1>
                !
                vm = gf_apply_metric(l_surface,l_potmix,fmpi,stars,atoms,cell,sphhar   &
                    &           ,sym,jspins,num_layers,fm1)
                                                                        
                !
                !     calculate the norm of fm1 and normalize vm it: vm = wfm1 /
                !     <fm1|w|fm1>
                !
                vmnorm = one/ddot(nmap,fm1,1,vm,1)
                CALL dscal(nmap,vmnorm,vm,1)
                                                                        
            ELSE IF (imix==7) THEN
                !
                !     ****************************************
                !     generalized anderson method
                !     ****************************************
                !
                !     calculate vm = alpha*wfm1 -\sum <fm1|w|vi> <fi1|w|vi><vi|
                !     convolute fm1 with the metrik and store in vm
                vm = gf_apply_metric(l_surface,l_potmix,fmpi,stars,atoms,cell,sphhar   &
                    &           ,sym,jspins,num_layers,fm1)
                !
                DO it = 2,iread
                    READ (59,REC = it-1) (ui(i),i=1,nmap),                   &
                        &              (vi(i),i = 1,nmap),dfivi
                    CALL daxpy(nmap,-am(it)*dfivi,vi,1,vm,1)
                END DO
                !
                vmnorm = ddot(nmap,fm1,1,vm,1)
                                                                        
                !     *   if (vmnorm<tol_10) stop
                CALL dscal(nmap,one/vmnorm,vm,1)
                !
                !     save dfivi(mit) for next iteration
                !
                dfivi = vmnorm
            END IF
            !
            !     ****************************************
            !
            !     write um,vm and dfivi on file broyd.?
            !
            npos = mit-1
            IF (mit>maxiter+1) THEN
                npos = MOD(mit-2,maxiter)+1
            ENDIF
            IF (fmpi%pe0) WRITE (59,REC = npos) (um(i),i = 1,nmap),(vm(i),i &
                &        = 1,nmap),dfivi
            !
            !     update rho(m+1)
            !
            !     calculate <fm|w|vm>
            !
            fmvm = ddot(nmap,vm,1,fm,1)
            CALL daxpy(nmap,one-fmvm,um,1,sm,1)
        END IF
        mit = mit + 1
        DEALLOCATE ( fm1,sm1,ui,um,vi,vm,am )
        CLOSE (57)
        CLOSE (59)
    END SUBROUTINE gf_broyden
                                                                        
!>
END
