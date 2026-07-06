!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_curvy2dprojector 
#include "juDFT_env.h"
      IMPLICIT NONE
      PRIVATE 
      COMPLEX, ALLOCATABLE, SAVE :: curvyproject(:,:,:)
      COMPLEX, ALLOCATABLE, SAVE :: curvystep(:,:,:,:)
      COMPLEX, ALLOCATABLE, SAVE :: circlestep(:,:,:)
      REAl, SAVE                 :: c
                                                      !remember the last
      INTEGER,SAVE               :: stored_layer = 0 
                                                      !layer calculated 
      PUBLIC :: gf_curvy2dprojector,gf_curvy2dealloc                   &
     &     ,gf_get_curvy2dproject                                       
      CONTAINS 
      !<-- S: gf_curvy2dinit(layer,lapw,cell,l_twothreeproj)            
                                                                        
!**************************************************                     
!     Calculate various step functions needed for                       
!     the use of curvy boundary surfaces.                               
!     This has to be done only once, at start up.                       
!     Frank Freimuth, October 2007                                      
!**************************************************                     
      SUBROUTINE gf_curvy2dinit(layer,lapw,cell,l_twothreeproj) 
                                                                        
                                                                        
      USE m_gf_types 
      USE hdf5 
      USE m_hdf_tools 
      use m_juDFT 
      USE m_gf_circlestep 
      USE m_gf_curvystep 
      USE m_gf_math 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)       :: layer 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      TYPE(t_cell),INTENT(IN)  :: cell 
      LOGICAL, INTENT(IN)      :: l_twothreeproj 
                                                                        
      !COMPLEX, ALLOCATABLE :: circlestep(:,:,:)
      !COMPLEX, ALLOCATABLE :: curvystep (:,:,:,:)
      INTEGER(HID_T)       :: fid,gid,varid,n,ffid,did 
      INTEGER              :: hdferr 
      INTEGER              :: ncutat,newat 
      REAL, ALLOCATABLE    :: cutpos(:,:),newpos(:,:),newrmt(:) 
      REAL, ALLOCATABLE    :: cutrmt(:) 
                                                                        
      INTEGER             :: k1max,k2max,k3max 
      INTEGER             :: dims(3) 
                                                                        
                                                                        
!kmax values (large enough for potential)                               
      k1max = ABS(lapw%g_max(1))+5
      k2max = ABS(lapw%g_max(2))+5
      k3max = MAXVAL(ABS(lapw%k%K3(:lapw%nvd,1)))+5
      CPP_juDFT_timestart("Curvy Projectors construction I")
      IF(.NOT.l_twothreeproj)k3max = 0 
#ifdef CPP_NEVER
      CALL io_hdfopen('gf_setup.hdf',H5F_ACC_RDWR_F,ffid,hdferr) 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      IF (io_groupexists(fid,"curvystep")) THEN 
         CALL io_gopen (fid,'curvystep',gid,hdferr) 
         CALL io_READ_att(gid,"Dims",dims) 
         CALL io_gclose(gid,hdferr) 
         IF (ALL(dims <= (/k1max,k2max,k3max/))) THEN 
            !Step function already calculated                           
            CALL io_gclose (fid, hdferr) 
            CALL io_hdfclose(ffid,hdferr) 
            RETURN 
         ENDIF 
         !present but wrong size!                                       
         WRITE(*,*) dims 
         WRITE(*,*) k1max,k2max,k3max 
         WRITE(*,*) "Curvystep with wrong dimensions -> recalculating" 
         CALL io_gdelete(fid,"curvystep",hdferr) 
      ENDIF 
#endif

                                                                        
      CALL io_hdfopen('gf_setup.hdf',H5F_ACC_RDONLY_F,ffid,hdferr)
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr)
!get information about cutting spheres                                  
      CALL io_gopen (fid,'stepinfo',gid,hdferr) 
      IF (hdferr <0 )  CALL juDFT_error                                   &
     &     ('Information about step functions not found')               
      CALL io_READ_att(gid, "c",c) 
      CALL io_READ_att(gid,"ncutat",ncutat) 
      CALL io_READ_att(gid,"nnewat",newat) 
      ALLOCATE( newpos(3,newat),cutpos(3,ncutat) ) 
      ALLOCATE( newrmt(newat),cutrmt(ncutat) ) 
      CALL io_READ_att(gid,"newpos",newpos) 
      CALL io_READ_att(gid,"newrmt",newrmt) 
      IF (ncutat>0) THEN 
         CALL io_READ_att(gid,"cutpos",cutpos) 
         CALL io_READ_att(gid,"cutrmt",cutrmt) 
      ENDIF 
      CALL io_gclose (gid, hdferr) 
      CALL io_gclose (fid, hdferr) 
                                                                        
!2d step function of the circles that the MT spheres cut out of         
!the planes at +-c/2.
      if (allocated(circlestep)) deallocate(circlestep)
      ALLOCATE( circlestep(-2*k1max:2*k1max,-2*k2max:2*k2max,2) ) 
      circlestep=cmplx(0.0,0.0) 
      IF (newat>0) CALL gf_circlestep(                                  &
     &                  cell%bmat,                                      &
     &                  newpos,newrmt,cell%amat(3,3),c,                 &
     &                  2*k1max,2*k2max,                                &
     &                  circlestep)                                     
      IF(ncutat>0) CALL gf_circlestep(                                  &
     &                  cell%bmat,                                      &
     &                  cutpos,cutrmt,cell%amat(3,3),c,                 &
     &                  2*k1max,2*k2max,                                &
     &                  circlestep)                                     
                                                                        
!2d step function of those parts of the MT sphere surfaces which        
!lie within the curvy boundary.
      if (allocated(curvystep)) deallocate(curvystep)
      ALLOCATE( curvystep(-2*k1max:2*k1max,-2*k2max:2*k2max,            &
     &                                        -k3max:k3max,2) )         
      curvystep=cmplx(0.0,0.0) 
      IF (newat>0) CALL gf_curvystep(                                   &
     &                  cell%bmat,                                      &
     &                  newpos,newrmt,cell%amat(3,3),c,                 &
     &                  2*k1max,2*k2max,k3max,                          &
     &                  curvystep)                                      
      IF(ncutat>0) CALL gf_curvystep(                                   &
     &                  cell%bmat,                                      &
     &                  cutpos,cutrmt,cell%amat(3,3),c,                 &
     &                  2*k1max,2*k2max,k3max,                          &
     &                  curvystep)                                      

#ifdef CPP_NEVER
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      CALL io_gcreate(fid,"curvystep",gid,hdferr) 
      CALL io_write_att(gid,"Dims",(/k1max,k2max,k3max/)) 
      CALL io_write_att(gid, "c",c) 
      CALL io_createvar(gid,"circlestep", H5T_NATIVE_DOUBLE,(/           &
     &     2,SIZE(circlestep,1),SIZE(circlestep,2),SIZE(circlestep,3)/)&
     &     ,did)                                                        
      CALL io_WRITE(did,(/-1,1,1,1/),(/1,SIZE(circlestep,1)             &
     &     ,SIZE(circlestep,2),SIZE(circlestep,3)/),circlestep)         
      CALL io_dclose(did,hdferr) 
                                                                        
      CALL io_createvar(gid,"curvystep", H5T_NATIVE_DOUBLE,(/            &
     &     2,SIZE(curvystep,1),SIZE(curvystep,2),SIZE(curvystep,3)     &
     &     ,SIZE(curvystep,4)/),did)                                    
      CALL io_WRITE(did,(/-1,1,1,1,1/),(/1,SIZE(curvystep,1)            &
     &     ,SIZE(curvystep,2),SIZE(curvystep,3),SIZE(curvystep,4)/)     &
     &     ,curvystep)                                                  
      CALL io_dclose(did,hdferr) 
                                                                        
      CALL io_gclose(gid,hdferr) 

                                                                        
      DEALLOCATE(curvystep) 
      DEALLOCATE(circlestep) 
#endif
      !CALL io_gclose(fid,hdferr)
      CALL io_hdfclose(ffid,hdferr)
      CPP_juDFT_timestop("Curvy Projectors construction I")

      END SUBROUTINE gf_curvy2dinit 
                                                                        
      !>                                                                
      !<-- S: gf_curvy2dprojector(layer,cell,lapw,l_twothreeproj)       
                                                                        
!**************************************************                     
!     Calculate the projector for curvy surfaces.                       
!     Frank Freimuth, October 2007                                      
!                                                                       
!  Note on noco:                                                        
!     in here nv(1) and nv2(1) are used, i.e. code is not suitable for l
!                                                                       
!**************************************************                     
      SUBROUTINE gf_curvy2dprojector(layer,                             &
     &     cell,lapw,l_twothreeproj)                                    
      USE m_gf_types 
      use m_juDFT 
      USE m_gf_math 
      USE hdf5 
      USE m_hdf_tools 
                                                                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)       :: layer 
      TYPE(t_cell),INTENT(IN)  :: cell 
      TYPE(t_lapw),INTENT(IN)  :: lapw 
      LOGICAL,INTENT(IN)       :: l_twothreeproj 
                                                                        
      !COMPLEX, ALLOCATABLE :: circlestep(:,:,:)
      !COMPLEX, ALLOCATABLE :: curvystep (:,:,:,:)
      COMPLEX              :: curvy2dproj_tmp(lapw%nv2(1),lapw%nv(1),2) 
      COMPLEX              :: exp1,norm,normcircle 
      INTEGER              :: n1,n2 
      INTEGER              :: k1max,k2max,k1,k2,k3max,k3,n3 
      INTEGER              :: k_old(3) 
      COMPLEX,ALLOCATABLE  :: curvyoverlaps(:,:,:) 
      !REAL                 :: c
      !INTEGER(HID_T)       :: fid,gid,ffid,did
      INTEGER              :: hdferr
      COMPLEX, ALLOCATABLE :: basisoverlaps(:,:,:) 
      COMPLEX, ALLOCATABLE :: basisoverlaps_inv(:,:,:) 


      IF (layer/=stored_layer) THEN 
         !try to read data from file
         !rewind(701+layer,err=100)
         !read(701+layer,err=100,end=100) n1,n2
         !if (n1==lapw%nv2(1).and.n2==lapw%nv(1)) then
         !   if(.not.allocated(curvyproject)) allocate(curvyproject(n1,n2,2))
        ! 	read(701+layer,err=100,end=100) curvyproject
        ! 	stored_layer=layer
         !	return
         !endif
         !reading failed
         CALL gf_curvy2dinit(layer,lapw,cell,l_twothreeproj)
         stored_layer = layer 
      ELSE 
         IF (ALLOCATED(curvyproject)) THEN 
            IF (SIZE(curvyproject,1) == lapw%nv2(1).AND.                 &
     &           SIZE(curvyproject,2) == lapw%nv(1)) THEN
               RETURN 
            ENDIF
         ENDIF
      ENDIF 

      CPP_juDFT_timestart("Curvy Projectors construction II")
                                                                        
      k1max=maxval(abs(lapw%kp%k1p(:lapw%nv2(1),1)))+2 
      k2max=maxval(abs(lapw%kp%k2p(:lapw%nv2(1),1)))+2 
      k3max = MAXVAL(ABS(lapw%k%k3(:lapw%nv(1),1))) 
                                                                        
                                                                        
                                                                        
!*****************************************************                  
!     Calculate overlap between curvy basis functions.                  
!*****************************************************                  
      ALLOCATE( curvyoverlaps(-2*k1max:2*k1max,-2*k2max:2*k2max,2) ) 
      curvyoverlaps(:,:,:)=cmplx(0.0,0.0) 
      norm=cmplx(cell%amat(3,3)/cell%omtil,0.0) 
      normcircle=cmplx(cell%amat(3,3)/cell%omtil,0.0) 
#ifdef CPP_NEVER
      CALL io_hdfopen('gf_setup.hdf',H5F_ACC_RDONLY_F,ffid,hdferr) 
      CALL io_gopen(ffid,io_layername(layer),fid,hdferr) 
      CALL io_gopen(fid,"curvystep",gid,hdferr) 
      CALL io_READ_att(gid, "c",c) 
      CALL io_READ_att(gid,"Dims",k_old) 
      IF (ANY((/k1max,k2max,k3max/)>k_old)) THEN 
         WRITE(*,*) "k_max:",k1max,k2max,k3max 
         WRITE(*,*) "File: ",k_old 
         CALL juDFT_error                                                 &
     &     ("Wrong dimension of curvystep-data")                        
      ENDIF
#endif
      k_old(1)=size(curvystep,1)/4
      k_old(2)=size(curvystep,2)/4
      k_old(3)=floor(size(curvystep,3)/2.)
      IF (ANY((/k1max,k2max,k3max/)>k_old)) THEN
         WRITE(*,*) "k_max:",k1max,k2max,k3max
         WRITE(*,*) "before: ",k_old
         CALL juDFT_error                                                 &
     &     ("Wrong dimension of curvystep-data")
      ENDIF
#ifdef CPP_NEVER
      !ALLOCATE(circlestep(-2*k_old(1):2*k_old(1),-2*k_old(2):2*k_old(2) &
     !&     ,2))
      !ALLOCATE(curvystep(-2*k_old(1):2*k_old(1),-2*k_old(2):2*k_old(2), &
     !&                   -k_old(3):k_old(3),2))
                                                                        
      CALL io_dopen(gid,"circlestep",did,hdferr) 
      CALL io_READ(did,(/-1,1,1,1/),(/1,SIZE(circlestep,1)              &
     &     ,SIZE(circlestep,2),SIZE(circlestep,3)/),circlestep)         
      CALL io_dclose(did,hdferr) 
                                                                        
      CALL io_dopen(gid,"curvystep",did,hdferr) 
      CALL io_read(did,(/-1,1,1,1,1/),(/1,SIZE(curvystep,1)             &
     &     ,SIZE(curvystep,2),SIZE(curvystep,3),SIZE(curvystep,4)/)     &
     &     ,curvystep)                                                  
      CALL io_dclose(did,hdferr) 
                                                                        
      CALL io_gclose(gid,hdferr) 
      CALL io_gclose(fid,hdferr) 
      CALL io_hdfclose(ffid,hdferr) 
#endif
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!generate the 2d-kronecker delta                                        
      curvyoverlaps(0,0,1)=cmplx(1.0,0.0) 
      curvyoverlaps(0,0,2)=cmplx(1.0,0.0) 
                                                                        
!subtract the circle                                                    
      curvyoverlaps(-2*k1max:2*k1max,-2*k2max:2*k2max,1:2)=             &
     &   curvyoverlaps(-2*k1max:2*k1max,-2*k2max:2*k2max,1:2)-          &
     &   circlestep(-2*k1max:2*k1max,-2*k2max:2*k2max,1:2)*norm         
                                                                        
!add the integrals over the sphere surfaces                             
      curvyoverlaps(-2*k1max:2*k1max,-2*k2max:2*k2max,1:2)=             &
     &   curvyoverlaps(-2*k1max:2*k1max,-2*k2max:2*k2max,1:2)+          &
     &   curvystep(-2*k1max:2*k1max,-2*k2max:2*k2max,0,1:2)*norm        
                                                                        
      ALLOCATE( basisoverlaps(lapw%nv2(1),lapw%nv2(1),2) ) 
      ALLOCATE( basisoverlaps_inv(lapw%nv2(1),lapw%nv2(1),2) ) 
                                                                        
                                                                        
      DO n2=1,lapw%nv2(1) 
        DO n1=1,lapw%nv2(1) 
           k1=lapw%kp%k1p(n2,1)-lapw%kp%k1p(n1,1) 
           k2=lapw%kp%k2p(n2,1)-lapw%kp%k2p(n1,1) 
           basisoverlaps(n1,n2,1)=curvyoverlaps(k1,k2,1) 
           basisoverlaps(n1,n2,2)=curvyoverlaps(k1,k2,2) 
        ENDDO 
      ENDDO 
      DEALLOCATE(curvyoverlaps) 
                                                                        
      basisoverlaps_inv(:,:,1)=mat_inverse(basisoverlaps(:,:,1)) 
      basisoverlaps_inv(:,:,2)=mat_inverse(basisoverlaps(:,:,2)) 
                                                                        
                                                                        
                                                                        
                                                                        
                                                                        
!***************************************************************        
!     Calculate projections of 3d basis onto 2d basis functions.        
!***************************************************************        
      IF(l_twothreeproj)THEN 
      if (allocated(curvyproject)) DEALLOCATE(curvyproject)
      ALLOCATE(curvyproject(lapw%nv2(1),lapw%nv(1),2)) 
      norm=cmplx(1.0/sqrt(cell%amat(3,3)),0.0) 
      normcircle=cmplx(sqrt(cell%amat(3,3))/cell%omtil,0.0) 
      curvyproject(:,:,:)=cmplx(0.0,0.0) 
      exp1=exp(cmplx(0.0,cell%bmat(3,3)*c/2.0)) 
                                                                        
!generate the 2d-kronecker delta                                        
      DO n1=1,lapw%nv(1) 
         curvyproject(lapw%k%kp(n1,1),n1,1)=                            &
     &          exp1**(-1.0*lapw%k%k3(n1,1))*norm                       
      ENDDO 
      DO n1=1,lapw%nv(1) 
         curvyproject(lapw%k%kp(n1,1),n1,2)=                            &
     &          exp1**( 1.0*lapw%k%k3(n1,1))*norm                       
      ENDDO 
                                                                        
!subtract the circle                                                    
      DO n2=1,lapw%nv(1) 
         k3=lapw%k%k3(n2,1) 
         DO n1=1,lapw%nv2(1) 
            k1=lapw%k%k1(n2,1)-lapw%kp%k1p(n1,1) 
            k2=lapw%k%k2(n2,1)-lapw%kp%k2p(n1,1) 
                                                                        
            curvyproject(n1,n2,1)=                                      &
     &          curvyproject(n1,n2,1)-                                  &
     &          exp1**(-k3)*normcircle*                                 &
     &          circlestep(k1,k2,1)                                     
               !n1                                                      
         ENDDO 
            !n2                                                         
      ENDDO 
                                                                        
      DO n2=1,lapw%nv(1) 
         k3=lapw%k%k3(n2,1) 
         DO n1=1,lapw%nv2(1) 
            k1=lapw%k%k1(n2,1)-lapw%kp%k1p(n1,1) 
            k2=lapw%k%k2(n2,1)-lapw%kp%k2p(n1,1) 
            curvyproject(n1,n2,2)=                                      &
     &          curvyproject(n1,n2,2)-                                  &
     &          exp1**(1.0*k3)*normcircle*                              &
     &          circlestep(k1,k2,2)                                     
               !n1                                                      
         ENDDO 
            !n2                                                         
      ENDDO 
                                                                        
!add the integrals over the sphere surfaces                             
                                                                        
      DO n2=1,lapw%nv(1) 
         k3=lapw%k%k3(n2,1) 
         DO n1=1,lapw%nv2(1) 
            k1=lapw%k%k1(n2,1)-lapw%kp%k1p(n1,1) 
            k2=lapw%k%k2(n2,1)-lapw%kp%k2p(n1,1) 
            curvyproject(n1,n2,1)=                                      &
     &          curvyproject(n1,n2,1)+                                  &
     &          exp1**(-1.0*k3)*normcircle*                             &
     &          curvystep(k1,k2,k3,1)                                   
               !n1                                                      
         ENDDO 
            !n2                                                         
      ENDDO 
                                                                        
      DO n2=1,lapw%nv(1) 
         k3=lapw%k%k3(n2,1) 
         DO n1=1,lapw%nv2(1) 
            k1=lapw%k%k1(n2,1)-lapw%kp%k1p(n1,1) 
            k2=lapw%k%k2(n2,1)-lapw%kp%k2p(n1,1) 
            curvyproject(n1,n2,2)=                                      &
     &          curvyproject(n1,n2,2)+                                  &
     &          exp1**(1.0*k3)*normcircle*                              &
     &          curvystep(k1,k2,k3,2)                                   
               !n1                                                      
         ENDDO 
            !n2                                                         
      ENDDO 
                                                                        
      curvy2dproj_tmp=curvyproject 
      curvyproject=cmplx(0.0,0.0) 
      DO n3=1,lapw%nv2(1) 
       DO n2=1,lapw%nv(1) 
        DO n1=1,lapw%nv2(1) 
           curvyproject(n1,n2,1)=curvyproject(n1,n2,1)+                 &
     &                       basisoverlaps_inv(n1,n3,1)*                &
     &                         curvy2dproj_tmp(n3,n2,1)                 
        ENDDO 
       ENDDO 
      ENDDO 
                                                                        
      DO n3=1,lapw%nv2(1) 
       DO n2=1,lapw%nv(1) 
        DO n1=1,lapw%nv2(1) 
           curvyproject(n1,n2,2)=curvyproject(n1,n2,2)+                 &
     &                       basisoverlaps_inv(n1,n3,2)*                &
     &                         curvy2dproj_tmp(n3,n2,2)                 
        ENDDO 
       ENDDO 
      ENDDO 
            !l_twothreeproj                                             
      ENDIF
      CPP_juDFT_timestop("Curvy Projectors construction II")

      !Save to SCRATCH file
      !open(701+layer,STATUS='SCRATCH',form="unformatted")
      !write(701+layer) lapw%nv2(1),lapw%nv(1)
      !write(701+layer) curvyproject
      END SUBROUTINE gf_curvy2dprojector 
                                                                        
      !>                                                                

      subroutine gf_curvy2dealloc()
      if (allocated(curvyproject)) deallocate(curvyproject)
      end subroutine

      !<-- S:gf_curvy2dproject(lapw,matrix,proj,transposed)             
      SUBROUTINE gf_get_curvy2dproject(lapw,proj,transposed) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-06-30) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_lapw),INTENT(IN)      :: lapw 
      COMPLEX,INTENT(OUT)          :: proj(:,:) 
      LOGICAL,INTENT(IN),OPTIONAL  :: transposed 
      !>                                                                
      !<-- Locals                                                       
      LOGICAL             :: tp 
      INTEGER             :: nv,nmat
      !>                                                                
      tp = .FALSE. 
      IF (PRESENT(transposed)) tp = transposed 
      proj(:,:) = 0.0


      IF (.NOT.tp) THEN 
         !<-- normal projector
         if (size(proj,2)>lapw%nmat_sphere) then
             nv=lapw%nv(1)
             nmat=lapw%nmat
         else
             nv=lapw%nv_sphere(1)
             nmat=lapw%nmat_sphere
         endif
         proj(1:lapw%nv2(1),1:nv) = curvyproject(1:lapw%nv2(1),1:nv,1)
         proj(lapw%nv2_tot+1:lapw%nv2_tot+lapw%nv2(1),1:nv) =   &
     &        curvyproject(1:lapw%nv2(1),1:nv,2)
                                                                        
         IF (lapw%nv2(1) /= lapw%nv2_tot) THEN 
            !noco                                                       
            proj(lapw%nv2_tot/2+1:lapw%nv2_tot/2+lapw%nv2(1),nmat/2+1:nmat/2+nv) =               &
     &           curvyproject(1:lapw%nv2(1),1:nv,1)
            proj(3*lapw%nv2_tot/2+1:3*lapw%nv2_tot/2+lapw%nv2(1)        &
     &           ,nmat/2+1:nmat/2+nv) =               &
     &           curvyproject(1:lapw%nv2(1),1:nv,2)
         ENDIF 
         !>                                                             
      ELSE 
         !<-- transposed version
         if (size(proj,1)>lapw%nmat_sphere) then
             nv=lapw%nv(1)
             nmat=lapw%nmat
         else
             nv=lapw%nv_sphere(1)
             nmat=lapw%nmat_sphere

         endif
         proj(1:nv,1:lapw%nv2(1)) = TRANSPOSE(CONJG(            &
     &        curvyproject(1:lapw%nv2(1),1:nv,1)))
         proj(1:nv,lapw%nv2_tot+1:lapw%nv2_tot+lapw%nv2(1)) =   &
     &        TRANSPOSE(CONJG(curvyproject(1:lapw%nv2(1),1:nv,2)&
     &        ))                                                        
          IF (lapw%nv2(1) /= lapw%nv2_tot) THEN 
             !noco                                                      
             proj(nmat/2+1:nmat/2+nv,lapw%nv2_tot/2   &
     &            +1:lapw%nv2_tot/2+lapw%nv2(1)) =                      &
     &            TRANSPOSE(CONJG(curvyproject(1:lapw%nv2(1),1:nv,1)))
             proj(nmat/2+1:nmat/2+nv,3*lapw%nv2_tot/2 &
     &            +1:3*lapw%nv2_tot/2+lapw%nv2(1)) =                    &
     &            TRANSPOSE(CONJG(curvyproject(1:lapw%nv2(1),1:nv,2)))
          ENDIF 
          !>                                                            
      ENDIF 
      END SUBROUTINE 
      !>                                                                
      END                                           
