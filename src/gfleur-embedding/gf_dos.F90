!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_dos 
      use m_juDFT 
      USE hdf5 
      USE m_hdf_tools 
      IMPLICIT NONE

      PRIVATE 
                                                        !norm of udot   
      REAL,ALLOCATABLE,SAVE    :: ddn(:,:,:,:)
                                                        !wigner matrices
      COMPLEX,ALLOCATABLE,SAVE :: d_wgn(:,:,:,:) 
                                                        !no of ops for  
      INTEGER,ALLOCATABLE,SAVE :: invarind(:,:) 
                                                        !each atom      
                                                        !sym ops for eac
      INTEGER,ALLOCATABLE,SAVE :: invarop(:,:,:) 
                                                                        
      INTEGER(HID_T),SAVE :: fid,ldos_id = 0,surfdos_id,intdos_id 
                                                                        
      INTEGER,SAVE             :: nx,ny,nz 
                                                                        
      PUBLIC gf_dos_init,gf_dos_mt,gf_writedos,gf_dos_plane,gf_dos_int 
      CONTAINS 
      !<-- S:gf_dos_plane(g,gfinp,lapw,cell,en,jspin,wk,nk)             
                                                                        
      SUBROUTINE gf_dos_plane(g,gfinp,lapw,cell,en,jspin,wk,nk) 
!*****************************************************************      
!  Calculate the DOS on some plane in the vacuum                        
!*****************************************************************      
#include "cpp_double.h"
      USE m_gf_types 
      USE m_gf_fft_singleton 
      IMPLICIT NONE 
!     Arguments                                                         
                                           ! Green function             
      COMPLEX, INTENT(IN)        :: g(:,:) 
      TYPE(t_gfinp),INTENT(IN)   :: gfinp 
      TYPE(t_cell),INTENT(IN)    :: cell 
      TYPE(t_lapw),INTENT(IN)    :: lapw 
                                                !energy,spin            
      INTEGER,INTENT(IN)         :: en,jspin,nk 
                                       !weight of kpt                   
      REAL   ,INTENT(IN)         :: wk 
                                                                        
                                                                        
                                   !position of plane                   
      REAL                :: z_dos 
      COMPLEX,ALLOCATABLE :: work(:,:,:,:),surf_dos(:,:) 
      INTEGER             :: n1,n2,k1_1,k1_2,k2_1,k2_2,in1,in2,z 
      COMPLEX             :: exp1,exp2,f1,f2 
                                                                        
                                                                        
      ALLOCATE(work(nx,ny,nx,ny),surf_dos(nx,ny)) 
                                                                        
      DO z=1,nz 
         z_dos=gfinp%z_min+(gfinp%z_max-gfinp%z_min)/(nz-1)*(z-1) 
         !<-- Put on FFT-grid                                           
         exp1 = EXP(CMPLX(0.0,cell%bmat(3,3)*z_dos)) 
         exp2 = EXP(CMPLX(0.0,-1.0*cell%bmat(3,3)*z_dos)) 
         work = 0.0 
         DO n1 = 1,lapw%nv_sphere(jspin)
            k1_1 = lapw%k%k1(n1,jspin) 
            IF (k1_1<0) THEN 
               k1_1 = nx+k1_1-1 
            ENDIF 
            k1_2 = lapw%k%k2(n1,jspin) 
            IF (k1_2<0) THEN 
               k1_2 = ny+k1_2-1 
            ENDIF 
            f1  = 1/cell%amat(3,3)*exp1**lapw%k%k3(n1,jspin) 
            DO n2 = 1,lapw%nv_sphere(jspin)
               k2_1 = lapw%k%k1(n2,jspin) 
               IF (k2_1<0) THEN 
                  k2_1 = nx+k2_1-1 
               ENDIF 
               k2_2 = lapw%k%k2(n2,jspin) 
               IF (k2_2<0) THEN 
                  k2_2 = ny+k2_2-1 
               ENDIF 
               f2  = exp2**lapw%k%k3(n2,jspin) 
               work(k1_1+1,k1_2+1,k2_1+1,k2_2+1) = work(k1_1+1,k1_2+1   &
     &              ,k2_1+1,k2_2+1)+g(n1,n2)*f1*f2                      
            ENDDO 
         ENDDO 
         !>                                                             
         work=fft(work,inv=.TRUE.) 
         !<-- Collect diagonal elements                                 
         surf_dos = 0.0 
         DO n1 = 1,nx 
            DO n2 = 1,ny 
               surf_dos(n1,n2) = surf_dos(n1,n2)+wk                     &
     &              *work(n1,n2,n1,n2)                                  
            ENDDO 
         ENDDO 
         CALL io_WRITE(surfdos_id,(/1,1,1,z,en,nk,jspin/),(/1,nx,ny,1,1 &
     &        ,1,1,1/),REAL(surf_dos))                                  
         CALL io_WRITE(surfdos_id,(/2,1,1,z,en,nk,jspin/),(/1,nx,ny,1,1 &
     &        ,1,1,1/),AIMAG(surf_dos))                                 
         !>                                                             
      ENDDO 
      DEALLOCATE(work,surf_dos) 
                                                                        
      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_dos_mt                                                  
                                                                        
      recursive SUBROUTINE gf_dos_mt(layer,atoms,gfinp,                           &
     &     gmat,sym,en,jspin_in,wk,nk,l_noco,lapw)
!*****************************************************************      
!     This subroutine constructs the density matrix for l<4 for all     
!     atom types.                                                       
!*****************************************************************      
#include "cpp_double.h"
      USE m_gf_types 
      use m_gf_ab_coef
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,INTENT(IN)         :: layer 
                                              ! The Green function      
      COMPLEX, INTENT(IN)        :: gmat(:,:) 
      TYPE(T_atoms),INTENT(IN)   :: atoms 
      TYPE(t_gfinp),INTENT(IN)   :: gfinp 
      TYPE(t_sym),INTENT(IN)     :: sym
      type(t_lapw),intent(in)    :: lapw
                                                             ! lapw info
                                                !energy,spin            
      INTEGER,INTENT(IN)         :: en,jspin_in,nk
                                       !weight of kpt                   
      REAL   ,INTENT(IN)         :: wk
      logical,intent(in)         :: l_noco
!     Locals
      INTEGER :: jspin,jspin1,jspin2
                                           !lm-loops                    
      INTEGER :: l,lm,lmp,m,mp 
                                           !atoms-loops                 
      INTEGER :: na,itype,nt 
                                           !result of sum over g,g'     
      COMPLEX :: uu,dd 
                                           !loop over rotations         
      INTEGER :: it 
                                           !index of rotation           
      INTEGER :: is,isi 
                                           !scaling for rotations       
      REAL    :: fac 
                                           !tmp for density matrix      
      COMPLEX :: n_tmp(-3:3,-3:3) 
                                           !tmp for wigner rot matrix   
      COMPLEX :: d_tmp(-3:3,-3:3) 
                                                     !tmps for matmul   
      COMPLEX :: n1_tmp(-3:3,-3:3),nr_tmp(-3:3,-3:3) 
                                                !tmp arrays for lapack-g
      COMPLEX,ALLOCATABLE ::vec1(:),vec2(:) 
      COMPLEX,ALLOCATABLE:: ldos_mat(:,:,:,:) 
      COMPLEX CPP_BLAS_cdotu 
      EXTERNAL CPP_BLAS_cdotu 

      ! In case of NOCO, call this subroutine recursivly
      if (l_noco) THEN
      	call gf_dos_mt(layer,atoms,gfinp,                           &
     &     gmat(:size(gmat,1)/2,:size(gmat,2)/2)                    &
     &     ,sym,en,1,wk,nk,.false.,lapw)
        call gf_dos_mt(layer,atoms,gfinp,                           &
     &     gmat(size(gmat,1)/2+1:,size(gmat,2)/2+1:)                    &
     &     ,sym,en,2,wk,nk,.false.,lapw)
        call gf_dos_mt(layer,atoms,gfinp,                           &
     &     gmat(:size(gmat,1)/2,size(gmat,2)/2+1:)                    &
     &     ,sym,en,3,wk,nk,.false.,lapw)
        return
      endif
      !<--Before first call, the wigner-matrixes and radial normalisatio
      if (jspin_in<3) THEN
        jspin=jspin_in
        jspin1=jspin_in
        jspin2=jspin_in
      else
        jspin=1
        jspin1=1
        jspin2=2
      endif
      !has to be initialized                                            
      IF (.NOT.ALLOCATED(d_wgn)) THEN 
         CALL juDFT_error("gf_mtdos called without Initialisation") 
      ENDIF 
                                                                        
      !>                                                                
                                                                        
      ALLOCATE(ldos_mat(-3:3,-3:3,0:3,atoms%ntype)) 
      ldos_mat=0.0 
      !<-- Loop over all atom types                                     
                                                                        
      ALLOCATE(vec1(lapw%nv_sphere(jspin)),vec2(lapw%nv_sphere(jspin)))
      nt = 0 
      DO itype=1,atoms%ntype 
         DO na = 1,atoms%neq(itype) 
            nt = nt + 1 
            DO l=0,3 
               n_tmp=0.0
               DO m=-l,l 
                  lm=(l*(l+1))+m 
                  DO mp=-l,l 
                     lmp=(l*(l+1))+mp 
                                                                        
                     !Use BLAS to construct the uu,dd

                     CALL CPP_BLAS_cgemv('n',lapw%nv_sphere(jspin)            &
     &                    ,lapw%nv_sphere(jspin),CMPLX(1.0,0.0),gmat             &
     &                    ,SIZE(gmat,1),CONJG(gf_ab_coef_vector(lmp,nt,1,jspin2))      &
     &                    ,1,CMPLX(0.0,0.0),vec1,1)                     
                     CALL CPP_BLAS_cgemv('n',lapw%nv_sphere(jspin)            &
     &                    ,lapw%nv_sphere(jspin),CMPLX(1.0,0.0),gmat             &
     &                    ,SIZE(gmat,1),CONJG(gf_ab_coef_vector(lmp,nt,2,jspin2))      &
     &                    ,1,CMPLX(0.0,0.0),vec2,1)                     
                     uu=(CPP_BLAS_cdotu(lapw%nv_sphere(jspin),gf_ab_coef_vector(lm,nt,1,jspin1),1,vec1,1))
                     dd=(CPP_BLAS_cdotu(lapw%nv_sphere(jspin),gf_ab_coef_vector(lm,nt,2,jspin1),1,vec2,1))
                     if (jspin_in==3) uu=uu*ddn(l,itype,layer,4)
                     n_tmp(m,mp) = (ddn(l,itype,layer,jspin_in)*dd+uu)           &
     &                    /atoms%neq(itype)                             
                  ENDDO 
               ENDDO 
               !<--some rotations have to be applied                    
                                                                        
               !compare n_mat of fleur                                  
               fac = 1.0/ ( invarind(nt,layer) * atoms%neq(itype) ) 
               DO it = 1, invarind(nt,layer) 
                  is = invarop(nt,it,layer) 
                  isi = sym%invtab(is) 
                  d_tmp(:,:) = CMPLX(0.0,0.0) 
                  IF (l==0) THEN 
                                    !no wigner matrix for l=0           
                     d_tmp(0,0)=1.0 
                  ELSE 
                     d_tmp(-l:l,-l:l)=d_wgn(-l:l,-l:l,l,isi) 
                  ENDIF 
                  nr_tmp = MATMUL( TRANSPOSE( CONJG(d_tmp) ) , n_tmp) 
                  n1_tmp =  MATMUL( nr_tmp, d_tmp ) 
                  ldos_mat(-l:l,-l:l,l,itype) =                         &
     &                 ldos_mat(-l:l,-l:l,l,itype)+                     &
     &                 wk*fac* n1_tmp(-l:l,-l:l)

               ENDDO 
                                                                        
               !>                                                       
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
      !>                                                                
                                                                        
      CALL io_WRITE(ldos_id,(/1,1,1,1,1,en,nk,jspin_in,layer/),(/1,7,7,4   &
     &     ,atoms%ntype,1,1,1,1/),REAL(ldos_mat))                       
      CALL io_WRITE(ldos_id,(/2,1,1,1,1,en,nk,jspin_in,layer/),(/1,7,7,4   &
     &     ,atoms%ntype,1,1,1,1/),AIMAG(ldos_mat))                      
                                                                        
      DEALLOCATE(vec1,vec2,ldos_mat) 
                                                                        
      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_dos_init                                                
                                                                        
      SUBROUTINE gf_dos_init(layers,gfinp,atoms,cell,sym,lapw,jspins    &
     &     ,nkpts,potential,enpara,l_noco)
!******************************************                             
!     initial PRIVATE arrays of the MODULE                              
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_types 
      USE m_gf_energies 
      USE m_intgr
      USE m_fleur_interface,ONLY:fleur_radfun,fleur_d_wigner 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_atoms),INTENT(IN)  :: atoms(:) 
      TYPE(t_gfinp),INTENT(IN)  :: gfinp 
      TYPE(t_sym),INTENT(IN)    :: sym 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_lapw),INTENT(IN)   :: lapw 
      INTEGER,INTENT(IN)        :: jspins,nkpts
      TYPE(t_potden),INTENT(IN) :: potential(:) 
      TYPE(t_enpara),INTENT(IN)    :: enpara(:)
      LOGICAL,INTENT(IN)           :: l_noco

      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER :: itype,l,layer,jspd,jspin
                                       !dummy args                      
      REAL    :: duds,dus,uds,us,wronk 
                                       !for radfun                      
      INTEGER :: noded,nodeu 
      REAL,ALLOCATABLE :: f(:,:,:),g(:,:,:)
      INTEGER          :: na,n,nt,jop,j3,j2,j1 
      REAL             :: gaminv(3),sr(3),s3 
      LOGICAL          :: firstcall 
      INTEGER(hid_t)   :: access_prp,access_mode 
      INTEGER          :: hdferr 
#ifdef CPP_HDFMPI
      INCLUDE 'mpif.h' 
      CALL h5pcreate_f(H5P_FILE_ACCESS_F, access_prp, hdferr) 
      CALL io_check("Could not create access properties",hdferr) 
      CALL h5pset_fapl_mpio_f(access_prp, MPI_COMM_WORLD,               &
     &     MPI_INFO_NULL,hdferr)                                        
      CALL io_check("Could not modify access properties",hdferr) 
#else                                                                   
      access_prp = H5P_DEFAULT_f 
#endif                                                                  
                                                                        
      !>                                                                
      firstcall=.FALSE. 
      jspd=jspins
      if (l_noco) jspd=3
         !<--Allocate module vars                                       
                                                                        
         IF(.NOT.ALLOCATED(ddn)) THEN 
            firstcall=.TRUE.
            IF (l_noco) THEN
                ALLOCATE(ddn(0:3,MAXVAL(atoms%ntype),layers%num_layers,4))
            ELSE
                ALLOCATE(ddn(0:3,MAXVAL(atoms%ntype),layers%num_layers,jspins))
            ENDIF
            ALLOCATE(d_wgn(-3:3,-3:3,3,sym%nop)) 
            ALLOCATE(invarind(MAXVAL(atoms%nat),layers%num_layers)) 
            ALLOCATE(invarop(MAXVAL(atoms%nat),sym%nop,layers%num_layers&
     &           ))                                                     
            CALL h5fcreate_f("gf_dos.hdf",H5F_ACC_TRUNC_F,fid, hdferr   &
     &           ,H5P_DEFAULT_f,access_prp)                             
            !create dataset                                             
            CALL io_createvar(fid,"bands",H5T_NATIVE_DOUBLE,            &
     &           (/2,7,7,4,MAXVAL(atoms%ntype),gf_noen(),nkpts,jspd   &
     &           ,layers%num_layers/),ldos_id)                          
            CALL io_createvar(fid,"intbands",H5T_NATIVE_DOUBLE,         &
     &           (/2,gf_noen(),nkpts,jspd                             &
     &           ,layers%num_layers/),intdos_id)                        
         IF (gfinp%l_surface) THEN 
            !prepare for planar d                                       
                                                                        
            nx = lapw%g_MAX(1)*2+1 
            ny = lapw%g_MAX(2)*2+1 
            nz = gfinp%nz_dos 
            !create_dataset                                             
            CALL io_createvar(fid,"surface",H5T_NATIVE_DOUBLE,          &
     &        (/2,nx,ny,nz,gf_noen(),nkpts,jspd/),surfdos_id)
         ENDIF 
         ENDIF 
                                                                        
         !>                                                             
                                                                        
      !<-- Set up radial functions to get ddn!                          

        DO layer = 1,layers%num_layers
             ALLOCATE(f(SIZE(atoms(layer)%rmsh,1),2,2),g(SIZE(atoms(layer     &
     &            )%rmsh,1),2,2))
            DO itype = 1,atoms(layer)%ntype
                DO l  = 0,3
                    DO jspin=1,jspins
                          CALL fleur_radfun(                                       &
     &                  l,enpara(layer)%el(l+1,itype,jspin),potential(layer &
     &                    )%mt(:,0,itype,jspin),itype,atoms(layer),f(:,:,jspin),g(: &
     &                   ,:,jspin),us,dus,uds,duds,ddn(l,itype,layer,jspin),nodeu,noded  &
     &                  ,wronk)
                    ENDDO
                    IF (l_noco) THEN
                         CALL intgr0(g(:,1,1)*g(:,1,2)+g(:,2,1)*g(:,2,2),atoms(layer)%rmsh(1,itype),  &
                         atoms(layer)%dx(itype),atoms(layer)%jri(itype),ddn(l,itype,layer,3))
                         CALL intgr0(f(:,1,1)*f(:,1,2)+f(:,2,1)*f(:,2,2),atoms(layer)%rmsh(1,itype),  &
                         atoms(layer)%dx(itype),atoms(layer)%jri(itype),ddn(l,itype,layer,4))
                    ENDIF
               ENDDO
             ENDDO
         DEALLOCATE (f,g) 
      ENDDO 
                                                                        
                                                                        
      !>                                                                
      IF (.NOT.firstcall) RETURN 
      !<-- set up wigner matrices                                       
                                                                        
      CALL fleur_d_wigner(                                              &
     &     sym%nop,sym%mrot,cell%bmat,3,                                &
     &     d_wgn)                                                       
                                                                        
      !>                                                                
      !<--local symmetry ops (see mapatom of FLEUR)                     
      invarind = 0 
      DO layer = 1,layers%num_layers 
         na    = 0 
         DO n  = 1,atoms(layer)%ntype 
            DO nt = 1,atoms(layer)%neq(n) 
               na = na+1 
               DO jop = 1,sym%nop 
                  gaminv = MATMUL(real(sym%mrot(:,:,jop)),atoms(layer)%taual(:&
     &                 ,na))+sym%tau(:,jop)-atoms(layer)%taual(:,na)    
                  DO j3 = -2,2 
                     sr(3) = gaminv(3) + REAL(j3) 
                     DO j2 = -2,2 
                        sr(2) = gaminv(2) + REAL(j2) 
                        DO j1 = -2,2 
                           sr(1) = gaminv(1) + REAL(j1) 
                                                          !distance in  
                           s3  = SQRT(dot_PRODUCT(sr,sr)) 
                                                          !interal units
                                                          !different fro
                           IF (s3<1.0E-6) THEN 
                              invarind(na,layer) = invarind(na,layer) + &
     &                             1                                    
                              invarop(na,invarind(na,layer),layer) = jop 
                           END IF 
                        END DO 
                     END DO 
                  ENDDO 
               ENDDO 
            ENDDO 
         ENDDO 
      ENDDO 
                                                                        
      !>                                                                
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S:gf_writedos                                                
                                                                        
      SUBROUTINE gf_writedos(mpi,layers,atoms,cell,gfinp,jspins,nkpts)
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      USE hdf5
      use m_hdf_tools
      USE m_gf_energies 
      USE m_gf_types 
      USE m_gf_math,ONLY:interpolate_analytic 
#include "juDFT_env.h"
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_mpi),INTENT(IN)  :: mpi
      TYPE(t_atoms),INTENT(IN)  :: atoms(:) 
      TYPE(t_cell),INTENT(IN)   :: cell(:) 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_gfinp),INTENT(IN)  :: gfinp 
      INTEGER,INTENT(IN)        :: jspins,nkpts 
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
                                                    !loops              
      INTEGER :: n,l,ll,jspin,layer 
                                                    !plot mode input    
      INTEGER :: proj,atom 
                                                    !gf_dos-file exists 
      LOGICAL :: l_gf_dos 
                                                    !no of LDOS         
      INTEGER :: ndos ,hdferr
                                                    !loop bounds        
      INTEGER :: l_max,l_min,n_max,n_min 
                                                    !flag for l-like    
      INTEGER,ALLOCATABLE :: l_dos(:,:,:) 
                                                    !flags for totdos   
      LOGICAL,ALLOCATABLE :: tot_dos(:,:) 
                                                    !collected raw dos  
      COMPLEX,ALLOCATABLE :: dos_raw(:,:,:) 
                                                      !diagonalized dos 
      COMPLEX,ALLOCATABLE :: dos_dia(:,:,:) 
                                                      !offdiagonal DOS  
      COMPLEX,ALLOCATABLE :: offdia(:,:) 
                                                    !energies for       
      REAL                :: xmin,xmax,dx,x 
                                                    !interpol           
                                                    !interpol. DOS      
      COMPLEX,ALLOCATABLE :: dos_inter(:,:) 
                                                                        
                                            !No of interpolation points 
      INTEGER,PARAMETER :: no_points = 1001 
      COMPLEX,ALLOCATABLE :: dos(:,:,:,:,:,:) 
      COMPLEX,ALLOCATABLE :: layer_dos(:,:,:) 
      REAL                :: energy 
      CHARACTER(LEN = 20) :: filename 
      NAMELIST /LDOS/ layer,atom,l,proj 
                                                                        
      !>                                                                
      ALLOCATE(l_dos(0:3,MAXVAL(atoms%ntype),layers%num_layers)) 
      ALLOCATE(tot_dos(MAXVAL(atoms%ntype),layers%num_layers)) 
                                                                        
                                                                        
      !<-- read gf_dos to determine what to write                       
                                                                        

      !first close the gf_setup.hdf file
      CPP_juDFT_timestart("Closing gf_dos.hdf")
      if (ldos_id.ne.0) then
         CALL io_dclose(ldos_id,n) 
         CALL io_dclose(intdos_id,n)
         IF (gfinp%l_surface) THEN 
            CALL io_dclose(surfdos_id,n) 
         ENDIF 
         CALL io_hdfclose(fid,n) 
      endif
      CPP_juDFT_timestop("Closing gf_dos.hdf")
      INQUIRE(FILE ="gf_dos",EXIST = l_gf_dos)
      l_dos= 0
      tot_dos=.FALSE.
      IF (.NOT.(l_gf_dos.and.mpi%pe0)) THEN
         RETURN 
      ELSE
         WRITE(*,*) "Writing DOS"
         !Open hdf-file and datasets
         CPP_juDFT_timestart("Opening gf_dos.hdf")
         call io_hdfopen("gf_dos.hdf",H5F_ACC_RDONLY_F,fid)
         call io_dopen(fid,"bands",ldos_id)

         CALL io_dopen(fid,"intbands",intdos_id,hdferr)
         call io_check("intdos?",hdferr)

         IF (gfinp%l_surface) THEN
            call io_dopen(fid,"surface",surfdos_id)
         endif
         CPP_juDFT_timestop("Opening gf_dos.hdf")
         OPEN(99,FILE="gf_dos",STATUS='old') 
         ndos=0 
         !<--loop over input until end of file                          
         DO 
            layer=1;atom = 0;l=0;proj=1 
            READ(99,LDOS,END=999,ERR=999) 
            IF (layer>layers%num_layers) THEN 
               layer = 1 

               WRITE(*,*) "Set layer = 1 in gf_dos!"
            ENDIF 
            IF (atom<1.OR.atom>atoms(layer)%ntype) THEN 
               n_min = 1 
               n_max=atoms(layer)%ntype 
            ELSE 
               n_min=atom 
               n_max=atom 
            ENDIF 
            DO atom=n_min,n_max 
               IF (l>3) THEN 
                  tot_dos(atom,layer) = .TRUE. 
                  ndos = ndos+1 
                  CYCLE 
               ENDIF 
               IF (l<0) THEN 
                  l_min=0 
                  l_max=3 
               ELSE 
                  l_min=l 
                  l_max=l 
               ENDIF 
               DO l=l_min,l_max 
                  !Add no of colums to dos-output                       
                  !Note: a later modification of the same l,n-mode will 
                  !result in a to large array, but should be corrected  
                  !later on                                             
                  IF (proj==1) THEN 
                     ndos=ndos+1 
                  ELSE IF(proj==2) THEN 
                     ndos=ndos+2*l+1 
                  ELSE IF(proj==3) THEN 
                     ndos=ndos+2*l+2 
                  ELSE IF(proj==4) THEN 
                     ndos=ndos+2*l+3 
                  ENDIF 
                  l_dos(l,atom,layer) = proj 
               ENDDO 
            ENDDO 
         ENDDO 
         !>                                                             
  999    CONTINUE 
         CLOSE(99) 
      ENDIF 
                                                                        
      !>                                                                
      !<-- construct the DOS                                            
      ALLOCATE(dos_raw(gf_noen(),jspins,ndos)) 
      ALLOCATE(DOS_dia(gf_noen(),jspins,0:7)) 
      ALLOCATE(offdia(gf_noen(),jspins)) 
      ALLOCATE(layer_dos(gf_noen(),jspins,layers%num_layers)) 
      layer_dos = 0.0 
      ndos = 0 
      DO layer = 1,layers%num_layers 
         ALLOCATE(dos(-3:3,-3:3,0:3,atoms(layer)%ntype,gf_noen()        &
     &        ,jspins))                                                 
         CALL priv_collectDOS(layer,dos,layer_dos(:,:,layer),gf_noen()  &
     &        ,nkpts,gfinp%l_band)                                      
         layer_dos(:,:,layer) = layer_dos(:,:,layer)*cell(layer)%volint 
         !WRITE(90+layer,*) dos(0,0,0,:,:,:)
         DO n = 1,atoms(layer)%ntype 
            IF (ALL(l_dos(:,n,layer)==0).AND..NOT.tot_dos(n,layer))     &
     &           CYCLE                                                  
               !<-- the total DOS for this atom                         
            IF (tot_dos(n,layer).OR.gfinp%l_intdos) THEN 
               CALL priv_totdos(n,dos_dia(:,:,0),dos) 
               IF (tot_dos(n,layer)) THEN 
                  ndos = ndos+1 
                  WRITE(6,'(a,i4,a,i3,a)') "Column "                    &
     &                 ,ndos," : atom =",n," totaldos"                  
                  dos_raw(:,:jspins,ndos) = dos_dia(:,:jspins,0) 
               ENDIF 
               layer_dos(:,:jspins,layer) = layer_dos(:,:jspin,layer)+dos_dia(:,:jspins,0)
            ENDIF 
            !>                                                          
            !<-- the l-decomposed dos                                   
            DO l = 0,3 
               IF (l_dos(l,n,layer)>0) THEN 
                  CALL priv_dosdia(n,l,dos_dia,offdia,dos               &
     &                 )                                                
                  IF (l_dos(l,n,layer) == 1.OR.l_dos(l,n,layer)>2)      &
     &                 THEN                                             
                     ndos = ndos+1 
                     WRITE(6,'(a,i4,a,i3,a,i1)')                        &
     &                    "Column ",ndos," : atom =",n,",l =",l         
                                                                        
                     dos_raw(:,:jspins,ndos) = dos_dia(:,:jspins,0) 
                  ENDIF 
                  IF (l_dos(l,n,layer)>1) THEN 
                     DO ll = 1,2*l+1 
                        ndos = ndos+1 
                        IF (jspin == 1) WRITE(6,'(a,i4,a,i3,a,i1,a)')   &
     &                       "Column ",ndos," : atom =",n,",l =",l      &
     &                       ,"sym"                                     
                        dos_raw(:,:jspins,ndos) = dos_dia(:,:jspins     &
     &                       ,ll)                                       
                     ENDDO 
                  ENDIF 
                  IF (l_dos(l,n,layer) == 4) THEN 
                     ndos           = ndos+1
                     IF (jspin == 1) WRITE(6,'(a,i4,a,i3,a,i1,a)') "Column ",ndos," : atom =",n,",l =",l,"Offdia"
                     dos_raw(:,:jspins,ndos) = offdia(:,:jspins) 
                  ENDIF 
               ENDIF 
            ENDDO 
         ENDDO 
         !>                                                             
         DEALLOCATE(dos) 
      ENDDO 
      DEALLOCATE(DOS_dia) 
                                                                        
      !>                                                                
      !<-- Writeout the gf_DOS.raw file                                 
                                                                        
      OPEN(99,FILE ="gf_DOS.raw",FORM ='formatted') 
      DO n = 1,gf_NoEn() 
         WRITE(99,'(99(f0.5,1x))') REAL(gf_Z(n,0)),                     &
     &        ((2*(jspin-1.5)*AIMAG(dos_raw(n,jspin,l)),l = 1,ndos)     &
     &        ,jspin                                      = 1,jspins)   
      ENDDO 
      CLOSE(99) 
                                                                        
      !>                                                                
                                                                        
      !<-- Writeout the gf_layerdos.raw file                            
                                                                        
      DO layer=1,layers%num_layers 
         WRITE(filename,"(a,i0.5)") "gf_layerdos.",layer 
         OPEN(99,FILE =filename,FORM ='formatted') 
         DO n = 1,gf_NoEn() 
            WRITE(99,'(99(f0.5,1x))') REAL(gf_Z(n,0)),                  &
     &           (2*(jspin-1.5)*(aimag(layer_dos(n,jspin,layer)))              &
     &           ,jspin = 1,jspins)                                     
         ENDDO 
         CLOSE(99) 
         !<-- Calculate the eigenvalue sum for the total energy         
         energy = AIMAG(SUM(gf_allz(layer)*layer_dos(:,1,layer))) 
         IF (jspins == 2) THEN 
            energy  = energy+AIMAG(SUM(gf_allz(layer)*layer_dos(:,2     &
     &           ,layer)))                                              
         ELSE 
            energy = energy*2 
         ENDIF 
         WRITE(6,*) "Eigenvalue sum:",energy,energy/cell(layer)%vol 
         !>                                                             
                                                                        
      ENDDO 
      !>                                                                
                                                                        
                                                                        
                                                                        
                                                                        
      !<-- energy-values for interpolation                              
      xmin = MINVAL(REAL(gf_allz(0)))-0.1 
      xmax = MAXVAL(REAL(gf_allz(0)))+0.1 
      dx = (xmax-xmin)/(1.*no_points) 
      !>                                                                
      !<-- calculate interpolated values                                
                                                                        
      OPEN(99,FILE ="gf_DOS.int",FORM ='formatted') 
      ALLOCATE(dos_inter(ndos,jspins)) 
      DO n = 1,no_points 
         x = xmin+n*dx 
         DO jspin = 1,jspins 
            DO l = 1,ndos 
               dos_inter(l,jspin) = interpolate_analytic(dos_raw(:,jspin&
     &              ,l),gf_ALLz(0),CMPLX(x,0.0))                        
            ENDDO 
         ENDDO 
         !OUTPUT                                                        
         WRITE(99,'(99(f0.6,1x))') x,((2*(jspin-1.5)*AIMAG(dos_inter(l  &
     &        ,jspin)),l = 1,ndos),jspin = 1,jspins)                    
      ENDDO 
      CLOSE(99) 
                                                                        
      !>                                                                
      DEALLOCATE(dos_inter,dos_raw) 
      DEALLOCATE(l_dos) 
      DEALLOCATE(tot_dos) 
      !<-- write planar_dos                                             
                                                                        
!      IF (gfinp%l_surface) THEN                                        
!         call priv_writesurface(gf_noen(),jspins,nkpts,gfinp%l_band)   
!      ENDIF                                                            
      CALL io_dclose(ldos_id,n) 
      CALL io_dclose(intdos_id,n)
      IF (gfinp%l_surface) THEN 
         CALL io_dclose(surfdos_id,n) 
      ENDIF 
      CALL io_hdfclose(fid,n) 
                                                                        
      !>                                                                
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S:gf_dos_int(layer,mpi,lapw,stars,jspin,omtil,G,en,nk)       
      RECURSIVE SUBROUTINE gf_dos_INT(gfinp,layer,mpi,lapw,stars,jspin,omtil,G,wk &
     &     ,en,nk,l_noco)
!-----------------------------------------------                        
!                                                                       
!           (last modified:10-03-04) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_constants,ONLY:pimach 
      USE m_gf_stepsanaly 
      USE m_hdf_tools 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      TYPE(t_gfinp),INTENT(IN) :: gfinp 
      INTEGER,INTENT(IN)      :: layer 
      TYPE(t_mpi),INTENT(IN)  :: mpi 
      TYPE(t_lapw),INTENT(IN) :: lapw 
      INTEGER, INTENT (IN)    :: jspin 
      REAL,    INTENT (IN)    :: omtil 
      TYPE(t_stars),INTENT(IN) :: stars 
      COMPLEX, INTENT (IN)    :: G(:,:) 
      REAL   ,INTENT(IN)      :: wk 
      INTEGER,INTENT(IN)      :: en,nk
      logical,intent(in)      :: l_noco
      !>                                                                
      !<-- Locals                                                       
      REAL    :: pi 
      INTEGER :: n1,n2,i1,i2,i3,in 
      COMPLEX :: dg,phase 
      COMPLEX,ALLOCATABLE :: rg(:,:,:) 
      COMPLEX,ALLOCATABLE :: qpw(:),qpw_w(:) 
      !>

      !In noco-case call subroutine recursively
      if (l_noco) THEN
        call gf_dos_int(gfinp,layer,mpi,lapw,stars,1,omtil,    &
            G(:size(G,1)/2,:size(g,2)/2),                      &
            wk,en,nk,.false.)
     call gf_dos_int(gfinp,layer,mpi,lapw,stars,2,omtil,    &
            G(size(G,1)/2+1:,size(g,2)/2+1:),                      &
            wk,en,nk,.false.)
      endif


      pi = pimach() 
      ALLOCATE(rg(-stars%mx1:stars%mx1,-stars%mx2:stars%mx2,            &
     &     -stars%mx3:stars%mx3))                                       
      ALLOCATE(qpw(stars%nq3),qpw_w(stars%nq3)) 
      qpw = 0.0 
      rg=CMPLX(0.0,0.0) 
      !<-- Calculate the reduced G-function                             
      DO n1=1,lapw%nv_sphere(jspin)
!         IF (ABS(lapw%k%k3(n1,jspin))>stars%mx3) CYCLE                 
         DO n2=1,lapw%nv_sphere(jspin)
!            IF (ABS(lapw%k%k3(n2,jspin))>stars%mx3) CYCLE              
            i1 = lapw%k%k1(n2,jspin) - lapw%k%k1(n1,jspin) 
            IF (iabs(i1)>stars%mx1) CYCLE 
            i2 = lapw%k%k2(n2,jspin) - lapw%k%k2(n1,jspin) 
            IF (iabs(i2)>stars%mx2) CYCLE 
            i3 = lapw%k%k3(n2,jspin) - lapw%k%k3(n1,jspin) 
            IF (iabs(i3)>stars%mx3) CYCLE 
            rg(i1,i2,i3)=rg(i1,i2,i3)+g(n2,n1) 
         ENDDO 
      ENDDO 
      !>                                                                
      !<--  Now calculate the density for all stars!                    
      DO i1=-stars%mx1,stars%mx1 
         DO i2=-stars%mx2,stars%mx2 
            DO i3=-stars%mx3,stars%mx3 
               in=stars%ig(i1,i2,i3) 
               IF (in==0) CYCLE 
               phase = stars%rgphs(i1,i2,i3)/ (stars%nstr(in)*omtil)/pi 
                                                            !factor 1/Pi
               dg=CMPLX(0.5,0.)*(rg(i1,i2,i3)-CONJG(rg(-i1,-i2,-i3))) 
               qpw(in)=qpw(in)+                                         &
     &              phase*dg                                            
               ! At this point one could also calculate the coef for -i1
               ! might be faster??? But doublecounting has to be avoided
            ENDDO 
         ENDDO 
      ENDDO 
      !>                                                                
      !<-- Convolute with step-function                                 
      CALL gf_initstepsanaly(stars,gfinp%napw(layer)) 
      CALL gf_gspaceconvolve(layer,stars,0.0,qpw,qpw_w)
                                                                        
      CALL io_WRITE(intdos_id,(/1,en,nk,jspin,layer/),(/1,1             &
     &     ,1,1,1/),wk*REAL(qpw_w(1)))                                  
      CALL io_WRITE(intdos_id,(/2,en,nk,jspin,layer/),(/1,1             &
     &     ,1,1,1/),wk*AIMAG(qpw_w(1)))                                 
      !>                                                                
                                                                        
      DEALLOCATE(qpw,qpw_w,rg) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- Private subroutines of the module                            
                                                                        
      !<-- S: priv_collectDOS                                           
                                                                        
      SUBROUTINE priv_collectDOS(layer,dos,int_dos,n_energies,nkpts     &
     &     ,l_band)                                                     
!******************************************                             
!     collects the DOS from all kpts stored in hdf-file                 
!     D. Wortmann                                                       
!******************************************                             
      USE m_gf_energies 
      use m_constants,only:pimach
      IMPLICIT NONE 
      INTEGER,INTENT(IN)     :: layer 
      COMPLEX,INTENT(OUT)    :: dos(:,:,:,:,:,:) 
      COMPLEX,INTENT(OUT)    :: int_dos(:,:) 
      INTEGER,INTENT(IN)     :: n_energies,nkpts 
      LOGICAL,INTENT(IN)     :: l_band 
                                                                        
      INTEGER               :: jsp,nk,en,l,m,nt,hdferr 
      REAL   ,ALLOCATABLE   :: tmp(:,:,:,:,:)
      complex,allocatable   :: band(:)
      CHARACTER(LEN = 15)   :: filename
      complex               :: c_tmp(n_energies)
      REAL                  :: s_tmp(2) 
                                                                        
      ALLOCATE(tmp(2,-3:3,-3:3,0:3,SIZE(dos,4))) 
      dos = 0 
      int_dos = 0.0 
                                                                        
      IF (ldos_id == 0) THEN 
         CALL io_hdfopen("gf_dos.hdf",H5F_ACC_RDONLY_F,fid,hdferr) 
         CALL io_dopen(fid,"bands",ldos_id,hdferr) 
         CALL io_dopen(fid,"intbands",intdos_id,hdferr) 
         IF (io_dataexists(fid,"surface")) CALL io_dopen(fid,"surface" &
     &        ,surfdos_id,hdferr)                                       
      ENDIF 
                                                                        
      IF (l_band) THEN 
         WRITE(filename,"(a,i0)") "gf_band_",layer 
         OPEN(99,FILE = filename,FORM ='formatted') 
         ALLOCATE(band(SIZE(dos,4))) 
      ENDIF 
      DO jsp = 1,SIZE(dos,6) 
         DO nk = 1,nkpts 
            call io_read(intdos_id,(/-1,1,nk,jsp,layer/),(/2,n_energies,1,1,1/),c_tmp)
            int_dos(:,jsp)=int_dos(:,jsp)+c_tmp
            DO en = 1,n_energies 
               !Read for this energy                                    
               CALL io_READ(ldos_id,(/1,1,1,1,1,en,nk,jsp,layer/),(/1,7 &
     &              ,7,4,SIZE(dos,4),1,1,1,1/),tmp(1,:,:,:,:))          
               CALL io_READ(ldos_id,(/2,1,1,1,1,en,nk,jsp,layer/),(/1,7 &
     &              ,7,4,SIZE(dos,4),1,1,1,1/),tmp(2,:,:,:,:))          

               !CALL io_READ(intdos_id,(/1,en,nk,jsp,layer/),(/2,1,1,1,1/),s_tmp)
               !int_dos(en,jsp) = int_dos(en,jsp)+CMPLX(s_tmp(1),s_tmp(2)&
     !&              )

               IF (l_band) THEN 
                                        !atoms%ntype                    
                  DO nt = 1,SIZE(dos,4) 
                     band(nt) = 0.0 
                                !sum over all l                         
                     DO l = 0,3 
                                   !trace over density matrix           
                        DO m =-l,l 
                           band(nt) = band(nt)+cmplx(tmp(1,m,m,l,nt),tmp(2,m,m,l,nt))
                                                                        
                        ENDDO 
                     ENDDO 
                  ENDDO 
                  WRITE(99,"(2i5,999(1x,f0.8))") jsp,nk,REAL(gf_Z(en,0))&
     &                 ,band/(cmplx(0,1.)*pimach())
               ENDIF 
               dos(:,:,:,:,en,jsp) = dos(:,:,:,:,en,jsp)+CMPLX(tmp(1,:,:&
     &              ,:,:),tmp(2,:,:,:,:))                               
            ENDDO 
            if (l_band) write(99,*)
         ENDDO 
      ENDDO 
      DEALLOCATE(tmp) 
      IF (l_band) THEN 
         DEALLOCATE(band) 
         CLOSE(99) 
      ENDIF 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: priv_collectDOS                                           
      SUBROUTINE priv_writesurface(n_energies,jspins,nkpts,l_band) 
!******************************************                             
!     collects the surface-DOS from all kpts stored in hdf-file         
!     and writes formatted files                                        
!     D. Wortmann                                                       
!******************************************                             
      USE m_gf_energies 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)     :: n_energies,jspins,nkpts 
      LOGICAL,INTENT(IN)     :: l_band 
                                                                        
      INTEGER               :: jsp,nk,en,l,ll,lll 
      REAL   ,ALLOCATABLE   :: dos(:,:,:),dosbox(:,:,:) 
                                                                        
      ALLOCATE(dos(nx,ny,nz),dosbox(nx,ny,nz)) 
                                                                        
      IF (l_band) THEN 
         OPEN(99,FILE ="gf_surface_band",FORM ='formatted') 
      ENDIF 
      OPEN(98,FILE ="gf_surface_dos",FORM ='formatted') 
      OPEN(97,FILE ="gf_zdos",FORM ='formatted') 
      DO jsp = 1,jspins 
         DO en = 1,n_energies 
            dos =0.0 
            DO nk = 1,nkpts 
               !Read for this energy                                    
               CALL io_read(surfdos_id,(/2,1,1,1,en,nk,jsp/),(/1,nx     &
     &              ,ny,nz,1,1,1,1/),dos)                               
               IF (l_band) THEN 
                  DO l = 1,nx 
                     DO ll = 1,ny 
                        DO lll = 1,nz 
                           WRITE(99,"(1i5,1x,f0.8,3i5,1x,f0.8)") jsp    &
     &                          ,REAL(gf_Z(en,0)),l,ll,lll,dosbox(l,ll  &
     &                          ,lll)                                   
                        ENDDO 
                     ENDDO 
                  ENDDO 
               ENDIF 
               dos = dosbox+dos 
                  !nk                                                   
            ENDDO 
            DO l = 1,nx 
               DO ll = 1,ny 
                  DO lll = 1,nz 
                     WRITE(98,"(i5,1x,f0.8,3i5,999(1x,f0.8))") jsp      &
     &                    ,REAL(gf_Z(en,0)),dos(l,ll,lll)               
                  ENDDO 
               ENDDO 
            ENDDO 
            DO lll = 1,nz 
               WRITE(97,"(i5,1x,f0.8,i5,999(1x,f0.8))") jsp,REAL(gf_Z(en&
     &              ,0)),lll,SUM(dos(:,:,lll))                          
            ENDDO 
         ENDDO 
      ENDDO 
      CLOSE(97) 
      CLOSE(97) 
      IF (l_band) CLOSE(99) 
      DEALLOCATE(dos,dosbox) 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: priv_totdos                                               
      SUBROUTINE priv_totdos(n,tot_dos,dos) 
!******************************************                             
!     Calculate total DOS for atom n                                    
!                          D. Wortmann                                  
!******************************************                             
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)  :: n 
      COMPLEX,INTENT(OUT) :: tot_dos(:,:) 
      COMPLEX,INTENT(IN)  :: dos(-3:,-3:,0:,:,:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER :: l,m,jspin 
      !>                                                                
      tot_dos = 0.0 
      DO jspin = 1,SIZE(tot_dos,2) 
                     !sum over all l                                    
         DO l  = 0,3 
                       !trace over density matrix                       
            DO m =-l,l 
               tot_dos(:,jspin) = tot_dos(:,jspin)+dos(m,m,l,n,:,jspin) 
            ENDDO 
         ENDDO 
      ENDDO 
      END SUBROUTINE 
      !>                                                                
      !<-- S: priv_dosdia                                               
                                                                        
      SUBROUTINE priv_dosdia(n,l,dos_dia,offdia,dos) 
!******************************************                             
!     Diagonalizes LDOS for give atom n and l                           
!                          D. Wortmann                                  
!******************************************                             
      USE m_gf_math 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN) :: n,l 
      COMPLEX,INTENT(OUT) :: dos_dia(:,:,0:) 
      COMPLEX,INTENT(OUT) :: offdia(:,:) 
      COMPLEX,INTENT(IN)  :: dos(-3:,-3:,0:,:,:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER :: m,mp,en,jspin 
      COMPLEX :: dos_tmp(-l:l,-l:l),dos_tmp1(-l:l,-l:l) 
      COMPLEX :: ew(-l:l),A(-l:l,-l:l),AINV(-l:l,-l:l) 
      !>                                                                
      dos_dia=0.0 
      offdia=0.0 
      !<--first create transformation matrix                            
      DO jspin = 1,SIZE(offdia,2) 
         dos_tmp = 0.0 
         DO m =-l,l 
            DO mp =-l,l 
               dos_tmp(m,mp) = SUM(dos(m,mp,l,n,:,jspin)) 
            ENDDO 
         ENDDO 
         !diagonalize dos_tmp                                           
         CALL eigenvalues(dos_tmp,ew,A) 
         IF (l == 2) THEN 
            A = 0.0 
            A(-2,-2) = 1 
            A(2,-2) = 1 
            A(-2,-1) =-1 
            A(2,-1) = 1 
            A(-1,0) = 1 
            A(1,0) = 1 
            A(-1,1) =-1 
            A(1,1) = 1 
            A(0,2) = SQRT(2.) 
            A      = A/SQRT(2.) 
         ENDIF 
         AINV=mat_inverse(A) 
                                                                        
         !>                                                             
         !<-- loop over energies and transform                          
         DO en = 1,SIZE(dos_dia,1) 
            dos_tmp = MATMUL(dos(-l:l,-l:l,l,n,en,jspin),A) 
            dos_tmp1 = MATMUL(AINV,dos_tmp) 
            DO m =-l,l 
               ! write(*,*) m,A(m,:),dos_tmp1(m,:)                      
               dos_dia(en,jspin,m+l+1) = dos_tmp1(m,m) 
               dos_dia(en,jspin,0)     = dos_dia(en,jspin,0)+dos_tmp1(m &
     &              ,m)                                                 
               !calcualate square of off-diagonal DOS                   
               DO mp =-l,l 
                  IF (m == mp) CYCLE 
                  offdia(en,jspin) = offdia(en,jspin)+CMPLX(0           &
     &                 ,AIMAG(dos_tmp1(m,mp))**2)                       
               ENDDO 
            ENDDO 
         ENDDO 
         !>                                                             
      ENDDO 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
                                                                        
      !>                                                                
      END                                           
