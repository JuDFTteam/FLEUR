!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_mix 
#ifdef CPP_MPI
      USE mpi
#endif
      USE m_constants, ONLY: oUnit
!-------------------------------------------------------------          
!Module provides an interface to the mixing part of FLEUR               
!as used in the GF-code                                                 
!------------------------------------------------------------           
      USE m_gf_iodop 
      IMPLICIT NONE
      integer,parameter:: nmz=250,nmzxy=100
      PRIVATE 
      PUBLIC gf_mix 
      CONTAINS 
      !<--S: gf_mix(atoms,cell,sphhar,stars,gfinp,mix,jspins,layers)    
                                                                        
      SUBROUTINE gf_mix(atoms,cell,sphhar,stars,gfinp,noco,mpi,sym,     &
     &     enpara,mix,jspins,layers)
!*****************************************************************      
! DESC:Subroutine to mix the potentials in a self-consistent embedding  
! calculation. Mostly an adapted version of mix.F                       
! Two different potential-files are needed!                             
!  gf_pot contains the old potential                                    
!  gf_pot_new contains the new potential                                
! after mixing the file gf_pot will be overwritten by the mix potential 
!                                                                       
!                          Daniel Wortmann, Mon May 27 11:39:02 2002    
!*****************************************************************      
      USE m_gf_types 
      USE m_stmix 
      use m_juDFT 
      USE m_gf_broyden 
      USE m_gf_potdis 
      USE m_gf_metric 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_sphhar),INTENT(IN) :: sphhar(:) 
      TYPE(t_stars),INTENT(IN)  :: stars(:) 
      TYPE(t_atoms),INTENT(IN)  :: atoms(:) 
      TYPE(t_mpi),INTENT(IN)    :: mpi 
      TYPE(t_cell),INTENT(IN)   :: cell(:) 
      TYPE(t_noco),INTENT(IN)   :: noco(:) 
      TYPE(t_sym),INTENT(IN)    :: sym 
      type(t_enpara),intent(in) :: enpara(:)
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      TYPE(t_gfmix),INTENT(IN)    :: mix 
      TYPE(t_layers),INTENT(IN) :: layers 
      INTEGER,INTENT(IN)        :: jspins
                                                                        
      !<-- locals                                                       
      INTEGER               :: fileid 
      INTEGER               :: len,layer,mmap,layer_loop,l,ierr,ds,dl
      REAL                  :: volume,distance 
      INTEGER,DIMENSION(layers%num_layers) :: datastart,dataend 
      REAL, ALLOCATABLE                    :: sm(:),fsm(:),metric(:),sm_old(:)
                                                                        
      !>                                                                
#ifdef CPP_MPI                                                          
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 
#endif                                                                  
                                                                        
      IF (gfinp%l_potmix) THEN 
         fileID = GF_POTFILE 
      ELSE 
         fileID = GF_CDNFILE 
      ENDIF 
                                                                        
      IF (.NOT.(mpi%k_kpts(1) == 1)) RETURN 
      IF (.NOT.(mpi%pe0)) RETURN 
                                                                        
      !<-- Calculate distances                                          
                                                                        
      distance = 0.0;volume = 0.0 
      DO layer = 1,layers%num_layers 
         CALL gf_potdis(jspins,atoms(layer),stars(layer)                &
     &        ,sphhar(layer),mpi,sym,cell(layer),gfinp%l_potmix         &
     &        ,layer,distance,volume,enpara(layer))
      ENDDO 
      WRITE(oUnit,*) "Total distance of all layers:",1000*SQRT(distance     &
     &     /volume)                                                     
                                                                        
      !>                                                                
                                                                        
      !<-- print info on mixing scheme                                  
                                                                        
      WRITE(16,*) "Mixing:",mix%imix 
      WRITE(16,*) "Alpha:",mix%alpha 
                                                                        
      IF (.NOT.ANY((/0,3,5,7/)  == mix%imix))                           &
     &     CALL juDFT_error('mix: imix =/= 0,3,5,7 ')                     
                                                                        
      !>                                                                
                                                                        
      !<-- Calculate length of mixing vector                            
                                                                        
      len = 0 
      DO layer = 1,layers%num_layers 
         len  = len+2*stars(layer)%ng3 +                                &
     &        atoms(layer)%ntype*(MAXVAL(sphhar(layer)%nlh)+1)          &
     &        *MAXVAL(atoms(layer)%jri)                                 
                                                                        
      ENDDO 
      if (gfinp%l_surface) len=len+nmz+2*nmz*stars(1)%ng2-1

      len = len*jspins 
                                                                        
      ALLOCATE (sm_old(len),sm(len),fsm(len),metric(len))
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- construct large mixing vector                                
      CALL  priv_makevector(fileid,layers,mpi,mix,sphhar,atoms,stars,cell,sym,jspins,&
     &     dataend,datastart,sm,fsm,gfinp%l_surface)

      !>
      sm_old=sm

                                                                        
      !<-- mixing of the densities                                      
                                                                        
      IF (mix%imix == 0) THEN 
         DO layer = 1,layers%num_layers 
            !straight mixing is done layer by layer                     
            mmap = dataend(layer)-datastart(layer)+1 
            CALL stmix(                                                 &
     &           mmap,jspins,0,                                         &
     &           mmap,mmap/jspins,jspins,noco(1)%l_noco,mix%alpha       &
     &           ,mix%spinf,fsm(datastart(layer):dataend(layer))        &
     &           ,sm(datastart(layer):dataend(layer)))                  
         ENDDO 
      ELSE 
         CALL gf_broyden("gf_broyd",gfinp%l_potmix,gfinp%l_surface,mpi,                 &
     &        mix%imix,mix%maxiter,mix%alpha,fsm,stars,atoms,sphhar     &
     &        ,cell,jspins,layers%num_layers,sm)                        
      END IF 
                                                                        
      !>                                                                
      CALL priv_decomposevector(fileid,layers,mpi,sphhar,atoms,stars    &
     &     ,gfinp,jspins,sm,gfinp%l_surface)
                                                                        
      !<-- Calculate the Distance                                       
      metric = gf_apply_metric(gfinp%l_surface,gfinp%l_potmix,mpi,stars,atoms,cell,sphhar   &
     &     ,jspins,layers%num_layers,fsm)
      write(oUnit,*) "Distance for all layers:"
      DO layer=1,layers%num_layers
          write(oUnit,"('Layer:',i5,' distance=',f0.8)") layer, 1000.*SQRT(ABS(dot_PRODUCT(metric(datastart(layer):dataend(layer)) &
     &     ,fsm(datastart(layer):dataend(layer))))/cell(layer)%vol)
          ds=datastart(layer)
          dl=dataend(layer)-datastart(layer)
          dl=dl/jspins
          write(oUnit,"('    Interstitial distance=',f0.8)")   &
                1000.*SQRT(ABS(dot_PRODUCT(metric(ds:ds+stars(layer)%ng3),fsm(ds:ds+stars(layer)%ng3)))/cell(layer)%volint)
          write(oUnit,"('    MT distance          =',f0.8)")   &
                1000.*SQRT(ABS(dot_PRODUCT(metric(ds+stars(layer)%ng3+1:ds+dl),fsm(ds+stars(layer)%ng3+1:ds+dl)))/(cell(layer)%vol-cell(layer)%volint))
          if (jspins==2) then
             ds=ds+dl+1
             write(oUnit,"('    Interstitial distance=',f0.8)")   &
                1000.*SQRT(ABS(dot_PRODUCT(metric(ds:ds+stars(layer)%ng3),fsm(ds:ds+stars(layer)%ng3)))/cell(layer)%volint)
             write(oUnit,"('    MT distance=          ',f0.8)")   &
                1000.*SQRT(ABS(dot_PRODUCT(metric(ds+stars(layer)%ng3+1:ds+dl),fsm(ds+stars(layer)%ng3+1:ds+dl)))/(cell(layer)%vol-cell(layer)%volint))
          endif
      enddo
      if (gfinp%l_surface) &
           write(oUnit,"('Vacuum:',5x,' distance=',f0.8)") 1000.*SQRT(ABS(dot_PRODUCT(metric(dataend(layers%num_layers:)) &
     &     ,fsm(dataend(layers%num_layers:)))))/cell(1)%area/10
      WRITE(oUnit,*) "Mix: total distance =", 1000.*SQRT(ABS(dot_PRODUCT(metric,fsm))/SUM(cell%vol))
      !>                                                                
      !<--Stepsize
      sm=sm-sm_old
      metric = gf_apply_metric(gfinp%l_surface,gfinp%l_potmix,mpi,stars,atoms,cell,sphhar   &
     &     ,jspins,layers%num_layers,sm)
      write(oUnit,*) "Stepsize for all layers:"
      DO layer=1,layers%num_layers
          write(oUnit,"('Layer:',i5,' step=',f0.8)") layer, 1000.*SQRT(ABS(dot_PRODUCT(metric(datastart(layer):dataend(layer)) &
     &     ,sm(datastart(layer):dataend(layer)))))/cell(layer)%vol
      enddo
      if (gfinp%l_surface) &
           write(oUnit,"('Vacuum:',5x,' step=',f0.8)") 1000.*SQRT(ABS(dot_PRODUCT(metric(dataend(layers%num_layers:)) &
     &     ,sm(dataend(layers%num_layers:)))))/cell(1)%area/10
      WRITE(oUnit,*) "Mix: total stepsize =", 1000.*SQRT(ABS(dot_PRODUCT(metric,sm))/SUM(cell%vol))

      DEALLOCATE (sm,fsm,metric,sm_old)
                                                                        
                                                                        
      RETURN 
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- S: priv_makevector(layers,sphhar,atoms,stars,jspins,old,datae
                                                                        
      SUBROUTINE priv_makevector(fileid,layers,mpi,mix,sphhar,atoms,stars,cell,sym   &
     &     ,jspins,dataend,datastart,sout,fsout,l_surface)
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_fft 
      USE m_gf_iodop 
      use m_gf_precond
      use m_gf_iodop
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)        :: fileid 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_mpi),INTENT(IN)    :: mpi 
      TYPE(t_sphhar),INTENT(IN) :: sphhar(:) 
      TYPE(t_atoms),INTENT(IN)  :: atoms(:) 
      TYPE(t_stars),INTENT(IN)  :: stars(:)
      type(t_cell),intent(in)   :: cell(:)
      type(t_sym),intent(in)    :: sym
      TYPE(t_gfmix),intent(in)    :: mix
      INTEGER,INTENT(IN)        :: jspins 
      INTEGER,INTENT(OUT)       :: dataend(:),datastart(:) 
      REAL   ,INTENT(OUT)       :: sout(:),fsout(:)
      logical,intent(in)        :: l_surface
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER             :: offset,layer,js,na,n,l,i 
      REAL   ,ALLOCATABLE :: rho(:,:,:,:,:),vz(:,:,:)
      COMPLEX,ALLOCATABLE :: vpw(:,:,:),vxy(:,:,:,:)
      LOGICAL             :: fixed(layers%num_layers+1),lexist
#ifdef CPP_MPI                                                          
#endif                                                                  
                                                                        
      !>                                                                
                                                                        
                                                                        
      !<-- check if some layers are kept fixed                          
                                                                        
      fixed = .FALSE. 
      OPEN(99,FILE ="fixed_layers") 
      DO 
         READ(99,*,END = 100,ERR=100) layer 
         fixed(layer) = .TRUE. 
         WRITE(*,*) "Fixed layer:",layer 
      ENDDO 
  100 CLOSE(99) 
      !>                                                                
                                                                        
      !<-- create a large vector                                        
                                                                        
      datastart = 0;dataend = 0;sout = 0.0 
      offset = 0 
      DO layer=1,layers%num_layers 
         datastart(layer) = offset+1 
         !<-- read charge                                               
                                                                        
         ALLOCATE(vpw(stars(layer)%ng3,jspins,2))
         ALLOCATE(rho(MAXVAL(atoms(layer)%jri),0:MAXVAL(sphhar(layer    &
     &        )%nlh),atoms(layer)%ntype,jspins,2))
         CALL gf_loddop(fileid,layer,jspins,                            &
     &        atoms(layer),stars(layer),sphhar(layer),                  &
     &        rho(:,:,:,:,1),vpw(:,:,1),old=.true.)
         CALL gf_loddop(fileid,layer,jspins,                            &
     &        atoms(layer),stars(layer),sphhar(layer),                  &
     &        rho(:,:,:,:,2),vpw(:,:,2),old=.false.)

         if (fixed(layer)) then
            rho(:,:,:,:,2)=0.0
            vpw(:,:,2)=0.0
         else
            rho(:,:,:,:,2) = rho(:,:,:,:,2)-rho(:,:,:,:,1)
            vpw(:,:,2)   = vpw(:,:,2)-vpw(:,:,1)
         endif

         !Now call preconditioner

         if (mix%precond>0) call gf_precond(layer,mpi,mix,sphhar(layer),atoms(layer),stars(layer),cell(layer),sym,jspins,vpw(:,:,2),rho(:,:,:,:,2))

         !>                                                             
         DO js = 1,jspins 
            !put plane-waves on real-space grid                         
            sout(offset+1:offset+stars(layer)%ng3) =                    &
     &           REAL(vpw(:stars(layer)%ng3,js,1))
            fsout(offset+1:offset+stars(layer)%ng3) =                    &
     &           REAL(vpw(:stars(layer)%ng3,js,2))
            offset = offset+stars(layer)%ng3 
            sout(offset+1:offset+stars(layer)%ng3) =                    &
     &           AIMAG(vpw(:stars(layer)%ng3,js,1))
            fsout(offset+1:offset+stars(layer)%ng3) =                    &
     &           AIMAG(vpw(:stars(layer)%ng3,js,2))
            offset = offset+stars(layer)%ng3 
            !put mt charge on grid                                      
            na = 1 
            DO n = 1,atoms(layer)%ntype 
               DO l = 0,sphhar(layer)%nlh(sym%ntypsy(na)) 
                  DO i = 1,atoms(layer)%jri(n) 
                     offset       = offset +1 
                     sout(offset) = rho(i,l,n,js,1)
                     fsout(offset) = rho(i,l,n,js,2)
                  END DO 
               END DO 
               na = na + atoms(layer)%neq(n) 
            END DO 
               !spins                                                   
         ENDDO 
         dataend(layer) = offset 
         DEALLOCATE(vpw,rho) 
            !layers                                                     
      ENDDO 

      !in case of a surface calculation we have to add the surface charge

      if (l_surface) then
         allocate(vz(nmz,jspins,2),vxy(nmzxy,stars(1)%ng2-1,jspins,2))
         call gf_iodop_readvacuum(GF_CDNFILE,vxy(:,:,:,1),vz(:,:,1),old=.true.)
         call gf_iodop_readvacuum(GF_CDNFILE,vxy(:,:,:,2),vz(:,:,2))
         if (fixed(layers%num_layers+1)) then
            vz(:,:,2)=0.0
            vxy(:,:,:,2)=0.0
         else
            vxy(:,:,:,2) = vxy(:,:,:,2)-vxy(:,:,:,1)
            vz(:,:,2)   = vz(:,:,2)-vz(:,:,1)
         endif
         DO js = 1,jspins
            sout(offset+1:offset+nmz) = vz(:,js,1)
            fsout(offset+1:offset+nmz) = vz(:,js,2)
            offset=offset+nmz
            DO n=1,stars(1)%ng2-1
                sout(offset+1:offset+nmzxy) = real(vxy(:,n,js,1))
                sout(offset+nmzxy+1:offset+2*nmzxy)= aimag(vxy(:,n,js,1))
                fsout(offset+1:offset+nmzxy) = real(vxy(:,n,js,2))
                fsout(offset+nmzxy+1:offset+2*nmzxy)= aimag(vxy(:,n,js,2))
                offset=offset+2*nmzxy
            enddo
        enddo
      endif

      !>                                                                
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: priv_decomposevector(layers,sphhar,atoms,stars,jspins,old,
                                                                        
      SUBROUTINE priv_decomposevector(fileid,layers,mpi,sphhar,atoms    &
     &     ,stars,gfinp,jspins,smix,l_surface)
!-----------------------------------------------                        
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_fft 
      USE m_gf_iodop 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)        :: fileid 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_mpi),INTENT(IN)    :: mpi 
      TYPE(t_sphhar),INTENT(IN) :: sphhar(:) 
      TYPE(t_atoms),INTENT(IN)  :: atoms(:) 
      TYPE(t_stars),INTENT(IN)  :: stars(:) 
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      INTEGER,INTENT(IN)        :: jspins 
      REAL   ,INTENT(IN)        :: smix(:) 
      logical,intent(in)        :: l_surface
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: offset,layer,js,na,n,l,i 
      REAL   ,ALLOCATABLE :: rho(:,:,:,:),vz(:,:)
      COMPLEX,ALLOCATABLE :: vpw(:,:),vxy(:,:,:)
      !>                                                                
                                                                        
                                                                        
      !<-- decompose large vector                                       
                                                                        
      offset = 0 
      DO layer = 1,layers%num_layers 
         ALLOCATE(vpw(stars(layer)%ng3,jspins)) 
         ALLOCATE(rho(MAXVAL(atoms(layer)%jri),0:MAXVAL(sphhar(layer    &
     &        )%nlh),atoms(layer)%ntype,jspins))                        
         DO js = 1,jspins 
            vpw(:stars(layer)%ng3,js) = CMPLX(smix(offset+1:offset      &
     &           +stars(layer)%ng3),smix(offset+1+stars(layer           &
     &           )%ng3:offset+2*stars(layer)%ng3))                      
            offset = offset+2*stars(layer)%ng3 
            !put mt charge on grid                                      
            na = 1 
            DO n = 1,atoms(layer)%ntype 
               DO l = 0,sphhar(layer)%nlh(sym%ntypsy(na)) 
                  DO i = 1,atoms(layer)%jri(n) 
                     offset = offset +1 
                     rho(i,l,n,js) = smix(offset) 
                  END DO 
               END DO 
               na = na + atoms(layer)%neq(n) 
            END DO 
               !spins                                                   
         ENDDO 
         CALL gf_renamepot(fileID,mpi%iodop_subcom,layer)
         CALL gf_wrtdop(fileid,layer,jspins,                            &
     &        gfinp,atoms(layer),stars(layer),sphhar(layer),            &
     &        rho,vpw,.FALSE.,mpi%self_subcom)                          
         DEALLOCATE(vpw,rho) 
            !layers                                                     
      ENDDO 
      if (l_surface) then
         allocate(vz(nmz,jspins),vxy(nmzxy,stars(1)%ng2-1,jspins))
         DO js = 1,jspins
            vz(:,js)=smix(offset+1:offset+nmz)
            offset=offset+nmz
            DO n=1,stars(1)%ng2-1
                vxy(:,n,js)=cmplx(smix(offset+1:offset+nmzxy),smix(offset+nmzxy+1:offset+2*nmzxy))
                offset=offset+2*nmzxy
            enddo
        enddo
        call gf_iodop_writevacuum(GF_CDNFILE,vxy(:,:,:),vz(:,:),mpi%self_subcom)
      endif

      !>                                                                
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
                                                                        
      END                                           
