!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_boundaryembpot 
          IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Save boundary value of Coulomb potential and pseudo charge to emb
!                 Daniel Wortmann, (08-07-10)                           
!-----------------------------------------------                        
      CONTAINS 
                                                                        
      !<-- S:gf_boundaryembpot(layers,stars,cell,gfinp)                 
      SUBROUTINE gf_boundaryembpot(layers,stars,cell,atoms,gfinp)
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-07-10) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_iodop 
      USE m_constants 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_stars),INTENT(IN)  :: stars(:) 
      TYPE(t_cell),INTENT(IN)   :: cell(:) 
      TYPE(t_atoms),INTENT(IN)  :: atoms(:)
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: pos 
      COMPLEX,ALLOCATABLE :: vb(:) 
      COMPLEX,ALLOCATABLE :: vpw(:) 
      COMPLEX,ALLOCATABLE :: psq(:,:) 
      INTEGER             :: side,layer 
      INTEGER             :: n2 
      REAL                :: dp(2),kg(2) 
      COMPLEX,ALLOCATABLE :: phase(:) 
      !>                                                                
                                                                        
      DO side = 1,2 
         IF (side == 2) THEN 
            layer = 1 
         ELSE 
            layer =layers%num_layers 
         ENDIF 
         !read pseudo charge and save it to embpot-file                 
         CALL gf_iodop_readpseudo(layer,psq) 
         !<-- shift charge by parallel vector                           
         ALLOCATE(phase(SIZE(psq,2))) 
         dp(1) = gfinp%dp1 
         dp(2) = gfinp%dp2 
         IF (side == 1) dp=-1*dp
         DO n2 = 1,stars(layer)%nq2 
            kg = stars(layer)%kv2(:,n2) 
            phase(n2) = EXP(CMPLX(0.0,2.*pi_const*dot_PRODUCT(dp,kg))) 
            psq(:,n2) =psq(:,n2)*phase(n2) 
         ENDDO 
         !>                                                             
         CALL priv_savepseudo(side,layers,atoms(layer),cell(layer),dp,psq)
         !read coulomb-pot                                              
         ALLOCATE(vpw(stars(layer)%nq3)) 
         CALL gf_lodcoul(GF_POTFILE,layer,vpw = vpw(:)) 
         ALLOCATE(vb(stars(layer)%nq2)) 
         CALL priv_calcboundary(side,layer,layers,stars(layer),vpw,pos  &
     &        ,vb)                                                      
         !vb = vb*phase
         CALL priv_saveboundary(side,layer,layers,vb,pos) 
         DEALLOCATE(vb,vpw,phase) 
      ENDDO 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:priv_calcboundary(side,layer,layers,stars,vpw,pos,vb)      
      SUBROUTINE priv_calcboundary(side,layer,layers,stars,vpw,pos,vb) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-07-10) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_embdesc 
      USE m_constants 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)        :: layer,side 
      TYPE(t_layers),INTENT(IN) :: layers 
      TYPE(t_stars),INTENT(IN)  :: stars 
      COMPLEX,INTENT(IN)        :: vpw(:) 
      REAL   ,INTENT(OUT)       :: pos 
      COMPLEX,INTENT(OUT)       :: vb(:) 
      !>                                                                
      !<-- Locals                                                       
      TYPE(t_embdesc) :: embdesc1,embdesc2 
      INTEGER         :: n,n2,index 
      COMPLEX         :: exp1 
      !>                                                                
      !<-- determine position of plane from embdesc                     
      CALL gf_readdesc_setup(layer,embdesc1,embdesc2) 
      IF (side == 1) THEN 
         IF (embdesc1%cut_atoms_out == 0) THEN 
            pos                                 = 0.0 
         ELSE 
            pos =                                                       &
     &           MAXVAL(ABS(embdesc1%atoms_pos(3,:embdesc1%cut_atoms_out&
     &           ,2)-embdesc1%atoms_rmt(:embdesc1%cut_atoms_out,2)))    
         ENDIF 
      ELSE 
         IF (embdesc2%cut_atoms_out == 0) THEN 
            pos = 0.0 
         ELSE 
            pos = MAXVAL(ABS(embdesc2%atoms_pos(3                       &
     &           ,:embdesc2%cut_atoms_out,2)                            &
     &           -embdesc2%atoms_rmt(:embdesc2%cut_atoms_out,2)))       
         ENDIF 
      ENDIF 
      CALL gf_deallocembdesc(embdesc1) 
      CALL gf_deallocembdesc(embdesc1) 
      !>                                                                
      write(*,*) "Warning matching plane for potential at c/2"
      pos=0.0
      vb = 0.0 
      IF (side == 1) THEN 
         exp1   = EXP(CMPLX(0.0,-2.*pi_const/layers%dt(layer)*(          &
     &        layers%c(layer)/2.-pos)))
      ELSE 
         exp1 = EXP(CMPLX(0.0,-2.*pi_const/layers%dt(layer)              &
     &        *(-layers%c(layer)/2.+pos)))
      ENDIF 
      write(*,*) side,exp1
      vb = 0.0 
      DO n2 = 1,stars%nq2 
         DO n =-stars%mx3,stars%mx3 
            index = stars%ig(stars%kv2(1,n2),stars%kv2(2,n2),n) 
            IF (index == 0) CYCLE 
            vb(n2) = vb(n2)+vpw(index)*exp1**n
         ENDDO 
      ENDDO 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:priv_saveboundary(side,layer,layers,vb,pos)                
      SUBROUTINE priv_saveboundary(side,layer,layers,vb,pos) 
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-07-10) D. Wortmann                        
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_hdf_tools 
      USE hdf5 
      USE m_gf_io2dmat 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)        :: layer,side 
      TYPE(t_layers),INTENT(IN) :: layers 
      COMPLEX,INTENT(IN)        :: vb(:) 
      REAL   ,INTENT(IN)        :: pos 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: hdferr 
      INTEGER(HID_T)      :: fid,gid 
      !>                                                                
      IF (side == 1) THEN 
         fid   = gf_io2dmatFID(IO2D_EMB,1,1,IO2D_WRITE) 
      ELSE 
         fid = gf_io2dmatFID(IO2D_EMB                                   &
     &        ,layers%num_layers,2,IO2D_WRITE)                          
      ENDIF 
      IF (fid ==-1) RETURN 
                                                                        
      CALL io_gopen(fid,"Boundary",gid,hdferr) 
      CALL io_WRITE_att(gid,"n2",SIZE(vb)) 
      CALL io_WRITE_att(gid,"bulk_real",REAL(vb)) 
      CALL io_write_att(gid,"bulk_imag",AIMAG(vb)) 
      CALL io_WRITE_att(gid,"pos",pos) 
      CALL io_gclose(gid,hdferr) 
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S:priv_savepseudo(side,layers,psq)                           
      SUBROUTINE priv_savepseudo(side,layers,atoms,cell,dp,psq)
!-----------------------------------------------                        
!                                                                       
!           (last modified:08-07-10) D. Wortmann                        
!-----------------------------------------------                        
      USE M_hdf_tools 
      USE hdf5 
      USE m_gf_io2dmat 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)        :: side 
      TYPE(t_layers),INTENT(IN) :: layers
      type(t_atoms),INTENT(IN)  :: atoms
      TYPE(t_cell),INTENT(IN)   :: cell
      REAL,INTENT(IN)           :: dp(2)
      COMPLEX,INTENT(IN)        :: psq(:,:) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      INTEGER(HID_T)      :: fid,gid,did 
      INTEGER             :: hdferr,layer 
      REAL                :: atpos(3,atoms%nat)
      REAL                :: rmt(atoms%nat)
      INTEGER             :: n,nn,na,natoms
                                                                        
      !>                                                                
                                                                        
      IF (side == 1) THEN 
         fid   = gf_io2dmatFID(IO2D_EMB,1,1,IO2D_WRITE) 
         layer=1 
      ELSE 
         fid = gf_io2dmatFID(IO2D_EMB,layers%num_layers,2,IO2D_WRITE) 
         layer=layers%num_layers 
      ENDIF 
                                                                        
      IF (fid ==-1) RETURN 
                                                                        
      CALL io_gcreate(fid,"Boundary",gid,hdferr) 
      CALL io_createvar(gid,"pseudo", H5T_NATIVE_DOUBLE,(/2,SIZE(psq,1) &
     &     ,SIZE(psq,2)/),did)                                          
      CALL io_WRITE(did,(/-1,1,1/),(/1,SIZE(psq,1),SIZE(psq,2)/),psq) 
      CALL io_dclose(did,hdferr) 
      CALL io_WRITE_att(gid,"c",layers%c(layer)) 
      CALL io_WRITE_att(gid,"dt",layers%dt(layer)) 

      !Save info on atoms
      na=0
      natoms=0
      DO n=1,atoms%ntype
        DO nn=1,atoms%neq(n)
           na=na+1
           if (abs(atoms%pos(3,na))-cell%c/2<atoms%rmt(n)) THEN
              natoms=natoms+1
              rmt(natoms)=atoms%rmt(n)
              atpos(1:2,natoms)=atoms%pos(1:2,na)+matmul(cell%amat(1:2,1:2),dp)
              atpos(3,natoms)=atoms%pos(3,na)
           ENDIF
        ENDDO
      ENDDO


      CALL io_createvar(gid,"rmt", H5T_NATIVE_DOUBLE,(/natoms/),did)
      CALL io_WRITE(did,(/1/),(/natoms/),rmt(:natoms))
      CALL io_dclose(did,hdferr)

      CALL io_createvar(gid,"atpos", H5T_NATIVE_DOUBLE,(/3,natoms/),did)
      CALL io_WRITE(did,(/1,1/),(/3,natoms/),atpos(:,:natoms))
      CALL io_dclose(did,hdferr)

      CALL io_gclose(gid,hdferr) 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      END                                           
