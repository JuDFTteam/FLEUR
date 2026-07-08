!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_charge 
#ifdef CPP_MPI
      USE mpi
#endif
      USE m_constants, ONLY: oUnit
          IMPLICIT NONE
!-----------------------------------------------                        
! DESC:Mix and possibly write out the charge                            
!                 Daniel Wortmann, (06-05-12)                           
!-----------------------------------------------                        
      CONTAINS 
                                                                        
      !<-- S: gf_charge(jspins,stars,atoms,cell,mix,gfinp,sphhar,sym,pwd
                                                                        
      SUBROUTINE gf_charge(jspins,fmpi,stars,atoms,cell,mix,gfinp,noco   &
     &     ,sphhar,sym,pwd,rho,layer,fix)                               
!-----------------------------------------------                        
!    finalize charge density by writing to file and mixing etc          
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_fft 
      USE m_gf_checkdop 
      USE m_gf_iodop 
      USE m_gf_types 
      use m_juDFT 
      use m_gf_vacuum_charge
      IMPLICIT NONE 
                                                                        
      !<-- Arguments                                                    
                                                                        
      INTEGER,INTENT(IN)        :: jspins 
      TYPE(t_gfmpi),INTENT(IN)    :: fmpi 
      TYPE(t_stars),INTENT(IN)  :: stars 
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_gfmix),INTENT(IN)    :: mix 
      TYPE(t_embinp),INTENT(IN)  :: gfinp 
      TYPE(t_sphhar),INTENT(IN) :: sphhar 
      TYPE(t_noco),INTENT(IN)   :: noco 
      TYPE(t_sym),INTENT(IN)    :: sym 
      COMPLEX,INTENT(INOUT)     :: pwd(:,:) 
      REAL   ,INTENT(INOUT)     :: rho(:,0:,:,:) 
      INTEGER,INTENT(IN)        :: layer 
      REAL   ,INTENT(IN)        :: fix(0:) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      REAL                :: q_r(27*stars%mx1*stars%mx2*stars%mx3) 
      REAL                :: q_MAX(27*stars%mx1*stars%mx2*stars%mx3) 
      INTEGER             :: jspin,i,ii,iii,i4 
      LOGICAL,POINTER     :: step1(:),step2(:),step3(:) 
      LOGICAL             :: l_nolimit,l_fix 
                                                                        
      !>                                                                
                                                                        

                                                                        
      !check charge                                                     
      CALL gf_checkdop(atoms,cell,stars,sphhar,sym,                     &
     &     .TRUE.,pwd,rho)                                              
                                                                        
                                                                        
      !<-- enforce charge-neutrality                                    
                                                                        
      INQUIRE(FILE ="qfix",EXIST = l_fix) 
      IF (l_fix) THEN 
         OPEN(99,FILE ="qfix") 
         loop:DO 
            READ(99,ERR = 100,END = 100) i 
            IF (i == 0) EXIT loop 
            IF (i<0) THEN 
               i = layer 
               EXIT loop 
            ENDIF 
            IF (layer == i) EXIT loop 
         ENDDO loop 
         WRITE(oUnit,*) "Qfix:",fix(i) 
         rho = fix(i)*rho 
         pwd = fix(i)*pwd 
  100    CLOSE(99) 
      ENDIF 
                                                                        
      !>                                                                
      !<-- save charge                                                  
      CALL gf_renamepot(gf_cdnfile,fmpi%iodop_subcom,layer)
      CALL gf_wrtdop(GF_CDNFILE,layer,jspins,                           &
     &     gfinp,atoms,stars,sphhar,                                    &
     &     rho(:,0:,:,:),pwd(:,:),.FALSE.,fmpi%iodop_subcom)             
      if (gfinp%l_surface.and.layer==1) call gf_vacuum_writecharge(fmpi)
      !>                                                                
                                                                        
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: gf_charge_sum(qtot_nuc,qtot_el)                           
      SUBROUTINE gf_charge_sum(surface,qtot_nuc,qtot_el)
!-----------------------------------------------                        
!     sum up the total charges in all layers                            
!     for the parallel version also collect from all PE                 
!           (last modified: 08-03-28) D. Wortmann                       
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments
      logical,intent(in)        :: surface
      REAL   ,INTENT(INOUT)     :: qtot_nuc(0:) 
      REAL   ,INTENT(INOUT)     :: qtot_el(0:) 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: q_nuc 
      REAL                :: q_el 
      INTEGER             :: rank,e
      !>                                                                
#ifdef CPP_MPI                                                          
#endif                                                                  
      q_nuc = SUM(qtot_nuc) 
      q_el = SUM(qtot_el) 
      rank=0
                                                                        
                                                                        
#ifdef CPP_MPI                                                          
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,e) 
      IF (rank == 0) THEN 
         CALL MPI_REDUCE(MPI_IN_PLACE,q_nuc,1,                          &
     &        MPI_DOUBLE_PRECISION,                                             &
     &        MPI_SUM,0,MPI_COMM_WORLD,e)                               
         CALL MPI_REDUCE(MPI_IN_PLACE,q_el,1,                           &
     &        MPI_DOUBLE_PRECISION,                                             &
     &        MPI_SUM,0,MPI_COMM_WORLD,e)                               
      ELSE 
         CALL MPI_REDUCE(q_nuc,MPI_IN_PLACE,1,                          &
     &        MPI_DOUBLE_PRECISION,                                             &
     &        MPI_SUM,0,MPI_COMM_WORLD,e)                               
         CALL MPI_REDUCE(q_el,MPI_IN_PLACE,1,                           &
     &        MPI_DOUBLE_PRECISION,                                             &
     &        MPI_SUM,0,MPI_COMM_WORLD,e)                               
      ENDIF 
      CALL MPI_BCAST(q_nuc,1,                                           &
     &     MPI_DOUBLE_PRECISION,                                                &
     &     0,MPI_COMM_WORLD,e)                                          
      CALL MPI_BCAST(q_el,1,                                            &
     &     MPI_DOUBLE_PRECISION,                                                &
     &     0,MPI_COMM_WORLD,e)                                          
#endif                                                                  
      qtot_nuc(0) = q_nuc 
      qtot_el(0)  = q_el 


      END SUBROUTINE 
      !>                                                                
      END                                           
