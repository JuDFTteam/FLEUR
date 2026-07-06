!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_psqpw 
      IMPLICIT NONE
!----------------------------------------------------------------       
!Module provides interface to FLEUR                                     
!EXTERNAL convn,convn_dim                                               
!USE psqpw                                                              
!----------------------------------------------------------------       
      PRIVATE 
      PUBLIC fleur_convn_setup,fleur_psqpw 
      CONTAINS 
      !<-- S:fleur_convn_setup(stars,atoms)                             
      SUBROUTINE fleur_convn_setup(stars,atoms) 
!-----------------------------------------------                        
!  Determines  convergence parameter for pseudo-charge                  
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_convn 
      USE m_convndim 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_stars),INTENT(IN)   :: stars 
      TYPE(t_atoms),INTENT(INOUT):: atoms 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: ncvd 
      !>                                                                
                                                                        
      CALL convn_DIM(                                                   &
     &     stars%gmax_inp,                                              &
     &     ncvd)                                                        
      ALLOCATE(atoms%ncv(atoms%ntype)) 
      CALL convn(                                                       &
     &     ncvd,atoms%ntype,stars%gmax_inp,atoms%rmt,                   &
     &     atoms%ncv)                                                   
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: fleur_psqpw(atoms,stars,sphhar,cell,qpw,rho,rht,psq)      
                                                                        
      SUBROUTINE fleur_psqpw(mpi_comm,atoms,stars,sphhar,cell,sym,qpw   &
     &     ,rho,psq,no_core)
!-----------------------------------------------                        
!calls psqpw from FLEUR                                                 
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_psqpw 
      USE m_constants 
      USE m_od_types
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)        :: mpi_comm 
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_stars),INTENT(IN)  :: stars 
      TYPE(t_sphhar),INTENT(IN) :: sphhar 
      TYPE(t_cell),INTENT(IN)   :: cell 
      TYPE(t_sym),INTENT(IN)    :: sym 
      REAL,INTENT(IN)           :: rho(:,0:,:,:) 
      COMPLEX,INTENT(IN)        :: qpw(:,:) 
                                                !pseudo-charge          
      COMPLEX,INTENT(OUT)       :: psq(:) 
      logical,intent(in),optional::no_core

                                                                        
      REAL,ALLOCATABLE :: rht(:,:,:) 
      real             :: zatom(size(atoms%zatom))
      !>                                                                
      TYPE (od_inp) :: odi 
      TYPE (od_sym) :: ods 
      INTEGER       :: irank,isize 
#ifdef CPP_MPI                                                          
      INCLUDE 'mpif.h' 
      INTEGER :: ierr(3) 
      CALL MPI_COMM_RANK (MPI_COMM,irank,ierr) 
      CALL MPI_COMM_SIZE (MPI_COMM,isize,ierr) 
#else                                                                   
      irank = 0;isize = 1 
#endif                                                                  
      odi%d1 = .FALSE. 
      zatom=atoms%zatom
      if (present(no_core)) then
        if (no_core) zatom=0.0
      endif
      ALLOCATE(rht(1,1,1)) 
      CALL psqpw(irank,isize,mpi_comm,                                  &
     &     atoms%ntype,SIZE(sphhar%clnu,3),SIZE(sphhar%clnu,2)-1        &
     &     ,stars%nq3,SIZE(atoms%rmsh,1),MAXVAL(atoms%lmax)             &
     &     ,SIZE(sphhar%clnu,1),1,atoms%ntype,stars%nq3,atoms%jri       &
     &     ,atoms%ntypsy,atoms%lmax,                                    &
     &     1,maxval(atoms%ncv)+1,                                       &
     &     atoms%rmt,atoms%rmsh,atoms%dx,zatom,atoms%volmts,      &
     &     stars%sk3,stars%nstr,cell%omtil,0.1                          &
     &     ,cell%z1, cell%area,.FALSE.,1,atoms%neq,1,0.1,sphhar%clnu    &
     &     ,sphhar%mlh,sphhar%nmem,sphhar%llh,sphhar%nlh,cell%bmat      &
     &     ,stars%kv3,stars%ig2,sym%nop,atoms%nat,sym%symor,sym%mrot    &
     &     ,sym%tau,atoms%taual,sym%invtab,atoms%ncv,qpw,rho,rht,0.0,0.0&
     &     ,odi,ods,.FALSE.,psq)                                        
      DEALLOCATE(rht) 
      END SUBROUTINE 
                                                                        
      !>                                                                
      END                                           
