!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_hsdata 
      use m_juDFT
!-----------------------------------------------                        
! DESC:Module to store, handle and communicate the FLEUR Hamiltonian    
!                 Daniel Wortmann, (05-08-24)                           
!-----------------------------------------------                        
#include "realcomplex.h"                                                
#include"cpp_double.h"                                                  
                                                                        
      use m_juDFT 
      IMPLICIT NONE
      PRIVATE 
      COMPLEX,ALLOCATABLE,TARGET,SAVE :: H(:),S(:)
      INTEGER,SAVE                        :: spin,k,stored_layer 
      PUBLIC gf_assignHS,gf_getlargeH_eS,gf_writeHS,gf_readHS           &
     &     ,gf_collectHS,gf_getHSaddr                                   
      CONTAINS 
      !<-- S: gf_getHSaddr(H,S)                                         
      SUBROUTINE gf_getHSaddr(layer,HP,SP) 
!-----------------------------------------------                        
!     get address of H and S                                            
!-----------------------------------------------                        
      IMPLICIT NONE 
      INTEGER,INTENT(IN)     :: layer 
      COMPLEX,POINTER:: hp(:),sp(:) 
                                                                        
      IF ((layer /= stored_layer).OR..NOT.allocated(h)) THEN
         IF (.NOT.gf_readHS(layer,spin,k)) THEN
             write(*,*) "Stored layer:",stored_layer
             write(*,*) "Required:",layer
             CALL juDFT_error(              &
     &        "Did not find the Data of the FLAPW-Hamiltonian")
         endif
      ENDIF 
      IF(.NOT.allocated(h)) CALL juDFT_error("gf_getHSaddr",calledby="gf_hsdata.F90")
      IF(.NOT.allocated(s)) CALL juDFT_error("gf_getHSaddr",calledby="gf_hsdata.F90")
      hp=>h 
      sp=>s 
      END SUBROUTINE 
                                                                        
      !<-- S: gf_assignHS(matsize,jspin,nk,H,S)                         
                                                                        
      SUBROUTINE gf_assignHS(matsize,jspin,nk,layer,Hp,Sp) 
!-----------------------------------------------                        
!     assignes the two pointers Hp,Sp to point to H and S               
!           (last modified: 05-08-24) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types,ONLY:t_LAPW 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)                  :: matsize,jspin,nk,layer 
      COMPLEX,POINTER             :: hp(:),sp(:) 
      !>                                                                
      INTEGER             :: err 
      !check if H,S can be reused                                       
      IF (allocated(H)) THEN
         IF (SIZE(H) /= matsize)  DEALLOCATE(H,s) 
      ENDIF 
      IF (.NOT.allocated(H)) THEN
         ALLOCATE(H(matsize),S(matsize),STAT=err) 
         IF (err /= 0) CALL juDFT_error                                   &
     &        ("Not enough memory for Hamiltonian")                     
      ENDIF 
      spin            = jspin 
      k               = nk 
      stored_layer    = layer 
      hp => H 
      sp => S 
      END SUBROUTINE 
                                                                        
      !>                                                                
      !<-- S: gf_getlargeH_eS(jspin,nk,H,S)                             
      SUBROUTINE gf_getlargeH_eS(layer,jspin,nk,en,Hinv) 
!-----------------------------------------------                        
!     Constructs the large 2D G^-1 = (e*S-H) matrix                     
!           (last modified: 05-08-24) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_energies,ONLY:gf_z 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)     :: layer,jspin,nk,en 
      COMPLEX,INTENT(OUT)    :: Hinv(:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: in,n1,n2 
      COMPLEX             :: z 
      !>                                                                
                                                                        
      !<-- Check if correct h,s in memory, otherwise try to read it     
                                                                        
      IF (.NOT.allocated(h).OR..NOT.(layer == stored_layer.AND.spin == &
     &     jspin.AND.k == nk)) THEN                                     
         IF (.NOT.gf_readHS(layer,jspin,nk)) CALL juDFT_error             &
     &        ("Did not find the Data of the FLAPW-Hamiltonian")        
      ENDIF 
      !>                                                                
      !<--Check if enough data is present                               
      IF ((SIZE(hinv,1)*(SIZE(hinv,1)+1))/2>SIZE(h)) CALL               &
     &     juDFT_error("Not enough data in FLAPW-Hamiltonian")            
      !>                                                                
!     The matrixes have been hermitian up to now which allowed to store 
!     them in packed storage. Now expand to full storage to use complex 
!     energies!                                                         
                                                                        
                       !The energy                                      
      z=gf_z(en,layer) 
              !the linear index of the matrix elements                  
      in = 0 
      DO n1  = 1, SIZE(Hinv,1) 
         DO n2 = 1, n1 - 1 
            in = in+1 
            Hinv(n1,n2)=z*s(in)-h(in) 
            Hinv(n2,n1) = z*CONJG(s(in))-CONJG(h(in)) 
         ENDDO 
         in=in+1 
         Hinv(n1,n2)=z*s(in)-h(in) 
      ENDDO 
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: gf_writeHS(savemem)                                       
      SUBROUTINE gf_writeHS(savemem) 
!-----------------------------------------------                        
!      writes the hamiltonian and Overlapp matrix to a unformated file  
!      The filename will contain the kpt&spin information               
!      so that all the different H and S used will be kept in the       
!      current directory.                                               
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      LOGICAL,INTENT(IN)     :: savemem 
      !<-- Locals                                                       
      INTEGER   :: n 
      !>                                                                
      OPEN(111,FILE = priv_makefilename(stored_layer,spin,k),FORM='unformatted',STATUS  ='REPLACE')
      WRITE(111) SIZE(h),stored_layer,spin,k 
      DO n = 1,SIZE(h) 
         WRITE(111) h(n) 
      ENDDO 
      DO n = 1,SIZE(h) 
         WRITE(111) s(n) 
      ENDDO 
      CLOSE(111) 
      !<--formatted output for debugging                                
      IF (.FALSE.) THEN 
         OPEN(111,FILE = priv_makefilename(stored_layer,spin,k),FORM='formatted',STATUS      ='REPLACE')
         WRITE(111,*) SIZE(h),stored_layer,spin,k 
         DO n = 1,SIZE(h) 
            WRITE(111,'(i6,4e15.5)') n,h(n),s(n) 
         ENDDO 
         CLOSE(111) 
      ENDIF 
      IF (savemem) THEN 
         DEALLOCATE(h,s) 
         !NULLIFY(h,s)
      ENDIF 
      END SUBROUTINE 
      !>                                                                
      !<-- F: gf_readHS(layer,jspin,nk)                                 
      FUNCTION gf_readHS(layer,jspin,nk) 
!-----------------------------------------------                        
!    reads the H and S matrices from the file                           
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: layer,jspin,nk 
      LOGICAL                :: gf_readHS 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n 
      !>                                                                
      OPEN(111,FILE = priv_makefilename(layer,jspin,nk),FORM='unformatted',STATUS  ='OLD',ERR = 100)
      READ(111,ERR = 100) n,stored_layer,spin,k 
      IF (spin/=jspin.OR.k/=nk.OR.layer/=stored_layer) CALL       &
     &     juDFT_error("Wrong data in FLEUR-Hamiltonian file")            
      !<-- The storage might have the wrong size or might not be aSSocia
      IF (.NOT.allocated(H)) THEN
         ALLOCATE(H(n),S(n)) 
      ELSEIF (n /= SIZE(h)) THEN 
         DEALLOCATE(H,S) 
         ALLOCATE(H(n),S(n)) 
      ENDIF 
      !>                                                                
      DO n = 1,SIZE(h) 
         READ(111,ERR = 100) h(n) 
      ENDDO 
      DO n = 1,SIZE(h) 
         READ(111,ERR = 100) s(n) 
      ENDDO 
      CLOSE(111,ERR = 100) 
      gf_readHS=.TRUE. 
      RETURN 
                        !something went wrong                           
  100 gf_readHS=.FALSE. 
                                                                        
      END FUNCTION 
      !>                                                                
      !<-- F: priv_makefilename(layer,jspin,nk)                         
      FUNCTION priv_makefilename(layer,jspin,nk) 
!-----------------------------------------------                        
!   generates a filename for the current jspin,nk                       
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)     :: layer,nk,jspin 
      CHARACTER(LEN = 15)       :: priv_makefilename 
      !>                                                                
#ifdef CPP_MPI                                                          
      WRITE(priv_makefilename,"(a,i4.0,i1,i4)") "HSmat",layer,jspin,nk 
#else                                                                   
      WRITE(priv_makefilename,"(a,i4.0,i1,i4)") "HSmat",layer 
#endif                                                                  
      END FUNCTION 
      !>                                                                
                                                                        
      !<-- S: gf_collectHS(mpi,nbasfcn,matsize)                         
      SUBROUTINE gf_collectHS(mpi,nbasfcn,matsize) 
!-----------------------------------------------                        
!  collects the data of H,S distributed over several PE such that all PE
!  have the full H,S afterwards                                         
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      use m_juDFT 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      TYPE(t_mpi),INTENT(IN) :: mpi 
      INTEGER,INTENT(IN)     :: nbasfcn,matsize 
      !>                                                                
      !<-- Locals                                                       
      COMPLEX,POINTER :: H_tmp(:),S_tmp(:) 
      INTEGER         :: n,i,start,lstart,lstop,pe,ierr
#ifndef CPP_MPI                                                         
      RETURN 
#else                                                                   
      INCLUDE 'mpif.h' 

      !>                                                                
      !<-- ALLOCATE large new H_tmp and S_tmp matrices                  
      ALLOCATE(H_tmp(matsize),S_tmp(matsize)) 
      !>                                                                
      !<--Loop over all columns                                         
      DO n = 1,nbasfcn 
         !<--Which PE holds this data?                                  
         pe = MOD(n,mpi%n_size)-1
         IF (pe<0) pe = mpi%n_size-1
         !>                                                             
         !<--Position in global array                                   
         start=((n-1)*n)/2+1 
         !>                                                             
         !<-- copy data from local to global array                      
         IF (pe == mpi%n_rank) THEN
            !<-- Find position in local array                           
            lstart = 1 
                                  !nrank=0,....      
            lstop  = mpi%n_rank+2
            DO WHILE(((lstop-lstart) /= n).AND. lstop <= SIZE(H)) 
               i      = lstop-lstart+mpi%n_size
               lstart = lstop 
               lstop  = lstop+i       
            ENDDO 
            IF (lstop>SIZE(h)) CALL juDFT_error                           &
     &           ("Distribution error in hsdata")                       
            !>                                                          
            h_tmp(start:start+n) = h(lstart:lstop) 
            S_tmp(start:start+n) = S(lstart:lstop) 
         ENDIF 
         !>                                                             
         !<-- communicate                                               
         !CALL MPI_BCAST( h_tmp(start:start+n),n,CPP_MPI_COMPLEX         &
     !&        ,pe,mpi%sub_comm,ierr)
         !CALL MPI_BCAST( S_tmp(start:start+n),n,CPP_MPI_COMPLEX         &
     !&        ,pe,mpi%sub_comm,ierr)
         !>                                                             
      ENDDO 
      !>                                                                
      DEALLOCATE(H,S) 
      call juDFT_error("Not implemented")
      !h => h_tmp
      !s => s_tmp
#endif                                                                  
      END SUBROUTINE 
      !>                                                                
      END                                           
