!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_cdntot 
          IMPLICIT NONE
!*****************************************************************      
! DESC:Calculate total charge in the system                             
!                          Daniel Wortmann, Sat Feb  7 21:47:56 2004    
!*****************************************************************      
      CONTAINS 
      !<-- S:                                                           
                                                                        
      SUBROUTINE gf_cdntot(layer,mpi,jspins,stars,cell,atoms,rho,qpw    &
     &     ,q_el,q_nuc)                                                 
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      USE m_constants, ONLY: pi_const, oUnit 
      USE m_intgr, ONLY : intgr3 
      USE m_gf_types 
      USE m_gf_stepsanaly 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      INTEGER,INTENT(IN)             :: layer 
      INTEGER, INTENT (IN)           :: jspins 
      TYPE(t_stars),INTENT(IN)       :: stars 
      TYPE(t_cell),INTENT(IN)        :: cell 
      TYPE(t_atoms),INTENT(IN)       :: atoms 
      TYPE(t_mpi),INTENT(IN)         :: mpi 
      COMPLEX,INTENT(IN)             :: qpw(:,:) 
      REAL   ,INTENT(IN)             :: rho(:,0:,:,:) 
      REAL   ,INTENT(OUT),OPTIONAL   :: q_nuc,q_el 
      !>                                                                
      !<-- Locals                                                       
                                                                        
                         !loops                                         
      INTEGER :: jspin,n 
                                                    !charges            
      REAL    :: q,qis,qistot,qtot,qmt(atoms%ntype) 
      REAL    :: w 
      COMPLEX :: warped(stars%ng3) 
      LOGICAL :: lexists 
      !>                                                                
                                                                        
      qtot=0.0 
      qistot=0.0 
      DO  jspin = 1,jspins 
         q = 0.E0 
         !<--mt charge                                                  
                                                                        
         DO n = 1,atoms%ntype 
            CALL intgr3(rho(:,0,n,jspin),atoms%rmsh(:,n)                &
     &           ,atoms%dx(n),atoms%jri(n),w)                           
            qmt(n) = w*SQRT(4*pi_const) 
            q = q + atoms%neq(n)*qmt(n) 
         ENDDO 
                                                                        
         !>                                                             
         !<--is region                                                  
                                                                        
         qis = 0. 
         CALL gf_initstepsanaly(stars,0) 
         CALL gf_gspaceconvolve(layer,stars,0.0,       &
     &        qpw(:,jspin),warped)                                      
         qis = real(warped(1))*cell%area*cell%amat(3,3) 
                                                                        
         !>                                                             
         qistot = qistot + qis 
         q = q + qis 
         WRITE (oUnit,FMT=8000) jspin,q,qis, (qmt(n),n=1,atoms%ntype) 
         qtot = qtot + q 
      ENDDO 
      !Sum up positive charge                                           
      q=0 
      DO n=1,atoms%ntype 
         q=q+atoms%zatom(n)*atoms%neq(n) 
      ENDDO 
                                                                        
      WRITE (oUnit,FMT=8020) qtot,q,q-qtot 
                                                                        
 8000 FORMAT (/,10x,'total charge for spin',i3,'=',f11.6,/,10x,         &
     &       'interst. charge =   ',f11.6,/,                            &
     &       (10x,'mt charge=          ',4f11.6,/))                     
 8020 FORMAT (/,10x,'total charge  =',f11.6,' pos. charge=',f11.6,      &
     &     ' difference=',f11.6)                                        
                                                                        
      IF (PRESENT(q_el)) q_el = qtot 
      IF (PRESENT(q_nuc)) q_nuc = q 
                                                                        
      INQUIRE(FILE ="gf_totalmix",EXIST= lexists) 
      IF (lexists) THEN 
         OPEN(99,FILE = "gf_totalmix",POSITION ="append") 
         WRITE(99,"(i5,2f30.20)") layer,qtot,q 
         CLOSE(99) 
      ENDIF 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
      END                                           
