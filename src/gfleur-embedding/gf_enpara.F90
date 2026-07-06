!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_enpara 
!----------------------------------------------------------------       
!Module provides interface to FLEUR                                     
!USE m_enpara,ONLY:r_enpara,w_enpara                                    
!----------------------------------------------------------------       
      use m_juDFT 
      IMPLICIT NONE 
      CONTAINS 
      SUBROUTINE gf_loadenpara(jspins,layers,atoms,enpara) 
!*****************************************************************      
!     DESC:loads enpara file                                            
!     Daniel Wortmann, Wed Aug 28 10:40:06 2002                         
!*****************************************************************      
      USE m_gf_types 
      USE m_fleur_INTERFACE 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)         :: jspins 
      TYPE(t_layers),INTENT(IN)  :: layers 
      TYPE(t_atoms),INTENT(IN)   :: atoms(:) 
      TYPE(t_enpara),INTENT(OUT) :: enpara(:) 
      !>                                                                
      LOGICAL              :: l_exist 
      CHARACTER(LEN = 200) :: line 
      INTEGER              :: l,n,i,ierr,zatom,lo,layer 
      INTEGER              :: ello_in(4),el_in(4) 
      INTEGER              :: ello(4,112),el(4,112) 
                                                                        
                                                                        
      INQUIRE(FILE ="enpara_atoms",EXIST = l_exist) 
      IF (.NOT.l_exist) THEN 
         CALL fleur_renpara(layers,atoms,jspins,enpara) 
      ELSE 
         !<-- read all lines in enpara_atoms                            
         OPEN(40,FILE="enpara_atoms",STATUS ="old") 
         DO 
            READ(40,"(a)",END = 99) line 
            IF (line(1:1) =="#") CYCLE 
            READ(line,*,ERR = 999) zatom,el_in,ello_in 
            el(:,zatom) = el_in 
            ello(:,zatom) = ello_in 
                                                                        
         ENDDO 
   99    CLOSE(40) 
         !>                                                             
         !<-- Allocate memory                                           
         DO layer = 1,layers%num_layers 
            ALLOCATE( enpara(layer)%lchange( 0:MAXVAL(atoms(layer)%lmax0&
     &           ),atoms(layer)%ntype,jspins)  )                        
            enpara(layer)%lchange = .FALSE. 
            ALLOCATE( enpara(layer)%llochg ( atoms(layer)%nlod,         &
     &           atoms(layer)%ntype,jspins)  )                          
            enpara(layer)%llochg = .FALSE. 
            ALLOCATE( enpara(layer)%el( 0:MAXVAL(atoms(layer)%lmax0),   &
     &           atoms(layer)%ntype,jspins)  )                          
            enpara(layer)%el =-999. 
            ALLOCATE( enpara(layer)%ello( atoms(layer)%nlod,            &
     &           atoms(layer)%ntype,jspins)  )                          
            enpara(layer)%ello =-999. 
            ALLOCATE( enpara(layer)%skiplo( atoms(layer)%ntype,jspins)  &
     &           )                                                      
            enpara(layer)%skiplo = 0 
         ENDDO 
         !>                                                             
         !<--assign enparas to atoms                                    
         WRITE(6,*) 
         WRITE(6,*) "Energy parameters from enpara_atoms" 
         WRITE(6,*) "Layer Atom   s  p  d  f  lo's" 
         DO l = 1,layers%num_layers 
            DO n=1,atoms(l)%ntype 
               enpara(l)%el(0:3,n,1)=el(:,nint(atoms(l)%zatom(n))) 
               enpara(l)%el(0:3,n,jspins)=el(:,nint(atoms(l)%zatom(n))) 
               enpara(l)%el(4:,n,:)=el(4,nint(atoms(l)%zatom(n))) 
               DO lo=1,atoms(l)%nlo(n) 
                  enpara(l)%ello(lo,n,:) = ello(atoms(l)%llo(lo,n)+1    &
     &                 ,nint(atoms(l)%zatom(n)))                        
               ENDDO 
               WRITE (6,"(2(i5,1x),4(i2,1x))")  l,n, (enpara(l)%el(i,n  &
     &              ,1),i = 0,3)                                        
               IF (atoms(l)%nlo(n) >= 1) THEN 
                  WRITE (6,"(20x,4(i2,1x))") (enpara(l                  &
     &                 )%ello(lo,n,1),lo = 1,atoms(l)%nlo(n))           
               ENDIF 
            ENDDO 
         ENDDO 
         !>                                                             
      ENDIF 
      RETURN 
  999 WRITE(*,*) "Error reading enpara_atoms:" 
      WRITE(*,*) line 
      CALL juDFT_error("Problem in enpara_atoms") 
      END SUBROUTINE 
      SUBROUTINE fleur_writeenpara(atoms,jspins,enpara) 
!*****************************************************************      
!     DESC:loads enpara file                                            
!     Daniel Wortmann, Wed Aug 28 10:40:06 2002                         
!*****************************************************************      
      USE m_gf_types 
      USE m_enpara,ONLY:w_enpara 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,INTENT(IN)         ::jspins 
      TYPE(t_atoms),INTENT(IN)   ::atoms 
      TYPE(t_enpara),INTENT(IN)  ::enpara 
                                                                        
      LOGICAL,PARAMETER::film=.TRUE. 
      INTEGER,PARAMETER::nw=1 
      INTEGER::jsp 
                                                                        
      OPEN(40,FILE="enpara") 
      !read the enpara-file                                             
      DO jsp=1,jspins 
         CALL w_enpara(                                                 &
     &        maxval(atoms%lmax),atoms%nlod,atoms%ntype,nw              &
     &        ,jsp,film,atoms%nlo,enpara%skiplo(:                       &
     &        ,jsp),enpara%ello(:,:,jsp),enpara%el(:,:,jsp)             &
     &        ,enpara%evac(:,jsp),enpara%lchange(:,:,jsp)               &
     &        ,enpara%llochg(:,:,jsp),enpara%lchg_v(1,jsp)              &
     &        ,enpara%enmix(jsp),16)                                    
      ENDDO 
      CLOSE(40) 
      END SUBROUTINE 
      END                                           
