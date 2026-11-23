!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_inp_describe 
      USE m_constants, ONLY: oUnit
      use m_juDFT 
      IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_inp_describe(layers,jspins,atoms,cell,gfinp) 
!*****************************************************************      
! DESC:Writes description of the setup to out                           
!                          Daniel Wortmann, Wed May  7 10:09:34 2003    
!*****************************************************************      
      USE m_gf_types 
      USE m_constants, ONLY: namat_const, oUnit 
      IMPLICIT NONE 
!     Arguments                                                         
      INTEGER,INTENT(IN)       ::jspins 
      TYPE(t_layers),INTENT(IN)::layers 
      TYPE(t_embinp),INTENT(IN) ::gfinp 
      TYPE(t_atoms),INTENT(IN) ::atoms(:) 
      TYPE(t_cell),INTENT(IN)  ::cell(:) 
                                                                        
      INTEGER:: n,na,ieq,ilo,i,j,layer 
                                                                        
                                                                        
      IF (jspins>1) WRITE(oUnit,'("Spin polarized calculation",/,/)') 
                                                                        
      IF (.NOT.gfinp%l_gmat) THEN 
         WRITE(oUnit,'(/,a30,/,/)')                                         &
     &     'WARNING! Green-function not calculated'                     
      ELSE 
         IF (gfinp%l_addemb) THEN 
            WRITE(oUnit,*) 'Embedded Green function calculation' 
            WRITE(oUnit,*) 
         ELSE 
            WRITE(oUnit,*) 'Green function with von-Neumann boundary ',     &
     &           'conditions will be calculated'                        
            WRITE(oUnit,*) 
         ENDIF 
      ENDIF 
      IF (gfinp%l_tmat) WRITE(oUnit,*) 'T-Matrix calculation' 
      IF (gfinp%l_CBS) THEN 
         WRITE(oUnit,*) 'CBS-calculation with the following parameters:' 
         WRITE(oUnit,*) 'eps_current:',gfinp%eps_current 
         WRITE(oUnit,*) 'eps_non_bloch:',gfinp%eps_non_bloch 
         WRITE(oUnit,*) 'CBS_bz:',gfinp%CBS_bz 
         WRITE(oUnit,*) 'dp:',gfinp%dp1,gfinp%dp2 
      ENDIF 
      WRITE(oUnit,*) 
      IF (gfinp%curr>0) WRITE(oUnit,*) 'Current will be calculated' 
      WRITE(oUnit,*) 
      WRITE(oUnit,*) 'Setup of Basis function:' 
      IF (gfinp%l_charge.OR.(gfinp%napw(1)<0)) THEN 
         WRITE(oUnit,*) 'SC-setup! Basis will be from spherical g-vectors' 
      ELSE 
         WRITE(oUnit,*)                                                     &
     &    'Basis functions will be from zylindrical set of g-vectors'   
         DO layer=1,layers%num_layers 
            WRITE(oUnit,*) "layer=",layer 
           WRITE(oUnit,*)                                                   &
     &        'z-Dependend basis functions in region I:',               &
     &              2*gfinp%napw(layer)+1                               
           IF (gfinp%npw==0) WRITE(oUnit,*)                                 &
     &        'Region II is treated analytically'                       
           IF (gfinp%npw<0) WRITE(oUnit,*)                                  &
     &        'Region II is treated analytically with shifted potential'
           IF (gfinp%npw>0) WRITE(oUnit,*)                                  &
     &        'z-Dependend basis functions in region II:',              &
     &          2*gfinp%npw+1                                           
         ENDDO 
      ENDIF 
      IF(layers%num_layers>1) WRITE(oUnit,*) 'Calculation with ',        &
     &                        layers%num_layers,' Layer'                
      DO layer=1,layers%num_layers 
        IF(layers%num_layers>1)WRITE(oUnit,*)'Layer=',layer 
        !The Cell                                                       
        WRITE(oUnit,'(/,/,"Unit Cell defined by:")') 
        WRITE (oUnit,FMT=8080) 
 8080   FORMAT (1x,'bravais matrices of real and reciprocal lattices',/) 
        DO 60 i = 1,3 
           WRITE (oUnit,8090)(cell(layer)%amat(i,j),j=1,3),                 &
     &                    (cell(layer)%bmat(i,j),j=1,3)                 
   60   CONTINUE 
 8090   FORMAT (3x,3f10.6,3x,3f10.6) 
        WRITE (oUnit,FMT=8100) cell(layer)%omtil,cell(layer)%vol,           &
     &                    cell(layer)%area                              
 8100   FORMAT (/,4x,'the volume of the unit cell omega-tilda=',f12.6,/,&
     &       10x,'the volume of the unit cell omega=',f12.6,/,2x,       &
     &       'the area of the two-dimensional unit cell=',f12.6)        
        !Now the atoms                                                  
        WRITE(oUnit,*) 
        na=0 
        WRITE(oUnit,*) 'No of atoms:',atoms(layer)%nat 
        WRITE(oUnit,*) 'No of types:',atoms(layer)%ntype 
        WRITE(oUnit,*) 'Atoms:' 
        DO n=1,atoms(layer)%ntype 
          WRITE(oUnit,*) '******TYPE:',n 
          WRITE(oUnit,'(a3,i3,3i5,2f10.6)') namat_const(atoms(layer)%nz(n)),      &
     &         atoms(layer)%nz(n),                                      &
     &         atoms(layer)%econf(n)%num_core_states,atoms(layer)%lmax(n),               &
     &         atoms(layer)%jri(n),atoms(layer)%rmt(n),                 &
     &         atoms(layer)%dx(n)                                       
          IF (atoms(layer)%lda_u(n)%l>=0) THEN 
            WRITE (oUnit,8180) atoms(layer)%lda_u(n)%l,                     &
     &            atoms(layer)%lda_u(n)%u,atoms(layer)%lda_u(n)%j       
 8180       FORMAT ('&ldaU l=',i1,',u=',f3.1,',j=',f3.1,' /') 
          ELSE 
            WRITE (oUnit,*) 'no LDAU+U   ' 
          ENDIF 
          IF (atoms(layer)%nlo(n)>0) THEN 
            WRITE (oUnit,9090) atoms(layer)%nlo(n),(atoms(layer)%llo(ilo,n),&
     &            ilo=1,atoms(layer)%nlo(n))                            
 9090       FORMAT ('Local Orbitals:nlo=',i2,',llo=',20i2) 
          ELSE 
            WRITE (oUnit,*) 'No local Orbitals' 
          ENDIF 
          !Loop over equivalent atoms                                   
          DO ieq=1,atoms(layer)%neq(n) 
            na = na + 1 
            WRITE (oUnit,9100) ieq,(atoms(layer)%taual(i,na),i=1,3) 
 9100       FORMAT ('Pos:',i2,'->',4f10.6) 
          ENDDO 
        ENDDO 
        IF (na/=atoms(layer)%nat)                                     &
     &     CALL juDFT_error('nat.ne.atoms(layer)%nat?')                   
            !layer                                                      
      ENDDO 
                                                                        
      END SUBROUTINE 
      END                                           
