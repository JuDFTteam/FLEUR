!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_kptsp1x1 
      use m_juDFT
      use m_juDFT 
      USE m_gf_types 
      USE m_constants 
      IMPLICIT NONE
      PRIVATE 
      INTEGER,PARAMETER :: maxngvec = 81 
      PUBLIC gf_kptsp1x1 
      CONTAINS 
      !<-- S:gf_kptsp1x1(cell,sym,kpts,gfinp)                           
                                                                        
      SUBROUTINE gf_kptsp1x1(cell,sym,kpts,gfinp) 
!-----------------------------------------------                        
!  Calculate the larger set of kpts needed in a                         
!  p(1x1) embedding potential calculation using the current             
!  kpts and the current information on the cell and symmetry.           
!  Additionally, a gf_setup.kpts.hdf file must be present for           
!  Symmetry and unit-cell information of the smaller system             
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_constants 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_cell),INTENT(IN)  :: cell 
      TYPE(t_sym),INTENT(IN)   :: sym 
      TYPE(t_kpts),INTENT(IN)  :: kpts 
      TYPE(t_embinp),INTENT(IN) :: gfinp 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: transmat(3,3) 
      REAL                :: gvec(3,maxngvec) 
      INTEGER             :: ngvec 
      TYPE(t_cell)        :: cell1 
      TYPE(t_sym)         :: sym1 
      REAL                :: kp(3) 
      REAL                :: kp_ALL(3,kpts%nkpt*sym%nop2*maxngvec) 
      INTEGER             :: kp_new_nk(kpts%nkpt*sym%nop2*maxngvec) 
      INTEGER             :: kp_ALL_nk(kpts%nkpt*sym%nop2*maxngvec) 
      INTEGER             :: n,nk,nk_all,nk_new,nk_new_all,nn 
      REAL                :: kp_all_new(3,kpts%nkpt*sym%nop2*maxngvec) 
      REAL                :: kp_new(3,kpts%nkpt*sym%nop2*maxngvec) 
      !>                                                                
                                                                        
      !<-- Read small setup                                             

      CALL priv_secondsetup(cell1,sym1) 

      !>                                                                
                                                                        
      !<-- Determine the transformation matrix and the set of g-vectors 
      !between the two setups                                           
      CALL priv_transsetup(gfinp%kpts,cell,cell1,transmat,gvec,ngvec) 
      !>                                                                
                                                                        
      !<-- Loop over all kpts and determine the corresponding kpts in   
      !small system                                                     
      !<-- Use symops to fill the BZ, check for duplicates              
      WRITE(*,*) "List of all possible old kpts:" 
      nk_all=0 
      DO nk = 1,kpts%nkpt 
         kp=kpts%bk(:,nk) 
!         kp = MATMUL(kpts%bk(:,nk),cell%amat)/2/pi_const               
         symloop:DO n = 1,sym%nop2 
            kp = MATMUL(REAL(sym%mrot(:,:,n)),kp) 
            !test if new kpts                                           
            DO nn = 1,nk_all 
               IF (ALL(kp_ALL(:,nn) == kp)) CYCLE symloop 
            ENDDO 
            nk_all = nk_all+1 
            kp_ALL(:,nk_all) = kp 
            kp_all_nk(nk_all) = nk 
            WRITE(*,"('n =',i5,' k = (',2(f10.5,1x),') irr.k =',i5)")   &
     &           nk_all,kp_ALL(:2,nk_all),kp_all_nk(nk_all)             
         ENDDO symloop 
      ENDDO 
      !>                                                                
      !<-- Use the translation g-vecs to make large kpoint-set          
      DO n=1,ngvec 
         DO nk=1,nk_all 
            kp_ALL(:,nk_all*(n-1)+nk) = kp_ALL(:,nk)+gvec(:,n) 
         ENDDO 
         kp_all_nk(nk_all*(n-1)+1:nk_all*n)=kp_all_nk(:nk_all) 
      ENDDO 
      nk_all= nk_all*ngvec 
      !>                                                                
                                                                        
      WRITE(*,*) "Transformation matrix" 
      WRITE(*,"(2(f10.5,1x))") transmat(1:2,1) 
      WRITE(*,"(2(f10.5,1x))")transmat(1:2,2) 
                                                                        
                                                                        
      !<-- Use the transformation matrix to transform to different      
      !unit-cell                                                        
      DO n=1,nk_all 
         kp_ALL(:,n)=matmul(transmat,kp_ALL(:,n)) 
      ENDDO 
      !>                                                                
                                                                        
      !<-- find new kpts                                                
      nk_new_all = 0 
      nk_new=0 
      nkloop:DO nk=1,nk_all 
         DO n=1,nk_new_all 
            IF (ALL(ABS(kp_ALL(:,nk)-kp_all_new(:,n))<1E-3)) THEN 
               WRITE(*,"('same:',4(f10.5,1x))") kp_ALL(1:2,nk)          &
     &              ,kp_all_new(1:2,n)                                  
               CYCLE nkloop 
            ENDIF 
            IF (ANY(ABS(kp_ALL(:,nk))>0.5)) THEN 
               WRITE(*,"('large:',4(f10.5,1x))") kp_ALL(1:2,nk) 
               CYCLE nkloop 
            ENDIF 
         ENDDO 
         !This is a new kpoint                                          
         nk_new=nk_new+1 
         kp_new(:,nk_new) = kp_ALL(:,nk) 
         kp_new_nk(nk_new) = kp_all_nk(nk) 
         WRITE(*,"('new:',i5,' k = (',2(f10.3,1x),')')") nk_new         &
     &        ,kp_ALL(1:2,nk)                                           
         !<-new kpoint and all symmetry equivalent                      
         DO n = 1,sym1%nop2 
            nk_new_all = nk_new_all+1 
            kp_ALL_new(:,nk_new_all) = MATMUL(REAL(sym1%mrot(:,:,n))    &
     &           ,kp_new(:,nk_new))                                     
!            WRITE(*,"('equivalent: k = (',2(f10.3,1x),')')")           
!     $           ,kp_ALL_new(1:2,nk_new_all)                           
         ENDDO 
         !>                                                             
      ENDDO nkloop 
      !>                                                                
      !<--Now write output                                              

      OPEN(99,FILE="kpts.info") 
      WRITE(99,"(2i5)") nk_new,kpts%nkpt 
      DO n=1,nk_new 
         WRITE(98,"(2i5)") n,kp_new_nk(n) 
      ENDDO 
      CLOSE(98) 
                                                                        
      OPEN(98,FILE="kpts.p1x1") 
      WRITE(98,"(i5,f10.5)") nk_new,1.0 
      DO n = 1,nk_new 
         WRITE(98,"(4f10.5)") kp_new(1:2,n),-1.0 
      ENDDO 
      CLOSE(98) 
       CALL juDFT_end("kpts-generated")

      !>                                                                
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: priv_secondsetup(cell1,sym1)                              

      SUBROUTINE priv_secondsetup(cell1,sym1) 
!-----------------------------------------------                        
!  reads the symmetry and the cell info from the gf_setup.kpts.hdf file 
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_iodop 
      IMPLICIT NONE 
      !<-- Arguments                                                    
      TYPE(t_cell),INTENT(OUT) :: cell1 
      TYPE(t_sym),INTENT(OUT) :: sym1 
      !>                                                                
      !<-- Locals                                                       
      LOGICAL       :: lexist 
      INTEGER       :: jspins 
      REAL          :: r1,r2,r3,r4 
      TYPE(t_atoms) :: atoms 
      !>                                                                
      INQUIRE(FILE ="gf_setup.kpts.hdf",EXIST = lexist) 
      IF (.NOT.lexist)                                                  &
     &     CALL juDFT_error                                               &
     &     ("You MUST give gf_setup.kpts.hdf-file of p(1x1) setup")     
       CALL juDFT_error("changed the code here",calledby="gf_kptsp1x1.F90")
!      CALL gf_readatt("gf_setup.kpts",jspins,r1,r2,r3,r4,atoms,cell1,sy
                                                                        
      !This is a memory leak because atoms is not properly deallocated  
                                                                        
      END SUBROUTINE 

      !>                                                                
                                                                        
      !<-- S: priv_transsetup(name,cell,cell1,transmat,gvec,ngvec)      
      SUBROUTINE priv_transsetup(name,cell,cell1,transmat,gvec,ngvec) 
!-----------------------------------------------                        
!  depending on the name of the transformation sets up the gvectors and 
!  the transformation matrix                                            
!                                                                       
!  At present the following transformations are known                   
!  - pnxm (n,m<9) -> a p(nxm) supercell                                 
!                                                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      IMPLICIT NONE 
      !<--Arguments                                                     
      CHARACTER(LEN=4),INTENT(IN) ::name 
      TYPE(t_cell),INTENT(IN)     :: cell,cell1 
      REAL   ,INTENT(OUT)         :: transmat(3,3) 
      REAL   ,INTENT(OUT)         :: gvec(3,MAXNGVEC) 
      INTEGER,INTENT(OUT)         :: ngvec 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n,m,i,ii 
      CHARACTER           :: mode 
      !>                                                                
      READ(name,"(a1,i1,1x,i1)") mode,n,m 
                                                                        
      WRITE(*,"('Generating kpts for:',a,'(',i1,'x',i1,')')") mode,n,m 
      WRITE(*,*) "G-Vectors:" 
      gvec = 0.0 
      IF (mode =="p") THEN 
         ngvec=0 
         DO i = 0,n-1 
            DO ii = 0,m-1 
               ngvec=ngvec+1 
               gvec(1:2,ngvec) = (/i,ii/) 
               WRITE(*,"(i3,'(',2(i5,1x),')')") ngvec,i,ii 
            ENDDO 
         ENDDO 
         IF (ngvec<2) THEN 
            WRITE(*,*) "Unknown mode:",name 
           CALL juDFT_error("gf_kptsp1x1") 
         ENDIF 
         transmat=0.0 
         transmat(1,1)=1/real(n) 
         transmat(2,2)=1/real(m) 
         IF (ABS(cell%area-cell1%area*n*m)>1E-5) THEN 
            WRITE(*,*) "Incorrect unit cell sizes" 
            WRITE(*,*) "Large cell:",cell%area 
            WRITE(*,*) "Small cell:",cell1%area 
            WRITE(*,*) "Scaling factor:",n*m 
            CALL juDFT_error("gf_kptsp1x1") 
         ENDIF 
      ELSE 
         WRITE(*,*) "Unknown mode:",name 
         CALL juDFT_error("gf_kptsp1x1") 
      ENDIF 
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
      END                                           
