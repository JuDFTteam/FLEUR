!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_kpts 
      USE m_constants, ONLY: oUnit
          IMPLICIT NONE
      CONTAINS 
      SUBROUTINE gf_loadkpts(cell,sym,l_pe0,kpts) 
!*****************************************************************      
! DESC:loads the kpts and enpara file                                   
!                          Daniel Wortmann, Wed Aug 28 10:40:06 2002    
!*****************************************************************      
      use m_juDFT 
      USE m_gf_types 
      IMPLICIT NONE 
!     Arguments                                                         
      LOGICAL,INTENT(IN)         ::l_pe0 
      TYPE(t_cell),INTENT(IN)    ::cell 
      TYPE(t_sym),INTENT(IN)     ::sym 
      TYPE(t_kpts),INTENT(OUT)   ::kpts 
      INTEGER                    :: i,nk 
      REAL                       :: fscale 
      CHARACTER                  :: testchar 
      !load kpts-file                                                   
      !Only Film-version is supported                                   
      !yuu feature not supported!                                       
      OPEN (41,FILE='kpts',FORM='formatted',STATUS='old') 
      READ (41,*) testchar 
      REWIND(41) 
      IF (testchar =="&") THEN 
         CALL priv_transformkpts(kpts,sym) 
      ELSE 
         READ (41,"(i5,f20.10)") kpts%nkpt,fscale 
         IF (fscale==0.0) fscale = 1.0 
         ALLOCATE(kpts%bk(3,kpts%nkpt)                                 &
     &        ,kpts%wtkpt(kpts%nkpt))                                 
         kpts%bk(3,:) = 0.0 
         DO nk = 1,kpts%nkpt 
            READ (41,"(4f10.5)") (kpts%bk(i,nk),i=1,2),kpts%wtkpt(nk) 
!         kpts%bk(:,nk)=matmul(cell%amat,kpts%bk(:,nk))/tpi/fscale      
            kpts%bk(:,nk) = kpts%bk(:,nk)/fscale 
         ENDDO 
      ENDIF 
      CLOSE (41) 
      WRITE(*,*) "Read ",kpts%nkpt," k-points" 
                                                                        
       !<-- output symmetry equivalent-kpoint-list,set weights          
       CALL priv_symkpt(sym,cell,kpts,l_pe0) 
       !>                                                               
       IF (ABS(SUM(kpts%wtkpt)-1.0)>0.01) THEN 
          WRITE(*,*) "Error: Total sum of kpts-weight incorrect" 
          WRITE(*,*) "Sum:",sum(kpts%wtkpt),' should be 1.0' 
          WRITE(*,*) "Set one weight to negative value to assign the"   &
     &         ," weights according to symmetry"                        
          CALL juDFT_error("gf_kpts") 
       ENDIF 
       RETURN 
                                                                        
      END SUBROUTINE 
                                                                        
      !<-- S: priv_symkpt(sym,cell,kpts)                                
      SUBROUTINE priv_symkpt(sym,cell,kpts,lpe0) 
!******************************************                             
!                                                                       
!                          D. Wortmann                                  
!******************************************                             
      use m_juDFT 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_sym),INTENT(IN)      :: sym 
      TYPE(t_cell),INTENT(IN)     :: cell 
      TYPE(t_kpts),INTENT(INOUT)  :: kpts 
      LOGICAL,INTENT(IN)          :: lpe0 
      !>                                                                
      !<--Locals                                                        
                                                                        
      INTEGER :: nk,n,i,nk1 
      INTEGER :: nkpt(kpts%nkpt) 
      LOGICAL :: new 
      REAL    :: bk(2),bk2(2,sym%nop,kpts%nkpt) 
                                                                        
      !>                                                                
                                                                        
                                                                        
      DO nk=1,kpts%nkpt 
         nkpt(nk)=0 
         DO n=1,sym%nop 
            bk=MATMUL(kpts%bk(1:2,nk),1.0*sym%mrot(1:2,1:2,n)) 
            new=.TRUE. 
                        !extra loop over kpts to check for double kpts  
            DO nk1=1,nk 
               DO i=1,nkpt(nk1) 
                  IF (ALL(bk(:)==bk2(:,i,nk1))) new=.FALSE. 
               ENDDO 
            ENDDO 
            IF (new) THEN 
               nkpt(nk)=nkpt(nk)+1 
               bk2(:,nkpt(nk),nk)=bk(:) 
            ENDIF 
         ENDDO 
      ENDDO 
      IF (ANY(kpts%wtkpt<0))  kpts%wtkpt=REAL(nkpt)/SUM(nkpt) 
      !if on pe0 write to output-file                                   
      IF(lpe0) THEN 
         WRITE(oUnit,8120) kpts%nkpt 
         WRITE(oUnit,8121) SUM(nkpt) 
 8120    FORMAT (1x,/,' number of k-points  =',i5) 
 8121    FORMAT (1x,'         ->in total BZ:',i6,/,'   No  Sym ',5x,    &
     &        'coordinates',t37,'weights')                              
         DO nk=1,kpts%nkpt 
            DO i=1,nkpt(nk) 
               bk(1:2) = MATMUL(bk2(1:2,i,nk),cell%bmat(1:2,1:2)) 
               WRITE(oUnit,'(i5,1x,i4,1x,3(f20.15,1x))')nk,i,bk(1),bk(2)    &
     &              ,kpts%wtkpt(nk)                                    
            ENDDO 
         ENDDO 
      ENDIF 
      !<--check that no kpts occurs double                              
      IF (ANY(nkpt==0)) THEN 
         WRITE(*,*) "Double k-point" 
      ENDIF 
      !>                                                                
                                                                        
                                                                        
      END SUBROUTINE 
      !>                                                                
                                                                        
      !<-- S: priv_transformkpts(kpts,sym)                              
      SUBROUTINE priv_transformkpts(kpts,sym) 
!-----------------------------------------------                        
!     transforms the kpts according to a transformation                 
!     specified at the beginning at the kpts-file                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_gf_math 
      IMPLICIT NONE 
      !<--Arguments                                                     
      TYPE(t_sym),INTENT(IN)     :: sym 
      TYPE(t_kpts),INTENT(INOUT) :: kpts 
      !>                                                                
      !<-- Locals                                                       
      REAL                :: trans(2,2) 
      REAL   ,ALLOCATABLE :: bk(:,:) 
      LOGICAL,ALLOCATABLE :: valid(:) 
      INTEGER             :: nx,ny 
      LOGICAL             :: invers 
      INTEGER             :: x,y,nkpts,nk,n,nn,nk1,nk_max 
      REAL                :: fscale,bk_test(2) 
                                                                        
      NAMELIST /rescale/nx,ny,trans,invers 
                                                                        
      !>                                                                
      READ (41,rescale) 
      IF (invers) trans=mat_inverse(trans) 
      READ (41,"(i5,f20.10)") nkpts,fscale 
      IF (fscale   == 0.0) fscale = 1.0 
                                                                        
      ALLOCATE(bk(2,nkpts*(2*nx+1)*(2*ny+1))) 
      ALLOCATE(valid(nkpts*(2*nx+1)*(2*ny+1))) 
      n = 1 
      DO nk = 1,nkpts 
         READ (41,"(4f10.5)") bk(:,n) 
         bk(:,n)   = bk(:,n)/fscale 
         nn = 1 
         DO x =-nx,nx 
            IF (x == 0) CYCLE 
            DO y =-ny,ny 
               IF (y == 0) CYCLE 
               bk(:,n+nn) = bk(:,n)+(/x,y/) 
               nn = nn+1 
            ENDDO 
         ENDDO 
         n = n+nn 
      ENDDO 
      nk_max = n-1 
      !OK now we have all possible k-point now transform them to new    
      !coordinates                                                      
      bk(:,:) = MATMUL(trans,bk) 
      !Move all points into first BZ                                    
      WHERE(ABS(bk)>0.5) 
         bk = bk-NINT(bk) 
      END WHERE 
      !Check for double k-points                                        
      valid = .TRUE. 
      WRITE(*,*) "Max no of k-points:",nk_max 
      DO nk = 1,nk_max 
         DO n = 1,sym%nop 
            bk_test = MATMUL(bk(1:2,nk),1.0*sym%mrot(1:2,1:2,n)) 
                             !extra loop over kpts to check for double k
            DO nk1  = 1,nk-1 
               IF (ALL(ABS(bk_test-bk(:,nk1))<1E-8)) valid(nk) = .FALSE. 
            ENDDO 
         ENDDO 
      ENDDO 
      !Now put kpts into type                                           
      kpts%nkpt=COUNT(valid(:nk_max)) 
      ALLOCATE(kpts%bk(3,kpts%nkpt),kpts%wtkpt(kpts%nkpt)) 
      nn = 1 
      DO nk = 1,nk_max 
         IF (valid(nk)) THEN 
            kpts%bk(:2,nn) = bk(:,nk) 
            kpts%bk(3,nn) = 0.0 
            kpts%wtkpt(nn) =-1. 
            nn=nn+1 
         ENDIF 
      ENDDO 
      DEALLOCATE(bk,valid) 
      END SUBROUTINE 
      !>                                                                
      END                                           
