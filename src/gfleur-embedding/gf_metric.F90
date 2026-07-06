!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_metric 
      IMPLICIT NONE
      integer,parameter:: nmz=250,nmzxy=100
      real,parameter   :: delz=0.1
      private
      public gf_apply_metric
!-----------------------------------------------                        
! DESC:The charge/potential metric used for mixing                      
!                 Daniel Wortmann, (08-01-03)                           
!-----------------------------------------------                        
      CONTAINS 
                                                                        
                                                                        
      !<-- F: gf_apply_metric(stars,atoms,sphhar,jspins,num_layers,a)   
      FUNCTION gf_apply_metric(l_surface,l_potmix,mpi,stars,atoms,cell,sphhar     &
     &     ,sym,jspins,num_layers,a)RESULT(b)                               
!-----------------------------------------------                        
!                                                                       
!             (last modified: 2004-00-00) D. Wortmann                   
!-----------------------------------------------                        
      USE m_Gf_Types 
      use m_metrz0
      use m_gf_stepsanaly, only: gf_gspaceconvolve, gf_initstepsanaly
      IMPLICIT NONE 
      !<-- Arguments                                                    
      LOGICAL,INTENT(IN)         :: l_potmix,l_surface
      TYPE(t_gfmpi),INTENT(IN)     :: mpi 
      TYPE(t_stars),INTENT(IN)   :: stars(:) 
      TYPE(t_atoms),INTENT(IN)   :: atoms(:) 
      TYPE(t_sphhar),INTENT(IN)  :: sphhar(:) 
      TYPE(t_cell),INTENT(IN)    :: cell(:)
      TYPE(t_sym),INTENT(IN)     :: sym(:)
      INTEGER,INTENT(IN)         :: jspins,num_layers 
      REAL   ,INTENT(IN)         :: a(:) 
      REAL                       :: b(size(a)) 
                                                                        
      !>                                                                
      !<-- Locals                                                       
      COMPLEX,ALLOCATABLE :: vpw(:),vpw_w(:) 
      INTEGER             :: imap,layer,js,na,n,l,j 
      REAL                :: dxn,dxn2,dxn4,dvol
      real                :: wght(nmz)
      !>                                                                
                                                                        
      ALLOCATE(vpw(MAXVAL(stars%ng3))) 
      ALLOCATE(vpw_w(MAXVAL(stars%ng3))) 
      !<-- apply the metric                                             
      imap = 0 
      DO layer = 1,num_layers 
         DO js = 1,jspins 
            !interstitial                                               
            vpw(:stars(layer)%ng3) = CMPLX(a(imap+1:imap+stars(layer    &
     &           )%ng3),a(imap+1+stars(layer)%ng3:imap+2*stars(layer    &
     &           )%ng3))                                                
            IF (stars(layer)%ng3<SIZE(vpw)) vpw(stars(layer)%ng3+1:) = 0.0
            CALL gf_initstepsanaly(stars(layer),0)
            CALL gf_gspaceconvolve(layer,stars(layer), &
     &           .0,vpw(:),vpw_w(:))
            b(imap+1:imap+stars(layer)%ng3) = REAL(vpw_w(:stars(layer   &
     &           )%ng3))*cell(layer)%omtil                              

            imap = imap+stars(layer)%ng3 
            b(imap+1:imap+stars(layer)%ng3) = AIMAG(vpw_w(:stars(layer  &
     &           )%ng3))*cell(layer)%omtil

            imap = imap+stars(layer)%ng3 
            !muffin-tin part                                            
            na = 1 
            DO  n = 1,atoms(layer)%ntype 
               dxn = atoms(layer)%neq(n)*atoms(layer)%dx(n)/3.0E0 
               dxn2 = 2.0E0 *dxn 
               dxn4 = 4.0E0 *dxn 
               DO l = 0,sphhar(layer)%nlh(sym(layer)%ntypsy(na)) 
                  imap = imap + 1 
                  b(imap) = a(imap)*dxn/atoms(layer)%rmsh(1,n) 
                                           !charge density mixing       

                  IF (.NOT.l_potmix) THEN 
                     DO j = 2,atoms(layer)%jri(n)-1,2 
                        imap = imap + 2 
                        b(imap-1) = a(imap-1)*  dxn4/atoms(layer        &
     &                       )%rmsh(j,n)                                
                        b(imap) = a(imap)*  dxn2/atoms(layer)%rmsh(j+1  &
     &                       ,n)                                        
                     ENDDO 
!     CHANGE JR 96/12/01                                                
!     take care when jri(n) is even                                     
                     imap = imap+1-MOD(atoms(layer)%jri(n),2) 
                     b(imap) = a(imap)*dxn/atoms(layer)%rmsh(atoms(layer&
     &                    )%jri(n),n)

                  ELSE 
!                                                                       
!     for the potential multiply by r^4                                 
!                                                                       
                     DO j = 2,atoms(layer)%jri(n)-1,2 
                        imap = imap + 2 
                        b(imap-1) = a(imap-1)*dxn4*atoms(layer)%rmsh(j,n&
     &                       )**3                                       
                        b(imap)   = a(imap)*dxn2*atoms(layer)%rmsh(j+1,n&
     &                       )**3                                       
                     ENDDO 
                     imap = imap+1-MOD(atoms(layer)%jri(n),2) 
                     b(imap) = a(imap)*dxn*atoms(layer)%rmsh(atoms(layer&
     &                    )%jri(n),n)**3                                
                  ENDIF 
               ENDDO 
               na = na + atoms(layer)%neq(n) 
            ENDDO 
               ! jspins                                                 
         ENDDO 
               ! layers                                                 
      ENDDO 
      !>

      !surface part
      if (l_surface) then
         dvol = cell(1)%area*delz
         DO js=1,jspins
             call metr_z0(nmz,wght)
             b(imap+1:imap+nmz)=a(imap+1:imap+nmz)*wght(:nmz)*dvol
             imap=imap+nmz
             call metr_z0(nmzxy,wght)
             DO n=1,stars(1)%ng2-1
                b(imap+1:imap+nmzxy)=a(imap+1:imap+nmzxy)*wght(:nmzxy)*dvol*stars(1)%nstr2(n+1)
                imap=imap+nmzxy
                b(imap+1:imap+nmzxy)=a(imap+1:imap+nmzxy)*wght(:nmzxy)*dvol*stars(1)%nstr2(n+1)
                imap=imap+nmzxy
             enddo
        enddo
      endif

      IF (imap+1<SIZE(b)) b(imap+1:) = 0.0 
      DEALLOCATE(vpw,vpw_w) 
      END FUNCTION 
      !>                                                                
      END                                           
