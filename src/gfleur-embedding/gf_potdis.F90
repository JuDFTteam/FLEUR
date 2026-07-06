!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_potdis 
          IMPLICIT NONE
          private
          public gf_potdis
      CONTAINS 
      SUBROUTINE gf_potdis(                                             &
     &                 jspins,atoms,stars,sphhar,mpi,sym,               &
     &                 cell,c_layer,l_pot,layer,distance,volume,enpara)
!*****************************************************************      
!     DESC:READ's in the old and new potential/charge-density and       
!     calculates the distances in MT+INT                                
!                                                                       
!     Daniel Wortmann, Fri Aug 23 15:18:32 2002                         
!*****************************************************************      
      USE m_gf_types 
      USE m_constants 
      USE m_gf_fft 
      USE m_gf_iodop 
      USE m_gf_stepsanaly 
      USE m_intgr, ONLY : intgr3 
      USE m_gf_plot 
      IMPLICIT NONE 
!     Arguments                                                         
      TYPE(t_atoms),INTENT(IN)  :: atoms 
      TYPE(t_stars),INTENT(IN)  :: stars 
      TYPE(t_cell),INTENT(IN)   :: cell
      !height of the physical region of the layer (was cell%c)
      REAL,INTENT(IN)           :: c_layer
      TYPE(t_sphhar),INTENT(IN) :: sphhar 
      TYPE(t_mpi),INTENT(IN)    :: mpi 
      TYPE(t_sym),INTENT(IN)    :: sym 
      INTEGER,INTENT(IN)        :: jspins 
      LOGICAL,INTENT(IN)        :: l_pot 
      INTEGER,INTENT(IN)        :: layer 
      REAL,INTENT(INOUT),OPTIONAL :: distance,volume 
      type(t_enpara),intent(in),optional :: enpara
                                                                        
!     Locals                                                            
      INTEGER :: js,na,n,lh,jspin 
      REAL    :: disint(3)
      REAL    :: dis(atoms%ntype,3) 
      REAL,ALLOCATABLE :: vr_old(:,:,:,:),vr_new(:,:,:,:),vr_dif(:,:,:,:&
     &     )                                                            
      REAL,ALLOCATABLE :: rh(:) 
      REAL             :: rh_s,totaldist 
      COMPLEX          :: vpw_old(stars%ng3,3) 
      COMPLEX          :: vpw_new(stars%ng3,3) 
      COMPLEX          :: vpw_dif(stars%ng3,3) 
      REAL             :: vpw_r(0:27*stars%mx1*stars%mx2*stars%mx3-1) 
      INTEGER          :: fileid 
      INTEGER          :: j,pl_mode 
                                                                        
                                                                        
                                                                        

                                                                        
      ALLOCATE(vr_old(MAXVAL(atoms%jri),0:MAXVAL(sphhar%nlh),atoms%ntype&
     &     ,jspins),vr_new(MAXVAL(atoms%jri),0:MAXVAL(sphhar%nlh)       &
     &     ,atoms%ntype,jspins),vr_dif(MAXVAL(atoms%jri)                &
     &     ,0:MAXVAL(sphhar%nlh),atoms%ntype,3),rh(MAXVAL(atoms%jri     &
     &     )))                                                          
                                                                        
      !<-- read the charges/potential                                   
                                                                        
      totaldist=0.0 
                                                                        
                                                                        
      IF (l_pot) THEN 
         fileid =GF_POTFILE 
         pl_mode = GF_PLOT_DIFFPOT 
      ELSE 
         fileid =GF_CDNFILE 
         pl_mode = GF_PLOT_DIFFRHO 
      ENDIF 
                                                                        
      CALL gf_loddop(fileid,layer,jspins,                               &
     &     atoms,stars,sphhar,                                          &
     &     vr_old,vpw_old,old=.true.)
      CALL gf_loddop(fileid,layer,jspins,                               &
     &     atoms,stars,sphhar,                                          &
     &     vr_new,vpw_new)
                                                                        
                                                                        
      !Differences for first spin                                       
      vr_dif(:,:,:,1) = vr_new(:,:,:,1)-vr_old(:,:,:,1) 
      vpw_dif(:,1) = vpw_new(:,1)-vpw_old(:,1) 
      IF (jspins>1) THEN 
         !Second spin                                                   
         vr_dif(:,:,:,2) = vr_new(:,:,:,2)-vr_old(:,:,:,2) 
         vpw_dif(:,2) = vpw_new(:,2)-vpw_old(:,2) 
         ! 'spin' potential difference                                  
         vr_dif(:,:,:,3) = vr_new(:,:,:,1)-vr_new(:,:,:,2)-vr_old(:,:,: &
     &        ,1)+vr_old(:,:,:,2)                                       
         vpw_dif(:,3) = vpw_new(:,1)-vpw_new(:,2)-vpw_old(:,1)          &
     &        +vpw_old(:,2)                                             
      END IF 

      if (.not.l_pot.and.present(enpara)) then
      ! Now create the l-dist
        CALL gf_loddop(GF_POTFILE,layer,jspins,                               &
                       atoms,stars,sphhar,                                          &
                          vr_new,vpw_new)
        call gf_project_differences(layer,atoms,enpara,vr_new,vr_dif(:,:,:,:min(jspins,2)))
      endif
                                                                        
      CALL gf_plot(layer,stars,cell,atoms,sym,jspins,vpw_dif            &
     &     ,pl_mode,sphhar,vr_dif)                                      
                                                                        
      !>                                                                
      !<-- init the stepfunctions                                       
      CALL gf_initstepsanaly(stars,0) 
                                                                        
      !>                                                                
                                                                        
      jspin = 1 
      IF (jspins>1) jspin = 3 
      DO  js              = 1,jspin 
         !<-- m.t. part                                                 
                                                                        
         na   = 1 
         DO  n = 1,atoms%ntype 
            DO j=1,atoms%jri(n) 
             IF(.NOT.l_pot)THEN 
               rh(j) = vr_dif(j,0,n,js)*vr_dif(j,0,n,js)/               &
     &                              (atoms%rmsh(j,n))**4
             ELSE
             	rh(j) =  vr_dif(j,0,n,js)*vr_dif(j,0,n,js)*fpi_const
!               rh(j) = vr_dif(j,0,n,js)*vr_dif(j,0,n,js)*               &
!     &                              (atoms%rmsh(j,n))**3*sqrt(fpi)
             ENDIF 
            ENDDO 
            DO lh = 1,sphhar%nlh(sym%ntypsy(na)) 
              IF(.NOT.l_pot)THEN 
                DO j=1,atoms%jri(n) 
                  rh(j) = rh(j) + vr_dif(j,lh,n,js)*vr_dif(j,lh,n,js)/  &
     &                               (atoms%rmsh(j,n))**4
                ENDDO 
              ELSE 
                DO j=1,atoms%jri(n) 
                  rh(j) = rh(j) + vr_dif(j,lh,n,js)*vr_dif(j,lh,n,js)*  &
     &                               (atoms%rmsh(j,n))**2               
                ENDDO 
              ENDIF 
            ENDDO 
            CALL intgr3(rh,atoms%rmsh(:,n),atoms%dx(n),atoms%jri(n&
     &           ),rh_s)                                                
            dis(n,js) = sqrt(atoms%neq(n)*rh_s/atoms%volmts(n)) 
            totaldist=totaldist+atoms%neq(n)*rh_s 
            na = na + atoms%neq(n) 
         ENDDO 
                                                                        
         !>                                                             
         !<-- interstitial part                                         
                                                                        
         CALL gf_fft3d(vpw_r,vpw_dif(:,js),stars,GF_FFT_TO_REAL_SPACE) 
         vpw_r = REAL(vpw_r)**2 
         CALL gf_fft3d(vpw_r,vpw_dif(:,js),stars,GF_FFT_TO_G_SPACE) 
                                                                        
         CALL gf_gspaceconvolve(layer,stars,0.0,       &
     &        vpw_dif(:,js),vpw_new(:,js))                              

         totaldist = totaldist+vpw_new(1,js)*cell%omtil

         disint(js) = SQRT(vpw_new(1,js)) 
                                                                        
         !>                                                             
      ENDDO 
      !<-- Output                                                       
                                                                        
      WRITE(oUnit,*) "Layer:",layer 
      IF (l_pot) THEN 
         WRITE(oUnit,'(a,a)')                                               &
     &        'Detailed Listing of Potential-Distances follows'         &
     &        ,'(in a.u.)'                                              
      ELSE 
         WRITE(oUnit,'(a,a)')                                               &
     &        'Detailed Listing of Charge-Density-Distances follows'    &
     &        ,'(in mBohr/a.u.)'                                        
      ENDIF 
      IF (jspins>1) THEN 
         WRITE(oUnit,'(a)')                                                 &
     &   'Atom:       First Spin      Second Spin     Magnetisation'    
      ELSE 
         WRITE(oUnit,'(a)')                                                 &
     &        'Atom:    Distance'                                       
      ENDIF 
      IF(l_pot) THEN 
         rh_s = 1.0 
      ELSE 
         rh_s = 1000.0 
      ENDIF 
      DO n = 1,atoms%ntype 
         WRITE(oUnit,8000) n,dis(n,1:jspin)*rh_s 
      ENDDO 
      WRITE(oUnit,8010) disint(1:jspin)*rh_s 
 8000 FORMAT ('----> ',i3,':',3f13.6) 
 8010 FORMAT ('----> INT:',3f13.6) 
                                                                        
                                                                        

                                                                        
      IF (PRESENT(distance)) distance   = distance+abs(totaldist)
      rh_s=cell%area*c_layer
      IF (PRESENT(volume)) volume = volume+rh_s
      IF (l_pot) THEN 
         WRITE(oUnit,*) layer," Layer Distance =",SQRT(abs(totaldist)/rh_s)
      ELSE 
         WRITE(oUnit,*) layer," Layer Distance =",1000*SQRT(abs(totaldist)/rh_s)
      ENDIF 
      !>                                                                
      DEALLOCATE(vr_old,vr_new,vr_dif,rh) 
      RETURN 
      END SUBROUTINE 

      subroutine gf_project_differences(layer,atoms,enpara,vr,drho)
      use m_gf_types
      use m_intgr, only: intgr3
        use m_radfun
      implicit none
      integer,intent(in)        :: layer
      type(t_atoms),intent(in)  :: atoms
      type(t_enpara),intent(in) :: enpara
      real,intent(in)           :: drho(:,0:,:,:)
      real,intent(in)           :: vr(:,:,:,:)


      integer   :: jspin,itype,l,jspins,noded,nodeu
      real      :: norm,norm1,us,dus,uds,duds,ddn,wronk
      TYPE(t_usdus) :: ud_loc
      real      :: l_proj(0:2)
      real,allocatable::f(:,:),ff(:,:),s(:),rho(:)


      jspins=size(drho,4)
      allocate(rho(size(drho,1)))
      allocate(s(size(drho,1)),ff(size(drho,1),2),f(size(drho,1),2))



      DO jspin=1,jspins
        DO itype=1,atoms%ntype
            rho=drho(:,0,itype,jspin)/atoms%rmsh(:,itype)**2
            CALL intgr3(abs(rho),atoms%rmsh(:,itype),atoms%dx(itype),atoms%jri(itype),norm1)
            DO l=0,2
                IF (.NOT.ALLOCATED(ud_loc%us)) CALL ud_loc%init(atoms,jspin)
                CALL radfun(l,itype,jspin,enpara%el0(l,itype,jspin),          &
                            vr(:,1,itype,jspin),atoms,f(:,:),ff(:,:),          &
                            ud_loc,nodeu,noded,wronk)
                 !construct the square of wavefunction
                 s=(f(:,1)*f(:,1)+f(:,2)*f(:,2))/atoms%rmsh(:,itype)**2
                 CALL intgr3(s,atoms%rmsh(:,itype),atoms%dx(itype),atoms%jri(itype),norm)
                 s=s/norm
                 CALL intgr3(s*rho,atoms%rmsh(:,itype),atoms%dx(itype),atoms%jri(itype),l_proj(l))
                 rho=rho-s*l_proj(l)
            enddo !l-loop
            CALL intgr3(abs(rho),atoms%rmsh(:,itype),atoms%dx(itype),atoms%jri(itype),norm)
            write(oUnit,"(a8,i3,i2,i5,5(e12.4,1x))") 'l-dist',layer,jspin,itype,l_proj,norm,norm1
            write(oUnit,"(a8,i3,i2,i5,4(f7.4,1x))") 'l-rdist',layer,jspin,itype,l_proj/sum(l_proj),norm1/norm

        enddo !atoms
      enddo !jspins

      end subroutine
      END                                           
