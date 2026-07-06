!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_fleur_hsetup 
      use m_juDFT
      !<-- some variables are created only once per spin and            
                                                                        
      !   saved for later use                                           
      !Put these into some types                                        
      USE m_gf_types,ONLY:t_tlmplm,t_raddata 
      IMPLICIT NONE
                                                                        
      !>                                                                
      CONTAINS 
      !<-- S: fleur_tlmplm(jspin,jspins,,,,,,,,                         
                                                                        
      SUBROUTINE fleur_tlmplm(jspin,jspins,sphhar,atoms,enpara,mpi,vr   &
     &     ,tlmplm_DATA,raddata,js,vs_mmp)                              
!-----------------------------------------------                        
!  interface to tlmplm from FLEUR                                       
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_tlmplm 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
                                                                        
      INTEGER,INTENT(IN)          :: jspin,jspins,js 
      COMPLEX,INTENT(IN)          :: vs_mmp(:,:,:) 
      TYPE(t_sphhar),INTENT(IN)   :: sphhar 
      TYPE(t_atoms),INTENT(IN)    :: atoms 
      TYPE(t_enpara),INTENT(IN)   :: enpara 
      TYPE(t_mpi),INTENT(IN)      :: mpi 
      REAL,    INTENT (IN)        :: vr(:,0:,:) 
                                                   !must be INOUT becaus
      TYPE(t_tlmplm),INTENT(INOUT)  :: tlmplm_DATA 
                                                   !allocated before!   
      TYPE(t_raddata),INTENT(INOUT) :: raddata 
                                                                        
      !>                                                                
      !<-- Locals                                                       
                                                                        
      INTEGER             :: lmd,lmplmd,loplod,irecl0 
                                                                        
      !>                                                                
      loplod = Atoms%nlod*(Atoms%nlod+1)/2 
      lmd = MAXVAL(Atoms%lmax0)*(MAXVAL(Atoms%lmax0)+2) 
      lmplmd = (lmd*(lmd+3))/2 
                                                                        
!      irecl0       = 8* (8* (lmplmd+1)+ (lmd+1)**2)
!      OPEN (28,FILE='tmat',ACCESS='direct',RECL=Irecl0)
!      IF ( Jspins==2 ) OPEN (38,FILE ='tmas',STATUS='unknown',        &
!     &        ACCESS='direct',RECL= Irecl0)
                                                                        
      CALL TLMPLM(SIZE(Sphhar%clnu,1)                                   &
     &     ,MAXVAL(Sphhar%nlh),MAXVAL(Atoms%ntypsy),Atoms%ntype         &
     &     ,MAXVAL(Atoms%jri),MAXVAL(Atoms%lmax0),Jspins,Atoms%ntype    &
     &     ,Atoms%dx,Atoms%rmsh,Atoms%jri,Atoms%lmax0,Atoms%ntypsy,1    &
     &     ,Atoms%nat,atoms%lnonsph,lmd,lmplmd,Sphhar%clnu,Sphhar%mlh   &
     &     ,Sphhar%nmem,Sphhar%llh,Sphhar%nlh,Atoms%neq,jspin,jspins,1,0       &
     &     ,mpi%Irank,SIZE(tlmplm_DATA%tuulo,3),SIZE(tlmplm_DATA%tuloulo&
     &     ,3),vr(1:,0:,1:),vr(1:,0,1:),1,1,RESHAPE(enpara%el,(/        &
     &      SIZE(enpara%el,1),SIZE(enpara%el,2),SIZE(enpara%el,3),1/)), &
     &     .FALSE.,Atoms%nlod,atoms%llod,loplod,enpara%ello,Atoms%llo   &
     &     ,Atoms%nlo,Atoms%lo1l,atoms%l_dulo,atoms%ulo_der,.FALSE.     &
     &     ,(SIZE(vs_mmp,1)-1)/2,SIZE(vs_mmp,3),atoms%lda_u%l,vs_mmp    &
     &     ,tlmplm_DATA,raddata)
!      CLOSE(28)
!      CLOSE(38)
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: fleur_usetup(jspins,atoms,mpi,sphhar,vr,vs_mmp)           
                                                                        
      SUBROUTINE fleur_usetup(jspins,atoms,mpi,sphhar,enpara,vr,vs_mmp,layer)
!-----------------------------------------------                        
!  interface to u_setup from fleur                                      
!           (last modified: 2004-00-00) D. Wortmann                     
!-----------------------------------------------                        
      USE m_usetup 
      USE m_gf_types 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)          :: jspins,layer
      TYPE(T_atoms),INTENT(IN)    :: atoms 
      TYPE(t_mpi),INTENT(IN)      :: mpi 
      TYPE(t_sphhar),INTENT(IN)   :: sphhar 
      TYPE(t_enpara),INTENT(IN)   :: enpara 
      REAL,    INTENT (IN)        :: vr(:,0:,:,:) 
      COMPLEX,INTENT (OUT)        ::vs_mmp(:,:,:,:) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER             :: n_u_l 
      INTEGER, PARAMETER  :: LMAXB = 3 
      REAL                :: e_ldau 
      !>                                                                
      n_u_l = atoms%n_u 
      CALL U_SETUP(Jspins,Atoms%ntype,MAXVAL(Atoms%lmax0),              &
     &           LMAXB,Maxval(Atoms%jri),Maxval(Sphhar%nlh),            &
     &           1,n_u_l,                                               &
     &           Jspins,atoms%lda_u,Atoms%neq,enpara%el(0:,1:,1:),vr,   &
     &           Atoms%jri,Atoms%dx,Atoms%rmsh,mpi%irank,               &
     &           vs_mmp,e_ldau,layer)
      IF (n_u_l /= atoms%n_u)  THEN
	       write(6,*) "LDA+U setup:",n_u_l,atoms%n_u
               CALL juDFT_error("Setup of LDA+U incorrect"&
     &     ,calledby ="fleur_hsetup.F90")
      ENDIF
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      !<-- S: fleur_hssphn(jspins,,,)                                   
                                                                        
      SUBROUTINE fleur_hssphn(jspins,jspin,atoms,sphhar,sym,cell        &
     &     ,mpi,soc,noco,lapw,enpara,bk,vr,vs_mmp,                      &
     &     tlmplm_data,raddata,                                         &
     &     h,s)                                                         
!-----------------------------------------------                        
!interface to hssphn from FLEUR                                         
!           (last modified: 05-03-24) D. Wortmann                       
!-----------------------------------------------                        
      USE m_gf_types 
      USE m_hsmt,ONLY    : hsmt
      USE m_od_types, ONLY : od_inp, od_sym 
      IMPLICIT NONE 
      !<--Arguments                                                     
      INTEGER,INTENT(IN)            :: jspins,jspin 
      TYPE(t_atoms),INTENT(IN)      :: atoms 
      TYPE(t_sphhar),INTENT(IN)     :: sphhar 
      TYPE(t_sym),INTENT(IN)        :: sym 
      TYPE(t_cell),INTENT(IN)       :: cell 
      TYPE(t_mpi),INTENT(IN)        :: mpi 
      TYPE(t_soc),INTENT(IN)        :: soc 
      TYPE(t_noco),INTENT(IN)       :: noco 
      TYPE(t_lapw),INTENT(INOUT)    :: lapw 
      TYPE(t_enpara),INTENT(IN)     :: enpara 
      REAL, INTENT(IN)              :: Bk(3) 
      REAL, INTENT(IN)              :: vr(:,0:,:,:) 
      COMPLEX,INTENT(IN)            :: vs_mmp(:,:,:,:) 
      TYPE(t_tlmplm),INTENT(INOUT)  :: tlmplm_DATA 
      TYPE(t_raddata),INTENT(INOUT) :: raddata 
      COMPLEX,INTENT(INOUT)         :: h(:),s(:) 
      INTEGER                       :: kveclo(atoms%nlotot) 
      !>                                                                
      !<-- Locals                                                       
      INTEGER,PARAMETER             :: LMAXB = 3 
      INTEGER                       :: loplod,matsize,lmd,lmplmd 
                                                 !no internat timing    
      REAL                          :: timedummy(4)
      LOGICAL                       :: soc_opt(atoms%ntype+2) 
!-odim                                                                  
      TYPE (od_inp) :: odi 
      TYPE (od_sym) :: ods 
!+odim                                                                  
      !>                                                                
      lmd = Maxval(Atoms%lmax0)*(Maxval(Atoms%lmax0)+2) 
      lmplmd = (lmd*(lmd+3))/2 
      loplod = Atoms%nlod*(Atoms%nlod+1)/2 
      Matsize = Lapw%nbasfcn*(Lapw%nbasfcn+1)/2 
                                                                        
      odi%d1 = .FALSE. 
      soc_opt = soc%soc 
                                                                        
      CALL hsmt(                                                              &
     &      lapw%nvd,maxval(atoms%lmax0),atoms%ntype,sphhar%nlhd,             &
     &      maxval(atoms%jri),atoms%nat,                                       &
     &      sym%nop,jspins,1,mpi%SUB_COMM,mpi%n_size,mpi%n_rank,               &
     &      atoms%nlotot,                                                      &
     &      jspin,jspins,1,atoms%ntype,matsize,mpi%irank,lmd,lmplmd,           &
     &      atoms%llod,atoms%nlod,loplod,atoms%n_u,lmaxb,                     &
     &      SIZE(tlmplm_data%tuulo,3),SIZE(tlmplm_data%tuloulo,3),0,          &
     &      .false.,noco%l_noco,noco%l_constr,noco%l_ss,soc%soc,cell%omtil,   &
     &      atoms%lnonsph,atoms%lda_u%l,atoms%lapw_l,atoms%invsat,atoms%nlo,    &
     &      atoms%llo,sym%invtab,atoms%ulo_der,lapw%nv,                       &
     &      sym%mrot,atoms%jri,lapw%k%k1,lapw%k%k2,lapw%k%k3,atoms%lmax,      &
     &      atoms%neq,atoms%ngopr,bk,lapw%k%rk,atoms%rmt,atoms%taual,    &
     &      atoms%rmsh,atoms%dx,vr,cell%bmat,enpara%el,enpara%ello,           &
     &      noco%alph,noco%beta,noco%qss,noco%b_con,atoms%l_dulo,soc_opt,     &
     &      vs_mmp,ods,odi,raddata,kveclo,h,s,tlmplm_data,lapw%nmat,          &
     &      timedummy)
                                                                        
      lapw%kveclo=kveclo 
                                                                        
      END SUBROUTINE 
                                                                        
      !>                                                                
                                                                        
      END                                           
